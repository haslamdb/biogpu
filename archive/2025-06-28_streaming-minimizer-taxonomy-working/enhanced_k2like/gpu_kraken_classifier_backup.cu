// gpu_kraken_classifier_paired.cuh
// GPU-accelerated Kraken2-style taxonomic classifier with PAIRED-END support
// Implements paired-end concordance, voting boost, and confidence scoring

#pragma once
#ifndef GPU_KRAKEN_CLASSIFIER_PAIRED_CUH
#define GPU_KRAKEN_CLASSIFIER_PAIRED_CUH

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>

// Paired read structure (similar to your existing profiler)
struct PairedRead {
    std::string read1;
    std::string read2;
    std::string read_id;
    bool is_paired;
    
    PairedRead(const std::string& r1, const std::string& r2 = "", 
               const std::string& id = "") 
        : read1(r1), read2(r2), read_id(id), is_paired(!r2.empty()) {}
};

// GPU Compact Hash Table (same as before)
struct GPUCompactHashTable {
    uint32_t* hash_cells;
    uint32_t table_size;
    uint32_t hash_mask;
    uint32_t lca_bits;
    uint32_t hash_bits;
};

// Classification parameters with paired-end options
struct ClassificationParams {
    int k = 35;
    int ell = 31;
    int spaces = 7;
    float confidence_threshold = 0.0f;
    bool use_spaced_seeds = true;
    int max_ambiguous_bases = 5;
    
    // Paired-end specific parameters
    bool use_paired_end_bonus = true;      // Boost for concordant pairs
    float paired_concordance_weight = 2.0f; // Weight multiplier for concordant pairs
    float min_pair_concordance = 0.5f;     // Minimum concordance to get bonus
    bool require_both_reads_classified = false; // Both reads must classify to get bonus
};

// Enhanced classification result with paired-end metrics
struct PairedReadClassification {
    uint32_t taxon_id;                 // Final classification (0 = unclassified)
    float confidence_score;            // Overall confidence [0,1]
    
    // Per-read metrics
    uint32_t read1_votes;              // Votes from read1
    uint32_t read2_votes;              // Votes from read2  
    uint32_t read1_kmers;              // Total k-mers in read1
    uint32_t read2_kmers;              // Total k-mers in read2
    
    // Paired-end metrics
    uint32_t concordant_votes;         // Votes where both reads agree
    float pair_concordance;            // Concordance rate [0,1]
    bool got_paired_bonus;             // Whether paired bonus was applied
    
    // Individual read classifications (for debugging)
    uint32_t read1_best_taxon;
    uint32_t read2_best_taxon;
    float read1_confidence;
    float read2_confidence;
};

// Taxonomy tree node
struct TaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    char name[64];
};

// GPU-accelerated Paired-End Kraken classifier
class PairedEndGPUKrakenClassifier {
private:
    // GPU data structures
    GPUCompactHashTable* d_hash_table;
    TaxonomyNode* d_taxonomy_tree;
    uint32_t* d_parent_lookup;
    uint32_t num_taxonomy_nodes;
    
    // Host taxonomy mapping
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    
    // Classification parameters
    ClassificationParams params;
    
    // GPU memory for paired-end batch processing
    char* d_reads_r1;
    char* d_reads_r2;
    uint32_t* d_read_offsets_r1;
    uint32_t* d_read_offsets_r2;
    uint32_t* d_read_lengths_r1;
    uint32_t* d_read_lengths_r2;
    bool* d_is_paired_flags;
    PairedReadClassification* d_results;
    
    // Paired-end voting arrays (per pair)
    uint32_t* d_pair_votes_r1;         // [pair_id * max_taxa + taxon_idx] = votes
    uint32_t* d_pair_votes_r2;
    uint32_t* d_concordant_votes;      // [pair_id * max_taxa + taxon_idx] = concordant votes
    
    // Configuration
    static const int DEFAULT_BATCH_SIZE = 10000;
    static const int MAX_READ_LENGTH = 1000;
    static const int THREADS_PER_BLOCK = 256;
    static const int MAX_TAXA_PER_PAIR = 128;  // Max taxa to track per pair
    
    bool database_loaded = false;
    
public:
    PairedEndGPUKrakenClassifier(const ClassificationParams& config = ClassificationParams());
    ~PairedEndGPUKrakenClassifier();
    
    // Database loading
    bool load_database(const std::string& database_directory);
    
    // Paired-end classification interface
    std::vector<PairedReadClassification> classify_paired_reads(
        const std::vector<PairedRead>& paired_reads
    );
    
    std::vector<PairedReadClassification> classify_paired_reads_batch(
        const std::vector<PairedRead>& paired_reads,
        int batch_size = DEFAULT_BATCH_SIZE
    );
    
    // Single-end compatibility (converts to unpaired PairedRead)
    std::vector<PairedReadClassification> classify_reads(
        const std::vector<std::string>& reads
    );
    
    // Configuration
    void set_confidence_threshold(float threshold) {
        params.confidence_threshold = threshold;
    }
    void set_paired_concordance_weight(float weight) {
        params.paired_concordance_weight = weight;
    }
    void enable_paired_end_bonus(bool enable) {
        params.use_paired_end_bonus = enable;
    }
    ClassificationParams get_params() const { return params; }
    
    // Utilities
    std::string get_taxon_name(uint32_t taxon_id) const;
    void print_database_stats() const;
    void print_paired_classification_stats(const std::vector<PairedReadClassification>& results) const;
    
    bool is_database_loaded() const { return database_loaded; }
    
private:
    // Internal methods
    bool allocate_gpu_memory(int max_pairs);
    void free_gpu_memory();
    
    bool load_hash_table(const std::string& hash_table_file);
    bool load_taxonomy_tree(const std::string& taxonomy_file);
    bool load_taxonomy_tree_from_db(const std::string& db_file);
    bool read_database_header(const std::string& db_file, ClassificationParams& db_params);
    
    void transfer_paired_reads_to_gpu(const std::vector<PairedRead>& paired_reads);
    void retrieve_paired_results_from_gpu(std::vector<PairedReadClassification>& results, int num_pairs);
};

// CUDA kernels for paired-end processing
__global__ void classify_paired_reads_kernel(
    const char* reads_r1,
    const char* reads_r2,
    const uint32_t* read_offsets_r1,
    const uint32_t* read_offsets_r2,
    const uint32_t* read_lengths_r1,
    const uint32_t* read_lengths_r2,
    const bool* is_paired_flags,
    const GPUCompactHashTable* hash_table,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup,
    uint32_t* pair_votes_r1,
    uint32_t* pair_votes_r2,
    uint32_t* concordant_votes,
    PairedReadClassification* results,
    int num_pairs,
    ClassificationParams params
);

__global__ void compute_paired_concordance_kernel(
    const uint32_t* pair_votes_r1,
    const uint32_t* pair_votes_r2,
    const bool* is_paired_flags,
    uint32_t* concordant_votes,
    PairedReadClassification* results,
    int num_pairs,
    ClassificationParams params
);

// Device functions
__device__ uint32_t lookup_lca_gpu(const GPUCompactHashTable* cht, uint64_t minimizer_hash);
__device__ uint64_t extract_minimizer_with_spaced_seeds(const char* sequence, int pos, int k, int ell, int spaces);
__device__ uint32_t find_paired_classification_with_confidence(
    const uint32_t* votes_r1,
    const uint32_t* votes_r2,
    const uint32_t* concordant_votes,
    int total_kmers_r1,
    int total_kmers_r2,
    float confidence_threshold,
    float paired_weight,
    bool use_paired_bonus,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup
);

// Utility functions
__device__ __host__ uint64_t hash_minimizer(const char* seq, int len);
__device__ __host__ uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces);
__device__ __host__ bool has_ambiguous_bases(const char* seq, int len);
__device__ __host__ uint32_t compute_compact_hash(uint64_t minimizer_hash);
__device__ __host__ uint32_t jenkins_hash(uint64_t key);

#endif // GPU_KRAKEN_CLASSIFIER_PAIRED_CUH

// ================================================================
// IMPLEMENTATION
// ================================================================

#ifndef GPU_KRAKEN_CLASSIFIER_HEADER_ONLY
// Implementation continues below - no need to include self
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <sstream>

#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// Constructor
PairedEndGPUKrakenClassifier::PairedEndGPUKrakenClassifier(const ClassificationParams& config)
    : params(config), d_hash_table(nullptr), d_taxonomy_tree(nullptr),
      d_parent_lookup(nullptr), num_taxonomy_nodes(0),
      d_reads_r1(nullptr), d_reads_r2(nullptr),
      d_read_offsets_r1(nullptr), d_read_offsets_r2(nullptr),
      d_read_lengths_r1(nullptr), d_read_lengths_r2(nullptr),
      d_is_paired_flags(nullptr), d_results(nullptr),
      d_pair_votes_r1(nullptr), d_pair_votes_r2(nullptr), d_concordant_votes(nullptr) {
    
    std::cout << "Initializing Paired-End GPU Kraken classifier..." << std::endl;
    std::cout << "Parameters: k=" << params.k << ", ell=" << params.ell 
              << ", spaces=" << params.spaces 
              << ", confidence=" << params.confidence_threshold << std::endl;
    std::cout << "Paired-end: bonus=" << (params.use_paired_end_bonus ? "ON" : "OFF")
              << ", weight=" << params.paired_concordance_weight << std::endl;
}

// Destructor (similar to single-end version)
PairedEndGPUKrakenClassifier::~PairedEndGPUKrakenClassifier() {
    free_gpu_memory();
    
    // Free database structures
    if (d_hash_table) {
        GPUCompactHashTable h_cht;
        CUDA_CHECK(cudaMemcpy(&h_cht, d_hash_table, sizeof(GPUCompactHashTable), 
                             cudaMemcpyDeviceToHost));
        if (h_cht.hash_cells) cudaFree(h_cht.hash_cells);
        cudaFree(d_hash_table);
    }
    if (d_taxonomy_tree) cudaFree(d_taxonomy_tree);
    if (d_parent_lookup) cudaFree(d_parent_lookup);
}

// Classify paired reads with batch processing
std::vector<PairedReadClassification> PairedEndGPUKrakenClassifier::classify_paired_reads_batch(
    const std::vector<PairedRead>& paired_reads,
    int batch_size) {
    
    if (!database_loaded) {
        throw std::runtime_error("Database not loaded. Call load_database() first.");
    }
    
    std::vector<PairedReadClassification> all_results;
    all_results.reserve(paired_reads.size());
    
    std::cout << "Classifying " << paired_reads.size() << " read pairs in batches of " 
              << batch_size << "..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (size_t batch_start = 0; batch_start < paired_reads.size(); batch_start += batch_size) {
        size_t batch_end = std::min(batch_start + batch_size, paired_reads.size());
        
        // Extract batch
        std::vector<PairedRead> batch_pairs(
            paired_reads.begin() + batch_start,
            paired_reads.begin() + batch_end
        );
        
        // Process batch
        auto batch_results = classify_paired_reads(batch_pairs);
        
        // Accumulate results
        all_results.insert(all_results.end(), 
                          batch_results.begin(), batch_results.end());
        
        // Report progress every 10 batches or at the end
        if ((batch_end - batch_size) % (batch_size * 10) == 0 || batch_end == paired_reads.size()) {
            std::cout << "Processed " << batch_end << "/" << paired_reads.size() 
                      << " read pairs..." << std::endl;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    double pairs_per_second = paired_reads.size() * 1000.0 / duration.count();
    std::cout << "Paired-end classification completed in " << duration.count() << " ms" << std::endl;
    std::cout << "Performance: " << std::fixed << std::setprecision(0) 
              << pairs_per_second << " pairs/second" << std::endl;
    
    return all_results;
}

// Classify single batch of paired reads
std::vector<PairedReadClassification> PairedEndGPUKrakenClassifier::classify_paired_reads(
    const std::vector<PairedRead>& paired_reads) {
    
    if (paired_reads.empty()) return {};
    
    int num_pairs = paired_reads.size();
    
    // Allocate GPU memory for this batch
    if (!allocate_gpu_memory(num_pairs)) {
        throw std::runtime_error("Failed to allocate GPU memory");
    }
    
    // Transfer paired reads to GPU
    transfer_paired_reads_to_gpu(paired_reads);
    
    // Reset voting arrays
    size_t vote_array_size = num_pairs * MAX_TAXA_PER_PAIR * sizeof(uint32_t);
    CUDA_CHECK(cudaMemset(d_pair_votes_r1, 0, vote_array_size));
    CUDA_CHECK(cudaMemset(d_pair_votes_r2, 0, vote_array_size));
    CUDA_CHECK(cudaMemset(d_concordant_votes, 0, vote_array_size));
    
    // Launch paired-end classification kernel
    int num_blocks = (num_pairs + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    
    classify_paired_reads_kernel<<<num_blocks, THREADS_PER_BLOCK>>>(
        d_reads_r1, d_reads_r2,
        d_read_offsets_r1, d_read_offsets_r2,
        d_read_lengths_r1, d_read_lengths_r2,
        d_is_paired_flags, d_hash_table, d_taxonomy_tree, d_parent_lookup,
        d_pair_votes_r1, d_pair_votes_r2, d_concordant_votes,
        d_results, num_pairs, params
    );
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Compute paired concordance and final classification
    compute_paired_concordance_kernel<<<num_blocks, THREADS_PER_BLOCK>>>(
        d_pair_votes_r1, d_pair_votes_r2, d_is_paired_flags,
        d_concordant_votes, d_results, num_pairs, params
    );
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Retrieve results
    std::vector<PairedReadClassification> results(num_pairs);
    retrieve_paired_results_from_gpu(results, num_pairs);
    
    return results;
}

// Allocate GPU memory for paired-end processing
bool PairedEndGPUKrakenClassifier::allocate_gpu_memory(int max_pairs) {
    // Free any existing memory
    free_gpu_memory();
    
    size_t read_buffer_size = max_pairs * MAX_READ_LENGTH;
    size_t vote_array_size = max_pairs * MAX_TAXA_PER_PAIR * sizeof(uint32_t);
    
    CUDA_CHECK(cudaMalloc(&d_reads_r1, read_buffer_size));
    CUDA_CHECK(cudaMalloc(&d_reads_r2, read_buffer_size));
    CUDA_CHECK(cudaMalloc(&d_read_offsets_r1, max_pairs * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_read_offsets_r2, max_pairs * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_read_lengths_r1, max_pairs * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_read_lengths_r2, max_pairs * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_is_paired_flags, max_pairs * sizeof(bool)));
    CUDA_CHECK(cudaMalloc(&d_results, max_pairs * sizeof(PairedReadClassification)));
    
    // Paired-end voting arrays
    CUDA_CHECK(cudaMalloc(&d_pair_votes_r1, vote_array_size));
    CUDA_CHECK(cudaMalloc(&d_pair_votes_r2, vote_array_size));
    CUDA_CHECK(cudaMalloc(&d_concordant_votes, vote_array_size));
    
    return true;
}

// Free GPU memory
void PairedEndGPUKrakenClassifier::free_gpu_memory() {
    if (d_reads_r1) { cudaFree(d_reads_r1); d_reads_r1 = nullptr; }
    if (d_reads_r2) { cudaFree(d_reads_r2); d_reads_r2 = nullptr; }
    if (d_read_offsets_r1) { cudaFree(d_read_offsets_r1); d_read_offsets_r1 = nullptr; }
    if (d_read_offsets_r2) { cudaFree(d_read_offsets_r2); d_read_offsets_r2 = nullptr; }
    if (d_read_lengths_r1) { cudaFree(d_read_lengths_r1); d_read_lengths_r1 = nullptr; }
    if (d_read_lengths_r2) { cudaFree(d_read_lengths_r2); d_read_lengths_r2 = nullptr; }
    if (d_is_paired_flags) { cudaFree(d_is_paired_flags); d_is_paired_flags = nullptr; }
    if (d_results) { cudaFree(d_results); d_results = nullptr; }
    if (d_pair_votes_r1) { cudaFree(d_pair_votes_r1); d_pair_votes_r1 = nullptr; }
    if (d_pair_votes_r2) { cudaFree(d_pair_votes_r2); d_pair_votes_r2 = nullptr; }
    if (d_concordant_votes) { cudaFree(d_concordant_votes); d_concordant_votes = nullptr; }
}

// Transfer paired reads to GPU
void PairedEndGPUKrakenClassifier::transfer_paired_reads_to_gpu(const std::vector<PairedRead>& paired_reads) {
    // Prepare data for both reads
    std::string concatenated_r1, concatenated_r2;
    std::vector<uint32_t> offsets_r1, offsets_r2;
    std::vector<uint32_t> lengths_r1, lengths_r2;
    std::vector<uint8_t> is_paired_flags;
    
    uint32_t current_offset_r1 = 0, current_offset_r2 = 0;
    
    for (const auto& pair : paired_reads) {
        // Read 1
        offsets_r1.push_back(current_offset_r1);
        lengths_r1.push_back(pair.read1.length());
        concatenated_r1 += pair.read1;
        current_offset_r1 += pair.read1.length();
        
        // Read 2 (or empty if single-end)
        offsets_r2.push_back(current_offset_r2);
        if (pair.is_paired && !pair.read2.empty()) {
            lengths_r2.push_back(pair.read2.length());
            concatenated_r2 += pair.read2;
            current_offset_r2 += pair.read2.length();
        } else {
            lengths_r2.push_back(0);
            // Add empty placeholder to maintain alignment
        }
        
        is_paired_flags.push_back(pair.is_paired);
    }
    
    // Transfer to GPU
    CUDA_CHECK(cudaMemcpy(d_reads_r1, concatenated_r1.c_str(),
                         concatenated_r1.length(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_reads_r2, concatenated_r2.c_str(),
                         concatenated_r2.length(), cudaMemcpyHostToDevice));
    
    CUDA_CHECK(cudaMemcpy(d_read_offsets_r1, offsets_r1.data(),
                         offsets_r1.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_read_offsets_r2, offsets_r2.data(),
                         offsets_r2.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    CUDA_CHECK(cudaMemcpy(d_read_lengths_r1, lengths_r1.data(),
                         lengths_r1.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_read_lengths_r2, lengths_r2.data(),
                         lengths_r2.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    CUDA_CHECK(cudaMemcpy(d_is_paired_flags, is_paired_flags.data(),
                         is_paired_flags.size() * sizeof(uint8_t), cudaMemcpyHostToDevice));
}

// Single-end compatibility wrapper
std::vector<PairedReadClassification> PairedEndGPUKrakenClassifier::classify_reads(
    const std::vector<std::string>& reads) {
    
    // Convert to paired format
    std::vector<PairedRead> paired_reads;
    paired_reads.reserve(reads.size());
    
    for (size_t i = 0; i < reads.size(); i++) {
        paired_reads.emplace_back(reads[i], "", "read_" + std::to_string(i));
    }
    
    return classify_paired_reads_batch(paired_reads);
}

// Print paired-end classification statistics
void PairedEndGPUKrakenClassifier::print_paired_classification_stats(
    const std::vector<PairedReadClassification>& results) const {
    
    if (results.empty()) {
        std::cout << "No classification results to analyze" << std::endl;
        return;
    }
    
    int classified = 0;
    int unclassified = 0;
    int paired_bonus_count = 0;
    float total_confidence = 0.0f;
    float total_concordance = 0.0f;
    int paired_count = 0;
    
    std::unordered_map<uint32_t, int> taxon_counts;
    
    for (const auto& result : results) {
        if (result.taxon_id > 0) {
            classified++;
            total_confidence += result.confidence_score;
            taxon_counts[result.taxon_id]++;
        } else {
            unclassified++;
        }
        
        if (result.got_paired_bonus) {
            paired_bonus_count++;
        }
        
        // Check if this was a paired read
        if (result.read2_kmers > 0) {
            paired_count++;
            total_concordance += result.pair_concordance;
        }
    }
    
    std::cout << "\n=== PAIRED-END CLASSIFICATION STATISTICS ===" << std::endl;
    std::cout << "Total read pairs: " << results.size() << std::endl;
    std::cout << "Paired reads: " << paired_count << std::endl;
    std::cout << "Single reads: " << (results.size() - paired_count) << std::endl;
    std::cout << "Classified: " << classified 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * classified / results.size()) << "%)" << std::endl;
    std::cout << "Unclassified: " << unclassified 
              << " (" << std::fixed << std::setprecision(1) 
              << (100.0 * unclassified / results.size()) << "%)" << std::endl;
    
    if (classified > 0) {
        std::cout << "Average confidence: " << std::fixed << std::setprecision(3) 
                  << (total_confidence / classified) << std::endl;
        std::cout << "Unique taxa detected: " << taxon_counts.size() << std::endl;
    }
    
    if (paired_count > 0) {
        std::cout << "Average pair concordance: " << std::fixed << std::setprecision(3) 
                  << (total_concordance / paired_count) << std::endl;
        std::cout << "Pairs with concordance bonus: " << paired_bonus_count 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * paired_bonus_count / paired_count) << "%)" << std::endl;
    }
}

// ================================================================
// CUDA KERNELS FOR PAIRED-END PROCESSING
// ================================================================

// Main paired-end classification kernel
__global__ void classify_paired_reads_kernel(
    const char* reads_r1,
    const char* reads_r2,
    const uint32_t* read_offsets_r1,
    const uint32_t* read_offsets_r2,
    const uint32_t* read_lengths_r1,
    const uint32_t* read_lengths_r2,
    const bool* is_paired_flags,
    const GPUCompactHashTable* hash_table,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup,
    uint32_t* pair_votes_r1,
    uint32_t* pair_votes_r2,
    uint32_t* concordant_votes,
    PairedReadClassification* results,
    int num_pairs,
    ClassificationParams params) {
    
    int pair_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_id >= num_pairs) return;
    
    const char* read1 = reads_r1 + read_offsets_r1[pair_id];
    const char* read2 = reads_r2 + read_offsets_r2[pair_id];
    uint32_t len1 = read_lengths_r1[pair_id];
    uint32_t len2 = read_lengths_r2[pair_id];
    bool is_paired = is_paired_flags[pair_id];
    
    // Initialize result
    PairedReadClassification& result = results[pair_id];
    result.taxon_id = 0;
    result.confidence_score = 0.0f;
    result.read1_votes = 0;
    result.read2_votes = 0;
    result.read1_kmers = 0;
    result.read2_kmers = 0;
    result.concordant_votes = 0;
    result.pair_concordance = 0.0f;
    result.got_paired_bonus = false;
    result.read1_best_taxon = 0;
    result.read2_best_taxon = 0;
    result.read1_confidence = 0.0f;
    result.read2_confidence = 0.0f;
    
    // Local voting arrays for this pair
    const int MAX_LOCAL_TAXA = 32;
    uint32_t local_taxon_ids_r1[MAX_LOCAL_TAXA];
    uint32_t local_vote_counts_r1[MAX_LOCAL_TAXA];
    uint32_t local_taxon_ids_r2[MAX_LOCAL_TAXA];
    uint32_t local_vote_counts_r2[MAX_LOCAL_TAXA];
    int num_taxa_r1 = 0, num_taxa_r2 = 0;
    
    // Process Read 1
    if (len1 >= params.k) {
        int total_kmers_r1 = len1 - params.k + 1;
        for (int i = 0; i < total_kmers_r1; i++) {
            if (has_ambiguous_bases(read1 + i, params.k)) continue;
            
            result.read1_kmers++;
            
            uint64_t minimizer = extract_minimizer_with_spaced_seeds(
                read1, i, params.k, params.ell, params.spaces
            );
            
            if (minimizer == UINT64_MAX) continue;
            
            uint32_t lca = lookup_lca_gpu(hash_table, minimizer);
            if (lca == 0) continue;
            
            // Add vote for read1
            bool found = false;
            for (int j = 0; j < num_taxa_r1; j++) {
                if (local_taxon_ids_r1[j] == lca) {
                    local_vote_counts_r1[j]++;
                    found = true;
                    break;
                }
            }
            
            if (!found && num_taxa_r1 < MAX_LOCAL_TAXA) {
                local_taxon_ids_r1[num_taxa_r1] = lca;
                local_vote_counts_r1[num_taxa_r1] = 1;
                num_taxa_r1++;
            }
        }
    }
    
    // Process Read 2 (if paired)
    if (is_paired && len2 >= params.k) {
        int total_kmers_r2 = len2 - params.k + 1;
        for (int i = 0; i < total_kmers_r2; i++) {
            if (has_ambiguous_bases(read2 + i, params.k)) continue;
            
            result.read2_kmers++;
            
            uint64_t minimizer = extract_minimizer_with_spaced_seeds(
                read2, i, params.k, params.ell, params.spaces
            );
            
            if (minimizer == UINT64_MAX) continue;
            
            uint32_t lca = lookup_lca_gpu(hash_table, minimizer);
            if (lca == 0) continue;
            
            // Add vote for read2
            bool found = false;
            for (int j = 0; j < num_taxa_r2; j++) {
                if (local_taxon_ids_r2[j] == lca) {
                    local_vote_counts_r2[j]++;
                    found = true;
                    break;
                }
            }
            
            if (!found && num_taxa_r2 < MAX_LOCAL_TAXA) {
                local_taxon_ids_r2[num_taxa_r2] = lca;
                local_vote_counts_r2[num_taxa_r2] = 1;
                num_taxa_r2++;
            }
        }
    }
    
    // Find best taxon for each read
    uint32_t best_taxon_r1 = 0, best_votes_r1 = 0;
    uint32_t best_taxon_r2 = 0, best_votes_r2 = 0;
    
    for (int i = 0; i < num_taxa_r1; i++) {
        if (local_vote_counts_r1[i] > best_votes_r1) {
            best_votes_r1 = local_vote_counts_r1[i];
            best_taxon_r1 = local_taxon_ids_r1[i];
        }
    }
    
    for (int i = 0; i < num_taxa_r2; i++) {
        if (local_vote_counts_r2[i] > best_votes_r2) {
            best_votes_r2 = local_vote_counts_r2[i];
            best_taxon_r2 = local_taxon_ids_r2[i];
        }
    }
    
    result.read1_best_taxon = best_taxon_r1;
    result.read2_best_taxon = best_taxon_r2;
    result.read1_votes = best_votes_r1;
    result.read2_votes = best_votes_r2;
    result.read1_confidence = (result.read1_kmers > 0) ? 
        (float)best_votes_r1 / result.read1_kmers : 0.0f;
    result.read2_confidence = (result.read2_kmers > 0) ? 
        (float)best_votes_r2 / result.read2_kmers : 0.0f;
    
    // Calculate concordance and final classification
    if (is_paired && best_taxon_r1 > 0 && best_taxon_r2 > 0) {
        // Check if both reads agree
        if (best_taxon_r1 == best_taxon_r2) {
            result.concordant_votes = best_votes_r1 + best_votes_r2;
            result.pair_concordance = 1.0f;  // Perfect concordance
            
            // Apply paired-end bonus
            if (params.use_paired_end_bonus && 
                result.pair_concordance >= params.min_pair_concordance) {
                result.taxon_id = best_taxon_r1;
                result.confidence_score = (result.read1_confidence + result.read2_confidence) * 
                                        params.paired_concordance_weight / 2.0f;
                result.got_paired_bonus = true;
            } else {
                result.taxon_id = best_taxon_r1;
                result.confidence_score = (result.read1_confidence + result.read2_confidence) / 2.0f;
            }
        } else {
            // Reads disagree - use better one or fall back to single-end logic
            if (result.read1_confidence > result.read2_confidence) {
                result.taxon_id = best_taxon_r1;
                result.confidence_score = result.read1_confidence;
            } else {
                result.taxon_id = best_taxon_r2;
                result.confidence_score = result.read2_confidence;
            }
            result.pair_concordance = 0.0f;
        }
    } else {
        // Single-end classification or only one read classified
        if (best_taxon_r1 > 0 && best_taxon_r2 == 0) {
            result.taxon_id = best_taxon_r1;
            result.confidence_score = result.read1_confidence;
        } else if (best_taxon_r2 > 0 && best_taxon_r1 == 0) {
            result.taxon_id = best_taxon_r2;
            result.confidence_score = result.read2_confidence;
        } else if (best_taxon_r1 > 0 && best_taxon_r2 > 0) {
            // Both classified but not paired, use better one
            if (result.read1_confidence > result.read2_confidence) {
                result.taxon_id = best_taxon_r1;
                result.confidence_score = result.read1_confidence;
            } else {
                result.taxon_id = best_taxon_r2;
                result.confidence_score = result.read2_confidence;
            }
        }
    }
    
    // Apply confidence threshold
    if (result.confidence_score < params.confidence_threshold) {
        result.taxon_id = 0;  // Unclassified
    }
}

// Compute paired concordance kernel (additional processing if needed)
__global__ void compute_paired_concordance_kernel(
    const uint32_t* pair_votes_r1,
    const uint32_t* pair_votes_r2,
    const bool* is_paired_flags,
    uint32_t* concordant_votes,
    PairedReadClassification* results,
    int num_pairs,
    ClassificationParams params) {
    
    int pair_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_id >= num_pairs) return;
    
    // This kernel can be used for additional paired-end processing
    // if needed, such as more sophisticated LCA computation across
    // read pairs. For now, the main kernel handles everything.
}

// Retrieve paired results from GPU
void PairedEndGPUKrakenClassifier::retrieve_paired_results_from_gpu(
    std::vector<PairedReadClassification>& results, int num_pairs) {
    
    CUDA_CHECK(cudaMemcpy(results.data(), d_results,
                         num_pairs * sizeof(PairedReadClassification),
                         cudaMemcpyDeviceToHost));
}

// Database loading methods
bool PairedEndGPUKrakenClassifier::load_database(const std::string& database_directory) {
    std::cout << "Loading database from " << database_directory << "..." << std::endl;
    
    // Load hash table
    std::string hash_table_file = database_directory + "/hash_table.k2d";
    if (!load_hash_table(hash_table_file)) {
        std::cerr << "Failed to load hash table from " << hash_table_file << std::endl;
        return false;
    }
    
    // Load taxonomy
    std::string taxonomy_file = database_directory + "/taxonomy.tsv";
    if (!load_taxonomy_tree(taxonomy_file)) {
        std::cerr << "Failed to load taxonomy from " << taxonomy_file << std::endl;
        return false;
    }
    
    database_loaded = true;
    std::cout << "Database loaded successfully!" << std::endl;
    return true;
}

bool PairedEndGPUKrakenClassifier::load_hash_table(const std::string& hash_table_file) {
    std::ifstream hash_in(hash_table_file, std::ios::binary);
    if (!hash_in.is_open()) {
        std::cerr << "Cannot open hash table file: " << hash_table_file << std::endl;
        return false;
    }
    
    // Read header
    uint64_t table_size, num_entries;
    hash_in.read(reinterpret_cast<char*>(&table_size), sizeof(uint64_t));
    hash_in.read(reinterpret_cast<char*>(&num_entries), sizeof(uint64_t));
    
    std::cout << "Hash table size: " << table_size << ", entries: " << num_entries << std::endl;
    
    // Create compact hash table on GPU
    GPUCompactHashTable h_cht;
    h_cht.table_size = table_size;
    h_cht.hash_mask = table_size - 1;
    h_cht.lca_bits = 20;  // Assuming 20 bits for LCA (up to 1M taxa)
    h_cht.hash_bits = 32 - h_cht.lca_bits;
    
    // Allocate hash table on GPU
    CUDA_CHECK(cudaMalloc(&h_cht.hash_cells, table_size * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemset(h_cht.hash_cells, 0, table_size * sizeof(uint32_t)));
    
    // Allocate and copy hash table structure to GPU
    CUDA_CHECK(cudaMalloc(&d_hash_table, sizeof(GPUCompactHashTable)));
    CUDA_CHECK(cudaMemcpy(d_hash_table, &h_cht, sizeof(GPUCompactHashTable), cudaMemcpyHostToDevice));
    
    // Read entries and build compact hash table on host first
    std::vector<uint32_t> host_hash_cells(table_size, 0);
    
    for (uint64_t i = 0; i < num_entries; i++) {
        uint64_t minimizer_hash;
        uint32_t lca_taxon, genome_count;
        float uniqueness_score;
        
        hash_in.read(reinterpret_cast<char*>(&minimizer_hash), sizeof(uint64_t));
        hash_in.read(reinterpret_cast<char*>(&lca_taxon), sizeof(uint32_t));
        hash_in.read(reinterpret_cast<char*>(&genome_count), sizeof(uint32_t));
        hash_in.read(reinterpret_cast<char*>(&uniqueness_score), sizeof(float));
        
        // Compute compact hash
        uint32_t compact_hash = compute_compact_hash(minimizer_hash);
        uint32_t pos = compact_hash & h_cht.hash_mask;
        
        // Linear probing to find empty slot
        int probes = 0;
        while (host_hash_cells[pos] != 0 && probes < 32) {
            pos = (pos + 1) & h_cht.hash_mask;
            probes++;
        }
        
        if (probes < 32) {
            // Store compact hash (upper bits) and LCA (lower bits)
            uint32_t stored_hash = compact_hash >> h_cht.lca_bits;
            uint32_t cell_value = (stored_hash << h_cht.lca_bits) | (lca_taxon & ((1U << h_cht.lca_bits) - 1));
            host_hash_cells[pos] = cell_value;
        }
    }
    
    // Copy hash table to GPU
    CUDA_CHECK(cudaMemcpy(h_cht.hash_cells, host_hash_cells.data(), 
                         table_size * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    hash_in.close();
    std::cout << "Hash table loaded: " << num_entries << " entries" << std::endl;
    return true;
}

bool PairedEndGPUKrakenClassifier::load_taxonomy_tree_from_db(const std::string& db_file) {
    std::ifstream in(db_file, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "Cannot open database file for taxonomy: " << db_file << std::endl;
        return false;
    }
    
    // Read header to get organism count
    uint32_t magic, version, k_size, m_size, num_organisms;
    uint64_t num_minimizer_hashes;
    
    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
    in.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
    in.read(reinterpret_cast<char*>(&num_organisms), sizeof(num_organisms));
    in.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
    
    std::cout << "Loading taxonomy for " << num_organisms << " organisms..." << std::endl;
    
    // Read organism metadata and build taxonomy mapping
    for (uint32_t i = 0; i < num_organisms; i++) {
        uint32_t taxonomy_id, taxon_level, minimizer_count;
        uint64_t genome_size;
        float gc_content;
        
        in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
        in.read(reinterpret_cast<char*>(&taxon_level), sizeof(taxon_level));
        in.read(reinterpret_cast<char*>(&gc_content), sizeof(gc_content));
        in.read(reinterpret_cast<char*>(&genome_size), sizeof(genome_size));
        in.read(reinterpret_cast<char*>(&minimizer_count), sizeof(minimizer_count));
        
        // Read organism name
        uint16_t name_length;
        in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
        std::string name(name_length, '\0');
        in.read(&name[0], name_length);
        
        // Read taxonomy path
        uint16_t taxonomy_length;
        in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
        std::string taxonomy_path(taxonomy_length, '\0');
        in.read(&taxonomy_path[0], taxonomy_length);
        
        // Store in taxonomy mappings
        taxon_names[taxonomy_id] = name;
        taxon_parents[taxonomy_id] = 0;  // Simple mapping for now
        
        // Create a simple taxonomy node
        TaxonomyNode node;
        node.taxon_id = taxonomy_id;
        node.parent_id = 0;
        node.rank = taxon_level;
        strncpy(node.name, name.c_str(), 63);
        node.name[63] = '\0';
    }
    
    in.close();
    
    // Create simple parent lookup (all point to root for now)
    std::vector<uint32_t> parent_lookup(1000000, 0);  // Support up to 1M taxa
    
    // Allocate taxonomy tree on GPU (simplified)
    std::vector<TaxonomyNode> host_taxonomy;
    for (const auto& [taxon_id, name] : taxon_names) {
        TaxonomyNode node;
        node.taxon_id = taxon_id;
        node.parent_id = 0;
        node.rank = 0;
        strncpy(node.name, name.c_str(), 63);
        node.name[63] = '\0';
        host_taxonomy.push_back(node);
    }
    
    num_taxonomy_nodes = host_taxonomy.size();
    
    if (num_taxonomy_nodes > 0) {
        CUDA_CHECK(cudaMalloc(&d_taxonomy_tree, num_taxonomy_nodes * sizeof(TaxonomyNode)));
        CUDA_CHECK(cudaMemcpy(d_taxonomy_tree, host_taxonomy.data(), 
                             num_taxonomy_nodes * sizeof(TaxonomyNode), cudaMemcpyHostToDevice));
        
        CUDA_CHECK(cudaMalloc(&d_parent_lookup, 1000000 * sizeof(uint32_t)));
        CUDA_CHECK(cudaMemcpy(d_parent_lookup, parent_lookup.data(), 
                             1000000 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    }
    
    std::cout << "Loaded taxonomy for " << taxon_names.size() << " organisms" << std::endl;
    return true;
}

bool PairedEndGPUKrakenClassifier::load_taxonomy_tree(const std::string& taxonomy_file) {
    std::ifstream tax_in(taxonomy_file);
    if (!tax_in.is_open()) {
        std::cerr << "Cannot open taxonomy file: " << taxonomy_file << std::endl;
        return false;
    }
    
    std::string line;
    // Skip header
    std::getline(tax_in, line);
    
    std::vector<TaxonomyNode> host_taxonomy;
    uint32_t max_taxon_id = 0;
    
    while (std::getline(tax_in, line)) {
        std::istringstream iss(line);
        uint32_t taxon_id, parent_id;
        std::string name;
        
        if (iss >> taxon_id) {
            iss.ignore(1); // skip tab
            std::getline(iss, name, '\t');
            iss >> parent_id;
            
            taxon_names[taxon_id] = name;
            taxon_parents[taxon_id] = parent_id;
            
            TaxonomyNode node;
            node.taxon_id = taxon_id;
            node.parent_id = parent_id;
            node.rank = 0;  // Not used for now
            strncpy(node.name, name.c_str(), 63);
            node.name[63] = '\0';
            
            host_taxonomy.push_back(node);
            max_taxon_id = std::max(max_taxon_id, taxon_id);
        }
    }
    
    tax_in.close();
    
    num_taxonomy_nodes = host_taxonomy.size();
    std::cout << "Loaded " << num_taxonomy_nodes << " taxonomy nodes, max ID: " << max_taxon_id << std::endl;
    
    // Allocate and copy taxonomy to GPU
    CUDA_CHECK(cudaMalloc(&d_taxonomy_tree, num_taxonomy_nodes * sizeof(TaxonomyNode)));
    CUDA_CHECK(cudaMemcpy(d_taxonomy_tree, host_taxonomy.data(), 
                         num_taxonomy_nodes * sizeof(TaxonomyNode), cudaMemcpyHostToDevice));
    
    // Create parent lookup table for fast access
    std::vector<uint32_t> parent_lookup(max_taxon_id + 1, 0);
    for (const auto& [taxon_id, parent_id] : taxon_parents) {
        if (taxon_id <= max_taxon_id) {
            parent_lookup[taxon_id] = parent_id;
        }
    }
    
    CUDA_CHECK(cudaMalloc(&d_parent_lookup, (max_taxon_id + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_parent_lookup, parent_lookup.data(), 
                         (max_taxon_id + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    return true;
}

// Utility device functions (same as single-end version)
__device__ uint32_t lookup_lca_gpu(const GPUCompactHashTable* cht, uint64_t minimizer_hash) {
    uint32_t compact_hash = compute_compact_hash(minimizer_hash);
    uint32_t pos = compact_hash & cht->hash_mask;
    uint32_t lca_mask = (1U << cht->lca_bits) - 1;
    
    for (int probe = 0; probe < 32; probe++) {
        uint32_t cell = cht->hash_cells[pos];
        if (cell == 0) return 0;
        
        uint32_t stored_hash = cell >> cht->lca_bits;
        uint32_t expected_hash = compact_hash >> cht->lca_bits;
        
        if (stored_hash == expected_hash) {
            return cell & lca_mask;
        }
        
        pos = (pos + 1) & cht->hash_mask;
    }
    
    return 0;
}

__device__ uint64_t extract_minimizer_with_spaced_seeds(
    const char* sequence, int pos, int k, int ell, int spaces) {
    
    uint64_t min_hash = UINT64_MAX;
    
    for (int i = 0; i <= k - ell; i++) {
        uint64_t hash = hash_minimizer(sequence + pos + i, ell);
        if (hash == UINT64_MAX) continue;
        
        if (spaces > 0) {
            hash = apply_spaced_seed_mask(hash, spaces);
        }
        
        if (hash < min_hash) {
            min_hash = hash;
        }
    }
    
    return min_hash;
}

// Utility functions (same as before)
__device__ __host__ uint64_t hash_minimizer(const char* seq, int len) {
    uint64_t hash = 0;
    for (int i = 0; i < len; i++) {
        int base;
        switch (seq[i]) {
            case 'A': case 'a': base = 0; break;
            case 'C': case 'c': base = 1; break;
            case 'G': case 'g': base = 2; break;
            case 'T': case 't': base = 3; break;
            default: return UINT64_MAX;
        }
        hash = (hash << 2) | base;
    }
    return hash;
}

__device__ __host__ uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces) {
    uint64_t masked = 0;
    int out_pos = 0;
    
    for (int i = 0; i < 32; i++) {
        if (i % (spaces + 1) != 0) {
            masked |= ((hash >> (i * 2)) & 3ULL) << (out_pos * 2);
            out_pos++;
        }
    }
    
    return masked;
}

__device__ __host__ bool has_ambiguous_bases(const char* seq, int len) {
    for (int i = 0; i < len; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return true;
        }
    }
    return false;
}

__device__ __host__ uint32_t compute_compact_hash(uint64_t minimizer_hash) {
    return jenkins_hash(minimizer_hash) & 0x7FFFFFFF;
}

__device__ __host__ uint32_t jenkins_hash(uint64_t key) {
    uint32_t hash = (uint32_t)(key ^ (key >> 32));
    hash += (hash << 10);
    hash ^= (hash >> 6);
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

#endif // GPU_KRAKEN_CLASSIFIER_HEADER_ONLY
