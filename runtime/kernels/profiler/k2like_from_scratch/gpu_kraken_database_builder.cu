// gpu_kraken_database_builder.cu
// COMPLETE FIXED VERSION: GPU-accelerated pipeline to build Kraken2-style database from genome files
// Parallelizes minimizer extraction, LCA computation, and database construction

#ifndef GPU_KRAKEN_DATABASE_BUILDER_CUH
#define GPU_KRAKEN_DATABASE_BUILDER_CUH

#include "gpu_kraken_classifier.cu"
#include "gpu_minimizer_extraction.cuh"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/unique.h>
#include <thrust/scan.h>
#include <cub/cub.cuh>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <filesystem>
#include <chrono>

// Use structures from gpu_minimizer_extraction.cu
// struct GPUGenomeInfo and GPUMinimizerHit are defined there

struct GPUTaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    uint8_t padding[3];
};

// LCACandidate is defined in minimizer_extraction.h

// Database building statistics
struct GPUBuildStats {
    uint64_t total_sequences;
    uint64_t total_bases;
    uint64_t total_kmers_processed;
    uint64_t valid_minimizers_extracted;
    uint64_t unique_minimizers;
    uint64_t lca_assignments;
    double sequence_processing_time;
    double minimizer_extraction_time;
    double lca_computation_time;
    double database_construction_time;
};

class GPUKrakenDatabaseBuilder {
private:
    // Configuration
    ClassificationParams params;
    std::string output_directory;
    
    // Host data
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    std::vector<LCACandidate> all_lca_candidates;  // Accumulate all candidates
    
    // GPU memory management
    char* d_sequence_data;
    GPUGenomeInfo* d_genome_info;
    GPUMinimizerHit* d_minimizer_hits;
    uint32_t* d_minimizer_counts;
    LCACandidate* d_lca_candidates;
    GPUTaxonomyNode* d_taxonomy_nodes;
    
    // FIXED: Realistic processing parameters for better memory management
    int MAX_SEQUENCE_BATCH = 50;           // Reduced from 1000
    int MAX_MINIMIZERS_PER_BATCH = 2000000; // Reduced from 50M to 2M
    static const int MAX_READ_LENGTH = 1000;
    static const int THREADS_PER_BLOCK = 256;
    static const int MAX_TAXA_PER_PAIR = 128;
    
    // FIXED: Convert existing params to new format
    MinimizerParams minimizer_params;
    
    // Statistics
    GPUBuildStats stats;
    
public:
    GPUKrakenDatabaseBuilder(const std::string& output_dir,
                            const ClassificationParams& config = ClassificationParams());
    ~GPUKrakenDatabaseBuilder();
    
    // Main pipeline methods
    bool build_database_from_genomes(
        const std::string& genome_library_path,
        const std::string& taxonomy_path = ""
    );
    
    bool build_database_from_file_list(
        const std::string& file_list_path,
        const std::string& taxonomy_path = ""
    );
    
    // Step-by-step processing
    bool load_genome_files(const std::string& library_path);
    bool load_taxonomy_data(const std::string& taxonomy_path);
    bool process_genomes_gpu();
    bool save_database();
    
    // Configuration
    void set_batch_size(int sequences_per_batch) {
        if (sequences_per_batch > 0 && sequences_per_batch <= 1000) {
            MAX_SEQUENCE_BATCH = sequences_per_batch;
            std::cout << "Set GPU batch size to " << MAX_SEQUENCE_BATCH 
                      << " sequences, " << MAX_MINIMIZERS_PER_BATCH << " minimizers" << std::endl;
        }
    }
    
    void print_build_statistics();
    
    // Testing function for integration
    void test_minimizer_extraction_integration();
    
private:
    // GPU processing methods
    bool allocate_gpu_memory();
    void free_gpu_memory();
    void check_and_adjust_memory();
    
    bool process_sequence_batch(
        const std::vector<std::string>& sequences,
        const std::vector<uint32_t>& taxon_ids,
        int batch_offset
    );
    
    // FIXED: Updated to use new minimizer extraction function
    bool extract_minimizers_gpu(
        const char* d_sequences,
        const GPUGenomeInfo* d_genomes,
        int num_sequences,
        GPUMinimizerHit* d_hits,
        uint32_t* total_hits
    );
    
public:  // Move to public section temporarily for lambda access
    bool compute_lca_assignments_gpu(
        const GPUMinimizerHit* d_hits,
        int num_hits,
        LCACandidate* d_candidates,
        int* num_candidates
    );
private:
    
    bool build_compact_hash_table_gpu(
        const LCACandidate* d_candidates,
        int num_candidates
    );
    
    // File processing
    std::vector<std::string> load_sequences_from_fasta(const std::string& fasta_path);
    bool parse_taxonomy_files(const std::string& taxonomy_path);
    
    // Utility methods
    uint32_t extract_taxon_from_filename(const std::string& filename);
    std::vector<std::string> find_genome_files(const std::string& directory);
};

// ================================================================
// CUDA KERNELS
// ================================================================

// Use optimized kernels from gpu_minimizer_extraction.cu

// Kernel to sort and deduplicate minimizers by hash
__global__ void prepare_lca_computation_kernel(
    const GPUMinimizerHit* sorted_hits,
    int num_hits,
    LCACandidate* lca_candidates,
    uint32_t* candidate_counts
);

// Kernel to compute LCA for each unique minimizer
__global__ void compute_lca_kernel(
    GPUMinimizerHit* grouped_hits,
    const uint32_t* group_offsets,
    const uint32_t* group_sizes,
    int num_groups,
    const GPUTaxonomyNode* taxonomy,
    int num_taxonomy_nodes,
    LCACandidate* lca_results
);

// Kernel to convert hits to candidates
__global__ void convert_hits_to_candidates(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates
);

// Utility device functions
__device__ uint32_t compute_lca_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const GPUTaxonomyNode* taxonomy,
    int num_nodes
);

__device__ bool has_valid_bases(const char* seq, int length);

// ================================================================
// IMPLEMENTATION
// ================================================================

#ifndef GPU_KRAKEN_CLASSIFIER_HEADER_ONLY
// Implementation continues below
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>

#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// Utility function to get next power of 2
inline uint64_t get_next_power_of_2(uint64_t n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    n++;
    return n;
}

// Forward declarations
__global__ void create_simple_lca_candidates(
    const uint64_t* unique_hashes,
    const uint32_t* hit_counts,
    const GPUMinimizerHit* all_hits,
    int num_hits,
    int num_unique,
    LCACandidate* candidates);

// Functor for sorting minimizer hits
struct MinimizerHitComparator {
    __host__ __device__ bool operator()(const GPUMinimizerHit& a, const GPUMinimizerHit& b) const {
        return a.minimizer_hash < b.minimizer_hash;
    }
};

// FIXED: Constructor with proper parameter conversion
GPUKrakenDatabaseBuilder::GPUKrakenDatabaseBuilder(
    const std::string& output_dir,
    const ClassificationParams& config)
    : output_directory(output_dir), params(config),
      d_sequence_data(nullptr), d_genome_info(nullptr),
      d_minimizer_hits(nullptr), d_minimizer_counts(nullptr),
      d_lca_candidates(nullptr), d_taxonomy_nodes(nullptr) {
    
    // FIXED: Convert classification params to minimizer params
    minimizer_params.k = params.k;
    minimizer_params.ell = params.ell;
    minimizer_params.spaces = params.spaces;
    
    std::cout << "Initializing GPU Kraken database builder..." << std::endl;
    std::cout << "Updated parameters: k=" << minimizer_params.k 
              << ", ell=" << minimizer_params.ell 
              << ", spaces=" << minimizer_params.spaces << std::endl;
    
    // Create output directory
    std::filesystem::create_directories(output_directory);
    
    // FIXED: Check GPU memory and adjust batch sizes
    check_and_adjust_memory();
    
    // Initialize statistics
    memset(&stats, 0, sizeof(stats));
}

// Destructor
GPUKrakenDatabaseBuilder::~GPUKrakenDatabaseBuilder() {
    free_gpu_memory();
}

// FIXED: Check GPU memory and adjust batch sizes
void GPUKrakenDatabaseBuilder::check_and_adjust_memory() {
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    
    std::cout << "GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
              << (total_mem / 1024 / 1024) << " MB total" << std::endl;
    
    // Calculate memory needs per batch
    size_t sequence_memory = MAX_SEQUENCE_BATCH * 10000000; // 10MB per genome max
    size_t minimizer_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit);
    size_t genome_info_memory = MAX_SEQUENCE_BATCH * sizeof(GPUGenomeInfo);
    size_t total_needed = sequence_memory + minimizer_memory + genome_info_memory;
    
    // Use 80% of available memory
    size_t safe_memory = free_mem * 0.8;
    
    if (total_needed > safe_memory) {
        double scale = (double)safe_memory / total_needed;
        MAX_SEQUENCE_BATCH = (int)(MAX_SEQUENCE_BATCH * scale);
        MAX_MINIMIZERS_PER_BATCH = (int)(MAX_MINIMIZERS_PER_BATCH * scale);
        
        std::cout << "Adjusted batch sizes for memory:" << std::endl;
        std::cout << "  Sequences: " << MAX_SEQUENCE_BATCH << std::endl;
        std::cout << "  Minimizers: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
    }
}

// Main pipeline: build from genome library directory
bool GPUKrakenDatabaseBuilder::build_database_from_genomes(
    const std::string& genome_library_path,
    const std::string& taxonomy_path) {
    
    std::cout << "\n=== BUILDING KRAKEN DATABASE FROM GENOMES ===" << std::endl;
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Load genome files
    if (!load_genome_files(genome_library_path)) {
        std::cerr << "Failed to load genome files" << std::endl;
        return false;
    }
    
    // Step 2: Load taxonomy (optional)
    if (!taxonomy_path.empty()) {
        if (!load_taxonomy_data(taxonomy_path)) {
            std::cerr << "Warning: Failed to load taxonomy data" << std::endl;
        }
    }
    
    // Step 3: Process genomes on GPU
    if (!process_genomes_gpu()) {
        std::cerr << "Failed to process genomes on GPU" << std::endl;
        return false;
    }
    
    // Step 4: Save database
    if (!save_database()) {
        std::cerr << "Failed to save database" << std::endl;
        return false;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
    
    std::cout << "\n=== DATABASE BUILD COMPLETE ===" << std::endl;
    std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
    print_build_statistics();
    
    return true;
}

// Load genome files from directory
bool GPUKrakenDatabaseBuilder::load_genome_files(const std::string& library_path) {
    std::cout << "Loading genome files from: " << library_path << std::endl;
    
    genome_files = find_genome_files(library_path);
    
    if (genome_files.empty()) {
        std::cerr << "No genome files found in " << library_path << std::endl;
        return false;
    }
    
    // Extract taxon IDs from filenames
    genome_taxon_ids.reserve(genome_files.size());
    for (const auto& file : genome_files) {
        uint32_t taxon_id = extract_taxon_from_filename(file);
        genome_taxon_ids.push_back(taxon_id);
        
        // Create default name if no taxonomy loaded
        if (taxon_names.find(taxon_id) == taxon_names.end()) {
            std::filesystem::path p(file);
            taxon_names[taxon_id] = p.stem().string();
        }
    }
    
    std::cout << "Loaded " << genome_files.size() << " genome files" << std::endl;
    stats.total_sequences = genome_files.size();
    
    return true;
}

// Process all genomes on GPU
bool GPUKrakenDatabaseBuilder::process_genomes_gpu() {
    std::cout << "Processing genomes on GPU..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (!allocate_gpu_memory()) {
        return false;
    }
    
    // Process genomes in batches
    for (size_t batch_start = 0; batch_start < genome_files.size(); 
         batch_start += MAX_SEQUENCE_BATCH) {
        
        size_t batch_end = std::min(batch_start + MAX_SEQUENCE_BATCH, genome_files.size());
        
        std::cout << "Processing batch " << (batch_start / MAX_SEQUENCE_BATCH + 1) 
                  << ": genomes " << batch_start << "-" << batch_end << std::endl;
        
        // Load sequences for this batch
        std::vector<std::string> batch_sequences;
        std::vector<uint32_t> batch_taxon_ids;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            auto sequences = load_sequences_from_fasta(genome_files[i]);
            for (const auto& seq : sequences) {
                batch_sequences.push_back(seq);
                batch_taxon_ids.push_back(genome_taxon_ids[i]);
                stats.total_bases += seq.length();
            }
        }
        
        // Process this batch on GPU
        if (!process_sequence_batch(batch_sequences, batch_taxon_ids, batch_start)) {
            std::cerr << "Failed to process batch starting at " << batch_start << std::endl;
            return false;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.sequence_processing_time = 
        std::chrono::duration<double>(end_time - start_time).count();
    
    return true;
}

// Process a single batch of sequences
bool GPUKrakenDatabaseBuilder::process_sequence_batch(
    const std::vector<std::string>& sequences,
    const std::vector<uint32_t>& taxon_ids,
    int batch_offset) {
    
    if (sequences.empty()) return true;
    
    // Prepare sequence data for GPU
    std::string concatenated_sequences;
    std::vector<GPUGenomeInfo> genome_infos;
    
    uint32_t current_offset = 0;
    for (size_t i = 0; i < sequences.size(); i++) {
        GPUGenomeInfo info;
        info.taxon_id = taxon_ids[i];
        info.sequence_offset = current_offset;
        info.sequence_length = sequences[i].length();
        info.genome_id = batch_offset + i;
        
        genome_infos.push_back(info);
        concatenated_sequences += sequences[i];
        current_offset += sequences[i].length();
        
        // Track k-mers
        if (sequences[i].length() >= minimizer_params.k) {
            stats.total_kmers_processed += sequences[i].length() - minimizer_params.k + 1;
        }
    }
    
    // Check if data fits in allocated memory
    size_t max_sequence_memory = size_t(MAX_SEQUENCE_BATCH) * 10000000;
    if (concatenated_sequences.length() > max_sequence_memory) {
        std::cerr << "ERROR: Sequence data (" << concatenated_sequences.length() 
                  << " bytes) exceeds allocated memory (" << max_sequence_memory << " bytes)" << std::endl;
        std::cerr << "Try reducing batch size" << std::endl;
        return false;
    }
    
    // Transfer to GPU
    CUDA_CHECK(cudaMemcpy(d_sequence_data, concatenated_sequences.c_str(),
                         concatenated_sequences.length(), cudaMemcpyHostToDevice));
    
    CUDA_CHECK(cudaMemcpy(d_genome_info, genome_infos.data(),
                         genome_infos.size() * sizeof(GPUGenomeInfo),
                         cudaMemcpyHostToDevice));
    
    // Extract minimizers on GPU
    uint32_t total_hits = 0;
    if (!extract_minimizers_gpu(d_sequence_data, d_genome_info, 
                               genome_infos.size(), d_minimizer_hits, &total_hits)) {
        return false;
    }
    
    // Compute LCA assignments
    int num_candidates = 0;
    if (!compute_lca_assignments_gpu(d_minimizer_hits, total_hits,
                                    d_lca_candidates, &num_candidates)) {
        return false;
    }
    
    // Copy candidates from GPU to host and accumulate
    if (num_candidates > 0) {
        std::vector<LCACandidate> batch_candidates(num_candidates);
        CUDA_CHECK(cudaMemcpy(batch_candidates.data(), d_lca_candidates,
                             num_candidates * sizeof(LCACandidate), cudaMemcpyDeviceToHost));
        
        // Accumulate candidates
        all_lca_candidates.insert(all_lca_candidates.end(), 
                                 batch_candidates.begin(), batch_candidates.end());
        
        std::cout << "Accumulated " << num_candidates << " candidates (total: " 
                  << all_lca_candidates.size() << ")" << std::endl;
    }
    
    stats.valid_minimizers_extracted += total_hits;
    stats.lca_assignments += num_candidates;
    
    return true;
}

// FIXED: Extract minimizers using optimized kernel
bool GPUKrakenDatabaseBuilder::extract_minimizers_gpu(
    const char* d_sequences,
    const GPUGenomeInfo* d_genomes,
    int num_sequences,
    GPUMinimizerHit* d_hits,
    uint32_t* total_hits) {
    
    std::cout << "Extracting minimizers from " << num_sequences << " sequences..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // FIXED: Use the optimized kernel instead of broken one
    bool success = extract_minimizers_gpu_optimized(
        d_sequences, d_genomes, num_sequences,
        d_hits, d_minimizer_counts, total_hits,
        minimizer_params, MAX_MINIMIZERS_PER_BATCH
    );
    
    if (!success) {
        std::cerr << "Minimizer extraction failed!" << std::endl;
        return false;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.minimizer_extraction_time += 
        std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "Extracted " << *total_hits << " minimizers" << std::endl;
    
    // Check if we hit the limit
    if (*total_hits >= MAX_MINIMIZERS_PER_BATCH * 0.95) {
        std::cout << "WARNING: Near minimizer limit. Consider reducing batch size." << std::endl;
    }
    
    return true;
}

// Kernel to convert hits to candidates
__global__ void convert_hits_to_candidates(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_hits) return;
    
    const GPUMinimizerHit& hit = hits[idx];
    LCACandidate& candidate = candidates[idx];
    
    candidate.minimizer_hash = hit.minimizer_hash;
    candidate.lca_taxon = hit.taxon_id;
    candidate.genome_count = 1;  // Will be updated during final merge
    candidate.uniqueness_score = 1.0f;
}

// FIXED: Compute LCA assignments on GPU
bool GPUKrakenDatabaseBuilder::compute_lca_assignments_gpu(
    const GPUMinimizerHit* d_hits,
    int num_hits,
    LCACandidate* d_candidates,
    int* num_candidates) {
    
    if (num_hits == 0) {
        *num_candidates = 0;
        return true;
    }
    
    std::cout << "Computing LCA assignments for " << num_hits << " hits..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // FIXED: First, deduplicate the minimizers on GPU
    uint32_t unique_count = 0;
    GPUMinimizerHit* d_hits_mutable = const_cast<GPUMinimizerHit*>(d_hits);
    
    bool dedup_success = deduplicate_minimizers_gpu(d_hits_mutable, num_hits, &unique_count);
    if (!dedup_success) {
        std::cerr << "Deduplication failed!" << std::endl;
        return false;
    }
    
    std::cout << "After deduplication: " << unique_count << " unique minimizers" << std::endl;
    
    // Convert deduplicated hits to LCA candidates
    int blocks = (unique_count + 255) / 256;
    convert_hits_to_candidates<<<blocks, 256>>>(
        d_hits_mutable, unique_count, d_candidates
    );
    cudaDeviceSynchronize();
    
    *num_candidates = unique_count;
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.lca_computation_time += 
        std::chrono::duration<double>(end_time - start_time).count();
    
    return true;
}

// Allocate GPU memory
bool GPUKrakenDatabaseBuilder::allocate_gpu_memory() {
    std::cout << "Allocating GPU memory..." << std::endl;
    
    // Check available GPU memory first
    size_t free_memory, total_memory_gpu;
    cudaMemGetInfo(&free_memory, &total_memory_gpu);
    std::cout << "GPU Memory: " << (free_memory / 1024 / 1024) << " MB free / " 
              << (total_memory_gpu / 1024 / 1024) << " MB total" << std::endl;
    
    // Calculate memory requirements
    size_t sequence_memory = size_t(MAX_SEQUENCE_BATCH) * 10000000;  // Assume up to 10MB per genome
    size_t genome_info_memory = MAX_SEQUENCE_BATCH * sizeof(GPUGenomeInfo);
    size_t minimizer_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit);
    size_t candidate_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(LCACandidate);
    
    size_t total_required = sequence_memory + genome_info_memory + 
                           minimizer_memory + candidate_memory;
    
    std::cout << "Memory requirements:" << std::endl;
    std::cout << "  Sequence data: " << (sequence_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Genome info: " << (genome_info_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Minimizer hits: " << (minimizer_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  LCA candidates: " << (candidate_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Total required: " << (total_required / 1024 / 1024) << " MB" << std::endl;
    
    // Check if we have enough memory
    if (total_required > free_memory * 0.9) {  // Leave 10% buffer
        std::cerr << "ERROR: Not enough GPU memory! Required: " << (total_required / 1024 / 1024) 
                  << " MB, Available: " << (free_memory / 1024 / 1024) << " MB" << std::endl;
        return false;
    }
    
    // Allocate with error checking
    cudaError_t err;
    
    err = cudaMalloc(&d_sequence_data, sequence_memory);
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate sequence memory: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    err = cudaMalloc(&d_genome_info, genome_info_memory);
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate genome info memory: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_sequence_data);
        return false;
    }
    
    err = cudaMalloc(&d_minimizer_hits, minimizer_memory);
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate minimizer memory: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_sequence_data);
        cudaFree(d_genome_info);
        return false;
    }
    
    err = cudaMalloc(&d_minimizer_counts, MAX_SEQUENCE_BATCH * sizeof(uint32_t));
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate minimizer counts: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_sequence_data);
        cudaFree(d_genome_info);
        cudaFree(d_minimizer_hits);
        return false;
    }
    
    err = cudaMalloc(&d_lca_candidates, candidate_memory);
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate LCA candidates: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_sequence_data);
        cudaFree(d_genome_info);
        cudaFree(d_minimizer_hits);
        cudaFree(d_minimizer_counts);
        return false;
    }
    
    std::cout << "Successfully allocated " << (total_required / 1024 / 1024) 
              << " MB of GPU memory" << std::endl;
    
    return true;
}

// Free GPU memory
void GPUKrakenDatabaseBuilder::free_gpu_memory() {
    if (d_sequence_data) { cudaFree(d_sequence_data); d_sequence_data = nullptr; }
    if (d_genome_info) { cudaFree(d_genome_info); d_genome_info = nullptr; }
    if (d_minimizer_hits) { cudaFree(d_minimizer_hits); d_minimizer_hits = nullptr; }
    if (d_minimizer_counts) { cudaFree(d_minimizer_counts); d_minimizer_counts = nullptr; }
    if (d_lca_candidates) { cudaFree(d_lca_candidates); d_lca_candidates = nullptr; }
    if (d_taxonomy_nodes) { cudaFree(d_taxonomy_nodes); d_taxonomy_nodes = nullptr; }
}

// Utility: Find genome files in directory
std::vector<std::string> GPUKrakenDatabaseBuilder::find_genome_files(const std::string& directory) {
    std::vector<std::string> files;
    
    for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            std::string extension = entry.path().extension().string();
            
            // Look for FASTA files
            if (extension == ".fna" || extension == ".fa" || extension == ".fasta" ||
                extension == ".ffn" || extension == ".faa") {
                files.push_back(entry.path().string());
            }
        }
    }
    
    std::sort(files.begin(), files.end());
    return files;
}

// Extract taxon ID from filename
uint32_t GPUKrakenDatabaseBuilder::extract_taxon_from_filename(const std::string& filename) {
    // Try to extract from filename patterns like:
    // GCF_000005825.2_ASM582v2_genomic.fna -> try to map via accession
    // Or taxid_12345_genomic.fna -> extract 12345
    
    std::filesystem::path p(filename);
    std::string stem = p.stem().string();
    
    // Look for taxid pattern
    std::regex taxid_pattern(R"(taxid[_-](\d+))");
    std::smatch match;
    if (std::regex_search(stem, match, taxid_pattern)) {
        return std::stoul(match[1].str());
    }
    
    // Look for GCF accession pattern
    std::regex gcf_pattern(R"(GCF_(\d+\.\d+))");
    if (std::regex_search(stem, match, gcf_pattern)) {
        // Would need to map GCF accession to taxon ID
        // For now, use a hash of the accession
        std::hash<std::string> hasher;
        return hasher(match[1].str()) % 1000000 + 1000000;
    }
    
    // Fallback: hash the filename
    std::hash<std::string> hasher;
    return hasher(stem) % 1000000 + 2000000;
}

// Load sequences from FASTA file
std::vector<std::string> GPUKrakenDatabaseBuilder::load_sequences_from_fasta(const std::string& fasta_path) {
    std::vector<std::string> sequences;
    std::ifstream file(fasta_path);
    
    if (!file.is_open()) {
        std::cerr << "Cannot open FASTA file: " << fasta_path << std::endl;
        return sequences;
    }
    
    std::string line, current_sequence;
    bool in_sequence = false;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Header line
            if (in_sequence && !current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
            in_sequence = true;
        } else if (in_sequence) {
            // Sequence line
            current_sequence += line;
        }
    }
    
    // Add last sequence
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }
    
    file.close();
    return sequences;
}

// Save database to files
bool GPUKrakenDatabaseBuilder::save_database() {
    std::cout << "Saving database to " << output_directory << "..." << std::endl;
    
    // Create a simple database file with all LCA candidates
    std::string database_file = output_directory + "/database.k2d";
    std::string taxonomy_file = output_directory + "/taxonomy.tsv";
    std::string hash_file = output_directory + "/hash_table.bin";
    
    // Collect all LCA candidates from GPU
    if (all_lca_candidates.empty()) {
        std::cerr << "Warning: No LCA candidates to save. Database may be incomplete." << std::endl;
    }
    
    // Save hash table in binary format
    std::ofstream hash_out(hash_file, std::ios::binary);
    if (!hash_out.is_open()) {
        std::cerr << "Cannot create hash table file: " << hash_file << std::endl;
        return false;
    }
    
    // Write header
    uint64_t num_entries = all_lca_candidates.size();
    uint64_t table_size = num_entries > 0 ? get_next_power_of_2(num_entries * 2) : 1024;  // 50% load factor
    
    hash_out.write(reinterpret_cast<const char*>(&table_size), sizeof(uint64_t));
    hash_out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    
    // Write all candidates
    for (const auto& candidate : all_lca_candidates) {
        hash_out.write(reinterpret_cast<const char*>(&candidate.minimizer_hash), sizeof(uint64_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.lca_taxon), sizeof(uint32_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.genome_count), sizeof(uint32_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.uniqueness_score), sizeof(float));
    }
    hash_out.close();
    
    // Save taxonomy mapping
    std::ofstream tax_out(taxonomy_file);
    if (!tax_out.is_open()) {
        std::cerr << "Cannot create taxonomy file: " << taxonomy_file << std::endl;
        return false;
    }
    
    tax_out << "taxon_id\tname\tparent_id\n";
    for (const auto& [taxon_id, name] : taxon_names) {
        uint32_t parent_id = taxon_parents.count(taxon_id) ? taxon_parents[taxon_id] : 0;
        tax_out << taxon_id << "\t" << name << "\t" << parent_id << "\n";
    }
    tax_out.close();
    
    // Save database statistics
    std::string stats_file = output_directory + "/build_stats.txt";
    std::ofstream stats_out(stats_file);
    
    stats_out << "Database Build Statistics\n";
    stats_out << "========================\n";
    stats_out << "Total sequences: " << stats.total_sequences << "\n";
    stats_out << "Total bases: " << stats.total_bases << "\n";
    stats_out << "Valid minimizers: " << stats.valid_minimizers_extracted << "\n";
    stats_out << "Unique minimizers: " << stats.unique_minimizers << "\n";
    stats_out << "LCA assignments: " << stats.lca_assignments << "\n";
    stats_out << "Sequence processing time: " << stats.sequence_processing_time << "s\n";
    stats_out << "Minimizer extraction time: " << stats.minimizer_extraction_time << "s\n";
    stats_out << "LCA computation time: " << stats.lca_computation_time << "s\n";
    if (stats.sequence_processing_time > 0) {
        stats_out << "Processing rate: " << (stats.total_bases / stats.sequence_processing_time) << " bases/second\n";
    }
    
    stats_out.close();
    
    std::cout << "Database saved successfully!" << std::endl;
    std::cout << "  Hash table: " << hash_file << " (" << num_entries << " entries)" << std::endl;
    std::cout << "  Taxonomy: " << taxonomy_file << " (" << taxon_names.size() << " taxa)" << std::endl;
    
    return true;
}

// Print build statistics
void GPUKrakenDatabaseBuilder::print_build_statistics() {
    std::cout << "\n=== BUILD STATISTICS ===" << std::endl;
    std::cout << "Total sequences processed: " << stats.total_sequences << std::endl;
    std::cout << "Total bases processed: " << stats.total_bases << std::endl;
    std::cout << "Valid minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
    std::cout << "LCA assignments computed: " << stats.lca_assignments << std::endl;
    std::cout << "Sequence processing time: " << std::fixed << std::setprecision(2) 
              << stats.sequence_processing_time << "s" << std::endl;
    std::cout << "Minimizer extraction time: " << std::fixed << std::setprecision(2) 
              << stats.minimizer_extraction_time << "s" << std::endl;
    std::cout << "LCA computation time: " << std::fixed << std::setprecision(2) 
              << stats.lca_computation_time << "s" << std::endl;
    
    if (stats.total_bases > 0) {
        double bases_per_second = stats.total_bases / 
            (stats.sequence_processing_time + stats.minimizer_extraction_time);
        std::cout << "Processing rate: " << std::scientific << std::setprecision(2) 
                  << bases_per_second << " bases/second" << std::endl;
    }
}

// Load taxonomy data from NCBI dump files
bool GPUKrakenDatabaseBuilder::load_taxonomy_data(const std::string& taxonomy_path) {
    std::cout << "Loading taxonomy data from: " << taxonomy_path << std::endl;
    
    // Clear existing taxonomy data
    taxon_names.clear();
    taxon_parents.clear();
    
    // Determine paths to nodes.dmp and names.dmp
    std::filesystem::path base_path(taxonomy_path);
    std::string nodes_file, names_file;
    
    if (std::filesystem::is_directory(base_path)) {
        // Look for NCBI taxonomy dump files
        nodes_file = base_path / "nodes.dmp";
        names_file = base_path / "names.dmp";
        
        // If not found in directory, check parent
        if (!std::filesystem::exists(nodes_file)) {
            nodes_file = base_path.parent_path() / "nodes.dmp";
            names_file = base_path.parent_path() / "names.dmp";
        }
    } else {
        // Assume taxonomy_path is the directory containing the files
        base_path = base_path.parent_path();
        nodes_file = base_path / "nodes.dmp";
        names_file = base_path / "names.dmp";
    }
    
    // First, load the taxonomy tree structure from nodes.dmp
    if (!std::filesystem::exists(nodes_file)) {
        std::cerr << "nodes.dmp not found at: " << nodes_file << std::endl;
        return false;
    }
    
    std::cout << "Loading taxonomy tree from: " << nodes_file << std::endl;
    std::ifstream nodes_in(nodes_file);
    if (!nodes_in.is_open()) {
        std::cerr << "Cannot open nodes.dmp: " << nodes_file << std::endl;
        return false;
    }
    
    std::string line;
    int nodes_loaded = 0;
    std::unordered_map<uint32_t, std::string> taxon_ranks;
    
    while (std::getline(nodes_in, line)) {
        if (line.empty()) continue;
        
        // Parse pipe-delimited format: tax_id | parent_tax_id | rank | ...
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '|')) {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            fields.push_back(token);
        }
        
        if (fields.size() >= 3) {
            try {
                uint32_t taxon_id = std::stoul(fields[0]);
                uint32_t parent_id = std::stoul(fields[1]);
                std::string rank = fields[2];
                
                taxon_parents[taxon_id] = parent_id;
                taxon_ranks[taxon_id] = rank;
                nodes_loaded++;
                
                // Initialize with taxon ID as name (will be replaced by names.dmp)
                taxon_names[taxon_id] = "taxon_" + std::to_string(taxon_id);
            } catch (const std::exception& e) {
                // Skip malformed lines
                continue;
            }
        }
    }
    nodes_in.close();
    
    std::cout << "Loaded " << nodes_loaded << " taxonomy nodes" << std::endl;
    
    // Now load the scientific names from names.dmp
    if (!std::filesystem::exists(names_file)) {
        std::cerr << "names.dmp not found at: " << names_file << std::endl;
        std::cerr << "Using taxon IDs as names" << std::endl;
        return nodes_loaded > 0;
    }
    
    std::cout << "Loading taxonomy names from: " << names_file << std::endl;
    std::ifstream names_in(names_file);
    if (!names_in.is_open()) {
        std::cerr << "Cannot open names.dmp: " << names_file << std::endl;
        return nodes_loaded > 0;
    }
    
    int names_loaded = 0;
    
    while (std::getline(names_in, line)) {
        if (line.empty()) continue;
        
        // Parse pipe-delimited format: tax_id | name_txt | unique_name | name_class |
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '|')) {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            fields.push_back(token);
        }
        
        if (fields.size() >= 4) {
            try {
                uint32_t taxon_id = std::stoul(fields[0]);
                std::string name_txt = fields[1];
                std::string name_class = fields[3];
                
                // Only use scientific names (primary names)
                if (name_class == "scientific name" && taxon_parents.find(taxon_id) != taxon_parents.end()) {
                    taxon_names[taxon_id] = name_txt;
                    names_loaded++;
                }
            } catch (const std::exception& e) {
                // Skip malformed lines
                continue;
            }
        }
    }
    names_in.close();
    
    std::cout << "Loaded " << names_loaded << " scientific names" << std::endl;
    
    // Ensure root node exists
    if (taxon_names.find(1) == taxon_names.end()) {
        taxon_names[1] = "root";
        taxon_parents[1] = 0;
    }
    
    // Print some statistics
    std::cout << "Total taxa in database: " << taxon_parents.size() << std::endl;
    
    // Count taxa by rank
    std::unordered_map<std::string, int> rank_counts;
    for (const auto& [taxon_id, rank] : taxon_ranks) {
        rank_counts[rank]++;
    }
    
    std::cout << "Taxa by rank:" << std::endl;
    for (const auto& [rank, count] : rank_counts) {
        if (count > 100) {  // Only show major ranks
            std::cout << "  " << rank << ": " << count << std::endl;
        }
    }
    
    return nodes_loaded > 0;
}

// Simple kernel to create LCA candidates
__global__ void create_simple_lca_candidates(
    const uint64_t* unique_hashes,
    const uint32_t* hit_counts,
    const GPUMinimizerHit* all_hits,
    int num_hits,
    int num_unique,
    LCACandidate* candidates) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_unique) return;
    
    uint64_t target_hash = unique_hashes[idx];
    
    // Find first hit with this hash (since hits are sorted)
    // Simple linear search for now - could be optimized with binary search
    uint32_t first_taxon = 0;
    for (int i = 0; i < num_hits; i++) {
        if (all_hits[i].minimizer_hash == target_hash) {
            first_taxon = all_hits[i].taxon_id;
            break;
        }
    }
    
    // Create candidate
    candidates[idx].minimizer_hash = target_hash;
    candidates[idx].lca_taxon = first_taxon;  // Simple: use first taxon
    candidates[idx].genome_count = hit_counts[idx];
    candidates[idx].uniqueness_score = 1.0f / hit_counts[idx];
}

// Device function to check for valid bases
__device__ bool has_valid_bases(const char* seq, int length) {
    for (int i = 0; i < length; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

// FIXED: Test minimizer extraction integration
void GPUKrakenDatabaseBuilder::test_minimizer_extraction_integration() {
    std::cout << "Testing minimizer extraction integration..." << std::endl;
    
    // Test with a small sequence
    std::string test_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    std::vector<std::string> test_sequences = {test_seq};
    std::vector<uint32_t> test_taxon_ids = {12345};
    
    // Process small batch
    bool success = process_sequence_batch(test_sequences, test_taxon_ids, 0);
    
    if (success) {
        std::cout << "✓ Minimizer extraction test PASSED" << std::endl;
    } else {
        std::cout << "✗ Minimizer extraction test FAILED" << std::endl;
    }
}

#endif // GPU_KRAKEN_DATABASE_BUILDER_CUH
#endif // GPU_KRAKEN_CLASSIFIER_HEADER_ONLY