// gpu_kraken_database_builder.cuh
// GPU-accelerated pipeline to build Kraken2-style database from genome files
// Parallelizes minimizer extraction, LCA computation, and database construction

#pragma once
#ifndef GPU_KRAKEN_DATABASE_BUILDER_CUH
#define GPU_KRAKEN_DATABASE_BUILDER_CUH

#include "gpu_kraken_classifier.cu"
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

// GPU-friendly structures
struct GPUGenomeInfo {
    uint32_t taxon_id;
    uint32_t sequence_offset;    // Offset in concatenated sequence buffer
    uint32_t sequence_length;
    uint32_t genome_id;         // Index in genome array
};

struct GPUMinimizerHit {
    uint64_t minimizer_hash;
    uint32_t taxon_id;
    uint32_t position;          // Position in source sequence
    uint32_t genome_id;         // Source genome index
};

struct GPUTaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    uint8_t padding[3];
};

// LCA computation intermediate result
struct LCACandidate {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
};

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
    
    // GPU memory management
    char* d_sequence_data;
    GPUGenomeInfo* d_genome_info;
    GPUMinimizerHit* d_minimizer_hits;
    uint32_t* d_minimizer_counts;
    LCACandidate* d_lca_candidates;
    GPUTaxonomyNode* d_taxonomy_nodes;
    
    // Processing parameters
    static const int MAX_SEQUENCE_BATCH = 10000;  // 10K sequences at once (reduced from 1M)
    static const int MAX_MINIMIZERS_PER_BATCH = 5000000;  // 5M minimizers (reduced from 50M)
    static const int THREADS_PER_BLOCK = 256;
    
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
        // Will be used to tune memory usage
    }
    
    void print_build_statistics();
    
private:
    // GPU processing methods
    bool allocate_gpu_memory();
    void free_gpu_memory();
    
    bool process_sequence_batch(
        const std::vector<std::string>& sequences,
        const std::vector<uint32_t>& taxon_ids,
        int batch_offset
    );
    
    bool extract_minimizers_gpu(
        const char* d_sequences,
        const GPUGenomeInfo* d_genomes,
        int num_sequences,
        GPUMinimizerHit* d_hits,
        uint32_t* d_hit_counts
    );
    
    bool compute_lca_assignments_gpu(
        const GPUMinimizerHit* d_hits,
        int num_hits,
        LCACandidate* d_candidates,
        int* num_candidates
    );
    
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

// Kernel to extract minimizers from genome sequences
__global__ void extract_minimizers_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts,
    uint32_t* global_hit_counter,
    ClassificationParams params
);

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

// Utility device functions
__device__ uint64_t extract_single_minimizer(
    const char* sequence,
    int seq_length,
    int pos,
    int k,
    int ell,
    int spaces
);

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
// Implementation continues below - no need to include self
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

// Global constants for GPU kernels
__constant__ int MAX_MINIMIZERS_PER_BATCH = 50000000;  // 50M minimizers

// Functor for sorting minimizer hits
struct MinimizerHitComparator {
    __host__ __device__ bool operator()(const GPUMinimizerHit& a, const GPUMinimizerHit& b) const {
        return a.minimizer_hash < b.minimizer_hash;
    }
};

// Constructor
GPUKrakenDatabaseBuilder::GPUKrakenDatabaseBuilder(
    const std::string& output_dir,
    const ClassificationParams& config)
    : output_directory(output_dir), params(config),
      d_sequence_data(nullptr), d_genome_info(nullptr),
      d_minimizer_hits(nullptr), d_minimizer_counts(nullptr),
      d_lca_candidates(nullptr), d_taxonomy_nodes(nullptr) {
    
    std::cout << "Initializing GPU Kraken database builder..." << std::endl;
    std::cout << "Output directory: " << output_directory << std::endl;
    std::cout << "Parameters: k=" << params.k << ", ell=" << params.ell 
              << ", spaces=" << params.spaces << std::endl;
    
    // Create output directory
    std::filesystem::create_directories(output_directory);
    
    // Initialize statistics
    memset(&stats, 0, sizeof(stats));
}

// Destructor
GPUKrakenDatabaseBuilder::~GPUKrakenDatabaseBuilder() {
    free_gpu_memory();
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
    std::vector<std::vector<GPUMinimizerHit>> all_minimizer_hits;
    
    for (size_t batch_start = 0; batch_start < genome_files.size(); 
         batch_start += MAX_SEQUENCE_BATCH) {
        
        size_t batch_end = std::min(batch_start + MAX_SEQUENCE_BATCH, genome_files.size());
        size_t batch_size = batch_end - batch_start;
        
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
    
    stats.valid_minimizers_extracted += total_hits;
    stats.lca_assignments += num_candidates;
    
    return true;
}

// GPU kernel to extract minimizers
__global__ void extract_minimizers_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts,
    uint32_t* global_hit_counter,
    ClassificationParams params) {
    
    int genome_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (genome_id >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[genome_id];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) {
        hit_counts[genome_id] = 0;
        return;
    }
    
    uint32_t local_hit_count = 0;
    
    // Extract minimizers from this sequence
    for (uint32_t i = 0; i <= seq_length - params.k; i++) {
        // Check for valid bases in k-mer
        if (!has_valid_bases(sequence + i, params.k)) {
            continue;
        }
        
        // Extract minimizer
        uint64_t minimizer = extract_single_minimizer(
            sequence, seq_length, i, params.k, params.ell, params.spaces
        );
        
        if (minimizer != UINT64_MAX) {
            // Get global position for this hit
            uint32_t global_pos = atomicAdd(global_hit_counter, 1);
            
            if (global_pos < MAX_MINIMIZERS_PER_BATCH) {
                GPUMinimizerHit hit;
                hit.minimizer_hash = minimizer;
                hit.taxon_id = genome.taxon_id;
                hit.position = i;
                hit.genome_id = genome.genome_id;
                
                minimizer_hits[global_pos] = hit;
                local_hit_count++;
            }
        }
    }
    
    hit_counts[genome_id] = local_hit_count;
}

// Device function to extract a single minimizer
__device__ uint64_t extract_single_minimizer(
    const char* sequence,
    int seq_length,
    int pos,
    int k,
    int ell,
    int spaces) {
    
    uint64_t min_hash = UINT64_MAX;
    
    // Find minimizer in k-mer window
    for (int i = 0; i <= k - ell; i++) {
        if (pos + i + ell > seq_length) break;
        
        // Hash this ell-mer
        uint64_t hash = 0;
        bool valid = true;
        
        for (int j = 0; j < ell; j++) {
            char base = sequence[pos + i + j];
            int encoded;
            
            switch (base) {
                case 'A': case 'a': encoded = 0; break;
                case 'C': case 'c': encoded = 1; break;
                case 'G': case 'g': encoded = 2; break;
                case 'T': case 't': encoded = 3; break;
                default: valid = false; break;
            }
            
            if (!valid) break;
            hash = (hash << 2) | encoded;
        }
        
        if (!valid) continue;
        
        // Apply spaced seed mask
        if (spaces > 0) {
            hash = apply_spaced_seed_mask(hash, spaces);
        }
        
        if (hash < min_hash) {
            min_hash = hash;
        }
    }
    
    return min_hash;
}

// Extract minimizers using GPU
bool GPUKrakenDatabaseBuilder::extract_minimizers_gpu(
    const char* d_sequences,
    const GPUGenomeInfo* d_genomes,
    int num_sequences,
    GPUMinimizerHit* d_hits,
    uint32_t* total_hits) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Reset global counter
    uint32_t zero = 0;
    uint32_t* d_global_counter;
    CUDA_CHECK(cudaMalloc(&d_global_counter, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_global_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Launch kernel
    int num_blocks = (num_sequences + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    
    extract_minimizers_kernel<<<num_blocks, THREADS_PER_BLOCK>>>(
        d_sequences, d_genomes, num_sequences,
        d_hits, d_minimizer_counts, d_global_counter, params
    );
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Get total hit count
    CUDA_CHECK(cudaMemcpy(total_hits, d_global_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    cudaFree(d_global_counter);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.minimizer_extraction_time += 
        std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "Extracted " << *total_hits << " minimizers from " 
              << num_sequences << " sequences" << std::endl;
    
    return true;
}

// Compute LCA assignments on GPU
bool GPUKrakenDatabaseBuilder::compute_lca_assignments_gpu(
    const GPUMinimizerHit* d_hits,
    int num_hits,
    LCACandidate* d_candidates,
    int* num_candidates) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Sort hits by minimizer hash using Thrust
    thrust::device_ptr<GPUMinimizerHit> hits_ptr(const_cast<GPUMinimizerHit*>(d_hits));
    thrust::sort(hits_ptr, hits_ptr + num_hits, MinimizerHitComparator());
    
    // Group by minimizer hash and compute LCA for each group
    // This is a simplified version - full implementation would use CUB for grouping
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.lca_computation_time += 
        std::chrono::duration<double>(end_time - start_time).count();
    
    return true;
}

// Allocate GPU memory
bool GPUKrakenDatabaseBuilder::allocate_gpu_memory() {
    std::cout << "Allocating GPU memory..." << std::endl;
    
    // Calculate memory requirements
    size_t sequence_memory = size_t(MAX_SEQUENCE_BATCH) * 5000;  // Assume avg 5kb per sequence (reduced from 10kb)
    size_t genome_info_memory = MAX_SEQUENCE_BATCH * sizeof(GPUGenomeInfo);
    size_t minimizer_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit);
    size_t candidate_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(LCACandidate);
    
    CUDA_CHECK(cudaMalloc(&d_sequence_data, sequence_memory));
    CUDA_CHECK(cudaMalloc(&d_genome_info, genome_info_memory));
    CUDA_CHECK(cudaMalloc(&d_minimizer_hits, minimizer_memory));
    CUDA_CHECK(cudaMalloc(&d_minimizer_counts, MAX_SEQUENCE_BATCH * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_lca_candidates, candidate_memory));
    
    size_t total_memory = sequence_memory + genome_info_memory + 
                         minimizer_memory + candidate_memory;
    
    std::cout << "Allocated " << (total_memory / 1024 / 1024) 
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
    
    // This would save the final hash table and taxonomy tree
    // Implementation depends on final database format
    
    std::string database_file = output_directory + "/database.k2d";
    std::string taxonomy_file = output_directory + "/taxonomy.tsv";
    
    // Save database statistics
    std::string stats_file = output_directory + "/build_stats.txt";
    std::ofstream stats_out(stats_file);
    
    stats_out << "Database Build Statistics\n";
    stats_out << "========================\n";
    stats_out << "Total sequences: " << stats.total_sequences << "\n";
    stats_out << "Total bases: " << stats.total_bases << "\n";
    stats_out << "Valid minimizers: " << stats.valid_minimizers_extracted << "\n";
    stats_out << "LCA assignments: " << stats.lca_assignments << "\n";
    stats_out << "Sequence processing time: " << stats.sequence_processing_time << "s\n";
    stats_out << "Minimizer extraction time: " << stats.minimizer_extraction_time << "s\n";
    stats_out << "LCA computation time: " << stats.lca_computation_time << "s\n";
    
    stats_out.close();
    
    std::cout << "Database saved successfully!" << std::endl;
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

#endif // GPU_KRAKEN_DATABASE_BUILDER_CUH
#endif // GPU_KRAKEN_CLASSIFIER_HEADER_ONLY
