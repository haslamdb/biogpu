// gpu_kraken_database_builder.cu
// Minimal implementation of GPU Kraken database builder

#include "gpu_kraken_database_builder.h"
#include "gpu_kraken_types.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <cuda_runtime.h>

#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// Constructor
GPUKrakenDatabaseBuilder::GPUKrakenDatabaseBuilder(const std::string& output_dir, const ClassificationParams& config)
    : output_directory(output_dir), params(config) {
    
    // Initialize member variables
    d_sequence_data = nullptr;
    d_genome_info = nullptr;
    d_minimizer_hits = nullptr;
    d_minimizer_counts = nullptr;
    d_global_hit_counter = nullptr;
    d_lca_candidates = nullptr;
    
    // Initialize minimizer params from classification params
    minimizer_params.k = params.k;
    minimizer_params.ell = params.ell;
    minimizer_params.spaces = params.spaces;
    minimizer_params.xor_mask = 0x3c8bfbb395c60474ULL;
    
    // Initialize statistics
    memset(&stats, 0, sizeof(stats));
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_directory);
    
    std::cout << "GPUKrakenDatabaseBuilder initialized" << std::endl;
}

// Destructor
GPUKrakenDatabaseBuilder::~GPUKrakenDatabaseBuilder() {
    free_gpu_memory();
}

// Initialize CUDA context
bool GPUKrakenDatabaseBuilder::initialize_cuda_context() {
    std::cout << "Initializing CUDA context..." << std::endl;
    
    // Initialize CUDA runtime
    cudaError_t err = cudaFree(0);
    if (err != cudaSuccess) {
        std::cerr << "Failed to initialize CUDA runtime: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    // Check for CUDA devices
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        std::cerr << "No CUDA devices found" << std::endl;
        return false;
    }
    
    std::cout << "CUDA context initialized successfully" << std::endl;
    return true;
}

// Main build pipeline
bool GPUKrakenDatabaseBuilder::build_database_from_genomes(const std::string& genome_library_path, const std::string& taxonomy_path) {
    std::cout << "\n=== Starting Database Build ===" << std::endl;
    
    if (!load_genome_files(genome_library_path)) {
        std::cerr << "Failed to load genome files" << std::endl;
        return false;
    }
    
    if (!taxonomy_path.empty() && !load_taxonomy_data(taxonomy_path)) {
        std::cerr << "Warning: Failed to load taxonomy data" << std::endl;
    }
    
    if (!process_genomes_gpu()) {
        std::cerr << "Failed to process genomes on GPU" << std::endl;
        return false;
    }
    
    if (!save_database()) {
        std::cerr << "Failed to save database" << std::endl;
        return false;
    }
    
    std::cout << "\n=== Database Build Complete ===" << std::endl;
    print_build_statistics();
    
    return true;
}

// Configuration methods
void GPUKrakenDatabaseBuilder::set_batch_size(int sequences_per_batch) {
    MAX_SEQUENCE_BATCH = sequences_per_batch;
}

void GPUKrakenDatabaseBuilder::set_minimizer_capacity(int minimizers_per_batch) {
    MAX_MINIMIZERS_PER_BATCH = minimizers_per_batch;
}

void GPUKrakenDatabaseBuilder::enable_auto_memory_scaling(bool enable, size_t memory_fraction) {
    auto_scale_memory = enable;
    max_gpu_memory_usage_fraction = memory_fraction;
}

// Print statistics
void GPUKrakenDatabaseBuilder::print_build_statistics() {
    std::cout << "\n=== Build Statistics ===" << std::endl;
    std::cout << "Total sequences: " << stats.total_sequences << std::endl;
    std::cout << "Total bases: " << stats.total_bases << std::endl;
    std::cout << "Total k-mers processed: " << stats.total_kmers_processed << std::endl;
    std::cout << "Valid minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
    std::cout << "Unique minimizers: " << stats.unique_minimizers << std::endl;
    std::cout << "LCA assignments: " << stats.lca_assignments << std::endl;
    
    if (total_minimizers_capped > 0) {
        std::cout << "\nWarning: " << total_minimizers_capped << " minimizers were dropped due to capacity limits" << std::endl;
        std::cout << "Batches that hit capacity limit: " << total_batches_capped << std::endl;
    }
}

// Internal methods - minimal implementations

bool GPUKrakenDatabaseBuilder::load_genome_files(const std::string& library_path) {
    std::cout << "Loading genome files from: " << library_path << std::endl;
    
    // For now, just create some dummy data for testing
    genome_files.push_back("test_genome_1.fna");
    genome_taxon_ids.push_back(1234);
    
    stats.total_sequences = 1;
    
    return true;
}

bool GPUKrakenDatabaseBuilder::load_taxonomy_data(const std::string& taxonomy_path) {
    std::cout << "Loading taxonomy data from: " << taxonomy_path << std::endl;
    
    // Minimal implementation - just add some test data
    taxon_names[1234] = "Test organism";
    taxon_parents[1234] = 1;
    
    return true;
}

bool GPUKrakenDatabaseBuilder::process_genomes_gpu() {
    std::cout << "Processing genomes on GPU..." << std::endl;
    
    if (!allocate_gpu_memory()) {
        return false;
    }
    
    // Process a dummy batch for testing
    std::vector<std::string> test_sequences = {"ACGTACGTACGTACGT"};
    std::vector<uint32_t> test_taxons = {1234};
    
    if (!process_sequence_batch(test_sequences, test_taxons, 0)) {
        return false;
    }
    
    return true;
}

bool GPUKrakenDatabaseBuilder::save_database() {
    std::cout << "Saving database to: " << output_directory << std::endl;
    
    // Create a dummy database file
    std::string db_file = output_directory + "/kraken_db.dat";
    std::ofstream out(db_file, std::ios::binary);
    if (!out) {
        std::cerr << "Failed to create database file" << std::endl;
        return false;
    }
    
    // Write header
    uint32_t version = 1;
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    
    // Write some statistics
    uint64_t num_minimizers = all_lca_candidates.size();
    out.write(reinterpret_cast<const char*>(&num_minimizers), sizeof(num_minimizers));
    
    out.close();
    std::cout << "Database saved successfully" << std::endl;
    
    return true;
}

bool GPUKrakenDatabaseBuilder::allocate_gpu_memory() {
    std::cout << "Allocating GPU memory..." << std::endl;
    
    size_t sequence_buffer_size = 1024 * 1024 * 100; // 100MB for sequences
    size_t genome_info_size = MAX_SEQUENCE_BATCH * sizeof(GPUGenomeInfo);
    size_t minimizer_hits_size = MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit);
    size_t minimizer_counts_size = MAX_SEQUENCE_BATCH * sizeof(uint32_t);
    size_t lca_candidates_size = MAX_MINIMIZERS_PER_BATCH * sizeof(LCACandidate);
    
    CUDA_CHECK(cudaMalloc(&d_sequence_data, sequence_buffer_size));
    CUDA_CHECK(cudaMalloc(&d_genome_info, genome_info_size));
    CUDA_CHECK(cudaMalloc(&d_minimizer_hits, minimizer_hits_size));
    CUDA_CHECK(cudaMalloc(&d_minimizer_counts, minimizer_counts_size));
    CUDA_CHECK(cudaMalloc(&d_global_hit_counter, sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_lca_candidates, lca_candidates_size));
    
    std::cout << "GPU memory allocated successfully" << std::endl;
    return true;
}

void GPUKrakenDatabaseBuilder::free_gpu_memory() {
    if (d_sequence_data) cudaFree(d_sequence_data);
    if (d_genome_info) cudaFree(d_genome_info);
    if (d_minimizer_hits) cudaFree(d_minimizer_hits);
    if (d_minimizer_counts) cudaFree(d_minimizer_counts);
    if (d_global_hit_counter) cudaFree(d_global_hit_counter);
    if (d_lca_candidates) cudaFree(d_lca_candidates);
    
    d_sequence_data = nullptr;
    d_genome_info = nullptr;
    d_minimizer_hits = nullptr;
    d_minimizer_counts = nullptr;
    d_global_hit_counter = nullptr;
    d_lca_candidates = nullptr;
}

void GPUKrakenDatabaseBuilder::check_and_adjust_memory_safe() {
    if (!auto_scale_memory) return;
    
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    
    size_t target_usage = (total_mem * max_gpu_memory_usage_fraction) / 100;
    size_t current_usage = total_mem - free_mem;
    
    if (current_usage > target_usage) {
        // Reduce batch sizes
        MAX_SEQUENCE_BATCH = std::max(1, MAX_SEQUENCE_BATCH / 2);
        MAX_MINIMIZERS_PER_BATCH = std::max(1000000, MAX_MINIMIZERS_PER_BATCH / 2);
        
        std::cout << "Adjusted batch sizes due to memory constraints:" << std::endl;
        std::cout << "  Sequences per batch: " << MAX_SEQUENCE_BATCH << std::endl;
        std::cout << "  Minimizers per batch: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
    }
}

void GPUKrakenDatabaseBuilder::validate_memory_settings() {
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    
    std::cout << "GPU Memory: " << free_mem / (1024*1024) << " MB free / " 
              << total_mem / (1024*1024) << " MB total" << std::endl;
}

bool GPUKrakenDatabaseBuilder::process_sequence_batch(const std::vector<std::string>& sequences, 
                                                     const std::vector<uint32_t>& taxon_ids, 
                                                     int batch_offset) {
    std::cout << "Processing batch of " << sequences.size() << " sequences" << std::endl;
    
    // For now, just update statistics
    for (const auto& seq : sequences) {
        stats.total_bases += seq.length();
        if (seq.length() >= minimizer_params.k) {
            stats.total_kmers_processed += seq.length() - minimizer_params.k + 1;
        }
    }
    
    // Create a dummy LCA candidate
    LCACandidate dummy_candidate;
    dummy_candidate.minimizer_hash = 0x1234567890ABCDEFULL;
    dummy_candidate.lca_taxon = 1234;
    dummy_candidate.genome_count = 1;
    dummy_candidate.uniqueness_score = 1.0f;
    
    all_lca_candidates.push_back(dummy_candidate);
    stats.valid_minimizers_extracted++;
    stats.lca_assignments++;
    
    return true;
}

bool GPUKrakenDatabaseBuilder::extract_minimizers_gpu(const char* d_sequences, 
                                                     const GPUGenomeInfo* d_genomes, 
                                                     int num_sequences, 
                                                     GPUMinimizerHit* d_hits, 
                                                     uint32_t* total_hits) {
    // Minimal implementation - just set some dummy values
    *total_hits = 1;
    return true;
}

bool GPUKrakenDatabaseBuilder::compute_lca_assignments_gpu(const GPUMinimizerHit* d_hits, 
                                                          int num_hits, 
                                                          LCACandidate* d_candidates, 
                                                          int* num_candidates) {
    // Minimal implementation
    *num_candidates = 1;
    return true;
}

std::vector<std::string> GPUKrakenDatabaseBuilder::load_sequences_from_fasta(const std::string& fasta_path) {
    std::vector<std::string> sequences;
    // Minimal implementation - return dummy sequence
    sequences.push_back("ACGTACGTACGTACGT");
    return sequences;
}

uint32_t GPUKrakenDatabaseBuilder::extract_taxon_from_filename(const std::string& filename) {
    // Minimal implementation
    return 1234;
}

uint32_t GPUKrakenDatabaseBuilder::compute_simple_lca_host(uint32_t taxon1, uint32_t taxon2) {
    // Minimal LCA implementation - just return the smaller taxon ID
    return std::min(taxon1, taxon2);
}