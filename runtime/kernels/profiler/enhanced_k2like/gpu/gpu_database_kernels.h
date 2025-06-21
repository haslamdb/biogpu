// gpu/gpu_database_kernels.h
// CUDA kernel declarations for GPU Kraken database builder
// Contains all GPU kernel interfaces and device function declarations

#ifndef GPU_DATABASE_KERNELS_H
#define GPU_DATABASE_KERNELS_H

#include <cuda_runtime.h>
#include <cstdint>

// Include common types
#include "../gpu_kraken_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// ===========================
// GPU Data Structures
// ===========================

// Note: GPUGenomeInfo, GPUMinimizerHit, LCACandidate and MinimizerParams 
// are all defined in gpu_kraken_types.h

// GPU taxonomy node structure
struct GPUTaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    uint8_t depth;
    uint16_t padding;
    
    __host__ __device__ GPUTaxonomyNode() : 
        taxon_id(0), parent_id(0), rank(0), depth(0), padding(0) {}
};

// GPU Compact Hash Table structure
struct GPUCompactHashTable {
    uint32_t* hash_cells;
    uint32_t hash_mask;
    uint32_t lca_bits;
    uint32_t max_taxon_id;
    
    __host__ __device__ GPUCompactHashTable() : 
        hash_cells(nullptr), hash_mask(0), lca_bits(0), max_taxon_id(0) {}
};

// ===========================
// Device Function Declarations
// ===========================

// Hashing functions
__device__ uint32_t jenkins_hash_gpu(uint64_t key);
__device__ uint32_t compute_compact_hash_gpu(uint64_t minimizer_hash);
__device__ uint32_t lookup_lca_gpu(const GPUCompactHashTable* cht, uint64_t minimizer_hash);
__device__ uint64_t murmur_hash3_64(const void* key, int len, uint32_t seed);
__device__ uint64_t compute_kmer_hash(const char* sequence, int k, uint64_t mask);
__device__ uint64_t compute_spaced_kmer_hash(const char* sequence, int k, int spaces, uint64_t mask);

// Minimizer extraction device functions
__device__ uint64_t extract_minimizer_sliding_window(const char* sequence, uint32_t start_pos, 
                                                    uint32_t k, uint32_t ell, uint32_t spaces, 
                                                    uint64_t xor_mask);
__device__ bool is_valid_minimizer(uint64_t minimizer_hash, uint64_t min_clear_hash);
__device__ bool has_valid_dna_bases(const char* sequence, int length);
__device__ bool has_valid_bases_device(const char* seq, int len);
__device__ bool validate_sequence_device(const char* sequence, int length);
__device__ void atomic_add_safe_device(uint32_t* address, uint32_t value);

// Sequence validation
__device__ bool validate_dna_sequence_gpu(const char* sequence, uint32_t length);
__device__ char complement_base(char base);
__device__ uint64_t reverse_complement_hash(uint64_t forward_hash, int k);

// LCA computation device functions
__device__ uint32_t compute_simple_lca_gpu(uint32_t taxon1, uint32_t taxon2);
__device__ uint32_t find_lca_simple_device(uint32_t taxon1, uint32_t taxon2,
                                          const uint32_t* parent_lookup,
                                          uint32_t max_taxon_id);
__device__ uint32_t compute_lca_from_taxonomy_gpu(uint32_t taxon1, uint32_t taxon2, 
                                                 const GPUTaxonomyNode* taxonomy_nodes, 
                                                 int num_nodes);

// Utility device functions
__device__ void atomic_add_float(float* address, float value);
__device__ uint32_t hash_combine(uint32_t a, uint32_t b);
__device__ bool sequences_equal(const char* seq1, const char* seq2, int length);

// ===========================
// CUDA Kernel Declarations
// ===========================

// Minimizer extraction kernels
__global__ void extract_minimizers_kraken2_improved_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts_per_genome,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    uint64_t min_clear_hash_value,
    uint64_t toggle_mask,
    int max_minimizers
);

__global__ void extract_minimizers_basic_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* total_hits,
    MinimizerParams params
);

__global__ void extract_minimizers_spaced_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts_per_genome,
    MinimizerParams params,
    int max_minimizers
);

// LCA computation kernels
__global__ void compute_lca_assignments_kernel(
    const GPUMinimizerHit* minimizer_hits,
    int num_hits,
    LCACandidate* lca_candidates,
    int* num_candidates
);

__global__ void convert_hits_to_candidates_fixed(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates
);

__global__ void merge_duplicate_minimizers_kernel(
    LCACandidate* candidates,
    int num_candidates,
    LCACandidate* merged_candidates,
    int* num_merged
);

// Advanced LCA computation with taxonomy
__global__ void compute_phylogenetic_lca_kernel(
    const GPUMinimizerHit* hits,
    int num_hits,
    const GPUTaxonomyNode* taxonomy_nodes,
    int num_taxonomy_nodes,
    LCACandidate* lca_candidates,
    int* num_candidates
);

// Memory and validation kernels
__global__ void initialize_gpu_memory_kernel(
    char* sequence_buffer,
    size_t sequence_buffer_size,
    GPUGenomeInfo* genome_info_buffer,
    size_t num_genomes,
    GPUMinimizerHit* minimizer_buffer,
    size_t num_minimizers
);

__global__ void validate_gpu_memory_kernel(
    const char* sequence_buffer,
    const GPUGenomeInfo* genome_info_buffer,
    int num_genomes,
    bool* validation_result
);

__global__ void clear_buffers_kernel(
    GPUMinimizerHit* minimizer_hits,
    uint32_t* counters,
    int buffer_size
);

// Quality control and filtering kernels
__global__ void filter_low_quality_minimizers_kernel(
    GPUMinimizerHit* minimizer_hits,
    int num_hits,
    float quality_threshold,
    GPUMinimizerHit* filtered_hits,
    int* num_filtered
);

__global__ void remove_duplicate_minimizers_kernel(
    GPUMinimizerHit* minimizer_hits,
    int num_hits,
    GPUMinimizerHit* unique_hits,
    int* num_unique
);

__global__ void calculate_uniqueness_scores_kernel(
    LCACandidate* candidates,
    int num_candidates,
    const uint32_t* genome_counts_per_taxon,
    int num_taxa
);

// Statistics and analysis kernels
__global__ void count_minimizers_per_taxon_kernel(
    const GPUMinimizerHit* hits,
    int num_hits,
    uint32_t* taxon_counts,
    int max_taxon_id
);

__global__ void calculate_coverage_statistics_kernel(
    const GPUGenomeInfo* genome_info,
    const uint32_t* hit_counts_per_genome,
    int num_genomes,
    float* coverage_stats
);

// ===========================
// Host-Side Kernel Launchers
// ===========================

// Host wrapper functions for kernel launches - matching implementation signatures
bool launch_minimizer_extraction_kernel(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint32_t* total_hits_output);

bool launch_improved_minimizer_kernel(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint64_t min_clear_hash_value,
    uint64_t toggle_mask,
    uint32_t* total_hits_output);

bool launch_lca_computation_kernel(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates,
    int* num_candidates,
    const uint32_t* parent_lookup = nullptr,
    uint32_t max_taxon_id = 0);

// Utility functions
LaunchConfig calculate_optimal_launch_config(int num_elements, int threads_per_block = 256, size_t shared_memory = 0);
bool check_kernel_execution_errors(const char* kernel_name);
void print_kernel_launch_info(const LaunchConfig& config, const char* kernel_name);
bool validate_kernel_parameters(const GPUBatchData& batch_data);

// Legacy signatures for compatibility
cudaError_t launch_memory_initialization_kernel(
    char* d_sequence_buffer,
    size_t sequence_buffer_size,
    GPUGenomeInfo* d_genome_info,
    size_t num_genomes,
    GPUMinimizerHit* d_minimizer_buffer,
    size_t num_minimizers,
    cudaStream_t stream = 0
);

cudaError_t launch_memory_validation_kernel(
    const char* d_sequence_buffer,
    const GPUGenomeInfo* d_genome_info,
    int num_genomes,
    bool* h_validation_result,
    cudaStream_t stream = 0
);

// Convenience host functions
bool extract_minimizers_gpu_wrapper(
    const char* d_sequences,
    const GPUGenomeInfo* d_genomes,
    int num_sequences,
    GPUMinimizerHit* d_hits,
    uint32_t* total_hits,
    const MinimizerParams& params,
    int max_minimizers
);

bool compute_lca_assignments_gpu_wrapper(
    const GPUMinimizerHit* d_hits,
    int num_hits,
    LCACandidate* d_candidates,
    int* num_candidates
);

// ===========================
// GPU Configuration and Limits
// ===========================

// GPU configuration constants
namespace GPUConfig {
    constexpr int MAX_THREADS_PER_BLOCK = 1024;
    constexpr int DEFAULT_THREADS_PER_BLOCK = 256;
    constexpr int MAX_BLOCKS_PER_GRID = 65535;
    constexpr int WARP_SIZE = 32;
    constexpr int MAX_SHARED_MEMORY_PER_BLOCK = 49152; // 48KB
    
    // Application-specific limits
    constexpr int MAX_SEQUENCE_LENGTH = 50000000;      // 50MB per sequence
    constexpr int MAX_GENOMES_PER_BATCH = 100;         // Maximum genomes per batch
    constexpr int MAX_MINIMIZERS_PER_GENOME = 1000000; // 1M minimizers per genome
    constexpr int MAX_K_VALUE = 31;                    // Maximum k-mer size
    constexpr int MIN_K_VALUE = 15;                    // Minimum k-mer size
}

// GPU memory management helpers
namespace GPUMemoryHelpers {
    // Calculate optimal grid dimensions
    dim3 calculate_grid_dimensions(int total_elements, int threads_per_block = GPUConfig::DEFAULT_THREADS_PER_BLOCK);
    
    // Calculate shared memory requirements
    size_t calculate_shared_memory_size(int threads_per_block, size_t per_thread_data);
    
    // Validate GPU capabilities
    bool validate_gpu_capabilities();
    bool check_compute_capability(int required_major, int required_minor);
    
    // Memory alignment helpers
    size_t align_memory_size(size_t size, size_t alignment = 256);
    bool is_memory_aligned(void* ptr, size_t alignment = 256);
}

// Error checking macros for kernels
#define CUDA_CHECK_KERNEL(call) do { \
    call; \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA kernel error at %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err)); \
        return false; \
    } \
} while(0)

#define CUDA_CHECK_LAUNCH(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA launch error at %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err)); \
        return err; \
    } \
} while(0)

// ===========================
// Performance Profiling Support
// ===========================

struct KernelPerformanceMetrics {
    float execution_time_ms;
    size_t memory_throughput_gb_per_s;
    int occupancy_percentage;
    size_t shared_memory_used;
    size_t global_memory_used;
};

// Kernel profiling functions
KernelPerformanceMetrics profile_minimizer_extraction_kernel(
    const char* d_sequence_data,
    const GPUGenomeInfo* d_genome_info,
    int num_genomes,
    const MinimizerParams& params
);

void print_kernel_performance_summary(const KernelPerformanceMetrics& metrics, const std::string& kernel_name);

#ifdef __cplusplus
}
#endif

#endif // GPU_DATABASE_KERNELS_H
