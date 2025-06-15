// gpu_minimizer_extraction.cuh
// Header file for GPU minimizer extraction functions
// This header should be included instead of the .cu file

#ifndef GPU_MINIMIZER_EXTRACTION_CUH
#define GPU_MINIMIZER_EXTRACTION_CUH

#include <cuda_runtime.h>
#include <cstdint>
#include "../../../include/biogpu/minimizer_extraction.h"

// Structure definitions (if not already in minimizer_extraction.h)
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

// Function declarations
bool extract_minimizers_gpu_optimized(
    const char* d_sequence_data,
    const GPUGenomeInfo* d_genome_info,
    int num_genomes,
    GPUMinimizerHit* d_minimizer_hits,
    uint32_t* d_hit_counts,
    uint32_t* total_hits,
    MinimizerParams params,
    int max_minimizers
);

bool deduplicate_minimizers_gpu(
    GPUMinimizerHit* d_minimizer_hits,
    uint32_t num_hits,
    uint32_t* final_count
);

void test_minimizer_extraction();

// Device function declarations
__device__ uint64_t encode_base(char base);
__device__ uint64_t hash_lmer(const char* sequence, int pos, int ell);
__device__ uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces, int ell);
__device__ uint64_t compute_canonical_minimizer(uint64_t hash, uint64_t xor_mask);
__device__ bool is_valid_sequence(const char* seq, int len);
__device__ uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    int kmer_pos, 
    int k, 
    int ell, 
    int spaces,
    uint64_t xor_mask
);

#endif // GPU_MINIMIZER_EXTRACTION_CUH