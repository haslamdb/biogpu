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
    uint32_t* final_count,
    uint32_t max_allocated_hits
);

void test_minimizer_extraction();

// Device function implementations - moved to header for cross-compilation unit access
__device__ inline uint64_t encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // Invalid base marker
    }
}

__device__ inline uint64_t hash_lmer(const char* sequence, int pos, int ell) {
    uint64_t hash = 0;
    for (int i = 0; i < ell; i++) {
        uint64_t base = encode_base(sequence[pos + i]);
        if (base == 4) return UINT64_MAX;  // Invalid sequence
        hash = (hash << 2) | base;
    }
    return hash;
}

__device__ inline uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces, int ell) {
    if (spaces == 0) return hash;
    
    uint64_t masked_hash = 0;
    int out_pos = 0;
    
    // Apply spaced seed pattern: keep every (spaces+1)th position
    for (int i = 0; i < ell; i++) {
        if (i % (spaces + 1) == 0) {
            uint64_t base = (hash >> (2 * (ell - 1 - i))) & 3;
            masked_hash = (masked_hash << 2) | base;
            out_pos++;
        }
    }
    
    return masked_hash;
}

__device__ inline uint64_t compute_canonical_minimizer(uint64_t hash, uint64_t xor_mask) {
    // Apply XOR shuffling to avoid bias toward low-complexity sequences
    return hash ^ xor_mask;
}

__device__ inline bool is_valid_sequence(const char* seq, int len) {
    for (int i = 0; i < len; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

// Sliding window minimizer extraction for a single k-mer - moved to header for cross-compilation unit access
__device__ inline uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    int kmer_pos, 
    int k, 
    int ell, 
    int spaces,
    uint64_t xor_mask
) {
    uint64_t min_hash = UINT64_MAX;
    
    // Slide window of size ell across the k-mer
    for (int i = 0; i <= k - ell; i++) {
        if (!is_valid_sequence(sequence + kmer_pos + i, ell)) {
            continue;
        }
        
        uint64_t lmer_hash = hash_lmer(sequence, kmer_pos + i, ell);
        if (lmer_hash == UINT64_MAX) continue;
        
        // Apply spaced seed mask
        uint64_t masked_hash = apply_spaced_seed_mask(lmer_hash, spaces, ell);
        
        // Apply XOR shuffling
        uint64_t canonical_hash = compute_canonical_minimizer(masked_hash, xor_mask);
        
        if (canonical_hash < min_hash) {
            min_hash = canonical_hash;
        }
    }
    
    return min_hash;
}

#endif // GPU_MINIMIZER_EXTRACTION_CUH