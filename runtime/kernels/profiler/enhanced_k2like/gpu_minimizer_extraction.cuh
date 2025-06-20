// gpu_minimizer_extraction.cuh
// Header file for GPU minimizer extraction functions
// This header should be included instead of the .cu file

#ifndef GPU_MINIMIZER_EXTRACTION_CUH
#define GPU_MINIMIZER_EXTRACTION_CUH

#include <cuda_runtime.h>
#include <cstdint>
#include "gpu_kraken_types.h"
#include "gpu/gpu_database_kernels.h"

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

extern "C" void test_minimizer_extraction();

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

// Compute reverse complement of a base
__device__ inline uint64_t complement_base(uint64_t base) {
    // A<->T (0<->3), C<->G (1<->2)
    return 3 - base;
}

// Proper MurmurHash3 implementation (matches Kraken2)
__device__ inline uint64_t murmur_hash3(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

// Compute canonical k-mer (minimum of forward and reverse complement)
__device__ inline uint64_t canonical_kmer(uint64_t kmer, int k) {
    uint64_t rev_comp = 0;
    uint64_t temp = kmer;
    
    // Compute reverse complement
    for (int i = 0; i < k; i++) {
        uint64_t base = temp & 3;
        rev_comp = (rev_comp << 2) | complement_base(base);
        temp >>= 2;
    }
    
    // Return minimum of forward and reverse complement
    return (kmer < rev_comp) ? kmer : rev_comp;
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

// Extract k-mer from sequence at given position
__device__ inline uint64_t extract_kmer(const char* sequence, int pos, int k) {
    uint64_t kmer = 0;
    for (int i = 0; i < k; i++) {
        uint64_t base = encode_base(sequence[pos + i]);
        if (base == 4) return UINT64_MAX; // Invalid base
        kmer = (kmer << 2) | base;
    }
    return kmer;
}

// Kraken2-style minimizer extraction structure
struct MinimizerWindow {
    uint64_t minimizer_hash;
    int minimizer_pos;
    bool valid;
};

// Extract minimizer for a k-mer using Kraken2's approach
__device__ inline MinimizerWindow extract_kmer_minimizer(
    const char* sequence, 
    int kmer_pos, 
    int k, 
    int ell,
    uint64_t xor_mask
) {
    MinimizerWindow result = {UINT64_MAX, -1, false};
    
    // Extract the k-mer
    uint64_t kmer = extract_kmer(sequence, kmer_pos, k);
    if (kmer == UINT64_MAX) return result;
    
    // Get canonical form
    uint64_t canon_kmer = canonical_kmer(kmer, k);
    
    // Find minimizer within this k-mer
    uint64_t min_hash = UINT64_MAX;
    int min_pos = -1;
    
    // Slide window of size ell across the k-mer
    for (int i = 0; i <= k - ell; i++) {
        uint64_t lmer = (canon_kmer >> (2 * (k - ell - i))) & ((1ULL << (2 * ell)) - 1);
        
        // Apply MurmurHash3
        uint64_t hash = murmur_hash3(lmer ^ xor_mask);
        
        if (hash < min_hash) {
            min_hash = hash;
            min_pos = i;
        }
    }
    
    result.minimizer_hash = min_hash;
    result.minimizer_pos = kmer_pos + min_pos;
    result.valid = true;
    
    return result;
}

// Compatibility wrapper for old code
__device__ inline uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    int kmer_pos, 
    int k, 
    int ell, 
    int spaces,
    uint64_t xor_mask
) {
    // Use the new implementation
    MinimizerWindow window = extract_kmer_minimizer(sequence, kmer_pos, k, ell, xor_mask);
    return window.valid ? window.minimizer_hash : UINT64_MAX;
}

#endif // GPU_MINIMIZER_EXTRACTION_CUH