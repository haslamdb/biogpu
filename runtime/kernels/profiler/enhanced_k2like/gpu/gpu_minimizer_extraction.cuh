// gpu_minimizer_extraction.cuh
// Device functions for minimizer extraction on GPU
// Provides core functionality for k-mer and minimizer extraction

#ifndef GPU_MINIMIZER_EXTRACTION_CUH
#define GPU_MINIMIZER_EXTRACTION_CUH

#include <cstdint>

// Helper functions for minimizer extraction
__device__ inline uint64_t encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // Invalid base marker
    }
}

__device__ inline uint64_t complement_base(uint64_t base) {
    return 3 - base;  // A<->T (0<->3), C<->G (1<->2)
}

// MurmurHash3 implementation (matches Kraken2)
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

// Kraken2-style minimizer extraction
__device__ inline uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    uint32_t kmer_pos,
    uint32_t k, 
    uint32_t ell, 
    uint32_t spaces, 
    uint64_t xor_mask
) {
    // Extract the k-mer
    uint64_t kmer = extract_kmer(sequence, kmer_pos, k);
    if (kmer == UINT64_MAX) return UINT64_MAX;
    
    // Get canonical form
    uint64_t canon_kmer = canonical_kmer(kmer, k);
    
    // Find minimizer within this k-mer
    uint64_t min_hash = UINT64_MAX;
    
    // Slide window of size ell across the k-mer
    for (uint32_t i = 0; i <= k - ell; i++) {
        uint64_t lmer = (canon_kmer >> (2 * (k - ell - i))) & ((1ULL << (2 * ell)) - 1);
        
        // Apply MurmurHash3 and XOR mask
        uint64_t hash = murmur_hash3(lmer ^ xor_mask);
        
        if (hash < min_hash) {
            min_hash = hash;
        }
    }
    
    return min_hash;
}

#endif // GPU_MINIMIZER_EXTRACTION_CUH