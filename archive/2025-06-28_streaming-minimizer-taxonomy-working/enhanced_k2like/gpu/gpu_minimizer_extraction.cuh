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

// Extract k-mer from sequence at given position with bounds checking
__device__ inline uint64_t extract_kmer(const char* sequence, int pos, int k, int seq_length) {
    // Bounds check
    if (pos < 0 || pos + k > seq_length) {
        return UINT64_MAX;  // Invalid position
    }
    
    uint64_t kmer = 0;
    for (int i = 0; i < k; i++) {
        uint64_t base = encode_base(sequence[pos + i]);
        if (base == 4) return UINT64_MAX; // Invalid base
        kmer = (kmer << 2) | base;
    }
    return kmer;
}

// Kraken2-style minimizer extraction with proper bounds checking
__device__ inline uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    uint32_t kmer_pos,
    uint32_t k, 
    uint32_t ell, 
    uint32_t spaces, 
    uint64_t xor_mask,
    uint32_t seq_length  // Add sequence length parameter
) {
    // Validate parameters
    if (k < ell || ell == 0) {
        return UINT64_MAX;  // Invalid parameters
    }
    
    // Check if we have enough sequence for a k-mer at this position
    if (kmer_pos + k > seq_length) {
        return UINT64_MAX;  // Not enough sequence
    }
    
    // Extract the k-mer with bounds checking
    uint64_t kmer = extract_kmer(sequence, kmer_pos, k, seq_length);
    if (kmer == UINT64_MAX) return UINT64_MAX;
    
    // Get canonical form
    uint64_t canon_kmer = canonical_kmer(kmer, k);
    
    // Find minimizer within this k-mer
    uint64_t min_hash = UINT64_MAX;
    
    // Slide window of size ell across the k-mer
    // Ensure we don't exceed the k-mer bounds
    uint32_t max_window_start = (k >= ell) ? (k - ell) : 0;
    
    for (uint32_t i = 0; i <= max_window_start; i++) {
        // Extract l-mer from the canonical k-mer
        // Shift right by appropriate amount to get the l-mer starting at position i
        uint64_t lmer = (canon_kmer >> (2 * (k - ell - i))) & ((1ULL << (2 * ell)) - 1);
        
        // Apply MurmurHash3 and XOR mask
        uint64_t hash = murmur_hash3(lmer ^ xor_mask);
        
        if (hash < min_hash) {
            min_hash = hash;
        }
    }
    
    return min_hash;
}

// Alternative version that takes sequence bounds into account
__device__ inline uint64_t extract_minimizer_safe(
    const char* sequence,
    uint32_t seq_start,     // Start position in global sequence buffer
    uint32_t seq_length,    // Length of this specific sequence
    uint32_t kmer_pos,      // Position within the sequence
    uint32_t k,
    uint32_t ell,
    uint32_t spaces,
    uint64_t xor_mask
) {
    // Ensure we're within bounds of this sequence
    if (kmer_pos >= seq_length || kmer_pos + k > seq_length) {
        return UINT64_MAX;
    }
    
    // Call the main function with proper bounds
    return extract_minimizer_sliding_window(
        sequence + seq_start,  // Offset to start of this sequence
        kmer_pos,
        k,
        ell,
        spaces,
        xor_mask,
        seq_length
    );
}

#endif // GPU_MINIMIZER_EXTRACTION_CUH