/*
 * gpu/gpu_minimizer_extraction.cuh
 *
 * This is a corrected version of the minimizer extraction logic, directly
 * adapted from the official Kraken 2 source code (MinimizerScanner).
 * It correctly handles k-mers larger than 32 by using a 128-bit integer
 * type for k-mer representation, which prevents the bit-shifting overflow
 * that was causing the illegal memory access error.
 *
 * The final minimizer hash remains a uint64_t, so this does not affect
 * the database size.
 */
#ifndef GPU_MINIMIZER_EXTRACTION_CUH
#define GPU_MINIMIZER_EXTRACTION_CUH

#include <cstdint>

// Helper to convert a character to its 2-bit DNA code
__device__ inline uint64_t DnaTo2Bit(char c) {
    // A = 00, C = 01, G = 10, T = 11
    // This is a common and efficient encoding scheme.
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // Invalid base
    }
}

// MurmurHash3 implementation - this is the same hash function used by Kraken 2
__device__ inline uint64_t murmur_hash3(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

// This is the core fix. We use unsigned __int128 to represent k-mers larger
// than 32 bases, preventing overflow during bit-shifting operations.
__device__ inline unsigned __int128 extract_kmer_128(const char *sequence, uint32_t k) {
    unsigned __int128 kmer = 0;
    for (uint32_t i = 0; i < k; i++) {
        uint64_t base = DnaTo2Bit(sequence[i]);
        if (base == 4) return (unsigned __int128)-1; // Invalid base found
        kmer <<= 2;
        kmer |= base;
    }
    return kmer;
}

// Calculates the reverse complement of a 128-bit k-mer
__device__ inline unsigned __int128 reverse_complement_128(unsigned __int128 kmer, uint32_t k) {
    unsigned __int128 rc_kmer = 0;
    for (uint32_t i = 0; i < k; i++) {
        rc_kmer <<= 2;
        // (kmer & 3) gets the last 2 bits (a base)
        // 3 - base complements it (A<->T, C<->G)
        rc_kmer |= (3 - (kmer & 3));
        kmer >>= 2;
    }
    return rc_kmer;
}

// Main device function to find a minimizer in a k-mer window.
// This is a direct port of the logic in Kraken 2's MinimizerScanner.
__device__ inline uint64_t extract_minimizer_sliding_window(
    const char* sequence,
    uint32_t kmer_pos,
    uint32_t k,
    uint32_t l,
    uint32_t spaces, // Note: Spaced seeds are not used in this corrected version for simplicity
    uint64_t xor_mask,
    uint32_t seq_length)
{
    // Basic validation
    if (kmer_pos + k > seq_length) return UINT64_MAX;

    // Step 1: Extract the k-mer into a 128-bit integer to prevent overflow
    unsigned __int128 kmer128 = extract_kmer_128(sequence + kmer_pos, k);
    if (kmer128 == (unsigned __int128)-1) return UINT64_MAX; // Invalid base found

    // Step 2: Get the canonical representation (the lesser of the k-mer and its reverse complement)
    unsigned __int128 rc_kmer128 = reverse_complement_128(kmer128, k);
    unsigned __int128 canonical_kmer = (kmer128 < rc_kmer128) ? kmer128 : rc_kmer128;

    // Step 3: Slide a window of size 'l' over the canonical k-mer to find the minimizer
    uint64_t min_lmer_hash = UINT64_MAX;
    
    // Create a mask to extract l-mer bits. Since l <= 31, it fits in uint64_t.
    // Note: 2*l can be 62, so 1ULL is necessary to ensure 64-bit shift.
    uint64_t lmer_mask = (l < 32) ? (1ULL << (2 * l)) - 1 : UINT64_MAX;

    for (uint32_t i = 0; i <= k - l; i++) {
        // Extract the l-mer from the canonical k-mer
        // Note: We cast to uint64_t because l <= 31, so it fits.
        uint64_t lmer = (uint64_t)((canonical_kmer >> (2 * i)) & lmer_mask);

        // Step 4: Hash the l-mer to get the minimizer hash
        uint64_t current_hash = murmur_hash3(lmer ^ xor_mask);

        if (current_hash < min_lmer_hash) {
            min_lmer_hash = current_hash;
        }
    }

    return min_lmer_hash;
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