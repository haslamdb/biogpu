// include/biogpu/minimizer_extraction.h
// Header file for GPU minimizer extraction functions

#ifndef BIOGPU_MINIMIZER_EXTRACTION_H
#define BIOGPU_MINIMIZER_EXTRACTION_H

#include <cuda_runtime.h>
#include <cstdint>

// Forward declarations
struct GPUMinimizerHit;
struct GPUGenomeInfo;

// Parameters for minimizer extraction
struct MinimizerParams {
    int k = 35;                               // k-mer length
    int ell = 31;                            // minimizer length  
    int spaces = 7;                          // spaced seed spacing
    uint64_t xor_mask = 0x123456789ABCDEFULL; // XOR shuffling constant
    
    MinimizerParams() = default;
    MinimizerParams(int k_val, int ell_val, int spaces_val) 
        : k(k_val), ell(ell_val), spaces(spaces_val) {}
};

// LCA candidate structure (for database building)
struct LCACandidate {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
};

// Host function declarations
extern "C" {
    // Main minimizer extraction function
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
    
    // Post-processing deduplication
    bool deduplicate_minimizers_gpu(
        GPUMinimizerHit* d_minimizer_hits,
        uint32_t num_hits,
        uint32_t* final_count
    );
    
    // Test and validation functions
    void test_minimizer_extraction();
    void validate_minimizer_algorithm(const char* test_sequence, int seq_len, MinimizerParams params);
}

// Device function declarations (for use in other CUDA files)
__device__ uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    int kmer_pos, 
    int k, 
    int ell, 
    int spaces,
    uint64_t xor_mask
);

__device__ uint64_t hash_lmer(const char* sequence, int pos, int ell);
__device__ uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces, int ell);
__device__ uint64_t compute_canonical_minimizer(uint64_t hash, uint64_t xor_mask);
__device__ bool is_valid_sequence(const char* seq, int len);

#endif // BIOGPU_MINIMIZER_EXTRACTION_H