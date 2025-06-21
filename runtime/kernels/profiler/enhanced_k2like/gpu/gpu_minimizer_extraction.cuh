// gpu_minimizer_extraction.cuh
// Device functions for minimizer extraction on GPU
// Provides core functionality for k-mer and minimizer extraction

#ifndef GPU_MINIMIZER_EXTRACTION_CUH
#define GPU_MINIMIZER_EXTRACTION_CUH

#include <cstdint>

// Device function for minimizer extraction
__device__ uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    uint32_t start_pos,
    uint32_t k, 
    uint32_t ell, 
    uint32_t spaces, 
    uint64_t xor_mask
) {
    // Note: Caller must ensure start_pos + k doesn't exceed sequence length
    // This is a simplified implementation - enhance for production use
    
    // Basic k-mer hashing (simplified)
    uint64_t hash = 0;
    for (uint32_t i = 0; i < k; i++) {
        char base = sequence[start_pos + i];
        
        // Handle both uppercase and lowercase
        if (base >= 'a' && base <= 'z') {
            base = base - 'a' + 'A';
        }
        
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            return UINT64_MAX;  // Invalid base
        }
        
        // Simple 2-bit encoding
        uint64_t base_val = 0;
        if (base == 'C') base_val = 1;
        else if (base == 'G') base_val = 2;
        else if (base == 'T') base_val = 3;
        
        hash = (hash << 2) | base_val;
    }
    
    return hash ^ xor_mask;
}

#endif // GPU_MINIMIZER_EXTRACTION_CUH