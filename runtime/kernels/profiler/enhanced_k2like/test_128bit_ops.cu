// test_128bit_ops.cu
// Test kernel to verify 128-bit integer operations and debug l-mer mask calculation

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>

// Test kernel for 128-bit integer operations
__global__ void test_128bit_operations(uint32_t k, uint32_t l) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        printf("=== Testing 128-bit Integer Operations ===\n");
        printf("k = %u, l = %u\n", k, l);
        
        // Test 1: Basic 128-bit integer creation and shifting
        unsigned __int128 test_val = 1;
        printf("\nTest 1: Basic 128-bit operations\n");
        printf("Initial value: 1\n");
        
        // Shift by 2*k bits (should work for k=35, which is 70 bits)
        test_val = test_val << (2 * k);
        printf("After shifting by %u bits: Success\n", 2 * k);
        
        // Test 2: L-mer mask calculation (the problematic part)
        printf("\nTest 2: L-mer mask calculation\n");
        printf("Calculating mask for l = %u (needs %u bits)\n", l, 2 * l);
        
        // Method 1: Direct calculation (problematic for l=31)
        if (2 * l < 64) {
            uint64_t lmer_mask_64 = (1ULL << (2 * l)) - 1;
            printf("64-bit mask calculation: 0x%016llx\n", lmer_mask_64);
        } else {
            printf("64-bit mask calculation would overflow (2*l = %u >= 64)\n", 2 * l);
        }
        
        // Method 2: Safe calculation using 128-bit integers
        unsigned __int128 lmer_mask_128;
        if (2 * l < 128) {
            lmer_mask_128 = ((unsigned __int128)1 << (2 * l)) - 1;
            // Print lower 64 bits
            uint64_t lower = (uint64_t)(lmer_mask_128 & 0xFFFFFFFFFFFFFFFFULL);
            uint64_t upper = (uint64_t)((lmer_mask_128 >> 64) & 0xFFFFFFFFFFFFFFFFULL);
            printf("128-bit mask calculation: 0x%016llx%016llx\n", upper, lower);
        } else {
            printf("Even 128-bit mask would overflow\n");
        }
        
        // Method 3: Alternative mask calculation for l=31
        if (l == 31) {
            // For l=31, we need a 62-bit mask
            // We can construct this safely
            uint64_t mask_lower = 0xFFFFFFFFFFFFFFFFULL >> 2; // 62 bits set
            printf("Alternative mask for l=31: 0x%016llx\n", mask_lower);
        }
        
        // Test 3: Verify mask extraction
        printf("\nTest 3: Testing mask extraction\n");
        unsigned __int128 test_kmer = 0x123456789ABCDEF0ULL;
        test_kmer = (test_kmer << 64) | 0xFEDCBA9876543210ULL;
        
        // Extract using safe mask
        if (l <= 32) {
            uint64_t safe_mask = (l == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * l)) - 1);
            uint64_t extracted = (uint64_t)(test_kmer & safe_mask);
            printf("Extracted l-mer with safe mask: 0x%016llx\n", extracted);
        }
        
        // Test 4: Memory alignment test
        printf("\nTest 4: Memory alignment\n");
        printf("sizeof(unsigned __int128) = %lu\n", sizeof(unsigned __int128));
        printf("alignof(unsigned __int128) = %lu\n", alignof(unsigned __int128));
        
        printf("\n=== Tests Complete ===\n");
    }
}

// Helper function to check CUDA errors
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        printf("CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
               cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

int main() {
    printf("Testing 128-bit integer operations on GPU\n");
    printf("Testing with k=35, l=31 (the problematic case)\n\n");
    
    // Get device properties
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("GPU: %s\n", prop.name);
    printf("Compute capability: %d.%d\n", prop.major, prop.minor);
    printf("Architecture: sm_%d%d\n\n", prop.major, prop.minor);
    
    // Launch test kernel
    test_128bit_operations<<<1, 1>>>(35, 31);
    
    // Check for kernel errors
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Also test with other values
    printf("\n\nTesting with k=31, l=15 (standard Kraken2 defaults)\n");
    test_128bit_operations<<<1, 1>>>(31, 15);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    
    printf("\nAll tests completed successfully!\n");
    return 0;
}