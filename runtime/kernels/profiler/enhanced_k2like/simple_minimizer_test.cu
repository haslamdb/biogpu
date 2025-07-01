// simple_minimizer_test.cu
// Simplified test to debug minimizer extraction kernel issues

#include <iostream>
#include <cuda_runtime.h>
#include <vector>
#include <iomanip>
#include "gpu/gpu_minimizer_extraction.cuh"
#include "gpu_kraken_types.h"

// Simple test kernel that extracts minimizers from a single sequence
__global__ void simple_minimizer_kernel(
    const char* sequence,
    int seq_length,
    uint64_t* minimizers,
    int* positions,
    int* num_minimizers,
    int k,
    int ell,
    uint64_t xor_mask) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Debug output from first thread
    if (tid == 0) {
        printf("Kernel started: seq_length=%d, k=%d, ell=%d\n", seq_length, k, ell);
        printf("First 50 chars: ");
        for (int i = 0; i < min(50, seq_length); i++) {
            printf("%c", sequence[i]);
        }
        printf("\n");
    }
    
    // Each thread processes one position
    if (tid < seq_length - k + 1) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, tid, k, ell, 0, xor_mask, seq_length
        );
        
        if (minimizer != UINT64_MAX) {
            // Use atomic to add to results
            int idx = atomicAdd(num_minimizers, 1);
            if (idx < 1000) { // Limit to prevent overflow
                minimizers[idx] = minimizer;
                positions[idx] = tid;
                
                // Debug output for first few minimizers
                if (idx < 5) {
                    printf("Thread %d: Found minimizer %llx at position %d\n", 
                           tid, minimizer, tid);
                }
            }
        }
    }
}

void test_sequence(const std::string& sequence, int k = 31, int ell = 31) {
    std::cout << "\n=== Testing sequence of length " << sequence.length() << " ===" << std::endl;
    std::cout << "Parameters: k=" << k << ", ell=" << ell << std::endl;
    
    // Allocate device memory
    char* d_sequence;
    uint64_t* d_minimizers;
    int* d_positions;
    int* d_num_minimizers;
    
    cudaMalloc(&d_sequence, sequence.length() + 1);
    cudaMalloc(&d_minimizers, 1000 * sizeof(uint64_t));
    cudaMalloc(&d_positions, 1000 * sizeof(int));
    cudaMalloc(&d_num_minimizers, sizeof(int));
    
    // Initialize
    cudaMemset(d_minimizers, 0, 1000 * sizeof(uint64_t));
    cudaMemset(d_positions, 0, 1000 * sizeof(int));
    cudaMemset(d_num_minimizers, 0, sizeof(int));
    
    // Copy sequence
    cudaMemcpy(d_sequence, sequence.c_str(), sequence.length() + 1, cudaMemcpyHostToDevice);
    
    // Launch kernel
    int threads = 256;
    int blocks = (sequence.length() - k + 1 + threads - 1) / threads;
    
    std::cout << "Launching kernel with " << blocks << " blocks, " << threads << " threads" << std::endl;
    
    simple_minimizer_kernel<<<blocks, threads>>>(
        d_sequence, sequence.length(),
        d_minimizers, d_positions, d_num_minimizers,
        k, ell, 0x3c8bfbb395c60474ULL
    );
    
    // Check for errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Kernel launch error: " << cudaGetErrorString(err) << std::endl;
        return;
    }
    
    cudaDeviceSynchronize();
    
    // Get results
    int num_minimizers;
    cudaMemcpy(&num_minimizers, d_num_minimizers, sizeof(int), cudaMemcpyDeviceToHost);
    
    std::cout << "Found " << num_minimizers << " minimizers" << std::endl;
    
    if (num_minimizers > 0 && num_minimizers < 1000) {
        std::vector<uint64_t> minimizers(num_minimizers);
        std::vector<int> positions(num_minimizers);
        
        cudaMemcpy(minimizers.data(), d_minimizers, num_minimizers * sizeof(uint64_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(positions.data(), d_positions, num_minimizers * sizeof(int), cudaMemcpyDeviceToHost);
        
        // Show first few
        std::cout << "\nFirst minimizers:" << std::endl;
        for (int i = 0; i < std::min(10, num_minimizers); i++) {
            std::cout << "  [" << i << "] pos=" << positions[i] 
                      << ", hash=" << std::hex << minimizers[i] << std::dec << std::endl;
        }
    }
    
    // Clean up
    cudaFree(d_sequence);
    cudaFree(d_minimizers);
    cudaFree(d_positions);
    cudaFree(d_num_minimizers);
}

int main() {
    // Check CUDA
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    
    // Test 1: Simple repetitive sequence
    std::string test1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    test_sequence(test1);
    
    // Test 2: Random-like sequence
    std::string test2 = "ACGTGCTAGCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
                        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
    test_sequence(test2);
    
    // Test 3: Poly-A sequence (low complexity)
    std::string test3 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    test_sequence(test3);
    
    // Test 4: Test with smaller k
    test_sequence(test1, 15, 15);
    
    return 0;
}
