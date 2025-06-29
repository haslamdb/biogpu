#include <iostream>
#include <cuda_runtime.h>
#include "gpu/gpu_minimizer_extraction.cuh"

__global__ void test_x_sequence_kernel(uint64_t* result) {
    // Test sequence with many x's
    const char* test_seq = "ATCGATCGATCGxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxATCGATCG";
    int seq_len = 94;  // Length of above string
    
    // Try to extract minimizer at position 10 (should hit x's)
    uint64_t minimizer = extract_minimizer_sliding_window(
        test_seq, 10, 31, 31, 7, 0x3c8bfbb395c60474ULL, seq_len
    );
    
    result[0] = minimizer;
    
    // Try at position 0 (should work)
    minimizer = extract_minimizer_sliding_window(
        test_seq, 0, 31, 31, 7, 0x3c8bfbb395c60474ULL, seq_len
    );
    
    result[1] = minimizer;
}

int main() {
    uint64_t* d_result;
    uint64_t h_result[2];
    
    cudaMalloc(&d_result, 2 * sizeof(uint64_t));
    
    test_x_sequence_kernel<<<1, 1>>>(d_result);
    cudaDeviceSynchronize();
    
    auto err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    cudaMemcpy(h_result, d_result, 2 * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
    std::cout << "Minimizer at pos 10 (with x's): " << std::hex << h_result[0] << std::endl;
    std::cout << "Minimizer at pos 0 (valid): " << std::hex << h_result[1] << std::endl;
    std::cout << "UINT64_MAX: " << std::hex << UINT64_MAX << std::endl;
    
    cudaFree(d_result);
    return 0;
}