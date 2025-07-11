#include <cuda_runtime.h>
#include <iostream>

__global__ void test_int128_kernel(bool* result) {
    // Test basic 128-bit arithmetic
    unsigned __int128 a = 1;
    a = a << 70;  // Shift by 70 bits (requires 128-bit)
    
    unsigned __int128 b = 1;
    b = b << 69;
    
    // Test comparison
    *result = (a > b);
}

int main() {
    bool* d_result;
    bool h_result = false;
    
    cudaMalloc(&d_result, sizeof(bool));
    
    test_int128_kernel<<<1, 1>>>(d_result);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Kernel launch error: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Kernel execution error: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    cudaMemcpy(&h_result, d_result, sizeof(bool), cudaMemcpyDeviceToHost);
    
    std::cout << "128-bit integer test " << (h_result ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Result: " << h_result << " (expected: 1)" << std::endl;
    
    cudaFree(d_result);
    
    return h_result ? 0 : 1;
}