#include <iostream>
#include <cuda_runtime.h>
#include <cstring>

__device__ __constant__ uint8_t nucleotide_map[256];

void init_nucleotide_map() {
    uint8_t h_map[256];
    memset(h_map, 255, 256);  // Invalid by default
    h_map['A'] = h_map['a'] = 0;
    h_map['C'] = h_map['c'] = 1;
    h_map['G'] = h_map['g'] = 2;
    h_map['T'] = h_map['t'] = 3;
    cudaMemcpyToSymbol(nucleotide_map, h_map, 256);
}

__global__ void test_kernel(char* test_seq, uint8_t* results) {
    int tid = threadIdx.x;
    if (tid < 4) {
        results[tid] = nucleotide_map[test_seq[tid]];
    }
}

int main() {
    init_nucleotide_map();
    
    char test_seq[] = "ACGT";
    char* d_test_seq;
    uint8_t* d_results;
    
    cudaMalloc(&d_test_seq, 4);
    cudaMalloc(&d_results, 4);
    
    cudaMemcpy(d_test_seq, test_seq, 4, cudaMemcpyHostToDevice);
    
    test_kernel<<<1, 4>>>(d_test_seq, d_results);
    
    uint8_t results[4];
    cudaMemcpy(results, d_results, 4, cudaMemcpyDeviceToHost);
    
    std::cout << "Nucleotide map test:\n";
    for (int i = 0; i < 4; i++) {
        std::cout << test_seq[i] << " -> " << (int)results[i] << "\n";
    }
    
    cudaFree(d_test_seq);
    cudaFree(d_results);
    
    return 0;
}