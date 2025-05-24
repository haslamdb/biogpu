#include <iostream>
#include <cuda_runtime.h>
#include <cstdlib>

#define SEQ_LEN 120
#define QRDR_INDEX 83
#define MUTANT_AA 'L'  // Simulated mutation: Ser83Leu

#define cudaCheckError(msg) { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA Error: " << msg << " - " << cudaGetErrorString(err) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
}

// GPU kernel
__global__ void detect_mutation(char* sequences, int* results, int num_seqs) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_seqs) {
        int offset = idx * SEQ_LEN;
        char aa = sequences[offset + QRDR_INDEX];
        results[idx] = (aa == MUTANT_AA) ? 1 : 0;
    }
}

// CPU fallback
void detect_mutation_cpu(const char* sequences, int* results, int num_seqs) {
    std::cout << "[INFO] Running CPU fallback version\n";
    for (int i = 0; i < num_seqs; ++i) {
        char aa = sequences[i * SEQ_LEN + QRDR_INDEX];
        results[i] = (aa == MUTANT_AA) ? 1 : 0;
    }
}

int main() {
    const int num_seqs = 100;
    char* h_sequences = new char[num_seqs * SEQ_LEN];
    int* h_results = new int[num_seqs];

    // Fill sequences with random amino acids
    for (int i = 0; i < num_seqs * SEQ_LEN; ++i) {
        h_sequences[i] = 'A' + rand() % 20;
    }

    // Inject mutations
    for (int i = 0; i < 10; ++i) {
        int idx = rand() % num_seqs;
        h_sequences[idx * SEQ_LEN + QRDR_INDEX] = MUTANT_AA;
    }

    // Detect CUDA device
    int device_count = 0;
    cudaError_t cudaStatus = cudaGetDeviceCount(&device_count);

    if (cudaStatus != cudaSuccess || device_count == 0) {
        std::cerr << "[WARNING] No CUDA-capable device detected. Falling back to CPU.\n";
        detect_mutation_cpu(h_sequences, h_results, num_seqs);
    } else {
        std::cout << "[INFO] CUDA-capable GPU found. Running GPU version.\n";

        // Allocate and copy
        char* d_sequences;
        int* d_results;
        cudaMalloc(&d_sequences, num_seqs * SEQ_LEN * sizeof(char));
        cudaCheckError("cudaMalloc d_sequences");

        cudaMalloc(&d_results, num_seqs * sizeof(int));
        cudaCheckError("cudaMalloc d_results");

        cudaMemcpy(d_sequences, h_sequences, num_seqs * SEQ_LEN * sizeof(char), cudaMemcpyHostToDevice);
        cudaCheckError("cudaMemcpy h_sequences");

        // Kernel launch
        int threadsPerBlock = 256;
        int blocksPerGrid = (num_seqs + threadsPerBlock - 1) / threadsPerBlock;
        detect_mutation<<<blocksPerGrid, threadsPerBlock>>>(d_sequences, d_results, num_seqs);
        cudaCheckError("kernel launch");

        cudaDeviceSynchronize();
        cudaCheckError("cudaDeviceSynchronize");

        // Copy back
        cudaMemcpy(h_results, d_results, num_seqs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaCheckError("cudaMemcpy results");

        // Cleanup
        cudaFree(d_sequences);
        cudaFree(d_results);
    }

    // Output
    std::cout << "Detected mutations at index " << QRDR_INDEX << ":\n";
    for (int i = 0; i < num_seqs; ++i) {
        if (h_results[i] == 1) {
            std::cout << "  Sequence " << i << ": mutation found\n";
        }
    }

    delete[] h_sequences;
    delete[] h_results;
    return 0;
}
