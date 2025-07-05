// runtime/common/gpu/gpu_sequence_buffer.cu
#include "gpu_sequence_buffer.h"
#include <iostream>
#include <cstring>
#include <cuda_runtime.h>

#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(error) << std::endl; \
            return false; \
        } \
    } while(0)

#define CUDA_CHECK_VOID(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(error) << std::endl; \
            return; \
        } \
    } while(0)

namespace BioGPU {

GPUSequenceBuffer::GPUSequenceBuffer(size_t max_seqs, size_t max_length) 
    : max_sequences(max_seqs), max_total_length(max_length),
      current_sequences(0), current_total_length(0) {
}

GPUSequenceBuffer::~GPUSequenceBuffer() {
    free();
}

bool GPUSequenceBuffer::allocate() {
    // Free any existing allocations
    free();
    
    // Allocate device memory
    CUDA_CHECK(cudaMalloc(&d_sequences, max_total_length * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&d_offsets, (max_sequences + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_lengths, max_sequences * sizeof(int)));
    
    // Allocate pinned host memory for fast transfers
    CUDA_CHECK(cudaMallocHost(&h_pinned_sequences, max_total_length * sizeof(char)));
    CUDA_CHECK(cudaMallocHost(&h_pinned_offsets, (max_sequences + 1) * sizeof(int)));
    CUDA_CHECK(cudaMallocHost(&h_pinned_lengths, max_sequences * sizeof(int)));
    
    // Initialize offsets
    CUDA_CHECK(cudaMemset(d_offsets, 0, sizeof(int)));
    h_pinned_offsets[0] = 0;
    
    return true;
}

void GPUSequenceBuffer::free() {
    // Free device memory
    if (d_sequences) {
        cudaFree(d_sequences);
        d_sequences = nullptr;
    }
    if (d_offsets) {
        cudaFree(d_offsets);
        d_offsets = nullptr;
    }
    if (d_lengths) {
        cudaFree(d_lengths);
        d_lengths = nullptr;
    }
    if (d_headers) {
        cudaFree(d_headers);
        d_headers = nullptr;
    }
    
    // Free pinned host memory
    if (h_pinned_sequences) {
        cudaFreeHost(h_pinned_sequences);
        h_pinned_sequences = nullptr;
    }
    if (h_pinned_offsets) {
        cudaFreeHost(h_pinned_offsets);
        h_pinned_offsets = nullptr;
    }
    if (h_pinned_lengths) {
        cudaFreeHost(h_pinned_lengths);
        h_pinned_lengths = nullptr;
    }
    
    current_sequences = 0;
    current_total_length = 0;
}

bool GPUSequenceBuffer::transferBatch(const SequenceBatch& batch) {
    if (batch.empty()) {
        return true;  // Nothing to transfer
    }
    
    // Check if batch fits
    size_t total_length = batch.getTotalBases();
    if (batch.size() > max_sequences || total_length > max_total_length) {
        std::cerr << "Error: Batch too large for buffer. "
                  << "Sequences: " << batch.size() << "/" << max_sequences
                  << ", Bases: " << total_length << "/" << max_total_length << std::endl;
        return false;
    }
    
    // Prepare flat format in pinned memory
    std::vector<char> flat_sequences;
    std::vector<int> offsets;
    const_cast<SequenceBatch&>(batch).prepareFlatFormat(flat_sequences, offsets);
    
    // Copy to pinned memory
    std::memcpy(h_pinned_sequences, flat_sequences.data(), flat_sequences.size());
    std::memcpy(h_pinned_offsets, offsets.data(), offsets.size() * sizeof(int));
    std::memcpy(h_pinned_lengths, batch.lengths.data(), batch.lengths.size() * sizeof(int));
    
    // Transfer to GPU
    CUDA_CHECK(cudaMemcpy(d_sequences, h_pinned_sequences, 
                         flat_sequences.size(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_pinned_offsets, 
                         offsets.size() * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_lengths, h_pinned_lengths, 
                         batch.lengths.size() * sizeof(int), cudaMemcpyHostToDevice));
    
    current_sequences = batch.size();
    current_total_length = flat_sequences.size();
    
    return true;
}

bool GPUSequenceBuffer::transferBatchAsync(const SequenceBatch& batch, cudaStream_t stream) {
    if (batch.empty()) {
        return true;  // Nothing to transfer
    }
    
    // Check if batch fits
    size_t total_length = batch.getTotalBases();
    if (batch.size() > max_sequences || total_length > max_total_length) {
        std::cerr << "Error: Batch too large for buffer. "
                  << "Sequences: " << batch.size() << "/" << max_sequences
                  << ", Bases: " << total_length << "/" << max_total_length << std::endl;
        return false;
    }
    
    // Prepare flat format in pinned memory
    std::vector<char> flat_sequences;
    std::vector<int> offsets;
    const_cast<SequenceBatch&>(batch).prepareFlatFormat(flat_sequences, offsets);
    
    // Copy to pinned memory
    std::memcpy(h_pinned_sequences, flat_sequences.data(), flat_sequences.size());
    std::memcpy(h_pinned_offsets, offsets.data(), offsets.size() * sizeof(int));
    std::memcpy(h_pinned_lengths, batch.lengths.data(), batch.lengths.size() * sizeof(int));
    
    // Async transfer to GPU
    CUDA_CHECK(cudaMemcpyAsync(d_sequences, h_pinned_sequences, 
                              flat_sequences.size(), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_offsets, h_pinned_offsets, 
                              offsets.size() * sizeof(int), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_lengths, h_pinned_lengths, 
                              batch.lengths.size() * sizeof(int), cudaMemcpyHostToDevice, stream));
    
    current_sequences = batch.size();
    current_total_length = flat_sequences.size();
    
    return true;
}

bool GPUSequenceBuffer::resize(size_t new_max_seqs, size_t new_max_length) {
    // If sizes are sufficient, no need to resize
    if (new_max_seqs <= max_sequences && new_max_length <= max_total_length) {
        return true;
    }
    
    // Free existing memory
    free();
    
    // Update sizes
    max_sequences = std::max(new_max_seqs, max_sequences);
    max_total_length = std::max(new_max_length, max_total_length);
    
    // Reallocate with new sizes
    return allocate();
}

void GPUSequenceBuffer::clear() {
    current_sequences = 0;
    current_total_length = 0;
    
    // Reset first offset to 0
    if (d_offsets) {
        cudaMemset(d_offsets, 0, sizeof(int));
    }
    if (h_pinned_offsets) {
        h_pinned_offsets[0] = 0;
    }
}

} // namespace BioGPU