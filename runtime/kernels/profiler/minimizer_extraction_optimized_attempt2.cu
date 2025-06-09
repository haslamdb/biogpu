// minimizer_extraction.cu
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdint>
#include <vector>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "minimizer_common.h"
#include "minimizer_extractor.h"

// Constants
#define BATCH_SIZE 10000
#define MAX_READ_LENGTH 300
#define MAX_MINIMIZERS_PER_READ 50

// Bit encoding for nucleotides
__device__ __constant__ uint8_t nucleotide_map[256];

// Host-side initialization
void init_nucleotide_map() {
    uint8_t h_map[256];
    memset(h_map, 255, 256);  // Invalid by default
    h_map['A'] = h_map['a'] = 0;
    h_map['C'] = h_map['c'] = 1;
    h_map['G'] = h_map['g'] = 2;
    h_map['T'] = h_map['t'] = 3;
    cudaError_t err = cudaMemcpyToSymbol(nucleotide_map, h_map, 256);
    if (err != cudaSuccess) {
        std::cerr << "Failed to initialize nucleotide map: " << cudaGetErrorString(err) << std::endl;
    }
}

// MurmurHash3 finalizer for better hash distribution
__device__ inline uint64_t murmur_hash3_finalizer(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccd;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53;
    key ^= key >> 33;
    return key;
}

// Reverse complement of a k-mer
__device__ uint64_t reverse_complement(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        uint8_t base = (kmer >> (2 * i)) & 3;
        uint8_t complement = 3 - base;  // A<->T, C<->G
        rc = (rc << 2) | complement;
    }
    return rc;
}

// GPU kernel for extracting minimizers from reads
__global__ void extract_minimizers_kernel(
    const char* sequences,           // Concatenated sequences
    const uint32_t* sequence_offsets,// Start position of each sequence
    const uint32_t* sequence_lengths,// Length of each sequence
    Minimizer* minimizers,           // Output minimizers
    uint32_t* minimizer_counts,      // Number of minimizers per read
    const int num_reads,
    const int k,                     // k-mer size
    const int m                      // minimizer window size
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const uint32_t seq_start = sequence_offsets[tid];
    const uint32_t seq_length = sequence_lengths[tid];
    
    if (seq_length < k) {
        minimizer_counts[tid] = 0;
        return;
    }
    
    // Local storage for window k-mers
    uint64_t window_kmers[32];  // Assuming m <= 32
    uint32_t window_positions[32];
    bool window_reverse[32];
    int window_size = 0;
    
    // Variables for rolling hash
    uint64_t forward_kmer = 0;
    uint64_t reverse_kmer = 0;
    const uint64_t kmer_mask = (1ULL << (2 * k)) - 1;
    
    // Initialize first k-mer
    int valid_bases = 0;
    for (int i = 0; i < k - 1; i++) {
        char base = sequences[seq_start + i];
        uint8_t encoded = nucleotide_map[base];
        
        if (encoded < 4) {
            forward_kmer = ((forward_kmer << 2) | encoded) & kmer_mask;
            reverse_kmer = (reverse_kmer >> 2) | ((uint64_t)(3 - encoded) << (2 * (k - 1)));
            valid_bases++;
        } else {
            // Invalid base, reset
            forward_kmer = 0;
            reverse_kmer = 0;
            valid_bases = 0;
        }
    }
    
    // Process sequence with sliding window
    int minimizer_count = 0;
    uint64_t prev_minimizer = UINT64_MAX;
    
    for (uint32_t pos = k - 1; pos < seq_length; pos++) {
        char base = sequences[seq_start + pos];
        uint8_t encoded = nucleotide_map[base];
        
        if (encoded < 4) {
            // Update k-mers
            forward_kmer = ((forward_kmer << 2) | encoded) & kmer_mask;
            reverse_kmer = (reverse_kmer >> 2) | ((uint64_t)(3 - encoded) << (2 * (k - 1)));
            valid_bases++;
            
            if (valid_bases >= k) {
                // Choose canonical k-mer (smaller of forward/reverse)
                uint64_t canonical_kmer = min(forward_kmer, reverse_kmer);
                uint64_t kmer_hash = murmur_hash3_finalizer(canonical_kmer);
                
                // Add to window
                window_kmers[window_size % m] = kmer_hash;
                window_positions[window_size % m] = pos - k + 1;
                window_reverse[window_size % m] = (reverse_kmer < forward_kmer);
                window_size++;
                
                // Find minimizer in current window
                if (window_size >= m) {
                    uint64_t min_hash = UINT64_MAX;
                    int min_idx = 0;
                    
                    int window_start = (window_size > m) ? (window_size - m) : 0;
                    for (int i = 0; i < m && window_start + i < window_size; i++) {
                        int idx = (window_start + i) % m;
                        if (window_kmers[idx] < min_hash) {
                            min_hash = window_kmers[idx];
                            min_idx = idx;
                        }
                    }
                    
                    // Only store if different from previous minimizer
                    if (min_hash != prev_minimizer && minimizer_count < MAX_MINIMIZERS_PER_READ) {
                        Minimizer* out_minimizer = &minimizers[tid * MAX_MINIMIZERS_PER_READ + minimizer_count];
                        out_minimizer->hash = min_hash;
                        out_minimizer->position = window_positions[min_idx];
                        out_minimizer->is_reverse = window_reverse[min_idx];
                        minimizer_count++;
                        prev_minimizer = min_hash;
                    }
                }
            }
        } else {
            // Invalid base, reset
            forward_kmer = 0;
            reverse_kmer = 0;
            valid_bases = 0;
            window_size = 0;
            prev_minimizer = UINT64_MAX;
        }
    }
    
    minimizer_counts[tid] = minimizer_count;
}

// Host wrapper class implementation
MinimizerExtractor::MinimizerExtractor(int k_size, int window_size) 
    : k(k_size), m(window_size), allocated_reads(0), allocated_sequence_length(0),
      d_sequences(nullptr), d_sequence_offsets(nullptr),
      d_sequence_lengths(nullptr), d_minimizers(nullptr),
      d_minimizer_counts(nullptr),
      h_sequences_pinned(nullptr), h_offsets_pinned(nullptr),
      h_lengths_pinned(nullptr), h_minimizers_pinned(nullptr),
      h_counts_pinned(nullptr) {
    
    init_nucleotide_map();
    
    // Create CUDA streams
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    
    // Pre-allocate pinned memory for typical batch
    size_t typical_batch = BATCH_SIZE;
    size_t typical_seq_length = typical_batch * 200;  // Assume 200bp average
    
    cudaMallocHost(&h_sequences_pinned, typical_seq_length);
    cudaMallocHost(&h_offsets_pinned, (typical_batch + 1) * sizeof(uint32_t));
    cudaMallocHost(&h_lengths_pinned, typical_batch * sizeof(uint32_t));
    cudaMallocHost(&h_minimizers_pinned, typical_batch * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
    cudaMallocHost(&h_counts_pinned, typical_batch * sizeof(uint32_t));
    
    allocated_sequence_length = typical_seq_length;
}

MinimizerExtractor::~MinimizerExtractor() {
    // Free device memory
    if (d_sequences) cudaFree(d_sequences);
    if (d_sequence_offsets) cudaFree(d_sequence_offsets);
    if (d_sequence_lengths) cudaFree(d_sequence_lengths);
    if (d_minimizers) cudaFree(d_minimizers);
    if (d_minimizer_counts) cudaFree(d_minimizer_counts);
    
    // Free pinned memory
    if (h_sequences_pinned) cudaFreeHost(h_sequences_pinned);
    if (h_offsets_pinned) cudaFreeHost(h_offsets_pinned);
    if (h_lengths_pinned) cudaFreeHost(h_lengths_pinned);
    if (h_minimizers_pinned) cudaFreeHost(h_minimizers_pinned);
    if (h_counts_pinned) cudaFreeHost(h_counts_pinned);
    
    // Destroy streams
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
}

void MinimizerExtractor::allocate_device_memory(size_t num_reads, size_t total_sequence_length) {
    // Only reallocate if we need more space
    if (num_reads > allocated_reads || d_sequences == nullptr) {
        // Free old memory
        if (d_sequences) cudaFree(d_sequences);
        if (d_sequence_offsets) cudaFree(d_sequence_offsets);
        if (d_sequence_lengths) cudaFree(d_sequence_lengths);
        if (d_minimizers) cudaFree(d_minimizers);
        if (d_minimizer_counts) cudaFree(d_minimizer_counts);
        
        // Allocate with some extra space to avoid frequent reallocation
        size_t alloc_reads = num_reads * 1.5;  // 50% extra
        size_t alloc_seq_length = total_sequence_length * 1.5;
        
        cudaMalloc(&d_sequences, alloc_seq_length);
        cudaMalloc(&d_sequence_offsets, (alloc_reads + 1) * sizeof(uint32_t));
        cudaMalloc(&d_sequence_lengths, alloc_reads * sizeof(uint32_t));
        cudaMalloc(&d_minimizers, alloc_reads * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
        cudaMalloc(&d_minimizer_counts, alloc_reads * sizeof(uint32_t));
        
        allocated_reads = alloc_reads;
        
        // Check for allocation errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA allocation failed: ") + cudaGetErrorString(error));
        }
    }
}

std::vector<std::vector<Minimizer>> MinimizerExtractor::extract_minimizers(
    const std::vector<std::string>& sequences
) {
    size_t num_reads = sequences.size();
    
    // Calculate total sequence length
    size_t total_length = 0;
    for (const auto& seq : sequences) {
        total_length += seq.length();
    }
    
    // Reallocate pinned memory if needed
    if (total_length > allocated_sequence_length) {
        if (h_sequences_pinned) cudaFreeHost(h_sequences_pinned);
        cudaMallocHost(&h_sequences_pinned, total_length * 1.5);
        allocated_sequence_length = total_length * 1.5;
    }
    
    // Reallocate other pinned buffers if needed
    if (num_reads > allocated_reads) {
        if (h_offsets_pinned) cudaFreeHost(h_offsets_pinned);
        if (h_lengths_pinned) cudaFreeHost(h_lengths_pinned);
        if (h_minimizers_pinned) cudaFreeHost(h_minimizers_pinned);
        if (h_counts_pinned) cudaFreeHost(h_counts_pinned);
        
        size_t alloc_reads = num_reads * 1.5;
        cudaMallocHost(&h_offsets_pinned, (alloc_reads + 1) * sizeof(uint32_t));
        cudaMallocHost(&h_lengths_pinned, alloc_reads * sizeof(uint32_t));
        cudaMallocHost(&h_minimizers_pinned, alloc_reads * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
        cudaMallocHost(&h_counts_pinned, alloc_reads * sizeof(uint32_t));
        // Note: allocated_reads is updated by allocate_device_memory
    }
    
    // Copy to pinned memory
    uint32_t current_offset = 0;
    for (size_t i = 0; i < num_reads; i++) {
        h_offsets_pinned[i] = current_offset;
        h_lengths_pinned[i] = sequences[i].length();
        memcpy(h_sequences_pinned + current_offset, 
               sequences[i].data(), 
               sequences[i].length());
        current_offset += sequences[i].length();
    }
    h_offsets_pinned[num_reads] = current_offset;
    
    // Allocate device memory
    allocate_device_memory(num_reads, current_offset);
    
    // Async copy to device
    cudaMemcpyAsync(d_sequences, h_sequences_pinned, 
                    current_offset, cudaMemcpyHostToDevice, stream1);
    cudaMemcpyAsync(d_sequence_offsets, h_offsets_pinned, 
                    (num_reads + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice, stream1);
    cudaMemcpyAsync(d_sequence_lengths, h_lengths_pinned, 
                    num_reads * sizeof(uint32_t), cudaMemcpyHostToDevice, stream1);
    
    // Wait for transfers to complete
    cudaStreamSynchronize(stream1);
    
    // Launch kernel with optimized configuration
    int block_size = 128;  // Reduced from 256 for better occupancy
    int num_blocks = (num_reads + block_size - 1) / block_size;
    
    // Use dynamic shared memory for better performance
    size_t shared_mem_size = 0;  // Can be tuned based on GPU
    
    extract_minimizers_kernel<<<num_blocks, block_size, shared_mem_size, stream1>>>(
        d_sequences, d_sequence_offsets, d_sequence_lengths,
        d_minimizers, d_minimizer_counts,
        num_reads, k, m
    );
    
    // Check for kernel errors
    cudaError_t kernel_error = cudaGetLastError();
    if (kernel_error != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel launch failed: ") + cudaGetErrorString(kernel_error));
    }
    
    // DEBUG: Check kernel execution
    cudaStreamSynchronize(stream1);
    cudaError_t kernel_sync_error = cudaGetLastError();
    if (kernel_sync_error != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel sync failed: ") + cudaGetErrorString(kernel_sync_error));
    }
    
    // Async copy results back
    cudaMemcpyAsync(h_minimizers_pinned, d_minimizers, 
                    num_reads * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer), 
                    cudaMemcpyDeviceToHost, stream2);
    cudaMemcpyAsync(h_counts_pinned, d_minimizer_counts, 
                    num_reads * sizeof(uint32_t), 
                    cudaMemcpyDeviceToHost, stream2);
    
    // Wait for all operations to complete
    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    
    // Check for execution errors
    cudaError_t sync_error = cudaGetLastError();
    if (sync_error != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA kernel execution failed: ") + cudaGetErrorString(sync_error));
    }
    
    // Organize results from pinned memory
    std::vector<std::vector<Minimizer>> results(num_reads);
    for (size_t i = 0; i < num_reads; i++) {
        results[i].reserve(h_counts_pinned[i]);
        for (uint32_t j = 0; j < h_counts_pinned[i]; j++) {
            results[i].push_back(h_minimizers_pinned[i * MAX_MINIMIZERS_PER_READ + j]);
        }
    }
    
    return results;
}