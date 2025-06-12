// minimizer_extraction.cu - Fixed version
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
    cudaMemcpyToSymbol(nucleotide_map, h_map, 256);
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
    
    // Local storage for window k-mers (circular buffer)
    uint64_t window_kmers[32];  // Assuming m <= 32
    uint32_t window_positions[32];
    bool window_reverse[32];
    
    // Variables for rolling hash
    uint64_t forward_kmer = 0;
    uint64_t reverse_kmer = 0;
    const uint64_t kmer_mask = (1ULL << (2 * k)) - 1;
    
    // Initialize first k-1 bases
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
    uint64_t prev_minimizer_hash = UINT64_MAX;
    uint32_t prev_minimizer_pos = UINT32_MAX;
    int window_start = 0;  // Start of valid window
    int window_end = 0;    // End of valid window (exclusive)
    
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
                uint32_t kmer_pos = pos - k + 1;
                
                // Add to circular buffer
                int buffer_idx = window_end % m;
                window_kmers[buffer_idx] = kmer_hash;
                window_positions[buffer_idx] = kmer_pos;
                window_reverse[buffer_idx] = (reverse_kmer < forward_kmer);
                window_end++;
                
                // Maintain window size
                if (window_end - window_start > m) {
                    window_start++;
                }
                
                // Find minimizer in current window
                if (window_end - window_start >= m || pos == seq_length - 1) {
                    uint64_t min_hash = UINT64_MAX;
                    int min_idx = -1;
                    uint32_t min_pos = 0;
                    bool min_reverse = false;
                    
                    // Search through valid window
                    for (int i = window_start; i < window_end; i++) {
                        int idx = i % m;
                        if (window_kmers[idx] < min_hash) {
                            min_hash = window_kmers[idx];
                            min_idx = idx;
                            min_pos = window_positions[idx];
                            min_reverse = window_reverse[idx];
                        }
                    }
                    
                    // Store minimizer if it's different from the previous one
                    if (min_idx != -1 && (min_hash != prev_minimizer_hash || min_pos != prev_minimizer_pos)) {
                        if (minimizer_count < MAX_MINIMIZERS_PER_READ) {
                            Minimizer* out_minimizer = &minimizers[tid * MAX_MINIMIZERS_PER_READ + minimizer_count];
                            out_minimizer->hash = min_hash;
                            out_minimizer->position = min_pos;
                            out_minimizer->is_reverse = min_reverse;
                            minimizer_count++;
                            
                            prev_minimizer_hash = min_hash;
                            prev_minimizer_pos = min_pos;
                        }
                    }
                }
            }
        } else {
            // Invalid base, reset everything
            forward_kmer = 0;
            reverse_kmer = 0;
            valid_bases = 0;
            window_start = window_end;  // Reset window
            prev_minimizer_hash = UINT64_MAX;
            prev_minimizer_pos = UINT32_MAX;
        }
    }
    
    minimizer_counts[tid] = minimizer_count;
}

// Host wrapper class implementation
MinimizerExtractor::MinimizerExtractor(int k_size, int window_size) 
    : k(k_size), m(window_size), allocated_reads(0), allocated_seq_length(0),
      allocated_pinned_reads(0), allocated_pinned_seq_length(0),
      d_sequences(nullptr), d_sequence_offsets(nullptr),
      d_sequence_lengths(nullptr), d_minimizers(nullptr),
      d_minimizer_counts(nullptr),
      h_sequences(nullptr), h_offsets(nullptr), h_lengths(nullptr),
      h_minimizers(nullptr), h_counts(nullptr) {
    
    // Initialize the nucleotide map once
    static bool map_initialized = false;
    if (!map_initialized) {
        init_nucleotide_map();
        map_initialized = true;
    }
    
    // Pre-allocate pinned memory for typical batch size
    size_t typical_batch = 10000;
    size_t typical_seq_length = typical_batch * 200;  // Assume 200bp average
    
    cudaMallocHost(&h_sequences, typical_seq_length);
    cudaMallocHost(&h_offsets, (typical_batch + 1) * sizeof(uint32_t));
    cudaMallocHost(&h_lengths, typical_batch * sizeof(uint32_t));
    cudaMallocHost(&h_minimizers, typical_batch * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
    cudaMallocHost(&h_counts, typical_batch * sizeof(uint32_t));
    
    allocated_pinned_reads = typical_batch;
    allocated_pinned_seq_length = typical_seq_length;
}

MinimizerExtractor::~MinimizerExtractor() {
    if (d_sequences) cudaFree(d_sequences);
    if (d_sequence_offsets) cudaFree(d_sequence_offsets);
    if (d_sequence_lengths) cudaFree(d_sequence_lengths);
    if (d_minimizers) cudaFree(d_minimizers);
    if (d_minimizer_counts) cudaFree(d_minimizer_counts);
    
    if (h_sequences) cudaFreeHost(h_sequences);
    if (h_offsets) cudaFreeHost(h_offsets);
    if (h_lengths) cudaFreeHost(h_lengths);
    if (h_minimizers) cudaFreeHost(h_minimizers);
    if (h_counts) cudaFreeHost(h_counts);
}

void MinimizerExtractor::allocate_device_memory(size_t num_reads, size_t total_sequence_length) {
    // Allocate device memory if needed
    if (num_reads > allocated_reads) {
        if (d_sequence_offsets) cudaFree(d_sequence_offsets);
        if (d_sequence_lengths) cudaFree(d_sequence_lengths);
        if (d_minimizers) cudaFree(d_minimizers);
        if (d_minimizer_counts) cudaFree(d_minimizer_counts);
        
        cudaMalloc(&d_sequence_offsets, (num_reads + 1) * sizeof(uint32_t));
        cudaMalloc(&d_sequence_lengths, num_reads * sizeof(uint32_t));
        cudaMalloc(&d_minimizers, num_reads * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
        cudaMalloc(&d_minimizer_counts, num_reads * sizeof(uint32_t));
        
        allocated_reads = num_reads;
    }
    
    if (total_sequence_length > allocated_seq_length) {
        if (d_sequences) cudaFree(d_sequences);
        cudaMalloc(&d_sequences, total_sequence_length);
        allocated_seq_length = total_sequence_length;
    }
}

std::vector<std::vector<Minimizer>> MinimizerExtractor::extract_minimizers(
    const std::vector<std::string>& sequences
) {
    size_t num_reads = sequences.size();
    if (num_reads == 0) return {};
    
    // Calculate total sequence length first
    size_t total_length = 0;
    for (size_t i = 0; i < num_reads; i++) {
        total_length += sequences[i].length();
    }
    
    // Reallocate pinned memory if needed
    if (num_reads > allocated_pinned_reads || total_length > allocated_pinned_seq_length) {
        // Free old pinned memory
        if (h_sequences) cudaFreeHost(h_sequences);
        if (h_offsets) cudaFreeHost(h_offsets);
        if (h_lengths) cudaFreeHost(h_lengths);
        if (h_minimizers) cudaFreeHost(h_minimizers);
        if (h_counts) cudaFreeHost(h_counts);
        
        // Allocate new pinned memory with some extra space
        size_t new_reads_capacity = std::max((size_t)(num_reads * 1.5), (size_t)10000);
        size_t new_seq_capacity = std::max((size_t)(total_length * 1.5), (size_t)(new_reads_capacity * 200));
        
        cudaMallocHost(&h_sequences, new_seq_capacity);
        cudaMallocHost(&h_offsets, (new_reads_capacity + 1) * sizeof(uint32_t));
        cudaMallocHost(&h_lengths, new_reads_capacity * sizeof(uint32_t));
        cudaMallocHost(&h_minimizers, new_reads_capacity * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
        cudaMallocHost(&h_counts, new_reads_capacity * sizeof(uint32_t));
        
        allocated_pinned_reads = new_reads_capacity;
        allocated_pinned_seq_length = new_seq_capacity;
    }
    
    // Prepare data for GPU
    total_length = 0;
    for (size_t i = 0; i < num_reads; i++) {
        h_offsets[i] = total_length;
        h_lengths[i] = sequences[i].length();
        
        // Copy sequence data
        memcpy(h_sequences + total_length, sequences[i].c_str(), sequences[i].length());
        total_length += sequences[i].length();
    }
    h_offsets[num_reads] = total_length;
    
    // Allocate device memory
    allocate_device_memory(num_reads, total_length);
    
    // Create CUDA stream for async operations
    cudaStream_t stream;
    cudaStreamCreate(&stream);
    
    // Copy to device asynchronously
    cudaMemcpyAsync(d_sequences, h_sequences, total_length, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_sequence_offsets, h_offsets, (num_reads + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_sequence_lengths, h_lengths, num_reads * sizeof(uint32_t), cudaMemcpyHostToDevice, stream);
    
    // Launch kernel
    int block_size = 256;
    int num_blocks = (num_reads + block_size - 1) / block_size;
    
    extract_minimizers_kernel<<<num_blocks, block_size, 0, stream>>>(
        d_sequences, d_sequence_offsets, d_sequence_lengths,
        d_minimizers, d_minimizer_counts,
        num_reads, k, m
    );
    
    // Copy results back asynchronously
    cudaMemcpyAsync(h_minimizers, d_minimizers, 
                    num_reads * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer), 
                    cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(h_counts, d_minimizer_counts, 
                    num_reads * sizeof(uint32_t), 
                    cudaMemcpyDeviceToHost, stream);
    
    // Wait for all operations to complete
    cudaStreamSynchronize(stream);
    cudaStreamDestroy(stream);
    
    // Check for errors
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(error) << std::endl;
        throw std::runtime_error("CUDA kernel execution failed");
    }
    
    // Organize results
    std::vector<std::vector<Minimizer>> results(num_reads);
    for (size_t i = 0; i < num_reads; i++) {
        results[i].reserve(h_counts[i]);
        for (uint32_t j = 0; j < h_counts[i]; j++) {
            results[i].push_back(h_minimizers[i * MAX_MINIMIZERS_PER_READ + j]);
        }
    }
    
    return results;
}