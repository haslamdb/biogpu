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
                    
                    // Look at the last m k-mers in the circular buffer
                    int actual_window_size = min(window_size, m);
                    for (int i = 0; i < actual_window_size; i++) {
                        if (window_kmers[i] < min_hash) {
                            min_hash = window_kmers[i];
                            min_idx = i;
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
    : k(k_size), m(window_size), allocated_reads(0), 
      d_sequences(nullptr), d_sequence_offsets(nullptr),
      d_sequence_lengths(nullptr), d_minimizers(nullptr),
      d_minimizer_counts(nullptr) {
    init_nucleotide_map();
}

MinimizerExtractor::~MinimizerExtractor() {
    if (d_sequences) cudaFree(d_sequences);
    if (d_sequence_offsets) cudaFree(d_sequence_offsets);
    if (d_sequence_lengths) cudaFree(d_sequence_lengths);
    if (d_minimizers) cudaFree(d_minimizers);
    if (d_minimizer_counts) cudaFree(d_minimizer_counts);
}

void MinimizerExtractor::allocate_device_memory(size_t num_reads, size_t total_sequence_length) {
    if (num_reads > allocated_reads) {
        if (d_sequences) cudaFree(d_sequences);
        if (d_sequence_offsets) cudaFree(d_sequence_offsets);
        if (d_sequence_lengths) cudaFree(d_sequence_lengths);
        if (d_minimizers) cudaFree(d_minimizers);
        if (d_minimizer_counts) cudaFree(d_minimizer_counts);
        
        cudaMalloc(&d_sequences, total_sequence_length);
        cudaMalloc(&d_sequence_offsets, (num_reads + 1) * sizeof(uint32_t));
        cudaMalloc(&d_sequence_lengths, num_reads * sizeof(uint32_t));
        cudaMalloc(&d_minimizers, num_reads * MAX_MINIMIZERS_PER_READ * sizeof(Minimizer));
        cudaMalloc(&d_minimizer_counts, num_reads * sizeof(uint32_t));
        
        allocated_reads = num_reads;
    }
}

std::vector<std::vector<Minimizer>> MinimizerExtractor::extract_minimizers(
    const std::vector<std::string>& sequences
) {
    size_t num_reads = sequences.size();
    
    // Prepare data for GPU
    std::vector<char> concatenated_sequences;
    std::vector<uint32_t> offsets(num_reads + 1);
    std::vector<uint32_t> lengths(num_reads);
    
    uint32_t current_offset = 0;
    for (size_t i = 0; i < num_reads; i++) {
        offsets[i] = current_offset;
        lengths[i] = sequences[i].length();
        concatenated_sequences.insert(concatenated_sequences.end(), 
                                    sequences[i].begin(), 
                                    sequences[i].end());
        current_offset += lengths[i];
    }
    offsets[num_reads] = current_offset;
    
    // Allocate device memory
    allocate_device_memory(num_reads, concatenated_sequences.size());
    
    // Copy to device
    cudaMemcpy(d_sequences, concatenated_sequences.data(), 
               concatenated_sequences.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sequence_offsets, offsets.data(), 
               offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sequence_lengths, lengths.data(), 
               lengths.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch kernel
    int block_size = 256;
    int num_blocks = (num_reads + block_size - 1) / block_size;
    
    extract_minimizers_kernel<<<num_blocks, block_size>>>(
        d_sequences, d_sequence_offsets, d_sequence_lengths,
        d_minimizers, d_minimizer_counts,
        num_reads, k, m
    );
    
    cudaDeviceSynchronize();
    
    // Copy results back
    std::vector<Minimizer> all_minimizers(num_reads * MAX_MINIMIZERS_PER_READ);
    std::vector<uint32_t> minimizer_counts(num_reads);
    
    cudaMemcpy(all_minimizers.data(), d_minimizers, 
               all_minimizers.size() * sizeof(Minimizer), cudaMemcpyDeviceToHost);
    cudaMemcpy(minimizer_counts.data(), d_minimizer_counts, 
               minimizer_counts.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    // Organize results
    std::vector<std::vector<Minimizer>> results(num_reads);
    for (size_t i = 0; i < num_reads; i++) {
        results[i].reserve(minimizer_counts[i]);
        for (uint32_t j = 0; j < minimizer_counts[i]; j++) {
            results[i].push_back(all_minimizers[i * MAX_MINIMIZERS_PER_READ + j]);
        }
    }
    
    return results;
}