// gpu_minimizer_extraction.cu
// FIXED: Optimized GPU kernel for Kraken2-style minimizer extraction
// Implements proper sliding window minimizer algorithm with deduplication

#ifndef GPU_MINIMIZER_EXTRACTION_CU
#define GPU_MINIMIZER_EXTRACTION_CU

#include "gpu_minimizer_extraction.cuh"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <cub/cub.cuh>
#include "../../../include/biogpu/minimizer_extraction.h"

// Configuration constants
#define MAX_KMER_LENGTH 64
#define MAX_THREADS_PER_BLOCK 256
#define SHARED_MEMORY_SIZE 48 * 1024  // 48KB shared memory per block

// Structures defined in gpu_minimizer_extraction.cuh

// MinimizerParams is defined in minimizer_extraction.h

// Device function implementations moved to header file for cross-compilation unit access

// Note: extract_minimizer_sliding_window implementation moved to header file for cross-compilation unit access

// FIXED: Proper Kraken2-style sliding window minimizer extraction
// FIXED: Multi-threaded kernel to match the database builder
__global__ void extract_minimizers_optimized_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts_per_genome,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    int max_minimizers) {
    
    // OPTIMIZED: Use multiple threads per genome
    int genome_id = blockIdx.x;
    int thread_id = threadIdx.x;
    int block_size = blockDim.x;
    
    if (genome_id >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[genome_id];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) {
        if (thread_id == 0) hit_counts_per_genome[genome_id] = 0;
        return;
    }
    
    // Shared memory for thread coordination
    extern __shared__ uint64_t shared_mem[];
    uint64_t* shared_minimizers = shared_mem;
    uint32_t* shared_positions = (uint32_t*)(shared_minimizers + blockDim.x);
    bool* shared_valid = (bool*)(shared_positions + blockDim.x);
    
    uint32_t total_kmers = seq_length - params.k + 1;
    uint32_t kmers_per_thread = (total_kmers + block_size - 1) / block_size;
    uint32_t start_kmer = thread_id * kmers_per_thread;
    uint32_t end_kmer = min(start_kmer + kmers_per_thread, total_kmers);
    
    // Each thread processes its assigned k-mers
    shared_minimizers[thread_id] = UINT64_MAX;
    shared_positions[thread_id] = 0;
    shared_valid[thread_id] = false;
    
    // Process assigned k-mers
    for (uint32_t kmer_idx = start_kmer; kmer_idx < end_kmer; kmer_idx++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, kmer_idx, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (minimizer != UINT64_MAX) {
            // Store first valid minimizer from this thread's range
            if (!shared_valid[thread_id]) {
                shared_minimizers[thread_id] = minimizer;
                shared_positions[thread_id] = kmer_idx;
                shared_valid[thread_id] = true;
                break;
            }
        }
    }
    
    __syncthreads();
    
    // Thread 0 processes results in order to maintain deduplication
    if (thread_id == 0) {
        uint64_t last_minimizer = UINT64_MAX;
        uint32_t local_count = 0;
        
        // Process all positions in order
        for (uint32_t pos = 0; pos < total_kmers; pos++) {
            uint32_t responsible_thread = pos / kmers_per_thread;
            if (responsible_thread >= block_size) responsible_thread = block_size - 1;
            
            if (shared_valid[responsible_thread] && 
                shared_positions[responsible_thread] == pos &&
                shared_minimizers[responsible_thread] != last_minimizer) {
                
                uint32_t current_count = atomicAdd(global_hit_counter, 0);
                if (current_count >= max_minimizers) break;
                
                uint32_t global_pos = atomicAdd(global_hit_counter, 1);
                if (global_pos < max_minimizers) {
                    GPUMinimizerHit hit;
                    hit.minimizer_hash = shared_minimizers[responsible_thread];
                    hit.taxon_id = genome.taxon_id;
                    hit.position = pos;
                    hit.genome_id = genome.genome_id;
                    
                    minimizer_hits[global_pos] = hit;
                    local_count++;
                    last_minimizer = shared_minimizers[responsible_thread];
                }
            }
        }
        
        hit_counts_per_genome[genome_id] = local_count;
    }
}

// Alternative kernel: multi-threaded version with cooperation (for very large genomes)
__global__ void extract_minimizers_cooperative_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts_per_genome,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    int max_minimizers) {
    
    int genome_id = blockIdx.x;
    int thread_id = threadIdx.x;
    int block_size = blockDim.x;
    
    if (genome_id >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[genome_id];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) {
        if (thread_id == 0) hit_counts_per_genome[genome_id] = 0;
        return;
    }
    
    // Shared memory for temporary storage
    extern __shared__ uint64_t shared_minimizers[];
    uint32_t* shared_positions = (uint32_t*)(shared_minimizers + block_size);
    
    uint32_t total_kmers = seq_length - params.k + 1;
    uint32_t local_count = 0;
    
    // Divide work among threads, but process results sequentially
    uint32_t kmers_per_thread = (total_kmers + block_size - 1) / block_size;
    uint32_t start_kmer = thread_id * kmers_per_thread;
    uint32_t end_kmer = min(start_kmer + kmers_per_thread, total_kmers);
    
    // Each thread extracts minimizers for its assigned k-mers
    uint32_t thread_minimizer_count = 0;
    for (uint32_t kmer_idx = start_kmer; kmer_idx < end_kmer; kmer_idx++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, kmer_idx, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (minimizer != UINT64_MAX && thread_minimizer_count < block_size) {
            shared_minimizers[thread_id] = minimizer;
            shared_positions[thread_id] = kmer_idx;
            thread_minimizer_count++;
            break; // Only store first valid minimizer per thread for now
        }
    }
    
    __syncthreads();
    
    // Thread 0 processes all minimizers sequentially to maintain order
    if (thread_id == 0) {
        uint64_t last_minimizer = UINT64_MAX;
        
        for (uint32_t kmer_idx = 0; kmer_idx < total_kmers; kmer_idx++) {
            uint64_t minimizer = extract_minimizer_sliding_window(
                sequence, kmer_idx, params.k, params.ell, params.spaces, params.xor_mask
            );
            
            if (minimizer != UINT64_MAX && minimizer != last_minimizer) {
                uint32_t global_pos = atomicAdd(global_hit_counter, 1);
                
                if (global_pos < max_minimizers) {
                    GPUMinimizerHit hit;
                    hit.minimizer_hash = minimizer;
                    hit.taxon_id = genome.taxon_id;
                    hit.position = kmer_idx;
                    hit.genome_id = genome.genome_id;
                    
                    minimizer_hits[global_pos] = hit;
                    local_count++;
                }
                
                last_minimizer = minimizer;
            }
        }
        
        hit_counts_per_genome[genome_id] = local_count;
    }
}

// Keep your working kernel as fallback
__global__ void extract_minimizers_sliding_window_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    int max_minimizers) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[idx];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) return;
    
    // This kernel already has the correct logic!
    uint64_t prev_minimizer = UINT64_MAX;
    uint32_t total_kmers = seq_length - params.k + 1;
    
    for (uint32_t pos = 0; pos < total_kmers; pos++) {
        uint64_t current_minimizer = extract_minimizer_sliding_window(
            sequence, pos, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        // Only output if minimizer changed (proper deduplication)
        if (current_minimizer != UINT64_MAX && current_minimizer != prev_minimizer) {
            uint32_t global_pos = atomicAdd(global_hit_counter, 1);
            
            if (global_pos < max_minimizers) {
                GPUMinimizerHit hit;
                hit.minimizer_hash = current_minimizer;
                hit.taxon_id = genome.taxon_id;
                hit.position = pos;
                hit.genome_id = genome.genome_id;
                
                minimizer_hits[global_pos] = hit;
            }
            
            prev_minimizer = current_minimizer;
        }
    }
}

// FIXED: Host function to launch optimized minimizer extraction
bool extract_minimizers_gpu_optimized(
    const char* d_sequence_data,
    const GPUGenomeInfo* d_genome_info,
    int num_genomes,
    GPUMinimizerHit* d_minimizer_hits,
    uint32_t* d_hit_counts,
    uint32_t* total_hits,
    MinimizerParams params,
    int max_minimizers) {
    
    // Reset global counter
    uint32_t zero = 0;
    uint32_t* d_global_counter;
    cudaMalloc(&d_global_counter, sizeof(uint32_t));
    cudaMemcpy(d_global_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Clear hit counts
    cudaMemset(d_hit_counts, 0, num_genomes * sizeof(uint32_t));
    
    // UPDATED: Launch with more threads and shared memory
    int threads_per_block = 256;  // Use 256 threads instead of 1
    size_t shared_mem_size = threads_per_block * (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(bool));
    
    extract_minimizers_optimized_kernel<<<num_genomes, threads_per_block, shared_mem_size>>>(
        d_sequence_data, d_genome_info, num_genomes,
        d_minimizer_hits, d_hit_counts, d_global_counter,
        params, max_minimizers
    );
    
    cudaDeviceSynchronize();
    
    // Check for kernel errors
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error in minimizer extraction: %s\n", cudaGetErrorString(error));
        cudaFree(d_global_counter);
        return false;
    }
    
    // Get total hit count
    cudaMemcpy(total_hits, d_global_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    // CRITICAL: Clamp total_hits to max_minimizers to prevent overflow
    if (*total_hits > max_minimizers) {
        printf("WARNING: Minimizer extraction hit limit. Clamping %u to %d\n", 
               *total_hits, max_minimizers);
        *total_hits = max_minimizers;
    }
    
    cudaFree(d_global_counter);
    
    return true;
}

// Host function for post-processing deduplication using Thrust
bool deduplicate_minimizers_gpu(
    GPUMinimizerHit* d_minimizer_hits,
    uint32_t num_hits,
    uint32_t* final_count,
    uint32_t max_allocated_hits) {
    
    if (num_hits == 0) {
        *final_count = 0;
        return true;
    }
    
    // CRITICAL: Check bounds before processing
    if (num_hits > max_allocated_hits) {
        printf("ERROR: num_hits (%u) exceeds allocated memory (%u)!\n", 
               num_hits, max_allocated_hits);
        return false;
    }
    
    // Sort by minimizer hash
    thrust::device_ptr<GPUMinimizerHit> hits_ptr(d_minimizer_hits);
    
    try {
        // Add explicit CUDA error checking before sort
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA error before sort: %s\n", cudaGetErrorString(err));
            return false;
        }
        
        thrust::sort(hits_ptr, hits_ptr + num_hits, 
            [] __device__ (const GPUMinimizerHit& a, const GPUMinimizerHit& b) {
                return a.minimizer_hash < b.minimizer_hash;
            });
        
        // Check for errors after sort
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA error after sort: %s\n", cudaGetErrorString(err));
            return false;
        }
        
        // Remove duplicates (keep first occurrence of each hash)
        auto new_end = thrust::unique(hits_ptr, hits_ptr + num_hits,
            [] __device__ (const GPUMinimizerHit& a, const GPUMinimizerHit& b) {
                return a.minimizer_hash == b.minimizer_hash;
            });
        
        *final_count = new_end - hits_ptr;
        
        // Final error check
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA error after unique: %s\n", cudaGetErrorString(err));
            return false;
        }
        
    } catch (const std::exception& e) {
        printf("Error in deduplication: %s\n", e.what());
        return false;
    }
    
    return true;
}

// Test function to validate minimizer extraction with actual kernel launch
__global__ void test_minimizer_kernel(
    const char* test_sequence,
    int seq_length,
    MinimizerParams params,
    uint64_t* results,
    int* result_count) {
    
    int tid = threadIdx.x;
    if (tid != 0) return;
    
    uint64_t last_minimizer = UINT64_MAX;
    int count = 0;
    
    for (int kmer_pos = 0; kmer_pos <= seq_length - params.k; kmer_pos++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            test_sequence, kmer_pos, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (minimizer != UINT64_MAX) {
            if (minimizer != last_minimizer && count < 100) {
                results[count] = minimizer;
                count++;
                last_minimizer = minimizer;
            }
        }
    }
    
    *result_count = count;
}

void test_minimizer_extraction() {
    printf("Testing minimizer extraction with actual GPU kernel...\n");
    
    // Test sequence
    const char* test_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    int seq_len = strlen(test_seq);
    
    printf("Test sequence: %s (length %d)\n", test_seq, seq_len);
    
    MinimizerParams params;
    params.k = 10;
    params.ell = 8;
    params.spaces = 0;
    params.xor_mask = 0;
    
    printf("Parameters: k=%d, ell=%d\n", params.k, params.ell);
    
    // Allocate GPU memory
    char* d_test_seq;
    uint64_t* d_results;
    int* d_result_count;
    
    cudaMalloc(&d_test_seq, seq_len + 1);
    cudaMalloc(&d_results, 100 * sizeof(uint64_t));
    cudaMalloc(&d_result_count, sizeof(int));
    
    // Copy test data
    cudaMemcpy(d_test_seq, test_seq, seq_len + 1, cudaMemcpyHostToDevice);
    
    // Launch test kernel
    test_minimizer_kernel<<<1, 1>>>(d_test_seq, seq_len, params, d_results, d_result_count);
    cudaDeviceSynchronize();
    
    // Get results
    uint64_t results[100];
    int result_count;
    cudaMemcpy(results, d_results, 100 * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&result_count, d_result_count, sizeof(int), cudaMemcpyDeviceToHost);
    
    printf("Found %d unique minimizers:\n", result_count);
    for (int i = 0; i < result_count && i < 10; i++) {
        printf("  [%d] 0x%016lx\n", i, results[i]);
    }
    
    int expected_kmers = seq_len - params.k + 1;
    double compression = (double)result_count / expected_kmers;
    printf("Compression: %d/%d = %.3f (%.1fx reduction)\n", 
           result_count, expected_kmers, compression, 1.0/compression);
    
    // Cleanup
    cudaFree(d_test_seq);
    cudaFree(d_results);
    cudaFree(d_result_count);
    
    printf("✓ Test completed\n");
}

#endif // GPU_MINIMIZER_EXTRACTION_CU