// gpu_minimizer_extraction.cu
// Optimized GPU kernel for Kraken2-style minimizer extraction
// Implements sliding window minimizer algorithm with proper deduplication

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <cub/cub.cuh>

// Configuration constants
#define MAX_KMER_LENGTH 64
#define MAX_THREADS_PER_BLOCK 256
#define SHARED_MEMORY_SIZE 48 * 1024  // 48KB shared memory per block

// Minimizer hit structure (same as your original)
struct GPUMinimizerHit {
    uint64_t minimizer_hash;
    uint32_t taxon_id;
    uint32_t position;
    uint32_t genome_id;
};

// Genome info structure
struct GPUGenomeInfo {
    uint32_t taxon_id;
    uint32_t sequence_offset;
    uint32_t sequence_length;
    uint32_t genome_id;
};

// Parameters for minimizer extraction
struct MinimizerParams {
    int k = 35;          // k-mer length
    int ell = 31;        // minimizer length
    int spaces = 7;      // spaced seed spacing
    uint64_t xor_mask = 0x123456789ABCDEFULL;  // XOR shuffling constant
};

// Device functions for minimizer computation
__device__ uint64_t encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // Invalid base marker
    }
}

__device__ uint64_t hash_lmer(const char* sequence, int pos, int ell) {
    uint64_t hash = 0;
    for (int i = 0; i < ell; i++) {
        uint64_t base = encode_base(sequence[pos + i]);
        if (base == 4) return UINT64_MAX;  // Invalid sequence
        hash = (hash << 2) | base;
    }
    return hash;
}

__device__ uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces, int ell) {
    if (spaces == 0) return hash;
    
    uint64_t masked_hash = 0;
    int out_pos = 0;
    
    // Apply spaced seed pattern: keep every (spaces+1)th position
    for (int i = 0; i < ell; i++) {
        if (i % (spaces + 1) == 0) {
            uint64_t base = (hash >> (2 * (ell - 1 - i))) & 3;
            masked_hash = (masked_hash << 2) | base;
            out_pos++;
        }
    }
    
    return masked_hash;
}

__device__ uint64_t compute_canonical_minimizer(uint64_t hash, uint64_t xor_mask) {
    // Apply XOR shuffling to avoid bias toward low-complexity sequences
    return hash ^ xor_mask;
}

__device__ bool is_valid_sequence(const char* seq, int len) {
    for (int i = 0; i < len; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

// Sliding window minimizer extraction for a single k-mer
__device__ uint64_t extract_minimizer_sliding_window(
    const char* sequence, 
    int kmer_pos, 
    int k, 
    int ell, 
    int spaces,
    uint64_t xor_mask) {
    
    uint64_t min_hash = UINT64_MAX;
    
    // Slide window of size ell across the k-mer
    for (int i = 0; i <= k - ell; i++) {
        if (!is_valid_sequence(sequence + kmer_pos + i, ell)) {
            continue;
        }
        
        uint64_t lmer_hash = hash_lmer(sequence, kmer_pos + i, ell);
        if (lmer_hash == UINT64_MAX) continue;
        
        // Apply spaced seed mask
        uint64_t masked_hash = apply_spaced_seed_mask(lmer_hash, spaces, ell);
        
        // Apply XOR shuffling
        uint64_t canonical_hash = compute_canonical_minimizer(masked_hash, xor_mask);
        
        if (canonical_hash < min_hash) {
            min_hash = canonical_hash;
        }
    }
    
    return min_hash;
}

// Optimized kernel: one block per genome, threads cooperate within block
__global__ void extract_minimizers_optimized_kernel(
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
    
    // Shared memory for collecting minimizers within block
    extern __shared__ uint64_t shared_minimizers[];
    uint32_t* shared_positions = (uint32_t*)(shared_minimizers + block_size);
    
    uint32_t local_minimizer_count = 0;
    uint32_t total_kmers = seq_length - params.k + 1;
    
    // Each thread processes a subset of k-mers using grid-stride
    for (uint32_t kmer_idx = thread_id; kmer_idx < total_kmers; kmer_idx += block_size) {
        
        // Extract minimizer for this k-mer
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, kmer_idx, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (minimizer != UINT64_MAX) {
            // Store in shared memory
            shared_minimizers[thread_id] = minimizer;
            shared_positions[thread_id] = kmer_idx;
            local_minimizer_count++;
        }
        
        __syncthreads();
        
        // Thread 0 collects and deduplicates minimizers from this iteration
        if (thread_id == 0) {
            // Simple deduplication within this block iteration
            for (int i = 0; i < block_size; i++) {
                if (i < total_kmers - (kmer_idx - thread_id) && shared_minimizers[i] != UINT64_MAX) {
                    
                    // Check for duplicates in this batch (simple O(n²) for small batches)
                    bool is_duplicate = false;
                    for (int j = 0; j < i; j++) {
                        if (shared_minimizers[j] == shared_minimizers[i]) {
                            is_duplicate = true;
                            break;
                        }
                    }
                    
                    if (!is_duplicate) {
                        // Get global position and write to output
                        uint32_t global_pos = atomicAdd(global_hit_counter, 1);
                        
                        if (global_pos < max_minimizers) {
                            GPUMinimizerHit hit;
                            hit.minimizer_hash = shared_minimizers[i];
                            hit.taxon_id = genome.taxon_id;
                            hit.position = shared_positions[i];
                            hit.genome_id = genome.genome_id;
                            
                            minimizer_hits[global_pos] = hit;
                        }
                    }
                }
            }
        }
        
        __syncthreads();
        
        // Clear shared memory for next iteration
        if (thread_id < block_size) {
            shared_minimizers[thread_id] = UINT64_MAX;
        }
        
        __syncthreads();
    }
    
    // Store total count for this genome (approximate)
    if (thread_id == 0) {
        hit_counts_per_genome[genome_id] = local_minimizer_count;
    }
}

// Alternative kernel: more sophisticated sliding window with deque-like behavior
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
    
    // Simplified sliding window minimizer computation
    // This processes one genome per thread, suitable for small genomes
    
    uint64_t prev_minimizer = UINT64_MAX;
    uint32_t total_kmers = seq_length - params.k + 1;
    
    for (uint32_t pos = 0; pos < total_kmers; pos++) {
        uint64_t current_minimizer = extract_minimizer_sliding_window(
            sequence, pos, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        // Only output if minimizer changed (simple deduplication)
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

// Host function to launch optimized minimizer extraction
extern "C" bool extract_minimizers_gpu_optimized(
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
    
    // Choose kernel based on number of genomes
    if (num_genomes <= 1024) {
        // Use one-block-per-genome approach for better load balancing
        int threads_per_block = min(MAX_THREADS_PER_BLOCK, 256);
        size_t shared_memory_size = 2 * threads_per_block * sizeof(uint64_t);
        
        extract_minimizers_optimized_kernel<<<num_genomes, threads_per_block, shared_memory_size>>>(
            d_sequence_data, d_genome_info, num_genomes,
            d_minimizer_hits, d_hit_counts, d_global_counter,
            params, max_minimizers
        );
    } else {
        // Use one-thread-per-genome for very large genome counts
        int threads_per_block = 256;
        int num_blocks = (num_genomes + threads_per_block - 1) / threads_per_block;
        
        extract_minimizers_sliding_window_kernel<<<num_blocks, threads_per_block>>>(
            d_sequence_data, d_genome_info, num_genomes,
            d_minimizer_hits, d_global_counter,
            params, max_minimizers
        );
    }
    
    cudaDeviceSynchronize();
    
    // Get total hit count
    cudaMemcpy(total_hits, d_global_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    cudaFree(d_global_counter);
    
    return true;
}

// Host function for post-processing deduplication using Thrust
extern "C" bool deduplicate_minimizers_gpu(
    GPUMinimizerHit* d_minimizer_hits,
    uint32_t num_hits,
    uint32_t* final_count) {
    
    if (num_hits == 0) {
        *final_count = 0;
        return true;
    }
    
    // Sort by minimizer hash
    thrust::device_ptr<GPUMinimizerHit> hits_ptr(d_minimizer_hits);
    thrust::sort(hits_ptr, hits_ptr + num_hits, 
        [] __device__ (const GPUMinimizerHit& a, const GPUMinimizerHit& b) {
            return a.minimizer_hash < b.minimizer_hash;
        });
    
    // Remove duplicates (keep first occurrence of each hash)
    auto new_end = thrust::unique(hits_ptr, hits_ptr + num_hits,
        [] __device__ (const GPUMinimizerHit& a, const GPUMinimizerHit& b) {
            return a.minimizer_hash == b.minimizer_hash;
        });
    
    *final_count = new_end - hits_ptr;
    
    return true;
}

// Test function to validate minimizer extraction
extern "C" void test_minimizer_extraction() {
    // Test with small example
    const char* test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    int seq_len = strlen(test_sequence);
    
    MinimizerParams params;
    params.k = 10;
    params.ell = 8;
    params.spaces = 2;
    
    printf("Testing minimizer extraction:\n");
    printf("Sequence: %s\n", test_sequence);
    printf("Parameters: k=%d, ell=%d, spaces=%d\n", params.k, params.ell, params.spaces);
    
    // Extract minimizers for each k-mer manually
    for (int pos = 0; pos <= seq_len - params.k; pos++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            test_sequence, pos, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        printf("K-mer at pos %d: minimizer = 0x%016lx\n", pos, minimizer);
    }
}