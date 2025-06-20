// gpu/gpu_database_kernels.cu
// CUDA kernel implementations for database building operations
// Extracted from monolithic gpu_kraken_database_builder.cu

#include "gpu_database_kernels.h"
#include "../gpu_kraken_types.h"
#include "gpu_minimizer_extraction.cuh"
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <vector>
#include <cstdio>
#include <algorithm>

// CUDA_CHECK_KERNEL macro is already defined in gpu_database_kernels.h

// ===========================
// Device Function Implementations
// ===========================

__device__ uint32_t jenkins_hash_gpu(uint64_t key) {
    uint32_t hash = (uint32_t)(key ^ (key >> 32));
    hash += (hash << 10);
    hash ^= (hash >> 6);
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

__device__ uint32_t compute_compact_hash_gpu(uint64_t minimizer_hash) {
    return jenkins_hash_gpu(minimizer_hash) & 0x7FFFFFFF;
}

__device__ uint32_t lookup_lca_gpu(const GPUCompactHashTable* cht, uint64_t minimizer_hash) {
    uint32_t compact_hash = compute_compact_hash_gpu(minimizer_hash);
    uint32_t pos = compact_hash & cht->hash_mask;
    uint32_t lca_mask = (1U << cht->lca_bits) - 1;
    
    for (int probe = 0; probe < 32; probe++) {
        uint32_t cell = cht->hash_cells[pos];
        if (cell == 0) return 0;
        
        uint32_t stored_hash = cell >> cht->lca_bits;
        uint32_t expected_hash = compact_hash >> cht->lca_bits;
        
        if (stored_hash == expected_hash) {
            return cell & lca_mask;
        }
        
        pos = (pos + 1) & cht->hash_mask;
    }
    
    return 0;
}

__device__ bool has_valid_bases_device(const char* seq, int len) {
    for (int i = 0; i < len; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

__device__ uint32_t find_lca_simple_device(uint32_t taxon1, uint32_t taxon2,
                                          const uint32_t* parent_lookup,
                                          uint32_t max_taxon_id) {
    if (taxon1 == taxon2) return taxon1;
    if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return 1;
    if (taxon1 == 0 || taxon2 == 0) return 1;
    
    // Simple LCA implementation
    uint32_t current1 = taxon1;
    uint32_t current2 = taxon2;
    
    for (int steps = 0; steps < 50; steps++) {
        if (current1 == current2) return current1;
        if (current1 == 0 || current1 == 1) break;
        if (current2 == 0 || current2 == 1) break;
        
        if (current1 > current2) {
            current1 = parent_lookup[current1];
        } else {
            current2 = parent_lookup[current2];
        }
    }
    
    return 1;  // Root fallback
}

__device__ uint32_t compute_simple_lca_gpu(uint32_t taxon1, uint32_t taxon2) {
    if (taxon1 == 0) return taxon2;
    if (taxon2 == 0) return taxon1;
    if (taxon1 == taxon2) return taxon1;
    
    // Simplified LCA - return smaller taxon ID
    return (taxon1 < taxon2) ? taxon1 : taxon2;
}

__device__ bool validate_sequence_device(const char* sequence, int length) {
    if (length <= 0) return false;
    return has_valid_bases_device(sequence, length);
}

__device__ void atomic_add_safe_device(uint32_t* address, uint32_t value) {
    atomicAdd(address, value);
}

// ===========================
// CUDA Kernel Implementations
// ===========================

__global__ void extract_minimizers_kraken2_improved_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts_per_genome,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    uint64_t min_clear_hash_value,
    uint64_t toggle_mask,
    int max_minimizers) {
    
    // Use multiple threads per genome for better performance
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
            // Apply subsampling if enabled
            if (min_clear_hash_value > 0 && minimizer < min_clear_hash_value) {
                continue;
            }
            
            // Apply toggle mask
            minimizer ^= toggle_mask;
            
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
                
                // Atomic add first, then check bounds to prevent overflow
                uint32_t global_pos = atomicAdd(global_hit_counter, 1);
                if (global_pos < max_minimizers) {
                    GPUMinimizerHit hit;
                    hit.minimizer_hash = shared_minimizers[responsible_thread];
                    // hit doesn't have taxon_id - stored in genome_info
                    hit.position = pos;
                    hit.genome_id = genome.genome_id;
                    
                    minimizer_hits[global_pos] = hit;
                    local_count++;
                    last_minimizer = shared_minimizers[responsible_thread];
                } else {
                    // Hit the limit, stop processing
                    break;
                }
            }
        }
        
        hit_counts_per_genome[genome_id] = local_count;
    }
}

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
                // hit doesn't have taxon_id - stored in genome_info
                hit.position = pos;
                hit.genome_id = genome.genome_id;
                
                minimizer_hits[global_pos] = hit;
            }
            
            prev_minimizer = current_minimizer;
        }
    }
}

__global__ void convert_hits_to_candidates_kernel(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_hits) return;
    
    const GPUMinimizerHit& hit = hits[idx];
    LCACandidate& candidate = candidates[idx];
    
    candidate.minimizer_hash = hit.minimizer_hash;
    candidate.lca_taxon = 0;  // TODO: Need to look up taxon_id from genome_info
    candidate.genome_count = 1;  // Will be updated during final merge
    candidate.uniqueness_score = 1.0f;
}

__global__ void merge_duplicate_minimizers_kernel(
    LCACandidate* candidates,
    int num_candidates,
    int* final_count) {
    
    // This is a simplified merge - in practice you'd use more sophisticated algorithms
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_candidates) return;
    
    // For now, just pass through - merging would be done on host with Thrust
    if (idx == 0) {
        *final_count = num_candidates;
    }
}

__global__ void compute_lca_for_minimizers_kernel(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates,
    const uint32_t* parent_lookup,
    uint32_t max_taxon_id) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_hits) return;
    
    const GPUMinimizerHit& hit = hits[idx];
    LCACandidate& candidate = candidates[idx];
    
    candidate.minimizer_hash = hit.minimizer_hash;
    candidate.lca_taxon = 0;  // TODO: Need to look up taxon_id from genome_info
    candidate.genome_count = 1;
    candidate.uniqueness_score = 1.0f;
    
    // For multiple hits with same minimizer, compute LCA
    // This is simplified - real implementation would group by minimizer first
}

__global__ void initialize_gpu_memory_kernel(
    char* sequence_buffer,
    GPUGenomeInfo* genome_buffer,
    GPUMinimizerHit* hit_buffer,
    size_t sequence_size,
    size_t genome_count,
    size_t hit_count) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Initialize sequence buffer
    if (idx < sequence_size) {
        sequence_buffer[idx] = 0;
    }
    
    // Initialize genome buffer
    if (idx < genome_count) {
        genome_buffer[idx].taxon_id = 0;
        genome_buffer[idx].sequence_offset = 0;
        genome_buffer[idx].sequence_length = 0;
        genome_buffer[idx].genome_id = 0;
    }
    
    // Initialize hit buffer
    if (idx < hit_count) {
        hit_buffer[idx].minimizer_hash = 0;
        hit_buffer[idx].genome_id = 0;
        hit_buffer[idx].position = 0;
        hit_buffer[idx].genome_id = 0;
    }
}

__global__ void validate_memory_integrity_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    bool* validation_results) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[idx];
    bool is_valid = true;
    
    // Validate genome info
    if (genome.sequence_length == 0 || genome.taxon_id == 0) {
        is_valid = false;
    }
    
    // Validate sequence data accessibility
    if (is_valid && genome.sequence_length > 0) {
        const char* sequence = sequence_data + genome.sequence_offset;
        // Simple validation - check first and last characters
        if (sequence[0] == 0 || sequence[genome.sequence_length - 1] == 0) {
            is_valid = false;
        }
    }
    
    validation_results[idx] = is_valid;
}

// ===========================
// Host Wrapper Functions
// ===========================

LaunchConfig calculate_optimal_launch_config(int num_elements, int threads_per_block, size_t shared_memory) {
    LaunchConfig config;
    
    config.threads_x = threads_per_block;
    config.blocks_x = (num_elements + threads_per_block - 1) / threads_per_block;
    config.shared_memory_bytes = shared_memory;
    config.stream = nullptr;  // Use default stream
    
    // Ensure grid size doesn't exceed limits
    if (config.blocks_x > 65535) {
        config.blocks_x = 65535;
    }
    
    return config;
}

bool launch_minimizer_extraction_kernel(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint32_t* total_hits_output) {
    
    // Calculate launch configuration
    LaunchConfig launch_config = calculate_optimal_launch_config(
        batch_data.max_genomes, 
        256,  // default threads per block
        0     // no shared memory needed
    );
    
    // Launch the kernel
    dim3 grid(launch_config.blocks_x, launch_config.blocks_y, launch_config.blocks_z);
    dim3 block(launch_config.threads_x, launch_config.threads_y, launch_config.threads_z);
    extract_minimizers_sliding_window_kernel<<<
        grid,
        block,
        launch_config.shared_memory_bytes,
        (cudaStream_t)launch_config.stream
    >>>(
        batch_data.d_sequence_data,
        batch_data.d_genome_info,
        batch_data.max_genomes,
        batch_data.d_minimizer_hits,
        batch_data.d_global_counter,
        params,
        batch_data.max_minimizers
    );
    
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
    
    // Get total hits
    CUDA_CHECK_KERNEL(cudaMemcpy(total_hits_output, batch_data.d_global_counter, 
                                 sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    return true;
}

bool launch_improved_minimizer_kernel(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint64_t min_clear_hash_value,
    uint64_t toggle_mask,
    uint32_t* total_hits_output) {
    
    // Use improved kernel with multiple threads per genome
    int threads_per_block = 256;
    size_t shared_mem_size = threads_per_block * (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(bool));
    
    extract_minimizers_kraken2_improved_kernel<<<
        batch_data.max_genomes, 
        threads_per_block, 
        shared_mem_size
    >>>(
        batch_data.d_sequence_data,
        batch_data.d_genome_info,
        batch_data.max_genomes,
        batch_data.d_minimizer_hits,
        batch_data.d_hit_counts,
        batch_data.d_global_counter,
        params,
        min_clear_hash_value,
        toggle_mask,
        batch_data.max_minimizers
    );
    
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
    
    // Get total hits and clamp to maximum
    CUDA_CHECK_KERNEL(cudaMemcpy(total_hits_output, batch_data.d_global_counter, 
                                 sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    if (*total_hits_output > batch_data.max_minimizers) {
        *total_hits_output = batch_data.max_minimizers;
    }
    
    return true;
}

bool launch_lca_computation_kernel(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates,
    int* num_candidates,
    const uint32_t* parent_lookup,
    uint32_t max_taxon_id) {
    
    if (num_hits <= 0) {
        *num_candidates = 0;
        return true;
    }
    
    LaunchConfig config = calculate_optimal_launch_config(num_hits, 256, 0);
    
    if (parent_lookup != nullptr) {
        // Use LCA computation with taxonomy
        dim3 grid(config.blocks_x, config.blocks_y, config.blocks_z);
        dim3 block(config.threads_x, config.threads_y, config.threads_z);
        compute_lca_for_minimizers_kernel<<<grid, block>>>(
            hits, num_hits, candidates, parent_lookup, max_taxon_id
        );
    } else {
        // Simple conversion without LCA computation
        dim3 grid(config.blocks_x, config.blocks_y, config.blocks_z);
        dim3 block(config.threads_x, config.threads_y, config.threads_z);
        convert_hits_to_candidates_kernel<<<grid, block>>>(
            hits, num_hits, candidates
        );
    }
    
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
    
    *num_candidates = num_hits;
    return true;
}

bool launch_memory_initialization_kernel_impl(
    char* sequence_buffer,
    size_t sequence_size,
    GPUGenomeInfo* genome_buffer,
    size_t genome_count,
    GPUMinimizerHit* hit_buffer,
    size_t hit_count) {
    
    size_t max_elements = std::max({sequence_size, genome_count, hit_count});
    LaunchConfig config = calculate_optimal_launch_config(max_elements, 256, 0);
    
    dim3 grid(config.blocks_x, config.blocks_y, config.blocks_z);
    dim3 block(config.threads_x, config.threads_y, config.threads_z);
    initialize_gpu_memory_kernel<<<grid, block>>>(
        sequence_buffer, genome_buffer, hit_buffer,
        sequence_size, genome_count, hit_count
    );
    
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
    return true;
}

bool launch_memory_validation_kernel_impl(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes) {
    
    bool* d_validation_results;
    CUDA_CHECK_KERNEL(cudaMalloc(&d_validation_results, num_genomes * sizeof(bool)));
    
    LaunchConfig config = calculate_optimal_launch_config(num_genomes, 256, 0);
    
    dim3 grid(config.blocks_x, config.blocks_y, config.blocks_z);
    dim3 block(config.threads_x, config.threads_y, config.threads_z);
    validate_memory_integrity_kernel<<<grid, block>>>(
        sequence_data, genome_info, num_genomes, d_validation_results
    );
    
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
    
    // Check results - use std::vector<uint8_t> instead of std::vector<bool>
    // because std::vector<bool> doesn't have data() method
    std::vector<uint8_t> host_results(num_genomes);
    CUDA_CHECK_KERNEL(cudaMemcpy(host_results.data(), d_validation_results, 
                                 num_genomes * sizeof(bool), cudaMemcpyDeviceToHost));
    
    cudaFree(d_validation_results);
    
    // Return true if all validations passed
    for (uint8_t result : host_results) {
        if (!result) return false;
    }
    
    return true;
}

bool check_kernel_execution_errors(const char* kernel_name) {
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        fprintf(stderr, "CUDA error in %s: %s\n", kernel_name, cudaGetErrorString(error));
        return false;
    }
    return true;
}

void print_kernel_launch_info(const LaunchConfig& config, const char* kernel_name) {
    printf("Launching %s with:\n", kernel_name);
    printf("  Grid size: (%d, %d, %d)\n", config.blocks_x, config.blocks_y, config.blocks_z);
    printf("  Block size: (%d, %d, %d)\n", config.threads_x, config.threads_y, config.threads_z);
    printf("  Shared memory: %zu bytes\n", config.shared_memory_bytes);
}

bool validate_kernel_parameters(const GPUBatchData& batch_data) {
    if (batch_data.d_sequence_data == nullptr) return false;
    if (batch_data.d_genome_info == nullptr) return false;
    if (batch_data.d_minimizer_hits == nullptr) return false;
    if (batch_data.max_genomes <= 0) return false;
    if (batch_data.max_minimizers <= 0) return false;
    
    return true;
}

// ===========================
// Legacy CUDA error-returning wrappers for compatibility
// ===========================

cudaError_t launch_memory_initialization_kernel(
    char* d_sequence_buffer,
    size_t sequence_buffer_size,
    GPUGenomeInfo* d_genome_info,
    size_t num_genomes,
    GPUMinimizerHit* d_minimizer_buffer,
    size_t num_minimizers,
    cudaStream_t stream) {
    
    bool result = launch_memory_initialization_kernel_impl(
        d_sequence_buffer,
        sequence_buffer_size,
        d_genome_info,
        num_genomes,
        d_minimizer_buffer,
        num_minimizers
    );
    
    return result ? cudaSuccess : cudaErrorUnknown;
}

cudaError_t launch_memory_validation_kernel(
    const char* d_sequence_buffer,
    const GPUGenomeInfo* d_genome_info,
    int num_genomes,
    bool* h_validation_result,
    cudaStream_t stream) {
    
    bool result = launch_memory_validation_kernel_impl(
        d_sequence_buffer,
        d_genome_info,
        num_genomes
    );
    
    if (h_validation_result) {
        *h_validation_result = result;
    }
    
    return cudaSuccess;
}
