// gpu/gpu_database_kernels.cu
// CUDA kernel implementations for database building operations
// Extracted from monolithic gpu_kraken_database_builder.cu

#include "gpu_database_kernels.h"
#include "../gpu_kraken_types.h"
#include "gpu_minimizer_extraction.cuh"
#include "../processing/contamination_detector.h"
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
// Feature Calculation Device Functions
// ===========================

__device__ uint8_t calculate_gc_category(const char* sequence, int k) {
    // Calculate GC content percentage and map to 3-bit category (0-7)
    int gc_count = 0;
    int valid_bases = 0;
    
    for (int i = 0; i < k; i++) {
        char c = sequence[i];
        // Convert to uppercase if needed
        if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        
        if (c == 'G' || c == 'C') {
            gc_count++;
            valid_bases++;
        } else if (c == 'A' || c == 'T') {
            valid_bases++;
        }
    }
    
    if (valid_bases == 0) return 0;
    
    // Calculate percentage and map to category
    float gc_percent = (float)gc_count / valid_bases * 100.0f;
    
    // Map GC% to 3-bit category (0-7)
    // 0: 0-20%, 1: 20-30%, 2: 30-40%, 3: 40-50%
    // 4: 50-60%, 5: 60-70%, 6: 70-80%, 7: 80-100%
    if (gc_percent < 20.0f) return 0;
    else if (gc_percent < 30.0f) return 1;
    else if (gc_percent < 40.0f) return 2;
    else if (gc_percent < 50.0f) return 3;
    else if (gc_percent < 60.0f) return 4;
    else if (gc_percent < 70.0f) return 5;
    else if (gc_percent < 80.0f) return 6;
    else return 7;
}

__device__ uint8_t calculate_sequence_complexity(const char* sequence, int k) {
    // Calculate sequence complexity using Shannon entropy approximation
    // Count frequency of each nucleotide
    int counts[4] = {0, 0, 0, 0}; // A, C, G, T
    int total = 0;
    
    for (int i = 0; i < k; i++) {
        char c = sequence[i];
        // Convert to uppercase if needed
        if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        
        switch (c) {
            case 'A': counts[0]++; total++; break;
            case 'C': counts[1]++; total++; break;
            case 'G': counts[2]++; total++; break;
            case 'T': counts[3]++; total++; break;
        }
    }
    
    if (total == 0) return 0;
    
    // Calculate Shannon entropy
    float entropy = 0.0f;
    for (int i = 0; i < 4; i++) {
        if (counts[i] > 0) {
            float p = (float)counts[i] / total;
            entropy -= p * log2f(p);
        }
    }
    
    // Map entropy to 3-bit score (0-7)
    // Max entropy for 4 symbols is 2.0
    // 0: very low complexity, 7: high complexity
    float normalized_entropy = entropy / 2.0f;
    uint8_t complexity_score = (uint8_t)(normalized_entropy * 7.0f + 0.5f);
    
    return min(complexity_score, (uint8_t)7);
}

__device__ bool check_position_clustering(uint32_t current_pos, uint32_t last_pos, int window_size) {
    // Simple clustering check: are minimizers within window_size of each other?
    // Note: This is a simplified version. In practice, you might want to track
    // multiple positions in shared memory for better clustering detection
    if (last_pos == UINT32_MAX) return false;
    
    int distance = abs((int)current_pos - (int)last_pos);
    return distance <= window_size;
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
        uint32_t last_position = UINT32_MAX;
        uint32_t local_count = 0;
        const int clustering_window = 100;  // Window size for position clustering detection
        
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
                    hit.taxon_id = static_cast<uint16_t>(genome.taxon_id);
                    hit.position = pos;
                    hit.genome_id = genome.genome_id;
                    hit.strand = MinimizerFlags::STRAND_FORWARD | MinimizerFlags::CLASSIFICATION_UNIQUE;  // Default to unique
                    
                    // Initialize ML weight to maximum (1.0 scaled to uint16_t)
                    hit.ml_weight = 65535;
                    
                    // Calculate features for the k-mer at this position
                    const char* kmer_seq = sequence + pos;
                    
                    // Calculate GC content category
                    uint8_t gc_category = calculate_gc_category(kmer_seq, params.k);
                    
                    // Calculate sequence complexity
                    uint8_t complexity_score = calculate_sequence_complexity(kmer_seq, params.k);
                    
                    // Check for position clustering
                    bool is_clustered = check_position_clustering(pos, last_position, clustering_window);
                    
                    // Encode features into feature_flags
                    uint16_t feature_flags = 0;
                    feature_flags = MinimizerFlags::set_gc_content_category(feature_flags, gc_category);
                    feature_flags = MinimizerFlags::set_complexity_score(feature_flags, complexity_score);
                    feature_flags = MinimizerFlags::set_position_bias(feature_flags, is_clustered);
                    
                    // Note: Contamination risk flag would typically be set based on 
                    // comparison with known contaminant databases, not implemented here
                    feature_flags = MinimizerFlags::set_contamination_risk(feature_flags, false);
                    
                    hit.feature_flags = feature_flags;
                    
                    minimizer_hits[global_pos] = hit;
                    local_count++;
                    last_minimizer = shared_minimizers[responsible_thread];
                    last_position = pos;
                } else {
                    // Hit the limit, stop processing
                    break;
                }
            }
        }
        
        hit_counts_per_genome[genome_id] = local_count;
    }
}

// Fixed kernel with global work distribution
__global__ void extract_minimizers_sliding_window_kernel_fixed(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    int max_minimizers,
    size_t sequence_buffer_size) {
    
    // Global thread ID
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total_threads = blockDim.x * gridDim.x;
    
    // Debug output from first thread
    if (tid == 0) {
        printf("KERNEL STARTED: num_genomes=%d, max_minimizers=%d\n", num_genomes, max_minimizers);
        printf("Block size: %d, Grid size: %d, Total threads: %d\n", blockDim.x, gridDim.x, total_threads);
    }
    
    // Calculate total work items (total k-mers across all genomes)
    uint32_t total_work = 0;
    for (int g = 0; g < num_genomes; g++) {
        if (genome_info[g].sequence_length >= params.k) {
            total_work += genome_info[g].sequence_length - params.k + 1;
        }
        // Debug first few genomes
        if (tid == 0 && g < 3) {
            printf("  Genome %d: length=%u, offset=%u\n", 
                   g, genome_info[g].sequence_length, genome_info[g].sequence_offset);
        }
    }
    
    if (tid == 0) {
        printf("Total k-mers to process: %u\n", total_work);
    }
    
    // Track minimizers found by this thread
    uint32_t thread_minimizer_count = 0;
    
    // Each thread processes multiple k-mers using stride pattern
    for (uint32_t work_idx = tid; work_idx < total_work; work_idx += total_threads) {
        // Find which genome and position this work item corresponds to
        uint32_t cumulative_kmers = 0;
        int genome_idx = -1;
        uint32_t local_kmer_idx = 0;
        
        for (int g = 0; g < num_genomes; g++) {
            if (genome_info[g].sequence_length < params.k) continue;
            
            uint32_t genome_kmers = genome_info[g].sequence_length - params.k + 1;
            if (work_idx < cumulative_kmers + genome_kmers) {
                genome_idx = g;
                local_kmer_idx = work_idx - cumulative_kmers;
                break;
            }
            cumulative_kmers += genome_kmers;
        }
        
        if (genome_idx < 0) continue;
        
        const GPUGenomeInfo& genome = genome_info[genome_idx];
        const char* sequence = sequence_data + genome.sequence_offset;
        
        // Multiple bounds checks for safety
        if (genome.sequence_offset >= sequence_buffer_size) {
            continue;
        }
        if (local_kmer_idx >= genome.sequence_length) {
            continue;
        }
        if (genome.sequence_offset + local_kmer_idx + params.k > sequence_buffer_size) {
            continue;
        }
        if (local_kmer_idx + params.k > genome.sequence_length) {
            continue;
        }
        
        // Debug first few k-mers
        if (work_idx < 3) {
            printf("Thread %d processing k-mer %u from genome %d at position %u\n", 
                   tid, work_idx, genome_idx, local_kmer_idx);
        }
        
        // Extract minimizer at this position
        uint64_t current_minimizer = extract_minimizer_sliding_window(
            sequence, local_kmer_idx, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (current_minimizer != UINT64_MAX) {
            // Simple deduplication - in production, would need better approach
            // For now, just output all valid minimizers
            uint32_t global_pos = atomicAdd(global_hit_counter, 1);
            
            if (global_pos < max_minimizers) {
                GPUMinimizerHit hit;
                hit.minimizer_hash = current_minimizer;
                hit.taxon_id = static_cast<uint16_t>(genome.taxon_id);
                hit.position = local_kmer_idx;
                hit.genome_id = genome_idx;
                hit.strand = MinimizerFlags::STRAND_FORWARD | MinimizerFlags::CLASSIFICATION_UNIQUE;
                hit.ml_weight = 65535;
                hit.feature_flags = 0; // Skip features for now
                
                minimizer_hits[global_pos] = hit;
                thread_minimizer_count++;
                
                // Debug output for first few minimizers
                if (global_pos < 3) {
                    printf("Found minimizer %llu at genome %d, position %u\n", 
                           current_minimizer, genome_idx, local_kmer_idx);
                }
            }
        }
    }
    
    // Debug output from a few threads
    if (tid < 5 && thread_minimizer_count > 0) {
        printf("Thread %d found %u minimizers\n", tid, thread_minimizer_count);
    }
    
    // Final debug output
    if (tid == 0) {
        printf("KERNEL COMPLETE\n");
    }
}

// Keep the old kernel for now but rename it
__global__ void extract_minimizers_sliding_window_kernel_old(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    int max_minimizers,
    size_t sequence_buffer_size) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Debug output from first thread
    if (idx == 0) {
        printf("OLD KERNEL: num_genomes=%d, max_minimizers=%d\n", num_genomes, max_minimizers);
    }
    
    if (idx >= num_genomes) return;
    
    // ... rest of old implementation
}

// Alternative: Multi-threaded per genome approach with shared memory
__global__ void extract_minimizers_multi_thread_per_genome_kernel(
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
    
    // Use one block per genome, multiple threads per block
    int genome_idx = blockIdx.x;
    int thread_id = threadIdx.x;
    int threads_per_block = blockDim.x;
    
    if (genome_idx >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[genome_idx];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) {
        if (thread_id == 0) hit_counts_per_genome[genome_idx] = 0;
        return;
    }
    
    // Shared memory for deduplication
    extern __shared__ uint64_t shared_mem[];
    uint64_t* last_minimizers = shared_mem;
    
    // Initialize shared memory
    if (thread_id < threads_per_block) {
        last_minimizers[thread_id] = UINT64_MAX;
    }
    __syncthreads();
    
    uint32_t total_kmers = seq_length - params.k + 1;
    uint32_t kmers_per_thread = (total_kmers + threads_per_block - 1) / threads_per_block;
    uint32_t start_pos = thread_id * kmers_per_thread;
    uint32_t end_pos = min(start_pos + kmers_per_thread, total_kmers);
    
    uint32_t local_count = 0;
    uint64_t prev_minimizer = (thread_id > 0) ? last_minimizers[thread_id - 1] : UINT64_MAX;
    
    for (uint32_t pos = start_pos; pos < end_pos; pos++) {
        uint64_t current_minimizer = extract_minimizer_sliding_window(
            sequence, pos, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (current_minimizer == UINT64_MAX) continue;
        
        // Apply subsampling if enabled
        if (min_clear_hash_value > 0 && current_minimizer < min_clear_hash_value) {
            continue;
        }
        
        // Apply toggle mask
        current_minimizer ^= toggle_mask;
        
        // Deduplication - only output if different from previous
        if (current_minimizer != prev_minimizer) {
            uint32_t global_pos = atomicAdd(global_hit_counter, 1);
            
            if (global_pos < max_minimizers) {
                GPUMinimizerHit hit;
                hit.minimizer_hash = current_minimizer;
                hit.taxon_id = static_cast<uint16_t>(genome.taxon_id);
                hit.position = pos;
                hit.genome_id = genome_idx;
                hit.strand = MinimizerFlags::STRAND_FORWARD | MinimizerFlags::CLASSIFICATION_UNIQUE;
                hit.ml_weight = 65535;
                
                // Calculate features safely
                const char* kmer_seq = sequence + pos;
                if (pos + params.k <= seq_length) {
                    uint8_t gc_category = calculate_gc_category(kmer_seq, params.k);
                    uint8_t complexity_score = calculate_sequence_complexity(kmer_seq, params.k);
                    
                    uint16_t feature_flags = 0;
                    feature_flags = MinimizerFlags::set_gc_content_category(feature_flags, gc_category);
                    feature_flags = MinimizerFlags::set_complexity_score(feature_flags, complexity_score);
                    hit.feature_flags = feature_flags;
                } else {
                    hit.feature_flags = 0;
                }
                
                minimizer_hits[global_pos] = hit;
                local_count++;
            }
            
            prev_minimizer = current_minimizer;
        }
    }
    
    // Store last minimizer for next thread
    if (thread_id < threads_per_block - 1) {
        last_minimizers[thread_id] = prev_minimizer;
    }
    
    // Sum up counts for this genome
    __shared__ uint32_t block_count;
    if (thread_id == 0) block_count = 0;
    __syncthreads();
    
    atomicAdd(&block_count, local_count);
    __syncthreads();
    
    if (thread_id == 0) {
        hit_counts_per_genome[genome_idx] = block_count;
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

// Option 2: Use the multi-thread per genome kernel
bool launch_improved_minimizer_kernel_fixed(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint64_t min_clear_hash_value,
    uint64_t toggle_mask,
    uint32_t* total_hits_output) {
    
    // Use one block per genome, multiple threads per block
    int threads_per_block = 256;
    int grid_size = batch_data.max_genomes;
    size_t shared_mem_size = threads_per_block * sizeof(uint64_t);
    
    std::cout << "Launching multi-thread kernel with:" << std::endl;
    std::cout << "  Grid size (genomes): " << grid_size << std::endl;
    std::cout << "  Block size (threads): " << threads_per_block << std::endl;
    std::cout << "  Shared memory: " << shared_mem_size << " bytes" << std::endl;
    
    // Allocate hit counts if not provided
    uint32_t* d_hit_counts = batch_data.d_hit_counts;
    if (!d_hit_counts) {
        cudaMalloc(&d_hit_counts, batch_data.max_genomes * sizeof(uint32_t));
        cudaMemset(d_hit_counts, 0, batch_data.max_genomes * sizeof(uint32_t));
    }
    
    // Clear any previous errors
    cudaGetLastError();
    
    extract_minimizers_multi_thread_per_genome_kernel<<<grid_size, threads_per_block, shared_mem_size>>>(
        batch_data.d_sequence_data,
        batch_data.d_genome_info,
        batch_data.max_genomes,
        batch_data.d_minimizer_hits,
        d_hit_counts,
        batch_data.d_global_counter,
        params,
        min_clear_hash_value,
        toggle_mask,
        batch_data.max_minimizers
    );
    
    cudaError_t launch_err = cudaGetLastError();
    if (launch_err != cudaSuccess) {
        std::cerr << "Kernel launch failed: " << cudaGetErrorString(launch_err) << std::endl;
        if (!batch_data.d_hit_counts && d_hit_counts) cudaFree(d_hit_counts);
        return false;
    }
    
    cudaDeviceSynchronize();
    
    // Get total hits
    cudaMemcpy(total_hits_output, batch_data.d_global_counter, 
               sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    // Clean up temporary allocation
    if (!batch_data.d_hit_counts && d_hit_counts) {
        cudaFree(d_hit_counts);
    }
    
    return true;
}

bool launch_minimizer_extraction_kernel(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint32_t* total_hits_output) {
    
    // Use the multi-thread per genome approach which is more stable
    // Pass 0 for min_clear_hash_value and toggle_mask to disable subsampling
    return launch_improved_minimizer_kernel_fixed(batch_data, params, 0, 0, total_hits_output);
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

// Kernel to classify minimizers as unique/canonical/redundant based on frequency
__global__ void classify_minimizers_by_frequency_kernel(
    GPUMinimizerHit* minimizer_hits,
    int num_hits,
    const uint32_t* minimizer_counts_per_species,  // Array indexed by [species_id * max_minimizers + minimizer_idx]
    uint32_t max_species_id,
    uint32_t canonical_threshold,    // e.g., minimizers appearing in >50% of genomes are canonical
    uint32_t redundant_threshold) {  // e.g., minimizers appearing in >80% of genomes are redundant
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_hits) return;
    
    GPUMinimizerHit& hit = minimizer_hits[idx];
    uint32_t species_id = hit.taxon_id;  // Assuming taxon_id represents species
    
    // For now, we'll use a simplified classification based on hash value
    // In a real implementation, you'd track frequency across genomes
    uint64_t hash_mod = hit.minimizer_hash % 100;
    
    // Clear existing classification bits
    hit.strand &= ~MinimizerFlags::CLASSIFICATION_MASK;
    
    // Classify based on simplified criteria
    if (hash_mod < 10) {
        // ~10% are redundant (very common)
        hit.strand |= MinimizerFlags::CLASSIFICATION_REDUNDANT;
    } else if (hash_mod < 40) {
        // ~30% are canonical (moderately common) 
        hit.strand |= MinimizerFlags::CLASSIFICATION_CANONICAL;
    } else {
        // ~60% remain unique
        hit.strand |= MinimizerFlags::CLASSIFICATION_UNIQUE;
    }
}

// Host wrapper for minimizer classification
bool classify_minimizers_by_frequency(
    GPUMinimizerHit* d_minimizer_hits,
    int num_hits,
    uint32_t canonical_threshold,
    uint32_t redundant_threshold) {
    
    if (num_hits <= 0) return true;
    
    LaunchConfig config = calculate_optimal_launch_config(num_hits, 256, 0);
    
    dim3 grid(config.blocks_x, config.blocks_y, config.blocks_z);
    dim3 block(config.threads_x, config.threads_y, config.threads_z);
    
    classify_minimizers_by_frequency_kernel<<<grid, block>>>(
        d_minimizer_hits,
        num_hits,
        nullptr,  // For now, pass nullptr for counts
        0,        // max_species_id
        canonical_threshold,
        redundant_threshold
    );
    
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
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

// ===========================
// Contamination Detection Integration
// ===========================

bool launch_minimizer_extraction_with_contamination_check(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint32_t* total_hits_output,
    const ContaminationConfig& contamination_config) {
    
    // First, extract minimizers as usual
    if (!launch_minimizer_extraction_kernel(batch_data, params, total_hits_output)) {
        return false;
    }
    
    // If no hits, nothing to check
    if (*total_hits_output == 0) {
        return true;
    }
    
    // Apply contamination detection
    ContaminationDetector detector(contamination_config);
    
    // Load contamination database if path is provided
    if (!contamination_config.contamination_db_path.empty()) {
        if (!detector.load_contamination_database(contamination_config.contamination_db_path)) {
            std::cerr << "Warning: Failed to load contamination database, using built-in patterns only" << std::endl;
        }
    }
    
    // Mark contamination in the extracted minimizers
    if (!detector.mark_contamination_in_batch(
            batch_data.d_minimizer_hits, 
            *total_hits_output, 
            ContaminationRiskLevel::LOW_RISK)) {
        std::cerr << "Warning: Contamination detection failed" << std::endl;
        // Non-fatal - continue with unmarked minimizers
    }
    
    // Get statistics
    const auto& stats = detector.get_statistics();
    std::cout << "Contamination detection complete:" << std::endl;
    std::cout << "  Human contamination: " << stats.human_contamination_count << std::endl;
    std::cout << "  Adapter contamination: " << stats.adapter_contamination_count << std::endl;
    std::cout << "  Low complexity: " << stats.low_complexity_count << std::endl;
    
    // Export contamination report if in debug mode
    if (contamination_config.contamination_db_path.find("debug") != std::string::npos) {
        detector.export_contamination_report("contamination_report.txt");
    }
    
    return true;
}
