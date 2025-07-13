// gpu/gpu_database_kernels.cu
//
// This file contains the final, corrected CUDA kernels for database building.
//
// FIX: The `extract_minimizers_stateful_kernel` has been completely rewritten
// to be a more direct and robust port of the official Kraken 2 MinimizerScanner
// algorithm. This resolves the illegal memory access error by correctly
// managing the sliding window state one base at a time.
//
// FIX: The LCA computation is a full, functional implementation using
// the Thrust library, removing all previous placeholders.

#include "gpu_database_kernels.h"
#include "../gpu_kraken_types.h"
#include "gpu_minimizer_extraction.cuh"
#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/iterator/constant_iterator.h>
#include <cstdio>

// ===========================
// Device Function Implementations
// ===========================

// Finds the Lowest Common Ancestor of two taxa by walking up the taxonomy tree.
__device__ uint32_t find_lca_simple_device(uint32_t taxon1, uint32_t taxon2,
                                          const uint32_t* parent_lookup,
                                          uint32_t max_taxon_id) {
    if (taxon1 == taxon2) return taxon1;
    if (taxon1 == 0 || taxon1 > max_taxon_id) return taxon2;
    if (taxon2 == 0 || taxon2 > max_taxon_id) return taxon1;

    uint32_t path1[64];
    uint32_t path2[64];
    int path1_len = 0;
    int path2_len = 0;

    uint32_t current = taxon1;
    while (current != 0 && current != 1 && path1_len < 64) {
        path1[path1_len++] = current;
        if (current > max_taxon_id) break;
        current = parent_lookup[current];
    }
    path1[path1_len++] = 1; // Add root

    current = taxon2;
    while (current != 0 && current != 1 && path2_len < 64) {
        path2[path2_len++] = current;
        if (current > max_taxon_id) break;
        current = parent_lookup[current];
    }
    path2[path2_len++] = 1; // Add root

    for (int i = 0; i < path1_len; i++) {
        for (int j = 0; j < path2_len; j++) {
            if (path1[i] == path2[j]) {
                return path1[i];
            }
        }
    }

    return 1; // Return root as the LCA if no other is found
}


// ===========================
// CUDA Kernel Implementations
// ===========================

// This kernel correctly implements the stateful, sliding-window minimizer
// extraction algorithm, inspired by the official Kraken 2 source code.
__global__ void extract_minimizers_stateful_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    int max_minimizers,
    uint64_t xor_mask,
    size_t sequence_buffer_size)
{
    // Add shared memory declaration
    extern __shared__ char shared_mem[];
    
    int genome_idx = blockIdx.x;
    if (genome_idx >= num_genomes) return;

    const GPUGenomeInfo& genome = genome_info[genome_idx];
    
    // Add validation
    if (genome.sequence_offset >= UINT32_MAX - genome.sequence_length) return; // Overflow check
    if (genome_idx >= num_genomes) return; // Double-check bounds
    
    // Bounds check for sequence access
    if (genome.sequence_offset + genome.sequence_length > sequence_buffer_size) return;
    
    const char* sequence = sequence_data + genome.sequence_offset;
    const uint32_t seq_length = genome.sequence_length;

    if (seq_length < params.k) return;

    const uint32_t window_len = params.k - params.ell + 1;
    if (window_len > 256) return;
    if (window_len == 0) return;  // Prevent division by zero

    // Use shared memory for arrays
    uint64_t* window_hashes = (uint64_t*)shared_mem;
    uint32_t* window_pos = (uint32_t*)(window_hashes + window_len);
    
    // Process sequence in chunks to avoid issues with very large sequences
    const uint32_t CHUNK_SIZE = 1000000; // 1MB chunks
    const uint32_t OVERLAP = params.k; // Overlap to ensure continuity
    
    unsigned __int128 current_lmer = 0;
    uint64_t lmer_mask = (params.ell < 64) ? ((unsigned __int128)1 << (2 * params.ell)) - 1 : (unsigned __int128)-1;
    uint64_t last_minimizer_hash = UINT64_MAX;
    
    for (uint32_t chunk_start = 0; chunk_start < seq_length; chunk_start += CHUNK_SIZE) {
        uint32_t chunk_end = min(chunk_start + CHUNK_SIZE + OVERLAP, seq_length);
        uint32_t chunk_length = chunk_end - chunk_start;
        
        // Skip if chunk is too small
        if (chunk_length < params.k) continue;
        
        // Reset deque for each chunk
        uint32_t deque_head = 0, deque_tail = 0;
        
        // Prime the first l-1 bases of the chunk
        uint32_t start_pos = (chunk_start == 0) ? 0 : chunk_start;
        for (uint32_t i = start_pos; i < start_pos + params.ell - 1 && i < chunk_end; i++) {
            uint64_t base = DnaTo2Bit(sequence[i]);
            if (base >= 4) {
                current_lmer = 0;
                continue;
            }
            current_lmer = (current_lmer << 2) | base;
        }

        // Initialize the first full window of l-mers
        for (uint32_t i = 0; i < window_len; i++) {
            uint32_t current_pos = start_pos + params.ell - 1 + i;
            if (current_pos >= chunk_end) break;
            uint64_t base = DnaTo2Bit(sequence[current_pos]);
            if (base >= 4) {
                deque_head = deque_tail = 0;
                continue;
            }
            current_lmer = ((current_lmer << 2) | base) & lmer_mask;

            unsigned __int128 rc_lmer = reverse_complement_128(current_lmer, params.ell);
            uint64_t canon_lmer = (uint64_t)(current_lmer < rc_lmer ? current_lmer : rc_lmer);
            uint64_t current_hash = murmur_hash3(canon_lmer ^ xor_mask);

            // **FIX**: Correct circular buffer indexing to prevent underflow
            while (deque_head != deque_tail) {
                uint32_t prev_idx = (deque_tail == 0) ? window_len - 1 : deque_tail - 1;
                if (window_hashes[prev_idx] <= current_hash) break;
                deque_tail = prev_idx;
            }
            window_hashes[deque_tail] = current_hash;
            window_pos[deque_tail] = current_pos - params.ell + 1;
            deque_tail = (deque_tail + 1) % window_len;
        }

        // Output the first minimizer if we have one
        if (deque_head != deque_tail && chunk_start == 0) {
            uint64_t minimizer_hash = window_hashes[deque_head];
            uint32_t write_idx = atomicAdd(global_hit_counter, 1);
            if (write_idx < max_minimizers) {
                // Use proper structure assignment to ensure write completes
                GPUMinimizerHit hit = {minimizer_hash, (uint32_t)genome_idx, window_pos[deque_head], (uint16_t)genome.taxon_id};
                minimizer_hits[write_idx] = hit;
                __threadfence(); // Ensure write is visible
            }
            last_minimizer_hash = minimizer_hash;
        }

        // Slide the window across the rest of the chunk
        uint32_t slide_start = max(start_pos + params.k, chunk_start);
        for (uint32_t i = slide_start; i < chunk_end; i++) {
            if (deque_head != deque_tail && window_pos[deque_head] == (i - params.k)) {
                deque_head = (deque_head + 1) % window_len;
            }

            uint64_t base = DnaTo2Bit(sequence[i]);
            if (base >= 4) {
                deque_head = deque_tail = 0;
                continue;
            }
            current_lmer = ((current_lmer << 2) | base) & lmer_mask;
            
            unsigned __int128 rc_lmer = reverse_complement_128(current_lmer, params.ell);
            uint64_t canon_lmer = (uint64_t)(current_lmer < rc_lmer ? current_lmer : rc_lmer);
            uint64_t current_hash = murmur_hash3(canon_lmer ^ xor_mask);

            // **FIX**: Correct circular buffer indexing to prevent underflow
            while (deque_head != deque_tail) {
                uint32_t prev_idx = (deque_tail == 0) ? window_len - 1 : deque_tail - 1;
                if (window_hashes[prev_idx] <= current_hash) break;
                deque_tail = prev_idx;
            }
            window_hashes[deque_tail] = current_hash;
            window_pos[deque_tail] = i - params.ell + 1;
            deque_tail = (deque_tail + 1) % window_len;

            uint64_t minimizer_hash = window_hashes[deque_head];
            if (minimizer_hash != last_minimizer_hash) {
                uint32_t write_idx = atomicAdd(global_hit_counter, 1);
                if (write_idx < max_minimizers) {
                    // Use proper structure assignment to ensure write completes
                    GPUMinimizerHit hit = {minimizer_hash, (uint32_t)genome_idx, window_pos[deque_head], (uint16_t)genome.taxon_id};
                    minimizer_hits[write_idx] = hit;
                    __threadfence(); // Ensure write is visible
                } else {
                    return; // Exit if we're out of space
                }
                last_minimizer_hash = minimizer_hash;
            }
        }
    } // End of chunk loop
}

// This kernel computes the final LCA for each unique minimizer.
__global__ void compute_lca_from_hits_kernel(
    const GPUMinimizerHit* sorted_hits,
    const uint32_t* run_starts,
    const uint32_t* run_ends,
    int num_unique_minimizers,
    LCACandidate* lca_candidates,
    const uint32_t* parent_lookup,
    uint32_t max_taxon_id)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_unique_minimizers) return;

    uint32_t start_idx = run_starts[tid];
    uint32_t end_idx = run_ends[tid];

    // Compute LCA for this run of identical minimizers
    uint32_t current_lca = sorted_hits[start_idx].taxon_id;
    for (uint32_t i = start_idx + 1; i < end_idx; i++) {
        current_lca = find_lca_simple_device(current_lca, sorted_hits[i].taxon_id, parent_lookup, max_taxon_id);
    }

    // Write the final LCA candidate
    lca_candidates[tid].minimizer_hash = sorted_hits[start_idx].minimizer_hash;
    lca_candidates[tid].lca_taxon = current_lca;
    lca_candidates[tid].genome_count = end_idx - start_idx;
    lca_candidates[tid].uniqueness_score = 1.0f / (float)(end_idx - start_idx);
}

// ===========================
// Host-Side Kernel Launchers
// ===========================

bool launch_minimizer_extraction_kernel(
    const GPUBatchData& batch_data,
    const MinimizerParams& params,
    uint32_t* total_hits_output)
{
    std::cout << "Launching stateful minimizer extraction kernel..." << std::endl;
    if (batch_data.max_genomes == 0) {
        *total_hits_output = 0;
        return true;
    }
    CUDA_CHECK_KERNEL(cudaMemset(batch_data.d_global_counter, 0, sizeof(uint32_t)));

    int threads_per_block = 1; // This kernel is serial per genome
    int grid_size = batch_data.max_genomes;
    uint64_t xor_mask = 0x3c8bfbb395c60474ULL;
    
    // Calculate shared memory size
    uint32_t window_len = params.k - params.ell + 1;
    size_t shared_mem_size = window_len * (sizeof(uint64_t) + sizeof(uint32_t));
    
    // Check if shared memory size is within limits
    if (shared_mem_size > 48 * 1024) { // 48KB limit
        std::cerr << "ERROR: Required shared memory " << shared_mem_size 
                  << " exceeds limit" << std::endl;
        return false;
    }

    extract_minimizers_stateful_kernel<<<grid_size, threads_per_block, shared_mem_size>>>(
        batch_data.d_sequence_data, batch_data.d_genome_info, batch_data.max_genomes,
        batch_data.d_minimizer_hits, batch_data.d_global_counter, params,
        batch_data.max_minimizers, xor_mask, batch_data.sequence_buffer_size
    );

    CUDA_CHECK_KERNEL(cudaGetLastError());
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());
    CUDA_CHECK_KERNEL(cudaMemcpy(total_hits_output, batch_data.d_global_counter,
                                 sizeof(uint32_t), cudaMemcpyDeviceToHost));

    if (*total_hits_output >= batch_data.max_minimizers) {
        std::cerr << "WARNING: Minimizer capacity reached. Found " << *total_hits_output
                  << " but buffer size is " << batch_data.max_minimizers << ". Clamping result." << std::endl;
        *total_hits_output = batch_data.max_minimizers;
    }
    std::cout << "✓ Kernel finished. Found " << *total_hits_output << " minimizers." << std::endl;
    return true;
}

bool launch_lca_computation_kernel(
    GPUMinimizerHit* hits, int num_hits,
    LCACandidate* candidates, int* num_candidates,
    const uint32_t* parent_lookup, uint32_t max_taxon_id)
{
    if (num_hits == 0) {
        *num_candidates = 0;
        return true;
    }
    std::cout << "Computing LCAs for " << num_hits << " minimizer hits..." << std::endl;

    thrust::device_ptr<GPUMinimizerHit> hits_ptr(hits);

    // 1. Sort hits by minimizer_hash to group identical minimizers together.
    thrust::sort(hits_ptr, hits_ptr + num_hits,
        [] __device__ (const GPUMinimizerHit& a, const GPUMinimizerHit& b) {
            return a.minimizer_hash < b.minimizer_hash;
        });
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());

    // 2. Find the start and end of each run of identical minimizer hashes.
    thrust::device_vector<uint64_t> keys(num_hits);
    thrust::transform(hits_ptr, hits_ptr + num_hits, keys.begin(),
        [] __device__ (const GPUMinimizerHit& hit) { return hit.minimizer_hash; });

    thrust::device_vector<uint32_t> run_starts(num_hits);
    thrust::device_vector<uint32_t> run_ends(num_hits);
    thrust::device_vector<uint64_t> unique_keys(num_hits);

    auto unique_end = thrust::reduce_by_key(
        keys.begin(), keys.end(),
        thrust::make_constant_iterator(1),
        unique_keys.begin(),
        run_ends.begin() // Use run_ends to store counts temporarily
    );

    int num_unique = unique_end.first - unique_keys.begin();
    *num_candidates = num_unique;
    std::cout << "Found " << num_unique << " unique minimizers." << std::endl;

    // Compute run_starts (exclusive scan) and run_ends (inclusive scan)
    thrust::exclusive_scan(run_ends.begin(), run_ends.begin() + num_unique, run_starts.begin());
    thrust::inclusive_scan(run_ends.begin(), run_ends.begin() + num_unique, run_ends.begin());

    // 3. Launch the kernel to compute LCA for each run.
    int threads = 256;
    int blocks = (num_unique + threads - 1) / threads;
    compute_lca_from_hits_kernel<<<blocks, threads>>>(
        hits,
        thrust::raw_pointer_cast(run_starts.data()),
        thrust::raw_pointer_cast(run_ends.data()),
        num_unique,
        candidates,
        parent_lookup,
        max_taxon_id
    );
    CUDA_CHECK_KERNEL(cudaDeviceSynchronize());

    std::cout << "✓ LCA computation completed." << std::endl;
    return true;
}

bool launch_memory_initialization_kernel_impl(
    char* sequence_buffer, size_t sequence_size,
    GPUGenomeInfo* genome_buffer, size_t genome_count,
    GPUMinimizerHit* hit_buffer, size_t hit_count) {

    if (sequence_buffer) CUDA_CHECK_KERNEL(cudaMemset(sequence_buffer, 0, sequence_size));
    if (genome_buffer) CUDA_CHECK_KERNEL(cudaMemset(genome_buffer, 0, genome_count * sizeof(GPUGenomeInfo)));
    if (hit_buffer) CUDA_CHECK_KERNEL(cudaMemset(hit_buffer, 0, hit_count * sizeof(GPUMinimizerHit)));
    return true;
}
