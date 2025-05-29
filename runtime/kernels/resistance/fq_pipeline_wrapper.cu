#include "fq_mutation_detector.cuh"
#include <cuda_runtime.h>
#include <iostream>

// Wrapper functions to call CUDA kernels from C++

extern "C" void launch_kmer_filter(
    const char* d_reads,
    const int* d_read_lengths,
    const int* d_read_offsets,
    FQMutationDetectorCUDA& detector,
    int num_reads,
    CandidateMatch* d_candidates,
    uint32_t* d_candidate_counts
) {
    dim3 block(256);
    dim3 grid((num_reads + block.x - 1) / block.x);
    
    kmer_filter_kernel<<<grid, block>>>(
        d_reads, d_read_lengths, d_read_offsets,
        detector.d_kmer_index, detector.d_kmer_sorted, detector.d_kmer_positions,
        detector.num_kmers, num_reads,
        d_candidates, d_candidate_counts, MAX_CANDIDATES_PER_READ
    );
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error in kmer_filter_kernel: " << cudaGetErrorString(err) << std::endl;
    }
}

extern "C" void launch_position_weighted_alignment(
    const char* d_reads,
    const int* d_read_lengths,
    const int* d_read_offsets,
    const CandidateMatch* d_candidates,
    const uint32_t* d_candidate_counts,
    FQMutationDetectorCUDA& detector,
    int num_reads,
    AlignmentResult* d_results,
    uint32_t* d_result_count
) {
    dim3 block(128);
    dim3 grid((num_reads + block.x - 1) / block.x);
    size_t shared_mem_size = 4096;
    
    position_weighted_alignment_kernel<<<grid, block, shared_mem_size>>>(
        d_reads, d_read_lengths, d_read_offsets,
        d_candidates, d_candidate_counts, MAX_CANDIDATES_PER_READ,
        detector.d_reference_sequences, detector.d_ref_lengths, detector.d_ref_offsets,
        detector.d_position_weights, detector.d_mutation_masks,
        detector.d_mutation_info, detector.d_mutation_counts,
        detector.align_params, num_reads,
        d_results, d_result_count
    );
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error in position_weighted_alignment_kernel: " << cudaGetErrorString(err) << std::endl;
    }
}