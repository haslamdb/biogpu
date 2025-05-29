#include "fq_mutation_detector.cuh"
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cooperative_groups.h>
#include <stdio.h>

namespace cg = cooperative_groups;

// Device functions
__device__ inline uint8_t encode_base(char base) {
    switch(base) {
        case 'A': case 'a': return BASE_A;
        case 'C': case 'c': return BASE_C;
        case 'G': case 'g': return BASE_G;
        case 'T': case 't': return BASE_T;
        default: return BASE_N;
    }
}

__device__ inline uint64_t encode_kmer(const char* seq, int pos) {
    uint64_t kmer = 0;
    for (int i = 0; i < KMER_LENGTH; i++) {
        uint8_t base = encode_base(seq[pos + i]);
        if (base == BASE_N) return UINT64_MAX; // Invalid k-mer
        kmer = (kmer << 2) | base;
    }
    return kmer;
}

// Stage 1: K-mer filtering kernel
__global__ void kmer_filter_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const KmerEntry* kmer_index,
    const uint64_t* kmer_sorted,
    const uint32_t* kmer_positions,
    const uint32_t num_kmers,
    const int num_reads,
    CandidateMatch* candidates,
    uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        
        if (read_len < KMER_LENGTH) continue;
        
        // Local candidate tracking
        CandidateMatch local_candidates[MAX_CANDIDATES_PER_READ];
        uint32_t local_count = 0;
        
        // Extract k-mers from read
        for (int pos = 0; pos <= read_len - KMER_LENGTH; pos++) {
            uint64_t kmer = encode_kmer(read, pos);
            if (kmer == UINT64_MAX) continue;
            
            // Binary search in sorted k-mer index
            uint32_t left = 0;
            uint32_t right = num_kmers;
            
            while (left < right) {
                uint32_t mid = (left + right) / 2;
                if (kmer_sorted[mid] < kmer) {
                    left = mid + 1;
                } else if (kmer_sorted[mid] > kmer) {
                    right = mid;
                } else {
                    // Found match - get all entries with this k-mer
                    uint32_t start = kmer_positions[mid];
                    uint32_t end = (mid + 1 < num_kmers) ? kmer_positions[mid + 1] : num_kmers;
                    
                    for (uint32_t i = start; i < end && local_count < MAX_CANDIDATES_PER_READ; i++) {
                        const KmerEntry& entry = kmer_index[i];
                        
                        // Check if we already have this candidate
                        bool found = false;
                        for (uint32_t j = 0; j < local_count; j++) {
                            if (local_candidates[j].gene_id == entry.gene_id &&
                                local_candidates[j].species_id == entry.species_id &&
                                local_candidates[j].seq_id == entry.seq_id) {
                                local_candidates[j].kmer_hits++;
                                found = true;
                                break;
                            }
                        }
                        
                        if (!found && local_count < MAX_CANDIDATES_PER_READ) {
                            local_candidates[local_count].gene_id = entry.gene_id;
                            local_candidates[local_count].species_id = entry.species_id;
                            local_candidates[local_count].seq_id = entry.seq_id;
                            local_candidates[local_count].kmer_hits = 1;
                            local_count++;
                        }
                    }
                    break;
                }
            }
        }
        
        // Write candidates to global memory
        uint32_t global_offset = read_idx * max_candidates_per_read;
        for (uint32_t i = 0; i < local_count; i++) {
            candidates[global_offset + i] = local_candidates[i];
        }
        candidate_counts[read_idx] = local_count;
    }
}

// Stage 2: Position-weighted alignment kernel
__global__ void position_weighted_alignment_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const CandidateMatch* candidates,
    const uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read,
    const char* reference_sequences,
    const uint32_t* ref_lengths,
    const uint32_t* ref_offsets,
    const float* position_weights,
    const uint8_t* mutation_masks,
    const MutationInfo* mutation_info,
    const uint32_t* mutation_counts,
    const AlignmentParams params,
    const int num_reads,
    AlignmentResult* results,
    uint32_t* result_count
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    // Shared memory for dynamic programming matrix
    extern __shared__ float shared_mem[];
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        uint32_t num_candidates = candidate_counts[read_idx];
        
        if (num_candidates == 0) continue;
        
        // Process each candidate for this read
        for (uint32_t cand_idx = 0; cand_idx < num_candidates; cand_idx++) {
            const CandidateMatch& candidate = candidates[read_idx * max_candidates_per_read + cand_idx];
            
            // Get reference sequence
            uint32_t ref_idx = candidate.seq_id;
            const char* ref_seq = reference_sequences + ref_offsets[ref_idx];
            const uint32_t ref_len = ref_lengths[ref_idx];
            
            // Get position weights and mutation mask for this sequence
            const float* seq_weights = position_weights + ref_offsets[ref_idx];
            const uint8_t* seq_mut_mask = mutation_masks + ref_offsets[ref_idx];
            
            // Sliding window alignment with position weights
            float best_score = -1000.0f;
            int best_pos = -1;
            int best_matches = 0;
            
            // Try different starting positions
            for (int start_pos = max(0, (int)ref_len - read_len - 50); 
                 start_pos < min((int)ref_len - read_len + 50, (int)ref_len); 
                 start_pos++) {
                
                float score = 0.0f;
                int matches = 0;
                
                // Simple scoring (full DP would be more complex)
                for (int i = 0; i < read_len && start_pos + i < ref_len; i++) {
                    uint8_t read_base = encode_base(read[i]);
                    uint8_t ref_base = encode_base(ref_seq[start_pos + i]);
                    
                    if (read_base == BASE_N || ref_base == BASE_N) continue;
                    
                    float pos_weight = seq_weights[start_pos + i];
                    
                    if (read_base == ref_base) {
                        score += params.match_score * pos_weight;
                        matches++;
                        
                        // Bonus for matching at mutation site
                        if (seq_mut_mask[start_pos + i]) {
                            score += params.mutation_match_bonus * pos_weight;
                        }
                    } else {
                        score += params.mismatch_score * pos_weight;
                        
                        // Penalty for mismatch at mutation site
                        if (seq_mut_mask[start_pos + i]) {
                            score -= params.mutation_mismatch_penalty * pos_weight;
                        }
                    }
                }
                
                if (score > best_score) {
                    best_score = score;
                    best_pos = start_pos;
                    best_matches = matches;
                }
            }
            
            // Check if alignment passes threshold
            float identity = (float)best_matches / read_len;
            if (best_score >= params.min_alignment_score && identity >= params.min_identity) {
                // Get result slot
                uint32_t result_idx = atomicAdd(result_count, 1);
                
                AlignmentResult& result = results[result_idx];
                result.read_id = read_idx;
                result.gene_id = candidate.gene_id;
                result.species_id = candidate.species_id;
                result.seq_id = candidate.seq_id;
                result.alignment_score = best_score;
                result.identity = identity;
                result.start_pos = best_pos;
                result.matches = best_matches;
                result.passes_threshold = true;
                
                // Check for mutations
                result.num_mutations_detected = 0;
                const MutationInfo* seq_mutations = mutation_info + ref_idx * MAX_MUTATIONS_PER_GENE;
                uint32_t num_mutations = mutation_counts[ref_idx];
                
                for (uint32_t m = 0; m < num_mutations && result.num_mutations_detected < MAX_MUTATIONS_PER_GENE; m++) {
                    const MutationInfo& mut = seq_mutations[m];
                    
                    // Check if mutation position is covered by alignment
                    if (best_pos <= mut.position && mut.position < best_pos + read_len) {
                        int read_pos = mut.position - best_pos;
                        uint8_t read_base = encode_base(read[read_pos]);
                        uint8_t ref_base = encode_base(ref_seq[mut.position]);
                        
                        // Store mutation detection result
                        result.mutations_detected[result.num_mutations_detected++] = 
                            (read_base == mut.mutant) ? 1 : 0;
                    }
                }
            }
        }
    }
}

// Implementation of FQMutationDetectorCUDA methods
FQMutationDetectorCUDA::FQMutationDetectorCUDA() {
        // Initialize default parameters
        align_params.match_score = 2.0f;
        align_params.mismatch_score = -1.0f;
        align_params.gap_open = -3.0f;
        align_params.gap_extend = -1.0f;
        align_params.mutation_match_bonus = 3.0f;
        align_params.mutation_mismatch_penalty = 2.0f;
        align_params.min_alignment_score = 50.0f;
        align_params.min_identity = 0.8f;
    }

FQMutationDetectorCUDA::~FQMutationDetectorCUDA() {
        // Free device memory
        if (d_kmer_index) cudaFree(d_kmer_index);
        if (d_kmer_sorted) cudaFree(d_kmer_sorted);
        if (d_kmer_positions) cudaFree(d_kmer_positions);
        if (d_reference_sequences) cudaFree(d_reference_sequences);
        if (d_ref_lengths) cudaFree(d_ref_lengths);
        if (d_ref_offsets) cudaFree(d_ref_offsets);
        if (d_position_weights) cudaFree(d_position_weights);
        if (d_mutation_masks) cudaFree(d_mutation_masks);
        if (d_mutation_info) cudaFree(d_mutation_info);
        if (d_mutation_counts) cudaFree(d_mutation_counts);
    }

void FQMutationDetectorCUDA::loadIndex(const char* index_path) {
        // TODO: Load index from HDF5 file
        // This would load the k-mer index, reference sequences, position weights, etc.
        printf("Loading index from: %s\n", index_path);
    }

void FQMutationDetectorCUDA::processReads(const char* r1_path, const char* r2_path, const char* output_path) {
        printf("Processing paired-end reads:\n");
        printf("  R1: %s\n", r1_path);
        printf("  R2: %s\n", r2_path);
        printf("  Output: %s\n", output_path);
        
        // TODO: Implement full pipeline
        // 1. Load reads from FASTQ files
        // 2. Transfer to GPU
        // 3. Run stage 1 k-mer filtering
        // 4. Run stage 2 alignment
        // 5. Merge results from R1 and R2
        // 6. Write results to output file
    }

// Main entry point for testing
extern "C" void run_fq_mutation_detection(
    const char* index_path,
    const char* r1_path,
    const char* r2_path,
    const char* output_path
) {
    FQMutationDetectorCUDA detector;
    detector.loadIndex(index_path);
    detector.processReads(r1_path, r2_path, output_path);
}