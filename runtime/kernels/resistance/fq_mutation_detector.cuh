#ifndef FQ_MUTATION_DETECTOR_CUH
#define FQ_MUTATION_DETECTOR_CUH

#include <cuda_runtime.h>
#include <stdint.h>

// Constants
#define KMER_LENGTH 15
#define MAX_READ_LENGTH 300
#define MAX_REF_LENGTH 2000
#define WARP_SIZE 32
#define BLOCK_SIZE 256
#define MAX_CANDIDATES_PER_READ 64
#define MAX_MUTATIONS_PER_GENE 20

// Base encoding
#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_T 3
#define BASE_N 4

// Scoring parameters
struct AlignmentParams {
    float match_score;
    float mismatch_score;
    float gap_open;
    float gap_extend;
    float mutation_match_bonus;
    float mutation_mismatch_penalty;
    float min_alignment_score;
    float min_identity;
};

// K-mer index structure
struct KmerEntry {
    uint64_t kmer;
    uint32_t gene_id;
    uint32_t species_id;
    uint32_t seq_id;
    uint16_t position;
};

// Mutation information
struct MutationInfo {
    uint16_t position;
    uint8_t wild_type;
    uint8_t mutant;
    uint16_t codon_position;
};

// Candidate match from stage 1
struct CandidateMatch {
    uint32_t gene_id;
    uint32_t species_id;
    uint32_t seq_id;
    uint16_t kmer_hits;
};

// Alignment result from stage 2
struct AlignmentResult {
    uint32_t read_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint32_t seq_id;
    float alignment_score;
    float identity;
    uint16_t start_pos;
    uint16_t matches;
    uint8_t mutations_detected[MAX_MUTATIONS_PER_GENE];
    uint8_t num_mutations_detected;
    bool passes_threshold;
};

// Forward declarations for kernels
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
);

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
);

// Host-side class declaration
class FQMutationDetectorCUDA {
public:
    // Device memory pointers
    KmerEntry* d_kmer_index;
    uint64_t* d_kmer_sorted;
    uint32_t* d_kmer_positions;
    char* d_reference_sequences;
    uint32_t* d_ref_lengths;
    uint32_t* d_ref_offsets;
    float* d_position_weights;
    uint8_t* d_mutation_masks;
    MutationInfo* d_mutation_info;
    uint32_t* d_mutation_counts;
    
    // Index metadata
    uint32_t num_kmers;
    uint32_t num_sequences;
    uint32_t total_ref_length;
    
    // Alignment parameters
    AlignmentParams align_params;

    FQMutationDetectorCUDA();
    ~FQMutationDetectorCUDA();
    
    void loadIndex(const char* index_path);
    void processReads(const char* r1_path, const char* r2_path, const char* output_path);
};

// Entry point
extern "C" void run_fq_mutation_detection(
    const char* index_path,
    const char* r1_path,
    const char* r2_path,
    const char* output_path
);

// Wrapper functions for calling kernels from C++
extern "C" {
    void launch_kmer_filter(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        FQMutationDetectorCUDA& detector,
        int num_reads,
        CandidateMatch* d_candidates,
        uint32_t* d_candidate_counts
    );

    void launch_position_weighted_alignment(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const CandidateMatch* d_candidates,
        const uint32_t* d_candidate_counts,
        FQMutationDetectorCUDA& detector,
        int num_reads,
        AlignmentResult* d_results,
        uint32_t* d_result_count
    );
}

#endif // FQ_MUTATION_DETECTOR_CUH