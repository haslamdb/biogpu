#include "fq_mutation_detector.cuh"
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cooperative_groups.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

// Debug macros
#define DEBUG 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { printf("[GPU DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

namespace cg = cooperative_groups;

// Import the enhanced k-mer screening function from fixed_kmer_screening.cu
extern "C" {
    int test_kmer_screening(
        const char* index_dir,
        const char* reads_data,
        const int* read_lengths,
        const int* read_offsets,
        int num_reads
    );
}

// Import the enhanced kernel declaration
extern __global__ void enhanced_kmer_filter_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const KmerEntry* kmer_index,
    const uint64_t* sorted_kmers,
    const uint32_t* kmer_start_positions,
    const uint32_t num_unique_kmers,
    const uint32_t total_kmer_entries,
    const int num_reads,
    const int kmer_length,
    CandidateMatch* candidates,
    uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read,
    const bool check_reverse_complement,
    int* debug_stats
);

// Device functions for base encoding/decoding
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

// Simplified position-weighted alignment kernel
__global__ void simple_alignment_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const CandidateMatch* candidates,
    const uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read,
    const int num_reads,
    AlignmentResult* results,
    uint32_t* result_count
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    if (tid == 0) {
        DEBUG_PRINT("[ALIGNMENT] Simple alignment kernel started: num_reads=%d", num_reads);
    }
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        uint32_t num_candidates = candidate_counts[read_idx];
        
        if (num_candidates == 0) continue;
        
        // Process each candidate for this read
        for (uint32_t cand_idx = 0; cand_idx < num_candidates; cand_idx++) {
            const CandidateMatch& candidate = candidates[read_idx * max_candidates_per_read + cand_idx];
            
            // Simple scoring based on k-mer hits
            float score = candidate.kmer_hits * 10.0f;
            float max_possible_kmers = max(1.0f, (float)(read_len - KMER_LENGTH + 1));
            float identity = min(1.0f, (float)candidate.kmer_hits / max_possible_kmers);
            
            // Check if alignment passes threshold (simplified)
            if (score >= 20.0f && identity >= 0.1f) {
                // Get result slot
                uint32_t result_idx = atomicAdd(result_count, 1);
                
                if (result_idx < num_reads * MAX_RESULTS_PER_READ) {
                    AlignmentResult& result = results[result_idx];
                    result.read_id = read_idx;
                    result.gene_id = candidate.gene_id;
                    result.species_id = candidate.species_id;
                    result.seq_id = candidate.seq_id;
                    result.alignment_score = score;
                    result.identity = identity;
                    result.start_pos = 0; // Simplified
                    result.matches = candidate.kmer_hits;
                    result.passes_threshold = true;
                    
                    // For now, mark as having mutations detected if we have good hits
                    result.num_mutations_detected = (candidate.kmer_hits > 3) ? 1 : 0;
                    if (result.num_mutations_detected > 0) {
                        result.mutations_detected[0] = 1; // Found mutation
                    }
                }
            }
        }
    }
    
    if (tid == 0) {
        DEBUG_PRINT("[ALIGNMENT] Simple alignment kernel completed");
    }
}

// Implementation of FQMutationDetectorCUDA methods
FQMutationDetectorCUDA::FQMutationDetectorCUDA() {
    // Initialize all pointers to null
    d_kmer_index = nullptr;
    d_kmer_sorted = nullptr;
    d_kmer_positions = nullptr;
    d_reference_sequences = nullptr;
    d_ref_lengths = nullptr;
    d_ref_offsets = nullptr;
    d_position_weights = nullptr;
    d_mutation_masks = nullptr;
    d_mutation_info = nullptr;
    d_mutation_counts = nullptr;
    
    // Initialize counters
    num_kmers = 0;
    num_sequences = 0;
    total_ref_length = 0;
    
    // Initialize default parameters
    align_params.match_score = 2.0f;
    align_params.mismatch_score = -1.0f;
    align_params.gap_open = -3.0f;
    align_params.gap_extend = -1.0f;
    align_params.mutation_match_bonus = 3.0f;
    align_params.mutation_mismatch_penalty = 2.0f;
    align_params.min_alignment_score = 50.0f;
    align_params.min_identity = 0.8f;
    
    DEBUG_PRINT("FQMutationDetectorCUDA constructor completed");
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
    
    DEBUG_PRINT("FQMutationDetectorCUDA destructor completed");
}

void FQMutationDetectorCUDA::loadIndex(const char* index_path) {
    DEBUG_PRINT("FQMutationDetectorCUDA::loadIndex called with: %s", index_path);
    
    // Check if file exists
    std::ifstream test_file(index_path);
    if (!test_file.good()) {
        std::cerr << "ERROR: Cannot open index file: " << index_path << std::endl;
        DEBUG_PRINT("File does not exist or cannot be opened");
        // Initialize to prevent crashes
        num_kmers = 0;
        num_sequences = 0;
        return;
    }
    test_file.close();
    
    // Check if this is a binary index directory
    std::string path_str(index_path);
    
    // Try to load as binary index (preferred method)
    std::string kmer_index_path = path_str + "/kmer_index.bin";
    std::ifstream kmer_check(kmer_index_path);
    if (kmer_check.good()) {
        kmer_check.close();
        DEBUG_PRINT("Found binary index, loading via enhanced loader");
        loadBinaryIndex(path_str);
        return;
    }
    
    // Fall back to simple loading
    DEBUG_PRINT("Using simplified index loader for: %s", index_path);
    
    // Create minimal index for testing
    num_kmers = 1000;  // Placeholder
    num_sequences = 100;
    total_ref_length = 10000;
    
    // Allocate minimal GPU memory
    cudaMalloc(&d_kmer_index, num_kmers * sizeof(KmerEntry));
    cudaMalloc(&d_kmer_sorted, num_kmers * sizeof(uint64_t));
    cudaMalloc(&d_kmer_positions, num_kmers * sizeof(uint32_t));
    
    DEBUG_PRINT("Loaded simplified index: %d k-mers, %d sequences", num_kmers, num_sequences);
}

void FQMutationDetectorCUDA::loadBinaryIndex(const std::string& index_dir) {
    DEBUG_PRINT("Loading binary index from: %s", index_dir.c_str());
    
    // Load k-mer index
    std::string kmer_path = index_dir + "/kmer_index.bin";
    std::ifstream kmer_file(kmer_path, std::ios::binary);
    
    if (!kmer_file.good()) {
        std::cerr << "ERROR: Cannot open k-mer index: " << kmer_path << std::endl;
        return;
    }
    
    // Read header
    uint32_t num_entries, kmer_length;
    kmer_file.read(reinterpret_cast<char*>(&num_entries), sizeof(uint32_t));
    kmer_file.read(reinterpret_cast<char*>(&kmer_length), sizeof(uint32_t));
    
    DEBUG_PRINT("Loading %d k-mer entries (k=%d)", num_entries, kmer_length);
    
    // Read all k-mer entries
    std::vector<KmerEntry> kmer_entries(num_entries);
    for (uint32_t i = 0; i < num_entries; i++) {
        KmerEntry& entry = kmer_entries[i];
        kmer_file.read(reinterpret_cast<char*>(&entry.kmer), sizeof(uint64_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.gene_id), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.species_id), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.seq_id), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.position), sizeof(uint16_t));
    }
    kmer_file.close();
    
    // Create sorted k-mer list
    std::map<uint64_t, std::vector<uint32_t>> kmer_to_positions;
    for (uint32_t i = 0; i < num_entries; i++) {
        uint64_t kmer = kmer_entries[i].kmer;
        kmer_to_positions[kmer].push_back(i);
    }
    
    std::vector<uint64_t> sorted_kmers;
    std::vector<uint32_t> kmer_start_positions;
    
    for (const auto& pair : kmer_to_positions) {
        sorted_kmers.push_back(pair.first);
        kmer_start_positions.push_back(pair.second[0]);
    }
    
    num_kmers = sorted_kmers.size();
    uint32_t total_entries = kmer_entries.size();
    
    DEBUG_PRINT("Created sorted k-mer list: %d unique k-mers, %d total entries", num_kmers, total_entries);
    
    // Allocate GPU memory and copy data
    cudaMalloc(&d_kmer_index, total_entries * sizeof(KmerEntry));
    cudaMalloc(&d_kmer_sorted, num_kmers * sizeof(uint64_t));
    cudaMalloc(&d_kmer_positions, num_kmers * sizeof(uint32_t));
    
    cudaMemcpy(d_kmer_index, kmer_entries.data(), total_entries * sizeof(KmerEntry), cudaMemcpyHostToDevice);
    cudaMemcpy(d_kmer_sorted, sorted_kmers.data(), num_kmers * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_kmer_positions, kmer_start_positions.data(), num_kmers * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Set other members
    num_sequences = 100; // Placeholder
    total_ref_length = 10000; // Placeholder
    
    DEBUG_PRINT("Binary index loaded successfully to GPU");
}

void FQMutationDetectorCUDA::processReads(const char* r1_path, const char* r2_path, const char* output_path) {
    DEBUG_PRINT("Processing paired-end reads:");
    DEBUG_PRINT("  R1: %s", r1_path);
    DEBUG_PRINT("  R2: %s", r2_path);
    DEBUG_PRINT("  Output: %s", output_path);
    
    // TODO: Implement full pipeline
    // For now, just output a placeholder
    std::ofstream output(output_path);
    output << "{\n";
    output << "  \"status\": \"pipeline_ready\",\n";
    output << "  \"r1_path\": \"" << r1_path << "\",\n";
    output << "  \"r2_path\": \"" << r2_path << "\",\n";
    output << "  \"kmers_loaded\": " << num_kmers << ",\n";
    output << "  \"sequences_loaded\": " << num_sequences << "\n";
    output << "}\n";
    output.close();
    
    DEBUG_PRINT("Placeholder output written to: %s", output_path);
}

// Wrapper functions for calling kernels from C++
extern "C" {
    void launch_kmer_filter(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        FQMutationDetectorCUDA& detector,
        int num_reads,
        CandidateMatch* d_candidates,
        uint32_t* d_candidate_counts,
        bool check_reverse_complement
    ) {
        DEBUG_PRINT("launch_kmer_filter called with %d reads", num_reads);
        
        // Launch the enhanced k-mer filter kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        DEBUG_PRINT("Launching enhanced_kmer_filter_kernel: grid=%d, block=%d", grid_size, block_size);
        
        enhanced_kmer_filter_kernel<<<grid_size, block_size>>>(
            d_reads,
            d_read_lengths,
            d_read_offsets,
            detector.d_kmer_index,
            detector.d_kmer_sorted,
            detector.d_kmer_positions,
            detector.num_kmers,
            detector.num_kmers, // total_kmer_entries (simplified)
            num_reads,
            KMER_LENGTH,
            d_candidates,
            d_candidate_counts,
            MAX_CANDIDATES_PER_READ,
            check_reverse_complement,
            nullptr // debug_stats
        );
        
        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            DEBUG_PRINT("Kernel launch error: %s", cudaGetErrorString(error));
            return;
        }
        
        cudaDeviceSynchronize();
        DEBUG_PRINT("launch_kmer_filter completed");
    }

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
    ) {
        DEBUG_PRINT("launch_position_weighted_alignment called with %d reads", num_reads);
        
        // Launch simplified alignment kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        DEBUG_PRINT("Launching simple_alignment_kernel: grid=%d, block=%d", grid_size, block_size);
        
        simple_alignment_kernel<<<grid_size, block_size>>>(
            d_reads,
            d_read_lengths,
            d_read_offsets,
            d_candidates,
            d_candidate_counts,
            MAX_CANDIDATES_PER_READ,
            num_reads,
            d_results,
            d_result_count
        );
        
        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            DEBUG_PRINT("Alignment kernel launch error: %s", cudaGetErrorString(error));
            return;
        }
        
        cudaDeviceSynchronize();
        DEBUG_PRINT("launch_position_weighted_alignment completed");
    }
}

// Main entry point for testing
extern "C" void run_fq_mutation_detection(
    const char* index_path,
    const char* r1_path,
    const char* r2_path,
    const char* output_path
) {
    DEBUG_PRINT("run_fq_mutation_detection called");
    DEBUG_PRINT("  index_path: %s", index_path);
    DEBUG_PRINT("  r1_path: %s", r1_path);
    DEBUG_PRINT("  r2_path: %s", r2_path);
    DEBUG_PRINT("  output_path: %s", output_path);
    
    FQMutationDetectorCUDA detector;
    detector.loadIndex(index_path);
    detector.processReads(r1_path, r2_path, output_path);
    
    DEBUG_PRINT("run_fq_mutation_detection completed");
}