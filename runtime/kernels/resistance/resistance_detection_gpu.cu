// resistance_detection_gpu.cu
// GPU-accelerated fluoroquinolone resistance detection pipeline
// Combines pileup-based variant calling with direct resistance matching

#ifndef RESISTANCE_DETECTION_GPU_CU
#define RESISTANCE_DETECTION_GPU_CU

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>
#include <cub/cub.cuh>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>

#include <map>

namespace cg = cooperative_groups;

// ============================================================================
// Core Data Structures
// ============================================================================

// Constants
#define MAX_GENES 10
#define MAX_POSITIONS_PER_GENE 1000
#define MAX_AMINO_ACIDS 20
#define MAX_RESISTANCE_POSITIONS 100
#define MIN_DEPTH 5
#define MIN_ALLELE_FREQ 0.1f
#define MAX_READS_PER_BATCH 100000
#define MAX_ALIGNMENTS_PER_READ 10
#define MINIMIZER_K 8  // For protein minimizers
#define MINIMIZER_W 4  // Window size

// Amino acid encoding
__device__ __constant__ char AA_TABLE[20] = {
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
};

// Known resistance mutations database
struct ResistanceMutation {
    uint32_t gene_id;
    uint16_t position;      // 1-based position
    char wildtype_aa;
    char resistant_aas[10]; // Up to 10 known resistant variants
    uint8_t num_resistant;
    float resistance_score;  // Clinical significance (0-1)
    char drug_affected[32]; // Drug name
};

// Protein alignment result
struct ProteinAlignment {
    uint32_t read_id;
    uint32_t gene_id;
    uint16_t ref_start;
    uint16_t query_start;
    uint16_t alignment_length;
    float alignment_score;
    float identity;
    // Mutation details
    uint8_t num_mutations;
    uint16_t mutation_positions[20];
    char ref_aas[20];
    char query_aas[20];
};

// Variant pileup structure
struct VariantPileup {
    // Use atomic counters for thread-safe updates
    // [gene][position][amino_acid] = count
    uint32_t* counts;  // Flattened 3D array
    uint32_t* depths;  // Total depth per position
    
    __device__ uint32_t* get_count_ptr(int gene, int pos, int aa_idx) const {
        return &counts[(gene * MAX_POSITIONS_PER_GENE + pos) * MAX_AMINO_ACIDS + aa_idx];
    }
    
    __device__ uint32_t* get_depth_ptr(int gene, int pos) const {
        return &depths[gene * MAX_POSITIONS_PER_GENE + pos];
    }
};

// Resistance call result
struct ResistanceCall {
    uint32_t gene_id;
    uint16_t position;
    char wildtype_aa;
    char observed_aa;
    float allele_frequency;
    uint32_t supporting_reads;
    uint32_t total_depth;
    float confidence_score;
    char drug_affected[32];
    bool is_resistance_mutation;
};

// Minimizer index for fast seeding
struct MinimizerIndex {
    uint64_t* minimizer_hashes;
    uint32_t* position_lists;
    uint32_t* list_starts;
    uint32_t* list_lengths;
    uint32_t num_minimizers;
};

// ============================================================================
// Device Functions
// ============================================================================

// Convert amino acid to index (host version)
inline int aa_to_index_host(char aa) {
    switch(aa) {
        case 'A': return 0;  case 'C': return 1;  case 'D': return 2;  case 'E': return 3;
        case 'F': return 4;  case 'G': return 5;  case 'H': return 6;  case 'I': return 7;
        case 'K': return 8;  case 'L': return 9;  case 'M': return 10; case 'N': return 11;
        case 'P': return 12; case 'Q': return 13; case 'R': return 14; case 'S': return 15;
        case 'T': return 16; case 'V': return 17; case 'W': return 18; case 'Y': return 19;
        default: return -1;
    }
}

// Convert amino acid to index (device version)
__device__ inline int aa_to_index(char aa) {
    switch(aa) {
        case 'A': return 0;  case 'C': return 1;  case 'D': return 2;  case 'E': return 3;
        case 'F': return 4;  case 'G': return 5;  case 'H': return 6;  case 'I': return 7;
        case 'K': return 8;  case 'L': return 9;  case 'M': return 10; case 'N': return 11;
        case 'P': return 12; case 'Q': return 13; case 'R': return 14; case 'S': return 15;
        case 'T': return 16; case 'V': return 17; case 'W': return 18; case 'Y': return 19;
        default: return -1;
    }
}

// Hash function for protein minimizers
__device__ inline uint64_t hash_protein_kmer(const char* seq, int pos, int k) {
    uint64_t hash = 0;
    const uint64_t prime = 31;
    
    for (int i = 0; i < k; i++) {
        int aa_idx = aa_to_index(seq[pos + i]);
        if (aa_idx < 0) return UINT64_MAX;  // Invalid
        hash = hash * prime + aa_idx;
    }
    
    return hash;
}

// Extract minimizers from a sequence
__device__ void extract_minimizers(
    const char* seq,
    int seq_len,
    uint64_t* minimizers,
    int* positions,
    int* num_minimizers,
    int k = MINIMIZER_K,
    int w = MINIMIZER_W
) {
    *num_minimizers = 0;
    
    for (int i = 0; i <= seq_len - k - w + 1; i++) {
        uint64_t min_hash = UINT64_MAX;
        int min_pos = -1;
        
        // Find minimum k-mer in window
        for (int j = 0; j < w; j++) {
            uint64_t hash = hash_protein_kmer(seq, i + j, k);
            if (hash < min_hash) {
                min_hash = hash;
                min_pos = i + j;
            }
        }
        
        if (min_pos >= 0 && (*num_minimizers == 0 || 
            minimizers[*num_minimizers - 1] != min_hash)) {
            minimizers[*num_minimizers] = min_hash;
            positions[*num_minimizers] = min_pos;
            (*num_minimizers)++;
            
            if (*num_minimizers >= 100) break;  // Limit
        }
    }
}

// Binary search in minimizer index
__device__ int find_minimizer_matches(
    uint64_t minimizer,
    const MinimizerIndex* index,
    uint32_t* match_positions,
    int max_matches
) {
    // Binary search for minimizer
    int left = 0;
    int right = index->num_minimizers - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        if (index->minimizer_hashes[mid] == minimizer) {
            // Found - get all positions
            uint32_t start = index->list_starts[mid];
            uint32_t length = index->list_lengths[mid];
            
            int num_copied = min((int)length, max_matches);
            for (int i = 0; i < num_copied; i++) {
                match_positions[i] = index->position_lists[start + i];
            }
            
            return num_copied;
        } else if (index->minimizer_hashes[mid] < minimizer) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    return 0;  // Not found
}

// Banded Smith-Waterman alignment
__device__ float banded_smith_waterman(
    const char* query,
    const char* ref,
    int query_len,
    int ref_len,
    int band_width,
    ProteinAlignment* result
) {
    // Simplified banded alignment for GPU
    const float MATCH = 5.0f;
    const float MISMATCH = -3.0f;
    const float GAP_OPEN = -8.0f;
    const float GAP_EXTEND = -2.0f;
    
    // Use shared memory for DP matrix
    extern __shared__ float dp_matrix[];
    
    float max_score = 0.0f;
    int max_i = 0, max_j = 0;
    
    // Initialize first row and column
    for (int i = 0; i <= band_width && i <= query_len; i++) {
        dp_matrix[i * (2 * band_width + 1)] = 0;
    }
    
    // Fill DP matrix with banding
    for (int i = 1; i <= query_len; i++) {
        int j_start = max(1, i - band_width);
        int j_end = min(ref_len, i + band_width);
        
        for (int j = j_start; j <= j_end; j++) {
            int band_j = j - i + band_width;
            
            float match_score = (query[i-1] == ref[j-1]) ? MATCH : MISMATCH;
            
            float diag = dp_matrix[(i-1) * (2 * band_width + 1) + band_j] + match_score;
            float up = dp_matrix[(i-1) * (2 * band_width + 1) + band_j + 1] + GAP_OPEN;
            float left = (band_j > 0) ? 
                dp_matrix[i * (2 * band_width + 1) + band_j - 1] + GAP_OPEN : 0;
            
            float score = fmaxf(0.0f, fmaxf(diag, fmaxf(up, left)));
            dp_matrix[i * (2 * band_width + 1) + band_j] = score;
            
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback to get alignment details
    result->alignment_score = max_score;
    result->num_mutations = 0;
    
    // Simplified traceback - just count mutations in aligned region
    int i = max_i, j = max_j;
    int align_len = 0;
    int matches = 0;
    
    while (i > 0 && j > 0 && align_len < 20) {
        if (query[i-1] != ref[j-1]) {
            result->mutation_positions[result->num_mutations] = j - 1;
            result->ref_aas[result->num_mutations] = ref[j-1];
            result->query_aas[result->num_mutations] = query[i-1];
            result->num_mutations++;
        } else {
            matches++;
        }
        i--; j--;
        align_len++;
    }
    
    result->ref_start = j;
    result->query_start = i;
    result->alignment_length = align_len;
    result->identity = (float)matches / align_len;
    
    return max_score;
}

// ============================================================================
// Main GPU Kernels
// ============================================================================

// Kernel 1: Seed and align proteins
__global__ void seed_and_align_kernel(
    const char* translated_reads,
    const int* read_lengths,
    const int* read_offsets,
    const char* reference_proteins,
    const int* ref_lengths,
    const int* ref_offsets,
    const MinimizerIndex* minimizer_index,
    const int num_reads,
    const int num_refs,
    ProteinAlignment* alignments,
    uint32_t* alignment_counts,
    const float min_score_threshold
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const char* read = translated_reads + read_offsets[tid];
    int read_len = read_lengths[tid];
    
    // Extract minimizers from read
    uint64_t read_minimizers[100];
    int minimizer_positions[100];
    int num_minimizers = 0;
    
    extract_minimizers(read, read_len, read_minimizers, minimizer_positions, 
                      &num_minimizers);
    
    // Find matches for each minimizer
    uint32_t match_positions[50];
    ProteinAlignment best_alignments[MAX_ALIGNMENTS_PER_READ];
    int num_alignments = 0;
    
    for (int i = 0; i < num_minimizers && num_alignments < MAX_ALIGNMENTS_PER_READ; i++) {
        int num_matches = find_minimizer_matches(
            read_minimizers[i], minimizer_index, match_positions, 50
        );
        
        // Try alignment at each match position
        for (int j = 0; j < num_matches; j++) {
            uint32_t ref_idx = match_positions[j] >> 16;
            uint32_t ref_pos = match_positions[j] & 0xFFFF;
            
            if (ref_idx >= num_refs) continue;
            
            const char* ref = reference_proteins + ref_offsets[ref_idx];
            int ref_len = ref_lengths[ref_idx];
            
            // Perform banded alignment
            ProteinAlignment temp_align;
            temp_align.read_id = tid;
            temp_align.gene_id = ref_idx;
            
            float score = banded_smith_waterman(
                read + minimizer_positions[i],
                ref + ref_pos,
                read_len - minimizer_positions[i],
                ref_len - ref_pos,
                20,  // band width
                &temp_align
            );
            
            if (score >= min_score_threshold) {
                best_alignments[num_alignments++] = temp_align;
            }
        }
    }
    
    // Store results
    alignment_counts[tid] = num_alignments;
    for (int i = 0; i < num_alignments; i++) {
        alignments[tid * MAX_ALIGNMENTS_PER_READ + i] = best_alignments[i];
    }
}

// Kernel 2: Build variant pileup from alignments
__global__ void build_pileup_kernel(
    const ProteinAlignment* alignments,
    const uint32_t* alignment_counts,
    const int num_reads,
    VariantPileup* pileup,
    const ResistanceMutation* known_mutations,
    const int num_known_mutations
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    uint32_t num_aligns = alignment_counts[tid];
    
    for (uint32_t i = 0; i < num_aligns; i++) {
        const ProteinAlignment& align = alignments[tid * MAX_ALIGNMENTS_PER_READ + i];
        
        // Process each mutation in the alignment
        for (int j = 0; j < align.num_mutations; j++) {
            uint16_t ref_pos = align.ref_start + align.mutation_positions[j];
            char query_aa = align.query_aas[j];
            int aa_idx = aa_to_index(query_aa);
            
            if (aa_idx >= 0) {
                // Atomic increment for this amino acid at this position
                atomicAdd(pileup->get_count_ptr(align.gene_id, ref_pos, aa_idx), 1);
                atomicAdd(pileup->get_depth_ptr(align.gene_id, ref_pos), 1);
            }
        }
        
        // Also count reference amino acids where there's no mutation
        for (int pos = 0; pos < align.alignment_length; pos++) {
            bool is_mutation_pos = false;
            for (int j = 0; j < align.num_mutations; j++) {
                if (align.mutation_positions[j] == pos) {
                    is_mutation_pos = true;
                    break;
                }
            }
            
            if (!is_mutation_pos) {
                // This position matches reference
                uint16_t ref_pos = align.ref_start + pos;
                
                // We'd need the reference sequence here to know which AA
                // For now, increment a "reference" counter
                atomicAdd(pileup->get_depth_ptr(align.gene_id, ref_pos), 1);
            }
        }
    }
}

// Kernel 3: Call variants from pileup
__global__ void call_variants_kernel(
    const VariantPileup* pileup,
    const ResistanceMutation* known_mutations,
    const int num_known_mutations,
    ResistanceCall* calls,
    uint32_t* num_calls,
    const float min_depth,
    const float min_allele_freq
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_known_mutations) return;
    
    const ResistanceMutation& mut = known_mutations[tid];
    
    // Get depth at this position
    uint32_t total_depth = *pileup->get_depth_ptr(mut.gene_id, mut.position - 1);
    
    if (total_depth < min_depth) return;
    
    // Check each possible resistant amino acid
    for (int i = 0; i < mut.num_resistant; i++) {
        char resistant_aa = mut.resistant_aas[i];
        int aa_idx = aa_to_index(resistant_aa);
        
        if (aa_idx < 0) continue;
        
        uint32_t alt_count = *pileup->get_count_ptr(mut.gene_id, mut.position - 1, aa_idx);
        float allele_freq = (float)alt_count / total_depth;
        
        if (allele_freq >= min_allele_freq) {
            // Add resistance call
            uint32_t call_idx = atomicAdd(num_calls, 1);
            
            ResistanceCall& call = calls[call_idx];
            call.gene_id = mut.gene_id;
            call.position = mut.position;
            call.wildtype_aa = mut.wildtype_aa;
            call.observed_aa = resistant_aa;
            call.allele_frequency = allele_freq;
            call.supporting_reads = alt_count;
            call.total_depth = total_depth;
            call.is_resistance_mutation = true;
            
            // Calculate confidence score
            float z = 1.96f;  // 95% confidence
            float p = allele_freq;
            float n = (float)total_depth;
            float stderr = sqrtf(p * (1.0f - p) / n);
            float lower_bound = fmaxf(0.0f, p - z * stderr);
            
            call.confidence_score = lower_bound * (1.0f - expf(-n / 100.0f));
            
            // Copy drug name
            for (int j = 0; j < 32; j++) {
                call.drug_affected[j] = mut.drug_affected[j];
            }
        }
    }
}

// ============================================================================
// Host Functions
// ============================================================================

class ResistanceDetectorGPU {
private:
    // Device memory
    char* d_translated_reads;
    int* d_read_lengths;
    int* d_read_offsets;
    
    char* d_reference_proteins;
    int* d_ref_lengths;
    int* d_ref_offsets;
    
    MinimizerIndex* d_minimizer_index;
    ResistanceMutation* d_known_mutations;
    
    ProteinAlignment* d_alignments;
    uint32_t* d_alignment_counts;
    
    VariantPileup* d_pileup;
    ResistanceCall* d_calls;
    uint32_t* d_num_calls;
    
    // Parameters
    int max_reads;
    int num_genes;
    int num_known_mutations;
    
public:
    ResistanceDetectorGPU(int max_reads_per_batch, int num_genes, int num_mutations) 
        : max_reads(max_reads_per_batch), num_genes(num_genes), 
          num_known_mutations(num_mutations) {
        
        // Allocate device memory
        size_t read_buffer_size = max_reads * 300;  // Max 300 AA per read
        cudaMalloc(&d_translated_reads, read_buffer_size);
        cudaMalloc(&d_read_lengths, max_reads * sizeof(int));
        cudaMalloc(&d_read_offsets, max_reads * sizeof(int));
        
        // Allocate alignment memory
        cudaMalloc(&d_alignments, max_reads * MAX_ALIGNMENTS_PER_READ * sizeof(ProteinAlignment));
        cudaMalloc(&d_alignment_counts, max_reads * sizeof(uint32_t));
        
        // Allocate pileup memory
        size_t pileup_size = num_genes * MAX_POSITIONS_PER_GENE * MAX_AMINO_ACIDS * sizeof(uint32_t);
        cudaMalloc(&d_pileup, sizeof(VariantPileup));
        
        VariantPileup h_pileup;
        cudaMalloc(&h_pileup.counts, pileup_size);
        cudaMalloc(&h_pileup.depths, num_genes * MAX_POSITIONS_PER_GENE * sizeof(uint32_t));
        cudaMemset(h_pileup.counts, 0, pileup_size);
        cudaMemset(h_pileup.depths, 0, num_genes * MAX_POSITIONS_PER_GENE * sizeof(uint32_t));
        
        cudaMemcpy(d_pileup, &h_pileup, sizeof(VariantPileup), cudaMemcpyHostToDevice);
        
        // Allocate results memory
        cudaMalloc(&d_calls, MAX_RESISTANCE_POSITIONS * 10 * sizeof(ResistanceCall));
        cudaMalloc(&d_num_calls, sizeof(uint32_t));
        
        d_reference_proteins = nullptr;
        d_minimizer_index = nullptr;
        d_known_mutations = nullptr;
    }
    
    ~ResistanceDetectorGPU() {
        cudaFree(d_translated_reads);
        cudaFree(d_read_lengths);
        cudaFree(d_read_offsets);
        cudaFree(d_alignments);
        cudaFree(d_alignment_counts);
        
        if (d_reference_proteins) cudaFree(d_reference_proteins);
        if (d_minimizer_index) {
            MinimizerIndex h_index;
            cudaMemcpy(&h_index, d_minimizer_index, sizeof(MinimizerIndex), cudaMemcpyDeviceToHost);
            cudaFree(h_index.minimizer_hashes);
            cudaFree(h_index.position_lists);
            cudaFree(h_index.list_starts);
            cudaFree(h_index.list_lengths);
            cudaFree(d_minimizer_index);
        }
        if (d_known_mutations) cudaFree(d_known_mutations);
        
        VariantPileup h_pileup;
        cudaMemcpy(&h_pileup, d_pileup, sizeof(VariantPileup), cudaMemcpyDeviceToHost);
        cudaFree(h_pileup.counts);
        cudaFree(h_pileup.depths);
        cudaFree(d_pileup);
        
        cudaFree(d_calls);
        cudaFree(d_num_calls);
    }
    
    void loadReferences(const std::vector<std::string>& sequences) {
        // Implementation would load reference sequences and build minimizer index
        // For brevity, showing the structure
        
        // Concatenate sequences
        std::string concat_refs;
        std::vector<int> lengths;
        std::vector<int> offsets;
        
        for (const auto& seq : sequences) {
            offsets.push_back(concat_refs.length());
            lengths.push_back(seq.length());
            concat_refs += seq;
        }
        
        // Copy to GPU
        cudaMalloc(&d_reference_proteins, concat_refs.length());
        cudaMalloc(&d_ref_lengths, lengths.size() * sizeof(int));
        cudaMalloc(&d_ref_offsets, offsets.size() * sizeof(int));
        
        cudaMemcpy(d_reference_proteins, concat_refs.c_str(), concat_refs.length(), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_ref_lengths, lengths.data(), lengths.size() * sizeof(int), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_ref_offsets, offsets.data(), offsets.size() * sizeof(int), 
                   cudaMemcpyHostToDevice);
        
        // Build minimizer index (simplified)
        buildMinimizerIndex(sequences);
    }
    
    void loadResistanceDatabase(const std::vector<ResistanceMutation>& mutations) {
        cudaMalloc(&d_known_mutations, mutations.size() * sizeof(ResistanceMutation));
        cudaMemcpy(d_known_mutations, mutations.data(), 
                   mutations.size() * sizeof(ResistanceMutation), 
                   cudaMemcpyHostToDevice);
    }
    
    std::vector<ResistanceCall> detectResistance(
        const std::vector<std::string>& translated_reads,
        float min_score = 50.0f,
        float min_depth = 5.0f,
        float min_allele_freq = 0.1f
    ) {
        // Prepare read data
        std::string concat_reads;
        std::vector<int> lengths;
        std::vector<int> offsets;
        
        for (const auto& read : translated_reads) {
            offsets.push_back(concat_reads.length());
            lengths.push_back(read.length());
            concat_reads += read;
        }
        
        int num_reads = translated_reads.size();
        
        // Copy to GPU
        cudaMemcpy(d_translated_reads, concat_reads.c_str(), concat_reads.length(), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, lengths.data(), num_reads * sizeof(int), 
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, offsets.data(), num_reads * sizeof(int), 
                   cudaMemcpyHostToDevice);
        
        // Reset counters
        cudaMemset(d_alignment_counts, 0, num_reads * sizeof(uint32_t));
        cudaMemset(d_num_calls, 0, sizeof(uint32_t));
        
        // Clear pileup
        VariantPileup h_pileup;
        cudaMemcpy(&h_pileup, d_pileup, sizeof(VariantPileup), cudaMemcpyDeviceToHost);
        size_t pileup_size = num_genes * MAX_POSITIONS_PER_GENE * MAX_AMINO_ACIDS * sizeof(uint32_t);
        cudaMemset(h_pileup.counts, 0, pileup_size);
        cudaMemset(h_pileup.depths, 0, num_genes * MAX_POSITIONS_PER_GENE * sizeof(uint32_t));
        
        // Launch kernels
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        // Stage 1: Seed and align
        size_t shared_mem_size = 100 * 100 * sizeof(float);  // For DP matrix
        seed_and_align_kernel<<<grid_size, block_size, shared_mem_size>>>(
            d_translated_reads, d_read_lengths, d_read_offsets,
            d_reference_proteins, d_ref_lengths, d_ref_offsets,
            d_minimizer_index, num_reads, num_genes,
            d_alignments, d_alignment_counts, min_score
        );
        
        cudaDeviceSynchronize();
        
        // Stage 2: Build pileup
        build_pileup_kernel<<<grid_size, block_size>>>(
            d_alignments, d_alignment_counts, num_reads,
            d_pileup, d_known_mutations, num_known_mutations
        );
        
        cudaDeviceSynchronize();
        
        // Stage 3: Call variants
        int var_grid_size = (num_known_mutations + block_size - 1) / block_size;
        call_variants_kernel<<<var_grid_size, block_size>>>(
            d_pileup, d_known_mutations, num_known_mutations,
            d_calls, d_num_calls, min_depth, min_allele_freq
        );
        
        cudaDeviceSynchronize();
        
        // Copy results back
        uint32_t num_calls;
        cudaMemcpy(&num_calls, d_num_calls, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::vector<ResistanceCall> results(num_calls);
        if (num_calls > 0) {
            cudaMemcpy(results.data(), d_calls, num_calls * sizeof(ResistanceCall), 
                       cudaMemcpyDeviceToHost);
        }
        
        return results;
    }
    
private:
    void buildMinimizerIndex(const std::vector<std::string>& sequences) {
        // Build minimizer index on host then copy to device
        // This is a simplified version - real implementation would be more sophisticated
        
        std::map<uint64_t, std::vector<uint32_t>> minimizer_map;
        
        for (size_t seq_idx = 0; seq_idx < sequences.size(); seq_idx++) {
            const std::string& seq = sequences[seq_idx];
            
            for (size_t i = 0; i <= seq.length() - MINIMIZER_K - MINIMIZER_W + 1; i++) {
                uint64_t min_hash = UINT64_MAX;
                int min_pos = -1;
                
                // Find minimizer in window
                for (int j = 0; j < MINIMIZER_W; j++) {
                    uint64_t hash = 0;
                    const uint64_t prime = 31;
                    
                    for (int k = 0; k < MINIMIZER_K; k++) {
                        int aa_idx = aa_to_index_host(seq[i + j + k]);
                        if (aa_idx < 0) break;
                        hash = hash * prime + aa_idx;
                    }
                    
                    if (hash < min_hash) {
                        min_hash = hash;
                        min_pos = i + j;
                    }
                }
                
                if (min_pos >= 0) {
                    uint32_t encoded = (seq_idx << 16) | (min_pos & 0xFFFF);
                    minimizer_map[min_hash].push_back(encoded);
                }
            }
        }
        
        // Convert to arrays
        std::vector<uint64_t> hashes;
        std::vector<uint32_t> all_positions;
        std::vector<uint32_t> starts;
        std::vector<uint32_t> lengths;
        
        for (const auto& pair : minimizer_map) {
            hashes.push_back(pair.first);
            starts.push_back(all_positions.size());
            lengths.push_back(pair.second.size());
            
            for (uint32_t pos : pair.second) {
                all_positions.push_back(pos);
            }
        }
        
        // Allocate and copy to GPU
        MinimizerIndex h_index;
        h_index.num_minimizers = hashes.size();
        
        cudaMalloc(&h_index.minimizer_hashes, hashes.size() * sizeof(uint64_t));
        cudaMalloc(&h_index.position_lists, all_positions.size() * sizeof(uint32_t));
        cudaMalloc(&h_index.list_starts, starts.size() * sizeof(uint32_t));
        cudaMalloc(&h_index.list_lengths, lengths.size() * sizeof(uint32_t));
        
        cudaMemcpy(h_index.minimizer_hashes, hashes.data(), 
                   hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_index.position_lists, all_positions.data(), 
                   all_positions.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_index.list_starts, starts.data(), 
                   starts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_index.list_lengths, lengths.data(), 
                   lengths.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        cudaMalloc(&d_minimizer_index, sizeof(MinimizerIndex));
        cudaMemcpy(d_minimizer_index, &h_index, sizeof(MinimizerIndex), cudaMemcpyHostToDevice);
    }
};

// C interface for integration
extern "C" {
    void* create_resistance_detector_gpu(int max_reads, int num_genes, int num_mutations) {
        return new ResistanceDetectorGPU(max_reads, num_genes, num_mutations);
    }
    
    void destroy_resistance_detector_gpu(void* detector) {
        delete static_cast<ResistanceDetectorGPU*>(detector);
    }
    
    int detect_resistance_gpu(
        void* detector,
        const char** translated_reads,
        int* read_lengths,
        int num_reads,
        void* results,
        float min_score,
        float min_depth,
        float min_allele_freq
    ) {
        ResistanceDetectorGPU* det = static_cast<ResistanceDetectorGPU*>(detector);
        
        std::vector<std::string> reads;
        for (int i = 0; i < num_reads; i++) {
            reads.emplace_back(translated_reads[i], read_lengths[i]);
        }
        
        auto calls = det->detectResistance(reads, min_score, min_depth, min_allele_freq);
        
        memcpy(results, calls.data(), calls.size() * sizeof(ResistanceCall));
        return calls.size();
    }
}

#endif // RESISTANCE_DETECTION_GPU_CU