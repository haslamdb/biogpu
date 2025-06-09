// fq_resistance_kernels.cu
// CUDA kernels for fluoroquinolone resistance detection
// Optimized for NVIDIA Titan Xp (Compute Capability 6.1)

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

// Device-compatible string functions
__device__ int device_strcmp(const char* s1, const char* s2) {
    while (*s1 && (*s1 == *s2)) {
        s1++;
        s2++;
    }
    return *(const unsigned char*)s1 - *(const unsigned char*)s2;
}

__device__ int device_strncmp(const char* s1, const char* s2, size_t n) {
    while (n && *s1 && (*s1 == *s2)) {
        s1++;
        s2++;
        n--;
    }
    if (n == 0) return 0;
    return *(const unsigned char*)s1 - *(const unsigned char*)s2;
}

__device__ void device_strcpy(char* dst, const char* src) {
    while ((*dst++ = *src++) != '\0');
}

// Constants for genetic code and QRDR mutations
__constant__ char GENETIC_CODE[64][4] = {
    "F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "*", "*", "C", "C", "*", "W",
    "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R",
    "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R",
    "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G"
};

// QRDR mutation database
struct QRDRMutation {
    char gene[8];           // e.g., "gyrA"
    int position;           // codon position
    char wild_type;         // wild type amino acid
    char mutant;            // resistant amino acid
    float resistance_score; // clinical significance (0-1)
};

__constant__ QRDRMutation QRDR_MUTATIONS[32] = {
    // E. coli mutations
    {"gyrA", 83, 'S', 'L', 0.9f},
    {"gyrA", 83, 'S', 'F', 0.95f},
    {"gyrA", 87, 'D', 'N', 0.7f},
    {"parC", 80, 'S', 'I', 0.6f},
    {"parC", 84, 'E', 'V', 0.4f},
    
    // P. aeruginosa mutations
    {"gyrA", 83, 'T', 'I', 0.9f},
    {"parC", 87, 'S', 'L', 0.7f},
    
    // S. aureus mutations (grlA = parC homolog)
    {"grlA", 80, 'S', 'F', 0.9f},
    
    // C. difficile mutations
    {"gyrA", 82, 'T', 'I', 0.95f},
    
    // Null terminator
    {"", 0, 0, 0, 0.0f}
};

// Result structure for each read
struct MutationResult {
    int read_id;
    char gene[8];
    int position;
    char wild_type;
    char mutant;
    float confidence;
    int read_position;  // where in read the mutation was found
};

// DNA sequence data structure
struct SequenceBatch {
    char* sequences;        // packed sequences
    int* seq_offsets;      // start position of each sequence
    int* seq_lengths;      // length of each sequence
    int num_sequences;
    int max_length;
};

// Convert DNA triplet to amino acid index
__device__ int codon_to_index(const char* codon) {
    int index = 0;
    for (int i = 0; i < 3; i++) {
        switch (codon[i]) {
            case 'A': case 'a': index += 0; break;
            case 'C': case 'c': index += 1; break;
            case 'G': case 'g': index += 2; break;
            case 'T': case 't': case 'U': case 'u': index += 3; break;
            default: return -1; // Invalid nucleotide
        }
        if (i < 2) index *= 4;
    }
    return index;
}

// Translate DNA codon to amino acid
__device__ char translate_codon(const char* codon) {
    int index = codon_to_index(codon);
    if (index < 0 || index >= 64) return 'X'; // Unknown
    return GENETIC_CODE[index][0];
}

// Check if a sequence matches a known gene pattern
__device__ bool matches_gene_pattern(const char* sequence, int length, const char* gene_name) {
    // Simplified gene matching - in production would use more sophisticated alignment
    
    if (device_strcmp(gene_name, "gyrA") == 0) {
        // Look for gyrA-specific k-mers around QRDR region
        const char* gyrA_pattern = "ATGGAT"; // Simplified pattern
        for (int i = 0; i <= length - 6; i++) {
            if (device_strncmp(sequence + i, gyrA_pattern, 6) == 0) {
                return true;
            }
        }
    } else if (device_strcmp(gene_name, "parC") == 0) {
        const char* parC_pattern = "ATGGAG"; // Simplified pattern
        for (int i = 0; i <= length - 6; i++) {
            if (device_strncmp(sequence + i, parC_pattern, 6) == 0) {
                return true;
            }
        }
    }
    
    return false;
}

// Main kernel for mutation detection
__global__ void detect_mutations_kernel(
    SequenceBatch sequences,
    MutationResult* results,
    int* result_count,
    int max_results
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    auto grid = cg::this_grid();
    
    // Shared memory for storing temporary results
    __shared__ MutationResult shared_results[256];
    __shared__ int shared_count;
    
    if (threadIdx.x == 0) {
        shared_count = 0;
    }
    __syncthreads();
    
    if (tid >= sequences.num_sequences) return;
    
    // Get sequence for this thread
    int seq_offset = sequences.seq_offsets[tid];
    int seq_length = sequences.seq_lengths[tid];
    const char* sequence = sequences.sequences + seq_offset;
    
    // Check each mutation in the database
    for (int mut_idx = 0; QRDR_MUTATIONS[mut_idx].position != 0; mut_idx++) {
        const QRDRMutation& mutation = QRDR_MUTATIONS[mut_idx];
        
        // Skip if sequence doesn't match the gene pattern
        if (!matches_gene_pattern(sequence, seq_length, mutation.gene)) {
            continue;
        }
        
        // Scan sequence for the mutation position
        // This is simplified - real implementation would do proper alignment
        for (int pos = 0; pos <= seq_length - 3; pos += 3) { // Step by codons
            if (pos + 3 > seq_length) break;
            
            // Translate codon
            char amino_acid = translate_codon(sequence + pos);
            
            // Check if this could be the mutation position
            // (In real implementation, would map to exact genomic coordinates)
            if (amino_acid == mutation.mutant) {
                // Found potential mutation
                int local_idx = atomicAdd(&shared_count, 1);
                
                if (local_idx < 256) { // Prevent overflow
                    MutationResult& result = shared_results[local_idx];
                    result.read_id = tid;
                    device_strcpy(result.gene, mutation.gene);
                    result.position = mutation.position;
                    result.wild_type = mutation.wild_type;
                    result.mutant = mutation.mutant;
                    result.confidence = mutation.resistance_score * 0.9f; // Adjust for read quality
                    result.read_position = pos;
                }
            }
        }
    }
    
    __syncthreads();
    
    // Copy shared results to global memory
    if (threadIdx.x == 0 && shared_count > 0) {
        int global_start = atomicAdd(result_count, shared_count);
        
        for (int i = 0; i < shared_count && global_start + i < max_results; i++) {
            results[global_start + i] = shared_results[i];
        }
    }
}

// Kernel for k-mer based gene classification
__global__ void classify_reads_kernel(
    SequenceBatch sequences,
    int* gene_assignments,  // Output: which gene each read maps to
    float* confidence_scores
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (tid >= sequences.num_sequences) return;
    
    int seq_offset = sequences.seq_offsets[tid];
    int seq_length = sequences.seq_lengths[tid];
    const char* sequence = sequences.sequences + seq_offset;
    
    // Initialize
    gene_assignments[tid] = -1; // No assignment
    confidence_scores[tid] = 0.0f;
    
    // Check against known gene patterns
    float best_score = 0.0f;
    int best_gene = -1;
    
    // Gene 0: gyrA
    if (matches_gene_pattern(sequence, seq_length, "gyrA")) {
        float score = 0.8f; // Base confidence
        if (score > best_score) {
            best_score = score;
            best_gene = 0;
        }
    }
    
    // Gene 1: parC
    if (matches_gene_pattern(sequence, seq_length, "parC")) {
        float score = 0.8f;
        if (score > best_score) {
            best_score = score;
            best_gene = 1;
        }
    }
    
    // Gene 2: grlA (S. aureus)
    if (matches_gene_pattern(sequence, seq_length, "grlA")) {
        float score = 0.8f;
        if (score > best_score) {
            best_score = score;
            best_gene = 2;
        }
    }
    
    gene_assignments[tid] = best_gene;
    confidence_scores[tid] = best_score;
}

// Host function to launch mutation detection
extern "C" {
    __attribute__((visibility("default")))
    int launch_mutation_detection(
        char* d_sequences,
        int* d_seq_offsets,
        int* d_seq_lengths,
        int num_sequences,
        void* d_results,
        int max_results
    ) {
        // Create sequence batch structure
        SequenceBatch sequences;
        sequences.sequences = d_sequences;
        sequences.seq_offsets = d_seq_offsets;
        sequences.seq_lengths = d_seq_lengths;
        sequences.num_sequences = num_sequences;
        
        // Allocate result counter
        int* d_result_count;
        cudaMalloc(&d_result_count, sizeof(int));
        cudaMemset(d_result_count, 0, sizeof(int));
        
        // Launch kernel
        int blockSize = 256;
        int gridSize = (num_sequences + blockSize - 1) / blockSize;
        
        detect_mutations_kernel<<<gridSize, blockSize>>>(
            sequences,
            (MutationResult*)d_results,
            d_result_count,
            max_results
        );
        
        // Get result count
        int result_count;
        cudaMemcpy(&result_count, d_result_count, sizeof(int), cudaMemcpyDeviceToHost);
        
        // Cleanup
        cudaFree(d_result_count);
        
        return result_count;
    }
    
    __attribute__((visibility("default")))
    int launch_read_classification(
        char* d_sequences,
        int* d_seq_offsets,
        int* d_seq_lengths,
        int num_sequences,
        int* d_gene_assignments,
        float* d_confidence_scores
    ) {
        // Create sequence batch structure
        SequenceBatch sequences;
        sequences.sequences = d_sequences;
        sequences.seq_offsets = d_seq_offsets;
        sequences.seq_lengths = d_seq_lengths;
        sequences.num_sequences = num_sequences;
        
        // Launch kernel
        int blockSize = 256;
        int gridSize = (num_sequences + blockSize - 1) / blockSize;
        
        classify_reads_kernel<<<gridSize, blockSize>>>(
            sequences,
            d_gene_assignments,
            d_confidence_scores
        );
        
        // Check for errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            return -1;
        }
        
        cudaDeviceSynchronize();
        return 0;
    }
} // extern "C"