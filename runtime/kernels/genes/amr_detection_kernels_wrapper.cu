#include "amr_detection_kernels.h"
#include <cuda_runtime.h>

// Genetic code table initialization
extern "C" void initializeGeneticCode() {
    // Standard genetic code
    const char h_genetic_code[64] = {
        'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',  // AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT
        'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',  // AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT
        'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',  // CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT
        'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',  // CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT
        'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',  // GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT
        'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',  // GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT
        '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S',  // TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT
        'W', 'C', '*', 'C', 'L', 'F', 'L', 'F'   // TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT
    };
    
    cudaMemcpyToSymbol(GENETIC_CODE, h_genetic_code, sizeof(h_genetic_code));
}

// Kernel launch wrappers
extern "C" void launch_build_bloom_filter_kernel(
    const char* amr_sequences,
    const uint32_t* sequence_offsets,
    const uint32_t* sequence_lengths,
    uint64_t* bloom_filter,
    const uint32_t num_sequences,
    const uint32_t bloom_size,
    const int k
) {
    dim3 blockSize(256);
    dim3 gridSize((num_sequences + blockSize.x - 1) / blockSize.x);
    
    build_bloom_filter_kernel<<<gridSize, blockSize>>>(
        amr_sequences, sequence_offsets, sequence_lengths,
        bloom_filter, num_sequences, bloom_size, k
    );
}

extern "C" void launch_generate_minimizers_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    Minimizer* minimizers,
    uint32_t* minimizer_counts,
    uint32_t* minimizer_offsets,
    const int num_reads,
    const int k,
    const int w
) {
    dim3 blockSize(256);
    dim3 gridSize((num_reads + blockSize.x - 1) / blockSize.x);
    
    generate_minimizers_kernel<<<gridSize, blockSize>>>(
        reads, read_offsets, read_lengths,
        minimizers, minimizer_counts, minimizer_offsets,
        num_reads, k, w
    );
}

extern "C" void launch_screen_minimizers_kernel(
    const Minimizer* minimizers,
    const uint32_t* minimizer_counts,
    const uint32_t* minimizer_offsets,
    const uint64_t* bloom_filter,
    const int bloom_size,
    bool* read_passes,
    const int num_reads
) {
    dim3 blockSize(256);
    dim3 gridSize((num_reads + blockSize.x - 1) / blockSize.x);
    
    screen_minimizers_kernel<<<gridSize, blockSize>>>(
        minimizers, minimizer_counts, minimizer_offsets,
        bloom_filter, read_passes, num_reads, bloom_size
    );
}

extern "C" void launch_translated_alignment_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    const bool* read_passes_filter,
    const char* amr_proteins,
    const uint32_t* protein_offsets,
    const uint32_t* protein_lengths,
    const AMRGeneEntry* gene_entries,
    AMRHit* hits,
    uint32_t* hit_counts,
    const int num_reads,
    const int num_proteins,
    const AMRDetectionConfig config
) {
    dim3 blockSize(128);
    dim3 gridSize((num_reads + blockSize.x - 1) / blockSize.x);
    size_t sharedMemSize = blockSize.x * 6 * 67 * sizeof(char); // 6 frames * max protein length * threads
    
    translated_alignment_kernel<<<gridSize, blockSize, sharedMemSize>>>(
        reads, read_offsets, read_lengths, read_passes_filter,
        amr_proteins, protein_offsets, protein_lengths, gene_entries,
        hits, hit_counts, num_reads, num_proteins, config
    );
}

extern "C" void launch_extend_alignments_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    const Minimizer* minimizers,
    const uint32_t* minimizer_offsets,
    const uint32_t* minimizer_counts,
    AMRHit* hits,
    uint32_t* hit_counts,
    const char* amr_proteins,
    const uint32_t* protein_offsets,
    const uint32_t* protein_lengths,
    const int num_reads
) {
    dim3 blockSize(256);
    dim3 gridSize((num_reads + blockSize.x - 1) / blockSize.x);
    
    extend_alignments_kernel<<<gridSize, blockSize>>>(
        reads, read_offsets, read_lengths,
        minimizers, minimizer_counts, minimizer_offsets,
        hits, hit_counts, amr_proteins, protein_offsets, protein_lengths,
        num_reads
    );
}

extern "C" void launch_update_coverage_stats_kernel(
    const AMRHit* hits,
    const uint32_t* hit_counts,
    AMRCoverageStats* coverage_stats,
    const int num_reads,
    const int num_genes
) {
    dim3 blockSize(256);
    dim3 gridSize((num_genes + blockSize.x - 1) / blockSize.x);
    
    update_coverage_stats_kernel<<<gridSize, blockSize>>>(
        hits, hit_counts, coverage_stats, num_reads, num_genes
    );
}

extern "C" void launch_finalize_coverage_stats_kernel(
    AMRCoverageStats* coverage_stats,
    const AMRGeneEntry* gene_entries,
    const int num_genes
) {
    dim3 blockSize(256);
    dim3 gridSize((num_genes + blockSize.x - 1) / blockSize.x);
    
    finalize_coverage_stats_kernel<<<gridSize, blockSize>>>(
        coverage_stats, gene_entries, num_genes
    );
}