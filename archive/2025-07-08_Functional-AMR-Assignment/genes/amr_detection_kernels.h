// amr_detection_kernels.h
#ifndef AMR_DETECTION_KERNELS_H
#define AMR_DETECTION_KERNELS_H

#include <cuda_runtime.h>
#include "amr_detection_pipeline.h"
#include "ncbi_amr_database_loader.h"

// CUDA kernel wrapper declarations
extern "C" {

void launch_build_bloom_filter_kernel(
    const char* amr_sequences,
    const uint32_t* sequence_offsets,
    const uint32_t* sequence_lengths,
    uint64_t* bloom_filter,
    const uint32_t num_sequences,
    const uint32_t bloom_size,
    const int k
);

void launch_generate_minimizers_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    Minimizer* minimizers,
    uint32_t* minimizer_counts,
    uint32_t* minimizer_offsets,
    const int num_reads,
    const int k,
    const int w
);

void launch_screen_minimizers_kernel(
    const Minimizer* minimizers,
    const uint32_t* minimizer_counts,
    const uint32_t* minimizer_offsets,
    const uint64_t* bloom_filter,
    const int bloom_size,
    bool* read_passes,
    const int num_reads
);

void launch_translated_alignment_kernel(
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
);

void launch_extend_alignments_kernel(
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
);

void launch_update_coverage_stats_kernel(
    const AMRHit* hits,
    const uint32_t* hit_counts,
    AMRCoverageStats* coverage_stats,
    const int num_reads,
    const int num_genes
);

void launch_finalize_coverage_stats_kernel(
    AMRCoverageStats* coverage_stats,
    const AMRGeneEntry* gene_entries,
    const int num_genes
);

// Initialize genetic code table
void initializeGeneticCode();

}

#endif // AMR_DETECTION_KERNELS_H