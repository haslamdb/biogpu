// amr_detection_kernels_wrapper.cu
#include "amr_detection_kernels.cu"

// GENETIC_CODE is already defined in amr_detection_kernels.cu

extern "C" {
    void initializeGeneticCode() {
        const char genetic_code_host[64] = {
            'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
            'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
            'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
            '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
        };
        cudaMemcpyToSymbol(GENETIC_CODE, genetic_code_host, sizeof(genetic_code_host));
    }
    void launch_build_bloom_filter_kernel(
        const char* amr_sequences,
        const uint32_t* sequence_offsets,
        const uint32_t* sequence_lengths,
        uint64_t* bloom_filter,
        const uint32_t num_sequences,
        const uint32_t bloom_size,
        const int k
    ) {
        int threads = 256;
        int blocks = (num_sequences + threads - 1) / threads;
        build_bloom_filter_kernel<<<blocks, threads>>>(
            amr_sequences, sequence_offsets, sequence_lengths,
            bloom_filter, num_sequences, bloom_size, k
        );
    }

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
    ) {
        int threads = 256;
        int blocks = (num_reads + threads - 1) / threads;
        generate_minimizers_kernel<<<blocks, threads>>>(
            reads, read_offsets, read_lengths,
            minimizers, minimizer_counts, minimizer_offsets,
            num_reads, k, w
        );
    }

    void launch_screen_minimizers_kernel(
        const Minimizer* minimizers,
        const uint32_t* minimizer_counts,
        const uint32_t* minimizer_offsets,
        const uint64_t* bloom_filter,
        const int bloom_size,
        bool* read_passes,
        const int num_reads
    ) {
        int threads = 256;
        int blocks = (num_reads + threads - 1) / threads;
        screen_minimizers_kernel<<<blocks, threads>>>(
            minimizers, minimizer_counts, minimizer_offsets,
            bloom_filter, read_passes, num_reads, bloom_size
        );
    }

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
    ) {
        int threads = 128;  // Less threads due to shared memory usage
        int blocks = (num_reads + threads - 1) / threads;
        
        // Calculate shared memory size for translations
        // Each thread needs space for 6-frame translation
        int max_read_length = config.max_read_length;
        int max_protein_len = (max_read_length / 3) + 1;
        size_t shared_mem_size = threads * max_protein_len * 6 * sizeof(char);
        
        translated_alignment_kernel<<<blocks, threads, shared_mem_size>>>(
            reads, read_offsets, read_lengths, read_passes_filter,
            amr_proteins, protein_offsets, protein_lengths,
            gene_entries, hits, hit_counts,
            num_reads, num_proteins, config
        );
    }

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
    ) {
        int threads = 256;
        int blocks = (num_reads + threads - 1) / threads;
        extend_alignments_kernel<<<blocks, threads>>>(
            reads, read_offsets, read_lengths,
            minimizers, minimizer_counts, minimizer_offsets,
            hits, hit_counts, amr_proteins, protein_offsets, protein_lengths, num_reads
        );
    }

    void launch_update_coverage_stats_kernel(
        const AMRHit* hits,
        const uint32_t* hit_counts,
        AMRCoverageStats* coverage_stats,
        const int num_reads,
        const int num_genes
    ) {
        int threads = 256;
        int blocks = (num_genes + threads - 1) / threads;
        update_coverage_stats_kernel<<<blocks, threads>>>(
            hits, hit_counts, coverage_stats,
            num_reads, num_genes
        );
    }

    void launch_finalize_coverage_stats_kernel(
        AMRCoverageStats* coverage_stats,
        const AMRGeneEntry* gene_entries,
        const int num_genes
    ) {
        int threads = 256;
        int blocks = (num_genes + threads - 1) / threads;
        finalize_coverage_stats_kernel<<<blocks, threads>>>(
            coverage_stats, gene_entries, num_genes
        );
    }
}