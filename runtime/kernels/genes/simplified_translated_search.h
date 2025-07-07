// simplified_translated_search.h
#ifndef SIMPLIFIED_TRANSLATED_SEARCH_H
#define SIMPLIFIED_TRANSLATED_SEARCH_H

#include <stdint.h>

// Simple structures for results
struct SimpleAlignment {
    uint32_t read_id;
    uint32_t protein_id;
    uint32_t gene_id;
    int8_t frame;
    uint16_t read_start;
    uint16_t read_end;
    uint16_t protein_start;
    uint16_t protein_end;
    float identity;
    bool valid;
};

struct SimpleCoverage {
    uint32_t gene_id;
    uint32_t total_reads;
    uint32_t covered_positions;
    uint16_t gene_length;
    float percent_coverage;
    float mean_depth;
    uint32_t* position_counts;
};

// C interface
extern "C" {
    // Create and destroy engine
    void* create_simple_translated_search(int batch_size);
    void destroy_simple_translated_search(void* engine);
    
    // Load proteins from FASTA
    int load_simple_proteins(void* engine, const char* fasta_path);
    
    // Search reads
    int search_simple_reads(
        void* engine,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        void* alignments_out,
        void* coverage_out
    );
}

#endif // SIMPLIFIED_TRANSLATED_SEARCH_H
