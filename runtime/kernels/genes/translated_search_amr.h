#ifndef TRANSLATED_SEARCH_AMR_H
#define TRANSLATED_SEARCH_AMR_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Create and destroy engine
void* create_translated_search_engine(int batch_size);
void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
void destroy_translated_search_engine(void* engine);

// Load protein database
int load_protein_database(void* engine, const char* db_path);

// Configure engine
void set_smith_waterman_enabled(void* engine, bool enabled);

// Search translated reads
int search_translated_reads(
    void* engine, 
    const char* d_reads, 
    const int* d_read_lengths,
    const int* d_read_offsets, 
    const bool* d_reads_to_process,
    int num_reads, 
    void* results, 
    uint32_t* result_counts
);

#ifdef __cplusplus
}
#endif

#endif // TRANSLATED_SEARCH_AMR_H
