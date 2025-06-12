// minimizer_extractor.h
#ifndef MINIMIZER_EXTRACTOR_H
#define MINIMIZER_EXTRACTOR_H

#include <vector>
#include <string>
#include "minimizer_common.h"

class MinimizerExtractor {
private:
    // Device memory pointers
    char* d_sequences;
    uint32_t* d_sequence_offsets;
    uint32_t* d_sequence_lengths;
    Minimizer* d_minimizers;
    uint32_t* d_minimizer_counts;
    
    // Pinned host memory for faster transfers
    char* h_sequences;
    uint32_t* h_offsets;
    uint32_t* h_lengths;
    Minimizer* h_minimizers;
    uint32_t* h_counts;
    
    // Parameters
    int k;
    int m;
    size_t allocated_reads;
    size_t allocated_seq_length;
    size_t allocated_pinned_reads;
    size_t allocated_pinned_seq_length;
    
    void allocate_device_memory(size_t num_reads, size_t total_sequence_length);
    
public:
    MinimizerExtractor(int k_size = 31, int window_size = 15);
    ~MinimizerExtractor();
    
    std::vector<std::vector<Minimizer>> extract_minimizers(
        const std::vector<std::string>& sequences
    );
};

#endif // MINIMIZER_EXTRACTOR_H