// minimizer_common.h
#ifndef MINIMIZER_COMMON_H
#define MINIMIZER_COMMON_H

#include <cstdint>
#include <vector>
#include <string>

// Structure to hold minimizer information
struct Minimizer {
    uint64_t hash;
    uint32_t position;
    bool is_reverse;
};

// Structure for a batch of reads
struct ReadBatch {
    std::vector<std::string> sequences;
    std::vector<std::string> headers;
    std::vector<std::string> qualities;
    
    void clear() {
        sequences.clear();
        headers.clear();
        qualities.clear();
    }
    
    size_t size() const {
        return sequences.size();
    }
    
    void reserve(size_t n) {
        sequences.reserve(n);
        headers.reserve(n);
        qualities.reserve(n);
    }
};

#endif // MINIMIZER_COMMON_H