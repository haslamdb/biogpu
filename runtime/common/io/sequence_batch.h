// runtime/common/io/sequence_batch.h
#ifndef SEQUENCE_BATCH_H
#define SEQUENCE_BATCH_H

#include <vector>
#include <string>
#include <algorithm>

namespace BioGPU {

struct SequenceBatch {
    std::vector<std::string> sequences;
    std::vector<std::string> headers;
    std::vector<std::string> qualities;
    std::vector<int> lengths;
    
    // Metadata
    std::string sample_name;
    bool is_paired_end = false;
    int batch_number = 0;
    int read_pair = 1;  // 1 for R1, 2 for R2
    
    // Clear batch data
    void clear() {
        sequences.clear();
        headers.clear();
        qualities.clear();
        lengths.clear();
    }
    
    // Reserve space
    void reserve(size_t n) {
        sequences.reserve(n);
        headers.reserve(n);
        qualities.reserve(n);
        lengths.reserve(n);
    }
    
    // Add a read
    void addRead(const std::string& header, const std::string& seq, 
                 const std::string& qual = "") {
        headers.push_back(header);
        sequences.push_back(seq);
        qualities.push_back(qual);
        lengths.push_back(seq.length());
    }
    
    // Get batch size
    size_t size() const { return sequences.size(); }
    bool empty() const { return sequences.empty(); }
    
    // GPU-friendly format
    void prepareFlatFormat(std::vector<char>& flat_sequences,
                          std::vector<int>& offsets) {
        flat_sequences.clear();
        offsets.clear();
        offsets.reserve(sequences.size() + 1);
        
        offsets.push_back(0);
        for (const auto& seq : sequences) {
            flat_sequences.insert(flat_sequences.end(), seq.begin(), seq.end());
            offsets.push_back(flat_sequences.size());
        }
    }
    
    // Statistics
    size_t getTotalBases() const {
        size_t total = 0;
        for (int len : lengths) total += len;
        return total;
    }
    
    size_t getMaxLength() const {
        if (lengths.empty()) return 0;
        return *std::max_element(lengths.begin(), lengths.end());
    }
    
    double getAverageLength() const {
        if (lengths.empty()) return 0;
        return (double)getTotalBases() / lengths.size();
    }
};

} // namespace BioGPU

#endif // SEQUENCE_BATCH_H