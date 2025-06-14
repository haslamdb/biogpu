#ifndef PAIRED_READ_COMMON_H
#define PAIRED_READ_COMMON_H

#include <string>
#include <vector>

// Common structures for paired-end read support across all profilers

struct PairedRead {
    std::string read1;
    std::string read2;      // Empty if single-end
    std::string read_id;
    bool is_paired;         // True if read2 is present
    
    // Constructor for single-end reads with just the read sequence
    explicit PairedRead(const std::string& single_read) 
        : read1(single_read), read_id(""), is_paired(false) {}
    
    // Constructor for paired-end reads
    PairedRead(const std::string& r1, const std::string& r2, const std::string& id = "")
        : read1(r1), read2(r2), read_id(id), is_paired(true) {}
};

struct OrganismProfile {
    uint32_t taxonomy_id;
    std::string name;
    std::string taxonomy_path;
    float abundance;
    float coverage_breadth;
    uint32_t unique_minimizers;
    uint32_t total_hits;
    uint32_t paired_hits;           // NEW: Number of read pairs mapping to this organism
    uint32_t single_hits;           // NEW: Number of single reads mapping
    float paired_concordance;       // NEW: Fraction of pairs where both reads map to same organism
    float confidence_score;
};

struct PairedMinimizerInfo {
    uint64_t minimizer_hash;
    uint32_t read_pair_id;          // Index into the paired reads array
    uint8_t read_number;            // 1 for read1, 2 for read2, 0 for single-end
    uint32_t position_in_read;      // Position of minimizer within the read
};

// Helper function to convert single reads to paired format
inline std::vector<PairedRead> convert_to_paired_format(const std::vector<std::string>& single_reads) {
    std::vector<PairedRead> paired_reads;
    paired_reads.reserve(single_reads.size());
    
    for (size_t i = 0; i < single_reads.size(); i++) {
        PairedRead pr(single_reads[i]);
        pr.read_id = "single_" + std::to_string(i);
        paired_reads.push_back(pr);
    }
    
    return paired_reads;
}

// Helper function to extract read sequences for processing
inline std::vector<std::string> extract_all_sequences(const std::vector<PairedRead>& paired_reads) {
    std::vector<std::string> all_sequences;
    
    for (const auto& pair : paired_reads) {
        all_sequences.push_back(pair.read1);
        if (pair.is_paired && !pair.read2.empty()) {
            all_sequences.push_back(pair.read2);
        }
    }
    
    return all_sequences;
}

// Structure for tracking paired-end concordance
struct PairConcordance {
    uint32_t read_pair_id;
    uint32_t organism_id_r1;        // Best matching organism for read1
    uint32_t organism_id_r2;        // Best matching organism for read2
    float score_r1;
    float score_r2;
    bool is_concordant;             // True if both reads map to same organism
};

#endif // PAIRED_READ_COMMON_H