// processing/minimizer_feature_extractor.h
// Two-pass minimizer feature extraction for enhanced ML-ready metadata

#ifndef MINIMIZER_FEATURE_EXTRACTOR_H
#define MINIMIZER_FEATURE_EXTRACTOR_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cstdint>
#include "../gpu_kraken_types.h"

// Statistics collected for each unique minimizer
struct MinimizerStatistics {
    uint64_t total_occurrences = 0;
    uint32_t unique_taxa = 0;
    std::unordered_map<uint32_t, uint32_t> taxon_occurrences;
    std::unordered_set<uint32_t> genome_taxa;
    uint64_t position_sum = 0;
    float average_position = 0.0f;
    
    // For GPU transfer
    MinimizerStatistics() = default;
};

// Feature extraction statistics
struct FeatureExtractionStats {
    size_t total_unique_minimizers = 0;
    size_t singleton_minimizers = 0;
    size_t rare_minimizers = 0;      // < 10 occurrences
    size_t common_minimizers = 0;    // > 50% of genomes
    size_t contamination_markers = 0;
    size_t total_genomes_processed = 0;
    double avg_taxonomic_spread = 0.0;
    double feature_computation_time = 0.0;
};

// Main feature extraction class
class MinimizerFeatureExtractor {
private:
    // Configuration
    size_t max_minimizers_;
    size_t max_genomes_;
    
    // Statistics storage
    std::unordered_map<uint64_t, MinimizerStatistics> minimizer_stats_map_;
    
    // GPU memory for second pass
    MinimizerStatistics* d_statistics_;
    uint64_t* d_unique_minimizers_;
    size_t num_unique_minimizers_;
    
    // Contamination patterns
    std::unordered_set<uint64_t> human_minimizer_hashes_;
    std::unordered_set<uint64_t> adapter_minimizer_hashes_;
    
    // Processing state
    size_t total_genomes_processed_;
    
public:
    explicit MinimizerFeatureExtractor(size_t max_minimizers = 10000000, 
                                      size_t max_genomes = 100000);
    ~MinimizerFeatureExtractor();
    
    // Two-pass processing
    bool process_first_pass(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits,
        const std::vector<uint32_t>& genome_taxon_ids
    );
    
    bool process_second_pass(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits
    );
    
    // Feature calculation methods
    float calculate_uniqueness_score(uint64_t minimizer_hash) const;
    float calculate_conservation_score(uint64_t minimizer_hash) const;
    bool identify_contamination_markers(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        std::vector<bool>& contamination_flags
    ) const;
    
    // Statistics and export
    FeatureExtractionStats get_statistics() const;
    void export_feature_statistics(const std::string& output_file) const;
    
    // Access to collected data
    const std::unordered_map<uint64_t, MinimizerStatistics>& get_minimizer_stats() const {
        return minimizer_stats_map_;
    }
    
    // Configuration
    void load_contamination_database(const std::string& contamination_db_path);
    void set_contamination_patterns(
        const std::unordered_set<uint64_t>& human_hashes,
        const std::unordered_set<uint64_t>& adapter_hashes
    );
    
private:
    // Internal processing methods
    void calculate_basic_features(GPUMinimizerHit* d_minimizer_hits, size_t num_hits);
    void mark_contamination_features(GPUMinimizerHit* d_minimizer_hits, size_t num_hits);
    void initialize_contamination_patterns();
    
    // Disable copy
    MinimizerFeatureExtractor(const MinimizerFeatureExtractor&) = delete;
    MinimizerFeatureExtractor& operator=(const MinimizerFeatureExtractor&) = delete;
};

// Helper namespace for feature extraction utilities
namespace FeatureExtractionUtils {
    // Run complete two-pass extraction
    bool run_two_pass_feature_extraction(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits,
        const std::vector<uint32_t>& genome_taxon_ids,
        size_t max_minimizers = 10000000
    );
    
    // Feature analysis utilities
    void print_feature_distribution(const std::vector<GPUMinimizerHit>& minimizer_hits);
    
    // Feature encoding/decoding helpers
    inline uint8_t get_gc_category(uint16_t feature_flags) {
        return feature_flags & 0x7;
    }
    
    inline uint8_t get_complexity_score(uint16_t feature_flags) {
        return (feature_flags >> 3) & 0x7;
    }
    
    inline bool has_position_bias(uint16_t feature_flags) {
        return (feature_flags & (1 << 6)) != 0;
    }
    
    inline bool is_contamination_risk(uint16_t feature_flags) {
        return (feature_flags & (1 << 7)) != 0;
    }
    
    inline uint8_t get_conservation_category(uint16_t feature_flags) {
        return (feature_flags >> 8) & 0x7;
    }
    
    inline bool is_unique_minimizer(uint16_t feature_flags) {
        return (feature_flags & (1 << 11)) != 0;
    }
    
    inline bool is_rare_minimizer(uint16_t feature_flags) {
        return (feature_flags & (1 << 12)) != 0;
    }
}

// Enhanced feature flags namespace for conservation and uniqueness
namespace EnhancedMinimizerFlags {
    // Existing flags (bits 0-7) are defined in MinimizerFlags namespace
    // Bits 0-2: GC content category
    // Bits 3-5: Sequence complexity score
    // Bit 6: Position bias indicator
    // Bit 7: Contamination risk flag
    
    // New flags for conservation and uniqueness
    // Bits 8-10: Conservation category (0-7)
    constexpr uint16_t CONSERVATION_MASK = 0x0700;
    constexpr uint16_t CONSERVATION_SHIFT = 8;
    
    // Bit 11: Unique minimizer flag (uniqueness > 0.9)
    constexpr uint16_t UNIQUE_MINIMIZER_FLAG = 0x0800;
    
    // Bit 12: Rare minimizer flag (< 5 occurrences)
    constexpr uint16_t RARE_MINIMIZER_FLAG = 0x1000;
    
    // Bits 13-15: Reserved for future use
    constexpr uint16_t RESERVED_MASK = 0xE000;
    
    // Helper functions for conservation category
    inline uint8_t get_conservation_category(uint16_t feature_flags) {
        return (feature_flags & CONSERVATION_MASK) >> CONSERVATION_SHIFT;
    }
    
    inline uint16_t set_conservation_category(uint16_t feature_flags, uint8_t category) {
        return (feature_flags & ~CONSERVATION_MASK) | ((category << CONSERVATION_SHIFT) & CONSERVATION_MASK);
    }
    
    // Helper functions for unique minimizer flag
    inline bool is_unique_minimizer(uint16_t feature_flags) {
        return (feature_flags & UNIQUE_MINIMIZER_FLAG) != 0;
    }
    
    inline uint16_t set_unique_minimizer(uint16_t feature_flags, bool is_unique) {
        if (is_unique) {
            return feature_flags | UNIQUE_MINIMIZER_FLAG;
        } else {
            return feature_flags & ~UNIQUE_MINIMIZER_FLAG;
        }
    }
    
    // Helper functions for rare minimizer flag
    inline bool is_rare_minimizer(uint16_t feature_flags) {
        return (feature_flags & RARE_MINIMIZER_FLAG) != 0;
    }
    
    inline uint16_t set_rare_minimizer(uint16_t feature_flags, bool is_rare) {
        if (is_rare) {
            return feature_flags | RARE_MINIMIZER_FLAG;
        } else {
            return feature_flags & ~RARE_MINIMIZER_FLAG;
        }
    }
}

#endif // MINIMIZER_FEATURE_EXTRACTOR_H