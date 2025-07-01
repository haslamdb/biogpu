// output/feature_exporter.h
// Feature export pipeline for ML training data generation

#ifndef FEATURE_EXPORTER_H
#define FEATURE_EXPORTER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <cstdint>
#include "../gpu_kraken_types.h"

// Forward declarations
class MinimizerFeatureExtractor;

// Feature export configuration
struct FeatureExportConfig {
    std::string output_directory;
    bool export_tsv = true;
    bool export_hdf5 = false;
    bool export_binary = false;
    
    // Feature selection
    bool include_sequence_features = true;
    bool include_taxonomic_distribution = true;
    bool include_positional_stats = true;
    bool include_cooccurrence = true;
    bool include_contamination = true;
    
    // Parameters
    int k_value = 31;
    int cooccurrence_window = 1000;  // bp window for co-occurrence
    int top_n_cooccurrences = 10;    // Top N co-occurring minimizers to export
    int min_occurrence_threshold = 1; // Minimum occurrences to export
    
    // Performance
    bool use_gpu_acceleration = true;
    size_t batch_size = 1000000;
};

// Decoded feature values for export
struct MinimizerFeatureValues {
    uint8_t gc_category;
    uint8_t complexity_score;
    uint8_t conservation_category;
    float uniqueness_score;
    bool has_position_bias;
    bool is_contamination;
    bool is_unique;
    bool is_rare;
};

// Aggregated features for ML training
struct AggregatedMinimizerFeatures {
    uint64_t minimizer_hash;
    uint32_t total_occurrences;
    uint32_t unique_taxa;
    
    // Positional statistics
    float mean_position;
    float position_std_dev;
    
    // Sequence features (averaged across occurrences)
    float mean_gc_content;
    float mean_complexity;
    
    // ML features
    float avg_uniqueness_score;
    uint8_t conservation_level;
    float contamination_frequency;
    
    // For GPU transfer
    AggregatedMinimizerFeatures() = default;
};

// Co-occurrence information
struct CooccurrencePattern {
    uint64_t neighbor_hash;
    float cooccurrence_score;
    uint32_t cooccurrence_count;
};

// Sparse co-occurrence matrix representation
using CooccurrenceMatrix = std::unordered_map<uint64_t, std::unordered_map<uint64_t, float>>;

// Main feature export class
class FeatureExporter {
private:
    FeatureExportConfig config_;
    
    // Aggregated feature storage
    std::vector<AggregatedMinimizerFeatures> aggregated_features_;
    std::unordered_map<uint64_t, std::unordered_map<uint32_t, uint32_t>> taxonomic_distributions_;
    CooccurrenceMatrix cooccurrence_patterns_;
    
    // Export statistics
    size_t total_minimizers_exported_;
    std::chrono::time_point<std::chrono::high_resolution_clock> export_start_time_;
    
public:
    explicit FeatureExporter(const FeatureExportConfig& config = FeatureExportConfig());
    
    // Main export method
    bool export_training_features(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        const MinimizerFeatureExtractor& feature_extractor,
        const std::unordered_map<uint32_t, std::string>& taxon_names
    );
    
    // Export format methods
    bool export_to_tsv(const std::unordered_map<uint32_t, std::string>& taxon_names);
    bool export_to_hdf5(const std::unordered_map<uint32_t, std::string>& taxon_names);
    bool export_to_binary();
    
    // Feature selection
    void select_features_by_occurrence(uint32_t min_occurrences);
    void select_features_by_uniqueness(float min_uniqueness);
    void filter_contaminated_features();
    
    // Statistics
    void generate_export_summary();
    size_t get_total_features_exported() const { return total_minimizers_exported_; }
    
    // Configuration
    void set_config(const FeatureExportConfig& config) { config_ = config; }
    const FeatureExportConfig& get_config() const { return config_; }
    
private:
    // Internal processing
    bool collect_and_aggregate_features(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        const MinimizerFeatureExtractor& feature_extractor
    );
    
    void calculate_cooccurrence_patterns(
        const std::vector<GPUMinimizerHit>& minimizer_hits
    );
    
    bool aggregate_features_cpu(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        const std::vector<uint64_t>& unique_minimizers,
        const MinimizerFeatureExtractor& feature_extractor
    );
    
    bool aggregate_features_gpu(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        const std::vector<uint64_t>& unique_minimizers
    );
    
    // Format-specific helpers
    void write_tsv_header(std::ofstream& out) const;
    void write_tsv_row(std::ofstream& out, const AggregatedMinimizerFeatures& features) const;
    
    // Disable copy
    FeatureExporter(const FeatureExporter&) = delete;
    FeatureExporter& operator=(const FeatureExporter&) = delete;
};

// Utility functions for feature export
namespace FeatureExportUtils {
    
    // Export features from a completed database
    bool export_features_for_ml_training(
        const std::string& database_path,
        const std::string& output_directory,
        const FeatureExportConfig& config = FeatureExportConfig()
    );
    
    // Feature analysis
    void print_feature_summary(const std::vector<AggregatedMinimizerFeatures>& features);
    
    // Feature encoding helpers
    inline std::string encode_gc_category(uint8_t category) {
        const char* categories[] = {
            "0-20%", "20-30%", "30-40%", "40-50%",
            "50-60%", "60-70%", "70-80%", "80-100%"
        };
        return (category < 8) ? categories[category] : "unknown";
    }
    
    inline std::string encode_complexity_level(uint8_t level) {
        const char* levels[] = {
            "very_low", "low", "medium_low", "medium",
            "medium_high", "high", "very_high", "complex"
        };
        return (level < 8) ? levels[level] : "unknown";
    }
    
    inline std::string encode_conservation_level(uint8_t level) {
        const char* levels[] = {
            "highly_specific", "species_specific", "genus_specific", "family_specific",
            "order_specific", "class_specific", "phylum_specific", "universal"
        };
        return (level < 8) ? levels[level] : "unknown";
    }
    
    // Feature normalization for ML
    inline float normalize_occurrence_count(uint32_t count, uint32_t max_count) {
        return (max_count > 0) ? (float)count / max_count : 0.0f;
    }
    
    inline float normalize_position(float position, float sequence_length) {
        return (sequence_length > 0) ? position / sequence_length : 0.5f;
    }
}

// Command-line integration helper
class FeatureExporterCLI {
public:
    static void add_command_line_options(/* command line parser */);
    static FeatureExportConfig parse_command_line_options(/* parsed args */);
    static bool run_feature_export(int argc, char** argv);
};

#endif // FEATURE_EXPORTER_H