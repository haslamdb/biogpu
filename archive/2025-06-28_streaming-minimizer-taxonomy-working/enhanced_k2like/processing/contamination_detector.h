// processing/contamination_detector.h
// Contamination detection for minimizers (human, adapter, low-complexity)

#ifndef CONTAMINATION_DETECTOR_H
#define CONTAMINATION_DETECTOR_H

#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <cmath>
#include "../gpu_kraken_types.h"

// Contamination risk levels
enum class ContaminationRiskLevel {
    NO_RISK = 0,
    LOW_RISK = 1,      // Low complexity, homopolymers
    MEDIUM_RISK = 2,   // Adapter sequences
    HIGH_RISK = 3      // Human sequences
};

// Contamination types
enum class ContaminationType {
    HUMAN = 1,
    ADAPTER = 2,
    LOW_COMPLEXITY = 4
};

// Configuration for contamination detection
struct ContaminationConfig {
    int k_value = 31;                           // K-mer size
    size_t max_contamination_patterns = 1000000; // Maximum patterns to store
    bool check_low_complexity = true;
    bool check_homopolymers = true;
    float homopolymer_threshold = 0.7f;         // 70% same base
    int min_complexity_score = 2;               // Minimum complexity bits
    std::string contamination_db_path;          // Path to pre-computed database
};

// Contamination detection statistics
struct ContaminationStats {
    size_t total_minimizers_checked = 0;
    size_t human_contamination_count = 0;
    size_t adapter_contamination_count = 0;
    size_t low_complexity_count = 0;
    size_t homopolymer_count = 0;
    size_t batches_processed = 0;
    double detection_time = 0.0;
};

// Main contamination detector class
class ContaminationDetector {
private:
    ContaminationConfig config_;
    ContaminationStats stats_;
    bool initialized_;
    
    // Host-side contamination patterns
    std::unordered_set<uint64_t> human_minimizer_hashes_;
    std::unordered_set<uint64_t> adapter_minimizer_hashes_;
    std::unordered_set<uint64_t> low_complexity_patterns_;
    
    // GPU-side sorted arrays for binary search
    uint64_t* d_human_hashes_;
    uint64_t* d_adapter_hashes_;
    uint64_t* d_combined_contamination_db_;
    size_t num_human_hashes_;
    size_t num_adapter_hashes_;
    size_t total_contamination_hashes_;
    
public:
    explicit ContaminationDetector(const ContaminationConfig& config = ContaminationConfig());
    ~ContaminationDetector();
    
    // Database loading and management
    bool load_contamination_database(const std::string& database_path);
    bool load_human_minimizers(const std::string& human_minimizer_file);
    bool load_adapter_sequences(const std::string& adapter_file);
    
    // Contamination checking methods
    ContaminationRiskLevel check_minimizer_contamination(uint64_t minimizer_hash) const;
    
    bool mark_contamination_in_batch(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits,
        ContaminationRiskLevel risk_threshold = ContaminationRiskLevel::LOW_RISK
    );
    
    // Batch checking for efficiency
    bool quick_contamination_check(
        const uint64_t* d_minimizer_hashes,
        uint8_t* d_contamination_flags,
        size_t num_minimizers
    );
    
    // Pattern management
    void add_custom_contamination_pattern(uint64_t minimizer_hash, ContaminationType type);
    void clear_contamination_patterns();
    
    // Statistics and reporting
    const ContaminationStats& get_statistics() const { return stats_; }
    void reset_statistics() { stats_ = ContaminationStats(); }
    bool export_contamination_report(const std::string& output_file) const;
    
    // Configuration
    void set_config(const ContaminationConfig& config) { config_ = config; }
    const ContaminationConfig& get_config() const { return config_; }
    
    // Database info
    size_t get_human_pattern_count() const { return human_minimizer_hashes_.size(); }
    size_t get_adapter_pattern_count() const { return adapter_minimizer_hashes_.size(); }
    size_t get_total_pattern_count() const { return total_contamination_hashes_; }
    
private:
    // Internal methods
    void initialize_builtin_patterns();
    bool update_gpu_databases();
    void add_adapter_pattern(const std::string& adapter_sequence);
    uint64_t compute_minimizer_from_sequence(const std::string& sequence, size_t pos) const;
    
    // Disable copy
    ContaminationDetector(const ContaminationDetector&) = delete;
    ContaminationDetector& operator=(const ContaminationDetector&) = delete;
};

// Utility functions for contamination detection
namespace ContaminationDetectionUtils {
    
    // Create contamination database from reference files
    bool create_contamination_database(
        const std::string& human_reference_path,
        const std::string& adapter_sequences_path,
        const std::string& output_database_path,
        int k_value = 31
    );
    
    // Load pre-computed human k-mer database
    bool load_human_kmer_database(
        const std::string& database_path,
        std::unordered_set<uint64_t>& human_kmers
    );
    
    // Common adapter sequences
    std::vector<std::string> get_common_adapter_sequences();
    
    // Analyze contamination patterns in results
    void analyze_contamination_patterns(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        const ContaminationDetector& detector
    );
    
    // Filter out contaminated minimizers
    bool filter_contaminated_minimizers(
        GPUMinimizerHit* d_minimizer_hits,
        size_t& num_hits,
        ContaminationRiskLevel threshold = ContaminationRiskLevel::MEDIUM_RISK
    );
    
    // Check if sequence is low complexity
    inline bool is_low_complexity_sequence(const std::string& sequence, float threshold = 1.5f) {
        if (sequence.empty()) return false;
        
        // Simple Shannon entropy calculation
        int counts[4] = {0, 0, 0, 0};
        for (char c : sequence) {
            switch (c) {
                case 'A': case 'a': counts[0]++; break;
                case 'C': case 'c': counts[1]++; break;
                case 'G': case 'g': counts[2]++; break;
                case 'T': case 't': counts[3]++; break;
            }
        }
        
        float entropy = 0.0f;
        int total = sequence.length();
        for (int i = 0; i < 4; i++) {
            if (counts[i] > 0) {
                float p = (float)counts[i] / total;
                entropy -= p * log2f(p);
            }
        }
        
        return entropy < threshold;
    }
    
    // Check if sequence is a homopolymer run
    inline bool is_homopolymer(const std::string& sequence, float threshold = 0.8f) {
        if (sequence.empty()) return false;
        
        char first_base = sequence[0];
        int same_count = 0;
        
        for (char c : sequence) {
            if (c == first_base) same_count++;
        }
        
        return (float)same_count / sequence.length() >= threshold;
    }
}

// Integration helper for database builder
class ContaminationDetectorIntegration {
public:
    static bool integrate_contamination_check(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits,
        const ContaminationConfig& config = ContaminationConfig()
    ) {
        ContaminationDetector detector(config);
        
        // Load default contamination database if available
        if (!config.contamination_db_path.empty()) {
            detector.load_contamination_database(config.contamination_db_path);
        }
        
        // Mark contamination
        return detector.mark_contamination_in_batch(
            d_minimizer_hits, 
            num_hits, 
            ContaminationRiskLevel::LOW_RISK
        );
    }
};

#endif // CONTAMINATION_DETECTOR_H