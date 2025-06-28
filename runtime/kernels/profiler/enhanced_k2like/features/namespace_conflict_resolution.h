// namespace_conflict_resolution.h
// Resolves conflicts between enhanced_minimizer_flags.h and minimizer_feature_extractor.h

#ifndef NAMESPACE_CONFLICT_RESOLUTION_H
#define NAMESPACE_CONFLICT_RESOLUTION_H

#include <cstdint>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

// Forward declare the conflicting functions to avoid redefinition
namespace MinimizerFlags {
    // These are already defined in minimizer_feature_extractor.h
    // We'll extend them rather than redefine them
    
    // Existing functions (already implemented):
    // uint8_t get_gc_content_category(uint32_t feature_flags);
    // uint8_t get_complexity_score(uint32_t feature_flags);  
    // bool has_position_bias(uint32_t feature_flags);
    // bool has_contamination_risk(uint32_t feature_flags);
    
    // NEW: Extended functions for uniqueness (non-conflicting names)
    inline uint8_t get_uniqueness_category_safe(uint32_t feature_flags) {
        return (feature_flags >> 8) & 0x7;  // Bits 8-10
    }
    
    inline uint32_t set_uniqueness_category_safe(uint32_t feature_flags, uint8_t category) {
        return (feature_flags & ~(0x7 << 8)) | ((category & 0x7) << 8);
    }
    
    inline bool is_unique_minimizer_safe(uint32_t feature_flags) {
        return (feature_flags & (1 << 11)) != 0;  // Bit 11
    }
    
    inline uint32_t set_unique_minimizer_safe(uint32_t feature_flags, bool is_unique) {
        if (is_unique) {
            return feature_flags | (1 << 11);
        } else {
            return feature_flags & ~(1 << 11);
        }
    }
    
    inline bool is_rare_minimizer_safe(uint32_t feature_flags) {
        return (feature_flags & (1 << 12)) != 0;  // Bit 12
    }
    
    inline uint32_t set_rare_minimizer_safe(uint32_t feature_flags, bool is_rare) {
        if (is_rare) {
            return feature_flags | (1 << 12);
        } else {
            return feature_flags & ~(1 << 12);
        }
    }
    
    inline bool is_reliable_minimizer_safe(uint32_t feature_flags) {
        return (feature_flags & (1 << 13)) != 0;  // Bit 13
    }
    
    inline uint32_t set_reliable_minimizer_safe(uint32_t feature_flags, bool is_reliable) {
        if (is_reliable) {
            return feature_flags | (1 << 13);
        } else {
            return feature_flags & ~(1 << 13);
        }
    }
    
    // Uniqueness category names
    inline const char* uniqueness_category_name(uint8_t category) {
        const char* names[] = {
            "Extremely Common", "Very Low", "Low", "Moderate",
            "High", "Very High", "Extremely High", "Singleton-like"
        };
        return (category < 8) ? names[category] : "Unknown";
    }
    
    // Convert ML weight to float uniqueness score
    inline float ml_weight_to_uniqueness(uint16_t ml_weight) {
        return (float)ml_weight / 65535.0f;
    }
    
    // Convert float uniqueness score to ML weight
    inline uint16_t uniqueness_to_ml_weight(float uniqueness_score) {
        uniqueness_score = std::max(0.0f, std::min(1.0f, uniqueness_score));
        return (uint16_t)(uniqueness_score * 65535.0f);
    }
}

// Utility namespace for quality analysis without conflicts
namespace MinimizerQuality {
    
    struct QualityMetrics {
        uint8_t uniqueness_category;
        float uniqueness_score;
        bool is_unique;
        bool is_rare;
        bool is_reliable;
        bool has_contamination;
        uint8_t gc_category;
        uint8_t complexity_score;
        float classification_weight;
    };
    
    inline QualityMetrics analyze_minimizer(uint32_t feature_flags, uint16_t ml_weight) {
        QualityMetrics metrics;
        
        metrics.uniqueness_category = MinimizerFlags::get_uniqueness_category_safe(feature_flags);
        metrics.uniqueness_score = MinimizerFlags::ml_weight_to_uniqueness(ml_weight);
        metrics.is_unique = MinimizerFlags::is_unique_minimizer_safe(feature_flags);
        metrics.is_rare = MinimizerFlags::is_rare_minimizer_safe(feature_flags);
        metrics.is_reliable = MinimizerFlags::is_reliable_minimizer_safe(feature_flags);
        
        // Use existing functions for established features
        metrics.has_contamination = (feature_flags & (1 << 7)) != 0;
        metrics.gc_category = feature_flags & 0x7;
        metrics.complexity_score = (feature_flags >> 3) & 0x7;
        
        // Calculate classification weight
        metrics.classification_weight = metrics.uniqueness_score;
        if (metrics.is_reliable) metrics.classification_weight *= 1.3f;
        if (metrics.has_contamination) metrics.classification_weight *= 0.2f;
        if (metrics.uniqueness_category == 0) metrics.classification_weight *= 0.1f;
        
        return metrics;
    }
    
    inline bool should_use_for_classification(const QualityMetrics& metrics, float threshold = 0.5f) {
        if (metrics.has_contamination) return false;
        if (metrics.uniqueness_category == 0) return false;  // Extremely common
        return metrics.classification_weight >= threshold;
    }
    
    inline void print_quality_summary(const std::vector<QualityMetrics>& metrics) {
        if (metrics.empty()) return;
        
        std::vector<int> uniqueness_dist(8, 0);
        int unique_count = 0, rare_count = 0, reliable_count = 0;
        float avg_score = 0.0f;
        
        for (const auto& m : metrics) {
            uniqueness_dist[m.uniqueness_category]++;
            if (m.is_unique) unique_count++;
            if (m.is_rare) rare_count++;
            if (m.is_reliable) reliable_count++;
            avg_score += m.uniqueness_score;
        }
        
        avg_score /= metrics.size();
        
        std::cout << "\n=== UNIQUENESS QUALITY SUMMARY ===" << std::endl;
        std::cout << "Total minimizers: " << metrics.size() << std::endl;
        std::cout << "Average uniqueness score: " << std::fixed << std::setprecision(3) << avg_score << std::endl;
        std::cout << "Unique minimizers (≥90%): " << unique_count << " (" 
                  << (100.0 * unique_count / metrics.size()) << "%)" << std::endl;
        std::cout << "Rare minimizers (≤3 occurrences): " << rare_count << " (" 
                  << (100.0 * rare_count / metrics.size()) << "%)" << std::endl;
        std::cout << "Reliable minimizers: " << reliable_count << " (" 
                  << (100.0 * reliable_count / metrics.size()) << "%)" << std::endl;
        
        std::cout << "\nUniqueness Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            if (uniqueness_dist[i] > 0) {
                std::cout << "  " << MinimizerFlags::uniqueness_category_name(i) 
                          << ": " << uniqueness_dist[i] << " (" 
                          << (100.0 * uniqueness_dist[i] / metrics.size()) << "%)" << std::endl;
            }
        }
    }
}

#endif // NAMESPACE_CONFLICT_RESOLUTION_H