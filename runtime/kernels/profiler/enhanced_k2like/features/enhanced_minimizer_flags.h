// enhanced_minimizer_flags.h
// Enhanced feature flags for minimizers with uniqueness and quality scoring
// Extends the existing MinimizerFlags namespace with new functionality

#ifndef ENHANCED_MINIMIZER_FLAGS_H
#define ENHANCED_MINIMIZER_FLAGS_H

#include <cstdint>
#include <string>

// Extended feature flags bit layout:
// Bits 0-2:   GC content category (0-7) [EXISTING]
// Bits 3-5:   Sequence complexity score (0-7) [EXISTING] 
// Bit 6:      Position bias indicator [EXISTING]
// Bit 7:      Contamination risk flag [EXISTING]
// Bits 8-10:  Uniqueness category (0-7) [NEW]
// Bit 11:     Unique minimizer flag (≥90% uniqueness) [NEW]
// Bit 12:     Rare minimizer flag (≤3 occurrences) [NEW]
// Bit 13:     Reliable minimizer flag (suitable for classification) [NEW]
// Bits 14-16: Co-occurrence score category (0-7) [RESERVED]
// Bit 17:     Position clustering indicator [RESERVED]
// Bits 18-20: Conservation level (species/genus/family/etc) [RESERVED]
// Bits 21-23: Phylogenetic spread category [RESERVED]
// Bits 24-26: ML confidence category [RESERVED]
// Bit 27:     High-confidence classification flag [RESERVED]
// Bits 28-30: Contamination type (human/adapter/low-complexity) [EXISTING]
// Bit 31:     Reserved for future use

namespace EnhancedMinimizerFlags {
    
    // ===========================
    // Uniqueness Scoring (NEW - Bits 8-10)
    // ===========================
    
    constexpr uint32_t UNIQUENESS_MASK = 0x00000700;
    constexpr uint32_t UNIQUENESS_SHIFT = 8;
    
    // Uniqueness categories
    enum class UniquenessCategory : uint8_t {
        EXTREMELY_COMMON = 0,   // 0.00-0.10 (likely false positives)
        VERY_LOW = 1,           // 0.10-0.30 (low specificity)
        LOW = 2,                // 0.30-0.50 (moderate specificity)
        MODERATE = 3,           // 0.50-0.70 (good specificity)
        HIGH = 4,               // 0.70-0.80 (high specificity) 
        VERY_HIGH = 5,          // 0.80-0.90 (very high specificity)
        EXTREMELY_HIGH = 6,     // 0.90-0.95 (excellent specificity)
        SINGLETON_LIKE = 7      // 0.95-1.00 (perfect specificity)
    };
    
    inline uint8_t get_uniqueness_category(uint32_t feature_flags) {
        return (feature_flags & UNIQUENESS_MASK) >> UNIQUENESS_SHIFT;
    }
    
    inline uint32_t set_uniqueness_category(uint32_t feature_flags, uint8_t category) {
        return (feature_flags & ~UNIQUENESS_MASK) | 
               ((static_cast<uint32_t>(category) << UNIQUENESS_SHIFT) & UNIQUENESS_MASK);
    }
    
    inline std::string uniqueness_category_to_string(uint8_t category) {
        const char* names[] = {
            "Extremely Common", "Very Low", "Low", "Moderate",
            "High", "Very High", "Extremely High", "Singleton-like"
        };
        return (category < 8) ? names[category] : "Unknown";
    }
    
    inline float uniqueness_category_to_score(uint8_t category) {
        const float scores[] = {0.05f, 0.20f, 0.40f, 0.60f, 0.75f, 0.85f, 0.925f, 0.975f};
        return (category < 8) ? scores[category] : 0.0f;
    }
    
    // ===========================
    // Quality Flags (NEW - Bits 11-13)
    // ===========================
    
    constexpr uint32_t UNIQUE_MINIMIZER_FLAG = 0x00000800;    // Bit 11
    constexpr uint32_t RARE_MINIMIZER_FLAG = 0x00001000;      // Bit 12  
    constexpr uint32_t RELIABLE_MINIMIZER_FLAG = 0x00002000;  // Bit 13
    
    inline bool is_unique_minimizer(uint32_t feature_flags) {
        return (feature_flags & UNIQUE_MINIMIZER_FLAG) != 0;
    }
    
    inline uint32_t set_unique_minimizer(uint32_t feature_flags, bool is_unique) {
        if (is_unique) {
            return feature_flags | UNIQUE_MINIMIZER_FLAG;
        } else {
            return feature_flags & ~UNIQUE_MINIMIZER_FLAG;
        }
    }
    
    inline bool is_rare_minimizer(uint32_t feature_flags) {
        return (feature_flags & RARE_MINIMIZER_FLAG) != 0;
    }
    
    inline uint32_t set_rare_minimizer(uint32_t feature_flags, bool is_rare) {
        if (is_rare) {
            return feature_flags | RARE_MINIMIZER_FLAG;
        } else {
            return feature_flags & ~RARE_MINIMIZER_FLAG;
        }
    }
    
    inline bool is_reliable_minimizer(uint32_t feature_flags) {
        return (feature_flags & RELIABLE_MINIMIZER_FLAG) != 0;
    }
    
    inline uint32_t set_reliable_minimizer(uint32_t feature_flags, bool is_reliable) {
        if (is_reliable) {
            return feature_flags | RELIABLE_MINIMIZER_FLAG;
        } else {
            return feature_flags & ~RELIABLE_MINIMIZER_FLAG;
        }
    }
    
    // ===========================
    // Reserved Future Features (Bits 14-27)
    // ===========================
    
    // Co-occurrence scoring (bits 14-16)
    constexpr uint32_t COOCCURRENCE_MASK = 0x0001C000;
    constexpr uint32_t COOCCURRENCE_SHIFT = 14;
    
    inline uint8_t get_cooccurrence_category(uint32_t feature_flags) {
        return (feature_flags & COOCCURRENCE_MASK) >> COOCCURRENCE_SHIFT;
    }
    
    inline uint32_t set_cooccurrence_category(uint32_t feature_flags, uint8_t category) {
        return (feature_flags & ~COOCCURRENCE_MASK) | 
               ((static_cast<uint32_t>(category) << COOCCURRENCE_SHIFT) & COOCCURRENCE_MASK);
    }
    
    // Position clustering (bit 17)
    constexpr uint32_t POSITION_CLUSTERING_FLAG = 0x00020000;
    
    inline bool has_position_clustering(uint32_t feature_flags) {
        return (feature_flags & POSITION_CLUSTERING_FLAG) != 0;
    }
    
    inline uint32_t set_position_clustering(uint32_t feature_flags, bool has_clustering) {
        if (has_clustering) {
            return feature_flags | POSITION_CLUSTERING_FLAG;
        } else {
            return feature_flags & ~POSITION_CLUSTERING_FLAG;
        }
    }
    
    // Conservation level (bits 18-20)
    constexpr uint32_t CONSERVATION_MASK = 0x001C0000;
    constexpr uint32_t CONSERVATION_SHIFT = 18;
    
    enum class ConservationLevel : uint8_t {
        SPECIES_SPECIFIC = 0,
        GENUS_SPECIFIC = 1,
        FAMILY_SPECIFIC = 2,
        ORDER_SPECIFIC = 3,
        CLASS_SPECIFIC = 4,
        PHYLUM_SPECIFIC = 5,
        KINGDOM_SPECIFIC = 6,
        UNIVERSAL = 7
    };
    
    inline uint8_t get_conservation_level(uint32_t feature_flags) {
        return (feature_flags & CONSERVATION_MASK) >> CONSERVATION_SHIFT;
    }
    
    inline uint32_t set_conservation_level(uint32_t feature_flags, uint8_t level) {
        return (feature_flags & ~CONSERVATION_MASK) | 
               ((static_cast<uint32_t>(level) << CONSERVATION_SHIFT) & CONSERVATION_MASK);
    }
    
    // Phylogenetic spread (bits 21-23)
    constexpr uint32_t PHYLO_SPREAD_MASK = 0x00E00000;
    constexpr uint32_t PHYLO_SPREAD_SHIFT = 21;
    
    inline uint8_t get_phylogenetic_spread(uint32_t feature_flags) {
        return (feature_flags & PHYLO_SPREAD_MASK) >> PHYLO_SPREAD_SHIFT;
    }
    
    inline uint32_t set_phylogenetic_spread(uint32_t feature_flags, uint8_t spread) {
        return (feature_flags & ~PHYLO_SPREAD_MASK) | 
               ((static_cast<uint32_t>(spread) << PHYLO_SPREAD_SHIFT) & PHYLO_SPREAD_MASK);
    }
    
    // ML confidence (bits 24-26)
    constexpr uint32_t ML_CONFIDENCE_MASK = 0x07000000;
    constexpr uint32_t ML_CONFIDENCE_SHIFT = 24;
    
    inline uint8_t get_ml_confidence(uint32_t feature_flags) {
        return (feature_flags & ML_CONFIDENCE_MASK) >> ML_CONFIDENCE_SHIFT;
    }
    
    inline uint32_t set_ml_confidence(uint32_t feature_flags, uint8_t confidence) {
        return (feature_flags & ~ML_CONFIDENCE_MASK) | 
               ((static_cast<uint32_t>(confidence) << ML_CONFIDENCE_SHIFT) & ML_CONFIDENCE_MASK);
    }
    
    // High-confidence classification (bit 27)
    constexpr uint32_t HIGH_CONFIDENCE_FLAG = 0x08000000;
    
    inline bool is_high_confidence(uint32_t feature_flags) {
        return (feature_flags & HIGH_CONFIDENCE_FLAG) != 0;
    }
    
    inline uint32_t set_high_confidence(uint32_t feature_flags, bool is_high_conf) {
        if (is_high_conf) {
            return feature_flags | HIGH_CONFIDENCE_FLAG;
        } else {
            return feature_flags & ~HIGH_CONFIDENCE_FLAG;
        }
    }
    
    // ===========================
    // Comprehensive Analysis Functions
    // ===========================
    
    struct MinimizerQuality {
        UniquenessCategory uniqueness;
        bool is_unique;
        bool is_rare;
        bool is_reliable;
        bool has_contamination_risk;
        uint8_t gc_category;
        uint8_t complexity_score;
        float estimated_specificity;
        float estimated_sensitivity;
    };
    
    inline MinimizerQuality analyze_minimizer_quality(uint32_t feature_flags, uint16_t ml_weight) {
        MinimizerQuality quality;
        
        quality.uniqueness = static_cast<UniquenessCategory>(get_uniqueness_category(feature_flags));
        quality.is_unique = is_unique_minimizer(feature_flags);
        quality.is_rare = is_rare_minimizer(feature_flags);
        quality.is_reliable = is_reliable_minimizer(feature_flags);
        quality.has_contamination_risk = (feature_flags & 0x80) != 0;  // Bit 7
        quality.gc_category = feature_flags & 0x7;                     // Bits 0-2
        quality.complexity_score = (feature_flags >> 3) & 0x7;         // Bits 3-5
        
        // Estimate specificity based on uniqueness
        quality.estimated_specificity = uniqueness_category_to_score(static_cast<uint8_t>(quality.uniqueness));
        
        // Estimate sensitivity (higher for more common minimizers, but lower specificity)
        quality.estimated_sensitivity = 1.0f - quality.estimated_specificity * 0.5f;
        
        // Adjust for contamination
        if (quality.has_contamination_risk) {
            quality.estimated_specificity *= 0.5f;  // Significant penalty
        }
        
        return quality;
    }
    
    inline std::string format_minimizer_quality(const MinimizerQuality& quality) {
        std::string result = "Quality{";
        result += "uniqueness=" + uniqueness_category_to_string(static_cast<uint8_t>(quality.uniqueness));
        result += ", unique=" + std::string(quality.is_unique ? "true" : "false");
        result += ", rare=" + std::string(quality.is_rare ? "true" : "false");
        result += ", reliable=" + std::string(quality.is_reliable ? "true" : "false");
        result += ", contamination=" + std::string(quality.has_contamination_risk ? "true" : "false");
        result += ", specificity=" + std::to_string(quality.estimated_specificity);
        result += "}";
        return result;
    }
    
    // Quality scoring for classification decisions
    inline float calculate_classification_weight(uint32_t feature_flags, uint16_t ml_weight) {
        MinimizerQuality quality = analyze_minimizer_quality(feature_flags, ml_weight);
        
        float weight = quality.estimated_specificity;
        
        // Boost weight for high-quality minimizers
        if (quality.is_reliable) weight *= 1.5f;
        if (quality.is_unique) weight *= 1.3f;
        
        // Penalize problematic minimizers
        if (quality.has_contamination_risk) weight *= 0.2f;
        if (quality.uniqueness == UniquenessCategory::EXTREMELY_COMMON) weight *= 0.1f;
        
        // Consider complexity and GC content
        if (quality.complexity_score < 2) weight *= 0.8f;  // Low complexity penalty
        if (quality.gc_category == 0 || quality.gc_category == 7) weight *= 0.9f;  // Extreme GC penalty
        
        return std::max(0.0f, std::min(1.0f, weight));
    }
    
    // Filter recommendation based on feature flags
    inline bool should_use_for_classification(uint32_t feature_flags, float quality_threshold = 0.5f) {
        MinimizerQuality quality = analyze_minimizer_quality(feature_flags, 0);
        
        // Immediate disqualifiers
        if (quality.has_contamination_risk) return false;
        if (quality.uniqueness == UniquenessCategory::EXTREMELY_COMMON) return false;
        
        // Quality-based decision
        float overall_quality = quality.estimated_specificity;
        if (quality.is_reliable) overall_quality += 0.2f;
        if (quality.is_unique) overall_quality += 0.1f;
        
        return overall_quality >= quality_threshold;
    }
    
    // ===========================
    // Debug and Visualization Functions
    // ===========================
    
    inline std::string format_feature_flags_detailed(uint32_t feature_flags) {
        std::string result = "FeatureFlags{\n";
        result += "  GC: " + std::to_string(feature_flags & 0x7) + "\n";
        result += "  Complexity: " + std::to_string((feature_flags >> 3) & 0x7) + "\n";
        result += "  Position bias: " + std::string((feature_flags & 0x40) ? "true" : "false") + "\n";
        result += "  Contamination: " + std::string((feature_flags & 0x80) ? "true" : "false") + "\n";
        result += "  Uniqueness: " + uniqueness_category_to_string(get_uniqueness_category(feature_flags)) + "\n";
        result += "  Unique flag: " + std::string(is_unique_minimizer(feature_flags) ? "true" : "false") + "\n";
        result += "  Rare flag: " + std::string(is_rare_minimizer(feature_flags) ? "true" : "false") + "\n";
        result += "  Reliable flag: " + std::string(is_reliable_minimizer(feature_flags) ? "true" : "false") + "\n";
        result += "}";
        return result;
    }
    
    inline void print_feature_flag_bits(uint32_t feature_flags) {
        std::cout << "Feature flags (0x" << std::hex << feature_flags << std::dec << "):" << std::endl;
        for (int i = 31; i >= 0; i--) {
            if (i == 31 || i == 27 || i == 23 || i == 19 || i == 15 || i == 11 || i == 7 || i == 3) {
                std::cout << " ";
            }
            std::cout << ((feature_flags >> i) & 1);
        }
        std::cout << std::endl;
        
        std::cout << "Bit ranges:" << std::endl;
        std::cout << "  31-28: Reserved | Contamination type" << std::endl;
        std::cout << "  27-24: High-conf | ML confidence" << std::endl;
        std::cout << "  23-20: Phylo spread | Conservation" << std::endl;
        std::cout << "  19-16: Position clustering | Co-occurrence" << std::endl;
        std::cout << "  15-12: Reliable | Rare | Unique" << std::endl;
        std::cout << "  11-8:  Uniqueness category" << std::endl;
        std::cout << "  7-4:   Contamination | Position bias" << std::endl;
        std::cout << "  3-0:   Complexity | GC content" << std::endl;
    }
}

#endif // ENHANCED_MINIMIZER_FLAGS_H
