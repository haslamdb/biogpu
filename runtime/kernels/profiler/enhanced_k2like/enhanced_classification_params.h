// enhanced_classification_params.h
// Enhanced parameters for multi-stage classification with phylogenetic validation

#ifndef ENHANCED_CLASSIFICATION_PARAMS_H
#define ENHANCED_CLASSIFICATION_PARAMS_H

#include "gpu_kraken_classifier.cu"  // Base ClassificationParams
#include "ncbi_taxonomy_loader.h"      // NCBI taxonomy loader
#include "gpu_phylogenetic_calculator.cuh"  // GPU phylogenetic functions
#include <vector>
#include <string>

namespace BioGPU {
namespace Enhanced {

// Phase 1 Enhanced Parameters (no coverage validation yet)
struct Phase1EnhancedParams : public ClassificationParams {
    // Multi-level filtering
    float primary_confidence_threshold = 0.1f;     // Initial filter
    float secondary_confidence_threshold = 0.3f;   // Stricter secondary filter
    
    // Basic k-mer requirements (not coverage-based)
    int min_kmers_for_classification = 10;         // Minimum k-mers needed
    
    // Phylogenetic consistency
    bool enable_phylogenetic_validation = true;
    float max_phylogenetic_distance = 0.5f;        // Maximum allowed inconsistency
    
    // Multi-minimizer validation
    bool require_multiple_minimizers = true;
    int min_distinct_minimizers = 3;               // Different minimizers required
    
    // Weighted phylogenetic parameters
    struct WeightedPhyloParams {
        bool enable_weighted_phylo = true;
        float distance_decay_rate = 2.0f;
        float same_genus_weight = 0.9f;
        float same_family_weight = 0.7f;
        float same_order_weight = 0.5f;
        float same_class_weight = 0.3f;
        float same_phylum_weight = 0.1f;
        enum WeightingMethod {
            EXPONENTIAL_DECAY,
            LINEAR_DECAY,
            RANK_BASED,
            HYBRID
        } weighting_method = HYBRID;
    } weighted_phylo;
    
    // Forward compatibility for post-hoc analysis
    struct ForwardCompatibility {
        bool track_kmer_positions = true;          // Enable for future post-hoc
        bool store_minimizer_details = true;       // Store extra minimizer info
        int max_kmers_to_track = 100;              // Per-read limit for position tracking
        bool enable_genome_mapping_prep = true;    // Prepare for multi-genome mapping
    } forward_compat;
};

// Enhanced result structure (Phase 1)
struct Phase1ClassificationResult {
    // Primary classification
    uint32_t taxon_id;
    float primary_confidence;
    float secondary_confidence;
    bool passed_stage1;
    bool passed_stage2;
    bool is_high_confidence_call;
    
    // K-mer metrics
    uint32_t total_kmers;
    uint32_t classified_kmers;
    uint32_t distinct_minimizers_found;
    
    // Phylogenetic validation
    float phylogenetic_consistency_score;
    bool passed_phylogenetic_filter;
    std::vector<uint32_t> conflicting_taxa;    // Host-side only
    
    // Forward compatibility: K-mer tracking for post-hoc analysis
    struct KmerTrackingData {
        std::vector<uint64_t> minimizer_hashes;   // Minimizers that hit
        std::vector<uint16_t> read_positions;     // Positions in read
        std::vector<uint32_t> contributing_taxa;  // Which taxa each contributed to
        std::vector<float> confidence_weights;    // Weight of each k-mer hit
        uint32_t num_tracked_kmers;               // How many we actually tracked
        bool tracking_overflow;                   // Whether we hit the limit
    } kmer_tracking;
    
    // Quality indicators
    bool has_taxonomic_conflicts;
    bool multiple_minimizer_requirement_met;
    float stage_agreement_score;                  // How much stages agree
    
    // Diagnostic information
    std::string classification_path;              // "Stage1->Stage2", etc.
    uint32_t primary_supporting_votes;
    uint32_t secondary_supporting_votes;
};

// GPU structure for k-mer tracking (forward compatibility)
struct GPUKmerTracker {
    uint64_t* minimizer_hashes;     // [read_id * max_kmers + kmer_idx]
    uint16_t* read_positions;       // Position in read where k-mer starts
    uint32_t* taxa_assignments;     // Which taxon each k-mer contributed to
    float* confidence_weights;      // Weight/contribution of each k-mer
    uint16_t* kmer_counts;          // Number of tracked k-mers per read
    bool* tracking_overflow_flags;  // Whether tracking hit limit per read
    uint32_t max_kmers_per_read;
    uint32_t max_reads;
};

// Phylogenetic-aware classification parameters
struct PhylogeneticClassificationParams {
    bool use_phylogenetic_validation = true;
    std::string taxonomy_nodes_path;     // Path to nodes.dmp
    std::string taxonomy_names_path;     // Path to names.dmp
    
    // Phylogenetic validation thresholds
    float max_phylogenetic_distance = 0.5f;
    float min_consistency_score = 0.7f;
    
    // Weighting method for phylogenetic consistency
    enum PhyloWeightingMethod {
        DISTANCE_BASED,     // Exponential decay based on phylogenetic distance
        RANK_BASED,         // Weight based on LCA taxonomic rank
        HYBRID             // Combination of distance and rank
    } weighting_method = HYBRID;
    
    float distance_decay_rate = 2.0f;   // For exponential decay weighting
    float distance_weight = 0.7f;       // For hybrid weighting
    float rank_weight = 0.3f;           // For hybrid weighting
};

} // namespace Enhanced
} // namespace BioGPU

#endif // ENHANCED_CLASSIFICATION_PARAMS_H