// phase1_enhanced_classification.cu
// Phase 1 implementation: Core enhancements without coverage validation
// Forward-compatible for future post-hoc positional analysis

#ifndef PHASE1_ENHANCED_CLASSIFICATION_CU
#define PHASE1_ENHANCED_CLASSIFICATION_CU

#include <cuda_runtime.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

// Forward declarations of types from gpu_kraken_classifier.cu
struct PairedRead {
    std::string read1;
    std::string read2;
    std::string read_id;
    bool is_paired;
    
    PairedRead(const std::string& r1, const std::string& r2 = "", 
               const std::string& id = "") 
        : read1(r1), read2(r2), read_id(id), is_paired(!r2.empty()) {}
};

struct GPUCompactHashTable {
    uint32_t* hash_cells;
    uint32_t table_size;
    uint32_t hash_mask;
    uint32_t lca_bits;
    uint32_t hash_bits;
};

struct ClassificationParams {
    int k = 35;
    int ell = 31;
    int spaces = 7;
    float confidence_threshold = 0.0f;
    bool use_spaced_seeds = true;
    int max_ambiguous_bases = 5;
    bool use_paired_end_bonus = true;
    float paired_concordance_weight = 2.0f;
    float min_pair_concordance = 0.5f;
    bool require_both_reads_classified = false;
};

struct PairedReadClassification {
    uint32_t taxon_id;
    float confidence_score;
    uint32_t read1_votes;
    uint32_t read2_votes;
    uint32_t read1_kmers;
    uint32_t read2_kmers;
    uint32_t concordant_votes;
    float pair_concordance;
    bool got_paired_bonus;
    uint32_t read1_best_taxon;
    uint32_t read2_best_taxon;
    float read1_confidence;
    float read2_confidence;
};

struct TaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    char name[64];
};

// Forward declaration of the classifier class
class PairedEndGPUKrakenClassifier;

#include "enhanced_classification_params.h"
#include "fast_enhanced_classifier.h"
#include "phase1_enhanced_classifier_with_phylo.h"
#include "../tools/compact_gpu_taxonomy.h"
#include "gpu_minimizer_extraction.cuh"

using namespace BioGPU::Enhanced;
using namespace BioGPU::CompactTaxonomy;

// Forward declarations of device functions
__device__ uint32_t find_lca_gpu_simple(
    uint32_t taxon1, 
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup,
    uint32_t max_taxon_id);

__device__ float calculate_weighted_phylo_consistency_gpu(
    const uint32_t* taxa,
    const uint32_t* votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup,
    Phase1EnhancedParams::WeightedPhyloParams params);

__device__ float calculate_simple_phylo_consistency_gpu(
    const uint32_t* taxa,
    const uint32_t* votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup);

__device__ float calculate_normalized_phylo_distance_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup);

__device__ bool are_taxa_related_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    int max_steps);

// Phase 1 Enhanced Classification Kernel
__global__ void phase1_enhanced_classification_kernel(
    const char* sequences,
    const uint32_t* read_offsets,
    const uint32_t* read_lengths,
    const GPUCompactHashTable* hash_table,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup,
    GPUKmerTracker* kmer_tracker,
    Phase1ClassificationResult* results,
    int num_reads,
    Phase1EnhancedParams params) {
    
    int read_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (read_id >= num_reads) return;
    
    const char* read_seq = sequences + read_offsets[read_id];
    uint32_t read_length = read_lengths[read_id];
    
    Phase1ClassificationResult& result = results[read_id];
    
    // Initialize result
    result.taxon_id = 0;
    result.primary_confidence = 0.0f;
    result.secondary_confidence = 0.0f;
    result.passed_stage1 = false;
    result.passed_stage2 = false;
    result.is_high_confidence_call = false;
    result.total_kmers = 0;
    result.classified_kmers = 0;
    result.distinct_minimizers_found = 0;
    result.phylogenetic_consistency_score = 0.0f;
    result.passed_phylogenetic_filter = false;
    result.has_taxonomic_conflicts = false;
    result.multiple_minimizer_requirement_met = false;
    
    if (read_length < params.k) return;
    
    // Local arrays for k-mer processing
    const int MAX_LOCAL_TAXA = 32;
    const int MAX_LOCAL_MINIMIZERS = 64;
    
    uint32_t hit_taxa[MAX_LOCAL_TAXA];
    uint32_t hit_votes[MAX_LOCAL_TAXA];
    uint64_t unique_minimizers[MAX_LOCAL_MINIMIZERS];
    uint16_t minimizer_positions[MAX_LOCAL_MINIMIZERS];
    uint32_t minimizer_taxa[MAX_LOCAL_MINIMIZERS];
    
    int num_taxa = 0;
    int num_unique_minimizers = 0;
    uint64_t last_minimizer = UINT64_MAX;
    
    // Forward compatibility: K-mer tracking setup
    uint64_t* tracked_minimizers = nullptr;
    uint16_t* tracked_positions = nullptr;
    uint32_t* tracked_taxa = nullptr;
    float* tracked_weights = nullptr;
    uint32_t max_tracking = params.forward_compat.max_kmers_to_track;
    uint32_t tracked_count = 0;
    bool tracking_enabled = params.forward_compat.track_kmer_positions;
    
    if (tracking_enabled && kmer_tracker) {
        tracked_minimizers = kmer_tracker->minimizer_hashes + (read_id * max_tracking);
        tracked_positions = kmer_tracker->read_positions + (read_id * max_tracking);
        tracked_taxa = kmer_tracker->taxa_assignments + (read_id * max_tracking);
        tracked_weights = kmer_tracker->confidence_weights + (read_id * max_tracking);
    }
    
    // ==================================================
    // STAGE 1: PRIMARY K-MER PROCESSING
    // ==================================================
    
    int total_kmers = read_length - params.k + 1;
    result.total_kmers = total_kmers;
    
    for (int pos = 0; pos < total_kmers; pos++) {
        // Extract minimizer
        uint64_t minimizer = extract_minimizer_sliding_window(
            read_seq, pos, params.k, params.ell, params.spaces, 0
        );
        
        if (minimizer == UINT64_MAX) continue;
        
        // Track distinct minimizers
        if (minimizer != last_minimizer) {
            if (num_unique_minimizers < MAX_LOCAL_MINIMIZERS) {
                unique_minimizers[num_unique_minimizers] = minimizer;
                minimizer_positions[num_unique_minimizers] = pos;
                num_unique_minimizers++;
            }
            last_minimizer = minimizer;
        }
        
        // Look up in hash table
        uint32_t lca = lookup_lca_gpu(hash_table, minimizer);
        if (lca == 0) continue;
        
        result.classified_kmers++;
        
        // Track for post-hoc analysis if enabled
        if (tracking_enabled && tracked_count < max_tracking) {
            tracked_minimizers[tracked_count] = minimizer;
            tracked_positions[tracked_count] = pos;
            tracked_taxa[tracked_count] = lca;
            tracked_weights[tracked_count] = 1.0f;  // Will be refined in Stage 2
            tracked_count++;
        }
        
        // Update taxon votes
        bool found = false;
        for (int i = 0; i < num_taxa; i++) {
            if (hit_taxa[i] == lca) {
                hit_votes[i]++;
                found = true;
                break;
            }
        }
        
        if (!found && num_taxa < MAX_LOCAL_TAXA) {
            hit_taxa[num_taxa] = lca;
            hit_votes[num_taxa] = 1;
            if (num_unique_minimizers <= MAX_LOCAL_MINIMIZERS) {
                minimizer_taxa[num_unique_minimizers - 1] = lca;
            }
            num_taxa++;
        }
    }
    
    result.distinct_minimizers_found = num_unique_minimizers;
    
    // Store tracking results
    if (tracking_enabled && kmer_tracker) {
        kmer_tracker->kmer_counts[read_id] = tracked_count;
        kmer_tracker->tracking_overflow_flags[read_id] = (tracked_count >= max_tracking);
    }
    
    if (result.classified_kmers == 0) return;
    
    // Find best taxon for Stage 1
    uint32_t best_taxon_stage1 = 0;
    uint32_t best_votes_stage1 = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (hit_votes[i] > best_votes_stage1) {
            best_votes_stage1 = hit_votes[i];
            best_taxon_stage1 = hit_taxa[i];
        }
    }
    
    result.primary_confidence = total_kmers > 0 ? (float)best_votes_stage1 / total_kmers : 0.0f;
    result.primary_supporting_votes = best_votes_stage1;
    
    // ==================================================
    // STAGE 1 VALIDATION
    // ==================================================
    
    // Check minimum k-mer requirement
    bool meets_kmer_threshold = best_votes_stage1 >= params.min_kmers_for_classification;
    
    // Check minimizer diversity requirement
    bool meets_minimizer_requirement = true;
    if (params.require_multiple_minimizers) {
        meets_minimizer_requirement = num_unique_minimizers >= params.min_distinct_minimizers;
    }
    result.multiple_minimizer_requirement_met = meets_minimizer_requirement;
    
    // Stage 1 decision
    result.passed_stage1 = (result.primary_confidence >= params.primary_confidence_threshold) &&
                           meets_kmer_threshold &&
                           meets_minimizer_requirement;
    
    if (!result.passed_stage1) {
        result.taxon_id = 0;  // Unclassified
        return;
    }
    
    // ==================================================
    // STAGE 2: ENHANCED VALIDATION
    // ==================================================
    
    // Phylogenetic consistency validation
    float phylo_score = 1.0f;
    bool passes_phylo = true;
    
    if (params.enable_phylogenetic_validation && num_taxa > 1) {
        if (params.weighted_phylo.enable_weighted_phylo) {
            // Weighted phylogenetic consistency
            phylo_score = calculate_weighted_phylo_consistency_gpu(
                hit_taxa, hit_votes, num_taxa, best_taxon_stage1,
                parent_lookup, depth_lookup, params.weighted_phylo
            );
        } else {
            // Simple phylogenetic consistency
            phylo_score = calculate_simple_phylo_consistency_gpu(
                hit_taxa, hit_votes, num_taxa, best_taxon_stage1,
                parent_lookup, depth_lookup
            );
        }
        
        passes_phylo = phylo_score >= (1.0f - params.max_phylogenetic_distance);
    }
    
    result.phylogenetic_consistency_score = phylo_score;
    result.passed_phylogenetic_filter = passes_phylo;
    
    // Detect taxonomic conflicts
    result.has_taxonomic_conflicts = (num_taxa > 3) && (phylo_score < 0.7f);
    
    // ==================================================
    // STAGE 2 COMPOSITE CONFIDENCE CALCULATION
    // ==================================================
    
    // Calculate enhanced confidence incorporating all factors
    float confidence_components[4];
    float component_weights[4];
    
    // Component 1: Primary confidence (k-mer vote ratio)
    confidence_components[0] = result.primary_confidence;
    component_weights[0] = 0.4f;
    
    // Component 2: Phylogenetic consistency
    confidence_components[1] = phylo_score;
    component_weights[1] = params.enable_phylogenetic_validation ? 0.3f : 0.0f;
    
    // Component 3: Minimizer diversity bonus
    float diversity_bonus = num_unique_minimizers >= params.min_distinct_minimizers ? 1.0f : 0.5f;
    confidence_components[2] = diversity_bonus;
    component_weights[2] = 0.2f;
    
    // Component 4: Taxonomic specificity (fewer conflicting taxa = higher confidence)
    float specificity_score = num_taxa <= 2 ? 1.0f : fmaxf(0.3f, 1.0f - (num_taxa - 2) * 0.1f);
    confidence_components[3] = specificity_score;
    component_weights[3] = 0.1f;
    
    // Normalize weights if phylogenetic validation is disabled
    if (!params.enable_phylogenetic_validation) {
        component_weights[0] = 0.5f;
        component_weights[2] = 0.3f;
        component_weights[3] = 0.2f;
    }
    
    // Calculate weighted composite confidence
    result.secondary_confidence = 0.0f;
    for (int i = 0; i < 4; i++) {
        result.secondary_confidence += confidence_components[i] * component_weights[i];
    }
    
    // Update k-mer weights for post-hoc analysis
    if (tracking_enabled && tracked_count > 0) {
        float kmer_weight = result.secondary_confidence / tracked_count;
        for (uint32_t i = 0; i < tracked_count; i++) {
            tracked_weights[i] = kmer_weight;
        }
    }
    
    // ==================================================
    // FINAL DECISION
    // ==================================================
    
    result.passed_stage2 = (result.secondary_confidence >= params.secondary_confidence_threshold) &&
                           passes_phylo;
    
    if (result.passed_stage2) {
        result.taxon_id = best_taxon_stage1;
        result.is_high_confidence_call = true;
        // result.classification_path = "Stage1->Stage2";
    } else if (result.passed_stage1) {
        result.taxon_id = best_taxon_stage1;
        result.is_high_confidence_call = false;
        // result.classification_path = "Stage1_only";
    } else {
        result.taxon_id = 0;  // Unclassified
        // result.classification_path = "Failed";
    }
    
    result.secondary_supporting_votes = best_votes_stage1;
    result.stage_agreement_score = result.passed_stage2 ? 1.0f : 
                                  (result.secondary_confidence / params.secondary_confidence_threshold);
}

// Device function for weighted phylogenetic consistency
__device__ float calculate_weighted_phylo_consistency_gpu(
    const uint32_t* taxa,
    const uint32_t* votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup,
    Phase1EnhancedParams::WeightedPhyloParams params) {
    
    if (num_taxa <= 1) return 1.0f;
    
    float total_weighted_score = 0.0f;
    uint32_t total_votes = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (votes[i] > 0) {
            // Calculate phylogenetic distance
            float distance = calculate_normalized_phylo_distance_gpu(
                primary_taxon, taxa[i], parent_lookup, depth_lookup
            );
            
            // Calculate weight based on distance
            float weight = 1.0f;
            switch (params.weighting_method) {
                case Phase1EnhancedParams::WeightedPhyloParams::EXPONENTIAL_DECAY:
                    weight = expf(-distance * params.distance_decay_rate);
                    break;
                case Phase1EnhancedParams::WeightedPhyloParams::LINEAR_DECAY:
                    weight = fmaxf(0.0f, 1.0f - distance);
                    break;
                case Phase1EnhancedParams::WeightedPhyloParams::HYBRID:
                default:
                    weight = 0.7f * expf(-distance * params.distance_decay_rate) + 
                            0.3f * fmaxf(0.0f, 1.0f - distance);
                    break;
            }
            
            total_weighted_score += weight * votes[i];
            total_votes += votes[i];
        }
    }
    
    return total_votes > 0 ? total_weighted_score / total_votes : 0.0f;
}

// Device function for simple phylogenetic consistency (fallback)
__device__ float calculate_simple_phylo_consistency_gpu(
    const uint32_t* taxa,
    const uint32_t* votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup) {
    
    if (num_taxa <= 1) return 1.0f;
    
    uint32_t consistent_votes = 0;
    uint32_t total_votes = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (votes[i] > 0) {
            total_votes += votes[i];
            
            // Check if taxa are phylogenetically related
            if (are_taxa_related_gpu(primary_taxon, taxa[i], parent_lookup, 10)) {
                consistent_votes += votes[i];
            }
        }
    }
    
    return total_votes > 0 ? (float)consistent_votes / total_votes : 0.0f;
}

// Device function to calculate normalized phylogenetic distance
__device__ float calculate_normalized_phylo_distance_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    const uint32_t* depth_lookup) {
    
    if (taxon1 == taxon2) return 0.0f;
    
    // Find LCA and calculate distance
    uint32_t lca = find_lca_gpu_simple(taxon1, taxon2, parent_lookup, depth_lookup, 1000000);
    
    uint32_t depth1 = depth_lookup[taxon1];
    uint32_t depth2 = depth_lookup[taxon2];
    uint32_t lca_depth = depth_lookup[lca];
    
    int total_distance = (depth1 - lca_depth) + (depth2 - lca_depth);
    
    // Normalize by maximum reasonable distance
    return fminf(1.0f, total_distance / 20.0f);
}

// Device function to check if taxa are related
__device__ bool are_taxa_related_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    int max_steps) {
    
    if (taxon1 == taxon2) return true;
    
    // Walk up from taxon1 to see if we find taxon2
    uint32_t current = taxon1;
    for (int i = 0; i < max_steps && current != 0 && current != 1; i++) {
        if (current == taxon2) return true;
        current = parent_lookup[current];
    }
    
    // Walk up from taxon2 to see if we find taxon1
    current = taxon2;
    for (int i = 0; i < max_steps && current != 0 && current != 1; i++) {
        if (current == taxon1) return true;
        current = parent_lookup[current];
    }
    
    return false;
}

// Host class for Phase 1 Enhanced Classification
// This works alongside the existing PairedEndGPUKrakenClassifier to add enhanced features
class Phase1EnhancedClassifier {
private:
    Phase1EnhancedParams enhanced_params;
    GPUKmerTracker* d_kmer_tracker;
    Phase1ClassificationResult* d_enhanced_results;
    int max_reads;
    bool tracking_enabled;
    
    // Core GPU data structures (will be shared/copied from base classifier)
    GPUCompactHashTable* d_hash_table;
    TaxonomyNode* d_taxonomy_tree;
    uint32_t* d_parent_lookup;
    uint32_t* d_depth_lookup;
    uint32_t num_taxonomy_nodes;
    
    // Read data on GPU
    char* d_sequences;
    uint32_t* d_read_offsets;
    uint32_t* d_read_lengths;
    size_t d_sequences_capacity;
    size_t max_sequence_length;
    
    // Base classifier instance for standard functionality
    PairedEndGPUKrakenClassifier* base_classifier;
    bool owns_base_classifier;
    
public:
    // Constructor that creates its own base classifier
    Phase1EnhancedClassifier(const Phase1EnhancedParams& config, int max_read_count = 50000) 
        : enhanced_params(config), max_reads(max_read_count), 
          d_kmer_tracker(nullptr), d_enhanced_results(nullptr), d_depth_lookup(nullptr),
          d_hash_table(nullptr), d_taxonomy_tree(nullptr), d_parent_lookup(nullptr),
          num_taxonomy_nodes(0), d_sequences(nullptr), d_read_offsets(nullptr),
          d_read_lengths(nullptr), d_sequences_capacity(0), max_sequence_length(1000),
          owns_base_classifier(true) {
        
        // Create base classifier with same basic params
        base_classifier = new PairedEndGPUKrakenClassifier(config);
        
        initialize();
    }
    
    // Constructor that uses an existing base classifier
    Phase1EnhancedClassifier(PairedEndGPUKrakenClassifier* existing_classifier,
                           const Phase1EnhancedParams& config, int max_read_count = 50000) 
        : enhanced_params(config), max_reads(max_read_count), 
          d_kmer_tracker(nullptr), d_enhanced_results(nullptr), d_depth_lookup(nullptr),
          d_hash_table(nullptr), d_taxonomy_tree(nullptr), d_parent_lookup(nullptr),
          num_taxonomy_nodes(0), d_sequences(nullptr), d_read_offsets(nullptr),
          d_read_lengths(nullptr), d_sequences_capacity(0), max_sequence_length(1000),
          base_classifier(existing_classifier), owns_base_classifier(false) {
        
        initialize();
    }
    
private:
    void initialize() {
        
        tracking_enabled = enhanced_params.forward_compat.track_kmer_positions;
        
        std::cout << "Phase 1 Enhanced Classifier initialized" << std::endl;
        std::cout << "  Two-stage filtering: " << enhanced_params.primary_confidence_threshold 
                  << " -> " << enhanced_params.secondary_confidence_threshold << std::endl;
        std::cout << "  Phylogenetic validation: " << (enhanced_params.enable_phylogenetic_validation ? "ON" : "OFF") << std::endl;
        std::cout << "  Weighted phylo: " << (enhanced_params.weighted_phylo.enable_weighted_phylo ? "ON" : "OFF") << std::endl;
        std::cout << "  K-mer tracking: " << (tracking_enabled ? "ON" : "OFF") << std::endl;
        
        allocate_gpu_memory();
    }
    
public:
    ~Phase1EnhancedClassifier() {
        free_enhanced_gpu_memory();
        if (owns_base_classifier && base_classifier) {
            delete base_classifier;
        }
    }
    
    // Load database through base classifier and setup enhanced structures
    bool load_database(const std::string& database_directory) {
        if (!base_classifier) return false;
        
        // Load database through base classifier
        if (!base_classifier->load_database(database_directory)) {
            std::cerr << "Failed to load database through base classifier" << std::endl;
            return false;
        }
        
        // Setup enhanced structures (we'll need to load additional data for depth lookup)
        std::string taxonomy_file = database_directory + "/taxonomy.tsv";
        if (!load_depth_lookup(taxonomy_file)) {
            std::cerr << "Warning: Could not load depth lookup table" << std::endl;
        }
        
        return true;
    }
    
    // Enhanced classification with all the new features
    std::vector<Phase1ClassificationResult> classify_enhanced(
        const std::vector<std::string>& reads) {
        
        if (reads.empty()) return {};
        
        int num_reads = std::min((int)reads.size(), max_reads);
        
        // First, get baseline classification from base classifier
        std::vector<PairedRead> paired_reads;
        paired_reads.reserve(num_reads);
        for (int i = 0; i < num_reads; i++) {
            paired_reads.emplace_back(reads[i], "", std::to_string(i));
        }
        
        // Get initial classifications
        auto base_results = base_classifier->classify_paired_reads(paired_reads);
        
        // Transfer reads to GPU for enhanced processing
        prepare_reads_batch(reads);
        
        // For now, we'll need to pass the database pointers through a setup method
        // This is a temporary solution until we can properly access base classifier's GPU data
        
        // Launch enhanced classification kernel
        int blocks = (num_reads + 255) / 256;
        
        if (d_hash_table && d_taxonomy_tree && d_parent_lookup) {
            phase1_enhanced_classification_kernel<<<blocks, 256>>>(
                d_sequences,
                d_read_offsets,
                d_read_lengths,
                d_hash_table,
                d_taxonomy_tree,
                d_parent_lookup,
                d_depth_lookup,
                d_kmer_tracker,
                d_enhanced_results,
                num_reads,
                enhanced_params
            );
        } else {
            std::cerr << "Warning: Enhanced classification requires database pointers to be set" << std::endl;
            // Fall back to converting base results
            return convert_base_to_enhanced_results(base_results);
        }
        
        cudaDeviceSynchronize();
        
        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in enhanced classification: " << cudaGetErrorString(error) << std::endl;
            return {};
        }
        
        // Retrieve results
        std::vector<Phase1ClassificationResult> results(num_reads);
        cudaMemcpy(results.data(), d_enhanced_results,
                   num_reads * sizeof(Phase1ClassificationResult),
                   cudaMemcpyDeviceToHost);
        
        // Retrieve k-mer tracking data if enabled
        if (tracking_enabled) {
            retrieve_kmer_tracking_data(results, num_reads);
        }
        
        return results;
    }
    
    void print_enhanced_statistics(const std::vector<Phase1ClassificationResult>& results) {
        if (results.empty()) return;
        
        int total = results.size();
        int stage1_passed = 0;
        int stage2_passed = 0;
        int high_confidence = 0;
        int phylo_passed = 0;
        int minimizer_req_met = 0;
        int has_conflicts = 0;
        
        float total_primary_conf = 0.0f;
        float total_secondary_conf = 0.0f;
        float total_phylo_score = 0.0f;
        
        for (const auto& result : results) {
            if (result.passed_stage1) stage1_passed++;
            if (result.passed_stage2) stage2_passed++;
            if (result.is_high_confidence_call) high_confidence++;
            if (result.passed_phylogenetic_filter) phylo_passed++;
            if (result.multiple_minimizer_requirement_met) minimizer_req_met++;
            if (result.has_taxonomic_conflicts) has_conflicts++;
            
            if (result.taxon_id > 0) {
                total_primary_conf += result.primary_confidence;
                total_secondary_conf += result.secondary_confidence;
                total_phylo_score += result.phylogenetic_consistency_score;
            }
        }
        
        int classified = stage2_passed > 0 ? stage2_passed : stage1_passed;
        
        std::cout << "\n=== PHASE 1 ENHANCED CLASSIFICATION STATISTICS ===" << std::endl;
        std::cout << "Total reads: " << total << std::endl;
        std::cout << "Stage 1 passed: " << stage1_passed << " (" 
                  << std::fixed << std::setprecision(1) << (100.0 * stage1_passed / total) << "%)" << std::endl;
        std::cout << "Stage 2 passed: " << stage2_passed << " (" 
                  << (100.0 * stage2_passed / total) << "%)" << std::endl;
        std::cout << "High confidence: " << high_confidence << " (" 
                  << (100.0 * high_confidence / total) << "%)" << std::endl;
        
        if (enhanced_params.enable_phylogenetic_validation) {
            std::cout << "Phylogenetic filter passed: " << phylo_passed << " (" 
                      << (100.0 * phylo_passed / total) << "%)" << std::endl;
        }
        
        std::cout << "Minimizer requirement met: " << minimizer_req_met << " (" 
                  << (100.0 * minimizer_req_met / total) << "%)" << std::endl;
        std::cout << "Taxonomic conflicts detected: " << has_conflicts << " (" 
                  << (100.0 * has_conflicts / total) << "%)" << std::endl;
        
        if (classified > 0) {
            std::cout << "\nAverage confidence scores:" << std::endl;
            std::cout << "  Primary: " << std::fixed << std::setprecision(3) 
                      << (total_primary_conf / classified) << std::endl;
            std::cout << "  Secondary: " << (total_secondary_conf / classified) << std::endl;
            if (enhanced_params.enable_phylogenetic_validation) {
                std::cout << "  Phylogenetic: " << (total_phylo_score / classified) << std::endl;
            }
        }
        
        // K-mer tracking statistics
        if (tracking_enabled) {
            int tracking_overflow_count = 0;
            uint32_t total_tracked_kmers = 0;
            
            for (const auto& result : results) {
                if (result.kmer_tracking.tracking_overflow) tracking_overflow_count++;
                total_tracked_kmers += result.kmer_tracking.num_tracked_kmers;
            }
            
            std::cout << "\nK-mer tracking statistics:" << std::endl;
            std::cout << "  Total k-mers tracked: " << total_tracked_kmers << std::endl;
            std::cout << "  Average per read: " << std::fixed << std::setprecision(1) 
                      << ((float)total_tracked_kmers / total) << std::endl;
            std::cout << "  Tracking overflow: " << tracking_overflow_count << " reads" << std::endl;
            std::cout << "  Ready for post-hoc analysis: " << (tracking_overflow_count == 0 ? "YES" : "PARTIAL") << std::endl;
        }
    }
    
    // Set database pointers (temporary solution for accessing base classifier's GPU data)
    void set_database_pointers(GPUCompactHashTable* hash_table, TaxonomyNode* taxonomy_tree,
                              uint32_t* parent_lookup, uint32_t num_nodes) {
        d_hash_table = hash_table;
        d_taxonomy_tree = taxonomy_tree;
        d_parent_lookup = parent_lookup;
        num_taxonomy_nodes = num_nodes;
    }
    
private:
    bool allocate_gpu_memory() {
        // Allocate results array
        cudaMalloc(&d_enhanced_results, max_reads * sizeof(Phase1ClassificationResult));
        
        // Allocate read data arrays
        d_sequences_capacity = max_reads * max_sequence_length;
        cudaMalloc(&d_sequences, d_sequences_capacity);
        cudaMalloc(&d_read_offsets, max_reads * sizeof(uint32_t));
        cudaMalloc(&d_read_lengths, max_reads * sizeof(uint32_t));
        
        // Allocate k-mer tracking if enabled
        if (tracking_enabled) {
            d_kmer_tracker = new GPUKmerTracker();
            
            uint32_t max_kmers_total = max_reads * enhanced_params.forward_compat.max_kmers_to_track;
            
            cudaMalloc(&d_kmer_tracker->minimizer_hashes, max_kmers_total * sizeof(uint64_t));
            cudaMalloc(&d_kmer_tracker->read_positions, max_kmers_total * sizeof(uint16_t));
            cudaMalloc(&d_kmer_tracker->taxa_assignments, max_kmers_total * sizeof(uint32_t));
            cudaMalloc(&d_kmer_tracker->confidence_weights, max_kmers_total * sizeof(float));
            cudaMalloc(&d_kmer_tracker->kmer_counts, max_reads * sizeof(uint16_t));
            cudaMalloc(&d_kmer_tracker->tracking_overflow_flags, max_reads * sizeof(bool));
            
            d_kmer_tracker->max_kmers_per_read = enhanced_params.forward_compat.max_kmers_to_track;
            d_kmer_tracker->max_reads = max_reads;
            
            std::cout << "K-mer tracking allocated: " << (max_kmers_total * 4 * sizeof(uint64_t) / 1024 / 1024) 
                      << " MB" << std::endl;
        }
        
        return true;
    }
    
    void free_enhanced_gpu_memory() {
        if (d_enhanced_results) { cudaFree(d_enhanced_results); d_enhanced_results = nullptr; }
        if (d_sequences) { cudaFree(d_sequences); d_sequences = nullptr; }
        if (d_read_offsets) { cudaFree(d_read_offsets); d_read_offsets = nullptr; }
        if (d_read_lengths) { cudaFree(d_read_lengths); d_read_lengths = nullptr; }
        if (d_depth_lookup) { cudaFree(d_depth_lookup); d_depth_lookup = nullptr; }
        
        if (d_kmer_tracker) {
            if (d_kmer_tracker->minimizer_hashes) cudaFree(d_kmer_tracker->minimizer_hashes);
            if (d_kmer_tracker->read_positions) cudaFree(d_kmer_tracker->read_positions);
            if (d_kmer_tracker->taxa_assignments) cudaFree(d_kmer_tracker->taxa_assignments);
            if (d_kmer_tracker->confidence_weights) cudaFree(d_kmer_tracker->confidence_weights);
            if (d_kmer_tracker->kmer_counts) cudaFree(d_kmer_tracker->kmer_counts);
            if (d_kmer_tracker->tracking_overflow_flags) cudaFree(d_kmer_tracker->tracking_overflow_flags);
            
            delete d_kmer_tracker;
            d_kmer_tracker = nullptr;
        }
    }
    
    void retrieve_kmer_tracking_data(std::vector<Phase1ClassificationResult>& results, int num_reads) {
        // Retrieve k-mer counts and overflow flags
        std::vector<uint16_t> kmer_counts(num_reads);
        std::vector<uint8_t> overflow_flags(num_reads);  // Use uint8_t instead of bool
        
        cudaMemcpy(kmer_counts.data(), d_kmer_tracker->kmer_counts,
                   num_reads * sizeof(uint16_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(overflow_flags.data(), d_kmer_tracker->tracking_overflow_flags,
                   num_reads * sizeof(bool), cudaMemcpyDeviceToHost);
        
        // Retrieve actual k-mer data
        uint32_t max_kmers_total = max_reads * enhanced_params.forward_compat.max_kmers_to_track;
        
        std::vector<uint64_t> all_minimizers(max_kmers_total);
        std::vector<uint16_t> all_positions(max_kmers_total);
        std::vector<uint32_t> all_taxa(max_kmers_total);
        std::vector<float> all_weights(max_kmers_total);
        
        cudaMemcpy(all_minimizers.data(), d_kmer_tracker->minimizer_hashes,
                   max_kmers_total * sizeof(uint64_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(all_positions.data(), d_kmer_tracker->read_positions,
                   max_kmers_total * sizeof(uint16_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(all_taxa.data(), d_kmer_tracker->taxa_assignments,
                   max_kmers_total * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(all_weights.data(), d_kmer_tracker->confidence_weights,
                   max_kmers_total * sizeof(float), cudaMemcpyDeviceToHost);
        
        // Populate result structures
        for (int read_id = 0; read_id < num_reads; read_id++) {
            auto& result = results[read_id];
            auto& tracking = result.kmer_tracking;
            
            tracking.num_tracked_kmers = kmer_counts[read_id];
            tracking.tracking_overflow = (bool)overflow_flags[read_id];
            
            // Extract k-mer data for this read
            uint32_t start_idx = read_id * enhanced_params.forward_compat.max_kmers_to_track;
            uint32_t count = tracking.num_tracked_kmers;
            
            tracking.minimizer_hashes.assign(
                all_minimizers.begin() + start_idx,
                all_minimizers.begin() + start_idx + count
            );
            tracking.read_positions.assign(
                all_positions.begin() + start_idx,
                all_positions.begin() + start_idx + count
            );
            tracking.contributing_taxa.assign(
                all_taxa.begin() + start_idx,
                all_taxa.begin() + start_idx + count
            );
            tracking.confidence_weights.assign(
                all_weights.begin() + start_idx,
                all_weights.begin() + start_idx + count
            );
        }
    }
    
    void prepare_reads_batch(const std::vector<std::string>& reads) {
        // Transfer reads to GPU
        std::vector<uint32_t> offsets(reads.size() + 1);
        std::vector<uint32_t> lengths(reads.size());
        
        // Calculate offsets and total size
        uint32_t total_size = 0;
        for (size_t i = 0; i < reads.size(); i++) {
            offsets[i] = total_size;
            lengths[i] = reads[i].length();
            total_size += lengths[i];
        }
        offsets[reads.size()] = total_size;
        
        // Concatenate all sequences
        std::string concatenated;
        concatenated.reserve(total_size);
        for (const auto& read : reads) {
            concatenated += read;
        }
        
        // Transfer to GPU
        cudaMemcpy(d_sequences, concatenated.c_str(), total_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, offsets.data(), offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, lengths.data(), lengths.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    }
    
    std::vector<Phase1ClassificationResult> convert_base_to_enhanced_results(
        const std::vector<PairedReadClassification>& base_results) {
        std::vector<Phase1ClassificationResult> enhanced_results;
        enhanced_results.reserve(base_results.size());
        
        for (const auto& base : base_results) {
            Phase1ClassificationResult enhanced;
            enhanced.taxon_id = base.taxon_id;
            enhanced.primary_confidence = base.confidence_score;
            enhanced.secondary_confidence = base.confidence_score;
            enhanced.passed_stage1 = base.taxon_id != 0;
            enhanced.passed_stage2 = base.taxon_id != 0 && base.confidence_score >= enhanced_params.secondary_confidence_threshold;
            enhanced.is_high_confidence_call = enhanced.passed_stage2;
            enhanced.total_kmers = base.read1_kmers;
            enhanced.classified_kmers = base.read1_votes;
            enhanced.distinct_minimizers_found = 0;  // Not available from base
            enhanced.phylogenetic_consistency_score = 1.0f;
            enhanced.passed_phylogenetic_filter = true;
            enhanced.has_taxonomic_conflicts = false;
            enhanced.multiple_minimizer_requirement_met = true;
            enhanced.stage_agreement_score = 1.0f;
            enhanced.primary_supporting_votes = base.read1_votes;
            enhanced.secondary_supporting_votes = base.read1_votes;
            
            enhanced_results.push_back(enhanced);
        }
        
        return enhanced_results;
    }
    
    bool load_depth_lookup(const std::string& taxonomy_file) {
        // This is a simplified version - you'll need to implement proper depth calculation
        // from the taxonomy file
        if (num_taxonomy_nodes == 0) return false;
        
        // Allocate depth lookup table
        std::vector<uint32_t> depths(num_taxonomy_nodes, 0);
        
        // TODO: Load taxonomy and calculate depths
        // For now, just initialize with dummy values
        for (uint32_t i = 0; i < num_taxonomy_nodes; i++) {
            depths[i] = i % 10;  // Placeholder
        }
        
        // Allocate and transfer to GPU
        cudaMalloc(&d_depth_lookup, num_taxonomy_nodes * sizeof(uint32_t));
        cudaMemcpy(d_depth_lookup, depths.data(), num_taxonomy_nodes * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        return true;
    }
};

// ================================================================
// MISSING KERNEL IMPLEMENTATIONS
// ================================================================

// Fast enhanced classification kernel (called from fast_enhanced_classifier.h)
__global__ void fast_enhanced_classification_kernel(
    const char* sequences,
    const uint32_t* read_offsets,
    const uint32_t* read_lengths,
    const GPUCompactHashTable* hash_table,
    const TaxonHashTable* taxonomy_table,
    const PhyloDistanceCache* distance_cache,
    GPUKmerTracker* kmer_tracker,
    Phase1ClassificationResult* results,
    int num_reads,
    FastEnhancedParams params) {
    
    int read_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (read_id >= num_reads) return;
    
    const char* sequence = sequences + read_offsets[read_id];
    uint32_t seq_length = read_lengths[read_id];
    
    // Initialize result using your existing structure
    Phase1ClassificationResult& result = results[read_id];
    result.taxon_id = 0;
    result.primary_confidence = 0.0f;
    result.secondary_confidence = 0.0f;
    result.passed_stage1 = false;
    result.passed_stage2 = false;
    result.is_high_confidence_call = false;
    result.total_kmers = 0;
    result.classified_kmers = 0;
    result.distinct_minimizers_found = 0;
    result.phylogenetic_consistency_score = 0.0f;
    result.passed_phylogenetic_filter = false;
    result.has_taxonomic_conflicts = false;
    result.multiple_minimizer_requirement_met = false;
    
    if (seq_length < params.k) return;
    
    // Use your existing logic but adapt parameters
    Phase1EnhancedParams adapted_params;
    adapted_params.k = params.k;
    adapted_params.ell = params.ell;
    adapted_params.spaces = params.spaces;
    adapted_params.primary_confidence_threshold = params.primary_confidence_threshold;
    adapted_params.secondary_confidence_threshold = params.secondary_confidence_threshold;
    adapted_params.enable_phylogenetic_validation = params.enable_phylogenetic_validation;
    adapted_params.min_distinct_minimizers = params.min_distinct_minimizers;
    adapted_params.forward_compat.track_kmer_positions = params.forward_compat.track_kmer_positions;
    adapted_params.forward_compat.max_kmers_to_track = params.forward_compat.max_kmers_to_track;
    
    // Reuse your existing kernel logic from phase1_enhanced_classification_kernel
    // but with compact taxonomy support
    
    // Local arrays for k-mer processing (from your implementation)
    const int MAX_LOCAL_TAXA = 32;
    uint32_t hit_taxa[MAX_LOCAL_TAXA];
    uint32_t hit_votes[MAX_LOCAL_TAXA];
    int num_taxa = 0;
    
    uint32_t total_kmers = seq_length - params.k + 1;
    result.total_kmers = total_kmers;
    
    // Extract minimizers using your existing logic
    for (int pos = 0; pos < total_kmers; pos++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, pos, params.k, params.ell, params.spaces, 0
        );
        
        if (minimizer == UINT64_MAX) continue;
        
        // Look up in hash table
        uint32_t lca = lookup_lca_gpu(hash_table, minimizer);
        if (lca == 0) continue;
        
        result.classified_kmers++;
        
        // Add vote (from your implementation)
        bool found = false;
        for (int i = 0; i < num_taxa; i++) {
            if (hit_taxa[i] == lca) {
                hit_votes[i]++;
                found = true;
                break;
            }
        }
        
        if (!found && num_taxa < MAX_LOCAL_TAXA) {
            hit_taxa[num_taxa] = lca;
            hit_votes[num_taxa] = 1;
            num_taxa++;
        }
    }
    
    // Find best classification (from your implementation)
    uint32_t best_taxon = 0;
    uint32_t best_votes = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (hit_votes[i] > best_votes) {
            best_votes = hit_votes[i];
            best_taxon = hit_taxa[i];
        }
    }
    
    // Apply your existing confidence calculation logic
    result.primary_confidence = total_kmers > 0 ? (float)best_votes / total_kmers : 0.0f;
    result.passed_stage1 = (result.primary_confidence >= params.primary_confidence_threshold);
    
    if (result.passed_stage1 && result.classified_kmers > 0) {
        result.secondary_confidence = (float)best_votes / result.classified_kmers;
        result.passed_stage2 = (result.secondary_confidence >= params.secondary_confidence_threshold);
        
        if (result.passed_stage2) {
            result.taxon_id = best_taxon;
            result.is_high_confidence_call = true;
        }
    }
}

// Phylogenetic enhanced classification kernel (called from phase1_enhanced_classifier_with_phylo.h)
__global__ void phase1_enhanced_classification_with_phylo_kernel(
    const char* sequences,
    const uint32_t* read_offsets,
    const uint32_t* read_lengths,
    const GPUCompactHashTable* hash_table,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    const uint8_t* rank_lookup,
    uint32_t max_taxon_id,
    GPUKmerTracker* kmer_tracker,
    Phase1ClassificationResult* results,
    int num_reads,
    Phase1EnhancedParams params,
    PhylogeneticClassificationParams phylo_params) {
    
    // This can delegate to your existing phase1_enhanced_classification_kernel
    // since it has the same core logic
    
    // Convert phylo_params to match your existing parameter structure
    Phase1EnhancedParams adapted_params = params;
    adapted_params.enable_phylogenetic_validation = phylo_params.use_phylogenetic_validation;
    adapted_params.max_phylogenetic_distance = phylo_params.max_phylogenetic_distance;
    
    // Call your existing kernel implementation
    phase1_enhanced_classification_kernel<<<1, 1>>>(
        sequences, read_offsets, read_lengths,
        hash_table, taxonomy_tree, parent_lookup,
        (uint32_t*)depth_lookup, // Cast for compatibility
        kmer_tracker, results, num_reads, adapted_params
    );
}

// Helper device functions that might be missing
__device__ uint32_t find_lca_gpu_simple(uint32_t taxon1, uint32_t taxon2, 
                                        const uint32_t* parent_lookup, 
                                        const uint32_t* depth_lookup, 
                                        uint32_t max_taxon_id) {
    if (taxon1 == taxon2) return taxon1;
    if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return 1;
    
    // Simple LCA implementation
    uint32_t current1 = taxon1;
    uint32_t current2 = taxon2;
    
    // Walk both up until they meet or reach root
    for (int steps = 0; steps < 50; steps++) {
        if (current1 == current2) return current1;
        if (current1 == 0 || current1 == 1) break;
        if (current2 == 0 || current2 == 1) break;
        
        if (current1 > current2) {
            current1 = parent_lookup[current1];
        } else {
            current2 = parent_lookup[current2];
        }
    }
    
    return 1; // Root fallback
}

#endif // PHASE1_ENHANCED_CLASSIFICATION_CU