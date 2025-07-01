// phase1_enhanced_classification_fixed.cu
// FIXED: Removed dynamic parallelism - all kernels launched from host only

#ifndef PHASE1_ENHANCED_CLASSIFICATION_FIXED_CU
#define PHASE1_ENHANCED_CLASSIFICATION_FIXED_CU

#include <cuda_runtime.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

// Include headers (not .cu files)
#include "gpu_kraken_classifier.h"
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

// FIXED: Main classification kernel - NO dynamic parallelism
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
        // Extract minimizer using the function from gpu_kraken_classifier.h
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
        
        // Look up in hash table using function from gpu_kraken_classifier.h
        uint32_t lca = lookup_lca_gpu_impl(hash_table, minimizer);
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
    } else if (result.passed_stage1) {
        result.taxon_id = best_taxon_stage1;
        result.is_high_confidence_call = false;
    } else {
        result.taxon_id = 0;  // Unclassified
    }
    
    result.secondary_supporting_votes = best_votes_stage1;
    result.stage_agreement_score = result.passed_stage2 ? 1.0f : 
                                  (result.secondary_confidence / params.secondary_confidence_threshold);
}

// FIXED: Fast enhanced classification kernel - NO dynamic parallelism
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
    
    // Convert FastEnhancedParams to Phase1EnhancedParams for compatibility
    Phase1EnhancedParams adapted_params;
    adapted_params.k = params.k;
    adapted_params.ell = params.ell;
    adapted_params.spaces = params.spaces;
    adapted_params.primary_confidence_threshold = params.primary_confidence_threshold;
    adapted_params.secondary_confidence_threshold = params.secondary_confidence_threshold;
    adapted_params.min_kmers_for_classification = params.min_kmers_for_classification;
    adapted_params.enable_phylogenetic_validation = params.enable_phylogenetic_validation;
    adapted_params.max_phylogenetic_distance = params.max_phylogenetic_distance;
    adapted_params.require_multiple_minimizers = params.require_multiple_minimizers;
    adapted_params.min_distinct_minimizers = params.min_distinct_minimizers;
    adapted_params.weighted_phylo = params.weighted_phylo;
    adapted_params.forward_compat = params.forward_compat;
    
    // Process single read using adapted parameters
    // (Same logic as main kernel but for single read)
    
    const char* read_seq = sequences + read_offsets[read_id];
    uint32_t read_length = read_lengths[read_id];
    
    Phase1ClassificationResult& result = results[read_id];
    
    // Initialize result
    result.taxon_id = 0;
    result.primary_confidence = 0.0f;
    result.secondary_confidence = 0.0f;
    result.passed_stage1 = false;
    result.passed_stage2 = false;
    
    if (read_length < adapted_params.k) return;
    
    // Simplified processing for fast mode
    int total_kmers = read_length - adapted_params.k + 1;
    int classified_kmers = 0;
    uint32_t best_taxon = 0;
    uint32_t best_votes = 0;
    
    const int MAX_TAXA = 16;  // Reduced for fast mode
    uint32_t taxa[MAX_TAXA];
    uint32_t votes[MAX_TAXA];
    int num_taxa = 0;
    
    for (int pos = 0; pos < total_kmers; pos++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            read_seq, pos, adapted_params.k, adapted_params.ell, adapted_params.spaces, 0
        );
        
        if (minimizer == UINT64_MAX) continue;
        
        uint32_t lca = lookup_lca_gpu_impl(hash_table, minimizer);
        if (lca == 0) continue;
        
        classified_kmers++;
        
        // Find or add taxon
        bool found = false;
        for (int i = 0; i < num_taxa; i++) {
            if (taxa[i] == lca) {
                votes[i]++;
                if (votes[i] > best_votes) {
                    best_votes = votes[i];
                    best_taxon = lca;
                }
                found = true;
                break;
            }
        }
        
        if (!found && num_taxa < MAX_TAXA) {
            taxa[num_taxa] = lca;
            votes[num_taxa] = 1;
            if (1 > best_votes) {
                best_votes = 1;
                best_taxon = lca;
            }
            num_taxa++;
        }
    }
    
    result.total_kmers = total_kmers;
    result.classified_kmers = classified_kmers;
    result.primary_confidence = total_kmers > 0 ? (float)best_votes / total_kmers : 0.0f;
    result.secondary_confidence = result.primary_confidence;  // Simplified for fast mode
    
    result.passed_stage1 = result.primary_confidence >= adapted_params.primary_confidence_threshold;
    result.passed_stage2 = result.secondary_confidence >= adapted_params.secondary_confidence_threshold;
    
    if (result.passed_stage1) {
        result.taxon_id = best_taxon;
        result.is_high_confidence_call = result.passed_stage2;
    }
}

// FIXED: Enhanced classification with phylogeny kernel - NO dynamic parallelism
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
    
    int read_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (read_id >= num_reads) return;
    
    // For now, use the standard enhanced classification
    // The phylogenetic enhancement would be added here
    
    const char* read_seq = sequences + read_offsets[read_id];
    uint32_t read_length = read_lengths[read_id];
    
    Phase1ClassificationResult& result = results[read_id];
    
    // Initialize result
    result.taxon_id = 0;
    result.primary_confidence = 0.0f;
    result.secondary_confidence = 0.0f;
    result.passed_stage1 = false;
    result.passed_stage2 = false;
    result.phylogenetic_consistency_score = 1.0f;
    result.passed_phylogenetic_filter = true;
    
    if (read_length < params.k) return;
    
    // Simplified implementation - would be enhanced with phylogenetic validation
    // For now, just mark as processed
    result.total_kmers = read_length - params.k + 1;
    result.classified_kmers = 0;  // Would be calculated
    result.primary_confidence = 0.0f;  // Would be calculated
    result.secondary_confidence = 0.0f;  // Would be calculated
}

// Device function implementations (same as before)
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
    int max_depth = max(depth1, depth2);
    
    return max_depth > 0 ? (float)total_distance / (2.0f * max_depth) : 1.0f;
}

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
    
    return 1;  // Root if no common ancestor found
}

#endif // PHASE1_ENHANCED_CLASSIFICATION_FIXED_CU
