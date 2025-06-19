// phase1_enhanced_classifier_with_phylo.h
// Enhanced Classification Host Class with NCBI Taxonomy Integration

#ifndef PHASE1_ENHANCED_CLASSIFIER_WITH_PHYLO_H
#define PHASE1_ENHANCED_CLASSIFIER_WITH_PHYLO_H

#include "enhanced_classification_params.h"
#include "ncbi_taxonomy_loader.h"
#include "gpu_phylogenetic_calculator.cuh"
#include <vector>
#include <iostream>
#include <iomanip>
#include <memory>

using namespace BioGPU::Enhanced;
using namespace BioGPU::Taxonomy;

// Forward declaration of CUDA kernel
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
    PhylogeneticClassificationParams phylo_params);

class Phase1EnhancedClassifierWithPhylogeny {
private:
    Phase1EnhancedParams params;
    PhylogeneticClassificationParams phylo_params;
    
    // GPU memory
    GPUKmerTracker* d_kmer_tracker;
    Phase1ClassificationResult* d_results;
    GPUPhylogeneticData d_phylo_data;
    
    // Host taxonomy data
    std::unique_ptr<NCBITaxonomyLoader> taxonomy_loader;
    
    // GPU data from base classifier (would be inherited)
    char* d_sequences;
    uint32_t* d_read_offsets;
    uint32_t* d_read_lengths;
    GPUCompactHashTable* d_hash_table;
    TaxonomyNode* d_taxonomy_tree;
    
    int max_reads;
    bool tracking_enabled;
    bool phylogeny_loaded;
    
public:
    Phase1EnhancedClassifierWithPhylogeny(
        const Phase1EnhancedParams& config, 
        const PhylogeneticClassificationParams& phylo_config,
        int max_read_count = 50000) 
        : params(config), phylo_params(phylo_config), max_reads(max_read_count), 
          d_kmer_tracker(nullptr), d_results(nullptr), phylogeny_loaded(false) {
        
        tracking_enabled = params.forward_compat.track_kmer_positions;
        
        std::cout << "Phase 1 Enhanced Classifier with Phylogeny initialized" << std::endl;
        std::cout << "  Two-stage filtering: " << params.primary_confidence_threshold 
                  << " -> " << params.secondary_confidence_threshold << std::endl;
        std::cout << "  Phylogenetic validation: " << (phylo_params.use_phylogenetic_validation ? "ON" : "OFF") << std::endl;
        std::cout << "  K-mer tracking: " << (tracking_enabled ? "ON" : "OFF") << std::endl;
        
        // Initialize phylogenetic data
        if (phylo_params.use_phylogenetic_validation) {
            load_phylogenetic_data();
        }
        
        allocate_gpu_memory();
    }
    
    ~Phase1EnhancedClassifierWithPhylogeny() {
        free_gpu_memory();
        if (phylogeny_loaded) {
            free_phylogenetic_data_gpu(&d_phylo_data);
        }
    }
    
    // Load NCBI taxonomy data
    bool load_phylogenetic_data() {
        if (phylo_params.taxonomy_nodes_path.empty() || phylo_params.taxonomy_names_path.empty()) {
            std::cerr << "Error: Taxonomy file paths not specified" << std::endl;
            std::cerr << "Please set taxonomy_nodes_path and taxonomy_names_path" << std::endl;
            return false;
        }
        
        std::cout << "Loading NCBI taxonomy for phylogenetic analysis..." << std::endl;
        
        taxonomy_loader = std::make_unique<NCBITaxonomyLoader>();
        
        if (!taxonomy_loader->load_ncbi_taxonomy(
                phylo_params.taxonomy_nodes_path,
                phylo_params.taxonomy_names_path)) {
            std::cerr << "Failed to load NCBI taxonomy files" << std::endl;
            return false;
        }
        
        // Copy phylogenetic data to GPU
        if (!copy_phylogenetic_data_to_gpu(
                taxonomy_loader->get_parent_lookup(),
                taxonomy_loader->get_depth_lookup(),
                taxonomy_loader->get_rank_lookup(),
                taxonomy_loader->get_max_taxon_id(),
                &d_phylo_data)) {
            std::cerr << "Failed to copy phylogenetic data to GPU" << std::endl;
            return false;
        }
        
        taxonomy_loader->print_statistics();
        phylogeny_loaded = true;
        return true;
    }
    
    // Enhanced classification with phylogenetic validation
    std::vector<Phase1ClassificationResult> classify_enhanced_with_phylogeny(
        const std::vector<std::string>& reads) {
        
        if (reads.empty()) return {};
        
        int num_reads = std::min((int)reads.size(), max_reads);
        
        // Transfer reads to GPU (reuse existing method from base class)
        transfer_reads_to_gpu(reads);
        
        // Launch enhanced classification kernel with phylogenetic data
        int blocks = (num_reads + 255) / 256;
        
        // Updated kernel call with phylogenetic parameters
        phase1_enhanced_classification_with_phylo_kernel<<<blocks, 256>>>(
            d_sequences,
            d_read_offsets,
            d_read_lengths,
            d_hash_table,
            d_taxonomy_tree,
            d_phylo_data.parent_lookup,
            d_phylo_data.depth_lookup,
            d_phylo_data.rank_lookup,
            d_phylo_data.max_taxon_id,
            d_kmer_tracker,
            d_results,
            num_reads,
            params,
            phylo_params
        );
        
        cudaDeviceSynchronize();
        
        // Check for errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error in enhanced classification: " << cudaGetErrorString(error) << std::endl;
            return {};
        }
        
        // Retrieve results
        std::vector<Phase1ClassificationResult> results(num_reads);
        cudaMemcpy(results.data(), d_results,
                   num_reads * sizeof(Phase1ClassificationResult),
                   cudaMemcpyDeviceToHost);
        
        // Retrieve k-mer tracking data if enabled
        if (tracking_enabled) {
            retrieve_kmer_tracking_data(results, num_reads);
        }
        
        // Enhance results with phylogenetic information
        enhance_results_with_lineage(results);
        
        return results;
    }
    
    // Add lineage information to results
    void enhance_results_with_lineage(std::vector<Phase1ClassificationResult>& results) {
        if (!phylogeny_loaded || !taxonomy_loader) return;
        
        for (auto& result : results) {
            if (result.taxon_id > 0) {
                // Add lineage information to classification_path
                std::string lineage = taxonomy_loader->get_lineage(result.taxon_id);
                if (!lineage.empty()) {
                    result.classification_path += " [" + lineage + "]";
                }
                
                // Get taxonomic rank information
                const NCBITaxonNode* taxon_info = taxonomy_loader->get_taxon_info(result.taxon_id);
                if (taxon_info) {
                    result.classification_path += " (rank: " + taxon_info->rank + ")";
                }
            }
        }
    }
    
    // Print enhanced statistics with phylogenetic information
    void print_enhanced_statistics_with_phylogeny(const std::vector<Phase1ClassificationResult>& results) {
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
        
        // Rank distribution
        std::unordered_map<std::string, int> rank_distribution;
        
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
                
                // Track rank distribution
                if (phylogeny_loaded && taxonomy_loader) {
                    const NCBITaxonNode* taxon_info = taxonomy_loader->get_taxon_info(result.taxon_id);
                    if (taxon_info) {
                        rank_distribution[taxon_info->rank]++;
                    }
                }
            }
        }
        
        int classified = stage2_passed > 0 ? stage2_passed : stage1_passed;
        
        std::cout << "\n=== ENHANCED CLASSIFICATION WITH PHYLOGENY STATISTICS ===" << std::endl;
        std::cout << "Total reads: " << total << std::endl;
        std::cout << "Stage 1 passed: " << stage1_passed << " (" 
                  << std::fixed << std::setprecision(1) << (100.0 * stage1_passed / total) << "%)" << std::endl;
        std::cout << "Stage 2 passed: " << stage2_passed << " (" 
                  << (100.0 * stage2_passed / total) << "%)" << std::endl;
        std::cout << "High confidence: " << high_confidence << " (" 
                  << (100.0 * high_confidence / total) << "%)" << std::endl;
        
        if (phylo_params.use_phylogenetic_validation) {
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
            if (phylo_params.use_phylogenetic_validation) {
                std::cout << "  Phylogenetic: " << (total_phylo_score / classified) << std::endl;
            }
        }
        
        // Print rank distribution
        if (!rank_distribution.empty()) {
            std::cout << "\nClassified reads by taxonomic rank:" << std::endl;
            for (const auto& [rank, count] : rank_distribution) {
                std::cout << "  " << rank << ": " << count << " (" 
                          << std::fixed << std::setprecision(1) 
                          << (100.0 * count / classified) << "%)" << std::endl;
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
    
    // Configuration methods
    void set_taxonomy_paths(const std::string& nodes_path, const std::string& names_path) {
        phylo_params.taxonomy_nodes_path = nodes_path;
        phylo_params.taxonomy_names_path = names_path;
    }
    
    void set_phylogenetic_weighting_method(PhylogeneticClassificationParams::PhyloWeightingMethod method) {
        phylo_params.weighting_method = method;
    }
    
    void set_phylogenetic_thresholds(float max_distance, float min_consistency) {
        phylo_params.max_phylogenetic_distance = max_distance;
        phylo_params.min_consistency_score = min_consistency;
    }
    
    // Load Kraken database from directory
    bool load_database(const std::string& db_dir) {
        std::cout << "Loading database from: " << db_dir << std::endl;
        
        // Load hash table
        std::string hash_table_file = db_dir + "/hash_table.k2d";
        if (!load_hash_table(hash_table_file)) {
            std::cerr << "Failed to load hash table from " << hash_table_file << std::endl;
            return false;
        }
        
        // Load phylogeny if not already loaded
        if (!phylogeny_loaded && !phylo_params.taxonomy_nodes_path.empty()) {
            if (!load_phylogenetic_data()) {
                std::cerr << "Failed to load phylogenetic data" << std::endl;
                return false;
            }
        }
        
        std::cout << "Database loaded successfully" << std::endl;
        return true;
    }
    
private:
    bool allocate_gpu_memory() {
        // Same as before - allocate results and k-mer tracking memory
        cudaError_t error = cudaMalloc(&d_results, max_reads * sizeof(Phase1ClassificationResult));
        if (error != cudaSuccess) {
            std::cerr << "Failed to allocate results memory: " << cudaGetErrorString(error) << std::endl;
            return false;
        }
        
        if (tracking_enabled) {
            d_kmer_tracker = new GPUKmerTracker();
            uint32_t max_kmers_total = max_reads * params.forward_compat.max_kmers_to_track;
            
            cudaMalloc(&d_kmer_tracker->minimizer_hashes, max_kmers_total * sizeof(uint64_t));
            cudaMalloc(&d_kmer_tracker->read_positions, max_kmers_total * sizeof(uint16_t));
            cudaMalloc(&d_kmer_tracker->taxa_assignments, max_kmers_total * sizeof(uint32_t));
            cudaMalloc(&d_kmer_tracker->confidence_weights, max_kmers_total * sizeof(float));
            cudaMalloc(&d_kmer_tracker->kmer_counts, max_reads * sizeof(uint16_t));
            cudaMalloc(&d_kmer_tracker->tracking_overflow_flags, max_reads * sizeof(bool));
            
            d_kmer_tracker->max_kmers_per_read = params.forward_compat.max_kmers_to_track;
            d_kmer_tracker->max_reads = max_reads;
        }
        
        return true;
    }
    
    void free_gpu_memory() {
        if (d_results) { 
            cudaFree(d_results); 
            d_results = nullptr; 
        }
        
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
        // Same implementation as before
        // [Implementation details omitted for brevity - same as previous version]
    }
    
    void transfer_reads_to_gpu(const std::vector<std::string>& reads) {
        // Placeholder - use existing read transfer logic from base class
        std::cout << "Transferring " << reads.size() << " reads to GPU..." << std::endl;
    }
    
    bool load_hash_table(const std::string& hash_table_file) {
        // Similar implementation to FastEnhancedClassifier
        std::ifstream hash_in(hash_table_file, std::ios::binary);
        if (!hash_in.is_open()) {
            std::cerr << "Cannot open hash table file: " << hash_table_file << std::endl;
            return false;
        }
        
        // Read header
        uint64_t table_size, num_entries;
        hash_in.read(reinterpret_cast<char*>(&table_size), sizeof(uint64_t));
        hash_in.read(reinterpret_cast<char*>(&num_entries), sizeof(uint64_t));
        
        std::cout << "Hash table size: " << table_size << ", entries: " << num_entries << std::endl;
        
        // Create compact hash table
        GPUCompactHashTable h_cht;
        h_cht.table_size = table_size;
        h_cht.hash_mask = table_size - 1;
        h_cht.lca_bits = 20;
        h_cht.hash_bits = 32 - h_cht.lca_bits;
        
        // Allocate on GPU
        cudaMalloc(&h_cht.hash_cells, table_size * sizeof(uint32_t));
        cudaMemset(h_cht.hash_cells, 0, table_size * sizeof(uint32_t));
        
        cudaMalloc(&d_hash_table, sizeof(GPUCompactHashTable));
        cudaMemcpy(d_hash_table, &h_cht, sizeof(GPUCompactHashTable), cudaMemcpyHostToDevice);
        
        // Read entries
        std::vector<uint32_t> host_hash_cells(table_size, 0);
        
        for (uint64_t i = 0; i < num_entries; i++) {
            uint64_t minimizer_hash;
            uint32_t lca_taxon, genome_count;
            float uniqueness_score;
            
            hash_in.read(reinterpret_cast<char*>(&minimizer_hash), sizeof(uint64_t));
            hash_in.read(reinterpret_cast<char*>(&lca_taxon), sizeof(uint32_t));
            hash_in.read(reinterpret_cast<char*>(&genome_count), sizeof(uint32_t));
            hash_in.read(reinterpret_cast<char*>(&uniqueness_score), sizeof(float));
            
            uint32_t pos = minimizer_hash & h_cht.hash_mask;
            host_hash_cells[pos] = lca_taxon;
        }
        
        cudaMemcpy(h_cht.hash_cells, host_hash_cells.data(), 
                   table_size * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        hash_in.close();
        return true;
    }
};

#endif // PHASE1_ENHANCED_CLASSIFIER_WITH_PHYLO_H