// taxonomy/taxonomy_processor.h
// NCBI taxonomy processing and phylogenetic calculations
// Handles taxonomy loading, LCA computation, and species tracking

#ifndef TAXONOMY_PROCESSOR_H
#define TAXONOMY_PROCESSOR_H

#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <memory>
#include <cstdint>
#include "../gpu_kraken_types.h"

// Forward declaration for compact taxonomy
namespace BioGPU {
namespace CompactTaxonomy {
// Dummy placeholder for CompactGPUTaxonomy to avoid incomplete type issues
class CompactGPUTaxonomy {
public:
    CompactGPUTaxonomy() = default;
    ~CompactGPUTaxonomy() = default;
};
}
}

// Enhanced NCBI taxonomy processor with phylogenetic capabilities
class EnhancedNCBITaxonomyProcessor {
private:
    std::unique_ptr<BioGPU::CompactTaxonomy::CompactGPUTaxonomy> compact_taxonomy;
    bool taxonomy_loaded;
    
    // Host-side lookup tables for database building
    std::unordered_map<uint32_t, uint32_t> parent_lookup;
    std::unordered_map<uint32_t, std::string> name_lookup;
    std::unordered_map<uint32_t, std::string> rank_lookup;
    std::unordered_map<uint32_t, uint8_t> depth_lookup;
    uint32_t max_taxon_id;
    
public:
    EnhancedNCBITaxonomyProcessor();
    ~EnhancedNCBITaxonomyProcessor();
    
    // Loading methods
    bool load_ncbi_taxonomy(const std::string& nodes_dmp_path, const std::string& names_dmp_path);
    bool load_from_compact_file(const std::string& compact_taxonomy_path);
    
    // Phylogenetic computation methods
    uint32_t compute_lca_of_species(const std::vector<uint32_t>& species_list);
    uint8_t calculate_distance_to_lca(uint32_t taxon, uint32_t lca);
    uint8_t calculate_phylogenetic_spread(const std::vector<uint32_t>& species_list, uint32_t lca);
    
    // Taxonomy information access
    std::string get_scientific_name(uint32_t taxon_id);
    std::string get_rank(uint32_t taxon_id);
    uint32_t get_parent(uint32_t taxon_id);
    uint8_t get_depth(uint32_t taxon_id);
    
    // Status and access
    bool is_loaded() const;
    BioGPU::CompactTaxonomy::CompactGPUTaxonomy* get_compact_taxonomy();
    const std::unordered_map<uint32_t, uint32_t>& get_parent_lookup() const { return parent_lookup; }
    const std::unordered_map<uint32_t, std::string>& get_name_lookup() const { return name_lookup; }
    
    // Validation and statistics
    size_t get_taxonomy_size() const;
    bool validate_taxon_id(uint32_t taxon_id) const;
    void print_taxonomy_statistics() const;
    
private:
    // Internal methods
    bool extract_host_lookup_tables();
    uint32_t find_lca_pair(uint32_t taxon1, uint32_t taxon2);
    std::vector<uint32_t> get_path_to_root(uint32_t taxon_id);
    uint8_t get_rank_weight(uint32_t taxon_id);
    bool parse_nodes_dmp(const std::string& nodes_file);
    bool parse_names_dmp(const std::string& names_file);
    void build_depth_lookup();
};

// Simple taxonomy processor for basic operations without compact taxonomy
class SimpleTaxonomyProcessor {
private:
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, std::string> taxon_ranks;
    bool taxonomy_loaded;
    
public:
    SimpleTaxonomyProcessor();
    
    // Basic loading from NCBI files
    bool load_ncbi_files(const std::string& nodes_dmp_path, const std::string& names_dmp_path);
    bool load_taxonomy_tsv(const std::string& taxonomy_tsv_path);
    
    // Basic LCA computation
    uint32_t compute_simple_lca(uint32_t taxon1, uint32_t taxon2);
    uint32_t compute_lca_of_list(const std::vector<uint32_t>& taxon_list);
    
    // Information access
    std::string get_name(uint32_t taxon_id) const;
    std::string get_rank(uint32_t taxon_id) const;
    uint32_t get_parent(uint32_t taxon_id) const;
    
    // Status
    bool is_loaded() const { return taxonomy_loaded; }
    size_t size() const { return taxon_parents.size(); }
    
    // Export functionality
    bool save_taxonomy_tsv(const std::string& output_path) const;
    void add_taxon(uint32_t taxon_id, uint32_t parent_id, const std::string& name, const std::string& rank = "no rank");
    
private:
    std::vector<uint32_t> get_lineage(uint32_t taxon_id) const;
    bool validate_taxonomy_consistency() const;
};

// Phylogenetic utilities and helper functions
namespace PhylogeneticUtils {
    // Distance calculations
    uint8_t calculate_taxonomic_distance(uint32_t taxon1, uint32_t taxon2, 
                                        const std::unordered_map<uint32_t, uint32_t>& parents);
    
    // Diversity metrics
    double calculate_phylogenetic_diversity(const std::vector<uint32_t>& species_list,
                                           const EnhancedNCBITaxonomyProcessor& taxonomy);
    
    uint8_t calculate_taxonomic_spread(const std::vector<uint32_t>& taxa_list, uint32_t lca,
                                      const std::unordered_map<uint32_t, uint32_t>& parents);
    
    // Validation utilities
    bool validate_species_list(const std::vector<uint32_t>& species_list);
    bool is_valid_taxon_id(uint32_t taxon_id);
    
    // Conversion utilities
    std::vector<uint32_t> extract_species_from_candidates(const std::vector<PhylogeneticLCACandidate>& candidates);
    std::unordered_map<uint32_t, uint16_t> count_genomes_per_species(const std::vector<PhylogeneticLCACandidate>& candidates);
    
    // File format utilities
    bool parse_taxonomy_line(const std::string& line, uint32_t& taxon_id, uint32_t& parent_id, 
                            std::string& name, std::string& rank);
    std::string format_taxonomy_line(uint32_t taxon_id, uint32_t parent_id, 
                                    const std::string& name, const std::string& rank);
}

// LCA computation algorithms
namespace LCAAlgorithms {
    // Simple LCA for two taxa
    uint32_t compute_lca_pair(uint32_t taxon1, uint32_t taxon2, 
                             const std::unordered_map<uint32_t, uint32_t>& parents);
    
    // Optimized LCA for multiple taxa
    uint32_t compute_lca_multiple(const std::vector<uint32_t>& taxa, 
                                 const std::unordered_map<uint32_t, uint32_t>& parents);
    
    // LCA with depth optimization
    uint32_t compute_lca_with_depths(const std::vector<uint32_t>& taxa,
                                    const std::unordered_map<uint32_t, uint32_t>& parents,
                                    const std::unordered_map<uint32_t, uint8_t>& depths);
    
    // Path-based LCA computation
    std::vector<uint32_t> get_path_to_root(uint32_t taxon_id, 
                                          const std::unordered_map<uint32_t, uint32_t>& parents);
    
    uint32_t find_common_ancestor_from_paths(const std::vector<std::vector<uint32_t>>& paths);
}

#endif // TAXONOMY_PROCESSOR_H
