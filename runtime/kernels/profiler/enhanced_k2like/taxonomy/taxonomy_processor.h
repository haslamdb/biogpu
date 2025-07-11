// taxonomy/taxonomy_processor.h
// Header for taxonomy processing and phylogenetic calculations

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

// Forward declarations
namespace BioGPU {
namespace CompactTaxonomy {

class CompactGPUTaxonomy {
public:
    explicit CompactGPUTaxonomy(bool enable_cache);
    ~CompactGPUTaxonomy();
    
    // Delete copy operations
    CompactGPUTaxonomy(const CompactGPUTaxonomy&) = delete;
    CompactGPUTaxonomy& operator=(const CompactGPUTaxonomy&) = delete;
    
    // Move operations
    CompactGPUTaxonomy(CompactGPUTaxonomy&&) = default;
    CompactGPUTaxonomy& operator=(CompactGPUTaxonomy&&) = default;
    
    bool build_from_ncbi_files(const std::string& nodes_file, const std::string& names_file);
    bool load_compact_taxonomy(const std::string& compact_file);
    bool save_compact_taxonomy(const std::string& compact_file) const;
    
    uint32_t get_parent(uint32_t taxon_id) const;
    std::string get_name(uint32_t taxon_id) const;
    std::string get_rank(uint32_t taxon_id) const;
    void get_all_taxon_ids(std::vector<uint32_t>& taxon_ids) const;
    
    bool build_gpu_structures();
    void free_gpu_memory();
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

} // namespace CompactTaxonomy
} // namespace BioGPU

// Main taxonomy processor class
class EnhancedNCBITaxonomyProcessor {
public:
    EnhancedNCBITaxonomyProcessor();
    ~EnhancedNCBITaxonomyProcessor();
    
    // Loading functions
    bool load_ncbi_taxonomy(const std::string& nodes_dmp_path, const std::string& names_dmp_path);
    bool load_from_compact_file(const std::string& compact_path);
    bool is_loaded() const { return taxonomy_loaded; }
    
    // LCA computation
    uint32_t compute_lca_of_species(const std::vector<uint32_t>& species_list);
    
    // Phylogenetic calculations
    double calculate_phylogenetic_spread(const std::vector<uint32_t>& species_list, uint32_t lca_taxon);
    uint8_t calculate_distance_to_lca(uint32_t species_taxon, uint32_t lca_taxon) const;
    
    // Accessors
    BioGPU::CompactTaxonomy::CompactGPUTaxonomy* get_compact_taxonomy();
    const std::unordered_map<uint32_t, uint32_t>& get_parent_map() const { return parents; }
    const std::unordered_map<uint32_t, std::string>& get_names_map() const { return names; }
    const std::unordered_map<uint32_t, std::string>& get_ranks_map() const { return ranks; }
    
private:
    bool extract_host_lookup_tables();
    void print_taxonomy_statistics();
    std::vector<uint32_t> get_path_to_root(uint32_t taxon_id) const;
    
    std::unique_ptr<BioGPU::CompactTaxonomy::CompactGPUTaxonomy> compact_taxonomy;
    std::unordered_map<uint32_t, uint32_t> parents;
    std::unordered_map<uint32_t, std::string> names;
    std::unordered_map<uint32_t, std::string> ranks;
    std::unordered_map<uint32_t, uint8_t> node_depths;
    bool taxonomy_loaded;
    uint32_t max_taxon_id;
};

// Utility namespaces
namespace TaxonomyUtils {
    bool save_taxonomy_tsv(const std::string& output_path, 
                          const std::unordered_map<uint32_t, uint32_t>& parents,
                          const std::unordered_map<uint32_t, std::string>& names,
                          const std::unordered_map<uint32_t, std::string>& ranks);
    
    std::string format_taxonomy_line(uint32_t taxon_id, uint32_t parent_id, 
                                    const std::string& name, const std::string& rank);
    
    bool validate_species_list(const std::vector<uint32_t>& species_list);
    bool is_valid_taxon_id(uint32_t taxon_id);
}

namespace LCAAlgorithms {
    uint32_t compute_lca_pair(uint32_t taxon1, uint32_t taxon2, 
                             const std::unordered_map<uint32_t, uint32_t>& parents);
    
    uint32_t compute_lca_multiple(const std::vector<uint32_t>& taxa, 
                                 const std::unordered_map<uint32_t, uint32_t>& parents);
    
    std::vector<uint32_t> get_path_to_root(uint32_t taxon_id, 
                                          const std::unordered_map<uint32_t, uint32_t>& parents);
    
    uint32_t compute_lca_with_depths(const std::vector<uint32_t>& taxa,
                                    const std::unordered_map<uint32_t, uint32_t>& parents,
                                    const std::unordered_map<uint32_t, uint8_t>& depths);
    
    uint32_t find_common_ancestor_from_paths(const std::vector<std::vector<uint32_t>>& paths);
}

namespace PhylogeneticUtils {
    uint8_t calculate_taxonomic_distance(uint32_t taxon1, uint32_t taxon2, 
                                        const std::unordered_map<uint32_t, uint32_t>& parents);
    
    double calculate_phylogenetic_diversity(const std::vector<uint32_t>& species_list,
                                           const EnhancedNCBITaxonomyProcessor& taxonomy);
    
    std::vector<uint32_t> extract_unique_species(const std::vector<uint32_t>& taxa_list);
    
    std::unordered_map<uint32_t, uint16_t> count_genomes_per_species(
        const std::vector<struct PhylogeneticLCACandidate>& candidates);
}