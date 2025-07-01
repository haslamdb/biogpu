// ncbi_taxonomy_loader.h
// Load and process NCBI taxonomy files for phylogenetic analysis

#ifndef NCBI_TAXONOMY_LOADER_H
#define NCBI_TAXONOMY_LOADER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace BioGPU {
namespace Taxonomy {

// NCBI taxonomy node structure
struct NCBITaxonNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    std::string rank;           // species, genus, family, etc.
    std::string name;           // scientific name
    uint8_t depth;              // distance from root
    
    // Rank encoding for fast comparison
    uint8_t rank_level;         // 0=root, 1=superkingdom, 2=phylum, etc.
};

// Taxonomic rank hierarchy (NCBI standard)
enum class TaxonomicRank : uint8_t {
    ROOT = 0,
    SUPERKINGDOM = 1,
    KINGDOM = 2,
    PHYLUM = 3,
    CLASS = 4,
    ORDER = 5,
    FAMILY = 6,
    GENUS = 7,
    SPECIES = 8,
    SUBSPECIES = 9,
    STRAIN = 10,
    NO_RANK = 255
};

class NCBITaxonomyLoader {
private:
    std::unordered_map<uint32_t, NCBITaxonNode> taxonomy_nodes;
    std::unordered_map<std::string, TaxonomicRank> rank_map;
    std::vector<uint32_t> parent_lookup;      // Fast parent lookup by taxon_id
    std::vector<uint8_t> depth_lookup;        // Fast depth lookup by taxon_id
    std::vector<uint8_t> rank_lookup;         // Fast rank lookup by taxon_id
    uint32_t max_taxon_id = 0;
    
    void initialize_rank_map() {
        rank_map["root"] = TaxonomicRank::ROOT;
        rank_map["superkingdom"] = TaxonomicRank::SUPERKINGDOM;
        rank_map["kingdom"] = TaxonomicRank::KINGDOM;
        rank_map["phylum"] = TaxonomicRank::PHYLUM;
        rank_map["class"] = TaxonomicRank::CLASS;
        rank_map["order"] = TaxonomicRank::ORDER;
        rank_map["family"] = TaxonomicRank::FAMILY;
        rank_map["genus"] = TaxonomicRank::GENUS;
        rank_map["species"] = TaxonomicRank::SPECIES;
        rank_map["subspecies"] = TaxonomicRank::SUBSPECIES;
        rank_map["strain"] = TaxonomicRank::STRAIN;
        rank_map["no rank"] = TaxonomicRank::NO_RANK;
    }
    
    TaxonomicRank string_to_rank(const std::string& rank_str) {
        auto it = rank_map.find(rank_str);
        return (it != rank_map.end()) ? it->second : TaxonomicRank::NO_RANK;
    }
    
public:
    NCBITaxonomyLoader() {
        initialize_rank_map();
    }
    
    // Load NCBI taxonomy from nodes.dmp and names.dmp
    bool load_ncbi_taxonomy(const std::string& nodes_dmp_path, 
                           const std::string& names_dmp_path) {
        
        std::cout << "Loading NCBI taxonomy..." << std::endl;
        
        // Load nodes.dmp first
        if (!load_nodes_dmp(nodes_dmp_path)) {
            return false;
        }
        
        // Load names.dmp
        if (!load_names_dmp(names_dmp_path)) {
            return false;
        }
        
        // Build lookup tables
        build_lookup_tables();
        
        // Calculate depths
        calculate_depths();
        
        std::cout << "Taxonomy loaded: " << taxonomy_nodes.size() << " nodes, "
                  << "max taxon ID: " << max_taxon_id << std::endl;
        
        return true;
    }
    
    // Load nodes.dmp (taxonomy structure)
    bool load_nodes_dmp(const std::string& nodes_dmp_path) {
        std::ifstream file(nodes_dmp_path);
        if (!file.is_open()) {
            std::cerr << "Cannot open nodes.dmp: " << nodes_dmp_path << std::endl;
            return false;
        }
        
        std::string line;
        int nodes_loaded = 0;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            // Parse nodes.dmp format: taxon_id | parent_id | rank | ...
            std::vector<std::string> fields = split_dmp_line(line);
            
            if (fields.size() >= 3) {
                try {
                    uint32_t taxon_id = std::stoul(fields[0]);
                    uint32_t parent_id = std::stoul(fields[1]);
                    std::string rank = fields[2];
                    
                    NCBITaxonNode node;
                    node.taxon_id = taxon_id;
                    node.parent_id = parent_id;
                    node.rank = rank;
                    node.rank_level = static_cast<uint8_t>(string_to_rank(rank));
                    node.depth = 0;  // Will be calculated later
                    
                    taxonomy_nodes[taxon_id] = node;
                    max_taxon_id = std::max(max_taxon_id, taxon_id);
                    nodes_loaded++;
                    
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing nodes.dmp line: " << line << std::endl;
                    continue;
                }
            }
        }
        
        file.close();
        std::cout << "Loaded " << nodes_loaded << " taxonomy nodes" << std::endl;
        return nodes_loaded > 0;
    }
    
    // Load names.dmp (scientific names)
    bool load_names_dmp(const std::string& names_dmp_path) {
        std::ifstream file(names_dmp_path);
        if (!file.is_open()) {
            std::cerr << "Cannot open names.dmp: " << names_dmp_path << std::endl;
            return false;
        }
        
        std::string line;
        int names_loaded = 0;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            // Parse names.dmp format: taxon_id | name | unique_name | name_class |
            std::vector<std::string> fields = split_dmp_line(line);
            
            if (fields.size() >= 4) {
                try {
                    uint32_t taxon_id = std::stoul(fields[0]);
                    std::string name = fields[1];
                    std::string name_class = fields[3];
                    
                    // Only use scientific names
                    if (name_class == "scientific name") {
                        auto it = taxonomy_nodes.find(taxon_id);
                        if (it != taxonomy_nodes.end()) {
                            it->second.name = name;
                            names_loaded++;
                        }
                    }
                    
                } catch (const std::exception& e) {
                    continue;
                }
            }
        }
        
        file.close();
        std::cout << "Loaded " << names_loaded << " scientific names" << std::endl;
        return names_loaded > 0;
    }
    
    // Build fast lookup tables for GPU
    void build_lookup_tables() {
        parent_lookup.resize(max_taxon_id + 1, 0);
        depth_lookup.resize(max_taxon_id + 1, 0);
        rank_lookup.resize(max_taxon_id + 1, static_cast<uint8_t>(TaxonomicRank::NO_RANK));
        
        for (const auto& [taxon_id, node] : taxonomy_nodes) {
            if (taxon_id <= max_taxon_id) {
                parent_lookup[taxon_id] = node.parent_id;
                rank_lookup[taxon_id] = node.rank_level;
            }
        }
        
        std::cout << "Built lookup tables for " << max_taxon_id + 1 << " taxa" << std::endl;
    }
    
    // Calculate depth from root for each taxon
    void calculate_depths() {
        // Use iterative approach to avoid stack overflow
        std::vector<bool> depth_calculated(max_taxon_id + 1, false);
        
        // Root node has depth 0
        if (taxonomy_nodes.find(1) != taxonomy_nodes.end()) {
            depth_lookup[1] = 0;
            depth_calculated[1] = true;
            taxonomy_nodes[1].depth = 0;
        }
        
        bool changed = true;
        int iterations = 0;
        
        while (changed && iterations < 50) {  // Prevent infinite loops
            changed = false;
            iterations++;
            
            for (auto& [taxon_id, node] : taxonomy_nodes) {
                if (!depth_calculated[taxon_id] && taxon_id <= max_taxon_id) {
                    uint32_t parent_id = node.parent_id;
                    
                    // Check if parent depth is known
                    if (parent_id <= max_taxon_id && depth_calculated[parent_id]) {
                        uint8_t depth = depth_lookup[parent_id] + 1;
                        depth_lookup[taxon_id] = depth;
                        node.depth = depth;
                        depth_calculated[taxon_id] = true;
                        changed = true;
                    }
                }
            }
        }
        
        std::cout << "Calculated depths in " << iterations << " iterations" << std::endl;
    }
    
    // Get parent lookup table for GPU
    const std::vector<uint32_t>& get_parent_lookup() const {
        return parent_lookup;
    }
    
    // Get depth lookup table for GPU  
    const std::vector<uint8_t>& get_depth_lookup() const {
        return depth_lookup;
    }
    
    // Get rank lookup table for GPU
    const std::vector<uint8_t>& get_rank_lookup() const {
        return rank_lookup;
    }
    
    // Get max taxon ID
    uint32_t get_max_taxon_id() const {
        return max_taxon_id;
    }
    
    // Find LCA of two taxa (host version)
    uint32_t find_lca(uint32_t taxon1, uint32_t taxon2) const {
        if (taxon1 == taxon2) return taxon1;
        if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return 1;  // Root
        
        // Get depths
        uint8_t depth1 = depth_lookup[taxon1];
        uint8_t depth2 = depth_lookup[taxon2];
        
        // Bring both to same depth
        uint32_t current1 = taxon1;
        uint32_t current2 = taxon2;
        
        while (depth1 > depth2 && current1 != 0) {
            current1 = parent_lookup[current1];
            depth1--;
        }
        
        while (depth2 > depth1 && current2 != 0) {
            current2 = parent_lookup[current2];
            depth2--;
        }
        
        // Find common ancestor
        for (int steps = 0; steps < 100; steps++) {
            if (current1 == current2) return current1;
            if (current1 == 0 || current2 == 0) break;
            
            current1 = parent_lookup[current1];
            current2 = parent_lookup[current2];
        }
        
        return 1;  // Root fallback
    }
    
    // Calculate phylogenetic distance between two taxa
    float calculate_phylogenetic_distance(uint32_t taxon1, uint32_t taxon2) const {
        if (taxon1 == taxon2) return 0.0f;
        if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return 1.0f;
        
        uint32_t lca = find_lca(taxon1, taxon2);
        
        uint8_t depth1 = depth_lookup[taxon1];
        uint8_t depth2 = depth_lookup[taxon2];
        uint8_t lca_depth = depth_lookup[lca];
        
        int total_distance = (depth1 - lca_depth) + (depth2 - lca_depth);
        
        // Normalize by maximum reasonable distance (20 levels)
        return std::min(1.0f, total_distance / 20.0f);
    }
    
    // Get taxon information
    const NCBITaxonNode* get_taxon_info(uint32_t taxon_id) const {
        auto it = taxonomy_nodes.find(taxon_id);
        return (it != taxonomy_nodes.end()) ? &it->second : nullptr;
    }
    
    // Get taxonomic lineage string
    std::string get_lineage(uint32_t taxon_id) const {
        std::vector<std::string> lineage_parts;
        uint32_t current = taxon_id;
        
        for (int steps = 0; steps < 50 && current != 0 && current != 1; steps++) {
            const NCBITaxonNode* node = get_taxon_info(current);
            if (node && !node->name.empty()) {
                lineage_parts.insert(lineage_parts.begin(), node->name);
            }
            
            if (current >= max_taxon_id) break;
            current = parent_lookup[current];
        }
        
        // Join with semicolons
        std::string result;
        for (size_t i = 0; i < lineage_parts.size(); i++) {
            if (i > 0) result += "; ";
            result += lineage_parts[i];
        }
        
        return result;
    }
    
    // Print taxonomy statistics
    void print_statistics() const {
        std::unordered_map<std::string, int> rank_counts;
        
        for (const auto& [taxon_id, node] : taxonomy_nodes) {
            rank_counts[node.rank]++;
        }
        
        std::cout << "\n=== TAXONOMY STATISTICS ===" << std::endl;
        std::cout << "Total taxa: " << taxonomy_nodes.size() << std::endl;
        std::cout << "Max taxon ID: " << max_taxon_id << std::endl;
        
        std::cout << "\nRank distribution:" << std::endl;
        for (const auto& [rank, count] : rank_counts) {
            std::cout << "  " << rank << ": " << count << std::endl;
        }
    }
    
private:
    // Parse DMP file line (pipe-delimited)
    std::vector<std::string> split_dmp_line(const std::string& line) {
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        
        while (std::getline(ss, field, '|')) {
            // Trim whitespace
            field.erase(0, field.find_first_not_of(" \t"));
            field.erase(field.find_last_not_of(" \t") + 1);
            fields.push_back(field);
        }
        
        return fields;
    }
};

} // namespace Taxonomy
} // namespace BioGPU

#endif // NCBI_TAXONOMY_LOADER_H