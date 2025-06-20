// taxonomy/taxonomy_processor.cu
// Implementation of taxonomy processing and phylogenetic calculations
// Extracted from monolithic database builder

#include "taxonomy_processor.h"
#include "../gpu_kraken_types.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>

// Include compact taxonomy if available
#ifdef COMPACT_TAXONOMY_AVAILABLE
#include "../tools/compact_gpu_taxonomy.h"
#endif

// SpeciesTrackingData, PhylogeneticLCACandidate, and ContributingTaxaArrays
// implementations are already defined as inline in gpu_kraken_types.h

// ===========================
// EnhancedNCBITaxonomyProcessor Implementation
// ===========================

EnhancedNCBITaxonomyProcessor::EnhancedNCBITaxonomyProcessor() 
    : taxonomy_loaded(false), max_taxon_id(0) {
}

EnhancedNCBITaxonomyProcessor::~EnhancedNCBITaxonomyProcessor() = default;

bool EnhancedNCBITaxonomyProcessor::load_ncbi_taxonomy(const std::string& nodes_dmp_path, const std::string& names_dmp_path) {
    std::cout << "Loading NCBI taxonomy from DMP files..." << std::endl;
    
#ifdef COMPACT_TAXONOMY_AVAILABLE
    compact_taxonomy = std::make_unique<BioGPU::CompactTaxonomy::CompactGPUTaxonomy>(false);
    
    if (!compact_taxonomy->build_from_ncbi_files(nodes_dmp_path, names_dmp_path)) {
        std::cerr << "Failed to build compact taxonomy from NCBI files" << std::endl;
        return false;
    }
    
    if (!extract_host_lookup_tables()) {
        std::cerr << "Failed to extract lookup tables from compact taxonomy" << std::endl;
        return false;
    }
#else
    // Fallback to direct parsing
    if (!parse_nodes_dmp(nodes_dmp_path) || !parse_names_dmp(names_dmp_path)) {
        return false;
    }
    build_depth_lookup();
#endif
    
    taxonomy_loaded = true;
    std::cout << "NCBI taxonomy loaded successfully" << std::endl;
    print_taxonomy_statistics();
    
    return true;
}

bool EnhancedNCBITaxonomyProcessor::load_from_compact_file(const std::string& compact_taxonomy_path) {
    std::cout << "Loading from pre-built compact taxonomy: " << compact_taxonomy_path << std::endl;
    
#ifdef COMPACT_TAXONOMY_AVAILABLE
    compact_taxonomy = std::make_unique<BioGPU::CompactTaxonomy::CompactGPUTaxonomy>(false);
    
    if (!compact_taxonomy->load_compact_taxonomy(compact_taxonomy_path)) {
        std::cerr << "Failed to load compact taxonomy file" << std::endl;
        return false;
    }
    
    if (!extract_host_lookup_tables()) {
        std::cerr << "Failed to extract lookup tables from compact taxonomy" << std::endl;
        return false;
    }
    
    taxonomy_loaded = true;
    std::cout << "Compact taxonomy loaded successfully" << std::endl;
    return true;
#else
    std::cerr << "Compact taxonomy not available in this build" << std::endl;
    return false;
#endif
}

uint32_t EnhancedNCBITaxonomyProcessor::compute_lca_of_species(const std::vector<uint32_t>& species_list) {
    if (!taxonomy_loaded || species_list.empty()) {
        return 1; // Root
    }
    
    if (species_list.size() == 1) {
        return species_list[0];
    }
    
    // Find LCA using proper taxonomy tree traversal
    uint32_t current_lca = species_list[0];
    for (size_t i = 1; i < species_list.size(); i++) {
        current_lca = find_lca_pair(current_lca, species_list[i]);
        if (current_lca == 1) break; // Already at root
    }
    
    return current_lca;
}

uint8_t EnhancedNCBITaxonomyProcessor::calculate_distance_to_lca(uint32_t taxon, uint32_t lca) {
    if (!taxonomy_loaded || taxon == lca) {
        return 0;
    }
    
    // Use depth lookup if available
    auto taxon_depth_it = depth_lookup.find(taxon);
    auto lca_depth_it = depth_lookup.find(lca);
    
    if (taxon_depth_it != depth_lookup.end() && lca_depth_it != depth_lookup.end()) {
        uint8_t taxon_depth = taxon_depth_it->second;
        uint8_t lca_depth = lca_depth_it->second;
        
        if (taxon_depth >= lca_depth) {
            return taxon_depth - lca_depth;
        }
    }
    
    // Fallback: count steps manually
    uint32_t current = taxon;
    uint8_t distance = 0;
    
    while (current != lca && current != 1 && distance < 50) {
        auto parent_it = parent_lookup.find(current);
        if (parent_it == parent_lookup.end()) break;
        
        current = parent_it->second;
        distance++;
        
        if (current == lca) break;
    }
    
    return (current == lca) ? distance : 255;
}

uint8_t EnhancedNCBITaxonomyProcessor::calculate_phylogenetic_spread(
    const std::vector<uint32_t>& species_list, uint32_t lca) {
    
    if (species_list.size() <= 1) return 0;
    
    std::vector<uint8_t> distances;
    for (uint32_t species : species_list) {
        distances.push_back(calculate_distance_to_lca(species, lca));
    }
    
    uint8_t max_dist = *std::max_element(distances.begin(), distances.end());
    uint8_t min_dist = *std::min_element(distances.begin(), distances.end());
    
    // Enhanced spread calculation using taxonomy ranks
    uint8_t range_spread = max_dist - min_dist;
    uint8_t diversity_factor = std::min((uint8_t)(species_list.size() / 5), (uint8_t)50);
    
    // Weight by taxonomic rank of LCA
    uint8_t rank_weight = get_rank_weight(lca);
    
    return std::min((uint8_t)255, (uint8_t)(range_spread * rank_weight + diversity_factor));
}

std::string EnhancedNCBITaxonomyProcessor::get_scientific_name(uint32_t taxon_id) {
    auto it = name_lookup.find(taxon_id);
    return (it != name_lookup.end()) ? it->second : ("taxon_" + std::to_string(taxon_id));
}

std::string EnhancedNCBITaxonomyProcessor::get_rank(uint32_t taxon_id) {
    auto it = rank_lookup.find(taxon_id);
    return (it != rank_lookup.end()) ? it->second : "no rank";
}

uint32_t EnhancedNCBITaxonomyProcessor::get_parent(uint32_t taxon_id) {
    auto it = parent_lookup.find(taxon_id);
    return (it != parent_lookup.end()) ? it->second : 0;
}

uint8_t EnhancedNCBITaxonomyProcessor::get_depth(uint32_t taxon_id) {
    auto it = depth_lookup.find(taxon_id);
    return (it != depth_lookup.end()) ? it->second : 0;
}

bool EnhancedNCBITaxonomyProcessor::is_loaded() const {
    return taxonomy_loaded;
}

BioGPU::CompactTaxonomy::CompactGPUTaxonomy* EnhancedNCBITaxonomyProcessor::get_compact_taxonomy() {
#ifdef COMPACT_TAXONOMY_AVAILABLE
    return compact_taxonomy.get();
#else
    return nullptr;
#endif
}

size_t EnhancedNCBITaxonomyProcessor::get_taxonomy_size() const {
    return parent_lookup.size();
}

bool EnhancedNCBITaxonomyProcessor::validate_taxon_id(uint32_t taxon_id) const {
    return parent_lookup.find(taxon_id) != parent_lookup.end();
}

void EnhancedNCBITaxonomyProcessor::print_taxonomy_statistics() const {
    std::cout << "\n=== TAXONOMY STATISTICS ===" << std::endl;
    std::cout << "Total taxa: " << parent_lookup.size() << std::endl;
    std::cout << "Named taxa: " << name_lookup.size() << std::endl;
    std::cout << "Ranked taxa: " << rank_lookup.size() << std::endl;
    std::cout << "Max taxon ID: " << max_taxon_id << std::endl;
    
    if (!depth_lookup.empty()) {
        auto max_depth_it = std::max_element(depth_lookup.begin(), depth_lookup.end(),
                                           [](const auto& a, const auto& b) { return a.second < b.second; });
        std::cout << "Maximum depth: " << (int)max_depth_it->second << std::endl;
    }
}

// Private methods
bool EnhancedNCBITaxonomyProcessor::extract_host_lookup_tables() {
#ifdef COMPACT_TAXONOMY_AVAILABLE
    if (!compact_taxonomy) {
        return false;
    }
    
    parent_lookup = compact_taxonomy->get_parent_lookup_map();
    name_lookup = compact_taxonomy->get_name_lookup_map();
    rank_lookup = compact_taxonomy->get_rank_lookup_map();
    depth_lookup = compact_taxonomy->get_depth_lookup_map();
    max_taxon_id = compact_taxonomy->get_max_taxon_id();
    
    std::cout << "Extracted host lookup tables from compact taxonomy" << std::endl;
    return true;
#else
    return false;
#endif
}

uint32_t EnhancedNCBITaxonomyProcessor::find_lca_pair(uint32_t taxon1, uint32_t taxon2) {
    if (taxon1 == taxon2) return taxon1;
    if (taxon1 == 1 || taxon2 == 1) return 1;
    
    // Get paths to root
    std::vector<uint32_t> path1 = get_path_to_root(taxon1);
    std::vector<uint32_t> path2 = get_path_to_root(taxon2);
    
    // Find first common ancestor
    std::set<uint32_t> ancestors1(path1.begin(), path1.end());
    
    for (uint32_t ancestor : path2) {
        if (ancestors1.find(ancestor) != ancestors1.end()) {
            return ancestor;
        }
    }
    
    return 1; // Root fallback
}

std::vector<uint32_t> EnhancedNCBITaxonomyProcessor::get_path_to_root(uint32_t taxon_id) {
    std::vector<uint32_t> path;
    uint32_t current = taxon_id;
    
    while (current != 1 && path.size() < 50) {
        path.push_back(current);
        auto parent_it = parent_lookup.find(current);
        if (parent_it == parent_lookup.end()) break;
        current = parent_it->second;
    }
    
    path.push_back(1); // Add root
    return path;
}

uint8_t EnhancedNCBITaxonomyProcessor::get_rank_weight(uint32_t taxon_id) {
    std::string rank = get_rank(taxon_id);
    
    // Weight factors based on taxonomic rank
    if (rank == "species" || rank == "subspecies") return 1;
    else if (rank == "genus") return 2;
    else if (rank == "family") return 3;
    else if (rank == "order") return 4;
    else if (rank == "class") return 5;
    else if (rank == "phylum") return 6;
    else if (rank == "kingdom" || rank == "superkingdom") return 7;
    else return 3; // Default for "no rank"
}

bool EnhancedNCBITaxonomyProcessor::parse_nodes_dmp(const std::string& nodes_file) {
    std::cout << "Parsing nodes.dmp: " << nodes_file << std::endl;
    
    std::ifstream file(nodes_file);
    if (!file.is_open()) {
        std::cerr << "Cannot open nodes.dmp: " << nodes_file << std::endl;
        return false;
    }
    
    std::string line;
    int nodes_loaded = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '|')) {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            fields.push_back(token);
        }
        
        if (fields.size() >= 3) {
            try {
                uint32_t taxon_id = std::stoul(fields[0]);
                uint32_t parent_id = std::stoul(fields[1]);
                std::string rank = fields[2];
                
                parent_lookup[taxon_id] = parent_id;
                rank_lookup[taxon_id] = rank;
                max_taxon_id = std::max(max_taxon_id, taxon_id);
                nodes_loaded++;
                
                // Default name
                name_lookup[taxon_id] = "taxon_" + std::to_string(taxon_id);
                
            } catch (const std::exception& e) {
                continue; // Skip invalid lines
            }
        }
    }
    
    file.close();
    std::cout << "Loaded " << nodes_loaded << " taxonomy nodes" << std::endl;
    return nodes_loaded > 0;
}

bool EnhancedNCBITaxonomyProcessor::parse_names_dmp(const std::string& names_file) {
    std::cout << "Parsing names.dmp: " << names_file << std::endl;
    
    std::ifstream file(names_file);
    if (!file.is_open()) {
        std::cerr << "Cannot open names.dmp: " << names_file << std::endl;
        return false;
    }
    
    std::string line;
    int names_loaded = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '|')) {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            fields.push_back(token);
        }
        
        if (fields.size() >= 4) {
            try {
                uint32_t taxon_id = std::stoul(fields[0]);
                std::string name_txt = fields[1];
                std::string name_class = fields[3];
                
                if (name_class == "scientific name" && parent_lookup.find(taxon_id) != parent_lookup.end()) {
                    name_lookup[taxon_id] = name_txt;
                    names_loaded++;
                }
            } catch (const std::exception& e) {
                continue; // Skip invalid lines
            }
        }
    }
    
    file.close();
    std::cout << "Loaded " << names_loaded << " scientific names" << std::endl;
    return true;
}

void EnhancedNCBITaxonomyProcessor::build_depth_lookup() {
    std::cout << "Building depth lookup table..." << std::endl;
    
    depth_lookup.clear();
    depth_lookup[1] = 0; // Root has depth 0
    
    // BFS to assign depths
    std::vector<uint32_t> current_level = {1};
    uint8_t current_depth = 0;
    
    while (!current_level.empty() && current_depth < 50) {
        std::vector<uint32_t> next_level;
        
        for (uint32_t taxon_id : current_level) {
            depth_lookup[taxon_id] = current_depth;
            
            // Find children
            for (const auto& [child_id, parent_id] : parent_lookup) {
                if (parent_id == taxon_id && depth_lookup.find(child_id) == depth_lookup.end()) {
                    next_level.push_back(child_id);
                }
            }
        }
        
        current_level = std::move(next_level);
        current_depth++;
    }
    
    std::cout << "Built depth lookup for " << depth_lookup.size() << " taxa" << std::endl;
}

// ===========================
// SimpleTaxonomyProcessor Implementation
// ===========================

SimpleTaxonomyProcessor::SimpleTaxonomyProcessor() : taxonomy_loaded(false) {
}

bool SimpleTaxonomyProcessor::load_ncbi_files(const std::string& nodes_dmp_path, const std::string& names_dmp_path) {
    std::cout << "Loading NCBI taxonomy with simple processor..." << std::endl;
    
    // Load nodes.dmp
    std::ifstream nodes_file(nodes_dmp_path);
    if (!nodes_file.is_open()) {
        std::cerr << "Cannot open nodes.dmp: " << nodes_dmp_path << std::endl;
        return false;
    }
    
    std::string line;
    while (std::getline(nodes_file, line)) {
        uint32_t taxon_id, parent_id;
        std::string rank;
        
        if (PhylogeneticUtils::parse_taxonomy_line(line, taxon_id, parent_id, taxon_names[taxon_id], rank)) {
            taxon_parents[taxon_id] = parent_id;
            taxon_ranks[taxon_id] = rank;
            taxon_names[taxon_id] = "taxon_" + std::to_string(taxon_id); // Default name
        }
    }
    nodes_file.close();
    
    // Load names.dmp if available
    std::ifstream names_file(names_dmp_path);
    if (names_file.is_open()) {
        while (std::getline(names_file, line)) {
            // Parse names.dmp format (simplified)
            std::istringstream iss(line);
            std::string taxon_str, name, unique_name, name_class;
            
            if (std::getline(iss, taxon_str, '|') &&
                std::getline(iss, name, '|') &&
                std::getline(iss, unique_name, '|') &&
                std::getline(iss, name_class, '|')) {
                
                // Trim whitespace
                name.erase(0, name.find_first_not_of(" \t"));
                name.erase(name.find_last_not_of(" \t") + 1);
                name_class.erase(0, name_class.find_first_not_of(" \t"));
                name_class.erase(name_class.find_last_not_of(" \t") + 1);
                
                if (name_class == "scientific name") {
                    try {
                        uint32_t taxon_id = std::stoul(taxon_str);
                        if (taxon_parents.find(taxon_id) != taxon_parents.end()) {
                            taxon_names[taxon_id] = name;
                        }
                    } catch (const std::exception&) {
                        continue;
                    }
                }
            }
        }
        names_file.close();
    }
    
    taxonomy_loaded = validate_taxonomy_consistency();
    std::cout << "Simple taxonomy loaded: " << taxon_parents.size() << " taxa" << std::endl;
    
    return taxonomy_loaded;
}

uint32_t SimpleTaxonomyProcessor::compute_simple_lca(uint32_t taxon1, uint32_t taxon2) {
    return LCAAlgorithms::compute_lca_pair(taxon1, taxon2, taxon_parents);
}

uint32_t SimpleTaxonomyProcessor::compute_lca_of_list(const std::vector<uint32_t>& taxon_list) {
    return LCAAlgorithms::compute_lca_multiple(taxon_list, taxon_parents);
}

std::string SimpleTaxonomyProcessor::get_name(uint32_t taxon_id) const {
    auto it = taxon_names.find(taxon_id);
    return (it != taxon_names.end()) ? it->second : ("taxon_" + std::to_string(taxon_id));
}

std::string SimpleTaxonomyProcessor::get_rank(uint32_t taxon_id) const {
    auto it = taxon_ranks.find(taxon_id);
    return (it != taxon_ranks.end()) ? it->second : "no rank";
}

uint32_t SimpleTaxonomyProcessor::get_parent(uint32_t taxon_id) const {
    auto it = taxon_parents.find(taxon_id);
    return (it != taxon_parents.end()) ? it->second : 0;
}

void SimpleTaxonomyProcessor::add_taxon(uint32_t taxon_id, uint32_t parent_id, const std::string& name, const std::string& rank) {
    taxon_parents[taxon_id] = parent_id;
    taxon_names[taxon_id] = name;
    taxon_ranks[taxon_id] = rank;
}

bool SimpleTaxonomyProcessor::save_taxonomy_tsv(const std::string& output_path) const {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        return false;
    }
    
    out << "taxon_id\tparent_id\tname\trank\n";
    for (const auto& [taxon_id, parent_id] : taxon_parents) {
        std::string name = get_name(taxon_id);
        std::string rank = get_rank(taxon_id);
        out << taxon_id << "\t" << parent_id << "\t" << name << "\t" << rank << "\n";
    }
    
    out.close();
    return true;
}

std::vector<uint32_t> SimpleTaxonomyProcessor::get_lineage(uint32_t taxon_id) const {
    return LCAAlgorithms::get_path_to_root(taxon_id, taxon_parents);
}

bool SimpleTaxonomyProcessor::validate_taxonomy_consistency() const {
    // Basic validation: ensure root exists and no cycles
    if (taxon_parents.find(1) == taxon_parents.end()) {
        std::cerr << "Warning: Root taxon (1) not found in taxonomy" << std::endl;
        return false;
    }
    
    // Check for basic consistency
    for (const auto& [taxon_id, parent_id] : taxon_parents) {
        if (taxon_id != 1 && parent_id != 0 && taxon_parents.find(parent_id) == taxon_parents.end()) {
            std::cerr << "Warning: Parent " << parent_id << " not found for taxon " << taxon_id << std::endl;
        }
    }
    
    return true;
}

// ===========================
// Utility Namespace Implementations
// ===========================

namespace PhylogeneticUtils {
    
    bool parse_taxonomy_line(const std::string& line, uint32_t& taxon_id, uint32_t& parent_id, 
                            std::string& name, std::string& rank) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '\t')) {
            fields.push_back(token);
        }
        
        if (fields.size() >= 4) {
            try {
                taxon_id = std::stoul(fields[0]);
                parent_id = std::stoul(fields[1]);
                name = fields[2];
                rank = fields[3];
                return true;
            } catch (const std::exception&) {
                return false;
            }
        }
        
        return false;
    }
    
    std::string format_taxonomy_line(uint32_t taxon_id, uint32_t parent_id, 
                                    const std::string& name, const std::string& rank) {
        return std::to_string(taxon_id) + "\t" + std::to_string(parent_id) + "\t" + name + "\t" + rank;
    }
    
    bool validate_species_list(const std::vector<uint32_t>& species_list) {
        if (species_list.empty()) return false;
        
        // Check for valid taxon IDs
        for (uint32_t taxon_id : species_list) {
            if (!is_valid_taxon_id(taxon_id)) {
                return false;
            }
        }
        
        return true;
    }
    
    bool is_valid_taxon_id(uint32_t taxon_id) {
        return taxon_id > 0 && taxon_id < 10000000; // Reasonable range
    }
}

namespace LCAAlgorithms {
    
    uint32_t compute_lca_pair(uint32_t taxon1, uint32_t taxon2, 
                             const std::unordered_map<uint32_t, uint32_t>& parents) {
        if (taxon1 == taxon2) return taxon1;
        if (taxon1 == 1 || taxon2 == 1) return 1;
        
        std::vector<uint32_t> path1 = get_path_to_root(taxon1, parents);
        std::vector<uint32_t> path2 = get_path_to_root(taxon2, parents);
        
        std::set<uint32_t> ancestors1(path1.begin(), path1.end());
        
        for (uint32_t ancestor : path2) {
            if (ancestors1.find(ancestor) != ancestors1.end()) {
                return ancestor;
            }
        }
        
        return 1; // Root fallback
    }
    
    uint32_t compute_lca_multiple(const std::vector<uint32_t>& taxa, 
                                 const std::unordered_map<uint32_t, uint32_t>& parents) {
        if (taxa.empty()) return 1;
        if (taxa.size() == 1) return taxa[0];
        
        uint32_t lca = taxa[0];
        for (size_t i = 1; i < taxa.size(); i++) {
            lca = compute_lca_pair(lca, taxa[i], parents);
            if (lca == 1) break; // Already at root
        }
        
        return lca;
    }
    
    std::vector<uint32_t> get_path_to_root(uint32_t taxon_id, 
                                          const std::unordered_map<uint32_t, uint32_t>& parents) {
        std::vector<uint32_t> path;
        uint32_t current = taxon_id;
        
        while (current != 1 && path.size() < 50) {
            path.push_back(current);
            auto parent_it = parents.find(current);
            if (parent_it == parents.end()) break;
            current = parent_it->second;
        }
        
        path.push_back(1); // Add root
        return path;
    }
}
