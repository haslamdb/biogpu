// compact_gpu_taxonomy.h
// Compact GPU taxonomy data structures and processing

#ifndef COMPACT_GPU_TAXONOMY_H
#define COMPACT_GPU_TAXONOMY_H

#include <cuda_runtime.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>

namespace BioGPU {
namespace CompactTaxonomy {

// Forward declarations
struct TaxonHashTable;
struct PhyloDistanceCache;

// Compact taxonomic node for GPU processing
struct CompactTaxonNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank_level;         // Encoded rank (0=root, 1=superkingdom, etc.)
    uint8_t depth_from_root;    // Distance from root node
    uint16_t reserved;          // Padding for alignment
};

// Hash table entry for GPU taxonomy lookup
struct TaxonHashEntry {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t depth;
    uint8_t rank_level;
    uint16_t name_offset;       // Offset into names array
};

// GPU-optimized hash table for taxonomy
struct TaxonHashTable {
    TaxonHashEntry* entries;    // Hash table entries
    uint32_t table_size;        // Size of hash table
    uint32_t num_entries;       // Number of valid entries
    uint32_t hash_mask;         // Mask for hash function
    char* names_data;           // Concatenated names string
    uint32_t names_data_size;   // Size of names data
};

// Phylogenetic distance cache for common LCA calculations
struct PhyloDistanceCache {
    uint32_t* taxon_pairs;      // [pair_index * 2 + {0,1}] = {taxon1, taxon2}
    uint8_t* distances;         // [pair_index] = distance
    uint32_t* lca_results;      // [pair_index] = LCA taxon
    uint32_t num_cached_pairs;  // Number of cached pairs
    uint32_t max_pairs;         // Maximum cache capacity
};

// Main compact taxonomy class
class CompactGPUTaxonomy {
private:
    // Host-side data
    std::unordered_map<uint32_t, CompactTaxonNode> taxonomy_nodes;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, std::string> taxon_ranks;
    std::unordered_map<uint32_t, uint32_t> parent_lookup;
    std::unordered_map<uint32_t, uint8_t> depth_lookup;
    std::unordered_map<uint32_t, uint8_t> rank_lookup;
    
    // GPU data
    TaxonHashTable* d_hash_table;
    PhyloDistanceCache* d_distance_cache;
    
    // Configuration
    bool use_distance_cache;
    uint32_t max_taxon_id;
    bool data_loaded;
    
public:
    CompactGPUTaxonomy(bool enable_cache = true) 
        : d_hash_table(nullptr), d_distance_cache(nullptr), 
          use_distance_cache(enable_cache), max_taxon_id(0), data_loaded(false) {
    }
    
    ~CompactGPUTaxonomy() {
        free_gpu_memory();
    }
    
    // Build compact taxonomy from NCBI files
    bool build_from_ncbi_files(const std::string& nodes_dmp_path, 
                               const std::string& names_dmp_path) {
        std::cout << "Building compact taxonomy from NCBI files..." << std::endl;
        
        if (!load_ncbi_nodes(nodes_dmp_path)) {
            std::cerr << "Failed to load nodes.dmp" << std::endl;
            return false;
        }
        
        if (!load_ncbi_names(names_dmp_path)) {
            std::cerr << "Failed to load names.dmp" << std::endl;
            return false;
        }
        
        if (!build_lookup_tables()) {
            std::cerr << "Failed to build lookup tables" << std::endl;
            return false;
        }
        
        if (!build_gpu_structures()) {
            std::cerr << "Failed to build GPU structures" << std::endl;
            return false;
        }
        
        data_loaded = true;
        std::cout << "Compact taxonomy built successfully" << std::endl;
        return true;
    }
    
    // Save compact taxonomy to binary file
    bool save_compact_taxonomy(const std::string& output_path) {
        if (!data_loaded) {
            std::cerr << "No taxonomy data to save" << std::endl;
            return false;
        }
        
        std::cout << "Saving compact taxonomy to: " << output_path << std::endl;
        
        std::ofstream file(output_path, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Cannot create output file: " << output_path << std::endl;
            return false;
        }
        
        // Write header
        uint32_t magic = 0x54415843; // "CTAX"
        uint32_t version = 1;
        file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
        file.write(reinterpret_cast<const char*>(&version), sizeof(version));
        file.write(reinterpret_cast<const char*>(&max_taxon_id), sizeof(max_taxon_id));
        
        uint32_t num_nodes = taxonomy_nodes.size();
        file.write(reinterpret_cast<const char*>(&num_nodes), sizeof(num_nodes));
        
        // Write taxonomy nodes
        for (const auto& [taxon_id, node] : taxonomy_nodes) {
            file.write(reinterpret_cast<const char*>(&node), sizeof(CompactTaxonNode));
        }
        
        // Write names
        uint32_t num_names = taxon_names.size();
        file.write(reinterpret_cast<const char*>(&num_names), sizeof(num_names));
        
        for (const auto& [taxon_id, name] : taxon_names) {
            file.write(reinterpret_cast<const char*>(&taxon_id), sizeof(taxon_id));
            uint32_t name_length = name.length();
            file.write(reinterpret_cast<const char*>(&name_length), sizeof(name_length));
            file.write(name.c_str(), name_length);
        }
        
        file.close();
        std::cout << "Compact taxonomy saved successfully" << std::endl;
        return true;
    }
    
    // Load compact taxonomy from binary file
    bool load_compact_taxonomy(const std::string& input_path) {
        std::cout << "Loading compact taxonomy from: " << input_path << std::endl;
        
        std::ifstream file(input_path, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Cannot open input file: " << input_path << std::endl;
            return false;
        }
        
        // Read header
        uint32_t magic, version;
        file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        file.read(reinterpret_cast<char*>(&version), sizeof(version));
        
        if (magic != 0x54415843) {
            std::cerr << "Invalid compact taxonomy file format" << std::endl;
            return false;
        }
        
        file.read(reinterpret_cast<char*>(&max_taxon_id), sizeof(max_taxon_id));
        
        // Read taxonomy nodes
        uint32_t num_nodes;
        file.read(reinterpret_cast<char*>(&num_nodes), sizeof(num_nodes));
        
        taxonomy_nodes.clear();
        for (uint32_t i = 0; i < num_nodes; i++) {
            CompactTaxonNode node;
            file.read(reinterpret_cast<char*>(&node), sizeof(CompactTaxonNode));
            taxonomy_nodes[node.taxon_id] = node;
        }
        
        // Read names
        uint32_t num_names;
        file.read(reinterpret_cast<char*>(&num_names), sizeof(num_names));
        
        taxon_names.clear();
        for (uint32_t i = 0; i < num_names; i++) {
            uint32_t taxon_id;
            file.read(reinterpret_cast<char*>(&taxon_id), sizeof(taxon_id));
            
            uint32_t name_length;
            file.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            
            std::string name(name_length, '\0');
            file.read(&name[0], name_length);
            
            taxon_names[taxon_id] = name;
        }
        
        file.close();
        
        // Rebuild lookup tables and GPU structures
        if (!build_lookup_tables()) {
            std::cerr << "Failed to rebuild lookup tables" << std::endl;
            return false;
        }
        
        if (!build_gpu_structures()) {
            std::cerr << "Failed to rebuild GPU structures" << std::endl;
            return false;
        }
        
        data_loaded = true;
        std::cout << "Compact taxonomy loaded successfully" << std::endl;
        return true;
    }
    
    // Accessor methods for database building integration
    std::unordered_map<uint32_t, uint32_t> get_parent_lookup_map() const {
        return parent_lookup;
    }
    
    std::unordered_map<uint32_t, std::string> get_name_lookup_map() const {
        return taxon_names;
    }
    
    std::unordered_map<uint32_t, std::string> get_rank_lookup_map() const {
        return taxon_ranks;
    }
    
    std::unordered_map<uint32_t, uint8_t> get_depth_lookup_map() const {
        return depth_lookup;
    }
    
    uint32_t get_max_taxon_id() const {
        return max_taxon_id;
    }
    
    // Get GPU structures for classification
    TaxonHashTable* get_gpu_hash_table() {
        return d_hash_table;
    }
    
    PhyloDistanceCache* get_gpu_distance_cache() {
        return d_distance_cache;
    }
    
    bool is_loaded() const {
        return data_loaded;
    }

private:
    bool load_ncbi_nodes(const std::string& nodes_path) {
        std::ifstream file(nodes_path);
        if (!file.is_open()) {
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            // Parse nodes.dmp format: taxid | parent_taxid | rank | ...
            std::vector<std::string> fields = split_dmp_line(line);
            
            if (fields.size() >= 3) {
                try {
                    uint32_t taxon_id = std::stoul(fields[0]);
                    uint32_t parent_id = std::stoul(fields[1]);
                    std::string rank = fields[2];
                    
                    CompactTaxonNode node;
                    node.taxon_id = taxon_id;
                    node.parent_id = parent_id;
                    node.rank_level = encode_rank(rank);
                    node.depth_from_root = 0; // Will be calculated later
                    
                    taxonomy_nodes[taxon_id] = node;
                    taxon_ranks[taxon_id] = rank;
                    max_taxon_id = std::max(max_taxon_id, taxon_id);
                    
                } catch (const std::exception& e) {
                    continue;
                }
            }
        }
        
        file.close();
        return !taxonomy_nodes.empty();
    }
    
    bool load_ncbi_names(const std::string& names_path) {
        std::ifstream file(names_path);
        if (!file.is_open()) {
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            // Parse names.dmp format: taxid | name | unique_name | name_class |
            std::vector<std::string> fields = split_dmp_line(line);
            
            if (fields.size() >= 4) {
                try {
                    uint32_t taxon_id = std::stoul(fields[0]);
                    std::string name = fields[1];
                    std::string name_class = fields[3];
                    
                    // Only store scientific names
                    if (name_class == "scientific name") {
                        taxon_names[taxon_id] = name;
                    }
                    
                } catch (const std::exception& e) {
                    continue;
                }
            }
        }
        
        file.close();
        return !taxon_names.empty();
    }
    
    bool build_lookup_tables() {
        // Build parent lookup
        parent_lookup.clear();
        depth_lookup.clear();
        rank_lookup.clear();
        
        for (const auto& [taxon_id, node] : taxonomy_nodes) {
            parent_lookup[taxon_id] = node.parent_id;
            rank_lookup[taxon_id] = node.rank_level;
        }
        
        // Calculate depths
        calculate_depths();
        
        return true;
    }
    
    void calculate_depths() {
        // Set root depth
        if (taxonomy_nodes.find(1) != taxonomy_nodes.end()) {
            depth_lookup[1] = 0;
            taxonomy_nodes[1].depth_from_root = 0;
        }
        
        // Iteratively calculate depths
        bool changed = true;
        int iterations = 0;
        
        while (changed && iterations < 50) {
            changed = false;
            iterations++;
            
            for (auto& [taxon_id, node] : taxonomy_nodes) {
                if (depth_lookup.find(taxon_id) == depth_lookup.end()) {
                    uint32_t parent_id = node.parent_id;
                    
                    if (depth_lookup.find(parent_id) != depth_lookup.end()) {
                        uint8_t depth = depth_lookup[parent_id] + 1;
                        depth_lookup[taxon_id] = depth;
                        node.depth_from_root = depth;
                        changed = true;
                    }
                }
            }
        }
    }
    
    bool build_gpu_structures() {
        // Allocate GPU hash table
        cudaError_t error = cudaMalloc(&d_hash_table, sizeof(TaxonHashTable));
        if (error != cudaSuccess) {
            return false;
        }
        
        // Build and copy hash table data
        // (Implementation details would go here)
        
        // Allocate distance cache if enabled
        if (use_distance_cache) {
            error = cudaMalloc(&d_distance_cache, sizeof(PhyloDistanceCache));
            if (error != cudaSuccess) {
                return false;
            }
            
            // Initialize cache
            // (Implementation details would go here)
        }
        
        return true;
    }
    
    void free_gpu_memory() {
        if (d_hash_table) {
            cudaFree(d_hash_table);
            d_hash_table = nullptr;
        }
        
        if (d_distance_cache) {
            cudaFree(d_distance_cache);
            d_distance_cache = nullptr;
        }
    }
    
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
    
    uint8_t encode_rank(const std::string& rank) {
        if (rank == "root") return 0;
        else if (rank == "superkingdom") return 1;
        else if (rank == "kingdom") return 2;
        else if (rank == "phylum") return 3;
        else if (rank == "class") return 4;
        else if (rank == "order") return 5;
        else if (rank == "family") return 6;
        else if (rank == "genus") return 7;
        else if (rank == "species") return 8;
        else if (rank == "subspecies") return 9;
        else return 255; // "no rank"
    }
};

} // namespace CompactTaxonomy
} // namespace BioGPU

#endif // COMPACT_GPU_TAXONOMY_H