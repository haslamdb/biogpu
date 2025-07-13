// fix_pipeline_resistance_loader.cpp
// Modified resistance database loader that doesn't expect drug/resistance info

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cstdint>

// Simplified resistance structure without drug info
struct CleanResistanceMutation {
    uint32_t gene_id;
    uint16_t position;
    char wildtype_aa;
    std::vector<char> mutant_aas;
};

class CleanResistanceLoader {
private:
    std::map<uint32_t, std::string> gene_id_to_name;
    std::map<uint32_t, std::string> species_id_to_name;
    std::vector<CleanResistanceMutation> mutations;
    
public:
    bool loadResistanceDatabase(const std::string& db_path) {
        std::string resistance_file = db_path + "/resistance_db.json";
        std::ifstream file(resistance_file);
        
        if (!file.is_open()) {
            std::cerr << "WARNING: Cannot open resistance database: " << resistance_file << std::endl;
            std::cerr << "Pipeline will continue without resistance mutation information" << std::endl;
            return false;
        }
        
        // Simple JSON parsing (would use proper JSON library in production)
        std::string line;
        bool in_gene_map = false;
        bool in_species_map = false;
        bool in_mutations = false;
        
        while (std::getline(file, line)) {
            // Detect sections
            if (line.find("\"gene_map\"") != std::string::npos) {
                in_gene_map = true;
                in_species_map = false;
                in_mutations = false;
                continue;
            }
            if (line.find("\"species_map\"") != std::string::npos) {
                in_gene_map = false;
                in_species_map = true;
                in_mutations = false;
                continue;
            }
            if (line.find("\"mutations\"") != std::string::npos) {
                in_gene_map = false;
                in_species_map = false;
                in_mutations = true;
                continue;
            }
            
            // Parse gene map entries
            if (in_gene_map && line.find("\":") != std::string::npos) {
                size_t quote1 = line.find("\"");
                size_t quote2 = line.find("\"", quote1 + 1);
                size_t quote3 = line.find("\"", quote2 + 1);
                size_t quote4 = line.find("\"", quote3 + 1);
                
                if (quote1 != std::string::npos && quote4 != std::string::npos) {
                    std::string id_str = line.substr(quote1 + 1, quote2 - quote1 - 1);
                    std::string name = line.substr(quote3 + 1, quote4 - quote3 - 1);
                    
                    try {
                        uint32_t id = std::stoul(id_str);
                        gene_id_to_name[id] = name;
                    } catch (...) {
                        // Skip malformed entries
                    }
                }
            }
            
            // Parse species map entries similarly
            if (in_species_map && line.find("\":") != std::string::npos) {
                size_t quote1 = line.find("\"");
                size_t quote2 = line.find("\"", quote1 + 1);
                size_t quote3 = line.find("\"", quote2 + 1);
                size_t quote4 = line.find("\"", quote3 + 1);
                
                if (quote1 != std::string::npos && quote4 != std::string::npos) {
                    std::string id_str = line.substr(quote1 + 1, quote2 - quote1 - 1);
                    std::string name = line.substr(quote3 + 1, quote4 - quote3 - 1);
                    
                    try {
                        uint32_t id = std::stoul(id_str);
                        species_id_to_name[id] = name;
                    } catch (...) {
                        // Skip malformed entries
                    }
                }
            }
            
            // Parse mutations
            if (in_mutations) {
                CleanResistanceMutation mut;
                
                // Look for gene_id
                size_t gene_id_pos = line.find("\"gene_id\":");
                if (gene_id_pos != std::string::npos) {
                    size_t start = line.find(":", gene_id_pos) + 1;
                    size_t end = line.find(",", start);
                    if (end == std::string::npos) end = line.find("}", start);
                    
                    std::string id_str = line.substr(start, end - start);
                    try {
                        mut.gene_id = std::stoul(id_str);
                    } catch (...) {
                        continue;
                    }
                }
                
                // Look for position
                size_t pos_pos = line.find("\"position\":");
                if (pos_pos != std::string::npos) {
                    size_t start = line.find(":", pos_pos) + 1;
                    size_t end = line.find(",", start);
                    if (end == std::string::npos) end = line.find("}", start);
                    
                    std::string pos_str = line.substr(start, end - start);
                    try {
                        mut.position = std::stoul(pos_str);
                    } catch (...) {
                        continue;
                    }
                }
                
                // Look for wildtype_aa
                size_t wt_pos = line.find("\"wildtype_aa\":");
                if (wt_pos != std::string::npos) {
                    size_t quote_start = line.find("\"", wt_pos + 14);
                    size_t quote_end = line.find("\"", quote_start + 1);
                    if (quote_start != std::string::npos && quote_end != std::string::npos) {
                        std::string wt = line.substr(quote_start + 1, quote_end - quote_start - 1);
                        if (!wt.empty()) {
                            mut.wildtype_aa = wt[0];
                        }
                    }
                }
                
                // For now, skip parsing mutant_aas array (complex)
                // In production, use proper JSON parser
                
                // If we have minimum required fields, add mutation
                if (mut.position > 0 && mut.wildtype_aa != 0) {
                    mutations.push_back(mut);
                }
            }
        }
        
        file.close();
        
        std::cout << "Loaded resistance database with " << gene_id_to_name.size() 
                  << " genes and " << mutations.size() << " mutation positions" << std::endl;
        
        return true;
    }
    
    std::string getGeneName(uint32_t gene_id) const {
        auto it = gene_id_to_name.find(gene_id);
        if (it != gene_id_to_name.end()) {
            return it->second;
        }
        return "Gene" + std::to_string(gene_id);
    }
    
    std::string getSpeciesName(uint32_t species_id) const {
        auto it = species_id_to_name.find(species_id);
        if (it != species_id_to_name.end()) {
            return it->second;
        }
        return "Species" + std::to_string(species_id);
    }
    
    const std::vector<CleanResistanceMutation>& getMutations() const {
        return mutations;
    }
    
    size_t getGeneCount() const {
        return gene_id_to_name.size();
    }
    
    size_t getSpeciesCount() const {
        return species_id_to_name.size();
    }
};

// Example usage in pipeline
void integrateWithPipeline(const std::string& resistance_db_path) {
    CleanResistanceLoader loader;
    
    if (loader.loadResistanceDatabase(resistance_db_path)) {
        std::cout << "Successfully loaded clean resistance database" << std::endl;
        std::cout << "Genes: " << loader.getGeneCount() << std::endl;
        std::cout << "Species: " << loader.getSpeciesCount() << std::endl;
        std::cout << "Known mutation positions: " << loader.getMutations().size() << std::endl;
        
        // Use the loader in your pipeline
        // The mutations only have position info, no drug/resistance data
    } else {
        std::cout << "Running pipeline without resistance mutation info" << std::endl;
    }
}