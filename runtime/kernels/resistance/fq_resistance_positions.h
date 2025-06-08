// fq_resistance_positions.h
// Dynamic FQ resistance position loader for QRDR mutations

#ifndef FQ_RESISTANCE_POSITIONS_H
#define FQ_RESISTANCE_POSITIONS_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

struct ResistanceMutation {
    std::string species;
    std::string gene;
    int position;  // 1-based position from CSV
    char wildtype_aa;
    char mutant_aa;
    std::string protein_id;
};

class FQResistanceDatabase {
private:
    // Map: gene -> set of QRDR positions (0-based)
    std::map<std::string, std::set<int>> qrdr_positions_by_gene;
    
    // Map: gene_position -> expected wildtype AA
    std::map<std::string, char> wildtype_aa_map;
    
    // All loaded mutations
    std::vector<ResistanceMutation> all_mutations;
    
public:
    FQResistanceDatabase() {}
    
    bool loadFromCSV(const std::string& csv_path) {
        std::ifstream file(csv_path);
        if (!file.good()) {
            std::cerr << "Error: Cannot open resistance mutation CSV: " << csv_path << std::endl;
            return false;
        }
        
        std::string line;
        bool first_line = true;
        
        while (std::getline(file, line)) {
            if (first_line) {
                first_line = false;
                continue; // Skip header
            }
            
            // Skip empty lines
            if (line.empty()) continue;
            
            // Parse CSV line (handle quoted fields)
            std::vector<std::string> fields;
            std::stringstream ss(line);
            std::string field;
            bool in_quotes = false;
            std::string current_field;
            
            for (char c : line) {
                if (c == '"') {
                    in_quotes = !in_quotes;
                } else if (c == ',' && !in_quotes) {
                    fields.push_back(current_field);
                    current_field.clear();
                } else {
                    current_field += c;
                }
            }
            fields.push_back(current_field); // Last field
            
            if (fields.size() >= 6) {
                ResistanceMutation mut;
                mut.species = fields[1];
                mut.gene = fields[2];
                
                // Parse position (handle potential errors)
                try {
                    mut.position = std::stoi(fields[3]);
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing position in line: " << line << std::endl;
                    continue;
                }
                
                // Make sure wildtype and mutant fields are not empty
                if (fields[4].empty() || fields[5].empty()) {
                    continue;
                }
                
                mut.wildtype_aa = fields[4][0];
                mut.mutant_aa = fields[5][0];
                mut.protein_id = fields.size() > 6 ? fields[6] : "";
                
                all_mutations.push_back(mut);
                
                // Store QRDR position (convert to 0-based)
                int pos_0based = mut.position - 1;
                qrdr_positions_by_gene[mut.gene].insert(pos_0based);
                
                // Store expected wildtype AA
                std::string key = mut.gene + "_" + std::to_string(pos_0based);
                wildtype_aa_map[key] = mut.wildtype_aa;
            }
        }
        
        file.close();
        
        std::cout << "Loaded " << all_mutations.size() << " FQ resistance mutations" << std::endl;
        std::cout << "QRDR positions by gene:" << std::endl;
        for (const auto& pair : qrdr_positions_by_gene) {
            std::cout << "  " << pair.first << ": ";
            for (int pos : pair.second) {
                std::cout << (pos + 1) << " ";  // Display as 1-based
            }
            std::cout << std::endl;
        }
        
        return true;
    }
    
    // Check if a position is a QRDR position for a given gene
    bool isQRDRPosition(const std::string& gene, int position_0based) const {
        auto it = qrdr_positions_by_gene.find(gene);
        if (it != qrdr_positions_by_gene.end()) {
            return it->second.find(position_0based) != it->second.end();
        }
        return false;
    }
    
    // Get expected wildtype AA at a position
    char getWildtypeAA(const std::string& gene, int position_0based) const {
        std::string key = gene + "_" + std::to_string(position_0based);
        auto it = wildtype_aa_map.find(key);
        if (it != wildtype_aa_map.end()) {
            return it->second;
        }
        return '?';
    }
    
    // Get all QRDR positions for a gene
    std::set<int> getQRDRPositions(const std::string& gene) const {
        auto it = qrdr_positions_by_gene.find(gene);
        if (it != qrdr_positions_by_gene.end()) {
            return it->second;
        }
        return std::set<int>();
    }
    
    // Map gene_id to gene name (based on the integrated database structure)
    std::string getGeneName(uint32_t gene_id) const {
        // This mapping should match the integrated database
        // Based on the protein_database.json, we have:
        // gene_list: ["grlA", "gyrA", "parC", "parE"]
        switch (gene_id) {
            case 0: return "grlA";
            case 1: return "gyrA";
            case 2: return "parC";
            case 3: return "parE";
            default: return "unknown";
        }
    }
};

// Global instance for use in both runtime and reporting
extern FQResistanceDatabase* g_fq_resistance_db;

#endif // FQ_RESISTANCE_POSITIONS_H