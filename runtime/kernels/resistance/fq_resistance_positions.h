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
    int position;
    char wildtype_aa;
    char mutant_aa;
    std::string protein_id;
};

// Concrete implementation of FQ resistance database
class FQResistanceDatabase {
private:
    std::map<std::string, std::set<int>> qrdr_positions_by_gene;
    std::map<std::string, char> wildtype_aa_map;
    std::vector<ResistanceMutation> all_mutations;
    
public:
    FQResistanceDatabase() {}
    
    bool loadFromCSV(const std::string& csv_path) {
        // Minimal implementation - just return true for now
        // In production, this would load the CSV data
        return true;
    }
    
    // Get wildtype amino acid at a position
    char getWildtypeAA(const std::string& gene, uint16_t position) const {
        // Check if we have wildtype info for this position
        std::string key = gene + "_" + std::to_string(position);
        auto it = wildtype_aa_map.find(key);
        if (it != wildtype_aa_map.end()) {
            return it->second;
        }
        return 'X'; // Unknown
    }
    
    // Check if a position is in QRDR region
    bool isQRDRPosition(const std::string& gene, uint16_t position) const {
        auto it = qrdr_positions_by_gene.find(gene);
        if (it != qrdr_positions_by_gene.end()) {
            // Convert to 0-based for internal storage
            return it->second.find(position - 1) != it->second.end();
        }
        return false;
    }
};

#endif // FQ_RESISTANCE_POSITIONS_H