// arg_database.h
#ifndef ARG_DATABASE_H
#define ARG_DATABASE_H

#include <string>
#include <vector>
#include <map>
#include <cstdint>

// Structure for ARG metadata
struct ARGEntry {
    uint32_t arg_id;
    std::string gene_name;
    std::string gene_symbol;
    std::string drug_class;
    std::string resistance_mechanism;
    std::string accession;
    uint32_t protein_id;  // Links to your protein database
    uint16_t length;
    float identity_threshold;  // Minimum identity for positive hit
};

// ARG database that integrates with your existing protein database
class ARGDatabase {
private:
    std::vector<ARGEntry> arg_entries;
    std::map<uint32_t, uint32_t> protein_to_arg;  // protein_id -> arg_id
    std::map<std::string, std::vector<uint32_t>> drug_class_args;
    std::map<std::string, std::vector<uint32_t>> mechanism_args;
    
public:
    bool loadFromNCBI(const std::string& fasta_path, const std::string& metadata_path);
    bool loadFromCARD(const std::string& card_json_path);
    
    const ARGEntry* getARGByProteinId(uint32_t protein_id) const;
    std::vector<uint32_t> getARGsByDrugClass(const std::string& drug_class) const;
    void exportForGPU(const std::string& output_path) const;
    
    // Integration with existing protein database
    bool integrateWithProteinDB(const std::string& protein_db_path);
};

#endif
