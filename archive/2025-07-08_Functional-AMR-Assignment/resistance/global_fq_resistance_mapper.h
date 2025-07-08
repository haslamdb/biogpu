// global_fq_resistance_mapper.h
// Unified global mapping system for fluoroquinolone resistance mutations
// This header provides a consistent interface for all pipeline components

#ifndef GLOBAL_FQ_RESISTANCE_MAPPER_H
#define GLOBAL_FQ_RESISTANCE_MAPPER_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cstdint>
#include <memory>

// Core resistance mutation structure used throughout the pipeline
struct FQResistanceMutation {
    std::string species;           // Normalized species name (underscores)
    std::string gene;              // Gene name (gyrA, gyrB, parC, parE, etc.)
    uint16_t position;             // 1-based amino acid position
    char wildtype_aa;              // Expected wildtype amino acid
    char mutant_aa;                // Resistant amino acid variant
    
    // Unique key for this mutation
    std::string getKey() const {
        return species + "_" + gene + "_" + std::to_string(position) + "_" + wildtype_aa + std::to_string(mutant_aa);
    }
    
    // Key for position only (to check all possible mutations at a position)
    std::string getPositionKey() const {
        return species + "_" + gene + "_" + std::to_string(position);
    }
};

// CUDA-compatible structure for GPU kernels
struct FQResistanceMutationGPU {
    uint32_t species_gene_id;      // Combined species-gene ID from your database
    uint16_t position;             // 1-based amino acid position
    char wildtype_aa;              // Expected wildtype amino acid
    char mutant_aas[8];            // Up to 8 resistant variants
    uint8_t num_mutants;           // Number of resistant variants
    float resistance_score;        // Clinical significance (reserved for future use)
};

// Main resistance mapper class
class GlobalFQResistanceMapper {
private:
    // Master mutation database loaded from CSV
    std::vector<FQResistanceMutation> all_mutations;
    
    // Lookup tables for fast access
    std::map<std::string, std::vector<FQResistanceMutation>> mutations_by_position;  // Key: species_gene_position
    std::map<std::string, std::set<char>> resistant_aas_by_position;                // Resistant AAs per position
    std::map<std::string, char> wildtype_by_position;                               // Wildtype AA per position
    
    // Species and gene name normalization maps
    std::map<std::string, std::string> species_normalization;
    std::map<std::string, std::string> gene_normalization;
    
    // ID mappings for database integration
    std::map<std::string, uint32_t> species_to_id;
    std::map<uint32_t, std::string> id_to_species;
    std::map<std::string, uint32_t> gene_to_id;
    std::map<uint32_t, std::string> id_to_gene;
    std::map<std::string, uint32_t> species_gene_to_id;  // Combined IDs
    std::map<uint32_t, std::pair<std::string, std::string>> id_to_species_gene;
    
    // GPU data (allocated on demand)
    FQResistanceMutationGPU* gpu_mutations;
    uint32_t num_gpu_mutations;
    bool gpu_data_ready;
    
    // Helper methods
    std::string normalizeSpeciesName(const std::string& species);
    std::string normalizeGeneName(const std::string& gene);
    void buildLookupTables();
    void initializeNormalizationMaps();
    
public:
    GlobalFQResistanceMapper();
    ~GlobalFQResistanceMapper();
    
    // Load resistance data from CSV file
    bool loadFromCSV(const std::string& csv_path);
    
    // Load hardcoded resistance data (used when CSV not provided)
    bool loadHardcodedData();
    
    // Load ID mappings from protein database metadata
    bool loadDatabaseMappings(const std::string& protein_db_path);
    
    // Check if a mutation is a known resistance mutation
    bool isResistanceMutation(const std::string& species, const std::string& gene, 
                             uint16_t position, char wildtype_aa, char observed_aa) const;
    
    // Get all resistant amino acids for a position
    std::vector<char> getResistantAAs(const std::string& species, const std::string& gene, 
                                     uint16_t position) const;
    
    // Get wildtype amino acid for a position
    char getWildtypeAA(const std::string& species, const std::string& gene, 
                       uint16_t position) const;
    
    // Check if position is a known resistance position
    bool isResistancePosition(const std::string& species, const std::string& gene, 
                             uint16_t position) const;
    
    // Get all mutations for a species-gene combination
    std::vector<FQResistanceMutation> getMutationsForGene(const std::string& species, 
                                                         const std::string& gene) const;
    
    // ID-based lookups for integration with existing pipeline
    bool isResistanceMutationByID(uint32_t species_gene_id, uint16_t position, 
                                  char wildtype_aa, char observed_aa) const;
    
    std::vector<char> getResistantAAsByID(uint32_t species_gene_id, uint16_t position) const;
    
    char getWildtypeAAByID(uint32_t species_gene_id, uint16_t position) const;
    
    // Convert between names and IDs
    uint32_t getSpeciesGeneID(const std::string& species, const std::string& gene) const;
    std::pair<std::string, std::string> getSpeciesGeneFromID(uint32_t species_gene_id) const;
    
    // GPU support
    bool prepareGPUData();  // Prepare data for GPU transfer
    FQResistanceMutationGPU* getGPUMutations() { return gpu_mutations; }
    uint32_t getNumGPUMutations() const { return num_gpu_mutations; }
    
    // Statistics and debugging
    void printSummary() const;
    void printMutationsForGene(const std::string& species, const std::string& gene) const;
    
    // Singleton pattern for global access
    static GlobalFQResistanceMapper& getInstance() {
        static GlobalFQResistanceMapper instance;
        return instance;
    }
};

// C interface for integration with existing code
extern "C" {
    // Initialize the global mapper
    int init_global_fq_mapper(const char* csv_path, const char* protein_db_path);
    
    // Check if a mutation is resistant
    bool is_fq_resistance_mutation(const char* species, const char* gene, 
                                  uint16_t position, char wildtype, char observed);
    
    // Get wildtype amino acid
    char get_fq_wildtype_aa(const char* species, const char* gene, uint16_t position);
    
    // ID-based queries
    bool is_fq_resistance_mutation_by_id(uint32_t species_gene_id, uint16_t position, 
                                        char wildtype, char observed);
    
    // Get GPU data pointer
    void* get_fq_gpu_mutations();
    uint32_t get_fq_gpu_mutation_count();
    
    // Cleanup
    void cleanup_global_fq_mapper();
}

#endif // GLOBAL_FQ_RESISTANCE_MAPPER_H