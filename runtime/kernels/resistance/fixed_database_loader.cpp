// fixed_database_loader.cpp
// Properly load and track species-gene mappings in the resistance detection pipeline

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <memory>
#include <json/json.h>  // You may need to install jsoncpp

// Structure to represent a species-gene protein
struct SpeciesGeneProtein {
    std::string species;
    std::string gene;
    int species_gene_id;
    std::string sequence;
    int length;
};

// Enhanced structure for protein matches that preserves species-gene info
struct EnhancedProteinMatch {
    uint32_t read_id;
    int8_t frame;
    
    // Species-gene specific identifiers
    uint32_t species_gene_id;  // Unique ID for species-gene combination
    std::string species_name;
    std::string gene_name;
    
    // Alignment details
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    
    // Mutation details
    uint8_t num_mutations;
    uint16_t mutation_positions[10];
    char ref_aas[10];
    char query_aas[10];
    
    std::string aligned_peptide;
};

class SpeciesGeneDatabaseLoader {
private:
    // Core data structures
    std::map<int, SpeciesGeneProtein> proteins_by_id;
    std::map<std::string, int> species_gene_to_id;
    std::map<int, std::pair<std::string, std::string>> id_to_species_gene;
    
    // For GPU transfer
    std::vector<std::string> concatenated_sequences;
    std::vector<int> sequence_lengths;
    std::vector<int> sequence_offsets;
    std::vector<int> species_gene_ids;
    
    // Statistics
    int total_sequences = 0;
    std::set<std::string> unique_species;
    std::set<std::string> unique_genes;
    
public:
    bool loadFromJSON(const std::string& json_path) {
        std::cout << "Loading species-specific protein database from: " << json_path << std::endl;
        
        std::ifstream file(json_path);
        if (!file.is_open()) {
            std::cerr << "ERROR: Cannot open database file: " << json_path << std::endl;
            return false;
        }
        
        Json::Value root;
        Json::Reader reader;
        
        if (!reader.parse(file, root)) {
            std::cerr << "ERROR: Failed to parse JSON: " << reader.getFormattedErrorMessages() << std::endl;
            return false;
        }
        
        // Load metadata
        const Json::Value& metadata = root["metadata"];
        std::cout << "Database version: " << metadata["version"].asString() << std::endl;
        std::cout << "Number of sequences: " << metadata["num_sequences"].asInt() << std::endl;
        std::cout << "Number of species: " << metadata["num_species"].asInt() << std::endl;
        std::cout << "Number of genes: " << metadata["num_genes"].asInt() << std::endl;
        
        // Load species-gene ID mapping
        const Json::Value& id_map = root["species_gene_id_map"];
        for (const auto& key : id_map.getMemberNames()) {
            int id = id_map[key].asInt();
            species_gene_to_id[key] = id;
            
            // Parse key to get species and gene
            size_t last_underscore = key.rfind('_');
            if (last_underscore != std::string::npos) {
                std::string species = key.substr(0, last_underscore);
                std::string gene = key.substr(last_underscore + 1);
                id_to_species_gene[id] = {species, gene};
            }
        }
        
        // Load sequences
        const Json::Value& sequences = root["sequences"];
        for (const auto& key : sequences.getMemberNames()) {
            const Json::Value& seq_data = sequences[key];
            
            SpeciesGeneProtein protein;
            protein.species = seq_data["species"].asString();
            protein.gene = seq_data["gene"].asString();
            protein.species_gene_id = seq_data["species_gene_id"].asInt();
            protein.sequence = seq_data["sequence"].asString();
            protein.length = seq_data["length"].asInt();
            
            proteins_by_id[protein.species_gene_id] = protein;
            
            unique_species.insert(protein.species);
            unique_genes.insert(protein.gene);
            total_sequences++;
        }
        
        std::cout << "Successfully loaded " << total_sequences << " protein sequences" << std::endl;
        std::cout << "Unique species: " << unique_species.size() << std::endl;
        std::cout << "Unique genes: " << unique_genes.size() << std::endl;
        
        // Prepare data for GPU
        prepareGPUData();
        
        return true;
    }
    
    bool loadFromTSV(const std::string& sequences_file, const std::string& mapping_file) {
        // Alternative loading from TSV files
        std::cout << "Loading from TSV files..." << std::endl;
        
        // Load ID mappings first
        std::ifstream map_file(mapping_file);
        if (!map_file.is_open()) {
            std::cerr << "ERROR: Cannot open mapping file: " << mapping_file << std::endl;
            return false;
        }
        
        std::string line;
        std::getline(map_file, line); // Skip header
        
        while (std::getline(map_file, line)) {
            std::istringstream iss(line);
            std::string species, gene;
            int id;
            
            if (iss >> species >> gene >> id) {
                std::string key = species + "_" + gene;
                species_gene_to_id[key] = id;
                id_to_species_gene[id] = {species, gene};
            }
        }
        
        std::cout << "Loaded " << species_gene_to_id.size() << " species-gene mappings" << std::endl;
        
        // Load sequences
        std::ifstream seq_file(sequences_file);
        if (!seq_file.is_open()) {
            std::cerr << "ERROR: Cannot open sequences file: " << sequences_file << std::endl;
            return false;
        }
        
        // Parse FASTA format
        std::string header, sequence;
        while (std::getline(seq_file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Process previous sequence if any
                if (!header.empty() && !sequence.empty()) {
                    processSequence(header, sequence);
                }
                
                header = line.substr(1);
                sequence.clear();
            } else {
                sequence += line;
            }
        }
        
        // Process last sequence
        if (!header.empty() && !sequence.empty()) {
            processSequence(header, sequence);
        }
        
        std::cout << "Loaded " << total_sequences << " sequences" << std::endl;
        
        prepareGPUData();
        return true;
    }
    
    // Get protein info by species-gene ID
    const SpeciesGeneProtein* getProtein(int species_gene_id) const {
        auto it = proteins_by_id.find(species_gene_id);
        if (it != proteins_by_id.end()) {
            return &it->second;
        }
        return nullptr;
    }
    
    // Get ID for species-gene combination
    int getSpeciesGeneId(const std::string& species, const std::string& gene) const {
        std::string key = species + "_" + gene;
        auto it = species_gene_to_id.find(key);
        if (it != species_gene_to_id.end()) {
            return it->second;
        }
        return -1;
    }
    
    // Get species and gene from ID
    std::pair<std::string, std::string> getSpeciesGene(int species_gene_id) const {
        auto it = id_to_species_gene.find(species_gene_id);
        if (it != id_to_species_gene.end()) {
            return it->second;
        }
        return {"", ""};
    }
    
    // Get GPU-ready data
    void getGPUData(const char** sequences, int* lengths, int* offsets, int* ids, int& num_sequences) {
        num_sequences = concatenated_sequences.size();
        
        for (int i = 0; i < num_sequences; i++) {
            sequences[i] = concatenated_sequences[i].c_str();
            lengths[i] = sequence_lengths[i];
            offsets[i] = sequence_offsets[i];
            ids[i] = species_gene_ids[i];
        }
    }
    
    // Convert protein match results preserving species-gene info
    std::vector<EnhancedProteinMatch> convertMatches(
        const std::vector<void*>& raw_matches,
        const std::vector<uint32_t>& match_counts
    ) {
        std::vector<EnhancedProteinMatch> results;
        
        for (size_t read_idx = 0; read_idx < match_counts.size(); read_idx++) {
            uint32_t count = match_counts[read_idx];
            
            for (uint32_t i = 0; i < count; i++) {
                // Convert raw match to enhanced match
                EnhancedProteinMatch match;
                
                // Extract basic match info from raw data
                // (Implementation depends on your raw match structure)
                
                // Lookup species and gene info
                auto species_gene = getSpeciesGene(match.species_gene_id);
                match.species_name = species_gene.first;
                match.gene_name = species_gene.second;
                
                results.push_back(match);
            }
        }
        
        return results;
    }
    
    void printSummary() const {
        std::cout << "\n=== Database Summary ===" << std::endl;
        std::cout << "Total sequences: " << total_sequences << std::endl;
        std::cout << "Species: ";
        for (const auto& species : unique_species) {
            std::cout << species << " ";
        }
        std::cout << std::endl;
        
        std::cout << "Genes: ";
        for (const auto& gene : unique_genes) {
            std::cout << gene << " ";
        }
        std::cout << std::endl;
        
        // Print sequence counts per gene
        std::map<std::string, int> gene_counts;
        for (const auto& pair : proteins_by_id) {
            gene_counts[pair.second.gene]++;
        }
        
        std::cout << "\nSequences per gene:" << std::endl;
        for (const auto& pair : gene_counts) {
            std::cout << "  " << pair.first << ": " << pair.second << std::endl;
        }
    }
    
private:
    void processSequence(const std::string& header, const std::string& sequence) {
        // Parse header to extract species and gene
        size_t last_underscore = header.rfind('_');
        if (last_underscore == std::string::npos) {
            std::cerr << "WARNING: Cannot parse header: " << header << std::endl;
            return;
        }
        
        std::string species = header.substr(0, last_underscore);
        std::string gene = header.substr(last_underscore + 1);
        
        // Get or create ID
        std::string key = species + "_" + gene;
        auto it = species_gene_to_id.find(key);
        if (it == species_gene_to_id.end()) {
            std::cerr << "WARNING: No ID mapping for " << key << std::endl;
            return;
        }
        
        int species_gene_id = it->second;
        
        // Create protein entry
        SpeciesGeneProtein protein;
        protein.species = species;
        protein.gene = gene;
        protein.species_gene_id = species_gene_id;
        protein.sequence = sequence;
        protein.length = sequence.length();
        
        proteins_by_id[species_gene_id] = protein;
        
        unique_species.insert(species);
        unique_genes.insert(gene);
        total_sequences++;
    }
    
    void prepareGPUData() {
        // Prepare concatenated sequences for GPU
        concatenated_sequences.clear();
        sequence_lengths.clear();
        sequence_offsets.clear();
        species_gene_ids.clear();
        
        int current_offset = 0;
        
        for (const auto& pair : proteins_by_id) {
            const SpeciesGeneProtein& protein = pair.second;
            
            concatenated_sequences.push_back(protein.sequence);
            sequence_lengths.push_back(protein.length);
            sequence_offsets.push_back(current_offset);
            species_gene_ids.push_back(protein.species_gene_id);
            
            current_offset += protein.length;
        }
    }
};

// Example usage in main pipeline
class UpdatedResistancePipeline {
private:
    SpeciesGeneDatabaseLoader db_loader;
    
public:
    bool initialize(const std::string& database_path) {
        // Try JSON first, then fall back to TSV
        if (database_path.find(".json") != std::string::npos) {
            return db_loader.loadFromJSON(database_path);
        } else {
            // Assume TSV format
            std::string seq_file = database_path + "/wildtype_proteins.fasta";
            std::string map_file = database_path + "/species_gene_to_id.tsv";
            return db_loader.loadFromTSV(seq_file, map_file);
        }
    }
    
    void processResults(const std::vector<EnhancedProteinMatch>& matches) {
        std::cout << "\n=== Resistance Detection Results ===" << std::endl;
        
        // Group by species and gene
        std::map<std::string, std::vector<const EnhancedProteinMatch*>> results_by_species_gene;
        
        for (const auto& match : matches) {
            std::string key = match.species_name + " " + match.gene_name;
            results_by_species_gene[key].push_back(&match);
        }
        
        // Report findings
        for (const auto& pair : results_by_species_gene) {
            std::cout << "\n" << pair.first << ":" << std::endl;
            
            for (const auto* match : pair.second) {
                std::cout << "  Read " << match->read_id << ": ";
                std::cout << "Score=" << match->alignment_score;
                std::cout << ", Identity=" << (match->identity * 100) << "%";
                
                if (match->num_mutations > 0) {
                    std::cout << ", Mutations: ";
                    for (int i = 0; i < match->num_mutations; i++) {
                        std::cout << match->ref_aas[i] << match->mutation_positions[i] 
                                 << match->query_aas[i];
                        if (i < match->num_mutations - 1) std::cout << ", ";
                    }
                }
                
                std::cout << std::endl;
            }
        }
    }
};

// Integration test
#ifdef STANDALONE_LOADER
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <database_path>" << std::endl;
        return 1;
    }
    
    UpdatedResistancePipeline pipeline;
    
    if (!pipeline.initialize(argv[1])) {
        std::cerr << "Failed to initialize pipeline" << std::endl;
        return 1;
    }
    
    // The pipeline would then process reads and use the species-gene specific database
    
    return 0;
}
#endif