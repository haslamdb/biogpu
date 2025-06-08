// fq_mutation_reporter.cpp
// Enhanced FQ resistance mutation reporter with detailed polymorphism analysis
// Generates comprehensive reports with species, gene, peptide sequences, and specific mutations

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <chrono>
// No external JSON library needed - we'll write JSON manually

// Structures matching the pipeline
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    uint8_t num_mutations;
    uint8_t mutation_positions[10];
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
    bool used_smith_waterman;
    char query_peptide[51];
    bool is_qrdr_alignment;
};

struct FQMutation {
    uint32_t species_id;
    std::string species_name;
    uint32_t gene_id;
    std::string gene_name;
    uint32_t protein_id;
    int position;           // Global position in protein
    char wildtype_aa;      // Reference amino acid
    char mutant_aa;        // Observed amino acid
    std::string mutation_code;  // e.g., "S87N"
    std::string peptide_sequence;  // Full translated peptide
    int peptide_start;     // Start position in protein
    int peptide_end;       // End position in protein
    float alignment_score;
    float identity;
    uint32_t read_id;
    int frame;
    int occurrences;       // Count of this specific mutation
    bool is_known_resistance;  // Whether this is a known resistance mutation
};

// Known FQ resistance mutations
const std::map<std::string, std::set<std::string>> KNOWN_FQ_MUTATIONS = {
    {"gyrA", {"S83L", "S83F", "S83A", "S83W", "S83Y", "D87N", "D87G", "D87Y", "D87A", "D87V"}},
    {"parC", {"S80I", "S80R", "E84K", "E84V", "E84G"}},
    {"grlA", {"S80F", "S80Y", "E84K", "E84G", "E84L"}},  // S. aureus equivalent
    {"grlB", {"E422D", "D432N"}},  // S. aureus ParE equivalent
    {"gyrB", {"D426N", "K447E"}},
    {"parE", {"D420N", "L445H"}}
};

class FQMutationReporter {
private:
    std::string output_path;
    std::vector<FQMutation> all_mutations;
    std::map<std::string, FQMutation> unique_mutations;  // Key: species_gene_mutation
    std::map<uint32_t, std::string> gene_names;
    std::map<uint32_t, std::string> species_names;
    
    // Statistics
    struct Stats {
        int total_protein_matches = 0;
        int total_mutations_found = 0;
        int unique_mutations = 0;
        int known_resistance_mutations = 0;
        std::map<std::string, int> mutations_by_species;
        std::map<std::string, int> mutations_by_gene;
        std::map<std::string, int> mutation_frequencies;
    } stats;

public:
    FQMutationReporter(const std::string& output) : output_path(output) {
        loadMetadataMappings();
    }
    
    void setGeneMapping(uint32_t id, const std::string& name) {
        gene_names[id] = name;
    }
    
    void setSpeciesMapping(uint32_t id, const std::string& name) {
        species_names[id] = name;
    }
    
    void processProteinMatch(const ProteinMatch& match) {
        stats.total_protein_matches++;
        
        if (match.num_mutations == 0) return;
        
        std::string species_name = getSpeciesName(match.species_id);
        std::string gene_name = getGeneName(match.gene_id);
        
        // Process each mutation
        for (int i = 0; i < match.num_mutations; i++) {
            FQMutation mut;
            mut.species_id = match.species_id;
            mut.species_name = species_name;
            mut.gene_id = match.gene_id;
            mut.gene_name = gene_name;
            mut.protein_id = match.protein_id;
            mut.position = match.ref_start + match.mutation_positions[i];
            mut.wildtype_aa = match.ref_aas[i];
            mut.mutant_aa = match.query_aas[i];
            mut.mutation_code = createMutationCode(mut.wildtype_aa, mut.position, mut.mutant_aa);
            mut.peptide_sequence = std::string(match.query_peptide);
            mut.peptide_start = match.ref_start;
            mut.peptide_end = match.ref_start + match.match_length - 1;
            mut.alignment_score = match.alignment_score;
            mut.identity = match.identity;
            mut.read_id = match.read_id;
            mut.frame = match.frame;
            mut.occurrences = 1;
            mut.is_known_resistance = isKnownResistanceMutation(gene_name, mut.mutation_code);
            
            // Add to collections
            std::string key = species_name + "_" + gene_name + "_" + mut.mutation_code;
            if (unique_mutations.find(key) != unique_mutations.end()) {
                unique_mutations[key].occurrences++;
            } else {
                unique_mutations[key] = mut;
                stats.unique_mutations++;
                if (mut.is_known_resistance) {
                    stats.known_resistance_mutations++;
                }
            }
            
            all_mutations.push_back(mut);
            stats.total_mutations_found++;
            
            // Update statistics
            stats.mutations_by_species[species_name]++;
            stats.mutations_by_gene[gene_name]++;
            stats.mutation_frequencies[mut.mutation_code]++;
        }
    }
    
    void generateReport() {
        // Generate both JSON and text reports
        generateJSONReport();
        generateTextReport();
        generateHTMLReport();
    }

private:
    void loadMetadataMappings() {
        // Try to load from integrated database metadata
        std::ifstream metadata_file("data/integrated_clean_db/protein/metadata.json");
        if (!metadata_file.good()) {
            // Try wildtype database
            metadata_file.open("data/wildtype_protein_db/metadata.json");
        }
        
        if (!metadata_file.good()) {
            std::cerr << "WARNING: Could not load metadata for species/gene names" << std::endl;
            // Use default mappings
            setupDefaultMappings();
            return;
        }
        
        // Simple JSON parsing for species_map and gene_map
        std::string line;
        bool in_species_map = false;
        bool in_gene_map = false;
        
        while (std::getline(metadata_file, line)) {
            // Check for section markers
            if (line.find("\"species_map\"") != std::string::npos) {
                in_species_map = true;
                in_gene_map = false;
                continue;
            }
            if (line.find("\"gene_map\"") != std::string::npos) {
                in_gene_map = true;
                in_species_map = false;
                continue;
            }
            
            // Parse "id": "name" entries
            if ((in_species_map || in_gene_map) && line.find(":") != std::string::npos) {
                size_t first_quote = line.find("\"");
                size_t second_quote = line.find("\"", first_quote + 1);
                size_t third_quote = line.find("\"", second_quote + 1);
                size_t fourth_quote = line.find("\"", third_quote + 1);
                
                if (first_quote != std::string::npos && fourth_quote != std::string::npos) {
                    std::string id_str = line.substr(first_quote + 1, second_quote - first_quote - 1);
                    std::string name = line.substr(third_quote + 1, fourth_quote - third_quote - 1);
                    
                    try {
                        uint32_t id = std::stoi(id_str);
                        if (in_species_map) {
                            species_names[id] = name;
                        } else if (in_gene_map) {
                            gene_names[id] = name;
                        }
                    } catch (...) {
                        // Skip invalid entries
                    }
                }
            }
            
            // Check for section end
            if ((in_species_map || in_gene_map) && line.find("}") != std::string::npos) {
                in_species_map = false;
                in_gene_map = false;
            }
        }
        
        metadata_file.close();
        
        if (species_names.empty() || gene_names.empty()) {
            setupDefaultMappings();
        }
    }
    
    void setupDefaultMappings() {
        // Only set minimal defaults as fallback
        // The actual mappings should come from the database
        std::cerr << "WARNING: Using fallback mappings - database metadata not loaded properly" << std::endl;
        gene_names[0] = "unknown_gene_0";
        gene_names[1] = "unknown_gene_1";
        species_names[0] = "unknown_species_0";
        species_names[1] = "unknown_species_1";
    }
    
    std::string getSpeciesName(uint32_t id) {
        if (species_names.find(id) != species_names.end()) {
            return species_names[id];
        }
        return "Unknown_species_" + std::to_string(id);
    }
    
    std::string getGeneName(uint32_t id) {
        if (gene_names.find(id) != gene_names.end()) {
            return gene_names[id];
        }
        return "Unknown_gene_" + std::to_string(id);
    }
    
    std::string createMutationCode(char wildtype, int position, char mutant) {
        std::stringstream ss;
        ss << wildtype << position << mutant;
        return ss.str();
    }
    
    bool isKnownResistanceMutation(const std::string& gene, const std::string& mutation) {
        auto gene_it = KNOWN_FQ_MUTATIONS.find(gene);
        if (gene_it != KNOWN_FQ_MUTATIONS.end()) {
            return gene_it->second.find(mutation) != gene_it->second.end();
        }
        return false;
    }
    
    void generateJSONReport() {
        std::ofstream json_file(output_path + "_fq_mutations.json");
        
        json_file << "{\n";
        
        // Summary statistics
        json_file << "  \"summary\": {\n";
        json_file << "    \"total_protein_matches\": " << stats.total_protein_matches << ",\n";
        json_file << "    \"total_mutations_found\": " << stats.total_mutations_found << ",\n";
        json_file << "    \"unique_mutations\": " << stats.unique_mutations << ",\n";
        json_file << "    \"known_resistance_mutations\": " << stats.known_resistance_mutations << "\n";
        json_file << "  },\n";
        
        // Detailed mutations
        json_file << "  \"mutations\": [\n";
        bool first = true;
        for (const auto& pair : unique_mutations) {
            const FQMutation& mut = pair.second;
            if (!first) json_file << ",\n";
            json_file << "    {\n";
            json_file << "      \"species\": \"" << mut.species_name << "\",\n";
            json_file << "      \"gene\": \"" << mut.gene_name << "\",\n";
            json_file << "      \"mutation\": \"" << mut.mutation_code << "\",\n";
            json_file << "      \"position\": " << mut.position << ",\n";
            json_file << "      \"wildtype_aa\": \"" << mut.wildtype_aa << "\",\n";
            json_file << "      \"mutant_aa\": \"" << mut.mutant_aa << "\",\n";
            json_file << "      \"peptide_sequence\": \"" << mut.peptide_sequence << "\",\n";
            json_file << "      \"peptide_region\": \"" << mut.peptide_start << "-" << mut.peptide_end << "\",\n";
            json_file << "      \"occurrences\": " << mut.occurrences << ",\n";
            json_file << "      \"is_known_resistance\": " << (mut.is_known_resistance ? "true" : "false") << ",\n";
            json_file << "      \"example_read_id\": " << mut.read_id << ",\n";
            json_file << "      \"example_frame\": " << static_cast<int>(mut.frame) << "\n";
            json_file << "    }";
            first = false;
        }
        json_file << "\n  ],\n";
        
        // Species breakdown
        json_file << "  \"mutations_by_species\": {\n";
        first = true;
        for (const auto& pair : stats.mutations_by_species) {
            if (!first) json_file << ",\n";
            json_file << "    \"" << pair.first << "\": " << pair.second;
            first = false;
        }
        json_file << "\n  },\n";
        
        // Gene breakdown
        json_file << "  \"mutations_by_gene\": {\n";
        first = true;
        for (const auto& pair : stats.mutations_by_gene) {
            if (!first) json_file << ",\n";
            json_file << "    \"" << pair.first << "\": " << pair.second;
            first = false;
        }
        json_file << "\n  }\n";
        
        json_file << "}\n";
        json_file.close();
    }
    
    void generateTextReport() {
        std::ofstream text_file(output_path + "_fq_mutations.txt");
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        
        text_file << "====================================\n";
        text_file << "FLUOROQUINOLONE RESISTANCE MUTATIONS\n";
        text_file << "====================================\n\n";
        text_file << "Generated: " << std::ctime(&time_t);
        text_file << "Total protein matches analyzed: " << stats.total_protein_matches << "\n";
        text_file << "Total mutations found: " << stats.total_mutations_found << "\n";
        text_file << "Unique mutations: " << stats.unique_mutations << "\n";
        text_file << "Known resistance mutations: " << stats.known_resistance_mutations << "\n\n";
        
        // Sort mutations by occurrence
        std::vector<std::pair<std::string, FQMutation>> sorted_mutations;
        for (const auto& pair : unique_mutations) {
            sorted_mutations.push_back(pair);
        }
        std::sort(sorted_mutations.begin(), sorted_mutations.end(),
                  [](const auto& a, const auto& b) {
                      return a.second.occurrences > b.second.occurrences;
                  });
        
        text_file << "DETAILED MUTATION LIST\n";
        text_file << "======================\n\n";
        
        for (const auto& pair : sorted_mutations) {
            const FQMutation& mut = pair.second;
            text_file << "Mutation: " << mut.mutation_code;
            if (mut.is_known_resistance) {
                text_file << " *** KNOWN RESISTANCE MUTATION ***";
            }
            text_file << "\n";
            text_file << "  Species: " << mut.species_name << "\n";
            text_file << "  Gene: " << mut.gene_name << "\n";
            text_file << "  Position: " << mut.position << "\n";
            text_file << "  Change: " << mut.wildtype_aa << " -> " << mut.mutant_aa << "\n";
            text_file << "  Occurrences: " << mut.occurrences << "\n";
            text_file << "  Peptide region: " << mut.peptide_start << "-" << mut.peptide_end << "\n";
            text_file << "  Peptide sequence: " << mut.peptide_sequence << "\n";
            text_file << "  Example read ID: " << mut.read_id << " (frame " << static_cast<int>(mut.frame) << ")\n\n";
        }
        
        text_file << "\nSPECIES SUMMARY\n";
        text_file << "===============\n";
        for (const auto& pair : stats.mutations_by_species) {
            text_file << pair.first << ": " << pair.second << " mutations\n";
        }
        
        text_file << "\nGENE SUMMARY\n";
        text_file << "============\n";
        for (const auto& pair : stats.mutations_by_gene) {
            text_file << pair.first << ": " << pair.second << " mutations\n";
        }
        
        text_file.close();
    }
    
    void generateHTMLReport() {
        std::ofstream html_file(output_path + "_fq_mutations.html");
        
        html_file << "<!DOCTYPE html>\n<html>\n<head>\n";
        html_file << "<title>FQ Resistance Mutations Report</title>\n";
        html_file << "<style>\n";
        html_file << "body { font-family: Arial, sans-serif; margin: 20px; }\n";
        html_file << "table { border-collapse: collapse; width: 100%; margin-top: 20px; }\n";
        html_file << "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n";
        html_file << "th { background-color: #4CAF50; color: white; }\n";
        html_file << "tr:nth-child(even) { background-color: #f2f2f2; }\n";
        html_file << ".known-resistance { background-color: #ffcccc; font-weight: bold; }\n";
        html_file << ".peptide { font-family: monospace; font-size: 12px; }\n";
        html_file << "</style>\n</head>\n<body>\n";
        
        html_file << "<h1>Fluoroquinolone Resistance Mutations Report</h1>\n";
        html_file << "<p>Generated: " << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << "</p>\n";
        
        html_file << "<h2>Summary</h2>\n";
        html_file << "<ul>\n";
        html_file << "<li>Total protein matches: " << stats.total_protein_matches << "</li>\n";
        html_file << "<li>Total mutations: " << stats.total_mutations_found << "</li>\n";
        html_file << "<li>Unique mutations: " << stats.unique_mutations << "</li>\n";
        html_file << "<li>Known resistance mutations: " << stats.known_resistance_mutations << "</li>\n";
        html_file << "</ul>\n";
        
        html_file << "<h2>Mutation Details</h2>\n";
        html_file << "<table>\n";
        html_file << "<tr><th>Mutation</th><th>Species</th><th>Gene</th><th>Position</th>";
        html_file << "<th>Change</th><th>Occurrences</th><th>Peptide Region</th><th>Peptide Sequence</th></tr>\n";
        
        // Sort by occurrence
        std::vector<std::pair<std::string, FQMutation>> sorted_mutations;
        for (const auto& pair : unique_mutations) {
            sorted_mutations.push_back(pair);
        }
        std::sort(sorted_mutations.begin(), sorted_mutations.end(),
                  [](const auto& a, const auto& b) {
                      return a.second.occurrences > b.second.occurrences;
                  });
        
        for (const auto& pair : sorted_mutations) {
            const FQMutation& mut = pair.second;
            html_file << "<tr";
            if (mut.is_known_resistance) {
                html_file << " class='known-resistance'";
            }
            html_file << ">";
            html_file << "<td>" << mut.mutation_code << "</td>";
            html_file << "<td>" << mut.species_name << "</td>";
            html_file << "<td>" << mut.gene_name << "</td>";
            html_file << "<td>" << mut.position << "</td>";
            html_file << "<td>" << mut.wildtype_aa << " → " << mut.mutant_aa << "</td>";
            html_file << "<td>" << mut.occurrences << "</td>";
            html_file << "<td>" << mut.peptide_start << "-" << mut.peptide_end << "</td>";
            html_file << "<td class='peptide'>" << mut.peptide_sequence << "</td>";
            html_file << "</tr>\n";
        }
        
        html_file << "</table>\n";
        html_file << "</body>\n</html>\n";
        html_file.close();
    }
};

// C interface for integration
extern "C" {
    void* create_fq_mutation_reporter(const char* output_path) {
        return new FQMutationReporter(std::string(output_path));
    }
    
    void destroy_fq_mutation_reporter(void* reporter) {
        if (reporter) {
            delete static_cast<FQMutationReporter*>(reporter);
        }
    }
    
    void set_fq_reporter_gene_mapping(void* reporter, uint32_t id, const char* name) {
        if (reporter && name) {
            FQMutationReporter* fq_reporter = static_cast<FQMutationReporter*>(reporter);
            fq_reporter->setGeneMapping(id, std::string(name));
        }
    }
    
    void set_fq_reporter_species_mapping(void* reporter, uint32_t id, const char* name) {
        if (reporter && name) {
            FQMutationReporter* fq_reporter = static_cast<FQMutationReporter*>(reporter);
            fq_reporter->setSpeciesMapping(id, std::string(name));
        }
    }
    
    void process_protein_match_for_fq_report(void* reporter, const ProteinMatch* match) {
        if (reporter && match) {
            static_cast<FQMutationReporter*>(reporter)->processProteinMatch(*match);
        }
    }
    
    void generate_fq_mutation_report(void* reporter) {
        if (reporter) {
            static_cast<FQMutationReporter*>(reporter)->generateReport();
        }
    }
}