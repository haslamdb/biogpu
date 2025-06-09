// pipeline_integration_example.cpp
// Example of how to integrate the global FQ resistance mapper into your pipeline

#include <iostream>
#include <string>
#include "global_fq_resistance_mapper.h"

// Example integration into clean_resistance_pipeline_main.cpp
class IntegratedCleanResistancePipeline {
private:
    // ... existing member variables ...
    
    // Add reference to global mapper
    GlobalFQResistanceMapper& fq_mapper;
    
public:
    IntegratedCleanResistancePipeline() 
        : fq_mapper(GlobalFQResistanceMapper::getInstance()) {
        
        // Initialize the global FQ mapper early in pipeline startup
        std::string fq_csv_path = "data/quinolone_resistance_mutation_table.csv";
        std::string protein_db_path = "data/integrated_clean_db/protein";
        
        int result = init_global_fq_mapper(fq_csv_path.c_str(), protein_db_path.c_str());
        if (result != 0) {
            std::cerr << "WARNING: Failed to initialize FQ resistance mapper" << std::endl;
            std::cerr << "FQ resistance detection will not be available" << std::endl;
        } else {
            std::cout << "FQ resistance mapper initialized successfully" << std::endl;
            fq_mapper.printSummary();
        }
    }
    
    // Modified protein match processing
    void processProteinMatch(const ProteinMatch& match) {
        // Get species and gene names from match
        auto species_gene = fq_mapper.getSpeciesGeneFromID(match.species_gene_id);
        if (species_gene.first.empty() || species_gene.second.empty()) {
            // Fallback to using individual IDs if combined ID doesn't work
            std::string species = getSpeciesName(match.species_id);
            std::string gene = getGeneName(match.gene_id);
            species_gene = {species, gene};
        }
        
        std::cout << "Processing match for " << species_gene.first 
                  << " " << species_gene.second << std::endl;
        
        // Check each mutation
        int fq_mutations_found = 0;
        for (int i = 0; i < match.num_mutations; i++) {
            uint16_t position = match.ref_start + match.mutation_positions[i] + 1; // 1-based
            char ref_aa = match.ref_aas[i];
            char query_aa = match.query_aas[i];
            
            // Check if this is a known FQ resistance mutation
            bool is_fq = fq_mapper.isResistanceMutation(
                species_gene.first, species_gene.second,
                position, ref_aa, query_aa
            );
            
            if (is_fq) {
                fq_mutations_found++;
                std::cout << "  FQ RESISTANCE DETECTED: " << species_gene.second 
                         << " " << ref_aa << position << query_aa << std::endl;
            } else {
                // Check if this is a resistance position with a different mutation
                char expected_wt = fq_mapper.getWildtypeAA(
                    species_gene.first, species_gene.second, position
                );
                
                if (expected_wt != 'X' && expected_wt != ref_aa) {
                    std::cout << "  Note: Position " << position 
                             << " expected wildtype " << expected_wt 
                             << " but reference has " << ref_aa << std::endl;
                }
            }
        }
        
        // Report results
        if (fq_mutations_found > 0) {
            stats.resistance_mutations += fq_mutations_found;
            
            // Add to JSON output
            json_output << "    {\n";
            json_output << "      \"type\": \"FQ_RESISTANCE\",\n";
            json_output << "      \"species\": \"" << species_gene.first << "\",\n";
            json_output << "      \"gene\": \"" << species_gene.second << "\",\n";
            json_output << "      \"mutations\": " << fq_mutations_found << ",\n";
            json_output << "      \"details\": [\n";
            
            // List specific mutations
            bool first = true;
            for (int i = 0; i < match.num_mutations; i++) {
                uint16_t position = match.ref_start + match.mutation_positions[i] + 1;
                char ref_aa = match.ref_aas[i];
                char query_aa = match.query_aas[i];
                
                if (fq_mapper.isResistanceMutation(species_gene.first, species_gene.second,
                                                  position, ref_aa, query_aa)) {
                    if (!first) json_output << ",\n";
                    json_output << "        {\"position\": " << position 
                               << ", \"change\": \"" << ref_aa << "->" << query_aa << "\"}";
                    first = false;
                }
            }
            
            json_output << "\n      ]\n";
            json_output << "    }";
        }
    }
    
    // Helper to get normalized species name
    std::string getSpeciesName(uint32_t species_id) {
        // First try the global mapper
        auto species_gene = fq_mapper.getSpeciesGeneFromID(species_id);
        if (!species_gene.first.empty()) {
            return species_gene.first;
        }
        
        // Fallback to existing mapping
        auto it = species_id_to_name.find(species_id);
        if (it != species_id_to_name.end()) {
            return it->second;
        }
        return "Species" + std::to_string(species_id);
    }
    
    // Helper to get normalized gene name
    std::string getGeneName(uint32_t gene_id) {
        // Use existing mapping
        auto it = gene_id_to_name.find(gene_id);
        if (it != gene_id_to_name.end()) {
            return it->second;
        }
        return "Gene" + std::to_string(gene_id);
    }
};

// Example for diagnostic report integration
void updateDiagnosticReportWithFQ(void* reporter, const ProteinMatch* match) {
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Get species and gene info
    auto species_gene = mapper.getSpeciesGeneFromID(match->species_gene_id);
    
    // Build detailed FQ resistance report
    std::stringstream fq_report;
    fq_report << "\nFQ RESISTANCE ANALYSIS:\n";
    
    bool has_fq_mutations = false;
    for (int i = 0; i < match->num_mutations; i++) {
        uint16_t position = match->ref_start + match->mutation_positions[i] + 1;
        
        if (mapper.isResistanceMutation(species_gene.first, species_gene.second,
                                       position, match->ref_aas[i], match->query_aas[i])) {
            has_fq_mutations = true;
            fq_report << "  - " << species_gene.second << " " 
                     << match->ref_aas[i] << position << match->query_aas[i]
                     << " [FLUOROQUINOLONE RESISTANCE]\n";
            
            // Get all known resistant variants at this position
            auto resistant_aas = mapper.getResistantAAs(species_gene.first, 
                                                       species_gene.second, position);
            fq_report << "    Known resistant variants at position " << position << ": ";
            for (char aa : resistant_aas) {
                fq_report << aa << " ";
            }
            fq_report << "\n";
        }
    }
    
    if (!has_fq_mutations) {
        fq_report << "  No known FQ resistance mutations detected\n";
    }
    
    // Add to diagnostic report
    // This would be integrated with your existing diagnostic reporter
    std::cout << fq_report.str();
}

// Main function showing complete integration
int main(int argc, char** argv) {
    // Initialize global FQ mapper once at program start
    const char* fq_csv = "data/quinolone_resistance_mutation_table.csv";
    const char* protein_db = "data/integrated_clean_db/protein";
    
    if (init_global_fq_mapper(fq_csv, protein_db) != 0) {
        std::cerr << "Failed to initialize FQ resistance database" << std::endl;
        return 1;
    }
    
    // Get instance for direct use
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Example: Check a specific mutation
    bool is_resistant = mapper.isResistanceMutation(
        "Escherichia_coli", "gyrA", 83, 'S', 'L'
    );
    std::cout << "E. coli gyrA S83L is " 
              << (is_resistant ? "RESISTANT" : "not resistant") << std::endl;
    
    // Example: Get all mutations for a gene
    auto ecoli_gyra_mutations = mapper.getMutationsForGene("Escherichia_coli", "gyrA");
    std::cout << "\nE. coli gyrA has " << ecoli_gyra_mutations.size() 
              << " known resistance mutations:" << std::endl;
    
    for (const auto& mut : ecoli_gyra_mutations) {
        std::cout << "  " << mut.wildtype_aa << mut.position << mut.mutant_aa << std::endl;
    }
    
    // Print summary
    mapper.printSummary();
    
    // Clean up (handled automatically by singleton destructor)
    return 0;
}