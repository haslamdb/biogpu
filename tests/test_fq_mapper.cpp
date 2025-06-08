// test_fq_mapper.cpp
// Test program for the global FQ resistance mapper

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include "global_fq_resistance_mapper.h"

void testSpecificMutations() {
    std::cout << "\n=== Testing Specific Known FQ Mutations ===" << std::endl;
    
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Test cases from the CSV
    struct TestCase {
        std::string species;
        std::string gene;
        uint16_t position;
        char wildtype;
        char mutant;
        bool expected_result;
    };
    
    std::vector<TestCase> test_cases = {
        // Known resistance mutations
        {"Escherichia_coli", "gyrA", 83, 'S', 'L', true},
        {"Escherichia_coli", "gyrA", 83, 'S', 'A', true},
        {"Escherichia_coli", "gyrA", 87, 'D', 'N', true},
        {"Escherichia_coli", "parC", 80, 'S', 'I', true},
        {"Escherichia_coli", "parE", 458, 'S', 'A', true},
        
        // Non-resistant mutations (wildtype)
        {"Escherichia_coli", "gyrA", 83, 'S', 'S', false},
        {"Escherichia_coli", "gyrA", 87, 'D', 'D', false},
        
        // Unknown mutations
        {"Escherichia_coli", "gyrA", 83, 'S', 'K', false},  // Not in CSV
        
        // Different species
        {"Staphylococcus_aureus", "gyrA", 84, 'S', 'L', true},
        {"Enterococcus_faecium", "gyrA", 83, 'S', 'I', true},
        {"Enterococcus_faecium", "parC", 80, 'S', 'I', true},
    };
    
    int passed = 0;
    int failed = 0;
    
    for (const auto& test : test_cases) {
        bool result = mapper.isResistanceMutation(
            test.species, test.gene, test.position, test.wildtype, test.mutant
        );
        
        std::cout << std::setw(25) << test.species 
                  << " " << std::setw(5) << test.gene 
                  << " " << test.wildtype << std::setw(3) << test.position << test.mutant
                  << " : " << (result ? "RESISTANT" : "not resistant")
                  << " [" << (result == test.expected_result ? "PASS" : "FAIL") << "]" << std::endl;
        
        if (result == test.expected_result) {
            passed++;
        } else {
            failed++;
        }
    }
    
    std::cout << "\nTest results: " << passed << " passed, " << failed << " failed" << std::endl;
}

void testSpeciesNormalization() {
    std::cout << "\n=== Testing Species Name Normalization ===" << std::endl;
    
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Test with both space and underscore versions
    std::vector<std::pair<std::string, std::string>> species_pairs = {
        {"Escherichia coli", "Escherichia_coli"},
        {"Staphylococcus aureus", "Staphylococcus_aureus"},
        {"Enterococcus faecium", "Enterococcus_faecium"},
    };
    
    for (const auto& pair : species_pairs) {
        // Test mutation lookup with space version
        bool space_result = mapper.isResistanceMutation(
            pair.first, "gyrA", 83, 'S', 'L'
        );
        
        // Test mutation lookup with underscore version
        bool underscore_result = mapper.isResistanceMutation(
            pair.second, "gyrA", 83, 'S', 'L'
        );
        
        std::cout << pair.first << " vs " << pair.second << ": "
                  << (space_result == underscore_result ? "NORMALIZED CORRECTLY" : "NORMALIZATION FAILED")
                  << std::endl;
    }
}

void testWildtypeRetrieval() {
    std::cout << "\n=== Testing Wildtype AA Retrieval ===" << std::endl;
    
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Test known positions
    struct WildtypeTest {
        std::string species;
        std::string gene;
        uint16_t position;
        char expected_wt;
    };
    
    std::vector<WildtypeTest> tests = {
        {"Escherichia_coli", "gyrA", 83, 'S'},
        {"Escherichia_coli", "gyrA", 87, 'D'},
        {"Escherichia_coli", "parC", 80, 'S'},
        {"Staphylococcus_aureus", "gyrA", 84, 'S'},
        {"Enterococcus_faecium", "parC", 84, 'E'},
    };
    
    for (const auto& test : tests) {
        char wt = mapper.getWildtypeAA(test.species, test.gene, test.position);
        std::cout << test.species << " " << test.gene << " position " << test.position
                  << ": wildtype = " << wt 
                  << " [" << (wt == test.expected_wt ? "CORRECT" : "INCORRECT") << "]" << std::endl;
    }
}

void testResistantAAsList() {
    std::cout << "\n=== Testing Resistant AAs Retrieval ===" << std::endl;
    
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Get all resistant AAs for E. coli gyrA position 83
    auto resistant_aas = mapper.getResistantAAs("Escherichia_coli", "gyrA", 83);
    
    std::cout << "E. coli gyrA position 83 resistant variants: ";
    for (char aa : resistant_aas) {
        std::cout << aa << " ";
    }
    std::cout << std::endl;
    
    // Should include at least L, A based on CSV
    bool has_L = std::find(resistant_aas.begin(), resistant_aas.end(), 'L') != resistant_aas.end();
    bool has_A = std::find(resistant_aas.begin(), resistant_aas.end(), 'A') != resistant_aas.end();
    
    if (has_L && has_A) {
        std::cout << "Correctly includes known resistant variants L and A" << std::endl;
    } else {
        std::cout << "ERROR: Missing expected resistant variants" << std::endl;
    }
}

void testGPUDataPreparation() {
    std::cout << "\n=== Testing GPU Data Preparation ===" << std::endl;
    
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    if (mapper.prepareGPUData()) {
        std::cout << "GPU data prepared successfully" << std::endl;
        std::cout << "Number of GPU mutations: " << mapper.getNumGPUMutations() << std::endl;
        
        void* gpu_ptr = mapper.getGPUMutations();
        std::cout << "GPU data pointer: " << (gpu_ptr ? "VALID" : "NULL") << std::endl;
    } else {
        std::cout << "ERROR: Failed to prepare GPU data" << std::endl;
    }
}

void printStatsBySpecies() {
    std::cout << "\n=== FQ Mutations by Species ===" << std::endl;
    
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    
    // Common species to check
    std::vector<std::string> species_list = {
        "Escherichia_coli",
        "Staphylococcus_aureus",
        "Enterococcus_faecium",
        "Pseudomonas_aeruginosa",
        "Klebsiella_pneumoniae"
    };
    
    for (const auto& species : species_list) {
        std::cout << "\n" << species << ":" << std::endl;
        
        std::vector<std::string> genes = {"gyrA", "gyrB", "parC", "parE"};
        for (const auto& gene : genes) {
            auto mutations = mapper.getMutationsForGene(species, gene);
            if (!mutations.empty()) {
                std::cout << "  " << gene << ": " << mutations.size() << " mutations" << std::endl;
                
                // Show first few
                int count = 0;
                for (const auto& mut : mutations) {
                    std::cout << "    " << mut.wildtype_aa << mut.position << mut.mutant_aa;
                    if (++count >= 3) {
                        if (mutations.size() > 3) std::cout << " ...";
                        break;
                    }
                    if (count < mutations.size()) std::cout << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fq_mutations.csv> <protein_db_path>" << std::endl;
        return 1;
    }
    
    std::string csv_path = argv[1];
    std::string protein_db_path = argv[2];
    
    std::cout << "=== FQ Resistance Mapper Test Suite ===" << std::endl;
    std::cout << "CSV file: " << csv_path << std::endl;
    std::cout << "Protein DB: " << protein_db_path << std::endl;
    
    // Initialize mapper
    if (init_global_fq_mapper(csv_path.c_str(), protein_db_path.c_str()) != 0) {
        std::cerr << "ERROR: Failed to initialize FQ mapper" << std::endl;
        return 1;
    }
    
    // Run tests
    testSpecificMutations();
    testSpeciesNormalization();
    testWildtypeRetrieval();
    testResistantAAsList();
    testGPUDataPreparation();
    printStatsBySpecies();
    
    // Print summary
    GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
    mapper.printSummary();
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    
    return 0;
}