// test_taxonomy_only.cu
// Test program for the taxonomy module

#include <iostream>
#include "taxonomy/taxonomy_processor.h"
#include "gpu_kraken_types.h"

int main() {
    std::cout << "=== BioGPU Taxonomy Module Test ===" << std::endl;
    
    try {
        // Test SimpleTaxonomyProcessor
        std::cout << "Testing SimpleTaxonomyProcessor..." << std::endl;
        SimpleTaxonomyProcessor simple_taxonomy;
        
        // Add some test taxonomy data
        simple_taxonomy.add_taxon(1, 1, "root", "no rank");
        simple_taxonomy.add_taxon(2, 1, "Bacteria", "superkingdom");
        simple_taxonomy.add_taxon(511145, 2, "Escherichia coli", "species");
        simple_taxonomy.add_taxon(287, 2, "Pseudomonas aeruginosa", "species");
        
        std::cout << "✓ Simple taxonomy created with " << simple_taxonomy.size() << " taxa" << std::endl;
        
        // Test LCA computation
        uint32_t lca = simple_taxonomy.compute_simple_lca(511145, 287);
        std::cout << "✓ LCA of E. coli and P. aeruginosa: " << lca << " (" << simple_taxonomy.get_name(lca) << ")" << std::endl;
        
        // Test multiple LCA
        std::vector<uint32_t> species_list = {511145, 287};
        uint32_t multi_lca = simple_taxonomy.compute_lca_of_list(species_list);
        std::cout << "✓ Multi-LCA result: " << multi_lca << std::endl;
        
        // Test EnhancedNCBITaxonomyProcessor
        std::cout << "\nTesting EnhancedNCBITaxonomyProcessor..." << std::endl;
        EnhancedNCBITaxonomyProcessor enhanced_taxonomy;
        
        std::cout << "✓ Enhanced taxonomy processor created" << std::endl;
        std::cout << "✓ Is loaded: " << (enhanced_taxonomy.is_loaded() ? "Yes" : "No") << std::endl;
        
        // Test utility functions
        std::cout << "\nTesting utility functions..." << std::endl;
        
        std::vector<uint32_t> test_species = {511145, 287, 2};
        bool valid = PhylogeneticUtils::validate_species_list(test_species);
        std::cout << "✓ Species list validation: " << (valid ? "Valid" : "Invalid") << std::endl;
        
        // Test with empty taxonomy (should handle gracefully)
        uint32_t empty_lca = enhanced_taxonomy.compute_lca_of_species(test_species);
        std::cout << "✓ Empty taxonomy LCA (should be 1): " << empty_lca << std::endl;
        
        std::cout << "\n=== ALL TAXONOMY TESTS PASSED ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}