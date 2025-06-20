#include "core/gpu_database_builder_core.h"
#include <iostream>

int main() {
    // Create builder with default configuration
    auto builder = DatabaseBuilderFactory::create_standard_builder("/output/path");
    
    // Build database from genome directory
    bool success = builder->build_database_from_genomes(
        "/path/to/genomes",     // Genome directory
        "/path/to/taxonomy"     // NCBI taxonomy (optional)
    );
    
    if (success) {
        builder->print_build_progress();
        std::cout << "Database built successfully!" << std::endl;
    }
    
    return success ? 0 : 1;
}