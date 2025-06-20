#include "core/gpu_database_builder_core.h"
#include <iostream>

int main() {
    // Create builder optimized for phylogenetic analysis
    auto builder = DatabaseBuilderFactory::create_phylogenetic_builder("/output/path");
    
    // Configure for enhanced processing
    builder->configure_phylogenetic_analysis(true);
    builder->configure_output_format(DatabaseFormat::ENHANCED_PHYLO, true);
    
    // Build from concatenated FNA with taxonomy
    bool success = builder->build_database_from_concatenated_fna(
        "/path/to/large_file.fna",          // Input FNA
        "/path/to/nodes.dmp",               // NCBI nodes
        "/path/to/names.dmp",               // NCBI names
        "/path/to/compact_taxonomy.bin"     // Pre-built compact taxonomy (optional)
    );
    
    if (success) {
        const auto& enhanced_stats = builder->get_enhanced_statistics();
        enhanced_stats.print_enhanced_stats();
    }
    
    return success ? 0 : 1;
}