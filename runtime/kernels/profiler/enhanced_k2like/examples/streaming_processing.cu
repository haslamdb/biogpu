#include "core/gpu_database_builder_core.h"

int main() {
    // Create streaming-optimized builder
    auto builder = DatabaseBuilderFactory::create_streaming_builder("/output/path");
    
    // Build from very large FNA file using streaming
    bool success = builder->build_database_from_streaming_fna(
        "/path/to/huge_database.fna",  // Large input file (100GB+)
        "/path/to/taxonomy"            // Optional taxonomy
    );
    
    return success ? 0 : 1;
}