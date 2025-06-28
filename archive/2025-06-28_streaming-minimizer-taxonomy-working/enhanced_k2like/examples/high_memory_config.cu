#include "core/gpu_database_builder_core.h"

int main() {
    // Create custom high-memory configuration
    DatabaseBuildConfig config = GPUKrakenDatabaseBuilder::create_high_memory_config();
    config.memory_config.minimizer_capacity = 25000000;  // 25M minimizers
    config.memory_config.sequence_batch_size = 50;       // Large batches
    config.memory_config.max_memory_fraction = 90;       // Use 90% of GPU memory
    
    auto builder = DatabaseBuilderFactory::create_custom_builder("/output/path", config);
    
    // Build with high capacity
    bool success = builder->build_database_from_genomes("/large/genome/collection");
    
    return success ? 0 : 1;
}