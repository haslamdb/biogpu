// test_fna_simple.cu
// Simple test to verify FNA processing and minimizer extraction works

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include "core/gpu_database_builder_core.h"
#include "gpu_kraken_types.h"

int main(int argc, char** argv) {
    std::cout << "=== Simple FNA Minimizer Extraction Test ===" << std::endl;
    
    // Configuration
    std::string fna_file = "data/type_strain_reference_genomes/library.fna";
    if (argc > 1) {
        fna_file = argv[1];
    }
    
    std::string output_dir = "test_output_fna";
    
    // Create output directory
    std::string mkdir_cmd = "mkdir -p " + output_dir;
    system(mkdir_cmd.c_str());
    
    // Setup configuration
    DatabaseBuildConfig config;
    config.k_value = 31;
    config.ell_value = 31;
    config.spaces_value = 7;
    
    // Memory configuration - key settings for our buffer management
    config.memory_config.minimizer_capacity = 10000000;            // 10M minimizers
    config.memory_config.sequence_batch_size = 10;                 // Process 10 genomes at a time
    config.memory_config.max_memory_fraction = 80;                 // Use 80% of GPU memory
    config.memory_config.reserved_memory_mb = 500;                 // Reserve 500MB
    
    config.enable_debug_mode = true;
    config.debug_output_dir = output_dir;
    
    // Create database builder
    std::cout << "\n1. Creating database builder..." << std::endl;
    GPUKrakenDatabaseBuilder builder(output_dir, config);
    
    // Test basic functionality first
    if (!builder.is_initialized()) {
        std::cerr << "Failed to initialize database builder" << std::endl;
        return 1;
    }
    std::cout << "✓ Builder initialized" << std::endl;
    
    // Build database from streaming FNA
    std::cout << "\n2. Building database from FNA file: " << fna_file << std::endl;
    std::cout << "   Using 50MB buffer with 90% threshold = " << (0.9 * 50) << "MB" << std::endl;
    std::cout << "   Batch size: " << config.memory_config.sequence_batch_size << " genomes" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    bool success = builder.build_database_from_streaming_fna(fna_file);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    if (!success) {
        std::cerr << "\n✗ Database build FAILED" << std::endl;
        return 1;
    }
    
    std::cout << "\n✓ Database build SUCCEEDED in " << duration.count() << " seconds" << std::endl;
    
    // Print statistics
    std::cout << "\n3. Build Statistics:" << std::endl;
    const auto& stats = builder.get_build_statistics();
    std::cout << "   Total sequences processed: " << stats.total_sequences << std::endl;
    std::cout << "   Total bases processed: " << stats.total_bases << std::endl;
    std::cout << "   Valid minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
    std::cout << "   Unique minimizers: " << stats.unique_minimizers << std::endl;
    std::cout << "   Build time: " << stats.build_time << " seconds" << std::endl;
    
    // Check if output files were created
    std::cout << "\n4. Checking output files..." << std::endl;
    std::string hash_file = output_dir + "/hash_table.k2d";
    std::string taxo_file = output_dir + "/taxonomy.tsv";
    std::string config_file = output_dir + "/config.txt";
    
    std::ifstream hf(hash_file), tf(taxo_file), cf(config_file);
    if (hf.good()) std::cout << "   ✓ " << hash_file << " created" << std::endl;
    else std::cout << "   ✗ " << hash_file << " NOT created" << std::endl;
    
    if (tf.good()) std::cout << "   ✓ " << taxo_file << " created" << std::endl;
    else std::cout << "   ✗ " << taxo_file << " NOT created" << std::endl;
    
    if (cf.good()) std::cout << "   ✓ " << config_file << " created" << std::endl;
    else std::cout << "   ✗ " << config_file << " NOT created" << std::endl;
    
    // Success criteria
    bool test_passed = success && 
                      stats.total_sequences > 0 && 
                      stats.valid_minimizers_extracted > 0;
    
    std::cout << "\n=== TEST " << (test_passed ? "PASSED" : "FAILED") << " ===" << std::endl;
    
    return test_passed ? 0 : 1;
}