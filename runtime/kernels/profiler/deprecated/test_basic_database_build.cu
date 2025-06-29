// test_basic_database_build.cu
// Test program to verify basic database building functionality

#include "core/gpu_database_builder_core.h"
#include <iostream>
#include <filesystem>
#include <chrono>

int main(int argc, char* argv[]) {
    std::cout << "=== GPU Kraken Database Builder Test ===" << std::endl;
    std::cout << "Testing basic database building functionality\n" << std::endl;
    
    // Set up paths
    std::string test_genomes_dir = "/home/david/Documents/Code/biogpu/data/type_strain_reference_genomes";
    std::string output_db_dir = "./test_db";
    
    // Verify input directory exists
    if (!std::filesystem::exists(test_genomes_dir)) {
        std::cerr << "Error: Test genomes directory not found: " << test_genomes_dir << std::endl;
        return 1;
    }
    
    // Count genome files
    int genome_count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(test_genomes_dir)) {
        if (entry.path().extension() == ".fna") {
            genome_count++;
        }
    }
    std::cout << "Found " << genome_count << " genome files in test directory" << std::endl;
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_db_dir);
    
    // 1. Create builder with simple config
    DatabaseBuildConfig config;
    config.k_value = 31;
    config.ell_value = 31;
    config.spaces_value = 0;  // No spaces for initial test
    config.memory_config.minimizer_capacity = 10000000;  // 10M minimizers
    config.memory_config.sequence_batch_size = 10;  // Process 10 genomes at a time
    config.enable_progress_reporting = true;
    config.validate_output = true;
    
    std::cout << "\nDatabase configuration:" << std::endl;
    std::cout << "  k-mer size: " << config.k_value << std::endl;
    std::cout << "  minimizer window: " << config.ell_value << std::endl;
    std::cout << "  output directory: " << output_db_dir << std::endl;
    
    try {
        // Create the database builder
        GPUKrakenDatabaseBuilder builder(output_db_dir, config);
        
        std::cout << "\nStarting database build..." << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // 2. Build from test genomes (no taxonomy for initial test)
        bool success = builder.build_database_from_genomes(
            test_genomes_dir,  // Directory with FASTA files
            ""  // No taxonomy for initial test
        );
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        if (success) {
            std::cout << "\n✅ Database built successfully!" << std::endl;
            std::cout << "Total build time: " << duration.count() << " seconds" << std::endl;
            
            // Print build statistics
            builder.print_build_progress();
            builder.print_memory_usage();
            
            // Check output files
            std::cout << "\nChecking output files:" << std::endl;
            std::vector<std::string> expected_files = {
                "hash_table.k2d",
                "taxonomy.tsv",
                "config.txt"
            };
            
            for (const auto& file : expected_files) {
                std::string full_path = output_db_dir + "/" + file;
                if (std::filesystem::exists(full_path)) {
                    auto file_size = std::filesystem::file_size(full_path);
                    std::cout << "  ✓ " << file << " (" << file_size << " bytes)" << std::endl;
                } else {
                    std::cout << "  ✗ " << file << " (missing)" << std::endl;
                }
            }
            
        } else {
            std::cerr << "\n❌ Database build failed!" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=== Test completed ===" << std::endl;
    return 0;
}