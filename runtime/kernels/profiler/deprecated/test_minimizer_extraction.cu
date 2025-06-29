// test_minimizer_extraction.cu
// Focused test for minimizer extraction functionality

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <filesystem>
#include "core/gpu_database_builder_core.h"
#include "gpu_kraken_types.h"

int main(int argc, char** argv) {
    std::cout << "=== Minimizer Extraction Test ===" << std::endl;
    
    // Configuration
    std::string genome_dir = "/home/david/Documents/Code/biogpu/data/type_strain_reference_genomes";
    std::string output_dir = "test_minimizer_output";
    
    // Create output directory
    std::filesystem::create_directories(output_dir);
    
    // Setup configuration with larger buffer to handle our genomes
    DatabaseBuildConfig config;
    config.k_value = 31;
    config.ell_value = 31;
    config.spaces_value = 7;
    
    // Set larger buffer and smaller batch size for successful processing
    config.memory_config.minimizer_capacity = 10000000;     // 10M minimizers
    config.memory_config.sequence_batch_size = 5;           // Process 5 genomes at a time
    config.memory_config.max_memory_fraction = 80;          // Use 80% of GPU memory
    config.memory_config.reserved_memory_mb = 500;          // Reserve 500MB
    
    // Enable debug mode for more output
    config.enable_debug_mode = true;
    config.debug_output_dir = output_dir;
    config.enable_progress_reporting = true;
    
    // Create database builder
    std::cout << "\n1. Creating database builder..." << std::endl;
    GPUKrakenDatabaseBuilder builder(output_dir, config);
    std::cout << "✓ Builder created" << std::endl;
    
    // Build database from first few genomes only
    std::cout << "\n2. Processing genome files..." << std::endl;
    std::cout << "   Batch size: " << config.memory_config.sequence_batch_size << " genomes" << std::endl;
    
    // Get list of FNA files
    std::vector<std::string> fna_files;
    for (const auto& entry : std::filesystem::directory_iterator(genome_dir)) {
        if (entry.path().extension() == ".fna") {
            fna_files.push_back(entry.path().string());
            if (fna_files.size() >= 5) break;  // Just process first 5 files
        }
    }
    
    std::cout << "   Found " << fna_files.size() << " FNA files to process" << std::endl;
    
    // Create a temporary file list
    std::string file_list_path = output_dir + "/test_files.txt";
    std::ofstream file_list(file_list_path);
    for (const auto& file : fna_files) {
        file_list << file << std::endl;
        std::cout << "   - " << std::filesystem::path(file).filename() << std::endl;
    }
    file_list.close();
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Build database from file list
    bool success = builder.build_database_from_file_list(file_list_path);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    if (!success) {
        std::cerr << "\n✗ Database build FAILED" << std::endl;
        return 1;
    }
    
    std::cout << "\n✓ Database build SUCCEEDED in " << duration.count() << " seconds" << std::endl;
    
    // Print statistics
    std::cout << "\n3. Minimizer Extraction Statistics:" << std::endl;
    const auto& stats = builder.get_build_statistics();
    std::cout << "   Total sequences processed: " << stats.total_sequences << std::endl;
    std::cout << "   Total bases processed: " << stats.total_bases << std::endl;
    std::cout << "   Valid minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
    std::cout << "   Unique minimizers: " << stats.unique_minimizers << std::endl;
    std::cout << "   Processing rate: " << stats.processing_rate_mb_per_sec << " MB/s" << std::endl;
    
    // Check minimizer extraction success
    if (stats.valid_minimizers_extracted > 0) {
        std::cout << "\n✓ MINIMIZER EXTRACTION SUCCESSFUL!" << std::endl;
        std::cout << "   Average minimizers per sequence: " 
                  << (stats.total_sequences > 0 ? stats.valid_minimizers_extracted / stats.total_sequences : 0) 
                  << std::endl;
        std::cout << "   Minimizer extraction rate: " 
                  << (stats.total_bases > 0 ? (double)stats.valid_minimizers_extracted / stats.total_bases * 1000 : 0) 
                  << " per kb" << std::endl;
    } else {
        std::cout << "\n✗ NO MINIMIZERS EXTRACTED" << std::endl;
    }
    
    // Check output files
    std::cout << "\n4. Output Files:" << std::endl;
    std::vector<std::string> expected_files = {
        output_dir + "/hash_table.k2d",
        output_dir + "/taxonomy.tsv", 
        output_dir + "/config.txt",
        output_dir + "/build_stats.txt"
    };
    
    for (const auto& file : expected_files) {
        if (std::filesystem::exists(file)) {
            auto size = std::filesystem::file_size(file);
            std::cout << "   ✓ " << std::filesystem::path(file).filename() 
                      << " (" << size << " bytes)" << std::endl;
        } else {
            std::cout << "   ✗ " << std::filesystem::path(file).filename() 
                      << " NOT FOUND" << std::endl;
        }
    }
    
    // Print first few lines of config file to verify parameters
    std::cout << "\n5. Database Configuration:" << std::endl;
    std::ifstream config_file(output_dir + "/config.txt");
    std::string line;
    int line_count = 0;
    while (std::getline(config_file, line) && line_count < 10) {
        std::cout << "   " << line << std::endl;
        line_count++;
    }
    
    // Overall test result
    bool test_passed = success && 
                      stats.total_sequences > 0 && 
                      stats.valid_minimizers_extracted > 0 &&
                      stats.unique_minimizers > 0;
    
    std::cout << "\n=== TEST " << (test_passed ? "PASSED" : "FAILED") << " ===" << std::endl;
    
    if (test_passed) {
        std::cout << "\nSummary:" << std::endl;
        std::cout << "- Successfully processed " << stats.total_sequences << " sequences" << std::endl;
        std::cout << "- Extracted " << stats.valid_minimizers_extracted << " minimizers" << std::endl;
        std::cout << "- Found " << stats.unique_minimizers << " unique minimizers" << std::endl;
        std::cout << "- Database files created successfully" << std::endl;
    }
    
    return test_passed ? 0 : 1;
}