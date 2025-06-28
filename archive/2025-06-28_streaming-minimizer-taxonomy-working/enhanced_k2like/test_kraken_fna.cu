// test_kraken_fna.cu
// Test program to verify Kraken-format FNA parsing

#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include "processing/genome_file_processor.h"
#include "gpu_kraken_types.h"

int main() {
    std::cout << "=== Testing Kraken-format FNA Parsing ===" << std::endl;
    
    // Configure processing
    FileProcessingConfig config;
    config.progress_reporting = true;
    config.max_file_count = 1000;  // Limit for testing
    config.progress_interval = 10;
    
    // Process the test file
    std::cout << "\nProcessing Kraken-format FNA file..." << std::endl;
    ConcatenatedFnaProcessor processor("test_kraken_format_larger.fna", config);
    
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    
    auto start = std::chrono::high_resolution_clock::now();
    bool success = processor.process_fna_file(genome_files, genome_taxon_ids, taxon_names, "./temp_genomes");
    auto end = std::chrono::high_resolution_clock::now();
    
    if (success) {
        std::cout << "\n✓ Processing completed successfully!" << std::endl;
        std::cout << "  Processing time: " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
                  << " ms" << std::endl;
        std::cout << "  Genomes found: " << genome_files.size() << std::endl;
        
        // Show unique taxon IDs
        std::set<uint32_t> unique_taxons(genome_taxon_ids.begin(), genome_taxon_ids.end());
        std::cout << "  Unique taxon IDs: " << unique_taxons.size() << std::endl;
        
        // Show first few taxon IDs
        std::cout << "  First 10 taxon IDs: ";
        size_t count = 0;
        for (uint32_t taxid : unique_taxons) {
            std::cout << taxid << " ";
            if (++count >= 10) break;
        }
        std::cout << std::endl;
        
        // Check a sample file for sequence integrity
        if (!genome_files.empty()) {
            std::cout << "\n  Checking first genome file..." << std::endl;
            std::ifstream check_file(genome_files[0]);
            if (check_file.is_open()) {
                std::string header, sequence;
                std::getline(check_file, header);
                
                // Read all sequence lines
                std::string line;
                size_t line_count = 0;
                while (std::getline(check_file, line)) {
                    sequence += line;
                    line_count++;
                }
                check_file.close();
                
                std::cout << "    File: " << genome_files[0] << std::endl;
                std::cout << "    Header: " << header.substr(0, 80) << "..." << std::endl;
                std::cout << "    Sequence length: " << sequence.length() << " bases" << std::endl;
                std::cout << "    Sequence lines: " << line_count << std::endl;
                
                // Check for whitespace
                bool has_whitespace = false;
                for (char c : sequence) {
                    if (std::isspace(c) || c == '\r' || c == '\n' || c == '\t') {
                        has_whitespace = true;
                        break;
                    }
                }
                std::cout << "    Contains whitespace: " << (has_whitespace ? "YES (ERROR!)" : "NO (good)") << std::endl;
                
                // Show first and last 30 bases
                if (sequence.length() >= 60) {
                    std::cout << "    First 30 bases: " << sequence.substr(0, 30) << std::endl;
                    std::cout << "    Last 30 bases:  " << sequence.substr(sequence.length() - 30, 30) << std::endl;
                }
            }
        }
        
        // Show species data statistics
        const auto& species_data = processor.get_species_data();
        std::cout << "\n  Species tracking statistics:" << std::endl;
        std::cout << "    Total species: " << species_data.total_species() << std::endl;
        std::cout << "    Total genomes tracked: " << species_data.total_genomes() << std::endl;
        
    } else {
        std::cerr << "✗ Processing failed!" << std::endl;
    }
    
    // Cleanup temporary directory
    try {
        std::filesystem::remove_all("./temp_genomes");
    } catch (...) {}
    
    return success ? 0 : 1;
}