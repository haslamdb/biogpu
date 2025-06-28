// test_fna_minimizer_extraction.cpp
// Test program to verify FNA file reading and minimizer extraction

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <string>
#include <cstring>

// Include your headers
#include "runtime/kernels/profiler/enhanced_k2like/core/gpu_database_builder_core.h"
#include "runtime/kernels/profiler/enhanced_k2like/gpu_kraken_types.h"

class FNAMinimizerTest {
private:
    std::string fna_file_path_;
    std::string output_dir_;
    
public:
    FNAMinimizerTest(const std::string& fna_file = "data/type_strain_reference_genomes/library.fna") 
        : fna_file_path_(fna_file), output_dir_("test_output") {
        
        // Create output directory
        std::filesystem::create_directories(output_dir_);
    }
    
    void run_all_tests() {
        std::cout << "\n=== FNA MINIMIZER EXTRACTION TEST SUITE ===\n" << std::endl;
        
        // Test 1: Verify FNA file exists and is readable
        if (!test_fna_file_readable()) {
            std::cerr << "FATAL: Cannot read FNA file: " << fna_file_path_ << std::endl;
            return;
        }
        
        // Test 2: Analyze FNA file structure
        test_analyze_fna_structure();
        
        // Test 3: Test the full pipeline with small subset
        test_minimizer_extraction_pipeline();
        
        // Test 4: Memory usage analysis
        test_memory_usage();
        
        // Test 5: Minimizer statistics
        test_minimizer_statistics();
        
        std::cout << "\n=== ALL TESTS COMPLETE ===\n" << std::endl;
    }
    
private:
    bool test_fna_file_readable() {
        std::cout << "Test 1: Checking FNA file accessibility..." << std::endl;
        
        if (!std::filesystem::exists(fna_file_path_)) {
            std::cerr << "  ERROR: File does not exist: " << fna_file_path_ << std::endl;
            return false;
        }
        
        auto file_size = std::filesystem::file_size(fna_file_path_);
        std::cout << "  ✓ File exists" << std::endl;
        std::cout << "  ✓ File size: " << (file_size / 1024 / 1024) << " MB" << std::endl;
        
        // Try to read first few lines
        std::ifstream file(fna_file_path_);
        if (!file.is_open()) {
            std::cerr << "  ERROR: Cannot open file for reading" << std::endl;
            return false;
        }
        
        std::string line;
        int lines_read = 0;
        while (std::getline(file, line) && lines_read < 5) {
            if (lines_read == 0) {
                std::cout << "  ✓ First header: " << line.substr(0, std::min(size_t(80), line.length())) 
                          << (line.length() > 80 ? "..." : "") << std::endl;
            }
            lines_read++;
        }
        
        file.close();
        return true;
    }
    
    void test_analyze_fna_structure() {
        std::cout << "\nTest 2: Analyzing FNA file structure..." << std::endl;
        
        std::ifstream file(fna_file_path_);
        if (!file.is_open()) return;
        
        int sequence_count = 0;
        int line_count = 0;
        size_t total_bases = 0;
        std::unordered_set<uint32_t> unique_taxids;
        std::vector<std::pair<uint32_t, std::string>> first_sequences;
        
        std::string line;
        uint32_t current_taxid = 0;
        std::string current_name;
        size_t current_seq_length = 0;
        
        while (std::getline(file, line)) {
            line_count++;
            
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Save previous sequence info
                if (current_taxid > 0 && sequence_count < 5) {
                    first_sequences.push_back({current_taxid, current_name + " (" + 
                                              std::to_string(current_seq_length) + " bp)"});
                }
                
                sequence_count++;
                current_seq_length = 0;
                
                // Parse header
                size_t taxid_pos = line.find("taxid|");
                if (taxid_pos != std::string::npos) {
                    size_t taxid_start = taxid_pos + 6;
                    size_t taxid_end = line.find('|', taxid_start);
                    if (taxid_end != std::string::npos) {
                        std::string taxid_str = line.substr(taxid_start, taxid_end - taxid_start);
                        try {
                            current_taxid = std::stoul(taxid_str);
                            unique_taxids.insert(current_taxid);
                            
                            // Extract species name
                            size_t name_start = line.find(' ', taxid_end);
                            if (name_start != std::string::npos) {
                                current_name = line.substr(name_start + 1);
                            }
                        } catch (...) {
                            current_taxid = 0;
                        }
                    }
                }
            } else {
                // Count bases in sequence line
                for (char c : line) {
                    if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
                        total_bases++;
                        current_seq_length++;
                    }
                }
            }
            
            // Progress report
            if (line_count % 1000000 == 0) {
                std::cout << "  Processed " << line_count << " lines..." << std::endl;
            }
        }
        
        // Save last sequence
        if (current_taxid > 0 && first_sequences.size() < 5) {
            first_sequences.push_back({current_taxid, current_name + " (" + 
                                      std::to_string(current_seq_length) + " bp)"});
        }
        
        file.close();
        
        // Report findings
        std::cout << "  ✓ Total sequences: " << sequence_count << std::endl;
        std::cout << "  ✓ Unique taxon IDs: " << unique_taxids.size() << std::endl;
        std::cout << "  ✓ Total bases: " << (total_bases / 1000000) << " Mbp" << std::endl;
        std::cout << "  ✓ Average sequence length: " << (total_bases / sequence_count) << " bp" << std::endl;
        
        std::cout << "\n  First few sequences:" << std::endl;
        for (const auto& [taxid, name] : first_sequences) {
            std::cout << "    Taxon " << taxid << ": " << name << std::endl;
        }
    }
    
    void test_minimizer_extraction_pipeline() {
        std::cout << "\nTest 3: Testing minimizer extraction pipeline..." << std::endl;
        
        // Create a small subset for testing
        std::string test_subset = output_dir_ + "/test_subset.fna";
        if (!create_test_subset(test_subset, 5)) {  // Extract first 5 genomes
            std::cerr << "  ERROR: Failed to create test subset" << std::endl;
            return;
        }
        
        std::cout << "  Created test subset with 5 genomes" << std::endl;
        
        // Configure database builder
        DatabaseBuildConfig config;
        config.k_value = 31;
        config.ell_value = 31;
        config.minimizer_capacity = 1000000;  // 1M minimizers for test
        config.sequence_batch_size = 3;       // Small batch to test batching
        config.enable_progress_reporting = true;
        config.enable_debug_mode = true;
        config.debug_output_dir = output_dir_;
        
        // Create builder
        GPUKrakenDatabaseBuilder builder(output_dir_ + "/test_db", config);
        
        std::cout << "\n  Running database build on subset..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        
        bool success = builder.build_database_from_concatenated_fna(
            test_subset,
            "",  // No taxonomy files for this test
            "",
            ""
        );
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end - start).count();
        
        if (success) {
            std::cout << "  ✓ Database build completed in " << std::fixed << std::setprecision(2) 
                      << duration << " seconds" << std::endl;
            
            // Get statistics
            const auto& stats = builder.get_build_statistics();
            std::cout << "  ✓ Sequences processed: " << stats.total_sequences << std::endl;
            std::cout << "  ✓ Minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
            std::cout << "  ✓ Unique minimizers: " << stats.unique_minimizers << std::endl;
            std::cout << "  ✓ Processing rate: " << (stats.total_bases / duration / 1e6) << " Mbp/s" << std::endl;
        } else {
            std::cerr << "  ✗ Database build failed!" << std::endl;
        }
    }
    
    void test_memory_usage() {
        std::cout << "\nTest 4: Analyzing memory requirements..." << std::endl;
        
        // Get file size
        auto file_size = std::filesystem::file_size(fna_file_path_);
        
        // Estimate based on file analysis
        size_t estimated_sequences = 20000;  // Rough estimate
        size_t avg_sequence_length = file_size / estimated_sequences;
        
        std::cout << "  File size: " << (file_size / 1024 / 1024) << " MB" << std::endl;
        std::cout << "  Estimated sequences: " << estimated_sequences << std::endl;
        std::cout << "  Average sequence length: " << avg_sequence_length << " bp" << std::endl;
        
        // Memory requirements
        size_t sequence_memory = file_size;
        size_t genome_info_memory = estimated_sequences * sizeof(GPUGenomeInfo);
        size_t minimizer_memory = estimated_sequences * 1000 * sizeof(GPUMinimizerHit); // ~1000 minimizers per genome
        size_t total_memory = sequence_memory + genome_info_memory + minimizer_memory;
        
        std::cout << "\n  Memory requirements:" << std::endl;
        std::cout << "    Sequence data: " << (sequence_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "    Genome info: " << (genome_info_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "    Minimizer hits: " << (minimizer_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "    Total estimated: " << (total_memory / 1024 / 1024) << " MB" << std::endl;
        
        // Get GPU memory info
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "\n  GPU memory available: " << (free_mem / 1024 / 1024) << " MB" << std::endl;
        
        if (total_memory > free_mem * 0.8) {
            std::cout << "  ⚠ WARNING: May need to process in multiple passes" << std::endl;
        } else {
            std::cout << "  ✓ Sufficient GPU memory for single-pass processing" << std::endl;
        }
    }
    
    void test_minimizer_statistics() {
        std::cout << "\nTest 5: Analyzing minimizer patterns..." << std::endl;
        
        // This would analyze the output from test 3
        std::string hash_table_file = output_dir_ + "/test_db/hash_table.k2d";
        if (std::filesystem::exists(hash_table_file)) {
            auto file_size = std::filesystem::file_size(hash_table_file);
            std::cout << "  ✓ Hash table created: " << (file_size / 1024) << " KB" << std::endl;
        }
        
        // Check for debug output
        std::string feature_file = output_dir_ + "/minimizer_features.tsv";
        if (std::filesystem::exists(feature_file)) {
            std::cout << "  ✓ Feature statistics exported" << std::endl;
            
            // Read first few lines
            std::ifstream features(feature_file);
            std::string line;
            int count = 0;
            while (std::getline(features, line) && count++ < 5) {
                if (count == 1) {
                    std::cout << "  Feature columns: " << line << std::endl;
                }
            }
            features.close();
        }
    }
    
    bool create_test_subset(const std::string& output_file, int num_genomes) {
        std::ifstream in(fna_file_path_);
        std::ofstream out(output_file);
        
        if (!in.is_open() || !out.is_open()) {
            return false;
        }
        
        std::string line;
        int genomes_written = 0;
        bool in_sequence = false;
        
        while (std::getline(in, line) && genomes_written < num_genomes) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (in_sequence) {
                    genomes_written++;
                }
                in_sequence = true;
            }
            
            if (genomes_written < num_genomes) {
                out << line << "\n";
            }
        }
        
        in.close();
        out.close();
        return true;
    }
};

// Main test function
int main(int argc, char** argv) {
    std::string fna_file = "data/type_strain_reference_genomes/library.fna";
    
    if (argc > 1) {
        fna_file = argv[1];
    }
    
    std::cout << "Testing FNA file: " << fna_file << std::endl;
    
    try {
        FNAMinimizerTest tester(fna_file);
        tester.run_all_tests();
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}