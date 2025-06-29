// test_database_serializer.cu
// Comprehensive test suite for database serializer module

#include <iostream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <random>
#include <chrono>
#include <cassert>
#include <cstring>

#include "output/database_serializer.h"
#include "gpu_kraken_types.h"

// Test helper functions
namespace TestHelpers {
    
    // Generate random hash value
    uint64_t generate_random_hash() {
        static std::random_device rd;
        static std::mt19937_64 gen(rd());
        static std::uniform_int_distribution<uint64_t> dist;
        return dist(gen);
    }
    
    // Generate test LCA candidates
    std::vector<LCACandidate> generate_test_candidates(size_t count) {
        std::vector<LCACandidate> candidates;
        candidates.reserve(count);
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> taxon_dist(1, 100000);
        std::uniform_int_distribution<uint32_t> count_dist(1, 100);
        
        for (size_t i = 0; i < count; i++) {
            LCACandidate candidate;
            candidate.minimizer_hash = generate_random_hash();
            candidate.lca_taxon = taxon_dist(gen);
            candidate.genome_count = count_dist(gen);
            candidate.uniqueness_score = 1.0f / candidate.genome_count;
            candidates.push_back(candidate);
        }
        
        return candidates;
    }
    
    // Generate test taxonomy
    std::pair<std::unordered_map<uint32_t, std::string>, 
              std::unordered_map<uint32_t, uint32_t>> 
    generate_test_taxonomy(size_t species_count) {
        std::unordered_map<uint32_t, std::string> taxon_names;
        std::unordered_map<uint32_t, uint32_t> taxon_parents;
        
        // Root taxon
        taxon_names[1] = "root";
        taxon_parents[1] = 0;
        
        // Generate species
        for (size_t i = 0; i < species_count; i++) {
            uint32_t taxon_id = 100 + i;
            taxon_names[taxon_id] = "Species_" + std::to_string(i);
            taxon_parents[taxon_id] = 1; // All species under root for simplicity
        }
        
        return {taxon_names, taxon_parents};
    }
    
    // Generate test phylogenetic candidates
    std::vector<PhylogeneticLCACandidate> generate_phylo_candidates(size_t count) {
        std::vector<PhylogeneticLCACandidate> candidates;
        candidates.reserve(count);
        
        auto basic_candidates = generate_test_candidates(count);
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> species_count_dist(1, 5);
        std::uniform_int_distribution<uint32_t> species_dist(100, 199);
        std::uniform_int_distribution<uint16_t> genome_count_dist(1, 10);
        
        for (const auto& basic : basic_candidates) {
            PhylogeneticLCACandidate phylo(basic);
            
            // Add contributing species
            size_t num_species = species_count_dist(gen);
            for (size_t i = 0; i < num_species; i++) {
                phylo.contributing_species.push_back(species_dist(gen));
                phylo.genome_counts_per_species.push_back(genome_count_dist(gen));
            }
            
            phylo.phylogenetic_spread = num_species * 10;
            phylo.max_phylogenetic_distance = num_species * 5;
            
            candidates.push_back(phylo);
        }
        
        return candidates;
    }
    
    // Create test directory
    std::string create_test_directory(const std::string& test_name) {
        std::string dir = "/tmp/kraken_test_" + test_name + "_" + 
                         std::to_string(std::chrono::system_clock::now().time_since_epoch().count());
        std::filesystem::create_directories(dir);
        return dir;
    }
    
    // Clean up test directory
    void cleanup_test_directory(const std::string& dir) {
        try {
            std::filesystem::remove_all(dir);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to clean up test directory: " << e.what() << std::endl;
        }
    }
}

// Test standard database serialization
void test_standard_database_serialization() {
    std::cout << "\n=== Testing Standard Database Serialization ===" << std::endl;
    
    // Create test directory
    std::string test_dir = TestHelpers::create_test_directory("standard");
    
    try {
        // Create serializer
        StandardDatabaseSerializer serializer(test_dir);
        
        // Generate test data
        auto candidates = TestHelpers::generate_test_candidates(1000);
        auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(100);
        
        // Create test stats
        GPUBuildStats stats;
        stats.total_sequences = 500;
        stats.total_bases = 5000000;
        stats.unique_minimizers = candidates.size();
        stats.total_kmers_processed = 10000000;
        stats.valid_minimizers_extracted = 1000000;
        
        // Test serialization
        bool success = serializer.save_standard_database(candidates, taxon_names, taxon_parents, stats);
        assert(success);
        std::cout << "✓ Standard database saved successfully" << std::endl;
        
        // Verify files exist
        assert(std::filesystem::exists(test_dir + "/hash_table.k2d"));
        assert(std::filesystem::exists(test_dir + "/taxonomy.tsv"));
        assert(std::filesystem::exists(test_dir + "/config.txt"));
        assert(std::filesystem::exists(test_dir + "/build_stats.txt"));
        std::cout << "✓ All expected files created" << std::endl;
        
        // Test validation
        assert(serializer.validate_database_consistency());
        std::cout << "✓ Database consistency validated" << std::endl;
        
        // Check file sizes
        size_t hash_size = std::filesystem::file_size(test_dir + "/hash_table.k2d");
        std::cout << "  Hash table size: " << (hash_size / 1024) << " KB" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        assert(false);
    }
    
    // Cleanup
    TestHelpers::cleanup_test_directory(test_dir);
    std::cout << "✓ Standard database serialization test passed" << std::endl;
}

// Test enhanced database serialization
void test_enhanced_database_serialization() {
    std::cout << "\n=== Testing Enhanced Database Serialization ===" << std::endl;
    
    std::string test_dir = TestHelpers::create_test_directory("enhanced");
    
    try {
        // Create serializer
        EnhancedDatabaseSerializer serializer(test_dir);
        
        // Generate test data
        auto phylo_candidates = TestHelpers::generate_phylo_candidates(500);
        auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(100);
        
        // Create contributing taxa arrays
        ContributingTaxaArrays taxa_arrays;
        for (const auto& candidate : phylo_candidates) {
            for (size_t i = 0; i < candidate.contributing_species.size(); i++) {
                taxa_arrays.add_entry(
                    candidate.contributing_species[i],
                    i * 10,  // phylogenetic distance
                    candidate.genome_counts_per_species[i]
                );
            }
        }
        
        // Create enhanced stats
        EnhancedBuildStats stats;
        stats.total_sequences = 250;
        stats.total_bases = 2500000;
        stats.unique_minimizers = phylo_candidates.size();
        stats.species_represented = 100;
        stats.minimizers_with_phylo_data = phylo_candidates.size();
        stats.contributing_taxa_array_size = taxa_arrays.total_entries();
        
        // Test serialization
        bool success = serializer.save_enhanced_database(
            phylo_candidates, taxa_arrays, taxon_names, stats
        );
        assert(success);
        std::cout << "✓ Enhanced database saved successfully" << std::endl;
        
        // Verify enhanced files exist
        assert(std::filesystem::exists(test_dir + "/enhanced_hash_table.k2d"));
        assert(std::filesystem::exists(test_dir + "/contributing_taxa.k2d"));
        assert(std::filesystem::exists(test_dir + "/enhanced_config.txt"));
        assert(std::filesystem::exists(test_dir + "/species_mapping.tsv"));
        assert(std::filesystem::exists(test_dir + "/phylogenetic_summary.txt"));
        std::cout << "✓ All enhanced files created" << std::endl;
        
        // Verify backward compatibility layer
        assert(std::filesystem::exists(test_dir + "/hash_table.k2d"));
        std::cout << "✓ Backward compatibility layer created" << std::endl;
        
        // Test validation
        assert(serializer.validate_enhanced_format());
        std::cout << "✓ Enhanced format validated" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        assert(false);
    }
    
    TestHelpers::cleanup_test_directory(test_dir);
    std::cout << "✓ Enhanced database serialization test passed" << std::endl;
}

// Test database validation utilities
void test_database_validation() {
    std::cout << "\n=== Testing Database Validation Utilities ===" << std::endl;
    
    std::string test_dir = TestHelpers::create_test_directory("validation");
    
    try {
        // Create a standard database first
        StandardDatabaseSerializer serializer(test_dir);
        auto candidates = TestHelpers::generate_test_candidates(100);
        auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(10);
        GPUBuildStats stats;
        
        serializer.save_standard_database(candidates, taxon_names, taxon_parents, stats);
        
        // Test validation functions
        assert(DatabaseValidation::validate_hash_table_integrity(test_dir + "/hash_table.k2d"));
        std::cout << "✓ Hash table integrity validated" << std::endl;
        
        assert(DatabaseValidation::validate_taxonomy_consistency(test_dir + "/taxonomy.tsv"));
        std::cout << "✓ Taxonomy consistency validated" << std::endl;
        
        assert(DatabaseValidation::check_database_completeness(test_dir));
        std::cout << "✓ Database completeness validated" << std::endl;
        
        // Test with missing file
        std::filesystem::remove(test_dir + "/config.txt");
        assert(!DatabaseValidation::check_database_completeness(test_dir));
        std::cout << "✓ Missing file detection works" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        assert(false);
    }
    
    TestHelpers::cleanup_test_directory(test_dir);
    std::cout << "✓ Database validation test passed" << std::endl;
}

// Test database optimization utilities
void test_database_optimization() {
    std::cout << "\n=== Testing Database Optimization Utilities ===" << std::endl;
    
    std::string test_dir = TestHelpers::create_test_directory("optimization");
    
    try {
        // Create a database
        StandardDatabaseSerializer serializer(test_dir);
        auto candidates = TestHelpers::generate_test_candidates(1000);
        auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(50);
        GPUBuildStats stats;
        
        serializer.save_standard_database(candidates, taxon_names, taxon_parents, stats);
        
        // Test size calculation
        size_t db_size = DatabaseOptimization::calculate_database_size(test_dir);
        assert(db_size > 0);
        std::cout << "✓ Database size calculated: " << (db_size / 1024) << " KB" << std::endl;
        
        // Test memory estimation
        double mem_usage = DatabaseOptimization::estimate_memory_usage(test_dir);
        assert(mem_usage > db_size);
        std::cout << "✓ Memory usage estimated: " << (mem_usage / 1024 / 1024) << " MB" << std::endl;
        
        // Test statistics printing
        std::cout << "\nDatabase Statistics:" << std::endl;
        DatabaseOptimization::print_database_statistics(test_dir);
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        assert(false);
    }
    
    TestHelpers::cleanup_test_directory(test_dir);
    std::cout << "✓ Database optimization test passed" << std::endl;
}

// Test format detection
void test_format_detection() {
    std::cout << "\n=== Testing Database Format Detection ===" << std::endl;
    
    std::string standard_dir = TestHelpers::create_test_directory("format_standard");
    std::string enhanced_dir = TestHelpers::create_test_directory("format_enhanced");
    
    try {
        // Create standard database
        {
            StandardDatabaseSerializer serializer(standard_dir);
            auto candidates = TestHelpers::generate_test_candidates(10);
            auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(5);
            GPUBuildStats stats;
            serializer.save_standard_database(candidates, taxon_names, taxon_parents, stats);
        }
        
        // Create enhanced database
        {
            EnhancedDatabaseSerializer serializer(enhanced_dir);
            auto phylo_candidates = TestHelpers::generate_phylo_candidates(10);
            auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(5);
            ContributingTaxaArrays taxa_arrays;
            EnhancedBuildStats stats;
            serializer.save_enhanced_database(phylo_candidates, taxa_arrays, taxon_names, stats);
        }
        
        // Test format detection
        DatabaseFormatConverter converter("", "");
        
        assert(converter.detect_database_format(standard_dir) == DatabaseFormat::STANDARD_KRAKEN2);
        std::cout << "✓ Standard format detected correctly" << std::endl;
        
        assert(converter.detect_database_format(enhanced_dir) == DatabaseFormat::ENHANCED_PHYLO);
        std::cout << "✓ Enhanced format detected correctly" << std::endl;
        
        assert(converter.is_valid_database_directory(standard_dir));
        assert(converter.is_valid_database_directory(enhanced_dir));
        assert(!converter.is_valid_database_directory("/nonexistent/directory"));
        std::cout << "✓ Directory validation works correctly" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        assert(false);
    }
    
    TestHelpers::cleanup_test_directory(standard_dir);
    TestHelpers::cleanup_test_directory(enhanced_dir);
    std::cout << "✓ Format detection test passed" << std::endl;
}

// Performance benchmarking test
void test_serialization_performance() {
    std::cout << "\n=== Testing Serialization Performance ===" << std::endl;
    
    std::string test_dir = TestHelpers::create_test_directory("performance");
    
    try {
        // Test with increasing sizes
        std::vector<size_t> test_sizes = {1000, 10000, 50000};
        
        for (size_t size : test_sizes) {
            std::cout << "\nTesting with " << size << " minimizers:" << std::endl;
            
            // Generate data
            auto start = std::chrono::high_resolution_clock::now();
            auto candidates = TestHelpers::generate_test_candidates(size);
            auto [taxon_names, taxon_parents] = TestHelpers::generate_test_taxonomy(size / 10);
            auto gen_end = std::chrono::high_resolution_clock::now();
            
            double gen_time = std::chrono::duration<double>(gen_end - start).count();
            std::cout << "  Data generation: " << gen_time << " seconds" << std::endl;
            
            // Time serialization
            StandardDatabaseSerializer serializer(test_dir);
            GPUBuildStats stats;
            stats.unique_minimizers = size;
            
            start = std::chrono::high_resolution_clock::now();
            bool success = serializer.save_standard_database(candidates, taxon_names, taxon_parents, stats);
            auto end = std::chrono::high_resolution_clock::now();
            
            assert(success);
            double save_time = std::chrono::duration<double>(end - start).count();
            std::cout << "  Serialization: " << save_time << " seconds" << std::endl;
            std::cout << "  Rate: " << (size / save_time) << " minimizers/second" << std::endl;
            
            // Clean up for next iteration
            std::filesystem::remove_all(test_dir);
            std::filesystem::create_directories(test_dir);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Performance test failed: " << e.what() << std::endl;
        assert(false);
    }
    
    TestHelpers::cleanup_test_directory(test_dir);
    std::cout << "\n✓ Performance test completed" << std::endl;
}

// Main test runner
int main(int argc, char* argv[]) {
    std::cout << "=== Database Serializer Test Suite ===" << std::endl;
    std::cout << "Running comprehensive tests..." << std::endl;
    
    try {
        // Run all tests
        test_standard_database_serialization();
        test_enhanced_database_serialization();
        test_database_validation();
        test_database_optimization();
        test_format_detection();
        test_serialization_performance();
        
        std::cout << "\n=== ALL TESTS PASSED ===" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n=== TEST SUITE FAILED ===" << std::endl;
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}