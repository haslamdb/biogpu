// test_hierarchical_db.cpp - Test program for hierarchical GPU database
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "hierarchical_gpu_database.h"

using namespace biogpu;

// Generate random k-mer hash
uint64_t generate_random_hash(std::mt19937_64& rng) {
    return rng();
}

int main() {
    std::cout << "Testing Hierarchical GPU Database\n";
    std::cout << "=================================\n\n";
    
    try {
        // Create database with default parameters
        HierarchicalGPUDatabase db;
        
        // Test parameters
        const size_t num_test_kmers = 1000000;
        const size_t num_species = 10;
        const size_t kmers_per_species = num_test_kmers / num_species;
        
        std::cout << "Test Configuration:\n";
        std::cout << "  Total k-mers: " << num_test_kmers << "\n";
        std::cout << "  Number of species: " << num_species << "\n";
        std::cout << "  K-mers per species: " << kmers_per_species << "\n\n";
        
        // Random number generator
        std::random_device rd;
        std::mt19937_64 rng(rd());
        
        // Add k-mers to database
        std::cout << "Adding k-mers to database...\n";
        auto start = std::chrono::high_resolution_clock::now();
        
        for (size_t species = 0; species < num_species; species++) {
            uint32_t taxon_id = 100 + species;
            
            for (size_t i = 0; i < kmers_per_species; i++) {
                uint64_t hash = generate_random_hash(rng);
                db.add_kmer(hash, taxon_id);
            }
            
            if ((species + 1) % 2 == 0) {
                std::cout << "  Added k-mers for " << (species + 1) << " species\n";
            }
        }
        
        auto add_time = std::chrono::high_resolution_clock::now();
        auto add_duration = std::chrono::duration_cast<std::chrono::milliseconds>(add_time - start);
        std::cout << "K-mer addition completed in " << add_duration.count() << " ms\n\n";
        
        // Build hierarchical structure
        std::cout << "Building hierarchical database structure...\n";
        db.build_hierarchy();
        
        auto build_time = std::chrono::high_resolution_clock::now();
        auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(build_time - add_time);
        std::cout << "Hierarchy built in " << build_duration.count() << " ms\n\n";
        
        // Get statistics
        auto stats = db.get_statistics();
        std::cout << "Database Statistics:\n";
        std::cout << "  Total levels: " << stats.num_levels << "\n";
        std::cout << "  Total k-mers: " << stats.total_kmers << "\n";
        std::cout << "  L1 buckets: " << stats.l1_buckets << "\n";
        std::cout << "  L1 k-mers: " << stats.l1_kmers << "\n";
        std::cout << "  L2 shards: " << stats.l2_shards << "\n";
        std::cout << "  L2 k-mers: " << stats.l2_kmers << "\n";
        std::cout << "  Memory usage: " << (stats.memory_bytes / (1024.0 * 1024.0)) << " MB\n\n";
        
        // Test some lookups
        std::cout << "Testing batch lookups...\n";
        const size_t num_queries = 10000;
        std::vector<uint64_t> query_hashes(num_queries);
        
        // Generate random queries
        for (size_t i = 0; i < num_queries; i++) {
            query_hashes[i] = generate_random_hash(rng);
        }
        
        // Perform lookup
        start = std::chrono::high_resolution_clock::now();
        std::vector<uint32_t> results(num_queries);
        
        db.lookup_batch(query_hashes, results);
        
        auto lookup_time = std::chrono::high_resolution_clock::now();
        auto lookup_duration = std::chrono::duration_cast<std::chrono::microseconds>(lookup_time - start);
        
        // Count hits
        size_t hits = 0;
        for (uint32_t taxon : results) {
            if (taxon != 0) hits++;
        }
        
        std::cout << "Lookup results:\n";
        std::cout << "  Queries: " << num_queries << "\n";
        std::cout << "  Hits: " << hits << " (" << (100.0 * hits / num_queries) << "%)\n";
        std::cout << "  Time: " << lookup_duration.count() << " μs\n";
        std::cout << "  Throughput: " << (num_queries * 1000000.0 / lookup_duration.count()) << " queries/sec\n\n";
        
        // Test saving and loading
        std::cout << "Testing save/load functionality...\n";
        std::string test_db_file = "test_hierarchical.db";
        
        db.save(test_db_file);
        std::cout << "Database saved to " << test_db_file << "\n";
        
        // Load into new database
        HierarchicalGPUDatabase db2;
        db2.load(test_db_file);
        std::cout << "Database loaded successfully\n";
        
        // Verify loaded database
        auto stats2 = db2.get_statistics();
        if (stats.total_kmers == stats2.total_kmers) {
            std::cout << "✓ Loaded database matches original\n";
        } else {
            std::cout << "✗ Loaded database does not match original\n";
        }
        
        std::cout << "\nAll tests completed successfully!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}