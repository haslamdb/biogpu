// test_compact_taxonomy.cpp
// Test program for compact GPU taxonomy

#include "compact_gpu_taxonomy.h"
#include <iostream>
#include <chrono>
#include <random>
#include <cassert>
#include <iomanip>

using namespace BioGPU::CompactTaxonomy;

// Test data - simplified taxonomy tree
void create_test_taxonomy(CompactGPUTaxonomy& taxonomy) {
    // Create a simple test taxonomy tree:
    // 1 (root)
    // ├── 2 (Bacteria)
    // │   ├── 10 (Proteobacteria)
    // │   │   ├── 20 (E. coli)
    // │   │   └── 21 (Salmonella)
    // │   └── 11 (Firmicutes)
    // │       ├── 22 (Bacillus)
    // │       └── 23 (Clostridium)
    // └── 3 (Archaea)
    //     └── 12 (Euryarchaeota)
    //         └── 24 (Methanobrevibacter)
    
    // Add taxa using internal methods (would normally load from files)
    taxonomy.add_test_taxon(1, 0, "root", "root", 0);
    taxonomy.add_test_taxon(2, 1, "Bacteria", "superkingdom", 1);
    taxonomy.add_test_taxon(3, 1, "Archaea", "superkingdom", 1);
    taxonomy.add_test_taxon(10, 2, "Proteobacteria", "phylum", 2);
    taxonomy.add_test_taxon(11, 2, "Firmicutes", "phylum", 2);
    taxonomy.add_test_taxon(12, 3, "Euryarchaeota", "phylum", 2);
    taxonomy.add_test_taxon(20, 10, "Escherichia coli", "species", 3);
    taxonomy.add_test_taxon(21, 10, "Salmonella enterica", "species", 3);
    taxonomy.add_test_taxon(22, 11, "Bacillus subtilis", "species", 3);
    taxonomy.add_test_taxon(23, 11, "Clostridium difficile", "species", 3);
    taxonomy.add_test_taxon(24, 12, "Methanobrevibacter smithii", "species", 3);
}

void test_basic_functionality() {
    std::cout << "\n=== Testing Basic Functionality ===" << std::endl;
    
    CompactGPUTaxonomy taxonomy(false); // No cache for basic test
    create_test_taxonomy(taxonomy);
    
    // Build GPU structures
    if (!taxonomy.build_test_taxonomy()) {
        std::cerr << "Failed to build test taxonomy" << std::endl;
        return;
    }
    
    // Test single lookups
    std::cout << "Testing single taxon lookups..." << std::endl;
    
    std::vector<uint32_t> test_taxons = {20, 21, 22, 23, 24};
    std::vector<uint32_t> expected_parents = {10, 10, 11, 11, 12};
    std::vector<uint8_t> expected_depths = {3, 3, 3, 3, 3};
    
    std::vector<uint32_t> parent_results;
    std::vector<uint8_t> depth_results;
    std::vector<uint8_t> rank_results;
    
    if (!taxonomy.lookup_taxons_gpu(test_taxons, parent_results, depth_results, rank_results)) {
        std::cerr << "GPU lookup failed" << std::endl;
        return;
    }
    
    // Verify results
    bool all_correct = true;
    for (size_t i = 0; i < test_taxons.size(); i++) {
        if (parent_results[i] != expected_parents[i]) {
            std::cerr << "Parent mismatch for taxon " << test_taxons[i] 
                      << ": got " << parent_results[i] 
                      << ", expected " << expected_parents[i] << std::endl;
            all_correct = false;
        }
        if (depth_results[i] != expected_depths[i]) {
            std::cerr << "Depth mismatch for taxon " << test_taxons[i] 
                      << ": got " << (int)depth_results[i] 
                      << ", expected " << (int)expected_depths[i] << std::endl;
            all_correct = false;
        }
    }
    
    if (all_correct) {
        std::cout << "✓ All single lookups passed" << std::endl;
    }
}

void test_lca_computation() {
    std::cout << "\n=== Testing LCA Computation ===" << std::endl;
    
    CompactGPUTaxonomy taxonomy(false);
    create_test_taxonomy(taxonomy);
    
    if (!taxonomy.build_test_taxonomy()) {
        std::cerr << "Failed to build test taxonomy" << std::endl;
        return;
    }
    
    // Test LCA pairs
    std::vector<std::pair<uint32_t, uint32_t>> test_pairs = {
        {20, 21},  // E. coli vs Salmonella -> Proteobacteria (10)
        {20, 22},  // E. coli vs Bacillus -> Bacteria (2)
        {20, 24},  // E. coli vs Methanobrevibacter -> root (1)
        {22, 23},  // Bacillus vs Clostridium -> Firmicutes (11)
        {20, 20}   // Same taxon -> itself (20)
    };
    
    std::vector<uint32_t> expected_lcas = {10, 2, 1, 11, 20};
    std::vector<uint32_t> lca_results;
    
    if (!taxonomy.compute_lca_batch_gpu(test_pairs, lca_results)) {
        std::cerr << "GPU LCA computation failed" << std::endl;
        return;
    }
    
    // Verify results
    bool all_correct = true;
    for (size_t i = 0; i < test_pairs.size(); i++) {
        if (lca_results[i] != expected_lcas[i]) {
            std::cerr << "LCA mismatch for pair (" << test_pairs[i].first 
                      << ", " << test_pairs[i].second << "): got " 
                      << lca_results[i] << ", expected " << expected_lcas[i] << std::endl;
            all_correct = false;
        }
    }
    
    if (all_correct) {
        std::cout << "✓ All LCA computations passed" << std::endl;
    }
}

void test_performance() {
    std::cout << "\n=== Testing Performance ===" << std::endl;
    
    CompactGPUTaxonomy taxonomy(true); // Enable cache
    
    // Generate larger test taxonomy
    std::cout << "Generating large test taxonomy..." << std::endl;
    taxonomy.add_test_taxon(1, 0, "root", "root", 0);
    
    // Add kingdoms
    for (int k = 0; k < 5; k++) {
        uint32_t kingdom_id = 2 + k;
        taxonomy.add_test_taxon(kingdom_id, 1, "Kingdom_" + std::to_string(k), "kingdom", 1);
        
        // Add phyla
        for (int p = 0; p < 10; p++) {
            uint32_t phylum_id = 100 + k * 10 + p;
            taxonomy.add_test_taxon(phylum_id, kingdom_id, "Phylum_" + std::to_string(phylum_id), "phylum", 2);
            
            // Add species
            for (int s = 0; s < 100; s++) {
                uint32_t species_id = 10000 + k * 1000 + p * 100 + s;
                taxonomy.add_test_taxon(species_id, phylum_id, "Species_" + std::to_string(species_id), "species", 3);
            }
        }
    }
    
    std::cout << "Built test taxonomy with " << taxonomy.get_test_size() << " taxa" << std::endl;
    
    if (!taxonomy.build_test_taxonomy()) {
        std::cerr << "Failed to build large test taxonomy" << std::endl;
        return;
    }
    
    // Generate random queries
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> species_dist(10000, 14999);
    
    const int num_queries = 10000;
    std::vector<uint32_t> query_taxons;
    std::vector<std::pair<uint32_t, uint32_t>> lca_pairs;
    
    for (int i = 0; i < num_queries; i++) {
        query_taxons.push_back(species_dist(gen));
        lca_pairs.push_back({species_dist(gen), species_dist(gen)});
    }
    
    // Benchmark lookups
    std::cout << "\nBenchmarking " << num_queries << " taxon lookups..." << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<uint32_t> parent_results;
    std::vector<uint8_t> depth_results;
    std::vector<uint8_t> rank_results;
    
    bool success = taxonomy.lookup_taxons_gpu(query_taxons, parent_results, depth_results, rank_results);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (success) {
        double throughput = (double)num_queries / duration.count() * 1000000;
        std::cout << "Lookup time: " << duration.count() << " μs" << std::endl;
        std::cout << "Throughput: " << std::fixed << std::setprecision(2) 
                  << throughput << " lookups/second" << std::endl;
    }
    
    // Benchmark LCA computation
    std::cout << "\nBenchmarking " << num_queries << " LCA computations..." << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    
    std::vector<uint32_t> lca_results;
    success = taxonomy.compute_lca_batch_gpu(lca_pairs, lca_results);
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (success) {
        double throughput = (double)num_queries / duration.count() * 1000000;
        std::cout << "LCA computation time: " << duration.count() << " μs" << std::endl;
        std::cout << "Throughput: " << std::fixed << std::setprecision(2) 
                  << throughput << " LCAs/second" << std::endl;
    }
}

void test_save_and_load() {
    std::cout << "\n=== Testing Save and Load ===" << std::endl;
    
    const std::string test_file = "test_compact_taxonomy.bin";
    
    // Create and save
    {
        CompactGPUTaxonomy taxonomy(false);
        create_test_taxonomy(taxonomy);
        
        if (!taxonomy.save_compact_taxonomy(test_file)) {
            std::cerr << "Failed to save taxonomy" << std::endl;
            return;
        }
        std::cout << "✓ Saved taxonomy to " << test_file << std::endl;
    }
    
    // Load and verify
    {
        CompactGPUTaxonomy taxonomy2(false);
        
        if (!taxonomy2.load_compact_taxonomy(test_file)) {
            std::cerr << "Failed to load taxonomy" << std::endl;
            return;
        }
        std::cout << "✓ Loaded taxonomy from " << test_file << std::endl;
        
        // Test a lookup to verify it works
        std::vector<uint32_t> test_taxons = {20};
        std::vector<uint32_t> parent_results;
        std::vector<uint8_t> depth_results;
        std::vector<uint8_t> rank_results;
        
        if (taxonomy2.lookup_taxons_gpu(test_taxons, parent_results, depth_results, rank_results)) {
            if (parent_results[0] == 10) {
                std::cout << "✓ Loaded taxonomy verified working" << std::endl;
            } else {
                std::cerr << "Loaded taxonomy gave wrong result" << std::endl;
            }
        }
    }
    
    // Cleanup
    std::remove(test_file.c_str());
}

int main() {
    std::cout << "=== Compact GPU Taxonomy Test Suite ===" << std::endl;
    
    // Check CUDA availability
    int device_count;
    cudaError_t error = cudaGetDeviceCount(&device_count);
    if (error != cudaSuccess || device_count == 0) {
        std::cerr << "No CUDA devices available" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
    
    // Run tests
    test_basic_functionality();
    test_lca_computation();
    test_performance();
    test_save_and_load();
    
    std::cout << "\n=== All tests completed ===" << std::endl;
    
    return 0;
}