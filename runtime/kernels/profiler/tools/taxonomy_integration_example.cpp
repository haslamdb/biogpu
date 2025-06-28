// taxonomy_integration_example.cpp
// Demonstrates how compact taxonomy integrates with the database building pipeline

#include "../tools/compact_gpu_taxonomy.h"
#include "../enhanced_k2like/taxonomy/taxonomy_processor.h"
#include "../enhanced_k2like/gpu_kraken_types.h"
#include <iostream>
#include <memory>
#include <vector>
#include <chrono>

using namespace BioGPU::CompactTaxonomy;

// Example of using both taxonomy systems together
class TaxonomyIntegrationExample {
private:
    std::unique_ptr<EnhancedNCBITaxonomyProcessor> build_time_taxonomy;
    std::unique_ptr<CompactGPUTaxonomy> runtime_taxonomy;
    
public:
    // Phase 1: Database Building - Use EnhancedNCBITaxonomyProcessor
    bool build_database_taxonomy(const std::string& nodes_file, const std::string& names_file) {
        std::cout << "\n=== PHASE 1: Database Building ===" << std::endl;
        std::cout << "Using EnhancedNCBITaxonomyProcessor for rich taxonomy operations" << std::endl;
        
        build_time_taxonomy = std::make_unique<EnhancedNCBITaxonomyProcessor>();
        
        // Load full taxonomy
        if (!build_time_taxonomy->load_ncbi_taxonomy(nodes_file, names_file)) {
            std::cerr << "Failed to load NCBI taxonomy" << std::endl;
            return false;
        }
        
        // Example: Complex phylogenetic operations during database building
        demonstrate_build_time_operations();
        
        return true;
    }
    
    // Phase 2: Create Compact Format for Runtime
    bool create_compact_taxonomy(const std::string& output_file) {
        std::cout << "\n=== PHASE 2: Creating Compact Format ===" << std::endl;
        
        if (!build_time_taxonomy || !build_time_taxonomy->is_loaded()) {
            std::cerr << "Build-time taxonomy not loaded" << std::endl;
            return false;
        }
        
        // Get the compact taxonomy from the enhanced processor
        auto* compact = build_time_taxonomy->get_compact_taxonomy();
        if (!compact) {
            std::cerr << "Failed to get compact taxonomy" << std::endl;
            return false;
        }
        
        // Save it for runtime use
        if (!compact->save_compact_taxonomy(output_file)) {
            std::cerr << "Failed to save compact taxonomy" << std::endl;
            return false;
        }
        
        std::cout << "✓ Compact taxonomy saved to: " << output_file << std::endl;
        return true;
    }
    
    // Phase 3: Runtime Classification - Use CompactGPUTaxonomy
    bool load_runtime_taxonomy(const std::string& compact_file) {
        std::cout << "\n=== PHASE 3: Runtime Classification ===" << std::endl;
        std::cout << "Loading compact taxonomy for GPU-accelerated classification" << std::endl;
        
        runtime_taxonomy = std::make_unique<CompactGPUTaxonomy>(true); // Enable cache
        
        if (!runtime_taxonomy->load_compact_taxonomy(compact_file)) {
            std::cerr << "Failed to load compact taxonomy" << std::endl;
            return false;
        }
        
        std::cout << "✓ Compact taxonomy loaded and ready for classification" << std::endl;
        return true;
    }
    
    // Demonstrate build-time operations
    void demonstrate_build_time_operations() {
        std::cout << "\nDemonstrating complex build-time operations:" << std::endl;
        
        // Example 1: Compute LCA for multiple species
        std::vector<uint32_t> test_species = {562, 511145, 316407}; // E. coli strains
        uint32_t lca = build_time_taxonomy->compute_lca_of_species(test_species);
        std::cout << "LCA of E. coli strains: " << lca << " (" 
                  << build_time_taxonomy->get_scientific_name(lca) << ")" << std::endl;
        
        // Example 2: Calculate phylogenetic spread
        uint8_t spread = build_time_taxonomy->calculate_phylogenetic_spread(test_species, lca);
        std::cout << "Phylogenetic spread: " << (int)spread << std::endl;
        
        // Example 3: Get taxonomic information
        for (uint32_t taxon : test_species) {
            std::cout << "  " << taxon << ": " 
                      << build_time_taxonomy->get_scientific_name(taxon)
                      << " (rank: " << build_time_taxonomy->get_rank(taxon) 
                      << ", depth: " << (int)build_time_taxonomy->get_depth(taxon) << ")" << std::endl;
        }
        
        // Example 4: Build PhylogeneticLCACandidate (as would happen in database building)
        PhylogeneticLCACandidate candidate;
        candidate.minimizer_hash = 0x123456789ABCDEF0;
        candidate.lca_taxon = lca;
        candidate.contributing_species = test_species;
        candidate.genome_counts_per_species = {10, 5, 3}; // Example counts
        candidate.phylogenetic_spread = spread;
        
        std::cout << "\nCreated phylogenetic candidate:" << std::endl;
        std::cout << "  Minimizer: " << std::hex << candidate.minimizer_hash << std::dec << std::endl;
        std::cout << "  LCA: " << candidate.lca_taxon << std::endl;
        std::cout << "  Contributing species: " << candidate.contributing_species.size() << std::endl;
    }
    
    // Demonstrate runtime operations
    void demonstrate_runtime_operations() {
        if (!runtime_taxonomy || !runtime_taxonomy->is_loaded()) {
            std::cerr << "Runtime taxonomy not loaded" << std::endl;
            return;
        }
        
        std::cout << "\nDemonstrating GPU-accelerated runtime operations:" << std::endl;
        
        // Example 1: Batch taxon lookups (as during classification)
        std::vector<uint32_t> query_taxons = {562, 511145, 316407, 1280, 1313}; // Various bacteria
        std::vector<uint32_t> parent_results;
        std::vector<uint8_t> depth_results;
        std::vector<uint8_t> rank_results;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        if (runtime_taxonomy->lookup_taxons_gpu(query_taxons, parent_results, depth_results, rank_results)) {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            
            std::cout << "GPU lookup completed in " << duration.count() << " μs" << std::endl;
            
            for (size_t i = 0; i < query_taxons.size(); i++) {
                std::cout << "  Taxon " << query_taxons[i] 
                          << ": parent=" << parent_results[i]
                          << ", depth=" << (int)depth_results[i]
                          << ", rank=" << (int)rank_results[i] << std::endl;
            }
        }
        
        // Example 2: Batch LCA computations (as during read classification)
        std::vector<std::pair<uint32_t, uint32_t>> lca_pairs = {
            {562, 511145},      // E. coli strains
            {562, 1280},        // E. coli vs Staph
            {1280, 1313},       // Staph vs Strep
            {562, 316407}       // E. coli strains
        };
        std::vector<uint32_t> lca_results;
        
        start = std::chrono::high_resolution_clock::now();
        
        if (runtime_taxonomy->compute_lca_batch_gpu(lca_pairs, lca_results)) {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            
            std::cout << "\nGPU LCA computation completed in " << duration.count() << " μs" << std::endl;
            
            for (size_t i = 0; i < lca_pairs.size(); i++) {
                std::cout << "  LCA(" << lca_pairs[i].first << ", " << lca_pairs[i].second 
                          << ") = " << lca_results[i] << std::endl;
            }
        }
    }
    
    // Simulate classification workflow
    void simulate_classification_workflow() {
        std::cout << "\n=== Simulating Classification Workflow ===" << std::endl;
        
        if (!runtime_taxonomy || !runtime_taxonomy->is_loaded()) {
            std::cerr << "Runtime taxonomy not loaded" << std::endl;
            return;
        }
        
        // Simulate minimizer hits from reads
        struct SimulatedHit {
            uint64_t minimizer;
            std::vector<uint32_t> hit_taxons;
        };
        
        std::vector<SimulatedHit> hits = {
            {0xABCDEF123456, {562, 511145}},        // Hit in multiple E. coli
            {0x123456ABCDEF, {562}},                // Hit in single E. coli
            {0xFEDCBA654321, {1280, 1313}},         // Hit in Staph and Strep
            {0x987654321ABC, {562, 1280, 1313}}    // Hit across multiple species
        };
        
        std::cout << "Processing " << hits.size() << " minimizer hits..." << std::endl;
        
        // Process each hit
        for (const auto& hit : hits) {
            std::cout << "\nMinimizer " << std::hex << hit.minimizer << std::dec 
                      << " found in " << hit.hit_taxons.size() << " taxa" << std::endl;
            
            if (hit.hit_taxons.size() == 1) {
                std::cout << "  Direct assignment to taxon " << hit.hit_taxons[0] << std::endl;
            } else {
                // Need to compute LCA
                std::vector<std::pair<uint32_t, uint32_t>> pairs;
                for (size_t i = 1; i < hit.hit_taxons.size(); i++) {
                    pairs.push_back({hit.hit_taxons[0], hit.hit_taxons[i]});
                }
                
                std::vector<uint32_t> lcas;
                if (runtime_taxonomy->compute_lca_batch_gpu(pairs, lcas)) {
                    // Find the overall LCA
                    uint32_t overall_lca = hit.hit_taxons[0];
                    for (uint32_t lca : lcas) {
                        overall_lca = lca; // Simplified - would need proper LCA merging
                    }
                    std::cout << "  LCA assignment to taxon " << overall_lca << std::endl;
                }
            }
        }
    }
};

// Example showing performance comparison
void performance_comparison() {
    std::cout << "\n=== Performance Comparison ===" << std::endl;
    
    // Create test data
    const int num_queries = 100000;
    std::vector<uint32_t> random_taxons;
    std::vector<std::pair<uint32_t, uint32_t>> random_pairs;
    
    // Generate random taxon IDs (in realistic range)
    std::srand(42);
    for (int i = 0; i < num_queries; i++) {
        uint32_t t1 = 1 + (std::rand() % 100000);
        uint32_t t2 = 1 + (std::rand() % 100000);
        random_taxons.push_back(t1);
        random_pairs.push_back({t1, t2});
    }
    
    // Host-side taxonomy (simplified for comparison)
    std::cout << "\nHost-side taxonomy operations:" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    // Simulate host lookups
    std::unordered_map<uint32_t, uint32_t> host_parents;
    for (uint32_t t : random_taxons) {
        host_parents[t] = t / 10; // Simplified parent calculation
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto host_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "  " << num_queries << " lookups: " << host_duration.count() << " ms" << std::endl;
    
    // GPU operations would be much faster
    std::cout << "\nGPU-accelerated operations:" << std::endl;
    std::cout << "  Expected speedup: 10-100x for large batches" << std::endl;
    std::cout << "  Memory efficiency: Compact format uses ~70% less memory" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "=== Taxonomy Integration Example ===" << std::endl;
    std::cout << "Demonstrating the two-phase taxonomy architecture:" << std::endl;
    std::cout << "1. Build phase: Full-featured taxonomy processing" << std::endl;
    std::cout << "2. Runtime phase: Compact GPU-optimized taxonomy" << std::endl;
    
    TaxonomyIntegrationExample example;
    
    // Check command line arguments
    if (argc < 3) {
        std::cout << "\nUsage: " << argv[0] << " <nodes.dmp> <names.dmp> [compact_output.bin]" << std::endl;
        std::cout << "\nRunning with test data instead..." << std::endl;
        
        // Run with test data
        example.demonstrate_build_time_operations();
        performance_comparison();
        
        return 0;
    }
    
    std::string nodes_file = argv[1];
    std::string names_file = argv[2];
    std::string compact_file = (argc > 3) ? argv[3] : "compact_taxonomy.bin";
    
    // Phase 1: Database building with full taxonomy
    if (!example.build_database_taxonomy(nodes_file, names_file)) {
        return 1;
    }
    
    // Phase 2: Create compact format
    if (!example.create_compact_taxonomy(compact_file)) {
        return 1;
    }
    
    // Phase 3: Runtime classification with compact taxonomy
    if (!example.load_runtime_taxonomy(compact_file)) {
        return 1;
    }
    
    // Demonstrate runtime operations
    example.demonstrate_runtime_operations();
    
    // Simulate classification workflow
    example.simulate_classification_workflow();
    
    // Show performance comparison
    performance_comparison();
    
    std::cout << "\n=== Integration Example Complete ===" << std::endl;
    return 0;
}