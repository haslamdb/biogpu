// test_cooccurrence_real_data.cu
// Test co-occurrence scoring with real genome data

#include <iostream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <cuda_runtime.h>

#include "core/gpu_database_builder_core.h"
#include "processing/genome_file_processor.h"
#include "taxonomy/taxonomy_processor.h"
#include "memory/gpu_memory_manager.h"
#include "gpu_kraken_types.h"

void print_cooccurrence_analysis(const std::vector<GPUMinimizerHit>& hits) {
    std::cout << "\n=== CO-OCCURRENCE ANALYSIS ===" << std::endl;
    
    // Count minimizers by co-occurrence category
    std::vector<size_t> category_counts(8, 0);
    size_t total_with_scores = 0;
    
    for (const auto& hit : hits) {
        uint8_t category = MinimizerFlags::get_cooccurrence_score(hit.feature_flags);
        if (category > 0 || MinimizerFlags::get_cooccurrence_score(hit.feature_flags) == 0) {
            category_counts[category]++;
            total_with_scores++;
        }
    }
    
    std::cout << "Total minimizers with co-occurrence scores: " << total_with_scores << std::endl;
    std::cout << "\nDistribution by category:" << std::endl;
    
    const char* category_names[] = {
        "Very low (< 0.1)",
        "Low (0.1-0.2)",
        "Low-medium (0.2-0.35)",
        "Medium (0.35-0.5)",
        "Medium-high (0.5-0.65)",
        "High (0.65-0.8)",
        "Very high (0.8-0.95)",
        "Extremely high (≥ 0.95)"
    };
    
    for (int i = 0; i < 8; i++) {
        double percent = (total_with_scores > 0) ? 
            (100.0 * category_counts[i] / total_with_scores) : 0.0;
        std::cout << "  Category " << i << " - " << category_names[i] << ": " 
                  << category_counts[i] << " (" << std::fixed << std::setprecision(1) 
                  << percent << "%)" << std::endl;
    }
    
    // Find some examples of high co-occurrence minimizers
    std::cout << "\nExamples of high co-occurrence minimizers (category >= 5):" << std::endl;
    int examples_shown = 0;
    for (const auto& hit : hits) {
        uint8_t category = MinimizerFlags::get_cooccurrence_score(hit.feature_flags);
        if (category >= 5 && examples_shown < 5) {
            std::cout << "  Hash: 0x" << std::hex << hit.minimizer_hash << std::dec
                      << ", Genome: " << hit.genome_id 
                      << ", Position: " << hit.position
                      << ", Category: " << (int)category << std::endl;
            examples_shown++;
        }
    }
}

int main(int argc, char** argv) {
    std::cout << "=== CO-OCCURRENCE SCORING TEST WITH REAL DATA ===" << std::endl;
    
    // Configuration
    std::string fna_file = "/home/david/Documents/Code/biogpu/data/test_50_genomes.fna";
    std::string output_dir = "./test_cooccurrence_output";
    std::string nodes_file = "/home/david/Documents/Code/biogpu/data/nodes.dmp";
    std::string names_file = "/home/david/Documents/Code/biogpu/data/names.dmp";
    
    // Check if files exist
    if (!std::filesystem::exists(fna_file)) {
        std::cerr << "Error: FNA file not found: " << fna_file << std::endl;
        return 1;
    }
    
    // Create output directory
    std::filesystem::create_directories(output_dir);
    
    // Step 1: Configure database builder with co-occurrence enabled
    std::cout << "\n1. Configuring database builder..." << std::endl;
    
    DatabaseBuildConfig config;
    config.k_value = 31;
    config.ell_value = 31;
    config.enable_cooccurrence_scoring = true;
    config.cooccurrence_window_size = 10;
    config.enable_debug_mode = true;
    config.enable_uniqueness_scoring = true;  // Also enable uniqueness for comparison
    
    std::cout << "  Co-occurrence scoring: ENABLED" << std::endl;
    std::cout << "  Window size: " << config.cooccurrence_window_size << std::endl;
    std::cout << "  Debug mode: ENABLED" << std::endl;
    
    // Step 2: Initialize database builder
    std::cout << "\n2. Initializing GPU Kraken database builder..." << std::endl;
    
    GPUKrakenDatabaseBuilder builder(output_dir, config);
    
    // Step 3: Build database from the concatenated FNA file
    std::cout << "\n3. Building database from concatenated FNA file..." << std::endl;
    std::cout << "  Using file: " << fna_file << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Use the concatenated FNA method with taxonomy files if available
    bool build_success = false;
    if (std::filesystem::exists(nodes_file) && std::filesystem::exists(names_file)) {
        std::cout << "  Taxonomy files found, building with taxonomy support" << std::endl;
        build_success = builder.build_database_from_concatenated_fna(
            fna_file, 
            nodes_file, 
            names_file
        );
    } else {
        std::cout << "  Building without taxonomy support" << std::endl;
        build_success = builder.build_database_from_concatenated_fna(fna_file);
    }
    
    if (!build_success) {
        std::cerr << "Database build failed" << std::endl;
        return 1;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "✓ Database built successfully in " << duration.count() << " seconds" << std::endl;
    
    // Step 4: Get and display build statistics
    std::cout << "\n4. Build Statistics:" << std::endl;
    
    auto stats = builder.get_build_statistics();
    stats.print_summary();
    
    // Get enhanced stats for co-occurrence data
    auto enhanced_stats = builder.get_enhanced_statistics();
    
    // Print uniqueness stats if available
    if (config.enable_uniqueness_scoring) {
        enhanced_stats.print_uniqueness_stats();
    }
    
    // Print co-occurrence stats
    enhanced_stats.print_cooccurrence_stats();
    
    // Step 5: Verify co-occurrence scores were computed
    std::cout << "\n5. Verification:" << std::endl;
    
    if (enhanced_stats.minimizers_with_cooccurrence_scores > 0) {
        std::cout << "✓ Co-occurrence scores successfully computed for " 
                  << enhanced_stats.minimizers_with_cooccurrence_scores << " minimizers" << std::endl;
        
        double high_cooc_percent = (enhanced_stats.minimizers_with_cooccurrence_scores > 0) ?
            (100.0 * enhanced_stats.high_cooccurrence_minimizers / 
             enhanced_stats.minimizers_with_cooccurrence_scores) : 0.0;
        
        double low_cooc_percent = (enhanced_stats.minimizers_with_cooccurrence_scores > 0) ?
            (100.0 * enhanced_stats.low_cooccurrence_minimizers / 
             enhanced_stats.minimizers_with_cooccurrence_scores) : 0.0;
        
        std::cout << "  High co-occurrence: " << std::fixed << std::setprecision(1) 
                  << high_cooc_percent << "%" << std::endl;
        std::cout << "  Low co-occurrence: " << low_cooc_percent << "%" << std::endl;
        
        // Database is automatically saved by the builder
        std::cout << "\n6. Database output:" << std::endl;
        std::cout << "✓ Database files saved to " << output_dir << std::endl;
        
        std::cout << "\n✅ CO-OCCURRENCE TEST PASSED" << std::endl;
        return 0;
    } else {
        std::cerr << "❌ ERROR: No co-occurrence scores were computed!" << std::endl;
        std::cerr << "Check that co-occurrence scoring is properly integrated." << std::endl;
        return 1;
    }
}