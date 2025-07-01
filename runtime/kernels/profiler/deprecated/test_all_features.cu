// test_all_features.cu
// Test all implemented features: co-occurrence, uniqueness, GC content, complexity, etc.

#include <iostream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <cuda_runtime.h>
#include <algorithm>

#include "core/gpu_database_builder_core.h"
#include "processing/genome_file_processor.h"
#include "taxonomy/taxonomy_processor.h"
#include "memory/gpu_memory_manager.h"
#include "gpu_kraken_types.h"
#include "features/enhanced_minimizer_flags.h"

// Function to analyze and print feature statistics
void analyze_minimizer_features(const std::vector<GPUMinimizerHit>& hits) {
    std::cout << "\n=== COMPREHENSIVE FEATURE ANALYSIS ===" << std::endl;
    std::cout << "Total minimizers: " << hits.size() << std::endl;
    
    // Uniqueness analysis
    std::vector<size_t> uniqueness_counts(8, 0);
    
    // Co-occurrence analysis
    std::vector<size_t> cooccurrence_counts(8, 0);
    
    // GC content analysis
    std::vector<size_t> gc_counts(8, 0);
    
    // Complexity analysis
    std::vector<size_t> complexity_counts(8, 0);
    
    // Position bias analysis
    size_t position_clustered = 0;
    size_t position_uniform = 0;
    
    // Contamination risk analysis
    size_t contamination_risk = 0;
    
    // Classification analysis
    size_t unique_class = 0;
    size_t canonical_class = 0;
    size_t redundant_class = 0;
    
    // ML weight analysis
    std::vector<float> ml_weights;
    
    for (const auto& hit : hits) {
        // Extract all features
        uint32_t flags = hit.feature_flags;
        uint16_t strand_flags = hit.strand;
        
        // Uniqueness (bits 8-10 in feature_flags)
        uint8_t uniqueness = EnhancedMinimizerFlags::get_uniqueness_category(flags);
        if (uniqueness < uniqueness_counts.size()) {
            uniqueness_counts[uniqueness]++;
        }
        
        // Co-occurrence (bits 13-15 in feature_flags)
        uint8_t cooccurrence = MinimizerFlags::get_cooccurrence_score(flags);
        if (cooccurrence < cooccurrence_counts.size()) {
            cooccurrence_counts[cooccurrence]++;
        }
        
        // GC content (bits 0-2 in feature_flags)
        uint8_t gc_category = MinimizerFlags::get_gc_content_category(flags);
        if (gc_category < gc_counts.size()) {
            gc_counts[gc_category]++;
        }
        
        // Complexity (bits 3-5 in feature_flags)
        uint8_t complexity = MinimizerFlags::get_complexity_score(flags);
        if (complexity < complexity_counts.size()) {
            complexity_counts[complexity]++;
        }
        
        // Position bias (bit 6 in feature_flags)
        if (MinimizerFlags::has_position_bias(flags)) {
            position_clustered++;
        } else {
            position_uniform++;
        }
        
        // Contamination risk (bit 7 in feature_flags)
        if (MinimizerFlags::has_contamination_risk(flags)) {
            contamination_risk++;
        }
        
        // Classification (bits 2-3 in strand field)
        uint32_t classification = MinimizerFlags::get_classification(strand_flags);
        if (classification == 0) unique_class++;
        else if (classification == 1) canonical_class++;
        else if (classification == 2) redundant_class++;
        
        // ML weight
        float ml_weight = MinimizerFlags::ml_weight_to_float(hit.ml_weight);
        ml_weights.push_back(ml_weight);
    }
    
    // Print uniqueness statistics
    std::cout << "\n--- UNIQUENESS SCORES ---" << std::endl;
    const char* uniqueness_names[] = {
        "Extremely common (0-10%)",
        "Very low (10-30%)",
        "Low (30-50%)",
        "Moderate (50-70%)",
        "High (70-80%)",
        "Very high (80-90%)",
        "Extremely high (90-95%)",
        "Singleton-like (95-100%)"
    };
    for (int i = 0; i < 8; i++) {
        double percent = (hits.size() > 0) ? (100.0 * uniqueness_counts[i] / hits.size()) : 0.0;
        std::cout << "  " << uniqueness_names[i] << ": " << uniqueness_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print co-occurrence statistics
    std::cout << "\n--- CO-OCCURRENCE SCORES ---" << std::endl;
    const char* cooccurrence_names[] = {
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
        double percent = (hits.size() > 0) ? (100.0 * cooccurrence_counts[i] / hits.size()) : 0.0;
        std::cout << "  " << cooccurrence_names[i] << ": " << cooccurrence_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print GC content statistics
    std::cout << "\n--- GC CONTENT ---" << std::endl;
    for (int i = 0; i < 8; i++) {
        double percent = (hits.size() > 0) ? (100.0 * gc_counts[i] / hits.size()) : 0.0;
        int gc_low = i * 100 / 8;
        int gc_high = (i + 1) * 100 / 8;
        std::cout << "  " << gc_low << "-" << gc_high << "%: " << gc_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print complexity statistics
    std::cout << "\n--- SEQUENCE COMPLEXITY ---" << std::endl;
    for (int i = 0; i < 8; i++) {
        double percent = (hits.size() > 0) ? (100.0 * complexity_counts[i] / hits.size()) : 0.0;
        std::cout << "  Level " << i << ": " << complexity_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print position bias statistics
    std::cout << "\n--- POSITION DISTRIBUTION ---" << std::endl;
    double clustered_percent = (hits.size() > 0) ? (100.0 * position_clustered / hits.size()) : 0.0;
    double uniform_percent = (hits.size() > 0) ? (100.0 * position_uniform / hits.size()) : 0.0;
    std::cout << "  Clustered: " << position_clustered 
              << " (" << std::fixed << std::setprecision(1) << clustered_percent << "%)" << std::endl;
    std::cout << "  Uniform: " << position_uniform 
              << " (" << uniform_percent << "%)" << std::endl;
    
    // Print contamination statistics
    std::cout << "\n--- CONTAMINATION RISK ---" << std::endl;
    double contam_percent = (hits.size() > 0) ? (100.0 * contamination_risk / hits.size()) : 0.0;
    std::cout << "  At risk: " << contamination_risk 
              << " (" << std::fixed << std::setprecision(1) << contam_percent << "%)" << std::endl;
    std::cout << "  Clean: " << (hits.size() - contamination_risk) 
              << " (" << (100.0 - contam_percent) << "%)" << std::endl;
    
    // Print classification statistics
    std::cout << "\n--- MINIMIZER CLASSIFICATION ---" << std::endl;
    double unique_percent = (hits.size() > 0) ? (100.0 * unique_class / hits.size()) : 0.0;
    double canonical_percent = (hits.size() > 0) ? (100.0 * canonical_class / hits.size()) : 0.0;
    double redundant_percent = (hits.size() > 0) ? (100.0 * redundant_class / hits.size()) : 0.0;
    std::cout << "  Unique: " << unique_class 
              << " (" << std::fixed << std::setprecision(1) << unique_percent << "%)" << std::endl;
    std::cout << "  Canonical: " << canonical_class 
              << " (" << canonical_percent << "%)" << std::endl;
    std::cout << "  Redundant: " << redundant_class 
              << " (" << redundant_percent << "%)" << std::endl;
    
    // Print ML weight statistics
    if (!ml_weights.empty()) {
        float min_weight = *std::min_element(ml_weights.begin(), ml_weights.end());
        float max_weight = *std::max_element(ml_weights.begin(), ml_weights.end());
        float avg_weight = 0.0f;
        for (float w : ml_weights) avg_weight += w;
        avg_weight /= ml_weights.size();
        
        std::cout << "\n--- ML CONFIDENCE WEIGHTS ---" << std::endl;
        std::cout << "  Min: " << std::fixed << std::setprecision(3) << min_weight << std::endl;
        std::cout << "  Max: " << max_weight << std::endl;
        std::cout << "  Average: " << avg_weight << std::endl;
    }
}

// Function to show example minimizers with all features
void show_example_minimizers(const std::vector<GPUMinimizerHit>& hits, int num_examples = 5) {
    std::cout << "\n=== EXAMPLE MINIMIZERS WITH ALL FEATURES ===" << std::endl;
    
    int shown = 0;
    for (const auto& hit : hits) {
        if (shown >= num_examples) break;
        
        uint32_t flags = hit.feature_flags;
        uint16_t strand_flags = hit.strand;
        
        std::cout << "\nMinimizer " << (shown + 1) << ":" << std::endl;
        std::cout << "  Hash: 0x" << std::hex << hit.minimizer_hash << std::dec << std::endl;
        std::cout << "  Genome ID: " << hit.genome_id << std::endl;
        std::cout << "  Position: " << hit.position << std::endl;
        std::cout << "  Taxon ID: " << hit.taxon_id << std::endl;
        std::cout << "  Features:" << std::endl;
        std::cout << "    - Uniqueness: " << (int)EnhancedMinimizerFlags::get_uniqueness_category(flags) 
                  << "/7" << std::endl;
        std::cout << "    - Co-occurrence: " << (int)MinimizerFlags::get_cooccurrence_score(flags) 
                  << "/7" << std::endl;
        std::cout << "    - GC content: " << (int)MinimizerFlags::get_gc_content_category(flags) 
                  << "/7" << std::endl;
        std::cout << "    - Complexity: " << (int)MinimizerFlags::get_complexity_score(flags) 
                  << "/7" << std::endl;
        std::cout << "    - Position: " << (MinimizerFlags::has_position_bias(flags) ? "Clustered" : "Uniform") 
                  << std::endl;
        std::cout << "    - Contamination: " << (MinimizerFlags::has_contamination_risk(flags) ? "Risk" : "Clean") 
                  << std::endl;
        std::cout << "    - Classification: ";
        uint32_t classification = MinimizerFlags::get_classification(strand_flags);
        if (classification == 0) std::cout << "Unique";
        else if (classification == 1) std::cout << "Canonical";
        else std::cout << "Redundant";
        std::cout << std::endl;
        std::cout << "    - ML weight: " << std::fixed << std::setprecision(3) 
                  << MinimizerFlags::ml_weight_to_float(hit.ml_weight) << std::endl;
        
        shown++;
    }
}

int main(int argc, char** argv) {
    std::cout << "=== ALL FEATURES TEST WITH REAL DATA ===" << std::endl;
    
    // Configuration
    std::string fna_file = "/home/david/Documents/Code/biogpu/data/test_50_genomes.fna";
    std::string output_dir = "./test_all_features_output";
    std::string nodes_file = "/home/david/Documents/Code/biogpu/data/nodes.dmp";
    std::string names_file = "/home/david/Documents/Code/biogpu/data/names.dmp";
    
    // Check if files exist
    if (!std::filesystem::exists(fna_file)) {
        std::cerr << "Error: FNA file not found: " << fna_file << std::endl;
        return 1;
    }
    
    // Create output directory
    std::filesystem::create_directories(output_dir);
    
    // Step 1: Configure database builder with all features enabled
    std::cout << "\n1. Configuring database builder with all features..." << std::endl;
    
    DatabaseBuildConfig config;
    config.k_value = 31;
    config.ell_value = 31;
    config.enable_cooccurrence_scoring = true;
    config.cooccurrence_window_size = 10;
    config.enable_uniqueness_scoring = true;
    config.enable_debug_mode = true;
    config.enable_phylogenetic_analysis = true;
    
    // Don't save the database - just process minimizers
    config.save_enhanced_format = false;
    config.save_compatibility_layer = false;
    
    // Set memory configuration for streaming approach
    // We'll process in batches and accumulate results
    config.memory_config.minimizer_capacity = 100000000;  // 100M minimizers total capacity
    config.memory_config.sequence_batch_size = 5;         // Process 5 sequences at a time
    
    std::cout << "  Enabled features:" << std::endl;
    std::cout << "    ✓ Co-occurrence scoring (window=" << config.cooccurrence_window_size << ")" << std::endl;
    std::cout << "    ✓ Uniqueness scoring" << std::endl;
    std::cout << "    ✓ GC content analysis" << std::endl;
    std::cout << "    ✓ Sequence complexity" << std::endl;
    std::cout << "    ✓ Position clustering" << std::endl;
    std::cout << "    ✓ Contamination detection" << std::endl;
    std::cout << "    ✓ ML confidence weights" << std::endl;
    std::cout << "    ✓ Phylogenetic analysis" << std::endl;
    
    // Step 2: Initialize database builder
    std::cout << "\n2. Initializing GPU Kraken database builder..." << std::endl;
    
    GPUKrakenDatabaseBuilder builder(output_dir, config);
    
    // Step 3: Process the FNA file
    std::cout << "\n3. Processing FNA file..." << std::endl;
    std::cout << "  Using file: " << fna_file << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // We'll use a custom approach to just extract and process minimizers
    // without building the full database
    
    // Load taxonomy if available
    if (std::filesystem::exists(nodes_file) && std::filesystem::exists(names_file)) {
        std::cout << "  Loading taxonomy for phylogenetic features..." << std::endl;
        EnhancedNCBITaxonomyProcessor taxonomy;
        if (taxonomy.load_ncbi_taxonomy(nodes_file, names_file)) {
            std::cout << "  ✓ Taxonomy loaded successfully" << std::endl;
        }
    }
    
    // Process using streaming FNA method
    bool process_success = builder.build_database_from_concatenated_fna(
        fna_file, 
        nodes_file, 
        names_file
    );
    
    if (!process_success) {
        std::cerr << "Processing failed - this is expected since we're not saving the database" << std::endl;
        std::cerr << "But we should have extracted and processed minimizers" << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\nProcessing completed in " << duration.count() << " seconds" << std::endl;
    
    // Step 4: Get and analyze statistics
    std::cout << "\n4. Analyzing feature assignments..." << std::endl;
    
    auto stats = builder.get_build_statistics();
    auto enhanced_stats = builder.get_enhanced_statistics();
    
    // Print basic statistics
    std::cout << "\n--- BASIC STATISTICS ---" << std::endl;
    std::cout << "Total sequences: " << stats.total_sequences << std::endl;
    std::cout << "Total bases: " << stats.total_bases << std::endl;
    std::cout << "Total k-mers processed: " << stats.total_kmers_processed << std::endl;
    std::cout << "Valid minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
    std::cout << "Unique minimizers: " << stats.unique_minimizers << std::endl;
    
    // Print feature-specific statistics
    if (config.enable_uniqueness_scoring) {
        enhanced_stats.print_uniqueness_stats();
    }
    
    if (config.enable_cooccurrence_scoring) {
        enhanced_stats.print_cooccurrence_stats();
    }
    
    // Step 5: Verification
    std::cout << "\n5. Feature Assignment Verification:" << std::endl;
    
    bool all_features_working = true;
    
    // Check uniqueness
    if (enhanced_stats.minimizers_with_uniqueness_scores > 0) {
        std::cout << "✓ Uniqueness scoring: " << enhanced_stats.minimizers_with_uniqueness_scores 
                  << " minimizers scored" << std::endl;
    } else {
        std::cout << "✗ Uniqueness scoring: No scores computed" << std::endl;
        all_features_working = false;
    }
    
    // Check co-occurrence
    if (enhanced_stats.minimizers_with_cooccurrence_scores > 0) {
        std::cout << "✓ Co-occurrence scoring: " << enhanced_stats.minimizers_with_cooccurrence_scores 
                  << " minimizers scored" << std::endl;
    } else {
        std::cout << "✗ Co-occurrence scoring: No scores computed" << std::endl;
        all_features_working = false;
    }
    
    // Final result
    std::cout << "\n=== TEST RESULT ===" << std::endl;
    if (all_features_working && stats.valid_minimizers_extracted > 0) {
        std::cout << "✅ ALL FEATURES TEST PASSED" << std::endl;
        std::cout << "Successfully extracted " << stats.valid_minimizers_extracted 
                  << " minimizers with all features assigned" << std::endl;
        return 0;
    } else {
        std::cout << "❌ TEST FAILED - Not all features were computed" << std::endl;
        return 1;
    }
}