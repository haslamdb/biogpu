// uniqueness_integration_patch.cu
// Integration patch to add uniqueness score calculation to the main database builder
// This shows how to modify gpu_database_builder_core.cu to include uniqueness scoring

// Add this include to the existing gpu_database_builder_core.cu
#include "uniqueness_score_implementation.cu"

// ===========================
// Modified process_accumulated_sequences() method
// Add this after the existing minimizer extraction in gpu_database_builder_core.cu
// ===========================

bool GPUKrakenDatabaseBuilder::process_accumulated_sequences() {
    std::cout << "Processing " << accumulated_genome_info_.size() << " accumulated sequences..." << std::endl;
    
    // [... existing code for minimizer extraction ...]
    
    if (total_hits_extracted > 0) {
        // EXISTING CODE: Process minimizers with feature extraction
        if (feature_extractor_) {
            std::cout << "Running feature extraction on minimizers..." << std::endl;
            
            // Create taxon IDs vector from genome info
            std::vector<uint32_t> genome_taxon_ids;
            for (const auto& info : accumulated_genome_info_) {
                genome_taxon_ids.push_back(info.taxon_id);
            }
            
            // MODIFIED: Use integrated uniqueness calculation
            std::cout << "Computing uniqueness scores..." << std::endl;
            if (!integrate_uniqueness_with_feature_extractor(
                    feature_extractor_.get(),
                    d_minimizer_hits, 
                    total_hits_extracted, 
                    genome_taxon_ids)) {
                std::cerr << "Uniqueness integration failed" << std::endl;
                return false;
            }
            
            std::cout << "✓ Feature extraction with uniqueness completed" << std::endl;
        } else {
            // Fallback: compute uniqueness without feature extractor
            std::cout << "Computing uniqueness scores (standalone)..." << std::endl;
            
            std::vector<uint32_t> genome_taxon_ids;
            for (const auto& info : accumulated_genome_info_) {
                genome_taxon_ids.push_back(info.taxon_id);
            }
            
            if (!compute_and_encode_uniqueness_scores(
                    d_minimizer_hits,
                    total_hits_extracted,
                    genome_taxon_ids,
                    accumulated_genome_info_.size())) {
                std::cerr << "Standalone uniqueness computation failed" << std::endl;
                return false;
            }
        }
        
        // NEW: Apply uniqueness filtering to improve quality
        if (config_.enable_uniqueness_filtering) {
            std::cout << "Applying uniqueness filtering..." << std::endl;
            size_t original_count = total_hits_extracted;
            
            if (!apply_uniqueness_filtering(
                    d_minimizer_hits,
                    total_hits_extracted,
                    config_.uniqueness_threshold,
                    config_.filter_extremely_common)) {
                std::cerr << "Uniqueness filtering failed" << std::endl;
                // Non-fatal - continue without filtering
            } else {
                std::cout << "Uniqueness filtering: " << original_count << " → " 
                          << total_hits_extracted << " minimizers" << std::endl;
                
                // Update statistics
                enhanced_stats_.minimizers_filtered_by_uniqueness = original_count - total_hits_extracted;
            }
        }
        
        // Copy minimizer hits back to host for further processing
        std::vector<GPUMinimizerHit> minimizer_hits(total_hits_extracted);
        cudaMemcpy(minimizer_hits.data(), d_minimizer_hits, 
                   total_hits_extracted * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
        
        // NEW: Print uniqueness quality summary
        if (config_.enable_debug_mode) {
            UniquenessUtils::print_minimizer_quality_summary(minimizer_hits);
        }
        
        // [... rest of existing code ...]
    }
    
    return true;
}

// ===========================
// Enhanced Configuration Structure
// Add these fields to DatabaseBuildConfig in gpu_database_builder_core.h
// ===========================

/*
struct DatabaseBuildConfig {
    // ... existing fields ...
    
    // NEW: Uniqueness scoring configuration
    bool enable_uniqueness_scoring = true;        // Enable uniqueness calculation
    bool enable_uniqueness_filtering = false;     // Enable filtering by uniqueness
    float uniqueness_threshold = 0.3f;            // Minimum uniqueness score to keep
    bool filter_extremely_common = true;          // Filter minimizers with uniqueness < 0.1
    
    // NEW: Quality control settings
    bool enable_quality_analysis = true;          // Print quality summaries
    float false_positive_reduction_target = 0.5f; // Target FP reduction (0.0-1.0)
};
*/

// ===========================
// Enhanced Statistics Structure  
// Add these fields to EnhancedBuildStats in gpu_kraken_types.h
// ===========================

/*
struct EnhancedBuildStats {
    // ... existing fields ...
    
    // NEW: Uniqueness-related statistics
    size_t minimizers_with_uniqueness_scores = 0;
    size_t minimizers_filtered_by_uniqueness = 0;
    size_t unique_minimizers_count = 0;           // ≥90% uniqueness
    size_t rare_minimizers_count = 0;             // ≤3 occurrences
    size_t reliable_minimizers_count = 0;         // Suitable for classification
    double average_uniqueness_score = 0.0;
    double estimated_false_positive_reduction = 0.0;
    
    void print_uniqueness_stats() const {
        std::cout << "\n=== UNIQUENESS STATISTICS ===" << std::endl;
        std::cout << "Minimizers with uniqueness scores: " << minimizers_with_uniqueness_scores << std::endl;
        std::cout << "Minimizers filtered by uniqueness: " << minimizers_filtered_by_uniqueness << std::endl;
        std::cout << "Unique minimizers (≥90%): " << unique_minimizers_count << std::endl;
        std::cout << "Rare minimizers (≤3 occurrences): " << rare_minimizers_count << std::endl;
        std::cout << "Reliable minimizers: " << reliable_minimizers_count << std::endl;
        std::cout << "Average uniqueness score: " << std::fixed << std::setprecision(3) 
                  << average_uniqueness_score << std::endl;
        std::cout << "Estimated false positive reduction: " << std::setprecision(1) 
                  << (estimated_false_positive_reduction * 100) << "%" << std::endl;
    }
};
*/

// ===========================
// Configuration Helper Functions
// Add these to DatabaseBuilderFactory in gpu_database_builder_core.h
// ===========================

DatabaseBuildConfig DatabaseBuilderFactory::create_high_precision_config(const std::string& output_dir) {
    DatabaseBuildConfig config = create_default_config();
    
    // Enable all quality features
    config.enable_uniqueness_scoring = true;
    config.enable_uniqueness_filtering = true;
    config.uniqueness_threshold = 0.5f;          // Stricter threshold
    config.filter_extremely_common = true;
    config.false_positive_reduction_target = 0.8f; // Aggressive FP reduction
    
    // Enable enhanced analysis
    config.enable_phylogenetic_analysis = true;
    config.enable_quality_analysis = true;
    config.enable_debug_mode = true;
    
    std::cout << "Created high-precision configuration for reduced false positives" << std::endl;
    return config;
}

DatabaseBuildConfig DatabaseBuilderFactory::create_sensitivity_config(const std::string& output_dir) {
    DatabaseBuildConfig config = create_default_config();
    
    // Balance between precision and sensitivity
    config.enable_uniqueness_scoring = true;
    config.enable_uniqueness_filtering = false;   // Don't filter - just score
    config.uniqueness_threshold = 0.1f;           // Very permissive
    config.filter_extremely_common = false;
    config.false_positive_reduction_target = 0.3f; // Moderate FP reduction
    
    std::cout << "Created sensitivity-optimized configuration" << std::endl;
    return config;
}

// ===========================
// Command Line Integration
// Add these options to any CLI parser
// ===========================

void add_uniqueness_command_line_options() {
    /*
    parser.add_option("--enable-uniqueness", 
                     "Enable uniqueness scoring (default: true)")
          .default_value(true);
    
    parser.add_option("--uniqueness-threshold", 
                     "Minimum uniqueness score for filtering (0.0-1.0)")
          .default_value(0.3f);
    
    parser.add_option("--filter-common", 
                     "Filter extremely common minimizers (default: true)")
          .default_value(true);
    
    parser.add_option("--target-fp-reduction", 
                     "Target false positive reduction (0.0-1.0)")
          .default_value(0.5f);
    */
}

// ===========================
// Usage Example
// ===========================

void example_usage_with_uniqueness() {
    // Create builder with uniqueness scoring enabled
    DatabaseBuildConfig config = DatabaseBuilderFactory::create_high_precision_config("/output");
    GPUKrakenDatabaseBuilder builder("/output", config);
    
    // Build database with uniqueness scoring
    if (builder.build_database_from_concatenated_fna(
            "microbial_genomes.fna",
            "taxonomy/nodes.dmp", 
            "taxonomy/names.dmp")) {
        
        std::cout << "Database built successfully with uniqueness scoring!" << std::endl;
        
        // Print enhanced statistics including uniqueness
        const auto& stats = builder.get_enhanced_statistics();
        stats.print_uniqueness_stats();
        
        // Estimate classification improvement
        std::cout << "\nExpected improvements:" << std::endl;
        std::cout << "- False positive reduction: " 
                  << (stats.estimated_false_positive_reduction * 100) << "%" << std::endl;
        std::cout << "- Reliable minimizers for classification: " 
                  << stats.reliable_minimizers_count << std::endl;
    }
}

// ===========================
// Testing and Validation
// ===========================

bool test_uniqueness_implementation() {
    std::cout << "Testing uniqueness score implementation..." << std::endl;
    
    // Create test minimizer hits
    std::vector<GPUMinimizerHit> test_hits = {
        {0x1111111111111111, 1000, 10, 0, 0, 0, 0, 0},  // Common minimizer
        {0x2222222222222222, 1001, 20, 0, 0, 0, 0, 0},  // Rare minimizer  
        {0x3333333333333333, 1002, 30, 0, 0, 0, 0, 0},  // Unique minimizer
    };
    
    // Simulate different occurrence counts
    std::vector<uint32_t> genome_taxon_ids = {1000, 1001, 1002};
    
    // Allocate GPU memory
    GPUMinimizerHit* d_hits;
    cudaMalloc(&d_hits, test_hits.size() * sizeof(GPUMinimizerHit));
    cudaMemcpy(d_hits, test_hits.data(), test_hits.size() * sizeof(GPUMinimizerHit), 
               cudaMemcpyHostToDevice);
    
    // Test uniqueness calculation
    bool success = compute_and_encode_uniqueness_scores(
        d_hits, test_hits.size(), genome_taxon_ids, 3.0f
    );
    
    if (success) {
        // Copy results back and verify
        cudaMemcpy(test_hits.data(), d_hits, test_hits.size() * sizeof(GPUMinimizerHit),
                   cudaMemcpyDeviceToHost);
        
        std::cout << "Test Results:" << std::endl;
        for (size_t i = 0; i < test_hits.size(); i++) {
            uint8_t uniqueness_cat = UniquenessUtils::get_uniqueness_category(test_hits[i].feature_flags);
            bool is_unique = UniquenessUtils::is_unique_minimizer(test_hits[i].feature_flags);
            bool is_reliable = UniquenessUtils::is_reliable_minimizer(test_hits[i].feature_flags);
            
            std::cout << "  Minimizer " << i << ": uniqueness=" << (int)uniqueness_cat 
                      << ", unique=" << is_unique << ", reliable=" << is_reliable << std::endl;
        }
    }
    
    cudaFree(d_hits);
    return success;
}
