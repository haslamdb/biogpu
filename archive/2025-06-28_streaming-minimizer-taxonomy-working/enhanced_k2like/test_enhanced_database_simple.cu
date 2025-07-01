// test_enhanced_database_simple.cu
// Simplified test for ML-enhanced database building
// Tests core functionality without all module dependencies

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <random>
#include <iomanip>
#include <cuda_runtime.h>
#include <cstring>

#include "gpu_kraken_types.h"
#include "output/database_serializer.h"

// Test configuration
struct TestConfig {
    std::string output_dir = "test_output_simple";
    int num_test_minimizers = 10000;
    bool verbose = true;
};

// Test results
struct TestResults {
    bool structure_size_test_passed = false;
    bool feature_encoding_test_passed = false;
    bool serialization_test_passed = false;
    bool ml_fields_test_passed = false;
    
    size_t old_structure_size = 0;
    size_t new_structure_size = 0;
    int gc_categories_tested = 0;
    int complexity_scores_tested = 0;
    int contamination_markers = 0;
    
    void print_summary() const {
        std::cout << "\n=== TEST RESULTS SUMMARY ===" << std::endl;
        std::cout << "Structure Size Test: " << (structure_size_test_passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << "Feature Encoding Test: " << (feature_encoding_test_passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << "Serialization Test: " << (serialization_test_passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << "ML Fields Test: " << (ml_fields_test_passed ? "PASSED" : "FAILED") << std::endl;
        
        std::cout << "\n=== SIZE COMPARISON ===" << std::endl;
        std::cout << "Old structure (20 bytes assumed): " << old_structure_size << " bytes" << std::endl;
        std::cout << "New structure (GPUMinimizerHit): " << new_structure_size << " bytes" << std::endl;
        std::cout << "Size increase: " << ((float)new_structure_size / old_structure_size - 1.0f) * 100 << "%" << std::endl;
        
        std::cout << "\n=== FEATURES TESTED ===" << std::endl;
        std::cout << "GC categories: " << gc_categories_tested << "/8" << std::endl;
        std::cout << "Complexity scores: " << complexity_scores_tested << "/8" << std::endl;
        std::cout << "Contamination markers: " << contamination_markers << std::endl;
    }
};

// Test 1: Structure Size Validation
bool test_structure_sizes(TestResults& results) {
    std::cout << "\n=== TEST 1: STRUCTURE SIZE VALIDATION ===" << std::endl;
    
    // Check GPUMinimizerHit size (should be 24 bytes)
    size_t gpu_hit_size = sizeof(GPUMinimizerHit);
    std::cout << "GPUMinimizerHit size: " << gpu_hit_size << " bytes" << std::endl;
    
    // Check field offsets
    std::cout << "Field offsets:" << std::endl;
    std::cout << "  minimizer_hash: " << offsetof(GPUMinimizerHit, minimizer_hash) << std::endl;
    std::cout << "  genome_id: " << offsetof(GPUMinimizerHit, genome_id) << std::endl;
    std::cout << "  position: " << offsetof(GPUMinimizerHit, position) << std::endl;
    std::cout << "  strand: " << offsetof(GPUMinimizerHit, strand) << std::endl;
    std::cout << "  taxon_id: " << offsetof(GPUMinimizerHit, taxon_id) << std::endl;
    std::cout << "  ml_weight: " << offsetof(GPUMinimizerHit, ml_weight) << std::endl;
    std::cout << "  feature_flags: " << offsetof(GPUMinimizerHit, feature_flags) << std::endl;
    
    // Check StreamlinedMinimizerMetadata size
    size_t metadata_size = sizeof(StreamlinedMinimizerMetadata);
    std::cout << "\nStreamlinedMinimizerMetadata size: " << metadata_size << " bytes" << std::endl;
    
    results.old_structure_size = 20;  // Assumed old size
    results.new_structure_size = gpu_hit_size;
    
    // Validate sizes
    if (gpu_hit_size != 24) {
        std::cerr << "ERROR: GPUMinimizerHit size is " << gpu_hit_size << ", expected 24" << std::endl;
        return false;
    }
    
    if (metadata_size != 32) {
        std::cerr << "ERROR: StreamlinedMinimizerMetadata size is " << metadata_size << ", expected 32" << std::endl;
        return false;
    }
    
    std::cout << "✓ Structure sizes validated" << std::endl;
    results.structure_size_test_passed = true;
    return true;
}

// Test 2: Feature Encoding
bool test_feature_encoding(TestResults& results) {
    std::cout << "\n=== TEST 2: FEATURE ENCODING TEST ===" << std::endl;
    
    // Test all GC categories
    for (uint8_t gc_cat = 0; gc_cat < 8; gc_cat++) {
        uint16_t flags = 0;
        flags = MinimizerFlags::set_gc_content_category(flags, gc_cat);
        uint8_t retrieved = MinimizerFlags::get_gc_content_category(flags);
        
        if (retrieved != gc_cat) {
            std::cerr << "ERROR: GC category encoding failed. Set " << (int)gc_cat 
                      << ", got " << (int)retrieved << std::endl;
            return false;
        }
        results.gc_categories_tested++;
    }
    
    // Test all complexity scores
    for (uint8_t complexity = 0; complexity < 8; complexity++) {
        uint16_t flags = 0;
        flags = MinimizerFlags::set_complexity_score(flags, complexity);
        uint8_t retrieved = MinimizerFlags::get_complexity_score(flags);
        
        if (retrieved != complexity) {
            std::cerr << "ERROR: Complexity score encoding failed. Set " << (int)complexity 
                      << ", got " << (int)retrieved << std::endl;
            return false;
        }
        results.complexity_scores_tested++;
    }
    
    // Test position bias flag
    uint16_t flags = 0;
    flags = MinimizerFlags::set_position_bias(flags, true);
    if (!MinimizerFlags::has_position_bias(flags)) {
        std::cerr << "ERROR: Position bias flag not set correctly" << std::endl;
        return false;
    }
    
    // Test contamination risk flag
    flags = 0;
    flags = MinimizerFlags::set_contamination_risk(flags, true);
    if (!MinimizerFlags::has_contamination_risk(flags)) {
        std::cerr << "ERROR: Contamination risk flag not set correctly" << std::endl;
        return false;
    }
    
    // Test combined flags
    flags = 0;
    flags = MinimizerFlags::set_gc_content_category(flags, 5);
    flags = MinimizerFlags::set_complexity_score(flags, 3);
    flags = MinimizerFlags::set_contamination_risk(flags, true);
    
    if (MinimizerFlags::get_gc_content_category(flags) != 5 ||
        MinimizerFlags::get_complexity_score(flags) != 3 ||
        !MinimizerFlags::has_contamination_risk(flags)) {
        std::cerr << "ERROR: Combined flag encoding failed" << std::endl;
        return false;
    }
    
    std::cout << "✓ Feature encoding test passed" << std::endl;
    results.feature_encoding_test_passed = true;
    return true;
}

// Test 3: ML Fields
bool test_ml_fields(TestResults& results) {
    std::cout << "\n=== TEST 3: ML FIELDS TEST ===" << std::endl;
    
    // Test ML weight encoding/decoding
    std::vector<float> test_confidences = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
    
    for (float confidence : test_confidences) {
        uint16_t encoded = MinimizerFlags::float_to_ml_weight(confidence);
        float decoded = MinimizerFlags::ml_weight_to_float(encoded);
        
        float error = std::abs(decoded - confidence);
        if (error > 0.0001f) {  // Allow small floating point error
            std::cerr << "ERROR: ML weight encoding failed. Original: " << confidence
                      << ", decoded: " << decoded << std::endl;
            return false;
        }
        
        std::cout << "  Confidence " << confidence << " -> " << encoded 
                  << " -> " << decoded << " ✓" << std::endl;
    }
    
    // Test StreamlinedMinimizerMetadata ML methods
    StreamlinedMinimizerMetadata metadata;
    metadata.set_ml_confidence(0.8f);
    float retrieved = metadata.get_ml_confidence();
    
    if (std::abs(retrieved - 0.8f) > 0.0001f) {
        std::cerr << "ERROR: StreamlinedMinimizerMetadata ML methods failed" << std::endl;
        return false;
    }
    
    std::cout << "✓ ML fields test passed" << std::endl;
    results.ml_fields_test_passed = true;
    return true;
}

// Test 4: Database Serialization
bool test_serialization(const TestConfig& config, TestResults& results) {
    std::cout << "\n=== TEST 4: DATABASE SERIALIZATION TEST ===" << std::endl;
    
    // Create output directory
    std::filesystem::create_directories(config.output_dir);
    
    // Create test data
    std::vector<PhylogeneticLCACandidate> test_candidates;
    
    for (int i = 0; i < config.num_test_minimizers; i++) {
        PhylogeneticLCACandidate candidate;
        candidate.minimizer_hash = 0x1000000000000000ULL + i;
        candidate.lca_taxon = 1000 + (i % 100);
        candidate.genome_count = 1 + (i % 10);
        candidate.uniqueness_score = (float)(i % 100) / 100.0f;
        candidate.phylogenetic_spread = i % 256;
        candidate.max_phylogenetic_distance = (i * 2) % 256;
        
        // Add some contributing species
        for (int j = 0; j < 3; j++) {
            candidate.contributing_species.push_back(1000 + j);
            candidate.genome_counts_per_species.push_back(1);
        }
        
        // Mark some as contamination
        if (i % 100 == 0) {
            results.contamination_markers++;
        }
        
        test_candidates.push_back(candidate);
    }
    
    // Test serialization
    try {
        EnhancedDatabaseSerializer serializer(config.output_dir);
        ContributingTaxaArrays dummy_taxa;
        std::unordered_map<uint32_t, std::string> taxon_names;
        EnhancedBuildStats stats;
        stats.total_sequences = 100;
        stats.total_bases = 1000000;
        stats.unique_minimizers = config.num_test_minimizers;
        
        if (!serializer.save_enhanced_database(test_candidates, dummy_taxa, taxon_names, stats)) {
            std::cerr << "Failed to save enhanced database" << std::endl;
            return false;
        }
        
        // Verify files
        std::vector<std::string> expected_files = {
            config.output_dir + "/enhanced_hash_table.bin",
            config.output_dir + "/ml_weights.bin",
            config.output_dir + "/feature_statistics.json",
            config.output_dir + "/contamination_markers.bin"
        };
        
        for (const auto& file : expected_files) {
            if (!std::filesystem::exists(file)) {
                std::cerr << "Missing file: " << file << std::endl;
                return false;
            }
            size_t file_size = std::filesystem::file_size(file);
            std::cout << "  Created: " << file << " (" << file_size << " bytes)" << std::endl;
        }
        
        // Read and verify version
        std::ifstream hash_file(config.output_dir + "/enhanced_hash_table.bin", std::ios::binary);
        uint64_t version;
        hash_file.read(reinterpret_cast<char*>(&version), sizeof(uint64_t));
        
        if (version != 3) {
            std::cerr << "ERROR: Database version is " << version << ", expected 3" << std::endl;
            return false;
        }
        
        std::cout << "  Database version: " << version << " ✓" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Serialization error: " << e.what() << std::endl;
        return false;
    }
    
    std::cout << "✓ Serialization test passed" << std::endl;
    results.serialization_test_passed = true;
    return true;
}

// Calculate GC content of a sequence
float calculate_gc_content(const std::string& seq) {
    int gc_count = 0;
    for (char c : seq) {
        if (c == 'G' || c == 'C' || c == 'g' || c == 'c') {
            gc_count++;
        }
    }
    return (float)gc_count / seq.length();
}

// Main test runner
int main(int argc, char* argv[]) {
    std::cout << "=== ENHANCED DATABASE BUILD SIMPLE TEST ===" << std::endl;
    std::cout << "Testing ML-enhanced 24-byte minimizer structure" << std::endl;
    
    TestConfig config;
    TestResults results;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--output-dir" && i + 1 < argc) {
            config.output_dir = argv[++i];
        } else if (arg == "--num-minimizers" && i + 1 < argc) {
            config.num_test_minimizers = std::stoi(argv[++i]);
        }
    }
    
    // Run tests
    bool all_passed = true;
    
    if (!test_structure_sizes(results)) {
        all_passed = false;
    }
    
    if (!test_feature_encoding(results)) {
        all_passed = false;
    }
    
    if (!test_ml_fields(results)) {
        all_passed = false;
    }
    
    if (!test_serialization(config, results)) {
        all_passed = false;
    }
    
    // Print summary
    results.print_summary();
    
    // Feature statistics demo
    std::cout << "\n=== FEATURE STATISTICS DEMO ===" << std::endl;
    
    // Example sequences with known properties
    std::vector<std::pair<std::string, std::string>> test_sequences = {
        {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "Low GC (0%)"},
        {"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", "High GC (100%)"},
        {"ATCGATCGATCGATCGATCGATCGATCGATCGAT", "Medium GC (50%)"},
        {"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", "Unknown bases"}
    };
    
    for (const auto& [seq, desc] : test_sequences) {
        float gc = calculate_gc_content(seq);
        uint8_t gc_category = (uint8_t)(gc * 8);
        if (gc_category > 7) gc_category = 7;
        
        std::cout << desc << ": GC=" << (gc * 100) << "%, Category=" << (int)gc_category << std::endl;
    }
    
    if (all_passed) {
        std::cout << "\n✅ ALL TESTS PASSED ✅" << std::endl;
        
        // Clean up test files
        if (std::filesystem::exists(config.output_dir)) {
            std::filesystem::remove_all(config.output_dir);
            std::cout << "Cleaned up test output directory" << std::endl;
        }
        
        return 0;
    } else {
        std::cout << "\n❌ SOME TESTS FAILED ❌" << std::endl;
        return 1;
    }
}