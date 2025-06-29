// test_ml_structures.cu
// Standalone test for ML-enhanced 24-byte minimizer structures
// No external dependencies - tests core type definitions only

#include <iostream>
#include <vector>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <filesystem>
#include "gpu_kraken_types.h"

// Simple test result tracking
struct TestResult {
    std::string test_name;
    bool passed;
    std::string details;
};

// Display test results
void print_results(const std::vector<TestResult>& results) {
    std::cout << "\n=== TEST RESULTS SUMMARY ===" << std::endl;
    int passed = 0;
    int total = results.size();
    
    for (const auto& result : results) {
        std::cout << std::setw(40) << std::left << result.test_name << ": ";
        if (result.passed) {
            std::cout << "✓ PASSED";
            passed++;
        } else {
            std::cout << "✗ FAILED";
        }
        if (!result.details.empty()) {
            std::cout << " (" << result.details << ")";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\nTotal: " << passed << "/" << total << " tests passed" << std::endl;
    
    if (passed == total) {
        std::cout << "\n✅ ALL TESTS PASSED ✅" << std::endl;
    } else {
        std::cout << "\n❌ SOME TESTS FAILED ❌" << std::endl;
    }
}

int main() {
    std::cout << "=== ML-ENHANCED MINIMIZER STRUCTURE TEST ===" << std::endl;
    std::cout << "Testing 24-byte GPUMinimizerHit structure with ML features\n" << std::endl;
    
    std::vector<TestResult> results;
    
    // Test 1: Structure Size
    {
        TestResult test{"Structure Size (GPUMinimizerHit)", false, ""};
        size_t actual_size = sizeof(GPUMinimizerHit);
        size_t expected_size = 24;
        
        test.passed = (actual_size == expected_size);
        test.details = std::to_string(actual_size) + " bytes";
        
        std::cout << "Test 1: Structure Size" << std::endl;
        std::cout << "  GPUMinimizerHit size: " << actual_size << " bytes (expected: " << expected_size << ")" << std::endl;
        std::cout << "  Field layout:" << std::endl;
        std::cout << "    minimizer_hash: offset " << offsetof(GPUMinimizerHit, minimizer_hash) << ", size " << sizeof(GPUMinimizerHit::minimizer_hash) << std::endl;
        std::cout << "    genome_id:      offset " << offsetof(GPUMinimizerHit, genome_id) << ", size " << sizeof(GPUMinimizerHit::genome_id) << std::endl;
        std::cout << "    position:       offset " << offsetof(GPUMinimizerHit, position) << ", size " << sizeof(GPUMinimizerHit::position) << std::endl;
        std::cout << "    strand:         offset " << offsetof(GPUMinimizerHit, strand) << ", size " << sizeof(GPUMinimizerHit::strand) << std::endl;
        std::cout << "    taxon_id:       offset " << offsetof(GPUMinimizerHit, taxon_id) << ", size " << sizeof(GPUMinimizerHit::taxon_id) << std::endl;
        std::cout << "    ml_weight:      offset " << offsetof(GPUMinimizerHit, ml_weight) << ", size " << sizeof(GPUMinimizerHit::ml_weight) << std::endl;
        std::cout << "    feature_flags:  offset " << offsetof(GPUMinimizerHit, feature_flags) << ", size " << sizeof(GPUMinimizerHit::feature_flags) << std::endl;
        
        results.push_back(test);
    }
    
    // Test 2: StreamlinedMinimizerMetadata Size
    {
        TestResult test{"Structure Size (StreamlinedMinimizerMetadata)", false, ""};
        size_t actual_size = sizeof(StreamlinedMinimizerMetadata);
        size_t expected_size = 32;
        
        test.passed = (actual_size == expected_size);
        test.details = std::to_string(actual_size) + " bytes";
        
        std::cout << "\nTest 2: StreamlinedMinimizerMetadata Size" << std::endl;
        std::cout << "  Size: " << actual_size << " bytes (expected: " << expected_size << ")" << std::endl;
        
        results.push_back(test);
    }
    
    // Test 3: ML Weight Encoding
    {
        TestResult test{"ML Weight Encoding", true, ""};
        std::cout << "\nTest 3: ML Weight Encoding/Decoding" << std::endl;
        
        std::vector<float> test_values = {0.0f, 0.1f, 0.5f, 0.9f, 1.0f};
        for (float original : test_values) {
            uint16_t encoded = MinimizerFlags::float_to_ml_weight(original);
            float decoded = MinimizerFlags::ml_weight_to_float(encoded);
            float error = std::abs(decoded - original);
            
            std::cout << "  " << std::fixed << std::setprecision(4) 
                      << original << " -> " << encoded << " -> " << decoded 
                      << " (error: " << error << ")" << std::endl;
            
            if (error > 0.0001f) {
                test.passed = false;
                test.details = "Excessive error in conversion";
            }
        }
        
        results.push_back(test);
    }
    
    // Test 4: GC Content Categories
    {
        TestResult test{"GC Content Categories", true, ""};
        std::cout << "\nTest 4: GC Content Category Encoding" << std::endl;
        
        for (uint8_t category = 0; category < 8; category++) {
            uint16_t flags = 0;
            flags = MinimizerFlags::set_gc_content_category(flags, category);
            uint8_t retrieved = MinimizerFlags::get_gc_content_category(flags);
            
            std::cout << "  Category " << (int)category << " -> flags: 0x" 
                      << std::hex << flags << std::dec 
                      << " -> retrieved: " << (int)retrieved << std::endl;
            
            if (retrieved != category) {
                test.passed = false;
                test.details = "Category encoding/decoding mismatch";
            }
        }
        
        results.push_back(test);
    }
    
    // Test 5: Complexity Scores
    {
        TestResult test{"Complexity Scores", true, ""};
        std::cout << "\nTest 5: Complexity Score Encoding" << std::endl;
        
        for (uint8_t score = 0; score < 8; score++) {
            uint16_t flags = 0;
            flags = MinimizerFlags::set_complexity_score(flags, score);
            uint8_t retrieved = MinimizerFlags::get_complexity_score(flags);
            
            std::cout << "  Score " << (int)score << " -> flags: 0x" 
                      << std::hex << flags << std::dec 
                      << " -> retrieved: " << (int)retrieved << std::endl;
            
            if (retrieved != score) {
                test.passed = false;
                test.details = "Score encoding/decoding mismatch";
            }
        }
        
        results.push_back(test);
    }
    
    // Test 6: Flag Combinations
    {
        TestResult test{"Combined Flag Encoding", true, ""};
        std::cout << "\nTest 6: Combined Flag Encoding" << std::endl;
        
        uint16_t flags = 0;
        flags = MinimizerFlags::set_gc_content_category(flags, 5);
        flags = MinimizerFlags::set_complexity_score(flags, 3);
        flags = MinimizerFlags::set_position_bias(flags, true);
        flags = MinimizerFlags::set_contamination_risk(flags, true);
        
        std::cout << "  Combined flags: 0x" << std::hex << flags << std::dec << std::endl;
        std::cout << "  GC category: " << (int)MinimizerFlags::get_gc_content_category(flags) << " (expected: 5)" << std::endl;
        std::cout << "  Complexity: " << (int)MinimizerFlags::get_complexity_score(flags) << " (expected: 3)" << std::endl;
        std::cout << "  Position bias: " << MinimizerFlags::has_position_bias(flags) << " (expected: 1)" << std::endl;
        std::cout << "  Contamination: " << MinimizerFlags::has_contamination_risk(flags) << " (expected: 1)" << std::endl;
        
        if (MinimizerFlags::get_gc_content_category(flags) != 5 ||
            MinimizerFlags::get_complexity_score(flags) != 3 ||
            !MinimizerFlags::has_position_bias(flags) ||
            !MinimizerFlags::has_contamination_risk(flags)) {
            test.passed = false;
            test.details = "Combined encoding failed";
        }
        
        results.push_back(test);
    }
    
    // Test 7: Database Version Check
    {
        TestResult test{"Database Version", true, ""};
        std::cout << "\nTest 7: Database Version Check" << std::endl;
        
        // Create a temporary file to test version writing
        std::string test_file = "test_version.bin";
        uint64_t version = 3;  // Version 3 for ML features
        
        // Write version
        std::ofstream out(test_file, std::ios::binary);
        out.write(reinterpret_cast<const char*>(&version), sizeof(uint64_t));
        out.close();
        
        // Read version back
        std::ifstream in(test_file, std::ios::binary);
        uint64_t read_version;
        in.read(reinterpret_cast<char*>(&read_version), sizeof(uint64_t));
        in.close();
        
        std::cout << "  Written version: " << version << std::endl;
        std::cout << "  Read version: " << read_version << std::endl;
        
        if (read_version != version) {
            test.passed = false;
            test.details = "Version mismatch";
        }
        
        // Clean up
        std::filesystem::remove(test_file);
        
        results.push_back(test);
    }
    
    // Test 8: Size Comparison
    {
        TestResult test{"Size Comparison (20 vs 24 byte)", true, ""};
        std::cout << "\nTest 8: Size Comparison" << std::endl;
        
        size_t old_size = 20;  // Assumed old structure size
        size_t new_size = sizeof(GPUMinimizerHit);
        float increase = ((float)new_size / old_size - 1.0f) * 100.0f;
        
        std::cout << "  Old structure: " << old_size << " bytes" << std::endl;
        std::cout << "  New structure: " << new_size << " bytes" << std::endl;
        std::cout << "  Size increase: " << std::fixed << std::setprecision(1) << increase << "%" << std::endl;
        
        test.details = "+" + std::to_string((int)increase) + "% size";
        results.push_back(test);
    }
    
    // Test 9: Feature Distribution Demo
    {
        TestResult test{"Feature Distribution Demo", true, ""};
        std::cout << "\nTest 9: Feature Distribution Statistics" << std::endl;
        
        // Simulate feature distribution
        std::vector<int> gc_dist(8, 0);
        std::vector<int> complexity_dist(8, 0);
        int contamination_count = 0;
        
        // Generate test data
        for (int i = 0; i < 1000; i++) {
            uint8_t gc_cat = i % 8;
            uint8_t complexity = (i / 8) % 8;
            
            gc_dist[gc_cat]++;
            complexity_dist[complexity]++;
            
            if (i % 50 == 0) {
                contamination_count++;
            }
        }
        
        std::cout << "  GC Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            std::cout << "    Category " << i << ": " << gc_dist[i] << " minimizers" << std::endl;
        }
        
        std::cout << "  Contamination markers: " << contamination_count << "/1000 ("
                  << (contamination_count * 100.0 / 1000.0) << "%)" << std::endl;
        
        results.push_back(test);
    }
    
    // Display final results
    print_results(results);
    
    // Return appropriate exit code
    for (const auto& result : results) {
        if (!result.passed) {
            return 1;
        }
    }
    
    return 0;
}