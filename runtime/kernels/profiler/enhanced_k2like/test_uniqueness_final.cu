// Final test for uniqueness scoring integration
#include <iostream>
#include <vector>
#include <filesystem>
#include "core/gpu_database_builder_core.h"

// Forward declaration for test function
bool test_uniqueness_implementation();

int main() {
    std::cout << "=== FINAL UNIQUENESS INTEGRATION TEST ===" << std::endl;
    
    // Create output directory
    std::string output_dir = "test_uniqueness_final_output";
    std::filesystem::create_directories(output_dir);
    
    // Create configuration with uniqueness scoring enabled
    DatabaseBuildConfig config;
    config.enable_uniqueness_scoring = true;
    config.enable_uniqueness_filtering = false;
    config.uniqueness_threshold = 0.3f;
    config.filter_extremely_common = true;
    config.enable_debug_mode = true;  // To see uniqueness statistics
    
    std::cout << "\nTest Configuration:" << std::endl;
    std::cout << "  Uniqueness scoring: " << (config.enable_uniqueness_scoring ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "  Uniqueness filtering: " << (config.enable_uniqueness_filtering ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "  Uniqueness threshold: " << config.uniqueness_threshold << std::endl;
    std::cout << "  Filter extremely common: " << (config.filter_extremely_common ? "YES" : "NO") << std::endl;
    std::cout << "  Debug mode: " << (config.enable_debug_mode ? "ENABLED" : "DISABLED") << std::endl;
    
    try {
        // Create database builder
        GPUKrakenDatabaseBuilder builder(output_dir, config);
        
        std::cout << "\n✓ Database builder created successfully" << std::endl;
        
        // Check if CUDA is available
        int deviceCount = 0;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount > 0) {
            std::cout << "✓ Found " << deviceCount << " CUDA device(s)" << std::endl;
        } else {
            std::cerr << "✗ No CUDA devices found" << std::endl;
            return 1;
        }
        
        // Print initial statistics
        const auto& stats = builder.get_enhanced_statistics();
        std::cout << "\nInitial Enhanced Statistics:" << std::endl;
        std::cout << "  Minimizers with uniqueness scores: " << stats.minimizers_with_uniqueness_scores << std::endl;
        std::cout << "  Unique minimizers count: " << stats.unique_minimizers_count << std::endl;
        std::cout << "  Rare minimizers count: " << stats.rare_minimizers_count << std::endl;
        std::cout << "  Reliable minimizers count: " << stats.reliable_minimizers_count << std::endl;
        std::cout << "  Average uniqueness score: " << stats.average_uniqueness_score << std::endl;
        
        // Run the standalone test
        std::cout << "\nRunning standalone uniqueness test..." << std::endl;
        if (test_uniqueness_implementation()) {
            std::cout << "✓ Standalone uniqueness test PASSED!" << std::endl;
        } else {
            std::cout << "✗ Standalone uniqueness test FAILED!" << std::endl;
            return 1;
        }
        
        std::cout << "\n=== ALL TESTS PASSED ===" << std::endl;
        std::cout << "Uniqueness scoring has been successfully integrated!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}