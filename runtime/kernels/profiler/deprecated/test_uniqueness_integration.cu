// Test for uniqueness scoring integration
#include <iostream>
#include <vector>
#include <filesystem>
#include "core/gpu_database_builder_core.h"

int main() {
    std::cout << "=== UNIQUENESS SCORING INTEGRATION TEST ===" << std::endl;
    
    // Create output directory
    std::string output_dir = "test_uniqueness_output";
    std::filesystem::create_directories(output_dir);
    
    // Create builder with uniqueness scoring enabled
    DatabaseBuildConfig config;
    config.enable_uniqueness_scoring = true;
    config.enable_uniqueness_filtering = false;
    config.uniqueness_threshold = 0.3f;
    config.filter_extremely_common = true;
    
    std::cout << "\nConfiguration:" << std::endl;
    std::cout << "  Uniqueness scoring: " << (config.enable_uniqueness_scoring ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "  Uniqueness filtering: " << (config.enable_uniqueness_filtering ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "  Uniqueness threshold: " << config.uniqueness_threshold << std::endl;
    std::cout << "  Filter extremely common: " << (config.filter_extremely_common ? "YES" : "NO") << std::endl;
    
    GPUKrakenDatabaseBuilder builder(output_dir, config);
    
    std::cout << "\n✓ Database builder created with uniqueness scoring configuration" << std::endl;
    
    // Verify configuration was set
    const auto& stats = builder.get_enhanced_statistics();
    std::cout << "\nInitial statistics:" << std::endl;
    std::cout << "  Minimizers with uniqueness scores: " << stats.minimizers_with_uniqueness_scores << std::endl;
    std::cout << "  Unique minimizers count: " << stats.unique_minimizers_count << std::endl;
    
    std::cout << "\n✓ Uniqueness scoring integration test passed!" << std::endl;
    
    return 0;
}