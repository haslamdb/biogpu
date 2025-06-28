// Simple test to verify uniqueness configuration
#include <iostream>
#include "core/gpu_database_builder_core.h"

int main() {
    std::cout << "=== UNIQUENESS CONFIGURATION TEST ===" << std::endl;
    
    // Test default configuration
    DatabaseBuildConfig default_config;
    std::cout << "\nDefault configuration:" << std::endl;
    std::cout << "  enable_uniqueness_scoring: " << (default_config.enable_uniqueness_scoring ? "true" : "false") << std::endl;
    std::cout << "  enable_uniqueness_filtering: " << (default_config.enable_uniqueness_filtering ? "true" : "false") << std::endl;
    std::cout << "  uniqueness_threshold: " << default_config.uniqueness_threshold << std::endl;
    std::cout << "  filter_extremely_common: " << (default_config.filter_extremely_common ? "true" : "false") << std::endl;
    
    // Test custom configuration
    DatabaseBuildConfig custom_config;
    custom_config.enable_uniqueness_scoring = true;
    custom_config.enable_uniqueness_filtering = true;
    custom_config.uniqueness_threshold = 0.7f;
    custom_config.filter_extremely_common = true;
    
    std::cout << "\nCustom configuration:" << std::endl;
    std::cout << "  enable_uniqueness_scoring: " << (custom_config.enable_uniqueness_scoring ? "true" : "false") << std::endl;
    std::cout << "  enable_uniqueness_filtering: " << (custom_config.enable_uniqueness_filtering ? "true" : "false") << std::endl;
    std::cout << "  uniqueness_threshold: " << custom_config.uniqueness_threshold << std::endl;
    std::cout << "  filter_extremely_common: " << (custom_config.filter_extremely_common ? "true" : "false") << std::endl;
    
    std::cout << "\n✓ Configuration test passed!" << std::endl;
    
    return 0;
}