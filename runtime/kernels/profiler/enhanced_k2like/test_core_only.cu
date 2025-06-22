// test_core_only.cu
// Test program for the core module

#include <iostream>
#include "core/gpu_database_builder_core.h"
#include "gpu_kraken_types.h"

int main() {
    std::cout << "=== BioGPU Core Module Test ===" << std::endl;
    
    try {
        // Test basic configuration
        DatabaseBuildConfig config = GPUKrakenDatabaseBuilder::create_default_config();
        std::cout << "✓ Default config created successfully" << std::endl;
        std::cout << "  k-value: " << config.k_value << std::endl;
        std::cout << "  ell-value: " << config.ell_value << std::endl;
        std::cout << "  minimizer capacity: " << config.memory_config.minimizer_capacity << std::endl;
        
        // Test builder creation
        std::cout << "\nTesting core builder creation..." << std::endl;
        GPUKrakenDatabaseBuilder builder("./test_output", config);
        std::cout << "✓ Core builder created successfully" << std::endl;
        
        // Test basic validation
        if (builder.validate_configuration()) {
            std::cout << "✓ Configuration validation passed" << std::endl;
        } else {
            std::cout << "✗ Configuration validation failed" << std::endl;
            return 1;
        }
        
        // Test factory patterns
        auto standard_builder = DatabaseBuilderFactory::create_standard_builder("./test_standard");
        if (standard_builder) {
            std::cout << "✓ Factory pattern works" << std::endl;
        }
        
        std::cout << "\n=== ALL CORE TESTS PASSED ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}