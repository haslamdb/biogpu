// test_memory_build.cu
// Simple test to verify gpu_memory_manager builds correctly

#include "memory/gpu_memory_manager.h"
#include <iostream>

int main() {
    std::cout << "Testing GPU Memory Manager build..." << std::endl;
    
    // Test instantiation with default config
    {
        MemoryConfig config;
        GPUMemoryManager manager(config);
        std::cout << "✓ GPUMemoryManager instantiated with default config" << std::endl;
    }
    
    // Test instantiation with custom config
    {
        MemoryConfig config;
        config.max_memory_fraction = 90;
        config.reserved_memory_mb = 256;
        GPUMemoryManager manager(config);
        std::cout << "✓ GPUMemoryManager instantiated with custom config" << std::endl;
    }
    
    // Test memory pool instantiation
    {
        GPUMemoryPool pool(512); // 512 MB
        std::cout << "✓ GPUMemoryPool instantiated" << std::endl;
    }
    
    std::cout << "All instantiation tests passed!" << std::endl;
    return 0;
}