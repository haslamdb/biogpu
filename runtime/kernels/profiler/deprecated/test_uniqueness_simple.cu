// Simple test for uniqueness scoring integration
#include <iostream>
#include <cuda_runtime.h>
#include "features/working_uniqueness_implementation.cu"

int main() {
    std::cout << "=== UNIQUENESS SCORING TEST ===" << std::endl;
    
    // Run the built-in test
    bool test_result = test_uniqueness_implementation();
    
    if (test_result) {
        std::cout << "\n✓ Uniqueness scoring test PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ Uniqueness scoring test FAILED!" << std::endl;
        return 1;
    }
}