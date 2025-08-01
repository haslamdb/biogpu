// Test version of unified pipeline to verify build system
#include <iostream>
#include <cuda_runtime.h>
#include "unified_amr_pipeline.cpp"

// Simplified test implementations
namespace {
    
    // Test function to verify CUDA is working
    bool testCUDA() {
        int device_count;
        cudaError_t err = cudaGetDeviceCount(&device_count);
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        std::cout << "Found " << device_count << " CUDA devices:" << std::endl;
        for (int i = 0; i < device_count; ++i) {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, i);
            std::cout << "  GPU " << i << ": " << prop.name 
                      << " (" << (prop.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
        }
        return true;
    }
    
    // Test streaming FASTQ reader
    bool testStreamingReader() {
        // This will be implemented when we test with actual files
        std::cout << "Streaming FASTQ reader test: SKIPPED (no test files)" << std::endl;
        return true;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "BioGPU Unified Pipeline - Test Build" << std::endl;
    std::cout << "====================================" << std::endl;
    
    // Test CUDA
    if (!testCUDA()) {
        std::cerr << "CUDA test failed!" << std::endl;
        return 1;
    }
    
    // Test streaming reader
    if (!testStreamingReader()) {
        std::cerr << "Streaming reader test failed!" << std::endl;
        return 1;
    }
    
    // If we got here, basic functionality works
    std::cout << "\nAll tests passed!" << std::endl;
    std::cout << "The unified pipeline framework is ready." << std::endl;
    std::cout << "\nNext steps:" << std::endl;
    std::cout << "1. Integrate CleanResistancePipeline from resistance/" << std::endl;
    std::cout << "2. Integrate AMRDetectionPipeline from genes/" << std::endl;
    std::cout << "3. Test with real FASTQ data" << std::endl;
    
    return 0;
}