// fastq_pipeline_debug.cpp
// Debug version to identify performance bottlenecks
#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <cuda_runtime.h>
#include "fastq_processing.h"

class Timer {
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start;
    std::string name;
public:
    Timer(const std::string& n) : name(n), start(Clock::now()) {}
    ~Timer() {
        auto end = Clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << name << ": " << duration.count() << " ms\n";
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fastq_file>\n";
        return 1;
    }
    
    std::string fastq_file = argv[1];
    
    // Check CUDA devices
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    std::cout << "CUDA Devices found: " << deviceCount << "\n";
    
    for (int i = 0; i < deviceCount; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        std::cout << "Device " << i << ": " << prop.name 
                  << " (Compute " << prop.major << "." << prop.minor << ")\n";
    }
    
    try {
        // Test small batch first
        std::cout << "\n=== Testing with small batch ===\n";
        {
            Timer t("Small batch test");
            biogpu::GPUMinimizerPipeline pipeline(31, 15, 1000, 1);  // Small batch, 1 thread
            
            size_t count = 0;
            pipeline.process_file(fastq_file, 
                [&count](const ReadBatch& batch, 
                        const std::vector<std::vector<Minimizer>>& minimizers) {
                    count += batch.size();
                    if (count >= 10000) return;  // Stop after 10k reads
                });
        }
        
        // Test file reading speed
        std::cout << "\n=== Testing file reading speed ===\n";
        {
            Timer t("Reading 100k reads");
            biogpu::FastqReader reader(fastq_file);
            ReadBatch batch;
            size_t total = 0;
            
            while (reader.read_batch(batch, 10000) && total < 100000) {
                total += batch.size();
            }
            std::cout << "Read " << total << " sequences\n";
        }
        
        // Test GPU processing speed
        std::cout << "\n=== Testing GPU processing speed ===\n";
        {
            // Generate test data
            std::vector<std::string> test_sequences;
            for (int i = 0; i < 10000; i++) {
                std::string seq(150, 'A');
                for (int j = 0; j < 150; j++) {
                    seq[j] = "ACGT"[rand() % 4];
                }
                test_sequences.push_back(seq);
            }
            
            Timer t("GPU processing 10k sequences");
            MinimizerExtractor extractor(31, 15);
            auto minimizers = extractor.extract_minimizers(test_sequences);
            std::cout << "Extracted minimizers from " << minimizers.size() << " sequences\n";
        }
        
        // Test with different batch sizes
        std::cout << "\n=== Testing different batch sizes ===\n";
        for (size_t batch_size : {1000, 5000, 10000, 20000}) {
            Timer t("Batch size " + std::to_string(batch_size));
            
            biogpu::GPUMinimizerPipeline pipeline(31, 15, batch_size, 1);
            size_t count = 0;
            
            pipeline.process_file(fastq_file, 
                [&count](const ReadBatch& batch, 
                        const std::vector<std::vector<Minimizer>>& minimizers) {
                    count += batch.size();
                    if (count >= 50000) return;  // Stop after 50k reads
                });
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}