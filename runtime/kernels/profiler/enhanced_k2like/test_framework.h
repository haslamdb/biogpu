// test_framework.h
// Modular test framework that can be extended without rewriting

#ifndef BIOGPU_TEST_FRAMEWORK_H
#define BIOGPU_TEST_FRAMEWORK_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional>
#include <memory>
#include <chrono>
#include <iomanip>
#include <cuda_runtime.h>

#include "gpu_kraken_types.h"
#include "processing/minimizer_feature_extractor.h"

// Base test configuration that can be extended
struct TestConfig {
    // File paths
    std::string fna_file = "/home/david/Documents/Code/biogpu/data/test_50_genomes.fna";
    std::string temp_dir = "/tmp/biogpu_test";
    std::string output_dir = "./test_output";
    
    // Processing parameters
    size_t batch_size = 5;
    size_t max_minimizers_per_batch = 10000000;
    
    // Minimizer parameters
    MinimizerParams params = {
        31,     // k
        31,     // ell
        7,      // spaces
        0x3c8bfbb395c60474ULL  // xor_mask
    };
    
    // Feature flags - enable/disable specific features
    bool enable_gc_content = true;
    bool enable_complexity = true;
    bool enable_uniqueness = true;
    bool enable_cooccurrence = true;
    bool enable_contamination = false;  // Add when ready
    bool enable_taxonomy = true;
    
    // Debug options
    bool verbose = true;
    bool save_intermediate_results = false;
    int num_examples_to_print = 10;
};

// Feature test interface
class IFeatureTest {
public:
    virtual ~IFeatureTest() = default;
    virtual std::string getName() const = 0;
    virtual bool isEnabled(const TestConfig& config) const = 0;
    virtual bool initialize() = 0;
    virtual bool processOnGPU(GPUMinimizerHit* d_hits, size_t num_hits, 
                             const std::vector<uint32_t>& taxon_ids) = 0;
    virtual void analyzeResults(const std::vector<GPUMinimizerHit>& hits) = 0;
    virtual void printSummary() const = 0;
};

// Test result collector
class TestResults {
private:
    struct FeatureStats {
        std::string name;
        bool passed = false;
        double processing_time_ms = 0.0;
        std::unordered_map<std::string, double> metrics;
        std::vector<std::string> warnings;
        std::vector<std::string> errors;
    };
    
    std::vector<FeatureStats> feature_results_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
    
public:
    void startTiming() {
        start_time_ = std::chrono::high_resolution_clock::now();
    }
    
    double getElapsedMs() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(now - start_time_).count();
    }
    
    void recordFeatureResult(const std::string& name, bool passed, double time_ms) {
        FeatureStats stats;
        stats.name = name;
        stats.passed = passed;
        stats.processing_time_ms = time_ms;
        feature_results_.push_back(stats);
    }
    
    void addMetric(const std::string& feature, const std::string& metric, double value) {
        for (auto& stats : feature_results_) {
            if (stats.name == feature) {
                stats.metrics[metric] = value;
                return;
            }
        }
    }
    
    void addWarning(const std::string& feature, const std::string& warning) {
        for (auto& stats : feature_results_) {
            if (stats.name == feature) {
                stats.warnings.push_back(warning);
                return;
            }
        }
    }
    
    void printSummary() const {
        std::cout << "\n=== TEST RESULTS SUMMARY ===" << std::endl;
        std::cout << std::string(80, '-') << std::endl;
        
        int passed = 0, failed = 0;
        double total_time = 0.0;
        
        for (const auto& stats : feature_results_) {
            std::cout << std::left << std::setw(30) << stats.name 
                      << " | " << std::setw(8) << (stats.passed ? "PASSED" : "FAILED")
                      << " | " << std::fixed << std::setprecision(1) 
                      << std::setw(10) << stats.processing_time_ms << " ms" << std::endl;
            
            if (stats.passed) passed++; else failed++;
            total_time += stats.processing_time_ms;
            
            // Print metrics
            for (const auto& [metric, value] : stats.metrics) {
                std::cout << "    " << metric << ": " << value << std::endl;
            }
            
            // Print warnings
            for (const auto& warning : stats.warnings) {
                std::cout << "    ⚠️  " << warning << std::endl;
            }
        }
        
        std::cout << std::string(80, '-') << std::endl;
        std::cout << "Total: " << passed << " passed, " << failed << " failed | "
                  << "Total time: " << total_time << " ms" << std::endl;
    }
};

// Main test framework
class BioGPUTestFramework {
private:
    TestConfig config_;
    std::vector<std::unique_ptr<IFeatureTest>> feature_tests_;
    std::unique_ptr<MinimizerFeatureExtractor> feature_extractor_;
    TestResults results_;
    
    // Stored data for analysis
    std::vector<GPUMinimizerHit> all_minimizers_;
    std::vector<uint32_t> all_taxon_ids_;
    
public:
    explicit BioGPUTestFramework(const TestConfig& config = TestConfig())
        : config_(config) {
        feature_extractor_ = std::make_unique<MinimizerFeatureExtractor>(100000000, 1000);
    }
    
    void registerFeatureTest(std::unique_ptr<IFeatureTest> test) {
        feature_tests_.push_back(std::move(test));
    }
    
    TestConfig& getConfig() { return config_; }
    MinimizerFeatureExtractor* getFeatureExtractor() { return feature_extractor_.get(); }
    TestResults& getResults() { return results_; }
    
    bool runAllTests() {
        std::cout << "=== BIOGPU MODULAR TEST FRAMEWORK ===" << std::endl;
        
        // Check prerequisites
        if (!checkPrerequisites()) {
            return false;
        }
        
        // Extract minimizers with basic features
        results_.startTiming();
        if (!extractMinimizers()) {
            return false;
        }
        
        // Run each registered feature test
        for (auto& test : feature_tests_) {
            if (!test->isEnabled(config_)) {
                std::cout << "\nSkipping " << test->getName() << " (disabled)" << std::endl;
                continue;
            }
            
            std::cout << "\n--- Running " << test->getName() << " ---" << std::endl;
            
            auto start = std::chrono::high_resolution_clock::now();
            
            if (!test->initialize()) {
                results_.recordFeatureResult(test->getName(), false, 0.0);
                results_.addWarning(test->getName(), "Initialization failed");
                continue;
            }
            
            // Process on GPU
            bool success = false;
            if (!all_minimizers_.empty()) {
                // Allocate GPU memory for minimizers
                GPUMinimizerHit* d_minimizers;
                cudaMalloc(&d_minimizers, all_minimizers_.size() * sizeof(GPUMinimizerHit));
                cudaMemcpy(d_minimizers, all_minimizers_.data(), 
                          all_minimizers_.size() * sizeof(GPUMinimizerHit), 
                          cudaMemcpyHostToDevice);
                
                success = test->processOnGPU(d_minimizers, all_minimizers_.size(), all_taxon_ids_);
                
                if (success) {
                    // Copy results back
                    cudaMemcpy(all_minimizers_.data(), d_minimizers,
                              all_minimizers_.size() * sizeof(GPUMinimizerHit),
                              cudaMemcpyDeviceToHost);
                }
                
                cudaFree(d_minimizers);
            }
            
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
            
            results_.recordFeatureResult(test->getName(), success, elapsed_ms);
            
            if (success) {
                test->analyzeResults(all_minimizers_);
                test->printSummary();
            }
        }
        
        // Print final analysis
        printFinalAnalysis();
        
        // Print results summary
        results_.printSummary();
        
        return true;
    }
    
    const std::vector<GPUMinimizerHit>& getMinimizers() const { return all_minimizers_; }
    
private:
    bool checkPrerequisites() {
        // Check GPU
        int device_count;
        cudaGetDeviceCount(&device_count);
        if (device_count == 0) {
            std::cerr << "No CUDA devices found!" << std::endl;
            return false;
        }
        
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        std::cout << "GPU: " << prop.name << std::endl;
        std::cout << "Memory: " << prop.totalGlobalMem / (1024*1024) << " MB" << std::endl;
        
        // Check input file
        if (!std::filesystem::exists(config_.fna_file)) {
            std::cerr << "Input file not found: " << config_.fna_file << std::endl;
            return false;
        }
        
        // Create directories
        std::filesystem::create_directories(config_.temp_dir);
        if (config_.save_intermediate_results) {
            std::filesystem::create_directories(config_.output_dir);
        }
        
        return true;
    }
    
    bool extractMinimizers() {
        std::cout << "\nExtracting minimizers from " << config_.fna_file << std::endl;
        
        // This is where you'd call your existing extraction code
        // For now, returning true to show the framework structure
        
        // In practice, you'd have something like:
        // StreamingFnaProcessor processor(config_.fna_file, config_.temp_dir, config_.batch_size);
        // ... process batches and fill all_minimizers_ ...
        
        std::cout << "Extracted " << all_minimizers_.size() << " minimizers" << std::endl;
        return true;
    }
    
    void printFinalAnalysis() {
        if (all_minimizers_.empty()) return;
        
        std::cout << "\n=== FINAL COMBINED ANALYSIS ===" << std::endl;
        
        // Print examples
        if (config_.num_examples_to_print > 0) {
            std::cout << "\nExample minimizers with all features:" << std::endl;
            int examples = std::min(config_.num_examples_to_print, 
                                   (int)all_minimizers_.size());
            
            for (int i = 0; i < examples; i++) {
                const auto& hit = all_minimizers_[i];
                std::cout << "\nMinimizer " << i << ":" << std::endl;
                printMinimizerDetails(hit);
            }
        }
    }
    
    void printMinimizerDetails(const GPUMinimizerHit& hit) {
        std::cout << "  Hash: 0x" << std::hex << hit.minimizer_hash << std::dec << std::endl;
        std::cout << "  Position: " << hit.position << std::endl;
        std::cout << "  Taxon ID: " << hit.taxon_id << std::endl;
        
        // Print all enabled features
        if (config_.enable_gc_content) {
            std::cout << "  GC category: " << (hit.feature_flags & 0x7) << "/7" << std::endl;
        }
        if (config_.enable_complexity) {
            std::cout << "  Complexity: " << ((hit.feature_flags >> 3) & 0x7) << "/7" << std::endl;
        }
        // Add more features as needed
    }
};

// Convenience macro for creating feature tests
#define DEFINE_FEATURE_TEST(ClassName, TestName, EnableFlag) \
class ClassName : public IFeatureTest { \
public: \
    std::string getName() const override { return TestName; } \
    bool isEnabled(const TestConfig& config) const override { return config.EnableFlag; }

#endif // BIOGPU_TEST_FRAMEWORK_H
