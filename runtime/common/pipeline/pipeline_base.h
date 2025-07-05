// runtime/common/pipeline/pipeline_base.h
#ifndef PIPELINE_BASE_H
#define PIPELINE_BASE_H

#include <memory>
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include "../io/sequence_batch.h"
#include "../gpu/gpu_sequence_buffer.h"

namespace BioGPU {

// Base configuration for all pipelines
struct PipelineBaseConfig {
    // GPU settings
    size_t max_sequences_per_batch = 100000;
    size_t max_bases_per_batch = 100000000;  // 100MB
    int device_id = 0;
    
    // Memory settings
    bool use_pinned_memory = true;
    size_t gpu_memory_limit = 0;  // 0 = auto-detect
    
    // Performance settings
    bool enable_profiling = false;
    int num_threads = 1;
    
    // Output settings
    std::string output_prefix = "";
    bool verbose = false;
};

// Base class for pipeline-specific results
class PipelineResultsBase {
public:
    virtual ~PipelineResultsBase() = default;
    
    // Common metrics
    size_t sequences_processed = 0;
    size_t bases_processed = 0;
    double processing_time_seconds = 0;
    
    // Virtual methods for pipeline-specific reporting
    virtual void writeReport(const std::string& filename) const = 0;
    virtual void writeSummary(std::ostream& out) const = 0;
};

// Abstract base class for all processing pipelines
class PipelineBase {
protected:
    PipelineBaseConfig config;
    std::string pipeline_name;
    
    // GPU resources
    cudaStream_t main_stream = nullptr;
    std::vector<cudaStream_t> worker_streams;
    
    // Memory pools
    void* d_workspace = nullptr;
    size_t workspace_size = 0;
    
    // Profiling
    cudaEvent_t start_event = nullptr;
    cudaEvent_t end_event = nullptr;
    std::vector<float> batch_times;
    
    // State
    bool initialized = false;
    size_t total_sequences = 0;
    size_t total_bases = 0;
    
public:
    PipelineBase(const std::string& name, const PipelineBaseConfig& cfg)
        : pipeline_name(name), config(cfg) {}
    
    virtual ~PipelineBase() {
        cleanup();
    }
    
    // Initialization and cleanup
    virtual bool initialize() {
        if (initialized) return true;
        
        // Set GPU device
        if (cudaSetDevice(config.device_id) != cudaSuccess) {
            return false;
        }
        
        // Create main stream
        if (cudaStreamCreate(&main_stream) != cudaSuccess) {
            return false;
        }
        
        // Create worker streams if needed
        worker_streams.resize(config.num_threads);
        for (auto& stream : worker_streams) {
            if (cudaStreamCreate(&stream) != cudaSuccess) {
                cleanup();
                return false;
            }
        }
        
        // Create profiling events
        if (config.enable_profiling) {
            cudaEventCreate(&start_event);
            cudaEventCreate(&end_event);
        }
        
        // Initialize pipeline-specific resources
        if (!initializePipeline()) {
            cleanup();
            return false;
        }
        
        initialized = true;
        return true;
    }
    
    virtual void cleanup() {
        if (!initialized) return;
        
        // Clean up pipeline-specific resources
        cleanupPipeline();
        
        // Free workspace
        if (d_workspace) {
            cudaFree(d_workspace);
            d_workspace = nullptr;
        }
        
        // Destroy streams
        if (main_stream) {
            cudaStreamDestroy(main_stream);
            main_stream = nullptr;
        }
        
        for (auto& stream : worker_streams) {
            if (stream) {
                cudaStreamDestroy(stream);
            }
        }
        worker_streams.clear();
        
        // Destroy events
        if (start_event) {
            cudaEventDestroy(start_event);
            start_event = nullptr;
        }
        if (end_event) {
            cudaEventDestroy(end_event);
            end_event = nullptr;
        }
        
        initialized = false;
    }
    
    // Main processing interface
    bool processBatch(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream = nullptr) {
        if (!initialized || !gpu_buffer) return false;
        
        // Use provided stream or default to main stream
        if (!stream) stream = main_stream;
        
        // Start profiling
        float elapsed_ms = 0;
        if (config.enable_profiling && start_event && end_event) {
            cudaEventRecord(start_event, stream);
        }
        
        // Call pipeline-specific processing
        bool success = processBatchImpl(gpu_buffer, stream);
        
        if (success) {
            // Update statistics
            total_sequences += gpu_buffer->getCurrentSequences();
            total_bases += gpu_buffer->getCurrentTotalLength();
        }
        
        // End profiling
        if (config.enable_profiling && start_event && end_event) {
            cudaEventRecord(end_event, stream);
            cudaEventSynchronize(end_event);
            cudaEventElapsedTime(&elapsed_ms, start_event, end_event);
            batch_times.push_back(elapsed_ms);
        }
        
        return success;
    }
    
    // Get results
    virtual std::unique_ptr<PipelineResultsBase> getResults() const = 0;
    
    // Configuration access
    const PipelineBaseConfig& getConfig() const { return config; }
    void setConfig(const PipelineBaseConfig& cfg) { 
        if (!initialized) config = cfg; 
    }
    
    // Pipeline information
    const std::string& getName() const { return pipeline_name; }
    bool isInitialized() const { return initialized; }
    
    // Statistics
    size_t getTotalSequences() const { return total_sequences; }
    size_t getTotalBases() const { return total_bases; }
    
    // Profiling data
    std::vector<float> getBatchTimes() const { return batch_times; }
    float getAverageBatchTime() const {
        if (batch_times.empty()) return 0;
        float sum = 0;
        for (float t : batch_times) sum += t;
        return sum / batch_times.size();
    }
    
protected:
    // Pipeline-specific implementation (must be overridden)
    virtual bool initializePipeline() = 0;
    virtual void cleanupPipeline() = 0;
    virtual bool processBatchImpl(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) = 0;
    
    // Helper methods for derived classes
    bool allocateWorkspace(size_t size) {
        if (d_workspace && workspace_size >= size) {
            return true;  // Already have enough space
        }
        
        // Free existing workspace
        if (d_workspace) {
            cudaFree(d_workspace);
            d_workspace = nullptr;
        }
        
        // Allocate new workspace
        if (cudaMalloc(&d_workspace, size) != cudaSuccess) {
            workspace_size = 0;
            return false;
        }
        
        workspace_size = size;
        return true;
    }
    
    // Logging helper
    void log(const std::string& message) const {
        if (config.verbose) {
            std::cerr << "[" << pipeline_name << "] " << message << std::endl;
        }
    }
    
    void logError(const std::string& message) const {
        std::cerr << "[" << pipeline_name << " ERROR] " << message << std::endl;
    }
};

} // namespace BioGPU

#endif // PIPELINE_BASE_H