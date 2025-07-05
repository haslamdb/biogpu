// runtime/common/pipeline/pipeline_coordinator.h
#ifndef PIPELINE_COORDINATOR_H
#define PIPELINE_COORDINATOR_H

#include <memory>
#include <vector>
#include <string>
#include <atomic>
#include <thread>
#include <cuda_runtime.h>
#include "../io/streaming_fastq_reader.h"
#include "../io/sequence_batch.h"
#include "../gpu/gpu_sequence_buffer.h"

namespace BioGPU {

// Forward declarations for pipeline-specific components
class ResistancePipeline;
class GenesPipeline;
class ProfilerPipeline;

struct PipelineConfig {
    // General settings
    size_t batch_size = 100000;
    size_t max_gpu_sequences = 500000;
    size_t max_gpu_bases = 500000000;  // 500MB
    int num_gpu_streams = 2;
    
    // Pipeline-specific flags
    bool enable_resistance = true;
    bool enable_genes = true;
    bool enable_profiler = true;
    
    // Output settings
    std::string output_dir = "output";
    bool write_intermediate_results = false;
    
    // Performance settings
    bool use_pinned_memory = true;
    bool overlap_compute_transfer = true;
    int num_worker_threads = 4;
};

struct PipelineResults {
    // Timing information
    double total_time_seconds = 0;
    double resistance_time_seconds = 0;
    double genes_time_seconds = 0;
    double profiler_time_seconds = 0;
    
    // Statistics
    size_t total_reads_processed = 0;
    size_t total_bases_processed = 0;
    size_t batches_processed = 0;
    
    // Pipeline-specific results (to be populated by individual pipelines)
    void* resistance_results = nullptr;
    void* genes_results = nullptr;
    void* profiler_results = nullptr;
};

class PipelineCoordinator {
private:
    PipelineConfig config;
    
    // Input reader
    std::unique_ptr<StreamingFastqReader> reader;
    
    // GPU resources
    std::vector<cudaStream_t> gpu_streams;
    std::vector<std::unique_ptr<GPUSequenceBuffer>> gpu_buffers;
    int current_buffer_idx = 0;
    
    // Pipeline instances
    std::unique_ptr<ResistancePipeline> resistance_pipeline;
    std::unique_ptr<GenesPipeline> genes_pipeline;
    std::unique_ptr<ProfilerPipeline> profiler_pipeline;
    
    // Processing state
    std::atomic<bool> processing_active{false};
    std::atomic<size_t> reads_processed{0};
    std::atomic<size_t> bases_processed{0};
    
    // Results
    PipelineResults results;
    
public:
    PipelineCoordinator(const PipelineConfig& cfg = PipelineConfig());
    ~PipelineCoordinator();
    
    // Initialize pipelines and GPU resources
    bool initialize();
    void cleanup();
    
    // Process single-end or paired-end FASTQ files
    bool processFastq(const std::string& fastq_file);
    bool processPairedFastq(const std::string& r1_file, const std::string& r2_file);
    
    // Process a list of samples from CSV
    bool processSampleList(const std::string& csv_file);
    
    // Get results
    const PipelineResults& getResults() const { return results; }
    
    // Configuration
    void setConfig(const PipelineConfig& cfg) { config = cfg; }
    const PipelineConfig& getConfig() const { return config; }
    
private:
    // Internal processing methods
    bool initializeGPUResources();
    void cleanupGPUResources();
    
    bool processBatch(std::shared_ptr<SequenceBatch> batch, int stream_idx);
    void processAllBatches();
    
    // Pipeline coordination
    bool runResistancePipeline(GPUSequenceBuffer* buffer, cudaStream_t stream);
    bool runGenesPipeline(GPUSequenceBuffer* buffer, cudaStream_t stream);
    bool runProfilerPipeline(GPUSequenceBuffer* buffer, cudaStream_t stream);
    
    // Utility methods
    void logProgress(size_t current_reads, size_t current_bases);
    void finalizeResults();
};

} // namespace BioGPU

#endif // PIPELINE_COORDINATOR_H