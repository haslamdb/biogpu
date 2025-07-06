#pragma once

#include "unified_pipeline_base.h"
#include "shared_kmer_generator.h"
#include "../common/io/sample_csv_parser.h"
#include "../common/gpu/gpu_sequence_buffer.h"

#include <memory>
#include <vector>
#include <chrono>

namespace BioGPU {
namespace Unified {

class UnifiedBioGPUPipeline {
public:
    UnifiedBioGPUPipeline();
    ~UnifiedBioGPUPipeline();
    
    // Initialize the pipeline with database paths
    bool initialize(const std::string& db_base_path);
    
    // Process a single sample
    void processSample(const SampleInfo& sample);
    
    // Process a batch of samples from CSV
    void processBatch(const std::string& csv_path);
    
    // Configuration
    void setConfig(const PipelineConfig& config) { config_ = config; }
    PipelineConfig& getConfig() { return config_; }
    
    // Module management
    void addModule(std::unique_ptr<AnalysisModule> module);
    void enableModule(const std::string& module_name, bool enable);
    
    // Report generation
    void generateCombinedReports();
    void generateMasterSummary();

private:
    // Configuration
    PipelineConfig config_;
    
    // Modules
    std::vector<std::unique_ptr<AnalysisModule>> modules_;
    
    // Shared components
    std::unique_ptr<SharedKmerGenerator> kmer_generator_;
    std::unique_ptr<GPUSequenceBuffer> gpu_buffer_;
    
    // CUDA resources
    cudaStream_t stream_io_;
    cudaStream_t stream_compute_;
    
    // Statistics
    struct PipelineStats {
        size_t total_samples = 0;
        size_t total_reads = 0;
        size_t total_bases = 0;
        std::chrono::duration<double> total_time{0};
        std::map<std::string, std::chrono::duration<double>> module_times;
    } stats_;
    
    // Helper methods
    void createOutputDirectories(const std::string& sample_name);
    void logProgress(const std::string& message);
    std::string getCurrentTimestamp();
    
    // Processing helpers
    void processBatchWithModules(
        SequenceBatch* batch,
        const std::string& sample_name);
    
    // Report helpers
    void writePipelineMetrics();
};

// Factory function for standard configuration
std::unique_ptr<UnifiedBioGPUPipeline> createStandardPipeline(
    bool enable_fq_resistance = true,
    bool enable_amr_genes = true,
    bool enable_taxonomy = true);

} // namespace Unified
} // namespace BioGPU