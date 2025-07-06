#pragma once

#include <memory>
#include <vector>
#include <string>
#include <cuda_runtime.h>

namespace BioGPU {
namespace Unified {

class GPUSequenceBuffer;
class SharedKmerGenerator;

class AnalysisModule {
public:
    virtual ~AnalysisModule() = default;
    
    virtual bool initialize(const std::string& db_path) = 0;
    
    virtual void processSequences(
        GPUSequenceBuffer* buffer, 
        SharedKmerGenerator* kmer_gen, 
        cudaStream_t stream) = 0;
    
    virtual void finalizeSample(const std::string& sample_name) = 0;
    
    virtual void generateReport(const std::string& output_path) = 0;
    
    virtual std::string getModuleName() const = 0;
    
    virtual bool isEnabled() const { return enabled_; }
    virtual void setEnabled(bool enabled) { enabled_ = enabled; }

protected:
    bool enabled_ = true;
};

struct PipelineConfig {
    bool enable_fq_resistance = true;
    bool enable_amr_genes = true;
    bool enable_taxonomy = true;
    int batch_size = 100000;
    std::string output_dir = "results";
    bool save_intermediate = false;
    bool verbose = false;
    
    int kmer_size_nucleotide = 31;
    int kmer_size_protein = 15;
    int minimizer_k = 15;
    int minimizer_w = 10;
    
    size_t gpu_memory_limit = 0;  // 0 = auto-detect
};

} // namespace Unified
} // namespace BioGPU