#pragma once

#include "../include/unified_pipeline_base.h"
#include "../../kernels/genes/amr_detection_pipeline_v2.h"
#include <memory>
#include <map>

namespace BioGPU {
namespace Unified {

class AMRGeneAdapter : public AnalysisModule {
public:
    AMRGeneAdapter();
    ~AMRGeneAdapter() override;
    
    bool initialize(const std::string& db_path) override;
    
    void processSequences(
        GPUSequenceBuffer* buffer,
        SharedKmerGenerator* kmer_gen,
        cudaStream_t stream) override;
    
    void finalizeSample(const std::string& sample_name) override;
    
    void generateReport(const std::string& output_path) override;
    
    std::string getModuleName() const override { return "AMR Gene Detection"; }
    
    // Configuration specific to AMR genes
    void setMinCoverage(float coverage) { min_coverage_ = coverage; }
    void setMinIdentity(float identity) { min_identity_ = identity; }

private:
    // Wrap the existing AMR detection pipeline
    std::unique_ptr<AMRDetectionPipelineV2> pipeline_;
    
    // Configuration
    float min_coverage_ = 0.8f;
    float min_identity_ = 0.9f;
    
    // Gene coverage tracking
    std::map<std::string, std::vector<uint32_t>> gene_coverage_;
    std::map<std::string, float> gene_rpkm_;
    
    // Accumulated results
    std::vector<AMRGeneResult> all_gene_results_;
    
    // Helper methods
    void calculateGeneCoverage(
        uint64_t* kmers,
        uint32_t kmer_count,
        cudaStream_t stream);
    
    void calculateGeneAbundance();
};

} // namespace Unified
} // namespace BioGPU