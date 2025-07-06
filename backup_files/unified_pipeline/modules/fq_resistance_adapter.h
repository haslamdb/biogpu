#pragma once

#include "../include/unified_pipeline_base.h"
#include "../../kernels/resistance/resistance_pipeline_v2.h"
#include <memory>

namespace BioGPU {
namespace Unified {

class FQResistanceAdapter : public AnalysisModule {
public:
    FQResistanceAdapter();
    ~FQResistanceAdapter() override;
    
    bool initialize(const std::string& db_path) override;
    
    void processSequences(
        GPUSequenceBuffer* buffer,
        SharedKmerGenerator* kmer_gen,
        cudaStream_t stream) override;
    
    void finalizeSample(const std::string& sample_name) override;
    
    void generateReport(const std::string& output_path) override;
    
    std::string getModuleName() const override { return "FQ Resistance Detection"; }
    
    // Configuration specific to FQ resistance
    void setMinAlleleDepth(uint32_t depth) { min_allele_depth_ = depth; }
    void setUseBloomFilter(bool use) { use_bloom_filter_ = use; }
    void setUseSmithWaterman(bool use) { use_smith_waterman_ = use; }

private:
    // Wrap the existing resistance pipeline
    std::unique_ptr<ResistancePipelineV2> pipeline_;
    
    // Configuration
    uint32_t min_allele_depth_ = 5;
    bool use_bloom_filter_ = true;
    bool use_smith_waterman_ = true;
    
    // Accumulated results across samples
    std::vector<SampleResistanceResult> all_results_;
    
    // Helper methods
    void convertKmersToResistanceFormat(
        uint64_t* unified_kmers,
        uint32_t* unified_counts,
        cudaStream_t stream);
};

} // namespace Unified
} // namespace BioGPU