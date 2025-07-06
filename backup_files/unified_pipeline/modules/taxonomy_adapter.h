#pragma once

#include "../include/unified_pipeline_base.h"
#include <memory>
#include <map>
#include <string>

// Forward declarations for Kraken-style components
struct KrakenDatabase;
struct TaxonomyNode;

namespace BioGPU {
namespace Unified {

class TaxonomyAdapter : public AnalysisModule {
public:
    TaxonomyAdapter();
    ~TaxonomyAdapter() override;
    
    bool initialize(const std::string& db_path) override;
    
    void processSequences(
        GPUSequenceBuffer* buffer,
        SharedKmerGenerator* kmer_gen,
        cudaStream_t stream) override;
    
    void finalizeSample(const std::string& sample_name) override;
    
    void generateReport(const std::string& output_path) override;
    
    std::string getModuleName() const override { return "Taxonomic Classification"; }
    
    // Configuration specific to taxonomy
    void setConfidenceThreshold(float threshold) { confidence_threshold_ = threshold; }
    void setMinHitGroups(int groups) { min_hit_groups_ = groups; }

private:
    // Kraken-style database
    std::unique_ptr<KrakenDatabase> kraken_db_;
    
    // Taxonomy tree
    std::map<uint32_t, TaxonomyNode> taxonomy_tree_;
    std::map<uint32_t, std::string> taxid_to_name_;
    
    // Classification results
    std::map<uint32_t, uint64_t> taxid_counts_;
    std::map<uint32_t, float> taxid_abundance_;
    
    // Configuration
    float confidence_threshold_ = 0.0f;
    int min_hit_groups_ = 2;
    
    // GPU memory for classification
    uint32_t* d_taxids_;
    float* d_scores_;
    
    // Helper methods
    void classifyMinimizers(
        uint64_t* minimizers,
        uint32_t minimizer_count,
        cudaStream_t stream);
    
    void resolveClassifications();
    void calculateAbundance();
    
    // Report generation
    void generateKrakenReport(const std::string& output_path);
    void generateBrackenReport(const std::string& output_path);
};

} // namespace Unified
} // namespace BioGPU