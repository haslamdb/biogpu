// runtime/kernels/genes/amr_detection_pipeline_v2.h
#ifndef AMR_DETECTION_PIPELINE_V2_H
#define AMR_DETECTION_PIPELINE_V2_H

#include "../../common/pipeline/pipeline_base.h"
#include "../../common/config/unified_config.h"
#include "../../common/gpu/gpu_sequence_buffer.h"
#include "ncbi_amr_database_loader.h"
#include <memory>
#include <vector>
#include <string>

namespace BioGPU {

// AMR-specific results
class AMRResults : public PipelineResultsBase {
public:
    struct AMRHit {
        uint32_t read_id;
        uint32_t gene_id;
        uint16_t ref_start;
        uint16_t ref_end;
        uint16_t read_start;
        uint16_t read_end;
        float identity;
        float coverage;
        int8_t frame;
        uint8_t num_mutations;
        bool is_complete_gene;
        std::string gene_name;
        std::string drug_class;
    };
    
    struct GeneCoverage {
        std::string gene_name;
        uint32_t gene_id;
        uint32_t total_reads;
        float mean_coverage;
        float coverage_uniformity;
        uint16_t covered_positions;
        uint16_t gene_length;
    };
    
    std::vector<AMRHit> hits;
    std::vector<GeneCoverage> gene_coverage;
    
    // High-confidence hits (>95% identity and coverage)
    std::vector<AMRHit> high_confidence_hits;
    
    // Override virtual methods
    void writeReport(const std::string& filename) const override;
    void writeSummary(std::ostream& out) const override;
    
    // AMR-specific reporting
    void writeClinicalReport(const std::string& filename) const;
    void writeDrugResistanceProfile(const std::string& filename) const;
};

// Structure to hold minimizer information
struct Minimizer {
    uint64_t hash;
    uint32_t pos;
    bool is_reverse;
};

// AMR Detection Pipeline - refactored to use shared components
class AMRDetectionPipeline : public PipelineBase {
private:
    // Configuration from unified config
    const GenesConfig& genes_config;
    
    // AMR Database
    std::unique_ptr<NCBIAMRDatabaseLoader> amr_db;
    
    // GPU memory for AMR-specific data
    Minimizer* d_minimizers = nullptr;
    uint32_t* d_minimizer_counts = nullptr;
    uint32_t* d_minimizer_offsets = nullptr;
    uint64_t* d_bloom_filter = nullptr;
    
    // GPU memory for results
    void* d_amr_hits = nullptr;  // AMRResults::AMRHit structs
    uint32_t* d_hit_counts = nullptr;
    uint32_t* d_coverage_stats = nullptr;
    
    // Bloom filter parameters
    size_t bloom_filter_size = 1ULL << 30;  // 1GB
    int bloom_k = 3;
    
    // Batch state
    size_t max_minimizers_per_batch = 0;
    size_t max_hits_per_batch = 0;
    
    // Results accumulator
    std::unique_ptr<AMRResults> accumulated_results;
    
public:
    AMRDetectionPipeline(const GenesConfig& config);
    ~AMRDetectionPipeline() = default;
    
    // Set database path
    void setDatabasePath(const std::string& nucl_path, const std::string& prot_path);
    
    // Get accumulated results
    std::unique_ptr<PipelineResultsBase> getResults() const override {
        auto results = std::make_unique<AMRResults>(*accumulated_results);
        results->sequences_processed = getTotalSequences();
        results->bases_processed = getTotalBases();
        results->processing_time_seconds = getAverageBatchTime() * getBatchTimes().size() / 1000.0;
        return results;
    }
    
protected:
    // Implement pure virtual methods from PipelineBase
    bool initializePipeline() override;
    void cleanupPipeline() override;
    bool processBatchImpl(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) override;
    
private:
    // AMR-specific processing steps
    bool loadAMRDatabase();
    bool buildBloomFilter();
    
    // GPU kernels wrapped for this pipeline
    void launchGenerateMinimizersKernel(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        int num_sequences,
        Minimizer* d_minimizers,
        uint32_t* d_minimizer_counts,
        uint32_t* d_minimizer_offsets,
        cudaStream_t stream);
    
    void launchScreenWithBloomFilterKernel(
        const Minimizer* d_minimizers,
        const uint32_t* d_minimizer_counts,
        const uint32_t* d_minimizer_offsets,
        const uint64_t* d_bloom_filter,
        int num_sequences,
        uint8_t* d_passed_filter,
        cudaStream_t stream);
    
    void launchTranslatedSearchKernel(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        const uint8_t* d_passed_filter,
        int num_sequences,
        const char* d_protein_db,
        const int* d_protein_offsets,
        int num_proteins,
        void* d_hits,
        uint32_t* d_hit_counts,
        cudaStream_t stream);
    
    // Result processing
    void processHits(cudaStream_t stream);
    void updateCoverageStats(cudaStream_t stream);
    
    // Helper methods
    size_t estimateMemoryRequirements(size_t num_sequences, size_t total_bases);
    void clearBatchResults();
};

} // namespace BioGPU

#endif // AMR_DETECTION_PIPELINE_V2_H