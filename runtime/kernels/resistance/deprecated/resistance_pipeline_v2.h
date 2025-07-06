// runtime/kernels/resistance/resistance_pipeline_v2.h
#ifndef RESISTANCE_PIPELINE_V2_H
#define RESISTANCE_PIPELINE_V2_H

#include "../../common/pipeline/pipeline_base.h"
#include "../../common/config/unified_config.h"
#include "../../common/gpu/gpu_sequence_buffer.h"
#include "../../common/gpu/bloom_filter.h"
#include "global_fq_resistance_mapper.h"
#include "fq_mutations_hardcoded.h"
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <atomic>

namespace BioGPU {

// Forward declarations
struct ResistanceMutation;
struct QRDRRegion;

// Resistance-specific results
class ResistanceResults : public PipelineResultsBase {
public:
    struct ResistanceHit {
        uint32_t read_id;
        std::string sample_name;
        std::string gene_name;
        std::string species;
        int position;
        char wildtype_aa;
        char mutant_aa;
        float confidence;
        bool is_qrdr;
        int read_support;
        std::string mutation_name;  // e.g., "GyrA_S83L"
    };
    
    struct AlleleFrequency {
        std::string gene_name;
        std::string species;
        int position;
        std::map<char, float> amino_acid_frequencies;
        char wildtype_aa;
        float wildtype_frequency;
        float total_coverage;
    };
    
    struct SampleResistanceProfile {
        std::string sample_name;
        std::vector<ResistanceHit> mutations;
        std::map<std::string, std::vector<AlleleFrequency>> allele_frequencies;  // gene -> frequencies
        bool has_fluoroquinolone_resistance;
        std::vector<std::string> resistance_mechanisms;
    };
    
    // Results storage
    std::map<std::string, SampleResistanceProfile> sample_profiles;
    std::vector<ResistanceHit> all_hits;
    
    // Override virtual methods
    void writeReport(const std::string& filename) const override;
    void writeSummary(std::ostream& out) const override;
    
    // Resistance-specific reporting
    void writeClinicalReport(const std::string& filename, const std::string& format = "html") const;
    void writeAlleleFrequencyReport(const std::string& filename) const;
    void writeJSONReport(const std::string& filename) const;
};

// Stage-specific match structures
struct BloomFilterMatch {
    uint32_t read_id;
    bool passed_r1;
    bool passed_r2;
};

struct NucleotideMatch {
    uint32_t read_id;
    uint32_t gene_id;
    int position;
    int mismatches;
    bool is_reverse;
};

struct ProteinMatch {
    uint32_t read_id;
    uint32_t gene_id;
    int protein_position;
    char amino_acid;
    float score;
    bool is_complete_codon;
    int frame;  // 1-3 for forward, -1 to -3 for reverse
};

// Allele frequency tracking
class AlleleFrequencyAnalyzer {
private:
    std::map<std::string, std::map<int, std::map<char, std::atomic<int>>>> position_counts;
    std::map<std::string, std::map<int, std::atomic<int>>> total_counts;
    mutable std::mutex mutex;
    
public:
    void addObservation(const std::string& gene, int position, char amino_acid);
    void addMultipleObservations(const std::string& gene, int position, char amino_acid, int count);
    ResistanceResults::AlleleFrequency calculateFrequency(const std::string& gene, int position, char wildtype_aa) const;
    std::vector<ResistanceResults::AlleleFrequency> getAllFrequencies(const std::string& gene) const;
    void clear();
};

// Fluoroquinolone Resistance Detection Pipeline - refactored to use shared components
class ResistancePipeline : public PipelineBase {
private:
    // Configuration from unified config
    const ResistanceConfig& resistance_config;
    
    // Batch processing parameters (maintaining compatibility)
    static constexpr int BATCH_SIZE = 10000;
    static constexpr int MAX_READ_LENGTH = 300;
    static constexpr int KMER_LENGTH = 15;
    
    // FQ Mutation Database
    std::unique_ptr<GlobalFQResistanceMapper> fq_mapper;
    std::vector<ResistanceMutation> mutation_database;
    std::map<std::string, QRDRRegion> qrdr_regions;
    
    // GPU memory for three-stage pipeline
    // Stage 1: Bloom filter (using shared component)
    std::unique_ptr<BloomFilter> bloom_filter;
    BloomFilterMatch* d_bloom_matches = nullptr;
    
    // Stage 2: Nucleotide k-mer matching
    uint64_t* d_kmer_table = nullptr;
    NucleotideMatch* d_nucleotide_matches = nullptr;
    uint32_t* d_nucleotide_match_counts = nullptr;
    
    // Stage 3: Translated protein search
    ProteinMatch* d_protein_matches = nullptr;
    uint32_t* d_protein_match_counts = nullptr;
    char* d_mutation_sequences = nullptr;  // Reference sequences for mutations
    
    // Result aggregation
    AlleleFrequencyAnalyzer allele_analyzer;
    std::unique_ptr<ResistanceResults> accumulated_results;
    
    // Current sample being processed
    std::string current_sample_name;
    
    // Performance metrics
    size_t stage1_passed = 0;
    size_t stage2_passed = 0;
    size_t stage3_passed = 0;
    
public:
    ResistancePipeline(const ResistanceConfig& config);
    ~ResistancePipeline() = default;
    
    // Set current sample name for batch processing
    void setCurrentSample(const std::string& sample_name) { 
        current_sample_name = sample_name; 
    }
    
    // Get accumulated results
    std::unique_ptr<PipelineResultsBase> getResults() const override {
        auto results = std::make_unique<ResistanceResults>(*accumulated_results);
        results->sequences_processed = getTotalSequences();
        results->bases_processed = getTotalBases();
        results->processing_time_seconds = getAverageBatchTime() * getBatchTimes().size() / 1000.0;
        return results;
    }
    
    // Additional methods for compatibility
    void calculateAlleleFrequencies();
    void generateClinicalReport(const std::string& output_path, const std::string& format = "html");
    
protected:
    // Implement pure virtual methods from PipelineBase
    bool initializePipeline() override;
    void cleanupPipeline() override;
    bool processBatchImpl(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) override;
    
private:
    // Three-stage pipeline implementation
    void stage1_bloomFilterScreening(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream);
    void stage2_nucleotideKmerMatching(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream);
    void stage3_translatedProteinSearch(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream);
    
    // Database initialization
    bool loadFQMutationDatabase();
    bool buildBloomFilterFromMutations();
    bool buildKmerTable();
    bool prepareProteinReferences();
    
    // QRDR detection logic
    void initializeQRDRRegions();
    bool isInQRDR(const std::string& gene, int position) const;
    
    // GPU kernel launchers
    void launchBloomFilterKernel(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        int num_sequences,
        BloomFilterMatch* d_matches,
        cudaStream_t stream);
    
    void launchNucleotideMatchingKernel(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        const BloomFilterMatch* d_bloom_matches,
        int num_sequences,
        const uint64_t* d_kmer_table,
        NucleotideMatch* d_matches,
        uint32_t* d_match_counts,
        cudaStream_t stream);
    
    void launchProteinTranslationKernel(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        const NucleotideMatch* d_nucleotide_matches,
        const uint32_t* d_match_counts,
        int num_sequences,
        const char* d_mutation_sequences,
        ProteinMatch* d_protein_matches,
        uint32_t* d_protein_match_counts,
        cudaStream_t stream);
    
    // Result processing
    void processProteinMatches(cudaStream_t stream);
    void aggregateMutationResults(const std::vector<ProteinMatch>& matches);
    void updateAlleleFrequencies(const std::vector<ProteinMatch>& matches);
    
    // Species-specific processing
    std::string normalizeSpeciesName(const std::string& species) const;
    std::string normalizeGeneName(const std::string& gene, const std::string& species) const;
    
    // Helper methods
    size_t estimateMemoryRequirements(size_t num_sequences);
    void clearBatchResults();
    ResistanceMutation getMutationInfo(uint32_t gene_id) const;
};

// QRDR Region definitions
struct QRDRRegion {
    std::string gene;
    std::string species;
    int start_position;
    int end_position;
    std::vector<int> critical_positions;  // Positions known to confer resistance
};

// Extended mutation structure for the pipeline
struct ResistanceMutation {
    uint32_t id;
    std::string gene_name;
    std::string species;
    int position;
    char wildtype_aa;
    char mutant_aa;
    std::string mutation_name;
    std::string drug_class;
    float resistance_level;  // 1.0 = high, 0.5 = moderate, 0.25 = low
    bool is_primary_mutation;
    std::string nucleotide_sequence;  // Reference sequence around mutation
    std::string protein_context;      // Protein sequence context
};

} // namespace BioGPU

#endif // RESISTANCE_PIPELINE_V2_H