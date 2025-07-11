// amr_detection_pipeline.h
#ifndef AMR_DETECTION_PIPELINE_H
#define AMR_DETECTION_PIPELINE_H

#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include "ncbi_amr_database_loader.h"
#include "sample_csv_parser.h"  // Reuse from FQ pipeline
// #include "bloom_filter.h"       // Reuse bloom filter

// Structure to track paired-end relationships
struct PairedReadInfo {
    uint32_t read_idx;      // Original read pair index
    bool is_read2;          // false for R1, true for R2
    uint32_t pair_offset;   // Offset to find the paired read results
};

// Structure to hold paired-end read data
struct ReadPairData {
    std::string read1_seq;
    std::string read2_seq;
    std::string read1_id;
    std::string read2_id;
    uint32_t pair_index;    // Original index in the batch
};

// Configuration for AMR detection
struct AMRDetectionConfig {
    // Bloom filter parameters
    bool use_bloom_filter = false;          // Disabled by default
    size_t bloom_filter_size = 1ULL << 30;  // 1GB bloom filter
    int bloom_k = 3;                        // Number of hash functions
    int kmer_length = 31;                   // DNA k-mer size for bloom filter
    
    // Minimizer parameters
    int minimizer_k = 15;                   // DNA minimizer k-mer size
    int minimizer_w = 10;                   // Window size
    
    // Protein search parameters
    int protein_kmer_size = 8;              // For protein alignment after translation
    
    // Alignment parameters (increased stringency)
    float min_identity = 0.95f;             // Minimum identity for AMR match (increased from 0.90f)
    float min_coverage = 0.90f;             // Minimum coverage of AMR gene (increased from 0.80f)
    int min_alignment_length = 75;          // Minimum alignment length (aa) (increased from 50)
    int band_width = 15;                    // Band width for banded SW
    
    // Paired-end concordance parameters
    float concordance_bonus = 2.0f;         // Bonus for concordant pairs
    float discord_penalty = 0.5f;           // Penalty for discordant pairs
    int expected_fragment_length = 500;     // Expected fragment length
    int fragment_length_stddev = 100;       // Standard deviation for fragment length
    float max_fragment_length_zscore = 3.0f; // Max z-score for valid fragments
    
    // Batch processing (for paired-end, this is per pair, so actual reads = 2x)
    int reads_per_batch = 100000;
    int max_read_length = 300;
    
    // Output options
    bool output_sam = false;
    bool output_hdf5 = true;
    std::string output_prefix = "amr_results";
    
    // Database paths
    std::string protein_db_path = "amr_protein_db";  // Default location
    
    // EM algorithm parameters
    bool use_em = false;               // Enable EM algorithm for multi-mapping reads
    int em_iterations = 10;            // Number of EM iterations
    float em_convergence = 0.001f;     // Convergence threshold
    float min_hit_coverage = -1.0f;    // Minimum coverage for high quality hits (-1 = no threshold)
};

// Structure to hold minimizer information
struct Minimizer {
    uint64_t hash;
    uint32_t pos;      // Position in read
    bool is_reverse;   // Strand
};

// AMR hit information
struct AMRHit {
    uint32_t read_id;
    uint32_t gene_id;
    uint16_t ref_start;     // Start position on reference protein
    uint16_t ref_end;       // End position on reference protein
    uint16_t read_start;    // Start position on read (translated)
    uint16_t read_end;      // End position on read (translated)
    float identity;
    float coverage;         // Fraction of AMR gene covered
    int8_t frame;          // Translation frame (-3 to +3, excluding 0)
    bool is_complete_gene;
    bool concordant;        // For paired-end reads, true if both reads hit same gene
    char gene_name[64];
    char drug_class[32];
    char gene_family[32];
    // Enhanced paired-end tracking
    uint32_t pair_id;       // Original read pair ID
    bool is_read2;          // false for R1, true for R2
    uint32_t mate_read_id;  // ID of the mate read
    float pair_score;       // Combined score if concordant
};

// Coverage tracking per AMR gene
struct AMRCoverageStats {
    uint32_t total_reads;           // Number of reads mapped to this gene
    uint32_t total_bases_mapped;    // Total bases from all reads mapped
    uint32_t covered_positions;     // Number of positions with at least 1x coverage
    uint32_t gene_length;           // Length of the gene (in amino acids for proteins)
    float percent_coverage;         // Percentage of gene covered (0-100)
    float mean_depth;              // Average depth across covered positions
    float mean_coverage;           // Average coverage across all positions
    float coverage_uniformity;     // Measure of how evenly coverage is distributed
    float rpkm;                    // Reads Per Kilobase per Million mapped reads
    float tpm;                     // Transcripts Per Million
    
    // Position-specific coverage tracking
    // Note: We'll allocate this dynamically based on gene length
    uint32_t* position_counts;     // Coverage depth at each position (allocated separately)
};

// Gene abundance information for clinical reporting
struct GeneAbundance {
    uint32_t gene_id;
    char gene_name[64];
    char drug_class[32];
    uint32_t read_count;
    float rpkm;
    float tpm;
    float coverage_depth;
    float coverage_breadth;  // Percentage of gene covered
};

// Structure to track alignments for a read pair
struct ReadPairAlignment {
    uint32_t pair_id;
    uint32_t r1_read_idx;
    uint32_t r2_read_idx;
    
    // R1 alignments
    std::vector<uint32_t> r1_gene_ids;
    std::vector<float> r1_scores;
    std::vector<int32_t> r1_positions;  // Position on reference
    std::vector<int8_t> r1_frames;      // Translation frame
    
    // R2 alignments  
    std::vector<uint32_t> r2_gene_ids;
    std::vector<float> r2_scores;
    std::vector<int32_t> r2_positions;  // Position on reference
    std::vector<int8_t> r2_frames;      // Translation frame
    
    // Concordant gene matches
    std::vector<uint32_t> concordant_genes;
    std::vector<float> concordant_scores;
    std::vector<int32_t> fragment_lengths;  // Fragment length for each concordant match
    std::vector<float> fragment_probabilities; // Probability based on fragment length distribution
    
    // Best concordant match
    int32_t best_fragment_length = -1;
    uint32_t best_gene_id = UINT32_MAX;
    float best_concordant_score = 0.0f;
    
    bool hasAlignment() const { 
        return !r1_gene_ids.empty() || !r2_gene_ids.empty(); 
    }
    
    bool isConcordant() const {
        return !concordant_genes.empty();
    }
    
    // Calculate fragment length from R1 and R2 positions
    int32_t calculateFragmentLength(int32_t r1_pos, int32_t r2_pos, 
                                   int32_t gene_length, int r1_len, int r2_len) const {
        // For protein alignments, positions are in amino acids
        // Convert to nucleotide positions and calculate fragment length
        int32_t r1_nt_pos = r1_pos * 3;
        int32_t r2_nt_pos = r2_pos * 3;
        int32_t gene_nt_length = gene_length * 3;
        
        // R2 should be downstream of R1 for proper pairs
        if (r2_nt_pos > r1_nt_pos) {
            // Fragment spans from R1 start to R2 end
            return (r2_nt_pos + r2_len) - r1_nt_pos;
        }
        return -1; // Invalid fragment
    }
};

class AMRDetectionPipeline {
public:
    // Enhanced EM data structures for Kallisto-style gene assignment
    struct ReadAssignment {
        uint32_t read_id;
        std::vector<uint32_t> candidate_genes;
        std::vector<float> alignment_scores;
        std::vector<float> assignment_probabilities;
        float total_score;
        bool high_quality;  // Filter for Smith-Waterman score >= 100
    };
    
    // Pair-centric assignment structure for paired-end EM
    struct PairedReadAssignment {
        uint32_t pair_id;
        uint32_t r1_read_id;
        uint32_t r2_read_id;
        
        // Compatible genes (where both reads can map consistently)
        std::vector<uint32_t> compatible_genes;
        
        // Individual hit information
        std::map<uint32_t, float> r1_gene_scores;  // gene_id -> score
        std::map<uint32_t, float> r2_gene_scores;  // gene_id -> score
        
        // Paired information
        std::map<uint32_t, float> concordance_scores;     // gene_id -> combined score
        std::map<uint32_t, int32_t> fragment_lengths;     // gene_id -> fragment length
        std::map<uint32_t, float> fragment_probabilities; // gene_id -> P(fragment|gene)
        
        // Assignment probabilities for each gene
        std::map<uint32_t, float> assignment_probabilities;
        
        // Total probability normalization factor
        float total_probability;
        
        bool hasR1Hit(uint32_t gene_id) const { return r1_gene_scores.count(gene_id) > 0; }
        bool hasR2Hit(uint32_t gene_id) const { return r2_gene_scores.count(gene_id) > 0; }
        bool isConcordant(uint32_t gene_id) const { return hasR1Hit(gene_id) && hasR2Hit(gene_id); }
        
        float getCombinedScore(uint32_t gene_id) const {
            float score = 0.0f;
            if (hasR1Hit(gene_id)) score += r1_gene_scores.at(gene_id);
            if (hasR2Hit(gene_id)) score += r2_gene_scores.at(gene_id);
            if (isConcordant(gene_id) && concordance_scores.count(gene_id)) {
                score = concordance_scores.at(gene_id);
            }
            return score;
        }
    };

    struct GeneAbundanceInfo {
        uint32_t gene_id;
        std::string gene_name;
        std::string gene_family;
        float effective_length;     // Length adjusted for mappability
        float abundance;           // Current abundance estimate (RPKM-like)
        float total_reads;         // Expected number of reads assigned
        float prior_weight;        // Prior probability based on known prevalence
    };

private:
    static constexpr size_t MAX_MATCHES_PER_READ = 32;  // Must match what's in the kernel
    
    // Configuration
    AMRDetectionConfig config;
    
    // Database
    std::unique_ptr<NCBIAMRDatabaseLoader> amr_db;
    
    // GPU memory for reads
    char* d_reads;
    int* d_read_offsets;
    int* d_read_lengths;
    uint32_t* d_read_ids;
    
    // GPU memory for minimizers
    Minimizer* d_minimizers;
    uint32_t* d_minimizer_counts;
    uint32_t* d_minimizer_offsets;
    
    // GPU memory for bloom filter
    uint64_t* d_bloom_filter;
    bool* d_read_passes_filter;  // Boolean array for bloom filter results
    
    // GPU memory for alignment results
    AMRHit* d_amr_hits;
    uint32_t* d_hit_counts;
    
    // Coverage statistics
    AMRCoverageStats* d_coverage_stats;
    std::vector<AMRCoverageStats> h_coverage_stats;
    
    // Batch processing
    int current_batch_size;
    uint64_t total_reads_processed;
    
    // Performance tracking
    std::chrono::steady_clock::time_point processing_start_time;
    uint64_t reads_processed_checkpoint;
    std::chrono::steady_clock::time_point last_performance_report;
    
    // Gene entries from database
    std::vector<AMRGeneEntry> gene_entries;
    
    // Translated search engine (reused across batches)
    void* translated_search_engine;
    bool search_engine_initialized;
    int engine_capacity;  // Maximum batch size for search engine
    
    // Paired-end read tracking
    std::vector<PairedReadInfo> paired_read_info;
    
    // Enhanced paired-end alignment tracking
    std::vector<ReadPairAlignment> read_pair_alignments;
    std::vector<ReadPairAlignment> current_pair_alignments;
    
public:
    AMRDetectionPipeline(const AMRDetectionConfig& cfg);
    ~AMRDetectionPipeline();
    
    // Initialize pipeline
    bool initialize(const std::string& amr_db_path);
    
    // Validation
    void validateGeneMappings();
    
    // Process a batch of reads
    void processBatch(const std::vector<std::string>& reads,
                     const std::vector<std::string>& read_ids);
    
    // Process a batch of paired-end reads (legacy - merges reads)
    void processBatchPaired(const std::vector<std::string>& reads1,
                           const std::vector<std::string>& reads2,
                           const std::vector<std::string>& read_ids);
    
    // Process a batch of paired-end reads (new - keeps R1/R2 separate)
    void processPairedBatch(const std::vector<ReadPairData>& read_pairs);
    
    // Main processing steps
    void buildBloomFilter();
    void generateMinimizers();
    void screenWithBloomFilter();
    void performTranslatedAlignment();
    void extendAlignments();
    void calculateCoverageStats();
    void updateCoverageStatsFromHits();  // Lightweight per-batch update
    void finalizeCoverageStats();        // Heavy calculation after all batches
    void calculateAbundanceMetrics();
    
    // EM algorithm methods
    void resolveAmbiguousAssignmentsEM();  // Main entry point for EM
    void runKallistoStyleEM();
    void buildGeneFamiliesMap();
    void buildSequenceSimilarityMatrix();
    void initializeGeneAbundances();
    void normalizeInitialAbundances();
    void updateAssignmentProbabilities();
    void updateGeneAbundances();
    void applyFamilyConstraints();
    
    // Coverage stats management
    void initializeCoverageStats();
    void freeCoverageStats();
    
    // Get results
    std::vector<AMRHit> getAMRHits();
    std::vector<AMRCoverageStats> getCoverageStats();
    
    // Methods for EM algorithm with accumulated hits
    void setAccumulatedHits(const std::vector<AMRHit>& hits);
    std::vector<AMRHit> getAllAccumulatedHits() const;
    void addBatchHits(const std::vector<AMRHit>& batch_hits);
    
    // Get database information
    uint32_t getNumGenes() const { 
        return amr_db ? amr_db->getNumGenes() : 0; 
    }
    
    AMRGeneEntry* getGPUGeneEntries() const { 
        return amr_db ? amr_db->getGPUGeneEntries() : nullptr; 
    }
    
    // Get configuration
    const AMRDetectionConfig& getConfig() const { return config; }
    
    // Get gene entries as vector (for report generation)
    std::vector<AMRGeneEntry> getGeneEntries() const {
        std::vector<AMRGeneEntry> entries;
        if (amr_db) {
            uint32_t num_genes = amr_db->getNumGenes();
            entries.resize(num_genes);
            
            AMRGeneEntry* gpu_entries = amr_db->getGPUGeneEntries();
            if (gpu_entries && num_genes > 0) {
                cudaMemcpy(entries.data(), gpu_entries, 
                          num_genes * sizeof(AMRGeneEntry), 
                          cudaMemcpyDeviceToHost);
            }
        }
        return entries;
    }
    
    // Get specific gene entry
    AMRGeneEntry getGeneEntry(uint32_t gene_id) const {
        AMRGeneEntry entry = {};
        if (amr_db && gene_id < amr_db->getNumGenes()) {
            AMRGeneEntry* gpu_entries = amr_db->getGPUGeneEntries();
            if (gpu_entries) {
                cudaMemcpy(&entry, &gpu_entries[gene_id], 
                          sizeof(AMRGeneEntry), 
                          cudaMemcpyDeviceToHost);
            }
        }
        return entry;
    }
    
    // Output methods
    void writeResults(const std::string& output_prefix);
    void generateClinicalReport(const std::string& output_file);
    void exportAbundanceTable(const std::string& output_file);
    
    // Clear accumulated results (useful between samples)
    void clearResults() {
        // Finalize stats before clearing if we have processed reads
        if (total_reads_processed > 0 && !h_coverage_stats.empty()) {
            finalizeCoverageStats();
        }
        
        // Clear host-side coverage statistics
        h_coverage_stats.clear();
        
        // Reset GPU hit counts to zero
        if (d_hit_counts) {
            cudaMemset(d_hit_counts, 0, config.reads_per_batch * sizeof(uint32_t));
        }
        
        // Reset coverage statistics on GPU
        if (d_coverage_stats && amr_db) {
            uint32_t num_genes = amr_db->getNumGenes();
            cudaMemset(d_coverage_stats, 0, num_genes * sizeof(AMRCoverageStats));
        }
        
        // Clear paired read info
        paired_read_info.clear();
        
        // Reset batch info
        current_batch_size = 0;
        total_reads_processed = 0;
        
        // Re-initialize coverage stats for next sample
        if (amr_db) {
            initializeCoverageStats();
        }
    }
    
private:
    // Memory management
    void allocateGPUMemory();
    void freeGPUMemory();
    void copyReadsToGPU(const std::vector<std::string>& reads);
    
    // Translated search engine management
    void resetTranslatedSearchEngine();
    
    // Paired-end concordance scoring (Note: ProteinMatch is defined in cpp file)
    void applyPairedConcordanceScoring(void* matches, uint32_t* match_counts, int num_reads);
    
    // Enhanced paired-end scoring integration
    void integratePairedEndScoring();
    
    // Fragment length distribution methods
    void estimateFragmentLengthDistribution(const std::vector<ReadPairAlignment>& alignments);
    float calculateFragmentLengthProbability(int32_t fragment_length) const;
    void updatePairAlignmentScores(std::vector<ReadPairAlignment>& alignments);
    void buildEnhancedPairAlignments(const std::vector<AMRHit>& batch_hits, size_t num_pairs);
    
    // Paired-end hit update methods
    void updateHitsWithPairedInfo(std::vector<AMRHit>& hits, const std::vector<ReadPairData>& read_pairs);
    void updateGPUHitsWithPairedInfo(const std::vector<AMRHit>& hits);
    
    
    // Gene family extraction methods
    std::string extractGeneFamily(const std::string& gene_name);
    
    // Hit filtering and assignment building
    float calculateAssignmentScore(const AMRHit& hit);
    bool isHighQualityHit(const AMRHit& hit);
    void buildReadAssignments(const std::vector<AMRHit>& all_hits);
    void buildPairedReadAssignments(const std::vector<AMRHit>& all_hits);
    
    // Paired-end EM methods
    void runPairedEndEM();
    void updatePairedAssignmentProbabilities();
    void updateGeneAbundancesFromPairs();
    float calculatePairProbability(const PairedReadAssignment& pair, uint32_t gene_id);
    
    // Sequence similarity calculation
    float calculateNameSimilarity(const std::string& name1, const std::string& name2);
    std::pair<std::string, std::string> parseGeneName(const std::string& name);
    float calculateCTXMSimilarity(const std::string& var1, const std::string& var2);
    float calculateTEMSimilarity(const std::string& var1, const std::string& var2);
    
    // Gene abundance initialization
    float getGeneFamilyPrior(const std::string& family, const std::string& gene_name);
    
    // Integration and reporting
    void updateCoverageStatsFromKallistoEM();
    void reportEMResults();
    void analyzeBetaLactamaseAssignments();

    // Add these as private member variables:
    std::map<std::string, std::vector<uint32_t>> gene_families_map;
    std::vector<std::vector<float>> gene_similarity_matrix;
    std::vector<ReadAssignment> read_assignments;
    std::vector<PairedReadAssignment> paired_read_assignments;  // For paired-end EM
    std::map<uint32_t, GeneAbundanceInfo> gene_abundance_info;
    std::vector<AMRHit> accumulated_hits;  // Store all hits across batches

    // Configuration constants for EM
    static constexpr float MIN_SMITH_WATERMAN_SCORE = 30.0f;  // Lowered for identity*length scoring
    static constexpr float EM_CONVERGENCE_THRESHOLD = 0.001f;
    static constexpr int MAX_EM_ITERATIONS = 100;
    static constexpr float SIMILARITY_WEIGHT = 0.3f;  // Weight for sequence similarity in EM
    
    // Fragment length distribution (learned from data)
    float estimated_fragment_mean = 500.0f;
    float estimated_fragment_std = 100.0f;
    bool fragment_dist_estimated = false;
};

#endif // AMR_DETECTION_PIPELINE_H
