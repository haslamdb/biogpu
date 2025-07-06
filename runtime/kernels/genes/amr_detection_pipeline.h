// amr_detection_pipeline.h
#ifndef AMR_DETECTION_PIPELINE_H
#define AMR_DETECTION_PIPELINE_H

#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <memory>
#include "ncbi_amr_database_loader.h"
#include "sample_csv_parser.h"  // Reuse from FQ pipeline
// #include "bloom_filter.h"       // Reuse bloom filter

// Configuration for AMR detection
struct AMRDetectionConfig {
    // Bloom filter parameters
    size_t bloom_filter_size = 1ULL << 30;  // 1GB bloom filter
    int bloom_k = 3;                        // Number of hash functions
    int kmer_length = 31;                   // DNA k-mer size for bloom filter
    
    // Minimizer parameters
    int minimizer_k = 15;                   // DNA minimizer k-mer size
    int minimizer_w = 10;                   // Window size
    
    // Protein search parameters
    int protein_kmer_size = 8;              // For protein alignment after translation
    
    // Alignment parameters
    float min_identity = 0.85f;             // Minimum identity for AMR match
    float min_coverage = 0.80f;             // Minimum coverage of AMR gene
    int min_alignment_length = 50;          // Minimum alignment length (aa)
    int band_width = 15;                    // Band width for banded SW
    
    // Batch processing
    int reads_per_batch = 100000;
    int max_read_length = 300;
    
    // Output options
    bool output_sam = false;
    bool output_hdf5 = true;
    std::string output_prefix = "amr_results";
    
    // Database paths
    std::string protein_db_path = "amr_protein_db";  // Default location
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
    uint8_t num_mutations;
    bool is_complete_gene;
    char gene_name[64];
    char drug_class[32];
};

// Coverage tracking per AMR gene
struct AMRCoverageStats {
    uint32_t total_reads;           // Number of reads mapped to this gene
    uint32_t total_bases_mapped;    // Total bases from all reads mapped
    uint32_t covered_positions;     // Number of positions with at least 1x coverage
    uint32_t gene_length;           // Length of the gene (in amino acids for proteins)
    float percent_coverage;         // Percentage of gene covered (0-100)
    float mean_depth;              // Average depth across covered positions
    float rpkm;                    // Reads Per Kilobase per Million mapped reads
    float tpm;                     // Transcripts Per Million
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

class AMRDetectionPipeline {
private:
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
    
    // GPU memory for alignment results
    AMRHit* d_amr_hits;
    uint32_t* d_hit_counts;
    
    // Coverage statistics
    AMRCoverageStats* d_coverage_stats;
    std::vector<AMRCoverageStats> h_coverage_stats;
    
    // Batch processing
    int current_batch_size;
    uint64_t total_reads_processed;
    
    // Gene entries from database
    std::vector<AMRGeneEntry> gene_entries;
    
public:
    AMRDetectionPipeline(const AMRDetectionConfig& cfg);
    ~AMRDetectionPipeline();
    
    // Initialize pipeline
    bool initialize(const std::string& amr_db_path);
    
    // Process a batch of reads
    void processBatch(const std::vector<std::string>& reads,
                     const std::vector<std::string>& read_ids);
    
    // Main processing steps
    void buildBloomFilter();
    void generateMinimizers();
    void screenWithBloomFilter();
    void performTranslatedAlignment();
    void extendAlignments();
    void calculateCoverageStats();
    void calculateAbundanceMetrics();
    
    // Get results
    std::vector<AMRHit> getAMRHits();
    std::vector<AMRCoverageStats> getCoverageStats();
    
    // Output methods
    void writeResults(const std::string& output_prefix);
    void generateClinicalReport(const std::string& output_file);
    void exportAbundanceTable(const std::string& output_file);
    
private:
    // Memory management
    void allocateGPUMemory();
    void freeGPUMemory();
    void copyReadsToGPU(const std::vector<std::string>& reads);
};

#endif // AMR_DETECTION_PIPELINE_H
