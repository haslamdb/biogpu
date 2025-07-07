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

// Structure to track paired-end relationships
struct PairedReadInfo {
    uint32_t read_idx;      // Original read pair index
    bool is_read2;          // false for R1, true for R2
    uint32_t pair_offset;   // Offset to find the paired read results
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
    
    // Alignment parameters
    float min_identity = 0.90f;             // Minimum identity for AMR match
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
    bool is_complete_gene;
    bool concordant;        // For paired-end reads, true if both reads hit same gene
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

class AMRDetectionPipeline {
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
    
    // Gene entries from database
    std::vector<AMRGeneEntry> gene_entries;
    
    // Translated search engine (reused across batches)
    void* translated_search_engine;
    bool search_engine_initialized;
    
    // Paired-end read tracking
    std::vector<PairedReadInfo> paired_read_info;
    
public:
    AMRDetectionPipeline(const AMRDetectionConfig& cfg);
    ~AMRDetectionPipeline();
    
    // Initialize pipeline
    bool initialize(const std::string& amr_db_path);
    
    // Process a batch of reads
    void processBatch(const std::vector<std::string>& reads,
                     const std::vector<std::string>& read_ids);
    
    // Process a batch of paired-end reads
    void processBatchPaired(const std::vector<std::string>& reads1,
                           const std::vector<std::string>& reads2,
                           const std::vector<std::string>& read_ids);
    
    // Main processing steps
    void buildBloomFilter();
    void generateMinimizers();
    void screenWithBloomFilter();
    void performTranslatedAlignment();
    void extendAlignments();
    void calculateCoverageStats();
    void calculateAbundanceMetrics();
    
    // Coverage stats management
    void initializeCoverageStats();
    void freeCoverageStats();
    
    // Get results
    std::vector<AMRHit> getAMRHits();
    std::vector<AMRCoverageStats> getCoverageStats();
    
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
};

#endif // AMR_DETECTION_PIPELINE_H
