// runtime/kernels/genes/hdf5_amr_writer.h
#ifndef HDF5_AMR_WRITER_H
#define HDF5_AMR_WRITER_H

#include <H5Cpp.h>
#include <vector>
#include <string>
#include <map>
#include "amr_detection_pipeline.h"

class HDF5AMRWriter {
private:
    H5::H5File* file;
    std::string filename;
    
    // Buffers for AMR hit data
    std::vector<uint32_t> read_ids;
    std::vector<uint32_t> gene_ids;
    std::vector<std::string> gene_names;
    std::vector<std::string> gene_families;
    std::vector<std::string> drug_classes;
    std::vector<float> identities;
    std::vector<float> coverages;
    std::vector<uint16_t> ref_starts;
    std::vector<uint16_t> ref_ends;
    std::vector<int8_t> frames;
    std::vector<bool> is_complete_genes;
    std::vector<bool> concordant_flags;
    
    // Coverage statistics buffers
    std::vector<uint32_t> cov_gene_ids;
    std::vector<std::string> cov_gene_names;
    std::vector<uint32_t> total_reads;
    std::vector<uint32_t> total_bases_mapped;
    std::vector<uint16_t> covered_positions;
    std::vector<uint16_t> gene_lengths;
    std::vector<float> percent_coverages;
    std::vector<float> mean_depths;
    std::vector<float> rpkms;
    std::vector<float> tpms;
    
    // Positional coverage data (for visualization)
    std::map<uint32_t, std::vector<uint32_t>> position_coverage_map;
    
    uint64_t total_hits_written;
    uint64_t total_genes_detected;
    
    // Control automatic flushing
    bool auto_flush_enabled;
    
public:
    HDF5AMRWriter(const std::string& output_path);
    ~HDF5AMRWriter();
    
    void initialize(const std::string& sample_name, 
                   const std::string& dna_fasta_path,
                   const std::string& protein_fasta_path);
    
    void addAMRHits(const std::vector<AMRHit>& hits);
    void addCoverageStats(const std::vector<AMRCoverageStats>& coverage_stats,
                         const std::vector<AMRGeneEntry>& gene_entries);
    void addPositionalCoverage(uint32_t gene_id, const uint32_t* position_counts, uint16_t gene_length);
    
    void flush();
    void finalize(const std::string& json_summary_path);
    
    // Control when flushing happens
    void setAutoFlush(bool enable) { auto_flush_enabled = enable; }
    void manualFlush() { flush(); }
    
private:
    void createDatasets();
    void createCoverageDatasets();
    void writeMetadata();
};

#endif // HDF5_AMR_WRITER_H
