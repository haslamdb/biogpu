#ifndef CLEAN_RESISTANCE_PIPELINE_H
#define CLEAN_RESISTANCE_PIPELINE_H

#include <string>
#include <memory>

// Forward declaration of the actual implementation
class CleanResistancePipelineImpl;

// Public interface for the resistance pipeline
class CleanResistancePipeline {
public:
    CleanResistancePipeline(bool use_bloom = true, bool use_sw = true, 
                           int min_allele_depth = 5, int min_report_depth = 0);
    ~CleanResistancePipeline();
    
    // Load required databases
    void loadDatabases(const std::string& nucleotide_index,
                      const std::string& protein_db,
                      const std::string& fq_csv_path = "");
    
    // Process paired-end reads
    void processReads(const std::string& read1_path,
                     const std::string& read2_path,
                     const std::string& output_prefix);
    
private:
    std::unique_ptr<CleanResistancePipelineImpl> pImpl;
};

#endif // CLEAN_RESISTANCE_PIPELINE_H