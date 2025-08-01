#ifndef CLEAN_RESISTANCE_PIPELINE_H
#define CLEAN_RESISTANCE_PIPELINE_H

#include <string>
#include <vector>
#include <memory>
#include <json/json.h>

// Forward declaration to avoid including the entire implementation
class CleanResistancePipeline {
public:
    CleanResistancePipeline(bool use_bloom = true, bool use_sw = true, 
                           int min_allele_depth = 5, int min_report_depth = 0);
    ~CleanResistancePipeline();
    
    // Load databases
    void loadDatabases(const std::string& nucleotide_index,
                      const std::string& protein_db,
                      const std::string& fq_csv_path = "");
    
    // Process reads
    void processReads(const std::string& read1_path,
                     const std::string& read2_path,
                     const std::string& output_prefix);
    
    // Process batch of reads (for unified pipeline)
    void processBatch(const std::vector<std::string>& read1_seqs,
                     const std::vector<std::string>& read2_seqs,
                     const std::vector<std::string>& read_ids,
                     int batch_number);
    
    // Get results
    Json::Value getResults() const;
    
    // Write output files
    void writeOutputFiles(const std::string& output_dir);
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

#endif // CLEAN_RESISTANCE_PIPELINE_H