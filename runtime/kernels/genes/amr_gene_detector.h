#ifndef AMR_GENE_DETECTOR_H
#define AMR_GENE_DETECTOR_H

#include "../shared/detector_interface.h"
#include <memory>

// Forward declarations
struct AMRPipelineData;

// Configuration specific to AMR gene detection
struct AMRGeneConfig : public DetectorConfig {
    std::string amr_db_path;
    std::string protein_db_path;
    bool enable_em_algorithm;
    int em_iterations;
    double min_hit_coverage;
    int kmer_length;
    bool merge_paired_reads;
    
    AMRGeneConfig() :
        DetectorConfig(),
        enable_em_algorithm(false),
        em_iterations(10),
        min_hit_coverage(0.80),
        kmer_length(31),
        merge_paired_reads(true) {}
};

class AMRGeneDetector : public DetectorInterface {
public:
    AMRGeneDetector();
    ~AMRGeneDetector() override;
    
    // Set extended configuration
    void setAMRConfig(const AMRGeneConfig& amr_config);
    
    // DetectorInterface implementation
    bool initialize(const DetectorConfig& config) override;
    void processBatch(const ReadBatch& batch) override;
    void finalize() override;
    void getResults(Json::Value& results) override;
    void writeOutputFiles(const std::string& output_dir) override;
    std::string getName() const override { return "AMRGeneDetector"; }
    
private:
    std::unique_ptr<AMRPipelineData> pipeline_data;
    AMRGeneConfig amr_config;
    
    // Internal methods
    bool loadAMRDatabases();
    void runEMAlgorithm();
    void generateAbundanceTables();
    void generateClinicalReport();
};

// Factory function
std::unique_ptr<DetectorInterface> createAMRGeneDetector();

#endif // AMR_GENE_DETECTOR_H