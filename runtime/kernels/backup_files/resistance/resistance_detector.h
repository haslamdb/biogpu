#ifndef RESISTANCE_DETECTOR_H
#define RESISTANCE_DETECTOR_H

#include "../shared/detector_interface.h"
#include <memory>

// Forward declarations for resistance pipeline internals
struct ResistancePipelineData;

// Configuration specific to resistance detection
struct ResistanceConfig : public DetectorConfig {
    std::string nucleotide_index_path;
    std::string protein_db_path;
    std::string fq_csv_path;  // Optional FQ mutations CSV
    bool enable_sw_alignment;
    int min_allele_depth;
    int min_report_depth;
    
    ResistanceConfig() : 
        DetectorConfig(),
        enable_sw_alignment(true),
        min_allele_depth(5),
        min_report_depth(0) {}
};

class ResistanceDetector : public DetectorInterface {
public:
    ResistanceDetector();
    ~ResistanceDetector() override;
    
    // Set extended configuration
    void setResistanceConfig(const ResistanceConfig& res_config);
    
    // DetectorInterface implementation
    bool initialize(const DetectorConfig& config) override;
    void processBatch(const ReadBatch& batch) override;
    void finalize() override;
    void getResults(Json::Value& results) override;
    void writeOutputFiles(const std::string& output_dir) override;
    std::string getName() const override { return "ResistanceDetector"; }
    
private:
    std::unique_ptr<ResistancePipelineData> pipeline_data;
    ResistanceConfig resistance_config;
    
    // Internal methods
    bool loadDatabases();
    void processResistanceMutations();
    void generateAlleleFrequencyReport();
};

// Factory function
std::unique_ptr<DetectorInterface> createResistanceDetector();

#endif // RESISTANCE_DETECTOR_H