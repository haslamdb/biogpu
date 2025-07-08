// runtime/kernels/genes/clinical_amr_report_generator.h
#ifndef CLINICAL_AMR_REPORT_GENERATOR_H
#define CLINICAL_AMR_REPORT_GENERATOR_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "amr_detection_pipeline.h"

class ClinicalAMRReportGenerator {
public:
    struct DrugClassSummary {
        std::string drug_class;
        std::vector<std::string> genes_detected;
        std::vector<std::string> high_confidence_genes;  // >90% coverage
        std::vector<std::string> moderate_confidence_genes;  // 50-90% coverage
        std::vector<std::string> low_confidence_genes;  // <50% coverage
        int total_reads;
        float max_tpm;
        std::string clinical_interpretation;
    };
    
    struct GeneSummary {
        std::string gene_name;
        std::string gene_family;
        std::string drug_class;
        uint32_t read_count;
        float percent_coverage;
        float mean_depth;
        float tpm;
        float rpkm;
        bool is_complete_gene;
        std::string confidence_level;
        std::vector<std::string> sample_hits;  // For multi-sample analysis
    };
    
    struct GeneFamilySummary {
        std::string gene_family;
        std::vector<std::string> gene_variants;  // e.g., blaKPC-2, blaKPC-3
        float total_tpm;
        float total_reads;
        float max_identity;
    };
    
private:
    std::string output_path;
    std::string sample_name;
    
    // Collected data
    std::map<std::string, DrugClassSummary> drug_class_summaries;
    std::map<std::string, GeneSummary> gene_summaries;
    std::map<std::string, GeneFamilySummary> gene_family_summaries;
    
    // Statistics
    int total_reads_processed;
    int reads_with_amr_hits;
    int total_genes_detected;
    int high_confidence_genes;
    
    // Clinical interpretations
    std::map<std::string, std::string> drug_class_interpretations;
    
public:
    ClinicalAMRReportGenerator(const std::string& output_path, const std::string& sample_name);
    
    void processAMRResults(const std::vector<AMRHit>& hits,
                          const std::vector<AMRCoverageStats>& coverage_stats,
                          const std::vector<AMRGeneEntry>& gene_entries,
                          int total_reads);
    
    void generateReports();
    
private:
    void initializeDrugClassInterpretations();
    void generateHTMLReport();
    void generateTextReport();
    void generateJSONReport();
    void generateTSVReports();
    
    std::string getConfidenceLevel(float coverage, float depth);
    std::string generateDrugClassInterpretation(const DrugClassSummary& summary);
    std::string getCurrentTimestamp();
};

#endif // CLINICAL_AMR_REPORT_GENERATOR_H
