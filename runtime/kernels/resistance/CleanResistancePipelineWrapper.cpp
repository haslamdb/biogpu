#include "CleanResistancePipeline.h"

// Include the actual implementation
#define private public  // Hack to access private members - remove in production
#include "clean_resistance_pipeline_main.cpp"
#undef private

// Implementation wrapper
class CleanResistancePipelineImpl {
public:
    // The actual pipeline from clean_resistance_pipeline_main.cpp
    std::unique_ptr<::CleanResistancePipeline> actual_pipeline;
    
    CleanResistancePipelineImpl(bool use_bloom, bool use_sw, 
                               int min_allele_depth, int min_report_depth) {
        actual_pipeline = std::make_unique<::CleanResistancePipeline>(
            use_bloom, use_sw, min_allele_depth, min_report_depth);
    }
};

// Public interface implementation
CleanResistancePipeline::CleanResistancePipeline(bool use_bloom, bool use_sw, 
                                               int min_allele_depth, int min_report_depth) 
    : pImpl(std::make_unique<CleanResistancePipelineImpl>(
        use_bloom, use_sw, min_allele_depth, min_report_depth)) {
}

CleanResistancePipeline::~CleanResistancePipeline() = default;

void CleanResistancePipeline::loadDatabases(const std::string& nucleotide_index,
                                           const std::string& protein_db,
                                           const std::string& fq_csv_path) {
    pImpl->actual_pipeline->loadDatabases(nucleotide_index, protein_db, fq_csv_path);
}

void CleanResistancePipeline::processReads(const std::string& read1_path,
                                          const std::string& read2_path,
                                          const std::string& output_prefix) {
    pImpl->actual_pipeline->processReads(read1_path, read2_path, output_prefix);
}