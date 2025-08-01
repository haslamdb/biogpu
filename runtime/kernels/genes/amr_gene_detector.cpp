#include "amr_gene_detector.h"
#include <iostream>
#include <fstream>
#include <filesystem>

// Include existing AMR pipeline headers
extern "C" {
    #include "amr_detection_kernels.h"
}
#include "amr_detection_pipeline.h"
#include "hdf5_amr_writer.h"
#include "clinical_amr_report_generator.h"

namespace fs = std::filesystem;

// Internal data structure to hold AMR pipeline state
struct AMRPipelineData {
    void* amr_pipeline;  // AMR pipeline handle
    void* amr_database;
    void* protein_database;
    
    // Results accumulation
    Json::Value amr_hits;
    Json::Value gene_abundances;
    Json::Value drug_resistance_profile;
    std::map<std::string, double> gene_coverage;
    
    AMRPipelineData() : 
        amr_pipeline(nullptr),
        amr_database(nullptr),
        protein_database(nullptr) {}
};

AMRGeneDetector::AMRGeneDetector() :
    pipeline_data(std::make_unique<AMRPipelineData>()) {
}

AMRGeneDetector::~AMRGeneDetector() {
    // Clean up resources
    if (pipeline_data->amr_pipeline) {
        // TODO: Call cleanup functions
        // destroy_amr_pipeline(pipeline_data->amr_pipeline);
    }
}

void AMRGeneDetector::setAMRConfig(const AMRGeneConfig& amr_config) {
    this->amr_config = amr_config;
    config = amr_config;  // Copy base config
}

bool AMRGeneDetector::initialize(const DetectorConfig& config) {
    this->config = config;
    
    reportProgress("amr_init", 0, "Initializing AMR gene detection pipeline");
    
    // Initialize AMR detection pipeline
    // TODO: Call actual initialization
    // AMRConfig pipeline_config;
    // pipeline_config.gpu_device = config.gpu_device;
    // pipeline_config.batch_size = config.batch_size;
    // pipeline_config.min_identity = config.min_identity;
    // pipeline_config.min_coverage = config.min_coverage;
    // pipeline_config.use_bloom_filter = config.use_bloom_filter;
    // 
    // pipeline_data->amr_pipeline = initialize_amr_pipeline(pipeline_config);
    
    if (!loadAMRDatabases()) {
        return false;
    }
    
    reportProgress("amr_init", 100, "AMR pipeline initialized");
    return true;
}

bool AMRGeneDetector::loadAMRDatabases() {
    reportProgress("amr_db_load", 0, "Loading AMR databases");
    
    // Load AMR reference database
    if (!amr_config.amr_db_path.empty()) {
        reportProgress("amr_db_load", 50, "Loading AMR reference database");
        // TODO: Call actual loading function
        // pipeline_data->amr_database = load_amr_database(amr_config.amr_db_path);
        
        if (!pipeline_data->amr_database) {
            std::cerr << "Failed to load AMR database from: " << amr_config.amr_db_path << std::endl;
            return false;
        }
    }
    
    // Load protein database if specified
    if (!amr_config.protein_db_path.empty()) {
        reportProgress("amr_db_load", 80, "Loading protein database");
        // TODO: Call actual loading function
        // pipeline_data->protein_database = load_protein_database(amr_config.protein_db_path);
    }
    
    reportProgress("amr_db_load", 100, "AMR databases loaded");
    return true;
}

void AMRGeneDetector::processBatch(const ReadBatch& batch) {
    reportProgress("amr_batch", 0, 
        "Processing batch of " + std::to_string(batch.num_reads) + " reads");
    
    // TODO: Convert ReadBatch to AMR pipeline format and process
    // AMRBatch amr_batch;
    // convert_to_amr_batch(batch, amr_batch, amr_config.merge_paired_reads);
    // 
    // AMRResults batch_results;
    // process_amr_batch(pipeline_data->amr_pipeline, amr_batch, batch_results);
    // 
    // // Accumulate results
    // accumulate_amr_results(pipeline_data->amr_hits, batch_results);
    
    reportProgress("amr_batch", 100, "Batch processing complete");
}

void AMRGeneDetector::finalize() {
    reportProgress("amr_finalize", 0, "Finalizing AMR gene detection");
    
    // Run EM algorithm if enabled
    if (amr_config.enable_em_algorithm) {
        reportProgress("amr_finalize", 30, "Running EM algorithm for multi-mapping reads");
        runEMAlgorithm();
    }
    
    // Generate abundance tables
    reportProgress("amr_finalize", 60, "Generating gene abundance tables");
    generateAbundanceTables();
    
    // Generate clinical report data
    reportProgress("amr_finalize", 90, "Generating clinical report");
    generateClinicalReport();
    
    reportProgress("amr_finalize", 100, "AMR detection finalized");
}

void AMRGeneDetector::runEMAlgorithm() {
    // TODO: Implement EM algorithm for resolving multi-mapping reads
    // This would refine gene abundance estimates
}

void AMRGeneDetector::generateAbundanceTables() {
    // TODO: Calculate normalized abundance values
    // Generate coverage statistics
    // Create gene family summaries
}

void AMRGeneDetector::generateClinicalReport() {
    // TODO: Analyze AMR genes for drug resistance implications
    // Map genes to drug classes
    // Generate resistance profile
}

void AMRGeneDetector::getResults(Json::Value& results) {
    results["amr_genes"] = pipeline_data->amr_hits;
    results["gene_abundances"] = pipeline_data->gene_abundances;
    results["drug_resistance"] = pipeline_data->drug_resistance_profile;
    
    Json::Value coverage_stats(Json::objectValue);
    for (const auto& [gene, coverage] : pipeline_data->gene_coverage) {
        coverage_stats[gene] = coverage;
    }
    results["coverage_statistics"] = coverage_stats;
    
    results["pipeline"] = "amr_genes";
    results["version"] = "1.0.0";
}

void AMRGeneDetector::writeOutputFiles(const std::string& output_dir) {
    reportProgress("amr_output", 0, "Writing AMR gene detection results");
    
    fs::path output_path(output_dir);
    fs::create_directories(output_path);
    
    // Write JSON results
    std::string json_path = (output_path / (config.sample_id + "_amr_genes.json")).string();
    std::ofstream json_file(json_path);
    Json::StreamWriterBuilder writer;
    std::unique_ptr<Json::StreamWriter> json_writer(writer.newStreamWriter());
    json_writer->write(pipeline_data->amr_hits, &json_file);
    json_file.close();
    
    // Write abundance TSV
    std::string tsv_path = (output_path / (config.sample_id + "_amr_abundance.tsv")).string();
    // TODO: Write TSV format abundance table
    
    // Write HDF5 results
    std::string hdf5_path = (output_path / (config.sample_id + "_amr_genes.h5")).string();
    // TODO: Call HDF5 writer
    
    // Write clinical HTML report
    std::string html_path = (output_path / (config.sample_id + "_amr_clinical_report.html")).string();
    // TODO: Generate HTML report
    
    reportProgress("amr_output", 100, "AMR results written successfully");
}

// Factory function
std::unique_ptr<DetectorInterface> createAMRGeneDetector() {
    return std::make_unique<AMRGeneDetector>();
}