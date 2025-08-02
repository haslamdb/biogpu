#include "resistance_detector.h"
#include <iostream>
#include <fstream>
#include <filesystem>

// Include existing resistance pipeline headers
#include "fq_mutation_detector.cuh"
#include "global_fq_resistance_mapper.h"
#include "hdf5_alignment_writer.h"
#include "clean_resistance_pipeline.h"
#include "diagnostic_report.h"

namespace fs = std::filesystem;

// Internal data structure to hold pipeline state
struct ResistancePipelineData {
    void* cuda_pipeline;  // CUDA pipeline handle
    void* nucleotide_index;
    void* protein_database;
    void* fq_resistance_data;
    
    // Results accumulation
    Json::Value mutation_results;
    std::vector<std::pair<std::string, double>> allele_frequencies;
    
    ResistancePipelineData() : 
        cuda_pipeline(nullptr),
        nucleotide_index(nullptr),
        protein_database(nullptr),
        fq_resistance_data(nullptr) {}
};

ResistanceDetector::ResistanceDetector() : 
    pipeline_data(std::make_unique<ResistancePipelineData>()) {
}

ResistanceDetector::~ResistanceDetector() {
    // Clean up CUDA resources
    if (pipeline_data->cuda_pipeline) {
        // TODO: Call cleanup functions from the resistance pipeline
        // destroy_resistance_pipeline(pipeline_data->cuda_pipeline);
    }
}

void ResistanceDetector::setResistanceConfig(const ResistanceConfig& res_config) {
    resistance_config = res_config;
    config = res_config;  // Copy base config
}

bool ResistanceDetector::initialize(const DetectorConfig& config) {
    this->config = config;
    
    reportProgress("resistance_init", 0, "Initializing resistance detection pipeline");
    
    // Initialize CUDA pipeline
    // TODO: Call actual initialization functions
    // pipeline_data->cuda_pipeline = initialize_resistance_pipeline(
    //     config.gpu_device,
    //     config.batch_size,
    //     config.use_bloom_filter,
    //     resistance_config.enable_sw_alignment
    // );
    
    if (!loadDatabases()) {
        return false;
    }
    
    reportProgress("resistance_init", 100, "Resistance pipeline initialized");
    return true;
}

bool ResistanceDetector::loadDatabases() {
    reportProgress("resistance_db_load", 0, "Loading resistance databases");
    
    // Load nucleotide k-mer index
    if (!resistance_config.nucleotide_index_path.empty()) {
        reportProgress("resistance_db_load", 33, "Loading nucleotide k-mer index");
        // TODO: Call actual loading function
        // pipeline_data->nucleotide_index = load_nucleotide_index(
        //     resistance_config.nucleotide_index_path
        // );
    }
    
    // Load protein database
    if (!resistance_config.protein_db_path.empty()) {
        reportProgress("resistance_db_load", 66, "Loading protein database");
        // TODO: Call actual loading function
        // pipeline_data->protein_database = load_protein_database(
        //     resistance_config.protein_db_path
        // );
    }
    
    // Load FQ resistance mutations if provided
    if (!resistance_config.fq_csv_path.empty()) {
        reportProgress("resistance_db_load", 90, "Loading FQ resistance mutations");
        // TODO: Call actual loading function
        // pipeline_data->fq_resistance_data = load_fq_mutations(
        //     resistance_config.fq_csv_path
        // );
    }
    
    reportProgress("resistance_db_load", 100, "Databases loaded successfully");
    return true;
}

void ResistanceDetector::processBatch(const ReadBatch& batch) {
    reportProgress("resistance_batch", 0, 
        "Processing batch of " + std::to_string(batch.num_reads) + " reads");
    
    // TODO: Convert ReadBatch to format expected by resistance pipeline
    // and call the actual processing function
    
    // Example pseudo-code:
    // ResistanceBatch res_batch;
    // convert_batch(batch, res_batch);
    // process_resistance_batch(pipeline_data->cuda_pipeline, res_batch);
    
    reportProgress("resistance_batch", 100, "Batch processing complete");
}

void ResistanceDetector::finalize() {
    reportProgress("resistance_finalize", 0, "Finalizing resistance detection");
    
    // Process accumulated results
    processResistanceMutations();
    
    // Generate allele frequency analysis
    generateAlleleFrequencyReport();
    
    reportProgress("resistance_finalize", 100, "Resistance detection finalized");
}

void ResistanceDetector::processResistanceMutations() {
    // TODO: Implement mutation processing logic
    // This would aggregate results from all batches
}

void ResistanceDetector::generateAlleleFrequencyReport() {
    // TODO: Implement allele frequency calculation
    // This would analyze depth and frequency of mutations
}

void ResistanceDetector::getResults(Json::Value& results) {
    results["resistance_mutations"] = pipeline_data->mutation_results;
    
    Json::Value allele_freq_array(Json::arrayValue);
    for (const auto& [mutation, frequency] : pipeline_data->allele_frequencies) {
        Json::Value entry;
        entry["mutation"] = mutation;
        entry["frequency"] = frequency;
        allele_freq_array.append(entry);
    }
    results["allele_frequencies"] = allele_freq_array;
    
    results["pipeline"] = "resistance";
    results["version"] = "1.0.0";
}

void ResistanceDetector::writeOutputFiles(const std::string& output_dir) {
    reportProgress("resistance_output", 0, "Writing resistance detection results");
    
    fs::path output_path(output_dir);
    fs::create_directories(output_path);
    
    // Write JSON results
    std::string json_path = (output_path / (config.sample_id + "_resistance.json")).string();
    std::ofstream json_file(json_path);
    Json::StreamWriterBuilder writer;
    std::unique_ptr<Json::StreamWriter> json_writer(writer.newStreamWriter());
    json_writer->write(pipeline_data->mutation_results, &json_file);
    json_file.close();
    
    // Write HDF5 results
    std::string hdf5_path = (output_path / (config.sample_id + "_resistance.h5")).string();
    // TODO: Call HDF5 writer
    
    // Write clinical report
    std::string html_path = (output_path / (config.sample_id + "_resistance_report.html")).string();
    // TODO: Generate HTML report
    
    reportProgress("resistance_output", 100, "Results written successfully");
}

// Factory function
std::unique_ptr<DetectorInterface> createResistanceDetector() {
    return std::make_unique<ResistanceDetector>();
}