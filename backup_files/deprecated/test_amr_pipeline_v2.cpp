// runtime/kernels/genes/test_amr_pipeline_v2.cpp
#include "amr_detection_pipeline_v2.h"
#include "../../common/io/streaming_fastq_reader.h"
#include "../../common/gpu/gpu_sequence_buffer.h"
#include "../../common/config/unified_config.h"
#include <iostream>
#include <chrono>

using namespace BioGPU;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fastq_file> [config_file]\n";
        return 1;
    }
    
    std::string fastq_file = argv[1];
    std::string config_file = argc > 2 ? argv[2] : "";
    
    // Load configuration
    UnifiedConfig config;
    if (!config_file.empty()) {
        auto loaded_config = ConfigLoader::loadFromFile(config_file);
        if (loaded_config) {
            config = *loaded_config;
            std::cout << "Loaded configuration from " << config_file << std::endl;
        } else {
            std::cerr << "Warning: Failed to load config file, using defaults\n";
            config = ConfigLoader::getDefault();
        }
    } else {
        config = ConfigLoader::getDefault();
        std::cout << "Using default configuration\n";
    }
    
    // Validate configuration
    std::vector<std::string> errors;
    if (!ConfigLoader::validate(config, errors)) {
        std::cerr << "Configuration errors:\n";
        for (const auto& error : errors) {
            std::cerr << "  - " << error << "\n";
        }
        return 1;
    }
    
    // Create AMR detection pipeline
    AMRDetectionPipeline amr_pipeline(config.genes_config);
    
    // Initialize pipeline
    std::cout << "Initializing AMR detection pipeline...\n";
    if (!amr_pipeline.initialize()) {
        std::cerr << "Failed to initialize AMR pipeline\n";
        return 1;
    }
    
    // Create streaming FASTQ reader
    StreamingFastqReader reader(config.pipeline_config.batch_size);
    if (!reader.open(fastq_file)) {
        std::cerr << "Failed to open FASTQ file: " << fastq_file << "\n";
        return 1;
    }
    
    // Create GPU sequence buffer
    GPUSequenceBuffer gpu_buffer(
        config.pipeline_config.max_gpu_sequences,
        config.pipeline_config.max_gpu_bases
    );
    if (!gpu_buffer.allocate()) {
        std::cerr << "Failed to allocate GPU buffer\n";
        return 1;
    }
    
    // Process batches
    std::cout << "Processing FASTQ file...\n";
    auto start_time = std::chrono::high_resolution_clock::now();
    
    size_t batch_count = 0;
    while (reader.hasNext()) {
        auto batch = reader.getNextBatch();
        if (!batch || batch->empty()) {
            break;
        }
        
        std::cout << "Processing batch " << ++batch_count 
                  << " with " << batch->size() << " sequences...\n";
        
        // Transfer batch to GPU
        if (!gpu_buffer.transferBatch(*batch)) {
            std::cerr << "Failed to transfer batch to GPU\n";
            continue;
        }
        
        // Process batch with AMR pipeline
        if (!amr_pipeline.processBatch(&gpu_buffer)) {
            std::cerr << "Failed to process batch\n";
            continue;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    // Get results
    auto results = amr_pipeline.getResults();
    auto* amr_results = dynamic_cast<AMRResults*>(results.get());
    
    if (amr_results) {
        std::cout << "\n=== Processing Complete ===\n";
        std::cout << "Total time: " << duration.count() << " seconds\n";
        std::cout << "Total reads: " << reader.getTotalReads() << "\n";
        std::cout << "Total bases: " << reader.getTotalBases() << "\n";
        std::cout << "Batches processed: " << batch_count << "\n";
        
        // Print summary
        std::cout << "\n";
        amr_results->writeSummary(std::cout);
        
        // Write detailed reports
        std::string output_prefix = config.pipeline_config.output_dir + "/amr_results";
        amr_results->writeReport(output_prefix + "_detailed.tsv");
        amr_results->writeClinicalReport(output_prefix + "_clinical_report.txt");
        
        std::cout << "\nResults written to:\n";
        std::cout << "  - " << output_prefix << "_detailed.tsv\n";
        std::cout << "  - " << output_prefix << "_clinical_report.txt\n";
        
        // Performance metrics
        if (config.genes_config.enable_profiling) {
            std::cout << "\n=== Performance Metrics ===\n";
            auto batch_times = amr_pipeline.getBatchTimes();
            if (!batch_times.empty()) {
                float total_gpu_time = 0;
                for (float t : batch_times) total_gpu_time += t;
                
                std::cout << "Average batch GPU time: " 
                          << amr_pipeline.getAverageBatchTime() << " ms\n";
                std::cout << "Total GPU time: " << total_gpu_time / 1000.0 << " seconds\n";
                std::cout << "GPU utilization: " 
                          << (total_gpu_time / 1000.0) / duration.count() * 100 << "%\n";
            }
        }
    }
    
    // Cleanup
    reader.close();
    gpu_buffer.free();
    amr_pipeline.cleanup();
    
    return 0;
}