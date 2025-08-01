#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <getopt.h>
#include <json/json.h>
#include <cuda_runtime.h>
#include <thread>
#include <future>
#include <fstream>
#include <gzstream.h>  // For reading gzipped FASTQ files
#include <zlib.h>

// Include the actual pipeline implementations
#include "resistance/CleanResistancePipeline.h"
#include "genes/amr_detection_pipeline.h"
#include "../common/io/streaming_fastq_reader.h"

namespace fs = std::filesystem;

// Progress reporting structure
struct ProgressReporter {
    std::string current_stage;
    int progress_percentage;
    std::string message;
    bool json_output;
    std::mutex output_mutex;
    
    ProgressReporter(bool use_json) : json_output(use_json), progress_percentage(0) {}
    
    void report(const std::string& stage, int percentage, const std::string& msg) {
        std::lock_guard<std::mutex> lock(output_mutex);
        current_stage = stage;
        progress_percentage = percentage;
        message = msg;
        
        if (json_output) {
            Json::Value root;
            root["stage"] = stage;
            root["progress"] = percentage;
            root["message"] = msg;
            Json::StreamWriterBuilder writer;
            std::unique_ptr<Json::StreamWriter> json_writer(writer.newStreamWriter());
            json_writer->write(root, &std::cout);
            std::cout << std::endl;
        } else {
            std::cerr << "[" << std::setw(3) << percentage << "%] " 
                      << stage << ": " << msg << std::endl;
        }
    }
};

// Command line options structure
struct PipelineOptions {
    std::string r1_path;
    std::string r2_path;
    std::string output_dir;
    std::string reference_db;
    std::string resistance_db;
    int gpu_device = 0;
    std::string sample_id;
    bool progress_json = false;
    int threads = 8;
    
    // Pipeline control
    bool run_resistance = true;
    bool run_amr_genes = true;
    
    // Multi-GPU support
    int resistance_gpu = 1;  // Default to A5000 (24GB)
    int amr_gpu = 0;        // Default to A6000 (48GB)
    bool use_multi_gpu = false;
    
    // Advanced options
    bool use_bloom_filter = false;
    double min_identity = 0.90;
    double min_coverage = 0.80;
    int batch_size = 50000;
    int min_allele_depth = 5;
    bool enable_em = false;
    int em_iterations = 10;
};

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " [OPTIONS]" << std::endl;
    std::cerr << "\nRequired arguments:" << std::endl;
    std::cerr << "  --r1 <path>              Forward reads (R1) FASTQ file" << std::endl;
    std::cerr << "  --r2 <path>              Reverse reads (R2) FASTQ file" << std::endl;
    std::cerr << "  --output-dir <path>      Output directory for results" << std::endl;
    std::cerr << "  --reference-db <path>    Microbial reference database path" << std::endl;
    std::cerr << "  --resistance-db <path>   Quinolone resistance database path" << std::endl;
    std::cerr << "  --sample-id <string>     Sample identifier" << std::endl;
    
    std::cerr << "\nOptional arguments:" << std::endl;
    std::cerr << "  --gpu-device <id>        GPU device ID (default: 0)" << std::endl;
    std::cerr << "  --use-multi-gpu          Enable multi-GPU mode" << std::endl;
    std::cerr << "  --resistance-gpu <id>    GPU for resistance pipeline" << std::endl;
    std::cerr << "  --amr-gpu <id>           GPU for AMR genes pipeline" << std::endl;
    std::cerr << "  --progress-json          Output progress in JSON format" << std::endl;
    std::cerr << "  --threads <num>          Number of CPU threads (default: 8)" << std::endl;
    std::cerr << "  --help                   Show this help message" << std::endl;
}

bool parse_arguments(int argc, char* argv[], PipelineOptions& options) {
    static struct option long_options[] = {
        {"r1", required_argument, 0, 0},
        {"r2", required_argument, 0, 0},
        {"output-dir", required_argument, 0, 0},
        {"reference-db", required_argument, 0, 0},
        {"resistance-db", required_argument, 0, 0},
        {"gpu-device", required_argument, 0, 0},
        {"use-multi-gpu", no_argument, 0, 0},
        {"resistance-gpu", required_argument, 0, 0},
        {"amr-gpu", required_argument, 0, 0},
        {"sample-id", required_argument, 0, 0},
        {"progress-json", no_argument, 0, 0},
        {"threads", required_argument, 0, 0},
        {"help", no_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
        if (c == 0) {
            std::string opt_name = long_options[option_index].name;
            
            if (opt_name == "r1") options.r1_path = optarg;
            else if (opt_name == "r2") options.r2_path = optarg;
            else if (opt_name == "output-dir") options.output_dir = optarg;
            else if (opt_name == "reference-db") options.reference_db = optarg;
            else if (opt_name == "resistance-db") options.resistance_db = optarg;
            else if (opt_name == "gpu-device") options.gpu_device = std::stoi(optarg);
            else if (opt_name == "use-multi-gpu") options.use_multi_gpu = true;
            else if (opt_name == "resistance-gpu") options.resistance_gpu = std::stoi(optarg);
            else if (opt_name == "amr-gpu") options.amr_gpu = std::stoi(optarg);
            else if (opt_name == "sample-id") options.sample_id = optarg;
            else if (opt_name == "progress-json") options.progress_json = true;
            else if (opt_name == "threads") options.threads = std::stoi(optarg);
            else if (opt_name == "help") {
                print_usage(argv[0]);
                return false;
            }
        }
    }
    
    // Validate required arguments
    if (options.r1_path.empty() || options.r2_path.empty() || 
        options.output_dir.empty() || options.reference_db.empty() || 
        options.resistance_db.empty() || options.sample_id.empty()) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        print_usage(argv[0]);
        return false;
    }
    
    // Create output directory
    fs::create_directories(options.output_dir);
    
    return true;
}

// Run resistance pipeline on specified GPU
void run_resistance_pipeline(const PipelineOptions& options, ProgressReporter& progress) {
    // Set GPU device for this thread
    cudaSetDevice(options.use_multi_gpu ? options.resistance_gpu : options.gpu_device);
    
    progress.report("resistance_init", 5, "Initializing resistance detection pipeline");
    
    try {
        // Create resistance pipeline instance
        CleanResistancePipeline pipeline(
            options.use_bloom_filter,    // use_bloom
            true,                        // use_sw (Smith-Waterman)
            options.min_allele_depth,    // min_allele_depth
            0                           // min_report_depth
        );
        
        // Load databases
        std::string nucleotide_index = options.reference_db + "/nucleotide_index";
        std::string protein_db = options.resistance_db + "/protein_db";
        std::string fq_csv = options.resistance_db + "/fq_mutations.csv";
        
        progress.report("resistance_init", 10, "Loading resistance databases");
        pipeline.loadDatabases(nucleotide_index, protein_db, fq_csv);
        
        // Process reads
        std::string output_prefix = options.output_dir + "/" + options.sample_id + "_resistance";
        progress.report("resistance_processing", 20, "Processing reads for resistance mutations");
        
        pipeline.processReads(options.r1_path, options.r2_path, output_prefix);
        
        progress.report("resistance_complete", 100, "Resistance detection complete");
        
    } catch (const std::exception& e) {
        progress.report("resistance_error", -1, std::string("Resistance pipeline error: ") + e.what());
        throw;
    }
}

// Run AMR gene detection pipeline on specified GPU
void run_amr_gene_pipeline(const PipelineOptions& options, ProgressReporter& progress) {
    // Set GPU device for this thread
    cudaSetDevice(options.use_multi_gpu ? options.amr_gpu : options.gpu_device);
    
    progress.report("amr_init", 5, "Initializing AMR gene detection pipeline");
    
    try {
        // Configure AMR detection
        AMRDetectionConfig config;
        config.use_bloom_filter = options.use_bloom_filter;
        config.min_identity = options.min_identity;
        config.min_coverage = options.min_coverage;
        config.reads_per_batch = options.batch_size;
        config.output_prefix = options.output_dir + "/" + options.sample_id;
        config.protein_db_path = options.reference_db + "/amr_protein_db";
        config.use_em = options.enable_em;
        config.em_iterations = options.em_iterations;
        
        // Create AMR detection pipeline
        AMRDetectionPipeline pipeline(config);
        
        // Initialize with database
        progress.report("amr_init", 10, "Loading AMR gene database");
        std::string amr_db_path = options.reference_db + "/AMR_CDS.fa," + 
                                  options.reference_db + "/AMRProt.fa";
        
        if (!pipeline.initialize(amr_db_path)) {
            throw std::runtime_error("Failed to initialize AMR detection pipeline");
        }
        
        // Process samples
        std::vector<SampleInfo> samples;
        samples.push_back({
            options.sample_id,
            "",  // base_path not used
            fs::path(options.r1_path).filename().string(),
            fs::path(options.r2_path).filename().string(),
            options.r1_path,
            options.r2_path
        });
        
        progress.report("amr_processing", 20, "Processing reads for AMR genes");
        pipeline.processSamples(samples, options.output_dir);
        
        progress.report("amr_complete", 100, "AMR gene detection complete");
        
    } catch (const std::exception& e) {
        progress.report("amr_error", -1, std::string("AMR pipeline error: ") + e.what());
        throw;
    }
}

int main(int argc, char* argv[]) {
    PipelineOptions options;
    
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    
    // Initialize progress reporter
    ProgressReporter progress(options.progress_json);
    
    progress.report("startup", 0, "Starting BioGPU unified pipeline");
    
    // Validate GPU configuration
    int device_count;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
        std::cerr << "Error: Failed to get GPU device count" << std::endl;
        return 1;
    }
    
    // Log GPU configuration
    if (!options.progress_json) {
        std::cout << "GPU Configuration:" << std::endl;
        for (int i = 0; i < device_count; ++i) {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, i);
            std::cout << "  GPU " << i << ": " << prop.name 
                      << " (" << (prop.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
        }
        
        if (options.use_multi_gpu) {
            std::cout << "\nMulti-GPU mode enabled:" << std::endl;
            std::cout << "  Resistance pipeline: GPU " << options.resistance_gpu << std::endl;
            std::cout << "  AMR genes pipeline: GPU " << options.amr_gpu << std::endl;
        }
    }
    
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (options.use_multi_gpu && options.run_resistance && options.run_amr_genes) {
            // Run both pipelines in parallel on different GPUs
            progress.report("parallel_processing", 10, "Running pipelines in parallel on multiple GPUs");
            
            auto resistance_future = std::async(std::launch::async, 
                run_resistance_pipeline, std::ref(options), std::ref(progress));
            
            auto amr_future = std::async(std::launch::async,
                run_amr_gene_pipeline, std::ref(options), std::ref(progress));
            
            // Wait for both to complete
            resistance_future.get();
            amr_future.get();
            
        } else {
            // Run sequentially on single GPU
            if (options.run_resistance) {
                run_resistance_pipeline(options, progress);
            }
            
            if (options.run_amr_genes) {
                run_amr_gene_pipeline(options, progress);
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Write unified summary
        Json::Value summary;
        summary["sample_id"] = options.sample_id;
        summary["processing_time_seconds"] = duration.count();
        summary["resistance_enabled"] = options.run_resistance;
        summary["amr_genes_enabled"] = options.run_amr_genes;
        summary["multi_gpu_used"] = options.use_multi_gpu;
        summary["output_directory"] = options.output_dir;
        
        std::string summary_path = fs::path(options.output_dir) / (options.sample_id + "_unified_summary.json");
        std::ofstream summary_file(summary_path);
        Json::StreamWriterBuilder writer;
        std::unique_ptr<Json::StreamWriter> json_writer(writer.newStreamWriter());
        json_writer->write(summary, &summary_file);
        summary_file.close();
        
        progress.report("complete", 100, 
            "Pipeline completed in " + std::to_string(duration.count()) + " seconds");
        
        return 0;
        
    } catch (const std::exception& e) {
        progress.report("error", -1, std::string("Pipeline failed: ") + e.what());
        return 1;
    }
}