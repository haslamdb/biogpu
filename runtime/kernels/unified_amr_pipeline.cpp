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

// Include unified batch processor
#include "shared/unified_batch_processor.h"
// Include actual pipeline implementations
#include "resistance/clean_resistance_pipeline_main.cpp"  // Contains CleanResistancePipeline class
#include "genes/amr_detection_pipeline.h"

namespace fs = std::filesystem;

// Progress reporting structure
struct ProgressReporter {
    std::string current_stage;
    int progress_percentage;
    std::string message;
    bool json_output;
    
    ProgressReporter(bool use_json) : json_output(use_json), progress_percentage(0) {}
    
    void report(const std::string& stage, int percentage, const std::string& msg) {
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
    int resistance_gpu = 1;  // Default to A5000 (24GB) - less memory intensive
    int amr_gpu = 0;        // Default to A6000 (48GB) - more memory intensive
    bool use_multi_gpu = false;
    
    // Advanced options
    bool use_bloom_filter = false;
    double min_identity = 0.90;
    double min_coverage = 0.80;
    int batch_size = 50000;  // Balanced between the two pipelines
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
    std::cerr << "  --gpu-device <id>        GPU device ID (default: 0) - single GPU mode" << std::endl;
    std::cerr << "  --use-multi-gpu          Enable multi-GPU mode (uses both GPUs)" << std::endl;
    std::cerr << "  --resistance-gpu <id>    GPU for resistance pipeline (default: 1 - A5000)" << std::endl;
    std::cerr << "  --amr-gpu <id>           GPU for AMR genes pipeline (default: 0 - A6000)" << std::endl;
    std::cerr << "  --progress-json          Output progress in JSON format" << std::endl;
    std::cerr << "  --threads <num>          Number of CPU threads (default: 8)" << std::endl;
    std::cerr << "  --use-bloom-filter       Enable bloom filter pre-screening" << std::endl;
    std::cerr << "  --batch-size <num>       Reads per batch (default: 50000)" << std::endl;
    std::cerr << "  --min-identity <float>   Minimum identity threshold (default: 0.90)" << std::endl;
    std::cerr << "  --min-coverage <float>   Minimum coverage threshold (default: 0.80)" << std::endl;
    std::cerr << "  --min-allele-depth <num> Minimum depth for allele frequency (default: 5)" << std::endl;
    std::cerr << "  --disable-resistance     Skip resistance mutation detection" << std::endl;
    std::cerr << "  --disable-amr-genes      Skip AMR gene detection" << std::endl;
    std::cerr << "  --enable-em              Enable EM algorithm for multi-mapping reads" << std::endl;
    std::cerr << "  --em-iterations <num>    Number of EM iterations (default: 10)" << std::endl;
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
        {"use-bloom-filter", no_argument, 0, 0},
        {"batch-size", required_argument, 0, 0},
        {"min-identity", required_argument, 0, 0},
        {"min-coverage", required_argument, 0, 0},
        {"min-allele-depth", required_argument, 0, 0},
        {"disable-resistance", no_argument, 0, 0},
        {"disable-amr-genes", no_argument, 0, 0},
        {"enable-em", no_argument, 0, 0},
        {"em-iterations", required_argument, 0, 0},
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
            else if (opt_name == "use-bloom-filter") options.use_bloom_filter = true;
            else if (opt_name == "batch-size") options.batch_size = std::stoi(optarg);
            else if (opt_name == "min-identity") options.min_identity = std::stod(optarg);
            else if (opt_name == "min-coverage") options.min_coverage = std::stod(optarg);
            else if (opt_name == "min-allele-depth") options.min_allele_depth = std::stoi(optarg);
            else if (opt_name == "disable-resistance") options.run_resistance = false;
            else if (opt_name == "disable-amr-genes") options.run_amr_genes = false;
            else if (opt_name == "enable-em") options.enable_em = true;
            else if (opt_name == "em-iterations") options.em_iterations = std::stoi(optarg);
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
    
    // Validate files exist
    if (!fs::exists(options.r1_path)) {
        std::cerr << "Error: R1 file does not exist: " << options.r1_path << std::endl;
        return false;
    }
    if (!fs::exists(options.r2_path)) {
        std::cerr << "Error: R2 file does not exist: " << options.r2_path << std::endl;
        return false;
    }
    if (!fs::exists(options.reference_db)) {
        std::cerr << "Error: Reference database does not exist: " << options.reference_db << std::endl;
        return false;
    }
    if (!fs::exists(options.resistance_db)) {
        std::cerr << "Error: Resistance database does not exist: " << options.resistance_db << std::endl;
        return false;
    }
    
    // Create output directory if it doesn't exist
    fs::create_directories(options.output_dir);
    
    return true;
}

// Unified FASTQ processing function using the batch processor
void process_fastq_unified(const PipelineOptions& options, ProgressReporter& progress) {
    using namespace BioGPU;
    
    progress.report("initialization", 5, "Initializing unified pipeline");
    
    // Create unified batch processor
    UnifiedBatchProcessor processor(options.batch_size, 10);  // 10 batches queue depth
    
    // Set up progress callback
    auto progress_callback = [&progress](const std::string& stage, int percentage, const std::string& message) {
        progress.report(stage, percentage, message);
    };
    
    // Create and configure resistance detector
    if (options.run_resistance) {
        progress.report("setup", 10, "Setting up resistance detector");
        
        auto resistance_detector = std::make_shared<ResistanceDetector>();
        
        // Configure resistance-specific settings
        ResistanceConfig res_config;
        res_config.gpu_device = options.use_multi_gpu ? options.resistance_gpu : options.gpu_device;
        res_config.batch_size = options.batch_size;
        res_config.use_bloom_filter = options.use_bloom_filter;
        res_config.min_identity = options.min_identity;
        res_config.min_coverage = options.min_coverage;
        res_config.threads = options.threads;
        res_config.output_dir = options.output_dir;
        res_config.sample_id = options.sample_id;
        res_config.nucleotide_index_path = options.reference_db + "/nucleotide_index";
        res_config.protein_db_path = options.resistance_db + "/protein_db";
        res_config.fq_csv_path = options.resistance_db + "/fq_mutations.csv";
        res_config.enable_sw_alignment = true;
        res_config.min_allele_depth = options.min_allele_depth;
        
        resistance_detector->setResistanceConfig(res_config);
        resistance_detector->setProgressCallback(progress_callback);
        processor.addDetector(resistance_detector);
    }
    
    // Create and configure AMR gene detector
    if (options.run_amr_genes) {
        progress.report("setup", 15, "Setting up AMR gene detector");
        
        auto amr_detector = std::make_shared<AMRGeneDetector>();
        
        // Configure AMR-specific settings
        AMRGeneConfig amr_config;
        amr_config.gpu_device = options.use_multi_gpu ? options.amr_gpu : options.gpu_device;
        amr_config.batch_size = options.batch_size;
        amr_config.use_bloom_filter = options.use_bloom_filter;
        amr_config.min_identity = options.min_identity;
        amr_config.min_coverage = options.min_coverage;
        amr_config.threads = options.threads;
        amr_config.output_dir = options.output_dir;
        amr_config.sample_id = options.sample_id;
        amr_config.amr_db_path = options.reference_db + "/amr_genes";
        amr_config.protein_db_path = options.reference_db + "/amr_proteins";
        amr_config.enable_em_algorithm = options.enable_em;
        amr_config.em_iterations = options.em_iterations;
        
        amr_detector->setAMRConfig(amr_config);
        amr_detector->setProgressCallback(progress_callback);
        processor.addDetector(amr_detector);
    }
    
    progress.report("processing", 20, "Starting FASTQ processing");
    
    // Process the FASTQ files
    bool success = false;
    if (!options.r2_path.empty()) {
        // Paired-end processing
        success = processor.processPairedEnd(options.r1_path, options.r2_path, progress_callback);
    } else {
        // Single-end processing
        success = processor.processSingleEnd(options.r1_path, progress_callback);
    }
    
    if (!success) {
        throw std::runtime_error("Failed to process FASTQ files");
    }
    
    progress.report("output", 90, "Writing output files");
    
    // Write all output files
    processor.writeAllOutputs(options.output_dir);
    
    // Get combined results
    Json::Value combined_results;
    processor.getAllResults(combined_results);
    
    // Write combined summary
    std::string summary_path = fs::path(options.output_dir) / (options.sample_id + "_unified_summary.json");
    std::ofstream summary_file(summary_path);
    Json::StreamWriterBuilder writer;
    std::unique_ptr<Json::StreamWriter> json_writer(writer.newStreamWriter());
    json_writer->write(combined_results, &summary_file);
    summary_file.close();
    
    progress.report("complete", 100, 
        "Pipeline completed. Processed " + std::to_string(processor.getTotalReadsProcessed()) + 
        " reads (" + std::to_string(processor.getTotalBasesProcessed()) + " bases)");
}

int main(int argc, char* argv[]) {
    PipelineOptions options;
    
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    
    // Validate GPU configuration
    int device_count;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
        std::cerr << "Error: Failed to get GPU device count" << std::endl;
        return 1;
    }
    
    if (options.use_multi_gpu) {
        // Multi-GPU mode validation
        if (options.resistance_gpu >= device_count || options.amr_gpu >= device_count) {
            std::cerr << "Error: Invalid GPU device IDs. Available GPUs: 0-" << (device_count-1) << std::endl;
            return 1;
        }
        if (options.resistance_gpu == options.amr_gpu) {
            std::cerr << "Warning: Both pipelines assigned to same GPU. Consider using different GPUs for better performance." << std::endl;
        }
        
        // Log GPU configuration
        cudaDeviceProp prop_resistance, prop_amr;
        cudaGetDeviceProperties(&prop_resistance, options.resistance_gpu);
        cudaGetDeviceProperties(&prop_amr, options.amr_gpu);
        
        if (!options.progress_json) {
            std::cerr << "Multi-GPU Configuration:" << std::endl;
            std::cerr << "  Resistance pipeline: GPU " << options.resistance_gpu 
                      << " (" << prop_resistance.name << ", " 
                      << (prop_resistance.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
            std::cerr << "  AMR genes pipeline: GPU " << options.amr_gpu 
                      << " (" << prop_amr.name << ", " 
                      << (prop_amr.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
        }
    } else {
        // Single GPU mode
        if (options.gpu_device >= device_count) {
            std::cerr << "Error: Invalid GPU device ID. Available GPUs: 0-" << (device_count-1) << std::endl;
            return 1;
        }
        if (cudaSetDevice(options.gpu_device) != cudaSuccess) {
            std::cerr << "Error: Failed to set GPU device " << options.gpu_device << std::endl;
            return 1;
        }
    }
    
    // Initialize progress reporter
    ProgressReporter progress(options.progress_json);
    
    progress.report("startup", 0, "Starting BioGPU unified pipeline");
    
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Run the unified processing
        process_fastq_unified(options, progress);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Output final summary
        Json::Value summary;
        summary["sample_id"] = options.sample_id;
        summary["processing_time_seconds"] = duration.count();
        summary["resistance_enabled"] = options.run_resistance;
        summary["amr_genes_enabled"] = options.run_amr_genes;
        summary["output_directory"] = options.output_dir;
        
        std::string summary_path = fs::path(options.output_dir) / (options.sample_id + "_summary.json");
        std::ofstream summary_file(summary_path);
        Json::StreamWriterBuilder writer;
        std::unique_ptr<Json::StreamWriter> json_writer(writer.newStreamWriter());
        json_writer->write(summary, &summary_file);
        summary_file.close();
        
        return 0;
        
    } catch (const std::exception& e) {
        progress.report("error", -1, std::string("Pipeline failed: ") + e.what());
        return 1;
    }
}