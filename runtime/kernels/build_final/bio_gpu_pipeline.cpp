#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono>
#include <thread>
#include <fstream>
#include <getopt.h>
#include <cuda_runtime.h>
#include <iomanip>
#include <ctime>

namespace fs = std::filesystem;

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
    bool use_multi_gpu = false;
    int resistance_gpu = 1;
    int amr_gpu = 0;
};

void print_usage(const char* program_name) {
    std::cerr << "BioGPU Unified Pipeline v2.0" << std::endl;
    std::cerr << "Usage: " << program_name << " [OPTIONS]" << std::endl;
    std::cerr << "\nRequired arguments:" << std::endl;
    std::cerr << "  --r1 <path>              Forward reads (R1) FASTQ file" << std::endl;
    std::cerr << "  --r2 <path>              Reverse reads (R2) FASTQ file" << std::endl;
    std::cerr << "  --output-dir <path>      Output directory for results" << std::endl;
    std::cerr << "  --reference-db <path>    Microbial reference database path" << std::endl;
    std::cerr << "  --resistance-db <path>   Quinolone resistance database path" << std::endl;
    std::cerr << "  --sample-id <string>     Sample identifier" << std::endl;
    std::cerr << "  --gpu-device <id>        GPU device ID (default: 0)" << std::endl;
    std::cerr << "\nOptional arguments:" << std::endl;
    std::cerr << "  --use-multi-gpu          Enable multi-GPU mode" << std::endl;
    std::cerr << "  --progress-json          Output progress in JSON format" << std::endl;
    std::cerr << "  --threads <num>          Number of threads (default: 8)" << std::endl;
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
        {"sample-id", required_argument, 0, 0},
        {"progress-json", no_argument, 0, 0},
        {"threads", required_argument, 0, 0},
        {"help", no_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        if (c == 'h' || (c == 0 && std::string(long_options[option_index].name) == "help")) {
            print_usage(argv[0]);
            return false;
        }
        
        if (c == 0) {
            std::string opt_name = long_options[option_index].name;
            if (opt_name == "r1") options.r1_path = optarg;
            else if (opt_name == "r2") options.r2_path = optarg;
            else if (opt_name == "output-dir") options.output_dir = optarg;
            else if (opt_name == "reference-db") options.reference_db = optarg;
            else if (opt_name == "resistance-db") options.resistance_db = optarg;
            else if (opt_name == "gpu-device") options.gpu_device = std::stoi(optarg);
            else if (opt_name == "use-multi-gpu") options.use_multi_gpu = true;
            else if (opt_name == "sample-id") options.sample_id = optarg;
            else if (opt_name == "progress-json") options.progress_json = true;
            else if (opt_name == "threads") options.threads = std::stoi(optarg);
        }
    }
    
    if (options.r1_path.empty() || options.r2_path.empty() || 
        options.output_dir.empty() || options.reference_db.empty() || 
        options.resistance_db.empty() || options.sample_id.empty()) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        print_usage(argv[0]);
        return false;
    }
    
    return true;
}

void report_progress(const std::string& stage, int percentage, 
                    const std::string& message, bool json_output) {
    if (json_output) {
        std::cout << "{\"stage\":\"" << stage << "\",\"progress\":" << percentage 
                  << ",\"message\":\"" << message << "\"}" << std::endl;
    } else {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::cerr << "[" << std::put_time(std::localtime(&time_t), "%H:%M:%S") << "] "
                  << "[" << std::setw(3) << percentage << "%] " 
                  << stage << ": " << message << std::endl;
    }
}

int main(int argc, char* argv[]) {
    PipelineOptions options;
    
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    
    // Create output directory
    fs::create_directories(options.output_dir);
    
    report_progress("startup", 0, "Starting BioGPU unified pipeline", options.progress_json);
    
    // Check GPUs
    int device_count;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
        report_progress("error", -1, "Failed to query CUDA devices", options.progress_json);
        return 1;
    }
    
    if (!options.progress_json) {
        std::cout << "\nGPU Configuration:" << std::endl;
        for (int i = 0; i < device_count; i++) {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, i);
            std::cout << "  GPU " << i << ": " << prop.name 
                      << " (" << (prop.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
        }
        std::cout << std::endl;
    }
    
    report_progress("initialization", 10, "Validating input files", options.progress_json);
    
    // Validate input files
    if (!fs::exists(options.r1_path)) {
        report_progress("error", -1, "R1 file not found: " + options.r1_path, options.progress_json);
        return 1;
    }
    if (!fs::exists(options.r2_path)) {
        report_progress("error", -1, "R2 file not found: " + options.r2_path, options.progress_json);
        return 1;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Simulate pipeline stages
    report_progress("database_loading", 20, "Loading reference databases", options.progress_json);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    
    if (options.use_multi_gpu) {
        report_progress("multi_gpu_init", 30, 
            "Initializing multi-GPU: Resistance on GPU " + std::to_string(options.resistance_gpu) + 
            ", AMR on GPU " + std::to_string(options.amr_gpu), options.progress_json);
    }
    
    report_progress("fastq_processing", 40, "Processing FASTQ files", options.progress_json);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    
    report_progress("resistance_detection", 60, "Running resistance mutation detection", options.progress_json);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    
    report_progress("amr_gene_detection", 80, "Running AMR gene detection", options.progress_json);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    
    report_progress("report_generation", 90, "Generating clinical reports", options.progress_json);
    
    // Write example output files
    std::string summary_file = options.output_dir + "/" + options.sample_id + "_summary.txt";
    std::ofstream summary(summary_file);
    summary << "BioGPU Unified Pipeline Results" << std::endl;
    summary << "===============================" << std::endl;
    summary << "Sample ID: " << options.sample_id << std::endl;
    summary << "Input R1: " << options.r1_path << std::endl;
    summary << "Input R2: " << options.r2_path << std::endl;
    summary << "Multi-GPU: " << (options.use_multi_gpu ? "Enabled" : "Disabled") << std::endl;
    summary << std::endl;
    summary << "Pipeline Status: Framework Integrated Successfully" << std::endl;
    summary << std::endl;
    summary << "Next Steps:" << std::endl;
    summary << "1. Link with resistance/clean_resistance_pipeline_main.cpp" << std::endl;
    summary << "2. Link with genes/amr_detection_pipeline.cpp" << std::endl;
    summary << "3. Implement streaming FASTQ reader integration" << std::endl;
    summary.close();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    report_progress("complete", 100, 
        "Pipeline completed in " + std::to_string(duration.count()) + " seconds", 
        options.progress_json);
    
    if (!options.progress_json) {
        std::cout << "\nOutput written to: " << options.output_dir << std::endl;
        std::cout << "Summary file: " << summary_file << std::endl;
    }
    
    return 0;
}
