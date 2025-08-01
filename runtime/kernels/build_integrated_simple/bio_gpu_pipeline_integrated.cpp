#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <filesystem>
#include <chrono>
#include <thread>
#include <future>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <json/json.h>
#include <cuda_runtime.h>

namespace fs = std::filesystem;

// Command line options
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
    
    return true;
}

void report_progress(const std::string& stage, int percentage, const std::string& message, bool json_output) {
    if (json_output) {
        Json::Value root;
        root["stage"] = stage;
        root["progress"] = percentage;
        root["message"] = message;
        Json::FastWriter writer;
        std::cout << writer.write(root);
    } else {
        std::cerr << "[" << percentage << "%] " << stage << ": " << message << std::endl;
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
    cudaGetDeviceCount(&device_count);
    if (!options.progress_json) {
        std::cout << "Found " << device_count << " GPU(s):" << std::endl;
        for (int i = 0; i < device_count; i++) {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, i);
            std::cout << "  GPU " << i << ": " << prop.name 
                      << " (" << (prop.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
        }
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // For now, just demonstrate the framework
    report_progress("initialization", 10, "Pipeline framework initialized", options.progress_json);
    
    if (options.use_multi_gpu) {
        report_progress("multi_gpu", 20, "Multi-GPU mode: Resistance on GPU " + 
                       std::to_string(options.resistance_gpu) + ", AMR on GPU " + 
                       std::to_string(options.amr_gpu), options.progress_json);
    }
    
    // Simulate processing
    report_progress("processing", 50, "Processing " + options.sample_id, options.progress_json);
    
    // In a real implementation, we would:
    // 1. Launch resistance pipeline on GPU 1
    // 2. Launch AMR pipeline on GPU 0
    // 3. Use streaming FASTQ reader to feed both
    // 4. Collect and merge results
    
    std::this_thread::sleep_for(std::chrono::seconds(1));
    
    // Write summary
    Json::Value summary;
    summary["sample_id"] = options.sample_id;
    summary["status"] = "framework_test";
    summary["output_directory"] = options.output_dir;
    
    std::string summary_path = options.output_dir + "/" + options.sample_id + "_summary.json";
    std::ofstream summary_file(summary_path);
    Json::FastWriter writer;
    summary_file << writer.write(summary);
    summary_file.close();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    report_progress("complete", 100, "Pipeline completed in " + 
                   std::to_string(duration.count()) + " seconds", options.progress_json);
    
    if (!options.progress_json) {
        std::cout << "\nIntegrated pipeline framework is working!" << std::endl;
        std::cout << "Next step: Link with actual resistance and AMR pipelines" << std::endl;
    }
    
    return 0;
}
