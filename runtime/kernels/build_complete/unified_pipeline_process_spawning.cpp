#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <chrono>
#include <thread>
#include <future>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include <cuda_runtime.h>
#include <zlib.h>
#include <cstring>
#include <iomanip>
#include <ctime>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

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
    bool run_resistance = true;
    bool run_amr = true;
    bool use_bloom_filter = false;
    float min_identity = 0.90f;
    float min_coverage = 0.80f;
    int batch_size = 50000;
};

// Progress reporter
class ProgressReporter {
private:
    bool json_output;
    std::mutex mutex;
    
public:
    ProgressReporter(bool use_json) : json_output(use_json) {}
    
    void report(const std::string& stage, int percentage, const std::string& message) {
        std::lock_guard<std::mutex> lock(mutex);
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
};

// Process spawner for running external pipelines
class ProcessSpawner {
private:
    std::string executable_path;
    std::vector<std::string> args;
    int gpu_device;
    bool capture_output;
    
public:
    ProcessSpawner(const std::string& path, int gpu = 0, bool capture = false) 
        : executable_path(path), gpu_device(gpu), capture_output(capture) {}
    
    void addArg(const std::string& arg) {
        args.push_back(arg);
    }
    
    void addArg(const std::string& flag, const std::string& value) {
        args.push_back(flag);
        args.push_back(value);
    }
    
    int run(std::string* output = nullptr) {
        // Prepare arguments
        std::vector<char*> c_args;
        c_args.push_back(const_cast<char*>(executable_path.c_str()));
        for (auto& arg : args) {
            c_args.push_back(const_cast<char*>(arg.c_str()));
        }
        c_args.push_back(nullptr);
        
        // Create pipes for output capture if needed
        int pipefd[2];
        if (capture_output && pipe(pipefd) == -1) {
            std::cerr << "Failed to create pipe" << std::endl;
            return -1;
        }
        
        pid_t pid = fork();
        if (pid == -1) {
            std::cerr << "Failed to fork process" << std::endl;
            return -1;
        }
        
        if (pid == 0) {
            // Child process
            // Set CUDA_VISIBLE_DEVICES
            std::string cuda_env = "CUDA_VISIBLE_DEVICES=" + std::to_string(gpu_device);
            putenv(const_cast<char*>(cuda_env.c_str()));
            
            // Redirect output if capturing
            if (capture_output) {
                close(pipefd[0]);  // Close read end
                dup2(pipefd[1], STDOUT_FILENO);
                dup2(pipefd[1], STDERR_FILENO);
                close(pipefd[1]);
            }
            
            // Execute the program
            execv(executable_path.c_str(), c_args.data());
            
            // If we get here, exec failed
            std::cerr << "Failed to execute: " << executable_path << std::endl;
            exit(1);
        }
        
        // Parent process
        if (capture_output) {
            close(pipefd[1]);  // Close write end
            
            // Read output
            if (output) {
                char buffer[4096];
                ssize_t n;
                while ((n = read(pipefd[0], buffer, sizeof(buffer) - 1)) > 0) {
                    buffer[n] = '\0';
                    *output += buffer;
                }
            }
            close(pipefd[0]);
        }
        
        // Wait for child to complete
        int status;
        waitpid(pid, &status, 0);
        
        if (WIFEXITED(status)) {
            return WEXITSTATUS(status);
        }
        return -1;
    }
};

void print_usage(const char* program_name) {
    std::cerr << "BioGPU Unified Pipeline v4.0 (Process Integration)" << std::endl;
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
    std::cerr << "  --use-bloom-filter       Enable bloom filter pre-screening" << std::endl;
    std::cerr << "  --batch-size <num>       Batch size (default: 50000)" << std::endl;
    std::cerr << "  --disable-resistance     Skip resistance detection" << std::endl;
    std::cerr << "  --disable-amr            Skip AMR gene detection" << std::endl;
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
        {"use-bloom-filter", no_argument, 0, 0},
        {"batch-size", required_argument, 0, 0},
        {"disable-resistance", no_argument, 0, 0},
        {"disable-amr", no_argument, 0, 0},
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
            else if (opt_name == "use-bloom-filter") options.use_bloom_filter = true;
            else if (opt_name == "batch-size") options.batch_size = std::stoi(optarg);
            else if (opt_name == "disable-resistance") options.run_resistance = false;
            else if (opt_name == "disable-amr") options.run_amr = false;
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

// Find executable in common locations
std::string find_executable(const std::string& name, const std::vector<std::string>& search_paths) {
    for (const auto& path : search_paths) {
        std::string full_path = path + "/" + name;
        if (fs::exists(full_path) && fs::is_regular_file(full_path)) {
            // Check if executable
            if (access(full_path.c_str(), X_OK) == 0) {
                return full_path;
            }
        }
    }
    return "";
}

// Run resistance pipeline
void run_resistance_pipeline(const PipelineOptions& options, ProgressReporter& progress) {
    if (options.use_multi_gpu) {
        cudaSetDevice(options.resistance_gpu);
    }
    
    progress.report("resistance_init", 5, "Initializing resistance detection");
    
    // Find resistance pipeline executable
    std::vector<std::string> search_paths = {
        "../resistance",
        "../resistance/build",
        "./resistance",
        "/usr/local/bin"
    };
    
    std::string exec_path = find_executable("clean_resistance_pipeline", search_paths);
    if (exec_path.empty()) {
        exec_path = find_executable("resistance_pipeline", search_paths);
    }
    
    if (exec_path.empty()) {
        progress.report("resistance_error", -1, "Resistance pipeline executable not found");
        return;
    }
    
    progress.report("resistance_found", 10, "Found resistance pipeline at: " + exec_path);
    
    // Prepare resistance output directory
    std::string resistance_output = options.output_dir + "/resistance";
    fs::create_directories(resistance_output);
    
    // Build command line for resistance pipeline
    // The resistance pipeline expects: <nucleotide_index> <protein_db> <reads_r1> <reads_r2> [options]
    ProcessSpawner resistance_proc(exec_path, options.resistance_gpu);
    
    // For now, use the fq_resistance_index as both nucleotide and protein db
    std::string resistance_index = "/home/david/projects/biogpu/data/fq_resistance_index";
    std::string protein_db = "/home/david/projects/biogpu/data/protein_resistance_db";
    
    // Add positional arguments first
    resistance_proc.addArg(resistance_index);
    resistance_proc.addArg(protein_db);
    resistance_proc.addArg(options.r1_path);
    resistance_proc.addArg(options.r2_path);
    
    // Add optional arguments
    resistance_proc.addArg("--output-prefix", resistance_output + "/" + options.sample_id);
    resistance_proc.addArg("--fq-csv", options.resistance_db);
    
    if (!options.use_bloom_filter) {
        resistance_proc.addArg("--no-bloom");
    }
    
    progress.report("resistance_processing", 50, "Running resistance mutation detection");
    
    // Execute resistance pipeline
    std::string output;
    int ret = resistance_proc.run(&output);
    
    if (ret == 0) {
        progress.report("resistance_complete", 100, "Resistance detection completed successfully");
    } else {
        progress.report("resistance_error", -1, "Resistance detection failed with code: " + std::to_string(ret));
    }
}

// Run AMR pipeline
void run_amr_pipeline(const PipelineOptions& options, ProgressReporter& progress) {
    if (options.use_multi_gpu) {
        cudaSetDevice(options.amr_gpu);
    }
    
    progress.report("amr_init", 5, "Initializing AMR gene detection");
    
    // Find AMR pipeline executable
    std::vector<std::string> search_paths = {
        "../genes",
        "../genes/build",
        "./genes",
        "/usr/local/bin"
    };
    
    std::string exec_path = find_executable("amr_detection", search_paths);
    
    if (exec_path.empty()) {
        progress.report("amr_error", -1, "AMR detection executable not found");
        return;
    }
    
    progress.report("amr_found", 10, "Found AMR pipeline at: " + exec_path);
    
    // Prepare AMR output directory
    std::string amr_output = options.output_dir + "/amr";
    fs::create_directories(amr_output);
    
    // Build command line for AMR pipeline
    // AMR detection expects: <amr_db_path> <input_csv> <output_dir> [options]
    // First create a temporary CSV file for the AMR pipeline
    std::string temp_csv = amr_output + "/" + options.sample_id + "_input.csv";
    std::ofstream csv_file(temp_csv);
    csv_file << "sample_name,fastq_path" << std::endl;
    csv_file << options.sample_id << "," << options.r1_path;
    if (!options.r2_path.empty()) {
        csv_file << "," << options.r2_path;
    }
    csv_file << std::endl;
    csv_file.close();
    
    ProcessSpawner amr_proc(exec_path, options.amr_gpu);
    // AMR pipeline expects database in format: dna.fasta,protein.fasta
    std::string amr_db_path = "/home/david/projects/biogpu/data/AMR_CDS.fa,/home/david/projects/biogpu/data/AMR_protein.fa";
    
    // Positional arguments
    amr_proc.addArg(amr_db_path);          // AMR database path
    amr_proc.addArg(temp_csv);              // Input CSV
    amr_proc.addArg(amr_output);            // Output directory
    
    // Optional arguments
    amr_proc.addArg("--min-identity");
    amr_proc.addArg(std::to_string(options.min_identity));
    amr_proc.addArg("--min-coverage");
    amr_proc.addArg(std::to_string(options.min_coverage));
    
    if (options.use_bloom_filter) {
        amr_proc.addArg("--use-bloom-filter");
    }
    
    progress.report("amr_processing", 50, "Running AMR gene detection");
    
    // Execute AMR pipeline
    std::string output;
    int ret = amr_proc.run(&output);
    
    if (ret == 0) {
        progress.report("amr_complete", 100, "AMR gene detection completed successfully");
    } else {
        progress.report("amr_error", -1, "AMR detection failed with code: " + std::to_string(ret));
    }
}

int main(int argc, char* argv[]) {
    PipelineOptions options;
    
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    
    // Create output directory
    fs::create_directories(options.output_dir);
    
    ProgressReporter progress(options.progress_json);
    progress.report("startup", 0, "Starting BioGPU unified pipeline");
    
    // Check GPUs
    int device_count;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
        progress.report("error", -1, "Failed to query CUDA devices");
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
    
    // Validate input files
    if (!fs::exists(options.r1_path) || !fs::exists(options.r2_path)) {
        progress.report("error", -1, "Input files not found");
        return 1;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        progress.report("initialization", 10, "Validating pipelines");
        
        // Run pipelines
        if (options.use_multi_gpu && options.run_resistance && options.run_amr) {
            progress.report("multi_gpu", 20, "Running pipelines in parallel on multiple GPUs");
            
            auto resistance_future = std::async(std::launch::async,
                run_resistance_pipeline, std::ref(options), std::ref(progress));
            
            auto amr_future = std::async(std::launch::async,
                run_amr_pipeline, std::ref(options), std::ref(progress));
            
            resistance_future.get();
            amr_future.get();
            
        } else {
            if (options.run_resistance) {
                run_resistance_pipeline(options, progress);
            }
            if (options.run_amr) {
                run_amr_pipeline(options, progress);
            }
        }
        
        // Generate unified summary
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::string summary_file = options.output_dir + "/" + options.sample_id + "_unified_summary.txt";
        std::ofstream summary(summary_file);
        summary << "BioGPU Unified Pipeline Summary" << std::endl;
        summary << "==============================" << std::endl;
        summary << "Sample ID: " << options.sample_id << std::endl;
        summary << "Processing time: " << duration.count() << " seconds" << std::endl;
        summary << "Multi-GPU: " << (options.use_multi_gpu ? "Enabled" : "Disabled") << std::endl;
        summary << "Resistance: " << (options.run_resistance ? "Completed" : "Skipped") << std::endl;
        summary << "AMR Genes: " << (options.run_amr ? "Completed" : "Skipped") << std::endl;
        summary << "\nResults:" << std::endl;
        
        // Check for resistance results
        if (options.run_resistance) {
            std::string res_clinical = options.output_dir + "/resistance/" + options.sample_id + "_clinical_report.txt";
            if (fs::exists(res_clinical)) {
                summary << "- Resistance clinical report: " << res_clinical << std::endl;
            }
        }
        
        // Check for AMR results
        if (options.run_amr) {
            std::string amr_tsv = options.output_dir + "/amr/" + options.sample_id + "_amr_genes.tsv";
            if (fs::exists(amr_tsv)) {
                summary << "- AMR genes report: " << amr_tsv << std::endl;
            }
        }
        
        summary.close();
        
        progress.report("complete", 100, "Pipeline completed in " + std::to_string(duration.count()) + " seconds");
        
    } catch (const std::exception& e) {
        progress.report("error", -1, std::string("Pipeline failed: ") + e.what());
        return 1;
    }
    
    return 0;
}