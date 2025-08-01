#include <iostream>
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include <getopt.h>

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
    std::cerr << "  --help                   Show this help message" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "BioGPU Unified Pipeline v1.0 (Minimal Build)" << std::endl;
    
    // Parse arguments
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"r1", required_argument, 0, 0},
        {"r2", required_argument, 0, 0},
        {"output-dir", required_argument, 0, 0},
        {"reference-db", required_argument, 0, 0},
        {"resistance-db", required_argument, 0, 0},
        {"sample-id", required_argument, 0, 0},
        {"gpu-device", required_argument, 0, 0},
        {"use-multi-gpu", no_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    std::string r1_path, r2_path, output_dir, reference_db, resistance_db, sample_id;
    int gpu_device = 0;
    bool use_multi_gpu = false;
    
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        if (c == 'h') {
            print_usage(argv[0]);
            return 0;
        }
        
        if (c == 0) {
            std::string opt_name = long_options[option_index].name;
            if (opt_name == "r1") r1_path = optarg;
            else if (opt_name == "r2") r2_path = optarg;
            else if (opt_name == "output-dir") output_dir = optarg;
            else if (opt_name == "reference-db") reference_db = optarg;
            else if (opt_name == "resistance-db") resistance_db = optarg;
            else if (opt_name == "sample-id") sample_id = optarg;
            else if (opt_name == "gpu-device") gpu_device = std::stoi(optarg);
            else if (opt_name == "use-multi-gpu") use_multi_gpu = true;
        }
    }
    
    // Check GPUs
    int device_count;
    cudaGetDeviceCount(&device_count);
    std::cout << "Found " << device_count << " GPU(s)" << std::endl;
    
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        std::cout << "  GPU " << i << ": " << prop.name 
                  << " (" << (prop.totalGlobalMem / (1024*1024*1024)) << " GB)" << std::endl;
    }
    
    std::cout << "\nPipeline Configuration:" << std::endl;
    std::cout << "  R1: " << r1_path << std::endl;
    std::cout << "  R2: " << r2_path << std::endl;
    std::cout << "  Output: " << output_dir << std::endl;
    std::cout << "  Sample ID: " << sample_id << std::endl;
    
    if (use_multi_gpu) {
        std::cout << "  Multi-GPU: ENABLED" << std::endl;
        std::cout << "    Resistance: GPU 1 (RTX A5000)" << std::endl;
        std::cout << "    AMR Genes: GPU 0 (RTX A6000)" << std::endl;
    } else {
        std::cout << "  GPU Device: " << gpu_device << std::endl;
    }
    
    std::cout << "\nThis is a minimal build for testing the framework." << std::endl;
    std::cout << "Full pipeline integration will be added in the next step." << std::endl;
    
    return 0;
}
