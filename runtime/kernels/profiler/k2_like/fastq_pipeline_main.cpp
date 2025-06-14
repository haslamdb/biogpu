// fastq_pipeline_main.cpp
// Main program for FASTQ processing pipeline
#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include "fastq_processing.h"

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <fastq_file> [options]\n";
    std::cerr << "\nOptions:\n";
    std::cerr << "  --threads <n>    Number of GPU processing threads (default: 2)\n";
    std::cerr << "  --batch <n>      Batch size for processing (default: 10000)\n";
    std::cerr << "  --k <n>          K-mer size (default: 31)\n";
    std::cerr << "  --m <n>          Minimizer window size (default: 15)\n";
    std::cerr << "  --quiet          Suppress progress messages\n";
    std::cerr << "  --help           Show this help message\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string fastq_file = argv[1];
    int num_threads = 2;
    size_t batch_size = 10000;
    int k = 31;
    int m = 15;
    bool quiet = false;
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::atoi(argv[++i]);
            if (num_threads < 1 || num_threads > 16) {
                std::cerr << "Error: threads must be between 1 and 16\n";
                return 1;
            }
        } else if (arg == "--batch" && i + 1 < argc) {
            batch_size = std::atoi(argv[++i]);
            if (batch_size < 100 || batch_size > 100000) {
                std::cerr << "Error: batch size must be between 100 and 100000\n";
                return 1;
            }
        } else if (arg == "--k" && i + 1 < argc) {
            k = std::atoi(argv[++i]);
            if (k < 15 || k > 31) {
                std::cerr << "Error: k must be between 15 and 31\n";
                return 1;
            }
        } else if (arg == "--m" && i + 1 < argc) {
            m = std::atoi(argv[++i]);
            if (m < 10 || m > 25) {
                std::cerr << "Error: m must be between 10 and 25\n";
                return 1;
            }
        } else if (arg == "--quiet") {
            quiet = true;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    try {
        // Create pipeline with specified parameters
        biogpu::GPUMinimizerPipeline pipeline(k, m, batch_size, num_threads);
        
        // Create statistics collector
        biogpu::MinimizerStats stats;
        
        // Process file
        if (!quiet) {
            std::cout << "Processing " << fastq_file << "...\n";
            std::cout << "Parameters:\n";
            std::cout << "  K-mer size: " << k << "\n";
            std::cout << "  Minimizer window: " << m << "\n";
            std::cout << "  Batch size: " << batch_size << "\n";
            std::cout << "  GPU threads: " << num_threads << "\n\n";
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Process with progress updates
        size_t processed_reads = 0;
        pipeline.process_file(fastq_file, 
            [&stats, &processed_reads, quiet](const ReadBatch& batch, 
                                             const std::vector<std::vector<Minimizer>>& minimizers) {
                stats.process_batch(batch, minimizers);
                processed_reads += batch.size();
                
                if (!quiet && processed_reads % 100000 == 0) {
                    std::cout << "Processed " << processed_reads << " reads...\n";
                }
            });
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        // Print results
        stats.print_summary();
        
        auto pipeline_stats = pipeline.get_statistics();
        std::cout << "\n=== Performance Metrics ===\n";
        std::cout << "Total processing time: " << duration.count() / 1000.0 << " seconds\n";
        std::cout << "Reads per second: " << pipeline_stats.total_reads * 1000.0 / duration.count() << "\n";
        std::cout << "Batches processed: " << pipeline_stats.total_batches << "\n";
        
        // Calculate bases processed
        size_t total_bases = 0;
        size_t avg_read_length = 150;  // Estimate, could calculate from actual data
        total_bases = pipeline_stats.total_reads * avg_read_length;
        std::cout << "Estimated throughput: " 
                  << (total_bases / 1000000.0) / (duration.count() / 1000.0) 
                  << " Mbases/second\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}