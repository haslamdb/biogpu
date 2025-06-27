// test_streaming_minimizer_extraction.cu
// Test program to combine StreamingFnaProcessor with GPU minimizer extraction
// This will help test and debug the minimizer extraction pipeline

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <filesystem>
#include <cuda_runtime.h>

// Include the necessary headers
#include "processing/genome_file_processor.h"
#include "gpu/gpu_database_kernels.h"
#include "gpu/gpu_minimizer_extraction.cuh"
#include "memory/gpu_memory_manager.h"
#include "gpu_kraken_types.h"

// Test configuration
struct TestConfig {
    std::string input_fna_path = "data/test_50_genomes.fna";
    std::string temp_dir = "/tmp/kraken_test";
    size_t batch_size = 10;  // Process 10 genomes at a time
    
    // Minimizer parameters
    uint32_t k_value = 31;
    uint32_t ell_value = 31;
    uint32_t spaces_value = 7;
    uint64_t xor_mask = 0x3c8bfbb395c60474ULL;
    
    // Memory configuration
    size_t max_minimizers_per_batch = 5000000;  // 5M minimizers per batch
    size_t sequence_buffer_size = 100 * 1024 * 1024;  // 100MB for sequences
};

// Helper function to print minimizer statistics
void print_minimizer_stats(const std::vector<GPUMinimizerHit>& minimizers) {
    if (minimizers.empty()) {
        std::cout << "No minimizers found!" << std::endl;
        return;
    }
    
    // Count unique minimizers
    std::unordered_set<uint64_t> unique_hashes;
    std::unordered_map<uint32_t, uint32_t> taxon_counts;
    
    for (const auto& hit : minimizers) {
        unique_hashes.insert(hit.minimizer_hash);
        taxon_counts[hit.taxon_id]++;
    }
    
    std::cout << "\nMinimizer Statistics:" << std::endl;
    std::cout << "  Total minimizers: " << minimizers.size() << std::endl;
    std::cout << "  Unique minimizers: " << unique_hashes.size() << std::endl;
    std::cout << "  Unique taxa: " << taxon_counts.size() << std::endl;
    
    // Show first few minimizers for debugging
    std::cout << "\nFirst 5 minimizers:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), minimizers.size()); i++) {
        const auto& hit = minimizers[i];
        std::cout << "  [" << i << "] hash=" << std::hex << hit.minimizer_hash << std::dec
                  << ", genome=" << hit.genome_id 
                  << ", pos=" << hit.position
                  << ", taxon=" << hit.taxon_id << std::endl;
    }
}

// Main test function
bool test_streaming_minimizer_extraction(const TestConfig& config) {
    std::cout << "=== Testing Streaming FNA Processing with GPU Minimizer Extraction ===" << std::endl;
    std::cout << "Input file: " << config.input_fna_path << std::endl;
    std::cout << "Batch size: " << config.batch_size << " genomes" << std::endl;
    std::cout << "K-mer size: " << config.k_value << std::endl;
    
    // Initialize GPU memory manager
    MemoryConfig mem_config;
    mem_config.minimizer_capacity = config.max_minimizers_per_batch;
    mem_config.sequence_batch_size = config.batch_size;
    
    GPUMemoryManager memory_manager(mem_config);
    if (!memory_manager.initialize()) {
        std::cerr << "Failed to initialize GPU memory manager" << std::endl;
        return false;
    }
    
    // Allocate GPU memory
    if (!memory_manager.allocate_sequence_memory(config.batch_size, config.sequence_buffer_size)) {
        std::cerr << "Failed to allocate sequence memory" << std::endl;
        return false;
    }
    
    if (!memory_manager.allocate_minimizer_memory(config.max_minimizers_per_batch)) {
        std::cerr << "Failed to allocate minimizer memory" << std::endl;
        return false;
    }
    
    // Initialize streaming processor
    StreamingFnaProcessor streaming_processor(
        config.input_fna_path, 
        config.temp_dir, 
        config.batch_size
    );
    
    // Process batches
    int batch_number = 0;
    size_t total_minimizers_extracted = 0;
    std::vector<GPUMinimizerHit> all_minimizers;  // Collect for analysis
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    std::vector<std::string> batch_files;
    std::vector<uint32_t> batch_taxons;
    
    while (streaming_processor.process_next_batch(batch_files, batch_taxons)) {
        batch_number++;
        std::cout << "\n--- Processing Batch " << batch_number << " ---" << std::endl;
        std::cout << "Genomes in batch: " << batch_files.size() << std::endl;
        
        // Load sequences from batch files
        std::vector<GPUGenomeInfo> genome_infos;
        std::string concatenated_sequences;
        size_t current_offset = 0;
        
        GenomeFileProcessor file_processor;
        
        for (size_t i = 0; i < batch_files.size(); i++) {
            std::cout << "Loading genome " << i << ": " << batch_files[i] 
                      << " (taxon=" << batch_taxons[i] << ")" << std::endl;
            
            auto sequences = file_processor.load_sequences_from_fasta(batch_files[i]);
            
            for (const auto& seq : sequences) {
                GPUGenomeInfo info;
                info.genome_id = genome_infos.size();
                info.sequence_offset = current_offset;
                info.sequence_length = seq.length();
                info.taxon_id = batch_taxons[i];
                info.minimizer_count = 0;
                
                genome_infos.push_back(info);
                concatenated_sequences += seq + '\0';
                current_offset += seq.length() + 1;
                
                std::cout << "  Sequence " << info.genome_id 
                          << ": length=" << seq.length() << " bp" << std::endl;
            }
        }
        
        if (genome_infos.empty()) {
            std::cout << "No sequences loaded from batch, skipping..." << std::endl;
            continue;
        }
        
        std::cout << "Total sequences in batch: " << genome_infos.size() << std::endl;
        std::cout << "Total sequence data: " << concatenated_sequences.length() << " bytes" << std::endl;
        
        // Copy data to GPU
        char* d_sequence_data = memory_manager.get_sequence_buffer();
        GPUGenomeInfo* d_genome_info = memory_manager.get_genome_info_buffer();
        GPUMinimizerHit* d_minimizer_hits = memory_manager.get_minimizer_buffer();
        uint32_t* d_hit_counter = memory_manager.get_global_counter();
        
        cudaMemcpy(d_sequence_data, concatenated_sequences.data(), 
                   concatenated_sequences.length(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_genome_info, genome_infos.data(), 
                   genome_infos.size() * sizeof(GPUGenomeInfo), cudaMemcpyHostToDevice);
        
        // Reset hit counter
        uint32_t zero = 0;
        cudaMemcpy(d_hit_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Launch minimizer extraction kernel
        MinimizerParams params;
        params.k = config.k_value;
        params.ell = config.ell_value;
        params.spaces = config.spaces_value;
        params.xor_mask = config.xor_mask;
        
        GPUBatchData batch_data;
        batch_data.d_sequence_data = d_sequence_data;
        batch_data.d_genome_info = d_genome_info;
        batch_data.d_minimizer_hits = d_minimizer_hits;
        batch_data.d_global_counter = d_hit_counter;
        batch_data.max_genomes = genome_infos.size();
        batch_data.max_minimizers = config.max_minimizers_per_batch;
        batch_data.sequence_buffer_size = concatenated_sequences.length();
        
        auto kernel_start = std::chrono::high_resolution_clock::now();
        
        uint32_t hits_extracted = 0;
        bool kernel_success = launch_improved_minimizer_kernel(
            batch_data, params, 0, 0, &hits_extracted
        );
        
        auto kernel_end = std::chrono::high_resolution_clock::now();
        auto kernel_duration = std::chrono::duration<double>(kernel_end - kernel_start).count();
        
        if (!kernel_success) {
            std::cerr << "Minimizer extraction kernel failed!" << std::endl;
            return false;
        }
        
        std::cout << "Minimizers extracted: " << hits_extracted 
                  << " (in " << std::fixed << std::setprecision(3) 
                  << kernel_duration << " seconds)" << std::endl;
        
        // Copy minimizers back to host for analysis
        if (hits_extracted > 0) {
            std::vector<GPUMinimizerHit> batch_minimizers(hits_extracted);
            cudaMemcpy(batch_minimizers.data(), d_minimizer_hits,
                       hits_extracted * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
            
            // Collect for final analysis
            all_minimizers.insert(all_minimizers.end(), 
                                 batch_minimizers.begin(), 
                                 batch_minimizers.end());
            
            total_minimizers_extracted += hits_extracted;
            
            // Print sample of minimizers from this batch
            std::cout << "\nSample minimizers from batch:" << std::endl;
            for (size_t i = 0; i < std::min(size_t(3), batch_minimizers.size()); i++) {
                const auto& hit = batch_minimizers[i];
                std::cout << "  Hash: " << std::hex << hit.minimizer_hash << std::dec
                          << ", Genome: " << hit.genome_id
                          << ", Position: " << hit.position
                          << ", Taxon: " << hit.taxon_id << std::endl;
            }
        }
        
        // Clean up temporary files
        for (const auto& file : batch_files) {
            std::filesystem::remove(file);
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration<double>(total_end - total_start).count();
    
    // Print final statistics
    std::cout << "\n=== FINAL STATISTICS ===" << std::endl;
    std::cout << "Total batches processed: " << batch_number << std::endl;
    std::cout << "Total genomes processed: " << streaming_processor.get_total_genomes() << std::endl;
    std::cout << "Total bases processed: " << streaming_processor.get_total_bases() << std::endl;
    std::cout << "Total minimizers extracted: " << total_minimizers_extracted << std::endl;
    std::cout << "Total processing time: " << std::fixed << std::setprecision(2) 
              << total_duration << " seconds" << std::endl;
    
    if (streaming_processor.get_total_bases() > 0) {
        double rate = streaming_processor.get_total_bases() / total_duration / 1024 / 1024;
        std::cout << "Processing rate: " << std::fixed << std::setprecision(2) 
                  << rate << " MB/s" << std::endl;
    }
    
    // Analyze all collected minimizers
    print_minimizer_stats(all_minimizers);
    
    // Clean up
    memory_manager.free_all_allocations();
    
    return true;
}

// Simple test kernel to verify GPU functionality
__global__ void test_minimizer_kernel(const char* sequence, int length, uint64_t* result) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        // Extract one test minimizer
        if (length >= 31) {
            *result = extract_minimizer_sliding_window(sequence, 0, 31, 31, 7, 0x3c8bfbb395c60474ULL, length);
        } else {
            *result = 0;
        }
    }
}

// Test basic minimizer extraction functionality
bool test_basic_minimizer_extraction() {
    std::cout << "\n=== Testing Basic Minimizer Extraction ===" << std::endl;
    
    // Test sequence
    std::string test_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    
    // Allocate GPU memory
    char* d_sequence;
    uint64_t* d_result;
    cudaMalloc(&d_sequence, test_seq.length() + 1);
    cudaMalloc(&d_result, sizeof(uint64_t));
    
    // Copy sequence to GPU
    cudaMemcpy(d_sequence, test_seq.c_str(), test_seq.length() + 1, cudaMemcpyHostToDevice);
    
    // Launch test kernel
    test_minimizer_kernel<<<1, 1>>>(d_sequence, test_seq.length(), d_result);
    cudaDeviceSynchronize();
    
    // Get result
    uint64_t result;
    cudaMemcpy(&result, d_result, sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
    std::cout << "Test sequence: " << test_seq.substr(0, 50) << "..." << std::endl;
    std::cout << "Minimizer hash: " << std::hex << result << std::dec << std::endl;
    
    // Clean up
    cudaFree(d_sequence);
    cudaFree(d_result);
    
    return result != 0 && result != UINT64_MAX;
}

int main(int argc, char** argv) {
    // Check CUDA availability
    int device_count;
    cudaError_t cuda_status = cudaGetDeviceCount(&device_count);
    if (cuda_status != cudaSuccess || device_count == 0) {
        std::cerr << "No CUDA devices available!" << std::endl;
        return 1;
    }
    
    // Print GPU info
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "Total memory: " << (prop.totalGlobalMem / 1024 / 1024) << " MB" << std::endl;
    
    // First test basic functionality
    if (!test_basic_minimizer_extraction()) {
        std::cerr << "Basic minimizer extraction test failed!" << std::endl;
        return 1;
    }
    std::cout << "✓ Basic minimizer extraction test passed" << std::endl;
    
    // Configure test
    TestConfig config;
    
    // Override with command line arguments if provided
    if (argc > 1) {
        config.input_fna_path = argv[1];
    }
    if (argc > 2) {
        config.batch_size = std::stoi(argv[2]);
    }
    if (argc > 3) {
        config.k_value = std::stoi(argv[3]);
    }
    
    // Create temp directory
    std::filesystem::create_directories(config.temp_dir);
    
    // Run the main test
    bool success = test_streaming_minimizer_extraction(config);
    
    // Clean up temp directory
    try {
        std::filesystem::remove_all(config.temp_dir);
    } catch (...) {
        // Ignore cleanup errors
    }
    
    if (success) {
        std::cout << "\n✓ All tests passed!" << std::endl;
        return 0;
    } else {
        std::cerr << "\n✗ Tests failed!" << std::endl;
        return 1;
    }
}