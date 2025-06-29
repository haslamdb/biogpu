// test_minimal_streaming.cu
// Minimal test for streaming FNA processor with GPU minimizer extraction
// Only uses implemented components

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cuda_runtime.h>
#include <filesystem>

// Include only what we need
#include "processing/genome_file_processor.h"
#include "gpu/gpu_minimizer_extraction.cuh"
#include "gpu_kraken_types.h"

// Simple test configuration
struct TestConfig {
    std::string input_fna_path;
    std::string temp_dir = "/tmp/kraken_test";
    size_t batch_size = 5;  // Small batch for testing
    
    // Minimizer parameters
    uint32_t k_value = 31;
    uint32_t ell_value = 31;
    uint32_t spaces_value = 7;
    uint64_t xor_mask = 0x3c8bfbb395c60474ULL;
};

// Simple kernel to extract minimizers without all the extra features
__global__ void simple_minimizer_kernel(
    const char* sequences,
    const uint32_t* seq_offsets,
    const uint32_t* seq_lengths,
    const uint32_t* taxon_ids,
    int num_sequences,
    uint64_t* minimizer_buffer,
    uint32_t* position_buffer,
    uint32_t* taxon_buffer,
    uint32_t* global_counter,
    uint32_t max_minimizers,
    uint32_t k,
    uint32_t ell,
    uint32_t spaces,
    uint64_t xor_mask) {
    
    int seq_idx = blockIdx.x;
    if (seq_idx >= num_sequences) return;
    
    uint32_t seq_start = seq_offsets[seq_idx];
    uint32_t seq_length = seq_lengths[seq_idx];
    uint32_t taxon_id = taxon_ids[seq_idx];
    
    if (seq_length < k) return;
    
    const char* sequence = sequences + seq_start;
    uint32_t num_kmers = seq_length - k + 1;
    
    // Each thread processes a portion of k-mers
    uint32_t tid = threadIdx.x;
    uint32_t threads_per_block = blockDim.x;
    
    for (uint32_t pos = tid; pos < num_kmers; pos += threads_per_block) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, pos, k, ell, spaces, xor_mask, seq_length
        );
        
        if (minimizer != UINT64_MAX) {
            uint32_t idx = atomicAdd(global_counter, 1);
            if (idx < max_minimizers) {
                minimizer_buffer[idx] = minimizer;
                position_buffer[idx] = pos;
                taxon_buffer[idx] = taxon_id;
            }
        }
    }
}

bool test_minimal_streaming(const TestConfig& config) {
    std::cout << "=== Minimal Streaming FNA Test ===" << std::endl;
    std::cout << "Input: " << config.input_fna_path << std::endl;
    std::cout << "Batch size: " << config.batch_size << std::endl;
    
    // Initialize streaming processor
    StreamingFnaProcessor processor(
        config.input_fna_path, 
        config.temp_dir, 
        config.batch_size
    );
    
    // Process batches
    int batch_num = 0;
    size_t total_minimizers = 0;
    
    std::vector<std::string> batch_files;
    std::vector<uint32_t> batch_taxons;
    
    while (processor.process_next_batch(batch_files, batch_taxons)) {
        batch_num++;
        std::cout << "\n--- Batch " << batch_num << " ---" << std::endl;
        std::cout << "Genomes: " << batch_files.size() << std::endl;
        
        // Load sequences from temp files
        std::vector<std::string> sequences;
        std::vector<uint32_t> seq_lengths;
        std::vector<uint32_t> seq_offsets;
        std::vector<uint32_t> taxon_ids;
        std::string concatenated;
        
        GenomeFileProcessor file_processor;
        
        for (size_t i = 0; i < batch_files.size(); i++) {
            auto seqs = file_processor.load_sequences_from_fasta(batch_files[i]);
            for (const auto& seq : seqs) {
                seq_offsets.push_back(concatenated.length());
                seq_lengths.push_back(seq.length());
                taxon_ids.push_back(batch_taxons[i]);
                concatenated += seq;
                
                std::cout << "  Genome " << i << ": " << seq.length() 
                          << " bp (taxon=" << batch_taxons[i] << ")" << std::endl;
            }
        }
        
        if (concatenated.empty()) continue;
        
        // Allocate GPU memory
        char* d_sequences;
        uint32_t* d_offsets;
        uint32_t* d_lengths;
        uint32_t* d_taxons;
        uint64_t* d_minimizers;
        uint32_t* d_positions;
        uint32_t* d_taxon_out;
        uint32_t* d_counter;
        
        size_t max_minimizers = concatenated.length();  // Worst case
        
        cudaMalloc(&d_sequences, concatenated.length());
        cudaMalloc(&d_offsets, seq_offsets.size() * sizeof(uint32_t));
        cudaMalloc(&d_lengths, seq_lengths.size() * sizeof(uint32_t));
        cudaMalloc(&d_taxons, taxon_ids.size() * sizeof(uint32_t));
        cudaMalloc(&d_minimizers, max_minimizers * sizeof(uint64_t));
        cudaMalloc(&d_positions, max_minimizers * sizeof(uint32_t));
        cudaMalloc(&d_taxon_out, max_minimizers * sizeof(uint32_t));
        cudaMalloc(&d_counter, sizeof(uint32_t));
        
        // Copy data to GPU
        cudaMemcpy(d_sequences, concatenated.data(), concatenated.length(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_offsets, seq_offsets.data(), seq_offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_lengths, seq_lengths.data(), seq_lengths.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_taxons, taxon_ids.data(), taxon_ids.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemset(d_counter, 0, sizeof(uint32_t));
        
        // Launch kernel
        int threads = 256;
        int blocks = seq_offsets.size();
        
        simple_minimizer_kernel<<<blocks, threads>>>(
            d_sequences, d_offsets, d_lengths, d_taxons,
            seq_offsets.size(),
            d_minimizers, d_positions, d_taxon_out, d_counter,
            max_minimizers,
            config.k_value, config.ell_value, config.spaces_value, config.xor_mask
        );
        
        cudaDeviceSynchronize();
        
        auto err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        // Get results
        uint32_t num_minimizers;
        cudaMemcpy(&num_minimizers, d_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::cout << "Minimizers extracted: " << num_minimizers << std::endl;
        total_minimizers += num_minimizers;
        
        // Sample first few minimizers
        if (num_minimizers > 0) {
            std::vector<uint64_t> sample_mins(std::min(5u, num_minimizers));
            std::vector<uint32_t> sample_pos(std::min(5u, num_minimizers));
            std::vector<uint32_t> sample_tax(std::min(5u, num_minimizers));
            
            cudaMemcpy(sample_mins.data(), d_minimizers, sample_mins.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(sample_pos.data(), d_positions, sample_pos.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(sample_tax.data(), d_taxon_out, sample_tax.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            
            std::cout << "Sample minimizers:" << std::endl;
            for (size_t i = 0; i < sample_mins.size(); i++) {
                std::cout << "  [" << i << "] hash=" << std::hex << sample_mins[i] << std::dec
                          << ", pos=" << sample_pos[i]
                          << ", taxon=" << sample_tax[i] << std::endl;
            }
        }
        
        // Cleanup GPU memory
        cudaFree(d_sequences);
        cudaFree(d_offsets);
        cudaFree(d_lengths);
        cudaFree(d_taxons);
        cudaFree(d_minimizers);
        cudaFree(d_positions);
        cudaFree(d_taxon_out);
        cudaFree(d_counter);
        
        // Clean up temp files
        for (const auto& file : batch_files) {
            std::filesystem::remove(file);
        }
    }
    
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total batches: " << batch_num << std::endl;
    std::cout << "Total genomes: " << processor.get_total_genomes() << std::endl;
    std::cout << "Total minimizers: " << total_minimizers << std::endl;
    
    return true;
}

int main(int argc, char** argv) {
    // Check CUDA
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        std::cerr << "No CUDA devices!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "GPU: " << prop.name << std::endl;
    
    // Setup test
    TestConfig config;
    config.input_fna_path = (argc > 1) ? argv[1] : "/home/david/Documents/Code/biogpu/data/ref_genomes.fna";
    
    // Create temp dir
    std::filesystem::create_directories(config.temp_dir);
    
    // Run test
    bool success = test_minimal_streaming(config);
    
    // Cleanup
    std::filesystem::remove_all(config.temp_dir);
    
    return success ? 0 : 1;
}