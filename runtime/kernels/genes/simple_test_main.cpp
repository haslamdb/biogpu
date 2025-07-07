// simple_test_main.cpp
// Test program for simplified translated search

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cuda_runtime.h>
#include "simplified_translated_search.h"

// Simple FASTQ reader
std::vector<std::pair<std::string, std::string>> readFastq(const std::string& filename, int max_reads = 1000) {
    std::vector<std::pair<std::string, std::string>> reads;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return reads;
    }
    
    std::string header, sequence, plus, quality;
    int count = 0;
    
    while (std::getline(file, header) && count < max_reads) {
        if (std::getline(file, sequence) && 
            std::getline(file, plus) && 
            std::getline(file, quality)) {
            
            if (!sequence.empty()) {
                std::string id = header.substr(1); // Remove '@'
                reads.push_back({id, sequence});
                count++;
            }
        }
    }
    
    file.close();
    return reads;
}

// Copy reads to GPU
bool copyReadsToGPU(const std::vector<std::pair<std::string, std::string>>& reads,
                   char*& d_reads, int*& d_read_lengths, int*& d_read_offsets) {
    
    // Concatenate reads
    std::vector<char> all_reads;
    std::vector<int> lengths;
    std::vector<int> offsets;
    
    for (const auto& read_pair : reads) {
        offsets.push_back(all_reads.size());
        lengths.push_back(read_pair.second.length());
        
        for (char c : read_pair.second) {
            all_reads.push_back(c);
        }
    }
    
    // Allocate GPU memory
    cudaError_t err;
    
    err = cudaMalloc(&d_reads, all_reads.size());
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for reads: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    err = cudaMalloc(&d_read_lengths, reads.size() * sizeof(int));
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for lengths: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    err = cudaMalloc(&d_read_offsets, reads.size() * sizeof(int));
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for offsets: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    // Copy data
    err = cudaMemcpy(d_reads, all_reads.data(), all_reads.size(), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "Failed to copy reads to GPU: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    err = cudaMemcpy(d_read_lengths, lengths.data(), lengths.size() * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "Failed to copy lengths to GPU: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    err = cudaMemcpy(d_read_offsets, offsets.data(), offsets.size() * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "Failed to copy offsets to GPU: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <protein_fasta> <reads_fastq> [output_prefix]" << std::endl;
        return 1;
    }
    
    std::string protein_fasta = argv[1];
    std::string reads_fastq = argv[2];
    std::string output_prefix = (argc > 3) ? argv[3] : "simple_results";
    
    std::cout << "=== Simple Translated Search Test ===" << std::endl;
    std::cout << "Protein database: " << protein_fasta << std::endl;
    std::cout << "Reads file: " << reads_fastq << std::endl;
    std::cout << "Output prefix: " << output_prefix << std::endl;
    
    // Create search engine
    void* engine = create_simple_translated_search(10000);
    if (!engine) {
        std::cerr << "Failed to create search engine" << std::endl;
        return 1;
    }
    
    // Load proteins
    std::cout << "\nLoading proteins..." << std::endl;
    if (load_simple_proteins(engine, protein_fasta.c_str()) != 0) {
        std::cerr << "Failed to load proteins" << std::endl;
        destroy_simple_translated_search(engine);
        return 1;
    }
    
    // Read FASTQ
    std::cout << "Reading FASTQ..." << std::endl;
    auto reads = readFastq(reads_fastq, 5000); // Limit for testing
    if (reads.empty()) {
        std::cerr << "No reads loaded" << std::endl;
        destroy_simple_translated_search(engine);
        return 1;
    }
    
    std::cout << "Loaded " << reads.size() << " reads" << std::endl;
    
    // Copy to GPU
    std::cout << "Copying reads to GPU..." << std::endl;
    char* d_reads = nullptr;
    int* d_read_lengths = nullptr;
    int* d_read_offsets = nullptr;
    
    if (!copyReadsToGPU(reads, d_reads, d_read_lengths, d_read_offsets)) {
        std::cerr << "Failed to copy reads to GPU" << std::endl;
        destroy_simple_translated_search(engine);
        return 1;
    }
    
    // Search
    std::cout << "Searching..." << std::endl;
    
    // Note: This is a simplified interface
    // In practice, you'd need to handle the results properly
    int result = search_simple_reads(
        engine,
        d_reads,
        d_read_lengths,
        d_read_offsets,
        reads.size(),
        nullptr, // alignments_out - would need proper handling
        nullptr  // coverage_out - would need proper handling
    );
    
    if (result == 0) {
        std::cout << "Search completed successfully!" << std::endl;
        // Results would be written here
    } else {
        std::cerr << "Search failed" << std::endl;
    }
    
    // Cleanup
    if (d_reads) cudaFree(d_reads);
    if (d_read_lengths) cudaFree(d_read_lengths);
    if (d_read_offsets) cudaFree(d_read_offsets);
    
    destroy_simple_translated_search(engine);
    
    return (result == 0) ? 0 : 1;
}
