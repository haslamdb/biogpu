// build_bloom_filter_simple.cpp
// Simplified version that processes existing kmer_index.bin files
// This allows reuse of already extracted k-mers from the pipelines

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cuda_runtime.h>

// External bloom filter functions
extern "C" {
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int save_bloom_filter(void* filter, const char* filename);
}

struct KmerIndexHeader {
    uint32_t kmer_length;
    uint32_t num_kmers;
    uint32_t num_sequences;
    uint32_t has_reverse_complement;
};

bool load_kmer_index(const std::string& index_path, std::vector<uint64_t>& kmers, int& kmer_length) {
    std::ifstream file(index_path, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open k-mer index: " << index_path << std::endl;
        return false;
    }
    
    // Read header
    KmerIndexHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    
    kmer_length = header.kmer_length;
    std::cout << "K-mer index info:" << std::endl;
    std::cout << "  K-mer length: " << header.kmer_length << std::endl;
    std::cout << "  Number of k-mers: " << header.num_kmers << std::endl;
    std::cout << "  Number of sequences: " << header.num_sequences << std::endl;
    std::cout << "  Has reverse complement: " << (header.has_reverse_complement ? "Yes" : "No") << std::endl;
    
    // Read k-mers
    kmers.resize(header.num_kmers);
    file.read(reinterpret_cast<char*>(kmers.data()), header.num_kmers * sizeof(uint64_t));
    
    file.close();
    return true;
}

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <kmer_index.bin> <output_bloom_filter.bin>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Examples:" << std::endl;
    std::cerr << "  # For FQ resistance genes:" << std::endl;
    std::cerr << "  " << program_name << " data/integrated_clean_db/nucleotide/kmer_index.bin fq_bloom_filter.bin" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  # For AMR genes (if k-mer index exists):" << std::endl;
    std::cerr << "  " << program_name << " data/amr_nucleotide_db/kmer_index.bin amr_bloom_filter.bin" << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string input_index = argv[1];
    std::string output_file = argv[2];
    
    // Load k-mer index
    std::vector<uint64_t> kmers;
    int kmer_length;
    
    if (!load_kmer_index(input_index, kmers, kmer_length)) {
        return 1;
    }
    
    if (kmers.empty()) {
        std::cerr << "Error: No k-mers loaded from index" << std::endl;
        return 1;
    }
    
    // Allocate GPU memory for k-mers
    uint64_t* d_kmers;
    size_t kmers_size = kmers.size() * sizeof(uint64_t);
    cudaError_t err = cudaMalloc(&d_kmers, kmers_size);
    if (err != cudaSuccess) {
        std::cerr << "Error: Failed to allocate GPU memory: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    // Copy k-mers to GPU
    err = cudaMemcpy(d_kmers, kmers.data(), kmers_size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "Error: Failed to copy k-mers to GPU: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_kmers);
        return 1;
    }
    
    // Create bloom filter
    void* bloom_filter = create_bloom_filter(kmer_length);
    if (!bloom_filter) {
        std::cerr << "Error: Failed to create bloom filter" << std::endl;
        cudaFree(d_kmers);
        return 1;
    }
    
    // Build bloom filter from k-mers
    std::cout << "\nBuilding bloom filter..." << std::endl;
    if (build_bloom_filter_from_index(bloom_filter, d_kmers, kmers.size()) != 0) {
        std::cerr << "Error: Failed to build bloom filter" << std::endl;
        destroy_bloom_filter(bloom_filter);
        cudaFree(d_kmers);
        return 1;
    }
    
    // Save bloom filter
    if (save_bloom_filter(bloom_filter, output_file.c_str()) != 0) {
        std::cerr << "Error: Failed to save bloom filter to " << output_file << std::endl;
        destroy_bloom_filter(bloom_filter);
        cudaFree(d_kmers);
        return 1;
    }
    
    std::cout << "\nSuccessfully saved bloom filter to: " << output_file << std::endl;
    
    // Cleanup
    destroy_bloom_filter(bloom_filter);
    cudaFree(d_kmers);
    
    return 0;
}