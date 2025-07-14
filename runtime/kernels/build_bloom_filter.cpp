// build_bloom_filter.cpp
// Builds a 15-mer nucleotide bloom filter from FASTA files
// Supports both single FASTA files and directories of sequences
// Includes forward and reverse complement k-mers

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <cuda_runtime.h>

// K-mer configuration
const int KMER_LENGTH = 15;
const uint64_t KMER_MASK = (1ULL << (2 * KMER_LENGTH)) - 1;

// Base encoding
const uint8_t BASE_A = 0;
const uint8_t BASE_C = 1;
const uint8_t BASE_G = 2;
const uint8_t BASE_T = 3;
const uint8_t BASE_INVALID = 255;

// External bloom filter functions
extern "C" {
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int save_bloom_filter(void* filter, const char* filename);
}

// Encode a nucleotide base
inline uint8_t encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return BASE_A;
        case 'C': case 'c': return BASE_C;
        case 'G': case 'g': return BASE_G;
        case 'T': case 't': case 'U': case 'u': return BASE_T;
        default: return BASE_INVALID;
    }
}

// Encode a k-mer to 64-bit integer
uint64_t encode_kmer(const std::string& sequence, size_t start) {
    uint64_t encoded = 0;
    
    for (size_t i = 0; i < KMER_LENGTH; i++) {
        uint8_t base = encode_base(sequence[start + i]);
        if (base == BASE_INVALID) {
            return UINT64_MAX; // Invalid k-mer
        }
        encoded = (encoded << 2) | base;
    }
    
    return encoded;
}

// Compute reverse complement of encoded k-mer
uint64_t reverse_complement_kmer(uint64_t kmer) {
    uint64_t rc = 0;
    
    for (int i = 0; i < KMER_LENGTH; i++) {
        uint8_t base = kmer & 3;
        uint8_t complement;
        
        switch (base) {
            case BASE_A: complement = BASE_T; break;
            case BASE_T: complement = BASE_A; break;
            case BASE_C: complement = BASE_G; break;
            case BASE_G: complement = BASE_C; break;
            default: return UINT64_MAX;
        }
        
        rc = (rc << 2) | complement;
        kmer >>= 2;
    }
    
    return rc;
}

// Extract k-mers from a sequence
void extract_kmers(const std::string& sequence, std::unordered_set<uint64_t>& kmers) {
    if (sequence.length() < KMER_LENGTH) {
        return;
    }
    
    for (size_t i = 0; i <= sequence.length() - KMER_LENGTH; i++) {
        uint64_t kmer = encode_kmer(sequence, i);
        if (kmer != UINT64_MAX) {
            kmers.insert(kmer);
            
            // Also add reverse complement
            uint64_t rc = reverse_complement_kmer(kmer);
            if (rc != UINT64_MAX) {
                kmers.insert(rc);
            }
        }
    }
}

// Read sequences from FASTA file
void read_fasta(const std::string& filename, std::unordered_set<uint64_t>& kmers) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    std::string line, sequence;
    int sequences_processed = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Process previous sequence if exists
            if (!sequence.empty()) {
                extract_kmers(sequence, kmers);
                sequences_processed++;
                sequence.clear();
            }
        } else {
            // Append to current sequence
            sequence += line;
        }
    }
    
    // Process last sequence
    if (!sequence.empty()) {
        extract_kmers(sequence, kmers);
        sequences_processed++;
    }
    
    file.close();
    std::cout << "Processed " << sequences_processed << " sequences from " << filename << std::endl;
}

// Read sequences from JSON files in FQ genes directory
void read_fq_genes_directory(const std::string& dir_path, std::unordered_set<uint64_t>& kmers) {
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << dir_path << std::endl;
        return;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        
        // Skip . and ..
        if (filename == "." || filename == "..") continue;
        
        std::string full_path = dir_path + "/" + filename;
        struct stat st;
        if (stat(full_path.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            // Process subdirectory
            read_fq_genes_directory(full_path, kmers);
        } else if (filename.find(".json") != std::string::npos) {
            // For JSON files, we need to parse them differently
            // For now, skip JSON parsing - in production, use a JSON library
            std::cout << "Note: JSON parsing not implemented. Skipping " << full_path << std::endl;
        } else if (filename.find(".fa") != std::string::npos || 
                   filename.find(".fasta") != std::string::npos ||
                   filename.find(".fna") != std::string::npos) {
            // Process FASTA files
            read_fasta(full_path, kmers);
        }
    }
    
    closedir(dir);
}

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <input> <output_bloom_filter.bin>" << std::endl;
    std::cerr << "  input: Can be either:" << std::endl;
    std::cerr << "    - A FASTA file (e.g., AMR_CDS.fa)" << std::endl;
    std::cerr << "    - A directory containing sequence files" << std::endl;
    std::cerr << "  output: Path to save the bloom filter binary" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Examples:" << std::endl;
    std::cerr << "  " << program_name << " ../../data/AMR_CDS.fa amr_bloom_filter.bin" << std::endl;
    std::cerr << "  " << program_name << " ../../data/fq_genes fq_bloom_filter.bin" << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string input_path = argv[1];
    std::string output_file = argv[2];
    
    // Check if input is file or directory
    struct stat st;
    if (stat(input_path.c_str(), &st) != 0) {
        std::cerr << "Error: Cannot access " << input_path << std::endl;
        return 1;
    }
    
    // Extract k-mers
    std::unordered_set<uint64_t> unique_kmers;
    
    if (S_ISDIR(st.st_mode)) {
        std::cout << "Processing directory: " << input_path << std::endl;
        read_fq_genes_directory(input_path, unique_kmers);
    } else {
        std::cout << "Processing file: " << input_path << std::endl;
        read_fasta(input_path, unique_kmers);
    }
    
    std::cout << "\nExtracted " << unique_kmers.size() << " unique " << KMER_LENGTH 
              << "-mers (including reverse complements)" << std::endl;
    
    if (unique_kmers.empty()) {
        std::cerr << "Error: No k-mers extracted" << std::endl;
        return 1;
    }
    
    // Convert to vector for GPU transfer
    std::vector<uint64_t> kmer_vector(unique_kmers.begin(), unique_kmers.end());
    
    // Allocate GPU memory for k-mers
    uint64_t* d_kmers;
    size_t kmers_size = kmer_vector.size() * sizeof(uint64_t);
    cudaError_t err = cudaMalloc(&d_kmers, kmers_size);
    if (err != cudaSuccess) {
        std::cerr << "Error: Failed to allocate GPU memory: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    // Copy k-mers to GPU
    err = cudaMemcpy(d_kmers, kmer_vector.data(), kmers_size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "Error: Failed to copy k-mers to GPU: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_kmers);
        return 1;
    }
    
    // Create bloom filter
    void* bloom_filter = create_bloom_filter(KMER_LENGTH);
    if (!bloom_filter) {
        std::cerr << "Error: Failed to create bloom filter" << std::endl;
        cudaFree(d_kmers);
        return 1;
    }
    
    // Build bloom filter from k-mers
    std::cout << "\nBuilding bloom filter..." << std::endl;
    if (build_bloom_filter_from_index(bloom_filter, d_kmers, kmer_vector.size()) != 0) {
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
    
    std::cout << "Successfully saved bloom filter to: " << output_file << std::endl;
    
    // Cleanup
    destroy_bloom_filter(bloom_filter);
    cudaFree(d_kmers);
    
    return 0;
}