// build_minimizer_bloom_filter.cpp
// Builds a minimizer-based bloom filter from nucleotide sequences
// Can combine multiple FASTA files (FQ resistance genes + AMR genes)
// Uses minimizers instead of all k-mers for better efficiency

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <deque>
#include <cuda_runtime.h>

// Minimizer configuration
const int KMER_LENGTH = 15;        // k-mer size
const int WINDOW_SIZE = 10;        // window size for minimizer selection
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

// Extract minimizers from a sequence using sliding window
void extract_minimizers(const std::string& sequence, std::unordered_set<uint64_t>& minimizers) {
    if (sequence.length() < KMER_LENGTH + WINDOW_SIZE - 1) {
        return;
    }
    
    // Use a deque to track k-mers in current window
    std::deque<std::pair<uint64_t, size_t>> window;
    
    // Process each window position
    for (size_t i = 0; i <= sequence.length() - KMER_LENGTH; i++) {
        uint64_t kmer = encode_kmer(sequence, i);
        if (kmer == UINT64_MAX) continue;
        
        // Remove k-mers that are outside the window
        while (!window.empty() && window.front().second <= i - WINDOW_SIZE) {
            window.pop_front();
        }
        
        // Remove k-mers that are larger than current k-mer
        while (!window.empty() && window.back().first > kmer) {
            window.pop_back();
        }
        
        // Add current k-mer
        window.push_back({kmer, i});
        
        // The front of the deque is the minimizer for this window
        if (i >= WINDOW_SIZE - 1) {
            uint64_t minimizer = window.front().first;
            minimizers.insert(minimizer);
            
            // Also add reverse complement
            uint64_t rc = reverse_complement_kmer(minimizer);
            if (rc != UINT64_MAX) {
                minimizers.insert(rc);
            }
        }
    }
}

// Read sequences from FASTA file and extract minimizers
void process_fasta(const std::string& filename, std::unordered_set<uint64_t>& minimizers) {
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
                extract_minimizers(sequence, minimizers);
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
        extract_minimizers(sequence, minimizers);
        sequences_processed++;
    }
    
    file.close();
    std::cout << "Processed " << sequences_processed << " sequences from " << filename << std::endl;
}

// Process directory of FASTA files
void process_directory(const std::string& dir_path, std::unordered_set<uint64_t>& minimizers) {
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << dir_path << std::endl;
        return;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        
        if (filename == "." || filename == "..") continue;
        
        std::string full_path = dir_path + "/" + filename;
        struct stat st;
        if (stat(full_path.c_str(), &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                process_directory(full_path, minimizers);
            } else if (filename.find(".fa") != std::string::npos || 
                       filename.find(".fasta") != std::string::npos ||
                       filename.find(".fna") != std::string::npos) {
                process_fasta(full_path, minimizers);
            }
        }
    }
    
    closedir(dir);
}

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <output_bloom_filter.bin> [input1] [input2] ..." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Arguments:" << std::endl;
    std::cerr << "  output:      Path to save the bloom filter binary" << std::endl;
    std::cerr << "  inputs:      Optional FASTA files or directories (see defaults below)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Default behavior (no inputs specified):" << std::endl;
    std::cerr << "  Builds from both:" << std::endl;
    std::cerr << "    - ../../data/resistance_db/reference_sequences.fasta (FQ resistance genes)" << std::endl;
    std::cerr << "    - ../../data/AMR_CDS.fa (AMR genes)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Examples:" << std::endl;
    std::cerr << "  # Use default sources (recommended):" << std::endl;
    std::cerr << "  " << program_name << " combined_bloom_filter.bin" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  # Custom sources:" << std::endl;
    std::cerr << "  " << program_name << " custom_bloom.bin /path/to/sequences.fa" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Configuration:" << std::endl;
    std::cerr << "  K-mer length: " << KMER_LENGTH << std::endl;
    std::cerr << "  Window size: " << WINDOW_SIZE << std::endl;
    std::cerr << "  Includes reverse complements: Yes" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string output_file = argv[1];
    
    // Extract minimizers from all input files/directories
    std::unordered_set<uint64_t> all_minimizers;
    
    // If no inputs provided, use default sources
    if (argc == 2) {
        std::cout << "No input files specified. Using default sources..." << std::endl;
        
        // Default FQ resistance sequences
        std::string fq_resistance_file = "../../data/resistance_db/reference_sequences.fasta";
        struct stat st;
        if (stat(fq_resistance_file.c_str(), &st) == 0) {
            std::cout << "\nProcessing FQ resistance genes: " << fq_resistance_file << std::endl;
            process_fasta(fq_resistance_file, all_minimizers);
        } else {
            std::cerr << "Warning: FQ resistance file not found: " << fq_resistance_file << std::endl;
        }
        
        // Default AMR genes
        std::string amr_file = "../../data/AMR_CDS.fa";
        if (stat(amr_file.c_str(), &st) == 0) {
            std::cout << "\nProcessing AMR genes: " << amr_file << std::endl;
            process_fasta(amr_file, all_minimizers);
        } else {
            std::cerr << "Warning: AMR genes file not found: " << amr_file << std::endl;
        }
        
        if (all_minimizers.empty()) {
            std::cerr << "Error: No default files found. Please specify input files." << std::endl;
            return 1;
        }
    } else {
        // Use provided input files
        for (int i = 2; i < argc; i++) {
            std::string input_path = argv[i];
            
            struct stat st;
            if (stat(input_path.c_str(), &st) != 0) {
                std::cerr << "Error: Cannot access " << input_path << std::endl;
                continue;
            }
            
            std::cout << "\nProcessing: " << input_path << std::endl;
            
            if (S_ISDIR(st.st_mode)) {
                process_directory(input_path, all_minimizers);
            } else {
                process_fasta(input_path, all_minimizers);
            }
        }
    }
    
    std::cout << "\nExtracted " << all_minimizers.size() << " unique minimizers "
              << "(k=" << KMER_LENGTH << ", w=" << WINDOW_SIZE << ", including RC)" << std::endl;
    
    if (all_minimizers.empty()) {
        std::cerr << "Error: No minimizers extracted" << std::endl;
        return 1;
    }
    
    // Convert to vector for GPU transfer
    std::vector<uint64_t> minimizer_vector(all_minimizers.begin(), all_minimizers.end());
    
    // Allocate GPU memory for minimizers
    uint64_t* d_minimizers;
    size_t minimizers_size = minimizer_vector.size() * sizeof(uint64_t);
    cudaError_t err = cudaMalloc(&d_minimizers, minimizers_size);
    if (err != cudaSuccess) {
        std::cerr << "Error: Failed to allocate GPU memory: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    // Copy minimizers to GPU
    err = cudaMemcpy(d_minimizers, minimizer_vector.data(), minimizers_size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "Error: Failed to copy minimizers to GPU: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_minimizers);
        return 1;
    }
    
    // Create bloom filter
    void* bloom_filter = create_bloom_filter(KMER_LENGTH);
    if (!bloom_filter) {
        std::cerr << "Error: Failed to create bloom filter" << std::endl;
        cudaFree(d_minimizers);
        return 1;
    }
    
    // Build bloom filter from minimizers
    std::cout << "\nBuilding bloom filter..." << std::endl;
    if (build_bloom_filter_from_index(bloom_filter, d_minimizers, minimizer_vector.size()) != 0) {
        std::cerr << "Error: Failed to build bloom filter" << std::endl;
        destroy_bloom_filter(bloom_filter);
        cudaFree(d_minimizers);
        return 1;
    }
    
    // Save bloom filter
    if (save_bloom_filter(bloom_filter, output_file.c_str()) != 0) {
        std::cerr << "Error: Failed to save bloom filter to " << output_file << std::endl;
        destroy_bloom_filter(bloom_filter);
        cudaFree(d_minimizers);
        return 1;
    }
    
    std::cout << "Successfully saved bloom filter to: " << output_file << std::endl;
    
    // Cleanup
    destroy_bloom_filter(bloom_filter);
    cudaFree(d_minimizers);
    
    return 0;
}