// test_minimizer.cpp
// Simple test program for minimizer extraction
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include "minimizer_common.h"

// Forward declaration of MinimizerExtractor class
class MinimizerExtractor {
public:
    MinimizerExtractor(int k_size = 31, int window_size = 15);
    ~MinimizerExtractor();
    std::vector<std::vector<Minimizer>> extract_minimizers(const std::vector<std::string>& sequences);
private:
    void allocate_device_memory(size_t num_reads, size_t total_sequence_length);
    int k;
    int m;
    size_t allocated_reads;
    void* d_sequences;
    void* d_sequence_offsets;
    void* d_sequence_lengths;
    void* d_minimizers;
    void* d_minimizer_counts;
};

// Simple FASTQ reader for testing
class SimpleFastqReader {
private:
    std::ifstream file;
    
public:
    SimpleFastqReader(const std::string& filename) : file(filename) {
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    }
    
    bool read_batch(ReadBatch& batch, size_t max_reads = 10000) {
        batch.clear();
        batch.reserve(max_reads);
        
        std::string line;
        size_t count = 0;
        
        while (count < max_reads && std::getline(file, line)) {
            if (line.empty() || line[0] != '@') continue;
            
            std::string header = line.substr(1);
            
            if (!std::getline(file, line)) break;
            std::string sequence = line;
            
            if (!std::getline(file, line)) break;  // + line
            
            if (!std::getline(file, line)) break;
            std::string quality = line;
            
            batch.headers.push_back(header);
            batch.sequences.push_back(sequence);
            batch.qualities.push_back(quality);
            count++;
        }
        
        return batch.size() > 0;
    }
};

// Generate test sequences if no file provided
std::vector<std::string> generate_test_sequences(size_t count = 1000) {
    std::vector<std::string> sequences;
    const char* bases = "ACGT";
    
    for (size_t i = 0; i < count; i++) {
        std::string seq;
        size_t length = 100 + (rand() % 150);  // 100-250bp
        
        for (size_t j = 0; j < length; j++) {
            seq += bases[rand() % 4];
        }
        
        sequences.push_back(seq);
    }
    
    return sequences;
}

int main(int argc, char** argv) {
    try {
        std::vector<std::string> sequences;
        
        if (argc == 2) {
            // Read from FASTQ file
            std::string filename = argv[1];
            std::cout << "Reading sequences from " << filename << "...\n";
            
            SimpleFastqReader reader(filename);
            ReadBatch batch;
            
            if (reader.read_batch(batch, 10000)) {
                sequences = batch.sequences;
                std::cout << "Read " << sequences.size() << " sequences\n";
            } else {
                std::cerr << "No sequences found in file\n";
                return 1;
            }
        } else {
            // Generate test sequences
            std::cout << "Generating 1000 test sequences...\n";
            sequences = generate_test_sequences(1000);
        }
        
        // Test minimizer extraction
        std::cout << "\nTesting minimizer extraction (k=31, m=15)...\n";
        
        MinimizerExtractor extractor(31, 15);
        
        auto start = std::chrono::high_resolution_clock::now();
        auto minimizers = extractor.extract_minimizers(sequences);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        // Print statistics
        size_t total_minimizers = 0;
        size_t max_minimizers = 0;
        size_t min_minimizers = SIZE_MAX;
        
        for (const auto& seq_minimizers : minimizers) {
            total_minimizers += seq_minimizers.size();
            max_minimizers = std::max(max_minimizers, seq_minimizers.size());
            min_minimizers = std::min(min_minimizers, seq_minimizers.size());
        }
        
        std::cout << "\nResults:\n";
        std::cout << "- Processed " << sequences.size() << " sequences\n";
        std::cout << "- Total minimizers: " << total_minimizers << "\n";
        std::cout << "- Average minimizers per sequence: " 
                  << (double)total_minimizers / sequences.size() << "\n";
        std::cout << "- Min minimizers in a sequence: " << min_minimizers << "\n";
        std::cout << "- Max minimizers in a sequence: " << max_minimizers << "\n";
        std::cout << "- Processing time: " << duration.count() << " ms\n";
        std::cout << "- Throughput: " << (sequences.size() * 1000.0 / duration.count()) 
                  << " sequences/second\n";
        
        // Print first few minimizers as examples
        if (minimizers.size() > 0 && minimizers[0].size() > 0) {
            std::cout << "\nFirst sequence minimizers (showing first 5):\n";
            for (size_t i = 0; i < std::min(size_t(5), minimizers[0].size()); i++) {
                std::cout << "  Hash: " << minimizers[0][i].hash
                          << ", Position: " << minimizers[0][i].position
                          << ", Reverse: " << (minimizers[0][i].is_reverse ? "yes" : "no")
                          << "\n";
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}