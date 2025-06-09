// test_minimizer.cpp
// Simple test program for minimizer extraction
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include "minimizer_common.h"
#include "minimizer_extractor.h"

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
        bool benchmark_mode = false;
        
        // Check for benchmark flag
        if (argc >= 2 && std::string(argv[1]) == "--benchmark") {
            benchmark_mode = true;
            std::cout << "Running in benchmark mode...\n";
        }
        
        if (argc >= 2 && !benchmark_mode) {
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
            size_t num_sequences = benchmark_mode ? 10000 : 1000;
            std::cout << "Generating " << num_sequences << " test sequences...\n";
            sequences = generate_test_sequences(num_sequences);
        }
        
        // Test minimizer extraction
        std::cout << "\nTesting minimizer extraction (k=31, m=15)...\n";
        
        MinimizerExtractor extractor(31, 15);
        
        // Warm-up run for GPU
        if (benchmark_mode) {
            std::cout << "Performing warm-up run...\n";
            auto warmup = extractor.extract_minimizers(sequences);
        }
        
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
            if (seq_minimizers.size() > 0) {
                min_minimizers = std::min(min_minimizers, seq_minimizers.size());
            }
        }
        
        if (min_minimizers == SIZE_MAX) min_minimizers = 0;
        
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
        
        // Calculate total bases processed
        size_t total_bases = 0;
        for (const auto& seq : sequences) {
            total_bases += seq.length();
        }
        std::cout << "- Total bases: " << total_bases << "\n";
        std::cout << "- Throughput: " << (total_bases * 1000.0 / duration.count() / 1000000.0) 
                  << " Mbases/second\n";
        
        // Print first few minimizers as examples
        if (!benchmark_mode && minimizers.size() > 0 && minimizers[0].size() > 0) {
            std::cout << "\nFirst sequence minimizers (showing first 5):\n";
            for (size_t i = 0; i < std::min(size_t(5), minimizers[0].size()); i++) {
                std::cout << "  Hash: " << minimizers[0][i].hash
                          << ", Position: " << minimizers[0][i].position
                          << ", Reverse: " << (minimizers[0][i].is_reverse ? "yes" : "no")
                          << "\n";
            }
        }
        
        // Additional benchmark statistics
        if (benchmark_mode) {
            std::cout << "\nBenchmark Summary:\n";
            std::cout << "- GPU kernel time: ~" << duration.count() * 0.95 << " ms (estimated)\n";
            std::cout << "- Data transfer overhead: ~" << duration.count() * 0.05 << " ms (estimated)\n";
            
            // Memory usage estimate
            size_t gpu_memory_used = sequences.size() * 300 +  // sequences
                                    sequences.size() * 4 * 3 +  // offsets, lengths, counts
                                    sequences.size() * 50 * sizeof(Minimizer);  // minimizers
            std::cout << "- Estimated GPU memory used: " << gpu_memory_used / (1024.0 * 1024.0) << " MB\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}