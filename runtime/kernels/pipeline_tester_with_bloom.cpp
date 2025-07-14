// pipeline_tester_with_bloom.cpp
//
// Enhanced pipeline tester with bloom filter pre-screening capability
// Combines sample sheet parsing, FASTQ reading, and bloom filter screening
//
// To Compile:
//   make pipeline_tester_bloom
//
// To Run:
//   ./pipeline_tester_bloom samples.csv                    # Without bloom filter (default)
//   ./pipeline_tester_bloom samples.csv --bloom-filter    # With bloom filter screening
//   ./pipeline_tester_bloom samples.csv --bloom-filter combined_bloom.bin  # Custom bloom filter
//

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <map>
#include <numeric>
#include <chrono>
#include <deque>
#include <cstring>
#include <cuda_runtime.h>

#include "sample_csv_parser.h"
#include "fastq_reader.h"

// Minimizer configuration (must match bloom filter)
const int KMER_LENGTH = 15;
const int WINDOW_SIZE = 10;
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
    int load_bloom_filter(void* filter, const char* filename);
    int bloom_filter_screen_reads_with_rc(void* filter, const char* d_reads, const int* d_read_lengths,
                                          const int* d_read_offsets, int num_reads, bool* d_read_passes,
                                          int* d_kmers_found, int min_kmers_threshold, bool check_rc,
                                          int* d_debug_stats);
}

// Test function signature
using TestFunc = std::function<bool(const BioGPU::SampleInfo&, bool, const std::string&)>;

struct TestCase {
    std::string name;
    TestFunc function;
};

// Global statistics
struct BloomStatistics {
    long total_reads = 0;
    long reads_passed = 0;
    long total_minimizers = 0;
    long minimizers_in_bloom = 0;
    double processing_time_ms = 0;
    
    void reset() {
        total_reads = 0;
        reads_passed = 0;
        total_minimizers = 0;
        minimizers_in_bloom = 0;
        processing_time_ms = 0;
    }
    
    void print() const {
        std::cout << "\n=== Bloom Filter Statistics ===" << std::endl;
        std::cout << "Total reads processed: " << total_reads << std::endl;
        std::cout << "Reads passed filter: " << reads_passed 
                  << " (" << (total_reads > 0 ? 100.0 * reads_passed / total_reads : 0) << "%)" << std::endl;
        std::cout << "Total minimizers extracted: " << total_minimizers << std::endl;
        std::cout << "Minimizers found in bloom filter: " << minimizers_in_bloom 
                  << " (" << (total_minimizers > 0 ? 100.0 * minimizers_in_bloom / total_minimizers : 0) << "%)" << std::endl;
        std::cout << "Processing time: " << processing_time_ms << " ms" << std::endl;
        std::cout << "==============================" << std::endl;
    }
};

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
            return UINT64_MAX;
        }
        encoded = (encoded << 2) | base;
    }
    
    return encoded;
}

// Extract minimizers from a sequence
int extract_minimizers_count(const std::string& sequence) {
    if (sequence.length() < KMER_LENGTH + WINDOW_SIZE - 1) {
        return 0;
    }
    
    int count = 0;
    std::deque<std::pair<uint64_t, size_t>> window;
    
    for (size_t i = 0; i <= sequence.length() - KMER_LENGTH; i++) {
        uint64_t kmer = encode_kmer(sequence, i);
        if (kmer == UINT64_MAX) continue;
        
        while (!window.empty() && window.front().second <= i - WINDOW_SIZE) {
            window.pop_front();
        }
        
        while (!window.empty() && window.back().first > kmer) {
            window.pop_back();
        }
        
        window.push_back({kmer, i});
        
        if (i >= WINDOW_SIZE - 1) {
            count++;
        }
    }
    
    return count;
}

// Process reads with bloom filter
bool process_batch_with_bloom(void* bloom_filter, 
                             const std::vector<FastqRecord>& batch_r1,
                             const std::vector<FastqRecord>& batch_r2,
                             BloomStatistics& stats) {
    if (batch_r1.empty()) return true;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int num_reads = batch_r1.size();
    stats.total_reads += num_reads;
    
    // Count minimizers
    for (int i = 0; i < num_reads; i++) {
        stats.total_minimizers += extract_minimizers_count(batch_r1[i].sequence);
        if (!batch_r2.empty()) {
            stats.total_minimizers += extract_minimizers_count(batch_r2[i].sequence);
        }
    }
    
    // Prepare batch for GPU
    ReadBatch gpu_batch_r1 = prepare_batch_for_gpu(batch_r1);
    ReadBatch gpu_batch_r2;
    if (!batch_r2.empty()) {
        gpu_batch_r2 = prepare_batch_for_gpu(batch_r2);
    }
    
    // Allocate GPU memory
    char* d_reads;
    int* d_read_lengths;
    int* d_read_offsets;
    bool* d_read_passes;
    int* d_kmers_found;
    
    int total_bases = gpu_batch_r1.total_bases + gpu_batch_r2.total_bases;
    int total_reads = num_reads * (batch_r2.empty() ? 1 : 2);
    
    cudaMalloc(&d_reads, total_bases);
    cudaMalloc(&d_read_lengths, total_reads * sizeof(int));
    cudaMalloc(&d_read_offsets, total_reads * sizeof(int));
    cudaMalloc(&d_read_passes, total_reads * sizeof(bool));
    cudaMalloc(&d_kmers_found, total_reads * sizeof(int));
    
    // Copy data to GPU
    if (!batch_r2.empty()) {
        // Interleave R1 and R2 reads
        std::vector<int> lengths, offsets;
        std::vector<char> sequences(total_bases);
        int current_offset = 0;
        
        for (int i = 0; i < num_reads; i++) {
            // R1
            lengths.push_back(gpu_batch_r1.lengths[i]);
            offsets.push_back(current_offset);
            memcpy(sequences.data() + current_offset, 
                   gpu_batch_r1.sequences + gpu_batch_r1.offsets[i], 
                   gpu_batch_r1.lengths[i]);
            current_offset += gpu_batch_r1.lengths[i];
            
            // R2
            lengths.push_back(gpu_batch_r2.lengths[i]);
            offsets.push_back(current_offset);
            memcpy(sequences.data() + current_offset, 
                   gpu_batch_r2.sequences + gpu_batch_r2.offsets[i], 
                   gpu_batch_r2.lengths[i]);
            current_offset += gpu_batch_r2.lengths[i];
        }
        
        cudaMemcpy(d_reads, sequences.data(), total_bases, cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, lengths.data(), total_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, offsets.data(), total_reads * sizeof(int), cudaMemcpyHostToDevice);
    } else {
        cudaMemcpy(d_reads, gpu_batch_r1.sequences, gpu_batch_r1.total_bases, cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, gpu_batch_r1.lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, gpu_batch_r1.offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice);
    }
    
    // Screen with bloom filter
    int min_kmers_threshold = 1;  // Very sensitive
    bool check_rc = true;  // Check reverse complement for R2
    
    bloom_filter_screen_reads_with_rc(bloom_filter, d_reads, d_read_lengths, d_read_offsets,
                                     total_reads, d_read_passes, d_kmers_found, 
                                     min_kmers_threshold, check_rc, nullptr);
    
    // Get results
    bool* h_read_passes = new bool[total_reads];
    std::vector<int> h_kmers_found(total_reads);
    cudaMemcpy(h_read_passes, d_read_passes, total_reads * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_kmers_found.data(), d_kmers_found, total_reads * sizeof(int), cudaMemcpyDeviceToHost);
    
    // Count passed reads and minimizers found
    if (!batch_r2.empty()) {
        // For paired-end, both reads must pass
        for (int i = 0; i < num_reads; i++) {
            if (h_read_passes[i*2] || h_read_passes[i*2+1]) {
                stats.reads_passed++;
            }
            stats.minimizers_in_bloom += h_kmers_found[i*2] + h_kmers_found[i*2+1];
        }
    } else {
        for (int i = 0; i < num_reads; i++) {
            if (h_read_passes[i]) {
                stats.reads_passed++;
            }
            stats.minimizers_in_bloom += h_kmers_found[i];
        }
    }
    
    // Cleanup
    cudaFree(d_reads);
    cudaFree(d_read_lengths);
    cudaFree(d_read_offsets);
    cudaFree(d_read_passes);
    cudaFree(d_kmers_found);
    
    delete[] h_read_passes;
    delete[] gpu_batch_r1.sequences;
    delete[] gpu_batch_r1.lengths;
    delete[] gpu_batch_r1.offsets;
    if (!batch_r2.empty()) {
        delete[] gpu_batch_r2.sequences;
        delete[] gpu_batch_r2.lengths;
        delete[] gpu_batch_r2.offsets;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.processing_time_ms += std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    return true;
}

// Enhanced FASTQ reader test with optional bloom filter
bool test_fastq_reader_with_bloom(const BioGPU::SampleInfo& sample, 
                                 bool use_bloom_filter, 
                                 const std::string& bloom_file) {
    std::cout << "  [RUNNING] test_fastq_reader_with_bloom for sample: " << sample.sample_name << std::endl;
    
    const int BATCH_SIZE = 1000;  // Larger batch for bloom filter efficiency
    long total_reads = 0;
    long total_bases = 0;
    bool success = true;
    
    BloomStatistics bloom_stats;
    void* bloom_filter = nullptr;
    
    // Load bloom filter if requested
    if (use_bloom_filter) {
        std::cout << "  [INFO] Loading bloom filter: " << bloom_file << std::endl;
        bloom_filter = create_bloom_filter(KMER_LENGTH);
        if (!bloom_filter || load_bloom_filter(bloom_filter, bloom_file.c_str()) != 0) {
            std::cerr << "  [ERROR] Failed to load bloom filter" << std::endl;
            return false;
        }
        std::cout << "  [INFO] Bloom filter loaded successfully" << std::endl;
    }
    
    if (sample.isPairedEnd()) {
        gzFile r1_file = gzopen(sample.read1_path.c_str(), "r");
        gzFile r2_file = gzopen(sample.read2_path.c_str(), "r");
        
        if (!r1_file || !r2_file) {
            std::cerr << "    [ERROR] Could not open paired-end files" << std::endl;
            if(r1_file) gzclose(r1_file);
            if(r2_file) gzclose(r2_file);
            if (bloom_filter) destroy_bloom_filter(bloom_filter);
            return false;
        }
        
        while (true) {
            std::vector<FastqRecord> batch_r1, batch_r2;
            if (!read_batch_paired(r1_file, r2_file, batch_r1, batch_r2, BATCH_SIZE)) {
                break;
            }
            if (batch_r1.empty()) break;
            
            total_reads += batch_r1.size();
            for(const auto& rec : batch_r1) total_bases += rec.sequence.length();
            for(const auto& rec : batch_r2) total_bases += rec.sequence.length();
            
            if (use_bloom_filter) {
                process_batch_with_bloom(bloom_filter, batch_r1, batch_r2, bloom_stats);
            }
        }
        
        gzclose(r1_file);
        gzclose(r2_file);
        
    } else {
        gzFile file = gzopen(sample.read1_path.c_str(), "r");
        if (!file) {
            std::cerr << "    [ERROR] Could not open single-end file" << std::endl;
            if (bloom_filter) destroy_bloom_filter(bloom_filter);
            return false;
        }
        
        while (true) {
            std::vector<FastqRecord> batch;
            if (!read_batch_single(file, batch, BATCH_SIZE)) {
                break;
            }
            if (batch.empty()) break;
            
            total_reads += batch.size();
            for(const auto& rec : batch) total_bases += rec.sequence.length();
            
            if (use_bloom_filter) {
                process_batch_with_bloom(bloom_filter, batch, {}, bloom_stats);
            }
        }
        gzclose(file);
    }
    
    // Print results
    if (total_reads > 0) {
        std::cout << "    [SUCCESS] Read " << total_reads << " records (" << total_bases << " bases)" << std::endl;
        
        if (use_bloom_filter) {
            bloom_stats.print();
        }
    } else {
        std::cerr << "    [FAIL] No reads were imported from the files" << std::endl;
        success = false;
    }
    
    if (bloom_filter) {
        destroy_bloom_filter(bloom_filter);
    }
    
    return success;
}

// Test runner class
class PipelineTester {
private:
    std::vector<TestCase> tests;
    BioGPU::SampleCSVParser parser;
    bool use_bloom_filter = false;
    std::string bloom_filter_file = "combined_bloom_filter.bin";
    
public:
    void register_test(const std::string& name, TestFunc func) {
        tests.push_back({name, func});
    }
    
    void set_bloom_filter(bool use_bloom, const std::string& bloom_file = "") {
        use_bloom_filter = use_bloom;
        if (!bloom_file.empty()) {
            bloom_filter_file = bloom_file;
        }
    }
    
    bool run(const std::string& csv_path) {
        std::cout << "--- Initializing Pipeline Test Runner ---" << std::endl;
        if (use_bloom_filter) {
            std::cout << "Bloom filter: ENABLED (file: " << bloom_filter_file << ")" << std::endl;
        } else {
            std::cout << "Bloom filter: DISABLED" << std::endl;
        }
        
        if (!parser.parseFile(csv_path, true)) {
            std::cerr << "[FATAL] Failed to parse sample CSV: " << csv_path << std::endl;
            return false;
        }
        parser.printSummary();
        
        int total_passed = 0;
        int total_failed = 0;
        
        for (size_t i = 0; i < parser.getSampleCount(); ++i) {
            const BioGPU::SampleInfo* sample = parser.getSample(i);
            if (!sample) continue;
            
            std::cout << "\n--- Testing Sample: " << sample->sample_name << " ---" << std::endl;
            
            if (!parser.validateSample(*sample)) {
                std::cerr << "  [SKIPPING] Sample files not found or invalid." << std::endl;
                continue;
            }
            
            for (const auto& test : tests) {
                std::cout << "-> Running Test: " << test.name << std::endl;
                bool result = test.function(*sample, use_bloom_filter, bloom_filter_file);
                if (result) {
                    total_passed++;
                } else {
                    total_failed++;
                }
            }
        }
        
        std::cout << "\n--- Test Suite Finished ---" << std::endl;
        std::cout << "Total Tests Passed: " << total_passed << std::endl;
        std::cout << "Total Tests Failed: " << total_failed << std::endl;
        std::cout << "--------------------------" << std::endl;
        
        return total_failed == 0;
    }
};

// Main function
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <samples.csv> [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --bloom-filter           Enable bloom filter screening (default: disabled)" << std::endl;
        std::cerr << "  --bloom-filter <file>    Enable bloom filter with custom file" << std::endl;
        return 1;
    }
    
    std::string csv_path = argv[1];
    bool use_bloom = false;
    std::string bloom_file = "combined_bloom_filter.bin";
    
    // Parse command line options
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--bloom-filter") {
            use_bloom = true;
            // Check if next argument is a file path
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                bloom_file = argv[++i];
            }
        }
    }
    
    PipelineTester tester;
    tester.set_bloom_filter(use_bloom, bloom_file);
    
    // Register tests
    tester.register_test("FASTQ Reader with Bloom Filter Test", test_fastq_reader_with_bloom);
    
    // Run the test suite
    bool success = tester.run(csv_path);
    
    return success ? 0 : 1;
}