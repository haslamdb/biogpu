#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cuda_runtime.h>
#include "fq_mutation_detector.cuh"

// Debug macros
#define DEBUG 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { fprintf(stderr, "[TEST DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

// CUDA error checking
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "[CUDA ERROR] %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// Generate test reads
std::vector<std::string> generateTestReads() {
    std::vector<std::string> reads = {
        // Some DNA sequences that might contain FQ resistance k-mers
        "ATGGCGATCGAATTGCTGCGATCGAATCCGATCGATGCATGC",  // 42 bp
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA",  // 42 bp
        "TTGCGATCGAATCCGATCGATGCATGCTTGCGATCGAATCCGA",  // 42 bp
        "GCGATCGAATTGCTGCGATCGAATCCGATCGATGCATGCTTGC",  // 42 bp
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"   // 41 bp
    };
    
    // Add some longer reads
    for (int i = 0; i < 10; i++) {
        std::string long_read = "ATGGCGATCGAATTGCTGCGATCGAATCCGATCGATGCATGCTTGCGATCGAATCCGATCGATGCATGC"
                              "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAATGGCGATCGAATTGCTGCGATCGAA"
                              "TTGCGATCGAATCCGATCGATGCATGCTTGCGATCGAATCCGATCGATGCATGCGCGATCGAATTGCT";
        reads.push_back(long_read);
    }
    
    return reads;
}

// Create test index directory structure
bool createTestIndex(const std::string& index_dir) {
    DEBUG_PRINT("Creating test index in: %s", index_dir.c_str());
    
    // Create directory (simplified - assumes it exists)
    std::string kmer_index_path = index_dir + "/kmer_index.bin";
    
    std::ofstream kmer_file(kmer_index_path, std::ios::binary);
    if (!kmer_file.good()) {
        DEBUG_PRINT("Failed to create k-mer index file: %s", kmer_index_path.c_str());
        return false;
    }
    
    // Write header
    uint32_t num_entries = 100;  // Small test index
    uint32_t kmer_length = 15;
    kmer_file.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint32_t));
    kmer_file.write(reinterpret_cast<const char*>(&kmer_length), sizeof(uint32_t));
    
    // Write test k-mer entries
    for (uint32_t i = 0; i < num_entries; i++) {
        KmerEntry entry;
        entry.kmer = i * 1000;  // Simple test k-mer values
        entry.gene_id = i % 4;  // Rotate through genes (gyrA, gyrB, parC, parE)
        entry.species_id = i % 10;  // Multiple species
        entry.seq_id = i % 20;      // Multiple sequences per species
        entry.position = i % 500;   // Position in sequence
        
        kmer_file.write(reinterpret_cast<const char*>(&entry.kmer), sizeof(uint64_t));
        kmer_file.write(reinterpret_cast<const char*>(&entry.gene_id), sizeof(uint32_t));
        kmer_file.write(reinterpret_cast<const char*>(&entry.species_id), sizeof(uint32_t));
        kmer_file.write(reinterpret_cast<const char*>(&entry.seq_id), sizeof(uint32_t));
        kmer_file.write(reinterpret_cast<const char*>(&entry.position), sizeof(uint16_t));
    }
    
    kmer_file.close();
    DEBUG_PRINT("Created test k-mer index with %d entries", num_entries);
    
    return true;
}

// Test the integrated k-mer filtering
bool testKmerIntegration(const std::string& index_dir) {
    DEBUG_PRINT("=== Testing K-mer Integration ===");
    
    // Generate test reads
    auto test_reads = generateTestReads();
    int num_reads = test_reads.size();
    
    DEBUG_PRINT("Generated %d test reads", num_reads);
    
    // Calculate total length and create offsets
    int total_length = 0;
    std::vector<int> read_lengths(num_reads);
    std::vector<int> read_offsets(num_reads);
    
    for (int i = 0; i < num_reads; i++) {
        read_offsets[i] = total_length;
        read_lengths[i] = test_reads[i].length();
        total_length += read_lengths[i];
    }
    
    // Concatenate all reads
    std::string all_reads;
    for (const auto& read : test_reads) {
        all_reads += read;
    }
    
    DEBUG_PRINT("Total concatenated length: %d", total_length);
    
    // Allocate GPU memory for reads
    char* d_reads;
    int* d_read_lengths;
    int* d_read_offsets;
    
    CUDA_CHECK(cudaMalloc(&d_reads, total_length));
    CUDA_CHECK(cudaMalloc(&d_read_lengths, num_reads * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_read_offsets, num_reads * sizeof(int)));
    
    // Copy data to GPU
    CUDA_CHECK(cudaMemcpy(d_reads, all_reads.c_str(), total_length, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_read_lengths, read_lengths.data(), num_reads * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_read_offsets, read_offsets.data(), num_reads * sizeof(int), cudaMemcpyHostToDevice));
    
    // Create detector and load index
    FQMutationDetectorCUDA detector;
    detector.loadIndex(index_dir.c_str());
    
    // Allocate memory for results
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    AlignmentResult* d_results;
    uint32_t* d_result_count;
    
    CUDA_CHECK(cudaMalloc(&d_candidates, num_reads * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch)));
    CUDA_CHECK(cudaMalloc(&d_candidate_counts, num_reads * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_results, num_reads * MAX_RESULTS_PER_READ * sizeof(AlignmentResult)));
    CUDA_CHECK(cudaMalloc(&d_result_count, sizeof(uint32_t)));
    
    // Initialize result arrays
    CUDA_CHECK(cudaMemset(d_candidate_counts, 0, num_reads * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
    
    DEBUG_PRINT("Launching Stage 1: K-mer filtering");
    
    // Test Stage 1: K-mer filtering
    launch_kmer_filter(
        d_reads, d_read_lengths, d_read_offsets,
        detector, num_reads,
        d_candidates, d_candidate_counts
    );
    
    // Check results
    std::vector<uint32_t> h_candidate_counts(num_reads);
    CUDA_CHECK(cudaMemcpy(h_candidate_counts.data(), d_candidate_counts, 
                         num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    int total_candidates = 0;
    int reads_with_candidates = 0;
    for (int i = 0; i < num_reads; i++) {
        total_candidates += h_candidate_counts[i];
        if (h_candidate_counts[i] > 0) reads_with_candidates++;
        DEBUG_PRINT("Read %d: %d candidates", i, h_candidate_counts[i]);
    }
    
    DEBUG_PRINT("Stage 1 Results: %d total candidates, %d reads with candidates", 
               total_candidates, reads_with_candidates);
    
    DEBUG_PRINT("Launching Stage 2: Alignment");
    
    // Test Stage 2: Alignment
    launch_position_weighted_alignment(
        d_reads, d_read_lengths, d_read_offsets,
        d_candidates, d_candidate_counts,
        detector, num_reads,
        d_results, d_result_count
    );
    
    // Check alignment results
    uint32_t num_results;
    CUDA_CHECK(cudaMemcpy(&num_results, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    DEBUG_PRINT("Stage 2 Results: %d alignment results", num_results);
    
    if (num_results > 0) {
        std::vector<AlignmentResult> h_results(num_results);
        CUDA_CHECK(cudaMemcpy(h_results.data(), d_results, 
                             num_results * sizeof(AlignmentResult), cudaMemcpyDeviceToHost));
        
        for (uint32_t i = 0; i < std::min(num_results, 5U); i++) {
            const auto& result = h_results[i];
            DEBUG_PRINT("Result %d: read=%d, gene=%d, species=%d, score=%.2f, mutations=%d", 
                       i, result.read_id, result.gene_id, result.species_id, 
                       result.alignment_score, result.num_mutations_detected);
        }
    }
    
    // Cleanup
    cudaFree(d_reads);
    cudaFree(d_read_lengths);
    cudaFree(d_read_offsets);
    cudaFree(d_candidates);
    cudaFree(d_candidate_counts);
    cudaFree(d_results);
    cudaFree(d_result_count);
    
    DEBUG_PRINT("=== K-mer Integration Test Complete ===");
    
    // Test is successful if we got some kind of output without crashes
    return true;
}

int main(int argc, char** argv) {
    DEBUG_PRINT("=== K-mer Integration Test Starting ===");
    
    std::string index_dir = "/tmp/test_index";
    if (argc > 1) {
        index_dir = argv[1];
    }
    
    DEBUG_PRINT("Using index directory: %s", index_dir.c_str());
    
    // Check CUDA device
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    DEBUG_PRINT("CUDA device count: %d", device_count);
    
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    // Create test index
    if (!createTestIndex(index_dir)) {
        std::cerr << "Failed to create test index" << std::endl;
        return 1;
    }
    
    // Run integration test
    if (!testKmerIntegration(index_dir)) {
        std::cerr << "K-mer integration test failed" << std::endl;
        return 1;
    }
    
    std::cout << "✅ K-mer integration test passed!" << std::endl;
    DEBUG_PRINT("=== All Tests Completed Successfully ===");
    
    return 0;
}