#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <zlib.h>
#include <hdf5.h>
#include <cuda_runtime.h>
#include "fq_mutation_detector.cuh"

// Debug macros
#define DEBUG 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { fprintf(stderr, "[DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); fflush(stderr); }

// CUDA error checking
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "[CUDA ERROR] %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// FASTQ record structure
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Read batch for GPU processing
struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

class FastqReader {
private:
    gzFile file;
    bool is_open;
    
public:
    FastqReader(const std::string& filename) : is_open(false) {
        DEBUG_PRINT("Opening FASTQ file: %s", filename.c_str());
        file = gzopen(filename.c_str(), "r");
        is_open = (file != NULL);
        if (!is_open) {
            DEBUG_PRINT("ERROR: Failed to open FASTQ file: %s", filename.c_str());
        } else {
            DEBUG_PRINT("Successfully opened FASTQ file");
        }
    }
    
    ~FastqReader() {
        if (is_open) {
            gzclose(file);
        }
    }
    
    bool readRecord(FastqRecord& record) {
        if (!is_open) return false;
        
        const int buffer_size = 1024;
        char buffer[buffer_size];
        
        // Read header
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.header = std::string(buffer);
        if (record.header.empty() || record.header[0] != '@') return false;
        
        // Read sequence
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.sequence = std::string(buffer);
        record.sequence.pop_back(); // Remove newline
        
        // Read plus line
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        
        // Read quality
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.quality = std::string(buffer);
        record.quality.pop_back(); // Remove newline
        
        return true;
    }
    
    bool isOpen() const { return is_open; }
};

class FQResistancePipeline {
private:
    // GPU memory for reads
    char* d_reads_r1;
    char* d_reads_r2;
    int* d_lengths_r1;
    int* d_lengths_r2;
    int* d_offsets_r1;
    int* d_offsets_r2;
    
    // Results
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    AlignmentResult* d_results;
    uint32_t* d_result_count;
    
    // Batch parameters
    const int batch_size = 10000;
    const int max_read_length = 300;
    
    FQMutationDetectorCUDA detector;
    
public:
    FQResistancePipeline() {
        DEBUG_PRINT("Initializing FQ Resistance Pipeline");
        DEBUG_PRINT("Batch size: %d, Max read length: %d", batch_size, max_read_length);
        
        // Allocate GPU memory for batch processing
        size_t read_buffer_size = batch_size * max_read_length;
        DEBUG_PRINT("Allocating GPU memory: read buffer size = %zu bytes", read_buffer_size);
        
        CUDA_CHECK(cudaMalloc(&d_reads_r1, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_reads_r2, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_lengths_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_lengths_r2, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r2, batch_size * sizeof(int)));
        
        // Allocate memory for results
        size_t candidates_size = batch_size * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch);
        size_t results_size = batch_size * MAX_CANDIDATES_PER_READ * sizeof(AlignmentResult);
        DEBUG_PRINT("Allocating GPU memory: candidates = %zu bytes, results = %zu bytes", 
                    candidates_size, results_size);
        
        CUDA_CHECK(cudaMalloc(&d_candidates, candidates_size));
        CUDA_CHECK(cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_results, results_size));
        CUDA_CHECK(cudaMalloc(&d_result_count, sizeof(uint32_t)));
        
        DEBUG_PRINT("GPU memory allocation completed successfully");
    }
    
    ~FQResistancePipeline() {
        cudaFree(d_reads_r1);
        cudaFree(d_reads_r2);
        cudaFree(d_lengths_r1);
        cudaFree(d_lengths_r2);
        cudaFree(d_offsets_r1);
        cudaFree(d_offsets_r2);
        cudaFree(d_candidates);
        cudaFree(d_candidate_counts);
        cudaFree(d_results);
        cudaFree(d_result_count);
    }
    
    void loadIndex(const std::string& index_path) {
        DEBUG_PRINT("Loading index from: %s", index_path.c_str());
        detector.loadIndex(index_path.c_str());
        DEBUG_PRINT("Index loading completed (check for errors above)");
    }
    
    ReadBatch prepareBatch(const std::vector<FastqRecord>& records) {
        ReadBatch batch;
        batch.num_reads = records.size();
        
        // Calculate total size needed
        batch.total_bases = 0;
        for (const auto& record : records) {
            batch.total_bases += record.sequence.length();
        }
        
        // Allocate host memory
        batch.sequences = new char[batch.total_bases];
        batch.lengths = new int[batch.num_reads];
        batch.offsets = new int[batch.num_reads];
        
        // Copy sequences
        int offset = 0;
        for (int i = 0; i < batch.num_reads; i++) {
            const std::string& seq = records[i].sequence;
            memcpy(batch.sequences + offset, seq.c_str(), seq.length());
            batch.lengths[i] = seq.length();
            batch.offsets[i] = offset;
            offset += seq.length();
        }
        
        return batch;
    }
    
    void processPairedReads(const std::string& r1_path, const std::string& r2_path, 
                           const std::string& output_path) {
        DEBUG_PRINT("Starting paired read processing");
        DEBUG_PRINT("R1: %s", r1_path.c_str());
        DEBUG_PRINT("R2: %s", r2_path.c_str());
        DEBUG_PRINT("Output: %s", output_path.c_str());
        
        FastqReader reader1(r1_path);
        FastqReader reader2(r2_path);
        
        if (!reader1.isOpen() || !reader2.isOpen()) {
            std::cerr << "ERROR: Failed to open input files" << std::endl;
            DEBUG_PRINT("Reader1 open: %s, Reader2 open: %s", 
                       reader1.isOpen() ? "YES" : "NO",
                       reader2.isOpen() ? "YES" : "NO");
            return;
        }
        
        DEBUG_PRINT("Both FASTQ files opened successfully");
        
        std::ofstream output(output_path);
        output << "{\n";
        output << "  \"sample\": \"" << r1_path << "\",\n";
        output << "  \"mutations\": [\n";
        
        int total_reads = 0;
        int total_mutations = 0;
        bool first_mutation = true;
        
        while (true) {
            // Read batch of paired reads
            std::vector<FastqRecord> batch_r1, batch_r2;
            
            DEBUG_PRINT("Reading batch %d (up to %d reads)", total_reads/batch_size + 1, batch_size);
            
            for (int i = 0; i < batch_size; i++) {
                FastqRecord rec1, rec2;
                bool got1 = reader1.readRecord(rec1);
                bool got2 = reader2.readRecord(rec2);
                
                if (!got1 || !got2) {
                    if (i == 0) {
                        DEBUG_PRINT("No more reads available");
                    } else {
                        DEBUG_PRINT("Batch partially filled: %d reads", i);
                    }
                    break;
                }
                
                batch_r1.push_back(rec1);
                batch_r2.push_back(rec2);
            }
            
            if (batch_r1.empty()) {
                DEBUG_PRINT("Empty batch, ending processing");
                break;
            }
            
            DEBUG_PRINT("Batch contains %zu read pairs", batch_r1.size());
            
            // Prepare batches
            ReadBatch batch1 = prepareBatch(batch_r1);
            ReadBatch batch2 = prepareBatch(batch_r2);
            
            // Transfer to GPU
            DEBUG_PRINT("Transferring batch to GPU: %d reads, %d total bases (R1)", 
                       batch1.num_reads, batch1.total_bases);
            
            CUDA_CHECK(cudaMemcpy(d_reads_r1, batch1.sequences, batch1.total_bases, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lengths_r1, batch1.lengths, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_offsets_r1, batch1.offsets, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            
            CUDA_CHECK(cudaMemcpy(d_reads_r2, batch2.sequences, batch2.total_bases, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lengths_r2, batch2.lengths, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_offsets_r2, batch2.offsets, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            
            // Reset result counter
            CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
            
            DEBUG_PRINT("GPU transfer completed");
            
            // Launch kernels for R1
            DEBUG_PRINT("Launching Stage 1: K-mer filtering for %d reads", batch1.num_reads);
            launch_kmer_filter(
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                detector, batch1.num_reads,
                d_candidates, d_candidate_counts
            );
            CUDA_CHECK(cudaDeviceSynchronize());
            DEBUG_PRINT("Stage 1 completed");
            
            // Check candidate counts
            std::vector<uint32_t> h_candidate_counts(batch1.num_reads);
            CUDA_CHECK(cudaMemcpy(h_candidate_counts.data(), d_candidate_counts, 
                                 batch1.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
            int total_candidates = 0;
            for (int i = 0; i < batch1.num_reads; i++) {
                total_candidates += h_candidate_counts[i];
            }
            DEBUG_PRINT("Total candidates found: %d", total_candidates);
            
            // Stage 2: Position-weighted alignment
            DEBUG_PRINT("Launching Stage 2: Position-weighted alignment");
            launch_position_weighted_alignment(
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                d_candidates, d_candidate_counts,
                detector, batch1.num_reads,
                d_results, d_result_count
            );
            CUDA_CHECK(cudaDeviceSynchronize());
            DEBUG_PRINT("Stage 2 completed");
            
            // Get results
            uint32_t num_results;
            CUDA_CHECK(cudaMemcpy(&num_results, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            DEBUG_PRINT("Number of alignment results: %u", num_results);
            
            if (num_results > 0) {
                std::vector<AlignmentResult> results(num_results);
                cudaMemcpy(results.data(), d_results, num_results * sizeof(AlignmentResult), 
                          cudaMemcpyDeviceToHost);
                
                // Output mutations found
                for (const auto& result : results) {
                    if (result.num_mutations_detected > 0) {
                        if (!first_mutation) output << ",\n";
                        output << "    {\n";
                        output << "      \"read_pair\": " << (result.read_id + total_reads) << ",\n";
                        output << "      \"gene_id\": " << result.gene_id << ",\n";
                        output << "      \"species_id\": " << result.species_id << ",\n";
                        output << "      \"alignment_score\": " << result.alignment_score << ",\n";
                        output << "      \"identity\": " << result.identity << ",\n";
                        output << "      \"mutations_detected\": " << (int)result.num_mutations_detected << "\n";
                        output << "    }";
                        first_mutation = false;
                        total_mutations++;
                    }
                }
            }
            
            // TODO: Process R2 reads similarly
            
            total_reads += batch1.num_reads;
            
            // Cleanup batch memory
            delete[] batch1.sequences;
            delete[] batch1.lengths;
            delete[] batch1.offsets;
            delete[] batch2.sequences;
            delete[] batch2.lengths;
            delete[] batch2.offsets;
            
            // Progress update
            if (total_reads % 100000 == 0) {
                std::cout << "Processed " << total_reads << " read pairs..." << std::endl;
            }
        }
        
        output << "\n  ],\n";
        output << "  \"summary\": {\n";
        output << "    \"total_reads\": " << total_reads << ",\n";
        output << "    \"mutations_found\": " << total_mutations << "\n";
        output << "  }\n";
        output << "}\n";
        output.close();
        
        std::cout << "Processing complete. Total reads: " << total_reads 
                  << ", Mutations found: " << total_mutations << std::endl;
    }
};

// Main entry point
int main(int argc, char** argv) {
    DEBUG_PRINT("=== FQ Pipeline GPU Starting ===");
    DEBUG_PRINT("Command line: argc=%d", argc);
    for (int i = 0; i < argc; i++) {
        DEBUG_PRINT("  argv[%d] = %s", i, argv[i]);
    }
    
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <index_path> <reads_R1.fastq.gz> <reads_R2.fastq.gz> <output.json>" << std::endl;
        return 1;
    }
    
    std::string index_path = argv[1];
    std::string r1_path = argv[2];
    std::string r2_path = argv[3];
    std::string output_path = argv[4];
    
    // Get input paths from user if files are on adjacent drive
    std::cout << "=== Fluoroquinolone Resistance Detection Pipeline ===" << std::endl;
    std::cout << "Index: " << index_path << std::endl;
    std::cout << "R1 reads: " << r1_path << std::endl;
    std::cout << "R2 reads: " << r2_path << std::endl;
    std::cout << "Output: " << output_path << std::endl;
    
    // Check input files exist
    DEBUG_PRINT("Checking input files...");
    std::ifstream idx_check(index_path);
    if (!idx_check.good()) {
        std::cerr << "ERROR: Index file not found: " << index_path << std::endl;
        DEBUG_PRINT("Index file check failed");
        return 1;
    }
    idx_check.close();
    DEBUG_PRINT("Index file exists");
    
    std::ifstream r1_check(r1_path);
    if (!r1_check.good()) {
        std::cerr << "ERROR: R1 file not found: " << r1_path << std::endl;
        DEBUG_PRINT("R1 file check failed");
        return 1;
    }
    r1_check.close();
    DEBUG_PRINT("R1 file exists");
    
    std::ifstream r2_check(r2_path);
    if (!r2_check.good()) {
        std::cerr << "ERROR: R2 file not found: " << r2_path << std::endl;
        DEBUG_PRINT("R2 file check failed");
        return 1;
    }
    r2_check.close();
    DEBUG_PRINT("R2 file exists");
    
    // Check CUDA device
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    DEBUG_PRINT("CUDA device count: %d", device_count);
    
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
    
    // Run pipeline
    try {
        DEBUG_PRINT("Creating pipeline instance");
        FQResistancePipeline pipeline;
        
        DEBUG_PRINT("Loading index");
        pipeline.loadIndex(index_path);
        
        DEBUG_PRINT("Starting paired read processing");
        pipeline.processPairedReads(r1_path, r2_path, output_path);
        
        DEBUG_PRINT("=== FQ Pipeline GPU Completed Successfully ===");
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Exception caught: " << e.what() << std::endl;
        DEBUG_PRINT("Exception: %s", e.what());
        return 1;
    } catch (...) {
        std::cerr << "ERROR: Unknown exception caught" << std::endl;
        DEBUG_PRINT("Unknown exception caught");
        return 1;
    }
    
    return 0;
}