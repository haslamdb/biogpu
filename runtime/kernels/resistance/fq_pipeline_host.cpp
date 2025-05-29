#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <zlib.h>
#include <hdf5.h>
#include <cuda_runtime.h>
#include "fq_mutation_detector.cuh"

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
        file = gzopen(filename.c_str(), "r");
        is_open = (file != NULL);
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
        // Allocate GPU memory for batch processing
        size_t read_buffer_size = batch_size * max_read_length;
        
        cudaMalloc(&d_reads_r1, read_buffer_size);
        cudaMalloc(&d_reads_r2, read_buffer_size);
        cudaMalloc(&d_lengths_r1, batch_size * sizeof(int));
        cudaMalloc(&d_lengths_r2, batch_size * sizeof(int));
        cudaMalloc(&d_offsets_r1, batch_size * sizeof(int));
        cudaMalloc(&d_offsets_r2, batch_size * sizeof(int));
        
        // Allocate memory for results
        cudaMalloc(&d_candidates, batch_size * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch));
        cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t));
        cudaMalloc(&d_results, batch_size * MAX_CANDIDATES_PER_READ * sizeof(AlignmentResult));
        cudaMalloc(&d_result_count, sizeof(uint32_t));
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
        detector.loadIndex(index_path.c_str());
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
        FastqReader reader1(r1_path);
        FastqReader reader2(r2_path);
        
        if (!reader1.isOpen() || !reader2.isOpen()) {
            std::cerr << "Failed to open input files" << std::endl;
            return;
        }
        
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
            
            for (int i = 0; i < batch_size; i++) {
                FastqRecord rec1, rec2;
                bool got1 = reader1.readRecord(rec1);
                bool got2 = reader2.readRecord(rec2);
                
                if (!got1 || !got2) break;
                
                batch_r1.push_back(rec1);
                batch_r2.push_back(rec2);
            }
            
            if (batch_r1.empty()) break;
            
            // Prepare batches
            ReadBatch batch1 = prepareBatch(batch_r1);
            ReadBatch batch2 = prepareBatch(batch_r2);
            
            // Transfer to GPU
            cudaMemcpy(d_reads_r1, batch1.sequences, batch1.total_bases, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lengths_r1, batch1.lengths, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_offsets_r1, batch1.offsets, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice);
            
            cudaMemcpy(d_reads_r2, batch2.sequences, batch2.total_bases, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lengths_r2, batch2.lengths, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_offsets_r2, batch2.offsets, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice);
            
            // Reset result counter
            cudaMemset(d_result_count, 0, sizeof(uint32_t));
            
            // Launch kernels for R1
            // Stage 1: K-mer filtering
            launch_kmer_filter(
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                detector, batch1.num_reads,
                d_candidates, d_candidate_counts
            );
            
            // Stage 2: Position-weighted alignment
            launch_position_weighted_alignment(
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                d_candidates, d_candidate_counts,
                detector, batch1.num_reads,
                d_results, d_result_count
            );
            
            // Get results
            uint32_t num_results;
            cudaMemcpy(&num_results, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost);
            
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
    
    // Check CUDA device
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
    
    // Run pipeline
    FQResistancePipeline pipeline;
    pipeline.loadIndex(index_path);
    pipeline.processPairedReads(r1_path, r2_path, output_path);
    
    return 0;
}