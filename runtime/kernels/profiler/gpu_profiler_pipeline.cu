// gpu_profiler_pipeline.cu - Fully GPU-accelerated profiling pipeline
#include <iostream>
#include <chrono>
#include <cuda_runtime.h>
#include "fastq_processing.h"
#include "minimizer_extractor.h"
#include "gpu_kmer_database.h"

using namespace biogpu;

class GPUProfilerPipeline {
private:
    std::unique_ptr<GPUKmerDatabase> database;
    std::unique_ptr<MinimizerExtractor> extractor;
    
    // GPU memory for batch processing
    uint64_t* d_minimizer_hashes = nullptr;
    uint32_t* d_minimizer_taxons = nullptr;
    uint32_t* d_minimizer_offsets = nullptr;
    uint32_t* d_minimizer_counts = nullptr;
    uint32_t* d_read_classifications = nullptr;
    float* d_confidence_scores = nullptr;
    
    size_t allocated_reads = 0;
    size_t allocated_minimizers = 0;
    
    void allocate_gpu_memory(size_t num_reads, size_t num_minimizers) {
        if (num_reads > allocated_reads) {
            if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
            if (d_minimizer_counts) cudaFree(d_minimizer_counts);
            if (d_read_classifications) cudaFree(d_read_classifications);
            if (d_confidence_scores) cudaFree(d_confidence_scores);
            
            cudaMalloc(&d_minimizer_offsets, (num_reads + 1) * sizeof(uint32_t));
            cudaMalloc(&d_minimizer_counts, num_reads * sizeof(uint32_t));
            cudaMalloc(&d_read_classifications, num_reads * sizeof(uint32_t));
            cudaMalloc(&d_confidence_scores, num_reads * sizeof(float));
            
            allocated_reads = num_reads;
        }
        
        if (num_minimizers > allocated_minimizers) {
            if (d_minimizer_hashes) cudaFree(d_minimizer_hashes);
            if (d_minimizer_taxons) cudaFree(d_minimizer_taxons);
            
            cudaMalloc(&d_minimizer_hashes, num_minimizers * sizeof(uint64_t));
            cudaMalloc(&d_minimizer_taxons, num_minimizers * sizeof(uint32_t));
            
            allocated_minimizers = num_minimizers;
        }
    }
    
public:
    GPUProfilerPipeline(const std::string& db_file, int k = 31, int m = 15) {
        database = std::make_unique<GPUKmerDatabase>();
        database->load_binary(db_file);
        
        extractor = std::make_unique<MinimizerExtractor>(k, m);
    }
    
    ~GPUProfilerPipeline() {
        if (d_minimizer_hashes) cudaFree(d_minimizer_hashes);
        if (d_minimizer_taxons) cudaFree(d_minimizer_taxons);
        if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
        if (d_minimizer_counts) cudaFree(d_minimizer_counts);
        if (d_read_classifications) cudaFree(d_read_classifications);
        if (d_confidence_scores) cudaFree(d_confidence_scores);
    }
    
    void process_paired_batch(const std::vector<std::string>& sequences1,
                             const std::vector<std::string>& sequences2,
                             std::vector<uint32_t>& classifications,
                             std::vector<float>& confidences) {
        // Process both forward and reverse reads together
        size_t num_pairs = sequences1.size();
        if (num_pairs == 0 || sequences1.size() != sequences2.size()) return;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Extract minimizers from both reads
        auto minimizer_results1 = extractor->extract_minimizers(sequences1);
        auto minimizer_results2 = extractor->extract_minimizers(sequences2);
        
        // Combine minimizers from paired reads
        std::vector<uint64_t> all_hashes;
        std::vector<uint32_t> offsets(num_pairs + 1);
        std::vector<uint32_t> counts(num_pairs);
        
        uint32_t current_offset = 0;
        for (size_t i = 0; i < num_pairs; i++) {
            offsets[i] = current_offset;
            
            // Add minimizers from both reads
            for (const auto& m : minimizer_results1[i]) {
                all_hashes.push_back(m.hash);
            }
            for (const auto& m : minimizer_results2[i]) {
                all_hashes.push_back(m.hash);
            }
            
            counts[i] = minimizer_results1[i].size() + minimizer_results2[i].size();
            current_offset += counts[i];
        }
        offsets[num_pairs] = current_offset;
        
        auto extract_time = std::chrono::high_resolution_clock::now();
        
        // Rest of processing is the same
        allocate_gpu_memory(num_pairs, all_hashes.size());
        
        cudaMemcpy(d_minimizer_hashes, all_hashes.data(), 
                   all_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_offsets, offsets.data(), 
                   offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_counts, counts.data(), 
                   counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        auto transfer_time = std::chrono::high_resolution_clock::now();
        
        database->lookup_batch_gpu(d_minimizer_hashes, d_minimizer_taxons, all_hashes.size());
        
        auto lookup_time = std::chrono::high_resolution_clock::now();
        
        int block_size = 256;
        int num_blocks = (num_pairs + block_size - 1) / block_size;
        
        classify_reads_kernel<<<num_blocks, block_size>>>(
            d_minimizer_taxons, d_minimizer_offsets, d_minimizer_counts,
            nullptr, d_read_classifications, d_confidence_scores, num_pairs
        );
        
        cudaDeviceSynchronize();
        auto classify_time = std::chrono::high_resolution_clock::now();
        
        classifications.resize(num_pairs);
        confidences.resize(num_pairs);
        
        cudaMemcpy(classifications.data(), d_read_classifications, 
                   num_pairs * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(confidences.data(), d_confidence_scores, 
                   num_pairs * sizeof(float), cudaMemcpyDeviceToHost);
        
        auto end = std::chrono::high_resolution_clock::now();
        
        auto extract_ms = std::chrono::duration_cast<std::chrono::microseconds>(extract_time - start).count();
        auto transfer_ms = std::chrono::duration_cast<std::chrono::microseconds>(transfer_time - extract_time).count();
        auto lookup_ms = std::chrono::duration_cast<std::chrono::microseconds>(lookup_time - transfer_time).count();
        auto classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(classify_time - lookup_time).count();
        auto return_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - classify_time).count();
        
        std::cout << "Paired batch timing (μs): extract=" << extract_ms 
                  << " transfer=" << transfer_ms
                  << " lookup=" << lookup_ms
                  << " classify=" << classify_ms
                  << " return=" << return_ms
                  << " total=" << extract_ms + transfer_ms + lookup_ms + classify_ms + return_ms
                  << " minimizers=" << all_hashes.size() << "\n";
    }
    
    void process_batch(const std::vector<std::string>& sequences,
                      std::vector<uint32_t>& classifications,
                      std::vector<float>& confidences) {
        size_t num_reads = sequences.size();
        if (num_reads == 0) return;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Step 1: Extract minimizers on GPU
        auto minimizer_results = extractor->extract_minimizers(sequences);
        
        // Prepare data for GPU processing
        std::vector<uint64_t> all_hashes;
        std::vector<uint32_t> offsets(num_reads + 1);
        std::vector<uint32_t> counts(num_reads);
        
        uint32_t current_offset = 0;
        for (size_t i = 0; i < num_reads; i++) {
            offsets[i] = current_offset;
            counts[i] = minimizer_results[i].size();
            
            for (const auto& m : minimizer_results[i]) {
                all_hashes.push_back(m.hash);  // Use the MurmurHash3 hash directly
            }
            
            current_offset += counts[i];
        }
        offsets[num_reads] = current_offset;
        
        auto extract_time = std::chrono::high_resolution_clock::now();
        
        // Allocate GPU memory
        allocate_gpu_memory(num_reads, all_hashes.size());
        
        // Transfer to GPU
        cudaMemcpy(d_minimizer_hashes, all_hashes.data(), 
                   all_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_offsets, offsets.data(), 
                   offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_counts, counts.data(), 
                   counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        auto transfer_time = std::chrono::high_resolution_clock::now();
        
        // Step 2: GPU database lookup
        database->lookup_batch_gpu(d_minimizer_hashes, d_minimizer_taxons, all_hashes.size());
        
        auto lookup_time = std::chrono::high_resolution_clock::now();
        
        // Step 3: GPU classification
        int block_size = 256;
        int num_blocks = (num_reads + block_size - 1) / block_size;
        
        classify_reads_kernel<<<num_blocks, block_size>>>(
            d_minimizer_taxons, d_minimizer_offsets, d_minimizer_counts,
            nullptr,  // Parent map not used in simple classification
            d_read_classifications, d_confidence_scores, num_reads
        );
        
        cudaDeviceSynchronize();
        auto classify_time = std::chrono::high_resolution_clock::now();
        
        // Transfer results back
        classifications.resize(num_reads);
        confidences.resize(num_reads);
        
        cudaMemcpy(classifications.data(), d_read_classifications, 
                   num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(confidences.data(), d_confidence_scores, 
                   num_reads * sizeof(float), cudaMemcpyDeviceToHost);
        
        auto end = std::chrono::high_resolution_clock::now();
        
        // Print timing breakdown
        auto extract_ms = std::chrono::duration_cast<std::chrono::microseconds>(extract_time - start).count();
        auto transfer_ms = std::chrono::duration_cast<std::chrono::microseconds>(transfer_time - extract_time).count();
        auto lookup_ms = std::chrono::duration_cast<std::chrono::microseconds>(lookup_time - transfer_time).count();
        auto classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(classify_time - lookup_time).count();
        auto return_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - classify_time).count();
        
        std::cout << "Batch timing (μs): extract=" << extract_ms 
                  << " transfer=" << transfer_ms
                  << " lookup=" << lookup_ms
                  << " classify=" << classify_ms
                  << " return=" << return_ms
                  << " total=" << extract_ms + transfer_ms + lookup_ms + classify_ms + return_ms
                  << " minimizers=" << all_hashes.size() << "\n";
    }
};

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <gpu_database.bin> <fastq_file> [fastq_file2]\n";
        std::cerr << "  Single-end mode: provide one FASTQ file\n";
        std::cerr << "  Paired-end mode: provide two FASTQ files (R1 and R2)\n";
        return 1;
    }
    
    std::string db_file = argv[1];
    std::string fastq_file1 = argv[2];
    std::string fastq_file2 = (argc == 4) ? argv[3] : "";
    bool paired_end = !fastq_file2.empty();
    
    try {
        // Initialize pipeline
        GPUProfilerPipeline pipeline(db_file);
        
        if (paired_end) {
            std::cout << "Running in paired-end mode\n";
            
            // Process paired files
            FastqReader reader1(fastq_file1);
            FastqReader reader2(fastq_file2);
            ReadBatch batch1, batch2;
            
            size_t total_pairs = 0;
            size_t classified_pairs = 0;
            
            auto start_time = std::chrono::high_resolution_clock::now();
            
            while (reader1.read_batch(batch1, 10000) && reader2.read_batch(batch2, 10000)) {
                if (batch1.size() != batch2.size()) {
                    std::cerr << "Warning: Batch size mismatch between paired files\n";
                    break;
                }
                
                std::vector<uint32_t> classifications;
                std::vector<float> confidences;
                
                pipeline.process_paired_batch(batch1.sequences, batch2.sequences, 
                                            classifications, confidences);
                
                // Count classified pairs
                for (size_t i = 0; i < classifications.size(); i++) {
                    if (classifications[i] > 0) {
                        classified_pairs++;
                    }
                }
                
                total_pairs += batch1.size();
                
                if (total_pairs % 100000 == 0) {
                    std::cout << "Processed " << total_pairs << " read pairs, "
                              << classified_pairs << " classified ("
                              << (100.0 * classified_pairs / total_pairs) << "%)\n";
                }
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            std::cout << "\n=== Final Results (Paired-End) ===\n";
            std::cout << "Total read pairs: " << total_pairs << "\n";
            std::cout << "Classified pairs: " << classified_pairs << " ("
                      << (100.0 * classified_pairs / total_pairs) << "%)\n";
            std::cout << "Processing time: " << duration.count() / 1000.0 << " seconds\n";
            std::cout << "Pairs per second: " << (total_pairs * 1000.0 / duration.count()) << "\n";
            
        } else {
            std::cout << "Running in single-end mode\n";
            
            // Process single file
            FastqReader reader(fastq_file1);
            ReadBatch batch;
            
            size_t total_reads = 0;
            size_t classified_reads = 0;
            
            auto start_time = std::chrono::high_resolution_clock::now();
            
            while (reader.read_batch(batch, 10000)) {
                std::vector<uint32_t> classifications;
                std::vector<float> confidences;
                
                pipeline.process_batch(batch.sequences, classifications, confidences);
                
                // Count classified reads
                for (size_t i = 0; i < classifications.size(); i++) {
                    if (classifications[i] > 0) {
                        classified_reads++;
                    }
                }
                
                total_reads += batch.size();
                
                if (total_reads % 100000 == 0) {
                    std::cout << "Processed " << total_reads << " reads, "
                              << classified_reads << " classified ("
                              << (100.0 * classified_reads / total_reads) << "%)\n";
                }
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            std::cout << "\n=== Final Results (Single-End) ===\n";
            std::cout << "Total reads: " << total_reads << "\n";
            std::cout << "Classified reads: " << classified_reads << " ("
                      << (100.0 * classified_reads / total_reads) << "%)\n";
            std::cout << "Processing time: " << duration.count() / 1000.0 << " seconds\n";
            std::cout << "Reads per second: " << (total_reads * 1000.0 / duration.count()) << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}