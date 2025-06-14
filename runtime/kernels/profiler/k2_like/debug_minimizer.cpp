// debug_minimizer.cpp - Debug program to test minimizer extraction
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include "minimizer_extractor.h"
#include "fastq_processing.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fastq_file>\n";
        return 1;
    }
    
    try {
        std::string fastq_file = argv[1];
        
        // Test 1: Direct minimizer extraction on a few reads
        std::cout << "=== Test 1: Direct Minimizer Extraction ===\n";
        {
            biogpu::FastqReader reader(fastq_file);
            ReadBatch batch;
            
            // Read just one batch
            if (reader.read_batch(batch, 100)) {
                std::cout << "Read " << batch.size() << " sequences\n";
                std::cout << "First sequence length: " << batch.sequences[0].length() << "\n";
                std::cout << "First 5 sequences:\n";
                
                for (size_t i = 0; i < std::min(size_t(5), batch.size()); i++) {
                    std::cout << "  Seq " << i << ": " << batch.sequences[i].substr(0, 50) << "...\n";
                }
                
                // Extract minimizers
                MinimizerExtractor extractor(31, 15);
                auto minimizers = extractor.extract_minimizers(batch.sequences);
                
                std::cout << "\nMinimizer counts for first 10 sequences:\n";
                for (size_t i = 0; i < std::min(size_t(10), minimizers.size()); i++) {
                    std::cout << "  Seq " << i << ": " << minimizers[i].size() << " minimizers\n";
                }
                
                // Calculate total
                size_t total_minimizers = 0;
                for (const auto& seq_mins : minimizers) {
                    total_minimizers += seq_mins.size();
                }
                
                std::cout << "\nTotal minimizers: " << total_minimizers << "\n";
                std::cout << "Average per sequence: " << (double)total_minimizers / batch.size() << "\n";
            }
        }
        
        // Test 2: Pipeline processing with limited reads
        std::cout << "\n=== Test 2: Pipeline Processing (1000 reads) ===\n";
        {
            biogpu::GPUMinimizerPipeline pipeline(31, 15, 100, 1);  // smaller batch, 1 GPU thread
            biogpu::MinimizerStats stats;
            
            size_t processed_reads = 0;
            size_t max_reads = 1000;
            
            auto start = std::chrono::high_resolution_clock::now();
            
            biogpu::FastqReader reader(fastq_file);
            ReadBatch batch;
            
            while (reader.read_batch(batch, 100) && processed_reads < max_reads) {
                auto minimizers = pipeline.process_batch(batch);
                stats.process_batch(batch, minimizers);
                
                processed_reads += batch.size();
                std::cout << "Processed " << processed_reads << " reads...\n";
                
                // Debug first batch
                if (processed_reads <= 100) {
                    size_t total_in_batch = 0;
                    for (const auto& seq_mins : minimizers) {
                        total_in_batch += seq_mins.size();
                    }
                    std::cout << "  Batch minimizers: " << total_in_batch 
                              << " (avg: " << (double)total_in_batch/batch.size() << " per read)\n";
                }
            }
            
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            
            std::cout << "\nProcessing time: " << duration.count() << " ms\n";
            std::cout << "Total reads: " << stats.get_total_reads() << "\n";
            std::cout << "Total minimizers: " << stats.get_total_minimizers() << "\n";
            std::cout << "Average minimizers per read: " << stats.get_average_minimizers_per_read() << "\n";
            
            // Don't call print_summary() to avoid hanging
            std::cout << "\n(Skipping detailed summary to avoid hanging)\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}