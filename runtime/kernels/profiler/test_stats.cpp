// test_stats.cpp - Test statistics collection issue
#include <iostream>
#include <vector>
#include <chrono>
#include "fastq_processing.h"

int main() {
    // Test the stats collection directly
    biogpu::MinimizerStats stats;
    
    // Simulate processing 100 batches of 10000 reads each
    for (int batch_num = 0; batch_num < 100; batch_num++) {
        ReadBatch batch;
        batch.reserve(10000);
        
        // Create fake batch
        for (int i = 0; i < 10000; i++) {
            batch.sequences.push_back("ACGTACGTACGT");
            batch.headers.push_back("read_" + std::to_string(batch_num * 10000 + i));
            batch.qualities.push_back("IIIIIIIIIIII");
        }
        
        // Create fake minimizers (about 14 per read)
        std::vector<std::vector<Minimizer>> minimizers;
        minimizers.reserve(10000);
        
        for (int i = 0; i < 10000; i++) {
            std::vector<Minimizer> read_mins;
            for (int j = 0; j < 14; j++) {
                Minimizer m;
                m.hash = (batch_num * 10000 + i) % 5000;  // Create some repeating hashes
                m.position = j * 10;
                m.is_reverse = false;
                read_mins.push_back(m);
            }
            minimizers.push_back(read_mins);
        }
        
        // Process batch
        auto start = std::chrono::high_resolution_clock::now();
        stats.process_batch(batch, minimizers);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        if (batch_num % 10 == 0) {
            std::cout << "Processed batch " << batch_num << " in " << duration.count() << " us\n";
            std::cout << "  Total reads so far: " << stats.get_total_reads() << "\n";
            std::cout << "  Total minimizers so far: " << stats.get_total_minimizers() << "\n";
            std::cout << "  Average: " << stats.get_average_minimizers_per_read() << "\n";
        }
    }
    
    std::cout << "\nFinal statistics:\n";
    std::cout << "Total reads: " << stats.get_total_reads() << "\n";
    std::cout << "Total minimizers: " << stats.get_total_minimizers() << "\n";
    std::cout << "Average: " << stats.get_average_minimizers_per_read() << "\n";
    
    std::cout << "\nTesting print_summary (this might hang)...\n";
    auto start = std::chrono::high_resolution_clock::now();
    stats.print_summary();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "print_summary took " << duration.count() << " ms\n";
    
    return 0;
}