// Test program for StreamingFnaProcessor
#include "processing/genome_file_processor.h"
#include <iostream>
#include <vector>
#include <filesystem>
#include <unistd.h>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fna_file>" << std::endl;
        return 1;
    }
    
    std::string fna_file = argv[1];
    std::string temp_dir = "/tmp/streaming_test_" + std::to_string(getpid());
    
    // Create temporary directory
    std::filesystem::create_directories(temp_dir);
    
    std::cout << "Testing StreamingFnaProcessor with: " << fna_file << std::endl;
    std::cout << "Temporary directory: " << temp_dir << std::endl;
    std::cout << "==================================================" << std::endl;
    
    // Create processor with batch size of 10
    StreamingFnaProcessor processor(fna_file, temp_dir, 10);
    
    std::vector<std::string> batch_files;
    std::vector<uint32_t> batch_taxons;
    
    int batch_num = 0;
    size_t total_genomes = 0;
    
    // Process all batches
    while (processor.process_next_batch(batch_files, batch_taxons)) {
        batch_num++;
        std::cout << "\nBatch " << batch_num << ":" << std::endl;
        std::cout << "  Genomes in batch: " << batch_files.size() << std::endl;
        
        // Print details of each genome in the batch
        for (size_t i = 0; i < batch_files.size(); i++) {
            std::cout << "  [" << (i+1) << "] Taxon: " << batch_taxons[i] 
                      << ", File: " << std::filesystem::path(batch_files[i]).filename() 
                      << std::endl;
            
            // Check file size
            auto file_size = std::filesystem::file_size(batch_files[i]);
            std::cout << "      Size: " << (file_size / 1024) << " KB" << std::endl;
        }
        
        total_genomes += batch_files.size();
    }
    
    std::cout << "\n==================================================" << std::endl;
    std::cout << "Processing complete!" << std::endl;
    std::cout << "Total batches: " << batch_num << std::endl;
    std::cout << "Total genomes processed: " << total_genomes << std::endl;
    std::cout << "Total genomes (from processor): " << processor.get_total_genomes() << std::endl;
    std::cout << "Total bases processed: " << (processor.get_total_bases() / (1024*1024)) << " MB" << std::endl;
    std::cout << "Sequences too large: " << processor.get_sequences_too_large() << std::endl;
    std::cout << "Sequences with invalid taxon: " << processor.get_sequences_with_invalid_taxon() << std::endl;
    
    // Note about filtered sequences
    std::cout << "\nNote: Sequences shorter than 1000bp are filtered out" << std::endl;
    
    // Check if more data available
    if (processor.has_more_data()) {
        std::cout << "WARNING: Processor reports more data available!" << std::endl;
    }
    
    // Cleanup
    std::cout << "\nKeeping temp directory for inspection. Remove manually: " << temp_dir << std::endl;
    // std::filesystem::remove_all(temp_dir);
    
    return 0;
}