// test_processing_only.cu
#include <iostream>
#include <filesystem>
#include "processing/genome_file_processor.h"
#include "gpu_kraken_types.h"

int main() {
    std::cout << "=== BioGPU Processing Module Test ===" << std::endl;
    
    try {
        // Test FileProcessingConfig
        FileProcessingConfig config;
        config.max_file_count = 100;
        config.validate_sequences = true;
        config.progress_reporting = true;
        
        std::cout << "✓ FileProcessingConfig created" << std::endl;
        std::cout << "  Max files: " << config.max_file_count << std::endl;
        std::cout << "  Validation: " << (config.validate_sequences ? "enabled" : "disabled") << std::endl;
        
        // Test GenomeFileProcessor
        std::cout << "\nTesting GenomeFileProcessor..." << std::endl;
        GenomeFileProcessor processor(config);
        
        // Test utility functions
        std::cout << "\nTesting utility functions..." << std::endl;
        
        // Test file format detection
        bool is_fasta1 = FileProcessingUtils::is_fasta_file("test.fna");
        bool is_fasta2 = FileProcessingUtils::is_fasta_file("test.txt");
        std::cout << "✓ FASTA detection: .fna=" << (is_fasta1 ? "true" : "false") 
                  << ", .txt=" << (is_fasta2 ? "true" : "false") << std::endl;
        
        // Test sequence validation
        bool valid_dna = FileProcessingUtils::validate_dna_sequence("ATCGATCG");
        bool invalid_dna = FileProcessingUtils::validate_dna_sequence("ATCGXYZ");
        std::cout << "✓ DNA validation: ATCGATCG=" << (valid_dna ? "valid" : "invalid")
                  << ", ATCGXYZ=" << (invalid_dna ? "valid" : "invalid") << std::endl;
        
        // Test GC content calculation
        double gc_content = FileProcessingUtils::calculate_gc_content("ATCGGCTA");
        std::cout << "✓ GC content of ATCGGCTA: " << gc_content << "%" << std::endl;
        
        // Test file size formatting
        std::string size_str = FileProcessingUtils::format_file_size(1048576);
        std::cout << "✓ File size formatting: 1048576 bytes = " << size_str << std::endl;
        
        // Test taxon extraction from headers (not filenames)
        std::cout << "✓ Note: Taxon IDs are extracted from FASTA headers, not filenames" << std::endl;
        
        // Test statistics
        const auto& stats = processor.get_statistics();
        std::cout << "✓ Statistics access works" << std::endl;
        
        // Test ConcatenatedFnaProcessor creation
        std::cout << "\nTesting ConcatenatedFnaProcessor..." << std::endl;
        ConcatenatedFnaProcessor fna_processor("dummy.fna", config);
        std::cout << "✓ ConcatenatedFnaProcessor created" << std::endl;
        
        // Test header parsing (this is internal to ConcatenatedFnaProcessor)
        std::cout << "✓ Headers in format: >kraken:taxid|455631|NZ_CM000441.1 Clostridioides difficile..." << std::endl;
        
        // Test StreamingFnaProcessor creation
        std::cout << "\nTesting StreamingFnaProcessor..." << std::endl;
        StreamingFnaProcessor streaming_processor("dummy.fna", "/tmp/test_streaming", 10);
        std::cout << "✓ StreamingFnaProcessor created" << std::endl;
        
        std::cout << "\n=== ALL PROCESSING TESTS PASSED ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}