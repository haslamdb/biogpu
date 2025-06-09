// test_simple_minimizer.cpp - Simple test to debug minimizer extraction
#include <iostream>
#include <vector>
#include <string>
#include "minimizer_extractor.h"

int main() {
    std::cout << "Creating MinimizerExtractor...\n";
    MinimizerExtractor extractor(31, 15);
    
    std::cout << "Creating test sequence...\n";
    std::vector<std::string> sequences = {
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    };
    
    std::cout << "Sequence length: " << sequences[0].length() << "\n";
    std::cout << "Extracting minimizers...\n";
    
    try {
        auto result = extractor.extract_minimizers(sequences);
        
        std::cout << "Extraction complete!\n";
        std::cout << "Number of sequences processed: " << result.size() << "\n";
        std::cout << "Minimizers found in first sequence: " << result[0].size() << "\n";
        
        if (result[0].size() > 0) {
            std::cout << "First few minimizers:\n";
            for (size_t i = 0; i < std::min(size_t(5), result[0].size()); i++) {
                std::cout << "  Hash: " << result[0][i].hash 
                          << ", Position: " << result[0][i].position 
                          << ", Reverse: " << (result[0][i].is_reverse ? "yes" : "no") << "\n";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}