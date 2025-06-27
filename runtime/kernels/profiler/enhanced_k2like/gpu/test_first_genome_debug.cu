#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "../gpu_kraken_types.h"
#include "gpu/gpu_minimizer_extraction.cuh"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fna_file>" << std::endl;
        return 1;
    }

    std::ifstream file(argv[1]);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << argv[1] << std::endl;
        return 1;
    }

    // Read first genome only
    std::string header, sequence;
    std::getline(file, header);
    std::string line;
    while (std::getline(file, line) && line[0] != '>') {
        sequence += line;
    }
    file.close();

    std::cout << "First genome header: " << header << std::endl;
    std::cout << "Sequence length: " << sequence.length() << std::endl;
    std::cout << "First 100 chars: " << sequence.substr(0, 100) << std::endl;

    // Check for invalid characters
    size_t invalid_count = 0;
    for (char c : sequence) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && 
            c != 'a' && c != 'c' && c != 'g' && c != 't' &&
            c != 'N' && c != 'n') {
            invalid_count++;
            if (invalid_count < 10) {
                std::cout << "Invalid character at position " << (&c - sequence.data()) 
                         << ": '" << c << "' (ASCII " << (int)c << ")" << std::endl;
            }
        }
    }
    std::cout << "Total invalid characters: " << invalid_count << std::endl;

    // Test minimizer extraction on CPU first
    MinimizerParams params;
    params.k = 31;
    params.ell = 31;
    params.spaces = 7;
    params.xor_mask = 0;

    std::cout << "\nTesting CPU extraction of first few k-mers..." << std::endl;
    for (int i = 0; i < 5 && i + params.k <= sequence.length(); i++) {
        std::string kmer = sequence.substr(i, params.k);
        std::cout << "K-mer " << i << ": " << kmer << std::endl;
    }

    return 0;
}