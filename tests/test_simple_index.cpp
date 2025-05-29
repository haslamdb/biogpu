#include <iostream>
#include <fstream>
#include "../runtime/kernels/resistance/fq_mutation_detector.cuh"

// Test program for simplified index
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <simple_index.h5>" << std::endl;
        return 1;
    }
    
    std::cout << "=== Testing Simplified FQ Index ===" << std::endl;
    std::cout << "Index file: " << argv[1] << std::endl;
    
    // Check file exists
    std::ifstream test(argv[1]);
    if (!test.good()) {
        std::cerr << "ERROR: Cannot open index file" << std::endl;
        return 1;
    }
    test.close();
    
    // Create detector and load index
    FQMutationDetectorCUDA detector;
    detector.loadIndex(argv[1]);
    
    std::cout << "\nIndex loaded successfully!" << std::endl;
    std::cout << "Number of k-mers: " << detector.num_kmers << std::endl;
    std::cout << "Number of sequences: " << detector.num_sequences << std::endl;
    std::cout << "Total reference length: " << detector.total_ref_length << std::endl;
    
    // Test k-mer search
    if (detector.num_kmers > 0) {
        std::cout << "\nTesting k-mer search..." << std::endl;
        
        // Create a test k-mer
        const char* test_seq = "AAAAAACAACCATCT";  // 15 bases
        uint64_t test_kmer = 0;
        for (int i = 0; i < 15; i++) {
            uint8_t base = 0;
            switch(test_seq[i]) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
            }
            test_kmer = (test_kmer << 2) | base;
        }
        
        std::cout << "Test k-mer: " << test_seq << std::endl;
        std::cout << "Encoded as: " << test_kmer << std::endl;
        
        // TODO: Test actual k-mer search on GPU
    }
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    return 0;
}