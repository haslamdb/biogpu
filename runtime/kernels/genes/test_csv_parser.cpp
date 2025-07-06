#include <iostream>
#include "sample_csv_parser.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv_file>" << std::endl;
        return 1;
    }
    
    BioGPU::SampleCSVParser parser;
    
    std::cout << "Testing CSV parser with: " << argv[1] << std::endl;
    
    if (!parser.parseFile(argv[1], false)) { // false = don't validate paths yet
        std::cerr << "Failed to parse CSV file" << std::endl;
        return 1;
    }
    
    std::cout << "\nParsed " << parser.getSampleCount() << " samples:" << std::endl;
    parser.printDetailed();
    
    // Test paired-end detection
    for (size_t i = 0; i < parser.getSampleCount(); i++) {
        const BioGPU::SampleInfo* sample = parser.getSample(i);
        if (sample) {
            std::cout << "\nSample " << i+1 << ": " << sample->sample_name << std::endl;
            std::cout << "  Paired-end: " << (sample->isPairedEnd() ? "Yes" : "No") << std::endl;
            std::cout << "  R1: " << sample->read1_path << std::endl;
            if (sample->isPairedEnd()) {
                std::cout << "  R2: " << sample->read2_path << std::endl;
            }
        }
    }
    
    return 0;
}