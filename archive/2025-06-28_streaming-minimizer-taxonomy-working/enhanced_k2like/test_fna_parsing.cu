// test_fna_parsing.cu
// Test program to verify FNA header parsing and sequence processing

#include <iostream>
#include <fstream>
#include <filesystem>
#include "processing/genome_file_processor.h"
#include "gpu_kraken_types.h"

int main() {
    std::cout << "=== Testing FNA Parsing ===" << std::endl;
    
    // Create test FNA file
    std::string test_fna = "test_concatenated.fna";
    std::ofstream outfile(test_fna);
    
    // Write test data with various whitespace scenarios
    outfile << ">kraken:taxid|455631|NZ_CM000441.1 Clostridioides difficile QCD-66c26 chromosome, whole genome shotgun sequence\n";
    outfile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    outfile << "CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n";
    outfile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n";
    outfile << "\n"; // Empty line
    outfile << ">kraken:taxid|562|NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome\n";
    outfile << "TTTAAAGGGCCCGGGATATATATATCGCGCGCGCGCGATATATATATGGGGCCCTTTAAA\n";
    outfile << "   AAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGG   \n"; // Leading/trailing spaces
    outfile << "\t\tTABBED\tSEQUENCE\tWITH\tTABS\t\t\n"; // Tabs
    outfile << ">kraken:taxid|1280|NZ_CP000253.1 Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete genome\n";
    outfile << "GGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT\r\n"; // Carriage return
    outfile << "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG\n";
    outfile.close();
    
    // Configure processing
    FileProcessingConfig config;
    config.progress_reporting = true;
    config.max_file_count = 100;
    
    // Process the file
    std::cout << "\nProcessing concatenated FNA file..." << std::endl;
    ConcatenatedFnaProcessor processor(test_fna, config);
    
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    
    bool success = processor.process_fna_file(genome_files, genome_taxon_ids, taxon_names, "./temp_genomes");
    
    if (success) {
        std::cout << "\n✓ Processing completed successfully!" << std::endl;
        std::cout << "  Genomes found: " << genome_files.size() << std::endl;
        std::cout << "  Taxon IDs extracted: ";
        for (size_t i = 0; i < genome_taxon_ids.size(); i++) {
            std::cout << genome_taxon_ids[i];
            if (i < genome_taxon_ids.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        
        std::cout << "\n  Taxon names:" << std::endl;
        for (const auto& [taxid, name] : taxon_names) {
            std::cout << "    " << taxid << " -> " << name << std::endl;
        }
        
        // Check the content of created files
        std::cout << "\n  Checking sequence integrity..." << std::endl;
        for (size_t i = 0; i < genome_files.size() && i < 3; i++) {
            std::ifstream check_file(genome_files[i]);
            std::string header, sequence;
            std::getline(check_file, header);
            
            // Read all sequence lines
            std::string line;
            while (std::getline(check_file, line)) {
                sequence += line;
            }
            check_file.close();
            
            std::cout << "    File " << i << ": " << genome_files[i] << std::endl;
            std::cout << "      Header: " << header.substr(0, 60) << "..." << std::endl;
            std::cout << "      Sequence length: " << sequence.length() << " bases" << std::endl;
            std::cout << "      First 60 bases: " << sequence.substr(0, 60) << std::endl;
            
            // Check for whitespace in sequence
            bool has_whitespace = false;
            for (char c : sequence) {
                if (std::isspace(c) || c == '\r' || c == '\n' || c == '\t') {
                    has_whitespace = true;
                    break;
                }
            }
            std::cout << "      Contains whitespace: " << (has_whitespace ? "YES (ERROR!)" : "NO (good)") << std::endl;
        }
        
    } else {
        std::cerr << "✗ Processing failed!" << std::endl;
    }
    
    // Cleanup
    std::filesystem::remove(test_fna);
    try {
        std::filesystem::remove_all("./temp_genomes");
    } catch (...) {}
    
    return success ? 0 : 1;
}