// build_db_from_kmers.cpp - Build GPU database from k-mer list
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../runtime/kernels/profiler/gpu_kmer_database.h"

using namespace biogpu;

// MurmurHash3 finalizer (matching the GPU implementation)
uint64_t murmur_hash3_finalizer(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccd;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53;
    key ^= key >> 33;
    return key;
}

// Convert k-mer string to 2-bit encoding
uint64_t encode_kmer(const std::string& kmer) {
    uint64_t encoded = 0;
    for (char c : kmer) {
        encoded <<= 2;
        switch (c) {
            case 'A': case 'a': encoded |= 0; break;
            case 'C': case 'c': encoded |= 1; break;
            case 'G': case 'g': encoded |= 2; break;
            case 'T': case 't': encoded |= 3; break;
            default: return 0; // Invalid k-mer
        }
    }
    return encoded;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <kmer_list.txt> <output_database.bin>\n";
        std::cerr << "\nExpected format of kmer_list.txt:\n";
        std::cerr << "# Comments start with #\n";
        std::cerr << "ACGTACGTACGTACGTACGTACGTACGTACG\t100\n";
        std::cerr << "TGCATGCATGCATGCATGCATGCATGCATGC\t101\n";
        std::cerr << "...\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
    std::cout << "Building GPU k-mer database from " << input_file << "...\n";
    
    // Create database
    GPUKmerDatabase database;
    
    // Read k-mer list
    std::ifstream infile(input_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open input file " << input_file << "\n";
        return 1;
    }
    
    std::string line;
    size_t line_num = 0;
    size_t kmers_added = 0;
    size_t max_taxon_id = 0;
    
    while (std::getline(infile, line)) {
        line_num++;
        
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // Parse line: kmer<tab>taxon_id
        std::istringstream iss(line);
        std::string kmer;
        uint32_t taxon_id;
        
        if (!(iss >> kmer >> taxon_id)) {
            std::cerr << "Warning: Invalid format at line " << line_num << "\n";
            continue;
        }
        
        // Encode k-mer and compute hash
        uint64_t encoded = encode_kmer(kmer);
        if (encoded == 0) {
            std::cerr << "Warning: Invalid k-mer at line " << line_num << ": " << kmer << "\n";
            continue;
        }
        
        uint64_t hash = murmur_hash3_finalizer(encoded);
        
        // Add to database
        database.add_kmer(hash, taxon_id);
        kmers_added++;
        max_taxon_id = std::max(max_taxon_id, (size_t)taxon_id);
        
        if (kmers_added % 100000 == 0) {
            std::cout << "  Processed " << kmers_added << " k-mers...\n";
        }
    }
    
    infile.close();
    
    std::cout << "\nDatabase statistics:\n";
    std::cout << "  Total k-mers: " << kmers_added << "\n";
    std::cout << "  Maximum taxon ID: " << max_taxon_id << "\n";
    
    // Build GPU hash table and save
    std::cout << "\nBuilding GPU hash table...\n";
    database.finalize_and_upload();
    
    std::cout << "Saving database to " << output_file << "...\n";
    database.save_binary(output_file);
    
    std::cout << "\nDatabase creation complete!\n";
    std::cout << "You can now use this database with:\n";
    std::cout << "  ./gpu_profiler_pipeline " << output_file << " <reads.fastq>\n";
    
    return 0;
}