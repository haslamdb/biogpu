// build_test_database.cpp - Build a test microbial database for GPU profiler
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <cstring>
#include "gpu_kmer_database.h"
#include "minimizer_common.h"

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

// Convert sequence to 2-bit encoding
uint64_t encode_kmer(const std::string& kmer) {
    uint64_t encoded = 0;
    for (char c : kmer) {
        encoded <<= 2;
        switch (c) {
            case 'A': case 'a': encoded |= 0; break;
            case 'C': case 'c': encoded |= 1; break;
            case 'G': case 'g': encoded |= 2; break;
            case 'T': case 't': encoded |= 3; break;
        }
    }
    return encoded;
}

// Reverse complement
std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.length(), 'N');
    for (size_t i = 0; i < seq.length(); i++) {
        switch (seq[seq.length() - 1 - i]) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            case 'a': rc[i] = 't'; break;
            case 't': rc[i] = 'a'; break;
            case 'c': rc[i] = 'g'; break;
            case 'g': rc[i] = 'c'; break;
        }
    }
    return rc;
}

// Get canonical k-mer (lexicographically smaller of forward/reverse)
std::string get_canonical_kmer(const std::string& kmer) {
    std::string rc = reverse_complement(kmer);
    return (kmer < rc) ? kmer : rc;
}

// Simulated organism data
struct TestOrganism {
    std::string name;
    uint32_t taxon_id;
    std::string marker_sequence;  // A characteristic sequence for this organism
};

// Create test organisms with known marker sequences
std::vector<TestOrganism> create_test_organisms() {
    return {
        // Common gut bacteria
        {"Escherichia_coli", 100, 
         "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
         "CGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGT"
         "AACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTT"},
        
        {"Klebsiella_pneumoniae", 101,
         "ATGAACATTAAAGGTCTGGTTGTTGCCGCTGCTGCCTTAGGTGGTGGCGCAACTGTCGCAGCAGAT"
         "ATTGCTCAGGATAAACTGGAAGAGAAACTGAACGCTGCGTTGCAGGCCTTTGATAAAAACAAAGAT"
         "GAAACCCGTGAAGGTAAAACCGATGCTCGCGCTTATCAGGCCTACAAGTTCAGCGTCGTTTCCGTC"},
        
        {"Enterococcus_faecalis", 102,
         "ATGACAGAAGAAAAATTTCATCTTGAAGATATTTTGAAGAAGCTCAAAGCTGAAACGGATAAAAAA"
         "GTCAATGAAAATGCACAAGAAGAAATTGAACGTTTAACGCTTGAAAAAGGTGTTTTAGATATTTAT"
         "GGTCCAACTCAAACAATGACTCGTAATGCTGAAAAAGCAATTGCGGATGAAATTGTCAGTAAAATT"},
        
        {"Staphylococcus_aureus", 103,
         "ATGAATATCAAAGAGCAAATTAAAGAATTAAAAGCAGAAGATGAAAAAGATAAAAAGAAAGATAAT"
         "GAACAAGCTCAAAAGAAACTAGGTGTTCTAGATATTGAAGTGACAAATAATGTAAAAGATAGTAAT"
         "AATCAAGATATTGATATCGTTAAAGAAGCAGAAAAGAAACAAGATATTAAAGCAGAAAAAGATGAA"},
        
        // Common respiratory pathogens
        {"Streptococcus_pneumoniae", 104,
         "ATGATTCGAGAAGATTTAGAAATTTATTTGAGTATTTCAGAACAAGATGGAAAAATTGAAATTGCC"
         "CGTCTTTCAGTAGATACAGCAGCTAAAGAAATTGAACGTCTTTCCCTAAGAGAAGATTTAAGAGAT"
         "GTTGAATTGGCTAAAGAAAAAGCTAAAGAAGCCTTGAAAGAATTGGATAGCCTGTATGCTGAAATT"},
        
        {"Haemophilus_influenzae", 105,
         "ATGGCAAAACACAATATCGTAGCGAACCTGCAAACCAAACAAGGTGAATTGCTGGCAGTTTTAGGT"
         "GGTAGCGGTAGCGGTTCTGCAGAAACACAAGATGAAACTCAAGAAGAAATCGCTGCGCTGCGTAAA"
         "CAAGCGATTAAAGACGCTCGTGAAGCGAAAACCGACGCGCGTGCGTATCAAGCTTACAAATTCAGC"},
        
        // Common UTI pathogens
        {"Proteus_mirabilis", 106,
         "ATGAGTAATAAAAAGCTGATTTACAGTGCAGGGATGTCAGGATTTGGTCGTGAAGCATTAGGTAAA"
         "GGGATGACGATGATGCCGGTATTAGCCGAAGATATTGCTCAAGATAAACTCGAGGAAAAACTTAAC"
         "GCCGCATTACAAGCCTTTGACAAGAACAAAGATGAAACACGTGAAGGCAAAACCGATGCCCGTGCT"},
        
        {"Pseudomonas_aeruginosa", 107,
         "ATGAGCAAATCCTTGTCGCGAATCGCCTTGCTGGCCGTCGCCCTGGCGGCCAACCTGGTGAGCAAA"
         "GATATTGCCCAGGGCAAACTGGAAGAGAAGCTGAACGCCGCCCTGCAGGCCTTCGACAAGAACAAG"
         "GATGAGACCCGCGAGGGTAAGACCGATGCCCGCGCCTATCAGGCCTACAAGTTCAGCGTGGTGTCC"},
        
        // Clostridium difficile (your research interest!)
        {"Clostridioides_difficile", 108,
         "ATGAAAGTAAAAGAATTTAAAGCAGAAGAAGAAATTGATGCTGAAAAAGCTAAAGCAGCTGTTGAA"
         "GCTATGGCTGCTGAAGATGCTAAAGCTGATATTGATGCTATTGAAGCTGATGAAGCTAAAGAAGAA"
         "GCTGATGAAGCTGTTGAAGCTAAAGCTGATGAAGCTATTAAAGCTGATGAAGCTGAAGATGCTAAA"}
    };
}

// Generate synthetic reads from an organism's sequence
std::vector<std::string> generate_reads_from_sequence(const std::string& sequence, 
                                                      int num_reads, 
                                                      int read_length,
                                                      double error_rate = 0.01) {
    std::vector<std::string> reads;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> pos_dist(0, sequence.length() - read_length);
    std::uniform_real_distribution<> error_dist(0.0, 1.0);
    std::uniform_int_distribution<> base_dist(0, 3);
    const char* bases = "ACGT";
    
    for (int i = 0; i < num_reads; i++) {
        int start_pos = pos_dist(gen);
        std::string read = sequence.substr(start_pos, read_length);
        
        // Add sequencing errors
        for (char& c : read) {
            if (error_dist(gen) < error_rate) {
                c = bases[base_dist(gen)];
            }
        }
        
        reads.push_back(read);
    }
    
    return reads;
}

// Extract all k-mers from a sequence
std::vector<std::pair<uint64_t, std::string>> extract_kmers(const std::string& sequence, int k) {
    std::vector<std::pair<uint64_t, std::string>> kmers;
    
    if (sequence.length() < k) return kmers;
    
    for (size_t i = 0; i <= sequence.length() - k; i++) {
        std::string kmer = sequence.substr(i, k);
        
        // Skip k-mers with N or other ambiguous bases
        bool valid = true;
        for (char c : kmer) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                valid = false;
                break;
            }
        }
        
        if (valid) {
            std::string canonical = get_canonical_kmer(kmer);
            uint64_t encoded = encode_kmer(canonical);
            uint64_t hash = murmur_hash3_finalizer(encoded);
            kmers.push_back({hash, canonical});
        }
    }
    
    return kmers;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <output_database.bin>\n";
        std::cerr << "\nThis tool creates a test microbial database with known organisms.\n";
        std::cerr << "The database can be used to test the GPU profiler pipeline.\n";
        return 1;
    }
    
    std::string output_file = argv[1];
    const int k = 31;  // k-mer size
    
    std::cout << "Building test microbial database...\n";
    std::cout << "K-mer size: " << k << "\n\n";
    
    // Create database
    GPUKmerDatabase database;
    
    // Get test organisms
    auto organisms = create_test_organisms();
    
    // Build k-mer database
    std::unordered_map<uint64_t, uint32_t> kmer_to_taxon;
    size_t total_kmers = 0;
    
    for (const auto& organism : organisms) {
        std::cout << "Processing " << organism.name 
                  << " (taxon " << organism.taxon_id << ")...\n";
        
        // Extract k-mers from marker sequence
        auto kmers = extract_kmers(organism.marker_sequence, k);
        
        // Also generate some synthetic genome sequence for more k-mers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> base_dist(0, 3);
        const char* bases = "ACGT";
        
        // Generate additional random sequence (simulating more of the genome)
        std::string additional_sequence;
        for (int i = 0; i < 5000; i++) {
            additional_sequence += bases[base_dist(gen)];
        }
        
        // Mix in some conserved regions from the marker
        for (int i = 0; i < 10; i++) {
            int insert_pos = (i * 500) % (additional_sequence.length() - organism.marker_sequence.length());
            additional_sequence.replace(insert_pos, 50, 
                                      organism.marker_sequence.substr(0, 50));
        }
        
        auto additional_kmers = extract_kmers(additional_sequence, k);
        
        // Add all k-mers to database
        size_t organism_kmers = 0;
        for (const auto& [hash, kmer] : kmers) {
            if (kmer_to_taxon.find(hash) == kmer_to_taxon.end()) {
                kmer_to_taxon[hash] = organism.taxon_id;
                database.add_kmer(hash, organism.taxon_id);
                organism_kmers++;
            }
        }
        
        // Add some of the additional k-mers (not all, to simulate incomplete database)
        for (size_t i = 0; i < additional_kmers.size(); i += 3) {  // Take every 3rd k-mer
            const auto& [hash, kmer] = additional_kmers[i];
            if (kmer_to_taxon.find(hash) == kmer_to_taxon.end()) {
                kmer_to_taxon[hash] = organism.taxon_id;
                database.add_kmer(hash, organism.taxon_id);
                organism_kmers++;
            }
        }
        
        std::cout << "  Added " << organism_kmers << " unique k-mers\n";
        total_kmers += organism_kmers;
    }
    
    std::cout << "\nTotal unique k-mers in database: " << total_kmers << "\n";
    
    // Finalize and save database
    std::cout << "\nFinalizing database and building GPU hash table...\n";
    database.finalize_and_upload();
    
    std::cout << "Saving database to " << output_file << "...\n";
    database.save_binary(output_file);
    
    std::cout << "\nDatabase creation complete!\n";
    std::cout << "Database size: " << database.get_database_size_bytes() / (1024.0 * 1024.0) << " MB\n";
    
    // Generate test FASTQ file
    std::cout << "\nGenerating test FASTQ file with mixed organisms...\n";
    std::string test_fastq = output_file + ".test.fastq";
    std::ofstream fastq_out(test_fastq);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> org_dist(0, organisms.size() - 1);
    
    int read_id = 0;
    for (int i = 0; i < organisms.size(); i++) {
        // Generate reads from each organism
        auto reads = generate_reads_from_sequence(organisms[i].marker_sequence, 100, 150);
        
        for (const auto& read : reads) {
            fastq_out << "@read_" << read_id++ << "_" << organisms[i].name << "\n";
            fastq_out << read << "\n";
            fastq_out << "+\n";
            fastq_out << std::string(read.length(), 'I') << "\n";  // Fake quality scores
        }
    }
    
    // Add some random mixed reads
    for (int i = 0; i < 500; i++) {
        int org_idx = org_dist(gen);
        auto reads = generate_reads_from_sequence(organisms[org_idx].marker_sequence, 1, 150);
        
        fastq_out << "@read_" << read_id++ << "_mixed\n";
        fastq_out << reads[0] << "\n";
        fastq_out << "+\n";
        fastq_out << std::string(reads[0].length(), 'I') << "\n";
    }
    
    fastq_out.close();
    std::cout << "Test FASTQ file created: " << test_fastq << "\n";
    std::cout << "Total reads: " << read_id << "\n";
    
    std::cout << "\nYou can now test the pipeline with:\n";
    std::cout << "./gpu_profiler_pipeline " << output_file << " " << test_fastq << "\n";
    
    return 0;
}