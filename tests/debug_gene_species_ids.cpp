#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cuda_runtime.h>
#include "fq_mutation_detector.cuh"

// Debug the k-mer index to understand gene/species ID assignments
void debugKmerIndex(const std::string& index_path) {
    std::cout << "=== Debugging K-mer Index Gene/Species IDs ===" << std::endl;
    
    // Load k-mer index
    std::string kmer_index_path = index_path + "/kmer_index.bin";
    std::ifstream kmer_file(kmer_index_path, std::ios::binary);
    
    if (!kmer_file.good()) {
        std::cerr << "ERROR: Cannot open k-mer index: " << kmer_index_path << std::endl;
        return;
    }
    
    // Read header
    uint32_t num_entries, kmer_length;
    kmer_file.read(reinterpret_cast<char*>(&num_entries), sizeof(uint32_t));
    kmer_file.read(reinterpret_cast<char*>(&kmer_length), sizeof(uint32_t));
    
    std::cout << "Index contains " << num_entries << " k-mer entries (k=" << kmer_length << ")" << std::endl;
    
    // Track unique gene and species IDs
    std::map<uint32_t, int> gene_id_counts;
    std::map<uint32_t, int> species_id_counts;
    std::map<std::pair<uint32_t, uint32_t>, int> gene_species_pairs;
    
    // Read and analyze all entries
    for (uint32_t i = 0; i < std::min(num_entries, 10000U); i++) { // Limit to first 10K for analysis
        KmerEntry entry;
        kmer_file.read(reinterpret_cast<char*>(&entry.kmer), sizeof(uint64_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.gene_id), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.species_id), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.seq_id), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&entry.position), sizeof(uint16_t));
        
        gene_id_counts[entry.gene_id]++;
        species_id_counts[entry.species_id]++;
        gene_species_pairs[{entry.gene_id, entry.species_id}]++;
        
        // Show first few entries
        if (i < 20) {
            std::cout << "Entry " << i << ": gene_id=" << entry.gene_id 
                     << ", species_id=" << entry.species_id 
                     << ", seq_id=" << entry.seq_id 
                     << ", kmer=" << entry.kmer << std::endl;
        }
    }
    kmer_file.close();
    
    // Analysis
    std::cout << "\n=== Gene ID Distribution ===" << std::endl;
    std::cout << "Unique gene IDs found: " << gene_id_counts.size() << std::endl;
    for (const auto& pair : gene_id_counts) {
        std::cout << "  Gene ID " << pair.first << ": " << pair.second << " k-mers" << std::endl;
    }
    
    std::cout << "\n=== Species ID Distribution ===" << std::endl;
    std::cout << "Unique species IDs found: " << species_id_counts.size() << std::endl;
    for (const auto& pair : species_id_counts) {
        std::cout << "  Species ID " << pair.first << ": " << pair.second << " k-mers" << std::endl;
    }
    
    std::cout << "\n=== Gene-Species Combinations ===" << std::endl;
    std::cout << "Unique gene-species pairs: " << gene_species_pairs.size() << std::endl;
    for (const auto& pair : gene_species_pairs) {
        std::cout << "  Gene " << pair.first.first << " + Species " << pair.first.second 
                 << ": " << pair.second << " k-mers" << std::endl;
    }
    
    // Check for the suspicious pattern
    if (gene_species_pairs.size() == 1 && gene_species_pairs.begin()->first.first == gene_species_pairs.begin()->first.second) {
        std::cout << "\n🚨 ISSUE DETECTED: All k-mers have gene_id == species_id == " 
                 << gene_species_pairs.begin()->first.first << std::endl;
        std::cout << "This suggests a bug in index creation where gene_id and species_id are being set to the same value." << std::endl;
    }
}

// Analyze what gene_id=6 and species_id=6 actually represent
void analyzeGeneSpeciesMapping() {
    std::cout << "\n=== Gene/Species ID Mapping Analysis ===" << std::endl;
    
    // Common FQ resistance gene mapping (if you have documentation)
    std::map<uint32_t, std::string> expected_genes = {
        {0, "gyrA"},
        {1, "gyrB"}, 
        {2, "parC"},
        {3, "parE"},
        {4, "qnrA"},
        {5, "qnrB"},
        {6, "qnrS"},  // This might be what you're seeing
        {7, "aac(6')-Ib-cr"}
    };
    
    std::map<uint32_t, std::string> expected_species = {
        {0, "E. coli"},
        {1, "K. pneumoniae"},
        {2, "P. aeruginosa"},
        {3, "S. aureus"},
        {4, "E. faecium"},   // VRE
        {5, "E. faecalis"},  // VRE
        {6, "C. difficile"}  // This might be what you're seeing
    };
    
    std::cout << "If gene_id=6 corresponds to: " << expected_genes[6] << std::endl;
    std::cout << "If species_id=6 corresponds to: " << expected_species[6] << std::endl;
    std::cout << "\nFor VRE12 sample, we'd expect:" << std::endl;
    std::cout << "  Species: E. faecium or E. faecalis (IDs 4 or 5)" << std::endl;
    std::cout << "  Genes: vanA, vanB resistance genes (not just FQ genes)" << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <index_directory>" << std::endl;
        std::cerr << "Example: " << argv[0] << " /data/fq_resistance_index" << std::endl;
        return 1;
    }
    
    std::string index_path = argv[1];
    
    debugKmerIndex(index_path);
    analyzeGeneSpeciesMapping();
    
    return 0;
}