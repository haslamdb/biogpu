#include <iostream>
#include <string>
#include <cstdint>
#include <fstream>

// MurmurHash3 finalizer (matching GPU implementation)
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
            default: return 0;
        }
    }
    return encoded;
}

// Compute reverse complement of encoded k-mer
uint64_t reverse_complement(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        uint64_t base = kmer & 3;
        rc = (rc << 2) | (3 - base);
        kmer >>= 2;
    }
    return rc;
}

// Get canonical k-mer
uint64_t get_canonical_kmer(const std::string& kmer_str) {
    uint64_t forward = encode_kmer(kmer_str);
    if (forward == 0) return 0;
    uint64_t reverse = reverse_complement(forward, kmer_str.length());
    return std::min(forward, reverse);
}

int main() {
    // Test with first few k-mers from the database
    std::ifstream file("data/test_kmers.txt");
    std::string line;
    int count = 0;
    
    std::cout << "Testing k-mer hash computation:\n";
    std::cout << "================================\n";
    
    while (std::getline(file, line) && count < 10) {
        if (line.empty() || line[0] == '#') continue;
        
        // Parse k-mer and taxon
        size_t tab_pos = line.find('\t');
        if (tab_pos == std::string::npos) continue;
        
        std::string kmer = line.substr(0, tab_pos);
        std::string taxon = line.substr(tab_pos + 1);
        
        // Compute hashes
        uint64_t forward = encode_kmer(kmer);
        uint64_t reverse = reverse_complement(forward, kmer.length());
        uint64_t canonical = get_canonical_kmer(kmer);
        uint64_t forward_hash = murmur_hash3_finalizer(forward);
        uint64_t canonical_hash = murmur_hash3_finalizer(canonical);
        
        std::cout << "K-mer " << count << ": " << kmer << " (taxon " << taxon << ")\n";
        std::cout << "  Forward:     " << forward << " -> hash: " << forward_hash << "\n";
        std::cout << "  Reverse:     " << reverse << "\n";
        std::cout << "  Canonical:   " << canonical << " -> hash: " << canonical_hash << "\n";
        std::cout << "  Uses:        " << (canonical == forward ? "forward" : "reverse") << "\n\n";
        
        count++;
    }
    
    return 0;
}