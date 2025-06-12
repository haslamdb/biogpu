#include <iostream>
#include <string>
#include <cstdint>

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

// Compute reverse complement of encoded k-mer
uint64_t reverse_complement(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        uint64_t base = kmer & 3;
        rc = (rc << 2) | (3 - base);  // Complement: A<->T (0<->3), C<->G (1<->2)
        kmer >>= 2;
    }
    return rc;
}

// Get canonical k-mer (minimum of forward and reverse complement)
uint64_t get_canonical_kmer(const std::string& kmer_str) {
    uint64_t forward = encode_kmer(kmer_str);
    if (forward == 0) return 0;  // Invalid k-mer
    
    uint64_t reverse = reverse_complement(forward, kmer_str.length());
    return std::min(forward, reverse);
}

// Decode k-mer for verification
std::string decode_kmer(uint64_t kmer, int k) {
    std::string result;
    for (int i = k-1; i >= 0; i--) {
        uint64_t base = (kmer >> (2*i)) & 3;
        result += "ACGT"[base];
    }
    return result;
}

int main() {
    // Test some k-mers and their reverse complements
    std::string test_kmers[] = {
        "ACGTACGT",
        "TGCATGCA",  // RC of first
        "AAAAAAAA",
        "TTTTTTTT",  // RC of third
        "ACGTGCAT",
        "ATGCACGT"   // RC of fifth
    };
    
    std::cout << "Testing canonical k-mer computation:\n";
    std::cout << "====================================\n";
    
    for (const auto& kmer : test_kmers) {
        uint64_t forward = encode_kmer(kmer);
        uint64_t reverse = reverse_complement(forward, kmer.length());
        uint64_t canonical = get_canonical_kmer(kmer);
        
        std::cout << "K-mer: " << kmer << "\n";
        std::cout << "  Forward encoded:  " << forward << "\n";
        std::cout << "  Reverse comp:     " << reverse << "\n";
        std::cout << "  Canonical:        " << canonical << "\n";
        std::cout << "  Canonical is:     " << (canonical == forward ? "forward" : "reverse") << "\n";
        std::cout << "  Decoded canon:    " << decode_kmer(canonical, kmer.length()) << "\n";
        std::cout << "\n";
    }
    
    return 0;
}