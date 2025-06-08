// bloom_cascade.cpp
#include <iostream>
#include <vector>
#include <string>
#include <bitset>
#include <functional>
#include <cmath>
#include <algorithm>
#include <random>
#include <cassert>

// MurmurHash3 for better distribution
uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const uint32_t r1 = 15;
    const uint32_t r2 = 13;
    const uint32_t m = 5;
    const uint32_t n = 0xe6546b64;
    
    uint32_t hash = seed;
    
    const int nblocks = len / 4;
    const uint32_t* blocks = (const uint32_t*) key;
    
    for (int i = 0; i < nblocks; i++) {
        uint32_t k = blocks[i];
        k *= c1;
        k = (k << r1) | (k >> (32 - r1));
        k *= c2;
        
        hash ^= k;
        hash = ((hash << r2) | (hash >> (32 - r2))) * m + n;
    }
    
    const uint8_t* tail = (const uint8_t*) (key + nblocks * 4);
    uint32_t k1 = 0;
    
    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1;
                k1 = (k1 << r1) | (k1 >> (32 - r1));
                k1 *= c2;
                hash ^= k1;
    }
    
    hash ^= len;
    hash ^= (hash >> 16);
    hash *= 0x85ebca6b;
    hash ^= (hash >> 13);
    hash *= 0xc2b2ae35;
    hash ^= (hash >> 16);
    
    return hash;
}

// Single Bloom filter implementation
class BloomFilter {
private:
    std::vector<bool> bits;
    size_t size;
    size_t num_hash_functions;
    size_t num_elements;
    std::vector<uint32_t> seeds;
    
public:
    BloomFilter(size_t expected_elements, double false_positive_rate) {
        // Calculate optimal size and number of hash functions
        size = -expected_elements * std::log(false_positive_rate) / (std::log(2) * std::log(2));
        size = ((size + 63) / 64) * 64;  // Round to multiple of 64 for alignment
        
        num_hash_functions = (size / expected_elements) * std::log(2);
        num_elements = 0;
        
        bits.resize(size, false);
        
        // Generate random seeds for hash functions
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> dis;
        
        for (size_t i = 0; i < num_hash_functions; i++) {
            seeds.push_back(dis(gen));
        }
    }
    
    void add(const std::string& item) {
        add(reinterpret_cast<const uint8_t*>(item.c_str()), item.length());
    }
    
    void add(const uint8_t* data, size_t len) {
        for (size_t i = 0; i < num_hash_functions; i++) {
            uint32_t hash = murmur3_32(data, len, seeds[i]);
            bits[hash % size] = true;
        }
        num_elements++;
    }
    
    bool contains(const std::string& item) const {
        return contains(reinterpret_cast<const uint8_t*>(item.c_str()), item.length());
    }
    
    bool contains(const uint8_t* data, size_t len) const {
        for (size_t i = 0; i < num_hash_functions; i++) {
            uint32_t hash = murmur3_32(data, len, seeds[i]);
            if (!bits[hash % size]) {
                return false;
            }
        }
        return true;
    }
    
    double get_false_positive_rate() const {
        size_t bits_set = std::count(bits.begin(), bits.end(), true);
        double ratio = (double)bits_set / size;
        return std::pow(ratio, num_hash_functions);
    }
    
    size_t get_size_bytes() const {
        return bits.size() / 8;
    }
};

// Hierarchical Bloom filter cascade for resistance detection
class ResistanceBloomCascade {
private:
    // Level 1: Genus-level resistance patterns
    BloomFilter genus_resistance;
    
    // Level 2: Species-specific mutations
    BloomFilter species_mutations;
    
    // Level 3: Exact mutation positions
    BloomFilter exact_positions;
    
    // Level 4: Co-occurring mutation patterns
    BloomFilter co_occurring_patterns;
    
    // Auxiliary data for verification
    struct MutationSignature {
        std::string gene;
        int position;
        char wild_type;
        char mutant;
        std::string organism;
        
        std::string to_string() const {
            return organism + ":" + gene + ":" + std::to_string(position) + 
                   ":" + wild_type + "->" + mutant;
        }
    };
    
    std::vector<MutationSignature> known_mutations;
    
public:
    ResistanceBloomCascade() 
        : genus_resistance(10000, 0.01),    // 10k patterns, 1% FPR
          species_mutations(50000, 0.001),  // 50k mutations, 0.1% FPR
          exact_positions(100000, 0.0001),  // 100k positions, 0.01% FPR
          co_occurring_patterns(20000, 0.01) // 20k patterns, 1% FPR
    {
        build_fluoroquinolone_cascade();
    }
    
    void build_fluoroquinolone_cascade() {
        // Level 1: Add genus-level patterns
        // These are broad patterns that indicate possible resistance
        genus_resistance.add("Enterobacteriaceae:fluoroquinolone_resistance");
        genus_resistance.add("Pseudomonas:fluoroquinolone_resistance");
        genus_resistance.add("Staphylococcus:fluoroquinolone_resistance");
        
        // Level 2: Add species-specific mutation patterns
        add_mutation("E.coli", "gyrA", 83, 'S', 'L');
        add_mutation("E.coli", "gyrA", 87, 'D', 'N');
        add_mutation("E.coli", "parC", 80, 'S', 'I');
        add_mutation("E.coli", "parC", 84, 'E', 'V');
        
        add_mutation("K.pneumoniae", "gyrA", 83, 'S', 'F');
        add_mutation("K.pneumoniae", "gyrA", 87, 'D', 'G');
        add_mutation("K.pneumoniae", "parC", 80, 'S', 'I');
        
        add_mutation("P.aeruginosa", "gyrA", 83, 'T', 'I');
        add_mutation("P.aeruginosa", "gyrA", 87, 'D', 'N');
        add_mutation("P.aeruginosa", "parC", 80, 'S', 'L');
        
        // Level 4: Add co-occurring patterns
        co_occurring_patterns.add("gyrA:S83L+parC:S80I");
        co_occurring_patterns.add("gyrA:S83F+gyrA:D87N");
        co_occurring_patterns.add("gyrA:T83I+parC:S80L");
    }
    
    void add_mutation(const std::string& organism, const std::string& gene, 
                     int position, char wild_type, char mutant) {
        MutationSignature sig{gene, position, wild_type, mutant, organism};
        known_mutations.push_back(sig);
        
        // Add to species level
        species_mutations.add(organism + ":" + gene + "_mutations");
        
        // Add to exact position level
        exact_positions.add(sig.to_string());
    }
    
    // Cascade screening result
    struct ScreeningResult {
        bool genus_hit;
        bool species_hit;
        bool exact_hit;
        bool co_occurring_hit;
        std::vector<std::string> potential_mutations;
        float confidence;
        
        bool has_resistance() const {
            return genus_hit || species_hit || exact_hit;
        }
    };
    
    // Screen a sequence through the cascade
    ScreeningResult screen_sequence(const std::string& sequence, 
                                   const std::string& organism_hint = "") {
        ScreeningResult result;
        result.genus_hit = false;
        result.species_hit = false;
        result.exact_hit = false;
        result.co_occurring_hit = false;
        result.confidence = 0.0;
        
        // Extract k-mers from sequence for checking
        const int kmer_size = 21;
        std::vector<std::string> kmers;
        
        for (size_t i = 0; i <= sequence.length() - kmer_size; i++) {
            kmers.push_back(sequence.substr(i, kmer_size));
        }
        
        // Level 1: Check genus patterns
        if (!organism_hint.empty()) {
            std::string genus = extract_genus(organism_hint);
            if (genus_resistance.contains(genus + ":fluoroquinolone_resistance")) {
                result.genus_hit = true;
                result.confidence += 0.1;
            }
        }
        
        // Level 2: Check species-specific patterns
        if (result.genus_hit && !organism_hint.empty()) {
            for (const auto& gene : {"gyrA", "gyrB", "parC", "parE"}) {
                std::string pattern = organism_hint + ":" + gene + "_mutations";
                if (species_mutations.contains(pattern)) {
                    result.species_hit = true;
                    result.confidence += 0.2;
                    break;
                }
            }
        }
        
        // Level 3: Check exact mutations (simplified - would use actual sequence alignment)
        if (result.species_hit) {
            for (const auto& mutation : known_mutations) {
                if (!organism_hint.empty() && mutation.organism != organism_hint) {
                    continue;
                }
                
                if (exact_positions.contains(mutation.to_string())) {
                    // In real implementation, would check if sequence contains mutation
                    result.exact_hit = true;
                    result.potential_mutations.push_back(mutation.to_string());
                    result.confidence += 0.3;
                }
            }
        }
        
        // Level 4: Check co-occurring patterns
        if (result.potential_mutations.size() >= 2) {
            // Check if any pairs of mutations are known to co-occur
            result.co_occurring_hit = true;
            result.confidence += 0.2;
        }
        
        // Normalize confidence
        result.confidence = std::min(result.confidence, 1.0f);
        
        return result;
    }
    
    // Batch screening for GPU simulation
    std::vector<ScreeningResult> batch_screen(const std::vector<std::string>& sequences,
                                            const std::vector<std::string>& organisms) {
        std::vector<ScreeningResult> results;
        results.reserve(sequences.size());
        
        #pragma omp parallel for
        for (size_t i = 0; i < sequences.size(); i++) {
            std::string org = (i < organisms.size()) ? organisms[i] : "";
            results.push_back(screen_sequence(sequences[i], org));
        }
        
        return results;
    }
    
    // Performance metrics
    void print_stats() const {
        std::cout << "Bloom Cascade Statistics:\n";
        std::cout << "Level 1 (Genus): " << genus_resistance.get_size_bytes() << " bytes, "
                  << "FPR: " << genus_resistance.get_false_positive_rate() << "\n";
        std::cout << "Level 2 (Species): " << species_mutations.get_size_bytes() << " bytes, "
                  << "FPR: " << species_mutations.get_false_positive_rate() << "\n";
        std::cout << "Level 3 (Exact): " << exact_positions.get_size_bytes() << " bytes, "
                  << "FPR: " << exact_positions.get_false_positive_rate() << "\n";
        std::cout << "Level 4 (Co-occur): " << co_occurring_patterns.get_size_bytes() << " bytes, "
                  << "FPR: " << co_occurring_patterns.get_false_positive_rate() << "\n";
        
        size_t total_size = genus_resistance.get_size_bytes() + 
                           species_mutations.get_size_bytes() +
                           exact_positions.get_size_bytes() +
                           co_occurring_patterns.get_size_bytes();
        std::cout << "Total size: " << total_size << " bytes (" 
                  << total_size / 1024.0 << " KB)\n";
    }
    
private:
    std::string extract_genus(const std::string& organism) {
        // Simple extraction - would be more sophisticated in practice
        if (organism.find("E.coli") != std::string::npos ||
            organism.find("K.pneumoniae") != std::string::npos) {
            return "Enterobacteriaceae";
        } else if (organism.find("P.aeruginosa") != std::string::npos) {
            return "Pseudomonas";
        } else if (organism.find("S.aureus") != std::string::npos) {
            return "Staphylococcus";
        }
        return "";
    }
};

// Example usage
#ifdef TEST
int main() {
    ResistanceBloomCascade cascade;
    
    // Print statistics
    cascade.print_stats();
    std::cout << "\n";
    
    // Test sequences (simplified - would be actual DNA sequences)
    std::vector<std::string> test_sequences = {
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",  // Mock sequence
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",  // Mock sequence
        "TACGTACGTACGTACGTACGTACGTACGTACGTACG"   // Mock sequence
    };
    
    std::vector<std::string> organisms = {
        "E.coli",
        "K.pneumoniae",
        "P.aeruginosa"
    };
    
    // Screen sequences
    auto results = cascade.batch_screen(test_sequences, organisms);
    
    // Print results
    for (size_t i = 0; i < results.size(); i++) {
        const auto& result = results[i];
        std::cout << "Sequence " << i << " (" << organisms[i] << "):\n";
        std::cout << "  Genus hit: " << (result.genus_hit ? "YES" : "NO") << "\n";
        std::cout << "  Species hit: " << (result.species_hit ? "YES" : "NO") << "\n";
        std::cout << "  Exact hit: " << (result.exact_hit ? "YES" : "NO") << "\n";
        std::cout << "  Co-occurring: " << (result.co_occurring_hit ? "YES" : "NO") << "\n";
        std::cout << "  Confidence: " << result.confidence << "\n";
        
        if (!result.potential_mutations.empty()) {
            std::cout << "  Potential mutations:\n";
            for (const auto& mut : result.potential_mutations) {
                std::cout << "    - " << mut << "\n";
            }
        }
        
        std::cout << "  Resistance detected: " 
                  << (result.has_resistance() ? "YES" : "NO") << "\n\n";
    }
    
    return 0;
}
#endif // TEST