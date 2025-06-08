// minimizer_classifier.cpp
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <cassert>

// Constants
const int KMER_SIZE = 21;
const int WINDOW_SIZE = 11;  // Window for minimizer selection
const uint64_t HASH_PRIME = 31;

// Data structures
struct MinimizerHit {
    uint64_t minimizer_hash;
    uint32_t genome_id;
    uint32_t position;
    uint16_t strand;  // 0 = forward, 1 = reverse
};

struct GenomeInfo {
    std::string name;
    std::string taxonomy;
    uint32_t genome_id;
    uint32_t length;
    std::vector<std::pair<uint32_t, uint32_t>> resistance_regions;  // Start, end positions
};

struct ClassificationResult {
    uint32_t read_id;
    uint32_t best_genome_id;
    float confidence;
    std::vector<uint32_t> multi_map_genomes;  // For reads mapping to multiple genomes
    bool has_resistance_region;
};

class MinimizerIndex {
private:
    std::unordered_multimap<uint64_t, MinimizerHit> index;
    std::vector<GenomeInfo> genome_info;
    
    // Hash function for k-mers
    uint64_t hash_kmer(const std::string& kmer) {
        uint64_t hash = 0;
        for (char c : kmer) {
            hash = hash * HASH_PRIME + encode_base(c);
        }
        return hash;
    }
    
    uint8_t encode_base(char base) {
        switch(base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0;  // Treat N as A for now
        }
    }
    
    std::string reverse_complement(const std::string& seq) {
        std::string rc(seq.length(), 'N');
        for (size_t i = 0; i < seq.length(); i++) {
            char c = seq[seq.length() - 1 - i];
            switch(c) {
                case 'A': rc[i] = 'T'; break;
                case 'T': rc[i] = 'A'; break;
                case 'C': rc[i] = 'G'; break;
                case 'G': rc[i] = 'C'; break;
                default: rc[i] = 'N';
            }
        }
        return rc;
    }

public:
    // Build index from genome sequences
    void build_index(const std::vector<std::pair<std::string, std::string>>& genomes) {
        uint32_t genome_id = 0;
        
        for (const auto& [name, sequence] : genomes) {
            GenomeInfo info;
            info.name = name;
            info.genome_id = genome_id;
            info.length = sequence.length();
            
            // Mark resistance regions (for now, mock data - would load from database)
            if (name.find("E.coli") != std::string::npos) {
                // Mock gyrA region
                info.resistance_regions.push_back({2337000, 2339000});
                // Mock parC region  
                info.resistance_regions.push_back({3176000, 3178000});
            }
            
            genome_info.push_back(info);
            
            // Extract minimizers using sliding window
            for (size_t i = 0; i <= sequence.length() - KMER_SIZE; i++) {
                std::string kmer = sequence.substr(i, KMER_SIZE);
                if (kmer.find('N') != std::string::npos) continue;  // Skip kmers with N
                
                // Get canonical k-mer (lexicographically smaller of forward/reverse)
                std::string rc_kmer = reverse_complement(kmer);
                bool use_forward = kmer < rc_kmer;
                uint64_t hash = hash_kmer(use_forward ? kmer : rc_kmer);
                
                // Check if this is a minimizer in its window
                bool is_minimizer = true;
                size_t window_start = (i >= WINDOW_SIZE/2) ? i - WINDOW_SIZE/2 : 0;
                size_t window_end = std::min(i + WINDOW_SIZE/2, sequence.length() - KMER_SIZE);
                
                for (size_t j = window_start; j <= window_end; j++) {
                    if (j == i) continue;
                    std::string window_kmer = sequence.substr(j, KMER_SIZE);
                    if (window_kmer.find('N') != std::string::npos) continue;
                    std::string window_rc = reverse_complement(window_kmer);
                    uint64_t window_hash = hash_kmer(window_kmer < window_rc ? window_kmer : window_rc);
                    if (window_hash < hash) {
                        is_minimizer = false;
                        break;
                    }
                }
                
                if (is_minimizer) {
                    MinimizerHit hit;
                    hit.minimizer_hash = hash;
                    hit.genome_id = genome_id;
                    hit.position = i;
                    hit.strand = use_forward ? 0 : 1;
                    index.insert({hash, hit});
                }
            }
            
            genome_id++;
        }
        
        std::cout << "Built index with " << index.size() << " minimizers from " 
                  << genomes.size() << " genomes\n";
    }
    
    // Classify a single read
    ClassificationResult classify_read(const std::string& read, uint32_t read_id) {
        ClassificationResult result;
        result.read_id = read_id;
        result.has_resistance_region = false;
        
        // Count hits per genome
        std::unordered_map<uint32_t, uint32_t> genome_hits;
        std::unordered_map<uint32_t, std::vector<uint32_t>> genome_positions;
        
        // Extract minimizers from read
        for (size_t i = 0; i <= read.length() - KMER_SIZE; i++) {
            std::string kmer = read.substr(i, KMER_SIZE);
            if (kmer.find('N') != std::string::npos) continue;
            
            std::string rc_kmer = reverse_complement(kmer);
            uint64_t hash = hash_kmer(kmer < rc_kmer ? kmer : rc_kmer);
            
            // Look up in index
            auto range = index.equal_range(hash);
            for (auto it = range.first; it != range.second; ++it) {
                genome_hits[it->second.genome_id]++;
                genome_positions[it->second.genome_id].push_back(it->second.position);
            }
        }
        
        // Find best match
        uint32_t max_hits = 0;
        for (const auto& [genome_id, hits] : genome_hits) {
            if (hits > max_hits) {
                max_hits = hits;
                result.best_genome_id = genome_id;
            }
        }
        
        // Calculate confidence
        uint32_t total_hits = 0;
        for (const auto& [genome_id, hits] : genome_hits) {
            total_hits += hits;
            if (hits >= max_hits * 0.8) {  // Within 80% of best hit
                result.multi_map_genomes.push_back(genome_id);
            }
        }
        
        result.confidence = (total_hits > 0) ? (float)max_hits / total_hits : 0.0;
        
        // Check if any hits are in resistance regions
        if (genome_positions.count(result.best_genome_id) > 0) {
            const auto& positions = genome_positions[result.best_genome_id];
            const auto& regions = genome_info[result.best_genome_id].resistance_regions;
            
            for (uint32_t pos : positions) {
                for (const auto& [start, end] : regions) {
                    if (pos >= start && pos <= end) {
                        result.has_resistance_region = true;
                        break;
                    }
                }
                if (result.has_resistance_region) break;
            }
        }
        
        return result;
    }
    
    // Batch classification (simulates GPU parallel processing)
    std::vector<ClassificationResult> classify_batch(const std::vector<std::string>& reads) {
        std::vector<ClassificationResult> results;
        results.reserve(reads.size());
        
        #pragma omp parallel for
        for (size_t i = 0; i < reads.size(); i++) {
            results.push_back(classify_read(reads[i], i));
        }
        
        return results;
    }
    
    // Get genome information
    const GenomeInfo& get_genome_info(uint32_t genome_id) const {
        return genome_info[genome_id];
    }
    
    // Save index to file
    void save_index(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        
        // Write genome info
        size_t num_genomes = genome_info.size();
        out.write(reinterpret_cast<const char*>(&num_genomes), sizeof(num_genomes));
        for (const auto& info : genome_info) {
            size_t name_len = info.name.length();
            out.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
            out.write(info.name.c_str(), name_len);
            out.write(reinterpret_cast<const char*>(&info.genome_id), sizeof(info.genome_id));
            out.write(reinterpret_cast<const char*>(&info.length), sizeof(info.length));
        }
        
        // Write index
        size_t index_size = index.size();
        out.write(reinterpret_cast<const char*>(&index_size), sizeof(index_size));
        for (const auto& [hash, hit] : index) {
            out.write(reinterpret_cast<const char*>(&hash), sizeof(hash));
            out.write(reinterpret_cast<const char*>(&hit), sizeof(hit));
        }
        
        out.close();
    }
};

// Example usage
#ifdef TEST
int main() {
    // Create test data
    std::vector<std::pair<std::string, std::string>> test_genomes = {
        {"E.coli_K12", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"},
        {"S.aureus_USA300", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"},
        {"P.aeruginosa_PAO1", "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"}
    };
    
    // Build index
    MinimizerIndex index;
    index.build_index(test_genomes);
    
    // Test classification
    std::vector<std::string> test_reads = {
        "ATCGATCGATCGATCGATCGATCGATCG",  // Should match E.coli
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTA",  // Should match S.aureus
        "TACGTACGTACGTACGTACGTACGTACG",  // Should match P.aeruginosa
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  // Should have low confidence
    };
    
    auto results = index.classify_batch(test_reads);
    
    // Print results
    for (size_t i = 0; i < results.size(); i++) {
        const auto& result = results[i];
        std::cout << "Read " << i << ": ";
        if (result.confidence > 0) {
            const auto& genome = index.get_genome_info(result.best_genome_id);
            std::cout << genome.name << " (confidence: " << result.confidence << ")";
            if (result.has_resistance_region) {
                std::cout << " [RESISTANCE REGION DETECTED]";
            }
            if (result.multi_map_genomes.size() > 1) {
                std::cout << " [MULTI-MAPPING]";
            }
        } else {
            std::cout << "UNCLASSIFIED";
        }
        std::cout << "\n";
    }
    
    return 0;
}
#endif // TEST