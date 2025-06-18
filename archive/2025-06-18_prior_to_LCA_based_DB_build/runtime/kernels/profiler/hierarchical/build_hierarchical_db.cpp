// build_hierarchical_db.cpp - Build hierarchical database from k-mer list
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <cstring>
#include "hierarchical_gpu_database.h"

using namespace biogpu;

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

// Configuration options
struct BuildConfig {
    size_t tier_size_mb = 512;           // Target size per tier
    size_t max_kmers_per_tier = 0;       // Auto-calculated from tier size
    bool sort_by_frequency = true;       // Put frequent k-mers in early tiers
    bool create_manifest = true;         // Create database manifest file
    size_t progress_interval = 1000000;  // Progress update interval
};

class HierarchicalDatabaseBuilder {
private:
    BuildConfig config;
    std::string output_prefix;
    
    // K-mer data with frequency counting
    std::map<uint64_t, std::pair<uint32_t, uint32_t>> kmer_data; // hash -> (taxon, frequency)
    
    // Tier information
    struct TierInfo {
        std::string filename;
        std::vector<std::pair<uint64_t, uint32_t>> entries;
        size_t estimated_size_bytes = 0;
        uint64_t min_hash = UINT64_MAX;
        uint64_t max_hash = 0;
    };
    
    std::vector<TierInfo> tiers;
    
public:
    HierarchicalDatabaseBuilder(const std::string& prefix, const BuildConfig& cfg = BuildConfig())
        : output_prefix(prefix), config(cfg) {
        
        // Calculate max k-mers per tier based on target size
        if (config.max_kmers_per_tier == 0) {
            size_t bytes_per_entry = sizeof(KmerEntry);  // 16 bytes
            config.max_kmers_per_tier = (config.tier_size_mb * 1024 * 1024) / bytes_per_entry;
        }
        
        std::cout << "Building hierarchical database:\n";
        std::cout << "  Output prefix: " << output_prefix << "\n";
        std::cout << "  Target tier size: " << config.tier_size_mb << " MB\n";
        std::cout << "  Max k-mers per tier: " << config.max_kmers_per_tier << "\n";
    }
    
    void read_kmer_file(const std::string& kmer_file) {
        std::cout << "Reading k-mer file: " << kmer_file << "\n";
        
        std::ifstream infile(kmer_file);
        if (!infile.is_open()) {
            throw std::runtime_error("Cannot open k-mer file: " + kmer_file);
        }
        
        std::string line;
        size_t line_num = 0;
        size_t kmers_read = 0;
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
            
            // Get canonical k-mer and compute hash
            uint64_t canonical = get_canonical_kmer(kmer);
            if (canonical == 0) {
                std::cerr << "Warning: Invalid k-mer at line " << line_num << ": " << kmer << "\n";
                continue;
            }
            
            uint64_t hash = murmur_hash3_finalizer(canonical);
            
            // Track frequency and taxon
            if (kmer_data.count(hash)) {
                kmer_data[hash].second++;  // Increment frequency
            } else {
                kmer_data[hash] = {taxon_id, 1};
            }
            
            kmers_read++;
            max_taxon_id = std::max(max_taxon_id, (size_t)taxon_id);
            
            if (kmers_read % config.progress_interval == 0) {
                std::cout << "  Read " << kmers_read << " k-mers...\n";
            }
        }
        
        std::cout << "K-mer reading complete:\n";
        std::cout << "  Total k-mers: " << kmers_read << "\n";
        std::cout << "  Unique k-mers: " << kmer_data.size() << "\n";
        std::cout << "  Maximum taxon ID: " << max_taxon_id << "\n";
    }
    
    void create_tiers() {
        std::cout << "\nCreating database tiers...\n";
        
        // Convert map to vector for sorting
        std::vector<std::tuple<uint64_t, uint32_t, uint32_t>> all_kmers;
        all_kmers.reserve(kmer_data.size());
        
        for (const auto& [hash, data] : kmer_data) {
            all_kmers.push_back({hash, data.first, data.second});
        }
        
        // Sort by frequency (most frequent first) if enabled
        if (config.sort_by_frequency) {
            std::cout << "Sorting k-mers by frequency...\n";
            std::sort(all_kmers.begin(), all_kmers.end(),
                     [](const auto& a, const auto& b) {
                         return std::get<2>(a) > std::get<2>(b);  // Sort by frequency descending
                     });
        } else {
            std::cout << "Sorting k-mers by hash...\n";
            std::sort(all_kmers.begin(), all_kmers.end());
        }
        
        // Partition into tiers
        size_t total_tiers = (all_kmers.size() + config.max_kmers_per_tier - 1) / config.max_kmers_per_tier;
        tiers.resize(total_tiers);
        
        std::cout << "Creating " << total_tiers << " tiers...\n";
        
        for (size_t i = 0; i < all_kmers.size(); i++) {
            size_t tier_idx = i / config.max_kmers_per_tier;
            
            auto& tier = tiers[tier_idx];
            auto& [hash, taxon, freq] = all_kmers[i];
            
            tier.entries.push_back({hash, taxon});
            tier.estimated_size_bytes += sizeof(KmerEntry);
            
            // Track hash range for each tier
            tier.min_hash = std::min(tier.min_hash, hash);
            tier.max_hash = std::max(tier.max_hash, hash);
        }
        
        // If not sorted by frequency, sort each tier by hash for efficient GPU lookup
        if (!config.sort_by_frequency) {
            for (auto& tier : tiers) {
                std::sort(tier.entries.begin(), tier.entries.end());
                if (!tier.entries.empty()) {
                    tier.min_hash = tier.entries.front().first;
                    tier.max_hash = tier.entries.back().first;
                }
            }
        }
        
        std::cout << "Tier creation complete:\n";
        for (size_t i = 0; i < tiers.size(); i++) {
            std::cout << "  Tier " << i << ": " << tiers[i].entries.size() 
                      << " k-mers, " << (tiers[i].estimated_size_bytes / 1024.0 / 1024.0) 
                      << " MB\n";
        }
    }
    
    void write_tier_files() {
        std::cout << "\nWriting tier files...\n";
        
        // Create output directory
        std::filesystem::create_directories(output_prefix);
        
        for (size_t i = 0; i < tiers.size(); i++) {
            auto& tier = tiers[i];
            tier.filename = output_prefix + "/tier_" + std::to_string(i) + ".bin";
            
            std::ofstream outfile(tier.filename, std::ios::binary);
            if (!outfile.is_open()) {
                throw std::runtime_error("Cannot create tier file: " + tier.filename);
            }
            
            // Write tier header
            struct TierHeader {
                uint64_t num_entries;
                uint64_t min_hash;
                uint64_t max_hash;
                uint64_t reserved[5];  // For future use
            } header;
            
            header.num_entries = tier.entries.size();
            header.min_hash = tier.min_hash;
            header.max_hash = tier.max_hash;
            memset(header.reserved, 0, sizeof(header.reserved));
            
            outfile.write(reinterpret_cast<const char*>(&header), sizeof(header));
            
            // Write k-mer entries
            for (const auto& [hash, taxon] : tier.entries) {
                KmerEntry entry = {hash, taxon, 0};
                outfile.write(reinterpret_cast<const char*>(&entry), sizeof(entry));
            }
            
            outfile.close();
            
            std::cout << "  Wrote " << tier.filename << " (" 
                      << (tier.estimated_size_bytes / 1024.0 / 1024.0) << " MB)\n";
        }
    }
    
    void create_manifest() {
        if (!config.create_manifest) return;
        
        std::cout << "\nCreating database manifest...\n";
        
        std::string manifest_file = output_prefix + "/manifest.json";
        std::ofstream manifest(manifest_file);
        
        manifest << "{\n";
        manifest << "  \"version\": \"1.0\",\n";
        manifest << "  \"type\": \"hierarchical_kmer_database\",\n";
        manifest << "  \"created\": \"" << std::time(nullptr) << "\",\n";
        manifest << "  \"total_kmers\": " << kmer_data.size() << ",\n";
        manifest << "  \"num_tiers\": " << tiers.size() << ",\n";
        manifest << "  \"tier_size_mb\": " << config.tier_size_mb << ",\n";
        manifest << "  \"sorted_by_frequency\": " << (config.sort_by_frequency ? "true" : "false") << ",\n";
        manifest << "  \"tiers\": [\n";
        
        for (size_t i = 0; i < tiers.size(); i++) {
            const auto& tier = tiers[i];
            manifest << "    {\n";
            manifest << "      \"tier_id\": " << i << ",\n";
            manifest << "      \"filename\": \"tier_" << i << ".bin\",\n";
            manifest << "      \"num_entries\": " << tier.entries.size() << ",\n";
            manifest << "      \"size_bytes\": " << tier.estimated_size_bytes << ",\n";
            manifest << "      \"min_hash\": " << tier.min_hash << ",\n";
            manifest << "      \"max_hash\": " << tier.max_hash << "\n";
            manifest << "    }";
            if (i < tiers.size() - 1) manifest << ",";
            manifest << "\n";
        }
        
        manifest << "  ]\n";
        manifest << "}\n";
        
        manifest.close();
        std::cout << "Created manifest: " << manifest_file << "\n";
    }
    
    void print_summary() {
        std::cout << "\n=== Database Build Summary ===\n";
        std::cout << "Output directory: " << output_prefix << "\n";
        std::cout << "Total unique k-mers: " << kmer_data.size() << "\n";
        std::cout << "Number of tiers: " << tiers.size() << "\n";
        
        size_t total_size = 0;
        for (const auto& tier : tiers) {
            total_size += tier.estimated_size_bytes;
        }
        
        std::cout << "Total database size: " << (total_size / 1024.0 / 1024.0) << " MB\n";
        std::cout << "Average tier size: " << (total_size / tiers.size() / 1024.0 / 1024.0) << " MB\n";
        
        // Memory efficiency analysis
        std::cout << "\n=== Memory Efficiency Analysis ===\n";
        std::cout << "For " << config.tier_size_mb << " MB GPU memory budget:\n";
        std::cout << "  Tiers that fit in memory: " << (config.tier_size_mb * 1024 * 1024) / (total_size / tiers.size()) << "\n";
        std::cout << "  Cache hit rate estimate: " << std::min(100.0, 100.0 * config.tier_size_mb / (total_size / 1024.0 / 1024.0)) << "%\n";
        
        std::cout << "\nUsage:\n";
        std::cout << "  Load in C++: HierarchicalGPUDatabase db; db.load_database(\"" << output_prefix << "\");\n";
        std::cout << "  With profiler: ./gpu_profiler_pipeline " << output_prefix << " reads.fastq\n";
    }
};

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <kmer_list.txt> <output_prefix> [options]\n";
    std::cerr << "\nOptions:\n";
    std::cerr << "  --tier-size <MB>      Target size per tier in MB (default: 512)\n";
    std::cerr << "  --sort-hash           Sort by hash instead of frequency\n";
    std::cerr << "  --no-manifest         Don't create manifest file\n";
    std::cerr << "  --progress <N>        Progress update interval (default: 1000000)\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << program_name << " database_kmers.txt pathogen_db --tier-size 256\n";
    std::cerr << "\nThis creates a hierarchical database optimized for streaming access\n";
    std::cerr << "with limited GPU memory. Each tier can be loaded independently.\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string kmer_file = argv[1];
    std::string output_prefix = argv[2];
    
    BuildConfig config;
    
    // Parse command line options
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--tier-size" && i + 1 < argc) {
            config.tier_size_mb = std::stoul(argv[++i]);
        } else if (arg == "--sort-hash") {
            config.sort_by_frequency = false;
        } else if (arg == "--no-manifest") {
            config.create_manifest = false;
        } else if (arg == "--progress" && i + 1 < argc) {
            config.progress_interval = std::stoul(argv[++i]);
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        HierarchicalDatabaseBuilder builder(output_prefix, config);
        
        // Step 1: Read k-mer data
        builder.read_kmer_file(kmer_file);
        
        // Step 2: Create tiers
        builder.create_tiers();
        
        // Step 3: Write tier files
        builder.write_tier_files();
        
        // Step 4: Create manifest
        builder.create_manifest();
        
        // Step 5: Summary
        builder.print_summary();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\nBuild completed in " << duration.count() << " seconds\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}