// hierarchical_gpu_database.h
#ifndef HIERARCHICAL_GPU_DATABASE_H
#define HIERARCHICAL_GPU_DATABASE_H

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <deque>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cuda_runtime.h>
#include "gpu_kmer_database.h"

namespace biogpu {

// Forward declarations
struct DatabaseTier;
class TierManager;

// Configuration for the hierarchical database
struct HierarchicalDBConfig {
    size_t max_gpu_memory_gb = 10;          // Maximum GPU memory to use
    size_t tier_size_mb = 512;              // Size of each tier in MB
    size_t cache_tiers = 3;                 // Number of tiers to keep cached
    bool use_lru_eviction = true;           // Use LRU for tier eviction
    bool preload_frequent_tiers = true;     // Preload most accessed tiers
};

// Statistics for performance monitoring
struct HierarchicalDBStats {
    uint64_t total_lookups = 0;
    uint64_t cache_hits = 0;
    uint64_t cache_misses = 0;
    uint64_t tier_loads = 0;
    uint64_t tier_evictions = 0;
    double avg_lookup_time_us = 0.0;
    
    // Additional stats for compatibility
    size_t num_levels = 1;
    size_t total_kmers = 0;
    size_t l1_buckets = 0;
    size_t l1_kmers = 0;
    size_t l2_shards = 0;
    size_t l2_kmers = 0;
    size_t memory_bytes = 0;
    
    double get_cache_hit_rate() const {
        return total_lookups > 0 ? (double)cache_hits / total_lookups : 0.0;
    }
};

// Individual database tier
struct DatabaseTier {
    std::string filename;
    uint64_t start_hash;           // First hash in this tier
    uint64_t end_hash;             // Last hash in this tier  
    size_t size_bytes;             // Size on disk
    size_t num_entries;            // Number of k-mer entries
    bool loaded = false;           // Currently in GPU memory
    uint64_t last_access_time = 0; // For LRU eviction
    uint64_t access_count = 0;     // For frequency-based loading
    
    // GPU memory for this tier
    KmerEntry* d_entries = nullptr;
    HashTableParams gpu_params;
    
    // For building - temporary host storage
    std::vector<std::pair<uint64_t, uint32_t>> host_entries;
    
    void clear_gpu_memory() {
        if (d_entries) {
            cudaFree(d_entries);
            d_entries = nullptr;
        }
        loaded = false;
    }
    
    ~DatabaseTier() {
        clear_gpu_memory();
    }
};

class HierarchicalGPUDatabase {
private:
    HierarchicalDBConfig config;
    std::vector<std::unique_ptr<DatabaseTier>> tiers;
    std::deque<size_t> lru_order;           // Indices of tiers in LRU order
    size_t current_gpu_usage = 0;           // Current GPU memory usage
    uint64_t access_counter = 0;            // Global access counter for LRU
    HierarchicalDBStats stats;
    bool sorted_by_frequency = true;        // How the database was built
    
    // GPU streams for async operations
    cudaStream_t load_stream;
    cudaStream_t lookup_stream;
    
    // Hash partitioning
    static const int HASH_PARTITION_BITS = 16;  // 64K partitions max
    
public:
    HierarchicalGPUDatabase(const HierarchicalDBConfig& cfg = HierarchicalDBConfig());
    ~HierarchicalGPUDatabase();
    
    // Building interface
    void start_building(const std::string& output_prefix);
    void add_kmer(uint64_t hash, uint32_t taxon_id);
    void finalize_build();
    
    // Loading interface  
    void load_database(const std::string& db_prefix);
    
    // Lookup interface
    void lookup_batch_gpu(const uint64_t* d_hashes, 
                         uint32_t* d_taxon_ids, 
                         size_t count,
                         cudaStream_t stream = 0);
    
    // Management
    void preload_frequent_tiers(size_t top_n = 5);
    void clear_cache();
    HierarchicalDBStats get_statistics() const { return stats; }
    void print_statistics() const;
    
    // Configuration
    void set_max_memory(size_t gb) { config.max_gpu_memory_gb = gb; }
    size_t get_tier_count() const { return tiers.size(); }
    
    // Additional methods for compatibility
    void save(const std::string& filename) const {
        // Simplified save implementation
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Cannot create database file: " + filename);
        }
        
        // Write header
        size_t num_tiers = tiers.size();
        size_t total_kmers = 0;
        for (const auto& tier : tiers) {
            total_kmers += tier->num_entries;
        }
        out.write(reinterpret_cast<const char*>(&num_tiers), sizeof(num_tiers));
        out.write(reinterpret_cast<const char*>(&total_kmers), sizeof(total_kmers));
        
        // Write tier information
        for (const auto& tier : tiers) {
            out.write(reinterpret_cast<const char*>(&tier->start_hash), sizeof(tier->start_hash));
            out.write(reinterpret_cast<const char*>(&tier->end_hash), sizeof(tier->end_hash));
            out.write(reinterpret_cast<const char*>(&tier->num_entries), sizeof(tier->num_entries));
        }
        
        out.close();
        std::cout << "Database saved to: " << filename << "\n";
    }
    
    void load(const std::string& filename) {
        load_database(filename);
    }
    
    void build_hierarchy() {
        finalize_build();
    }
    
    // CPU batch lookup for testing
    void lookup_batch(const std::vector<uint64_t>& hashes, std::vector<uint32_t>& results) {
        results.resize(hashes.size());
        
        // Allocate GPU memory for batch
        uint64_t* d_hashes;
        uint32_t* d_results;
        cudaMalloc(&d_hashes, hashes.size() * sizeof(uint64_t));
        cudaMalloc(&d_results, results.size() * sizeof(uint32_t));
        
        // Copy to GPU
        cudaMemcpy(d_hashes, hashes.data(), hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        
        // Do lookup
        lookup_batch_gpu(d_hashes, d_results, hashes.size());
        
        // Copy back
        cudaMemcpy(results.data(), d_results, results.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        // Cleanup
        cudaFree(d_hashes);
        cudaFree(d_results);
    }
    
private:
    // Tier management
    int find_tier_for_hash(uint64_t hash) const;
    void ensure_tier_loaded(size_t tier_idx);
    void load_tier(size_t tier_idx);
    void evict_lru_tier();
    void update_lru(size_t tier_idx);
    
    // Building helpers
    void partition_kmers_by_hash(const std::string& output_prefix);
    void build_tier_files(const std::string& output_prefix);
    
    // Memory management
    size_t estimate_tier_gpu_memory(const DatabaseTier& tier) const;
    bool has_memory_for_tier(size_t tier_idx) const;
};

// Implementation details

HierarchicalGPUDatabase::HierarchicalGPUDatabase(const HierarchicalDBConfig& cfg) 
    : config(cfg) {
    
    // Create CUDA streams
    cudaStreamCreate(&load_stream);
    cudaStreamCreate(&lookup_stream);
    
    std::cout << "Hierarchical database initialized:\n";
    std::cout << "  Max GPU memory: " << config.max_gpu_memory_gb << " GB\n";
    std::cout << "  Tier size: " << config.tier_size_mb << " MB\n";
    std::cout << "  Cache tiers: " << config.cache_tiers << "\n";
}

HierarchicalGPUDatabase::~HierarchicalGPUDatabase() {
    clear_cache();
    
    if (load_stream) cudaStreamDestroy(load_stream);
    if (lookup_stream) cudaStreamDestroy(lookup_stream);
}

void HierarchicalGPUDatabase::add_kmer(uint64_t hash, uint32_t taxon_id) {
    // Find which tier this k-mer belongs to based on hash
    uint64_t partition = hash >> (64 - HASH_PARTITION_BITS);
    
    // Ensure we have enough tiers
    size_t tier_idx = partition % std::max(size_t(1), tiers.size());
    
    if (tier_idx >= tiers.size()) {
        // Create new tier
        auto tier = std::make_unique<DatabaseTier>();
        tier->start_hash = partition << (64 - HASH_PARTITION_BITS);
        tier->end_hash = ((partition + 1) << (64 - HASH_PARTITION_BITS)) - 1;
        tiers.push_back(std::move(tier));
        tier_idx = tiers.size() - 1;
    }
    
    // Add to appropriate tier
    tiers[tier_idx]->host_entries.push_back({hash, taxon_id});
    tiers[tier_idx]->num_entries++;
}

void HierarchicalGPUDatabase::finalize_build() {
    std::cout << "Finalizing hierarchical database with " << tiers.size() << " tiers\n";
    
    // Sort each tier by hash for efficient GPU lookup
    stats.total_kmers = 0;
    stats.memory_bytes = 0;
    
    for (auto& tier : tiers) {
        std::sort(tier->host_entries.begin(), tier->host_entries.end());
        tier->size_bytes = tier->host_entries.size() * sizeof(KmerEntry);
        
        stats.total_kmers += tier->num_entries;
        stats.memory_bytes += tier->size_bytes;
        
        std::cout << "Tier with " << tier->num_entries << " entries ("
                  << tier->size_bytes / 1024 / 1024 << " MB)\n";
    }
    
    // Update stats
    stats.num_levels = 1;
    stats.l1_buckets = tiers.size();
    stats.l1_kmers = stats.total_kmers;
    stats.l2_shards = 0;
    stats.l2_kmers = 0;
    
    // Preload frequent tiers if enabled
    if (config.preload_frequent_tiers) {
        preload_frequent_tiers(config.cache_tiers);
    }
}

int HierarchicalGPUDatabase::find_tier_for_hash(uint64_t hash) const {
    // Use binary search since tiers are sorted by hash range
    int left = 0, right = tiers.size() - 1;
    
    while (left <= right) {
        int mid = left + (right - left) / 2;
        
        if (hash >= tiers[mid]->start_hash && hash <= tiers[mid]->end_hash) {
            return mid;
        } else if (hash < tiers[mid]->start_hash) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    
    return -1;  // Not found
}

void HierarchicalGPUDatabase::ensure_tier_loaded(size_t tier_idx) {
    auto& tier = tiers[tier_idx];
    
    if (tier->loaded) {
        // Update access statistics
        tier->access_count++;
        tier->last_access_time = ++access_counter;
        update_lru(tier_idx);
        stats.cache_hits++;
        return;
    }
    
    stats.cache_misses++;
    
    // Check if we need to evict tiers first
    while (!has_memory_for_tier(tier_idx)) {
        evict_lru_tier();
    }
    
    // Load the tier
    load_tier(tier_idx);
    update_lru(tier_idx);
}

void HierarchicalGPUDatabase::load_tier(size_t tier_idx) {
    auto& tier = tiers[tier_idx];
    
    std::cout << "Loading tier " << tier_idx << " (" 
              << tier->size_bytes / 1024 / 1024 << " MB, "
              << tier->num_entries << " entries)\n";
    
    // Allocate GPU memory
    size_t gpu_size = tier->num_entries * sizeof(KmerEntry);
    cudaMalloc(&tier->d_entries, gpu_size);
    
    // Copy data to GPU (async if we have host data)
    if (!tier->host_entries.empty()) {
        // Convert host entries to GPU format
        std::vector<KmerEntry> gpu_entries;
        gpu_entries.reserve(tier->host_entries.size());
        
        for (const auto& [hash, taxon] : tier->host_entries) {
            gpu_entries.push_back({hash, taxon, 0});
        }
        
        cudaMemcpyAsync(tier->d_entries, gpu_entries.data(),
                       gpu_size, cudaMemcpyHostToDevice, load_stream);
    } else {
        // Load from file
        std::ifstream file(tier->filename, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open tier file: " + tier->filename);
        }
        
        // Skip the tier header (64 bytes)
        struct TierHeader {
            uint64_t num_entries;
            uint64_t min_hash;
            uint64_t max_hash;
            uint64_t reserved[5];
        } header;
        
        file.read(reinterpret_cast<char*>(&header), sizeof(header));
        
        // Read k-mer entries
        std::vector<KmerEntry> file_entries(tier->num_entries);
        file.read(reinterpret_cast<char*>(file_entries.data()), gpu_size);
        file.close();
        
        cudaMemcpyAsync(tier->d_entries, file_entries.data(),
                       gpu_size, cudaMemcpyHostToDevice, load_stream);
        
        // Debug: show first few entries
        static int debug_tier_count = 0;
        if (debug_tier_count++ < 3) {
            std::cout << "  First few entries in tier " << tier_idx << ":\n";
            for (size_t i = 0; i < std::min(size_t(5), tier->num_entries); i++) {
                std::cout << "    Entry[" << i << "]: hash=" << file_entries[i].hash
                         << ", taxon=" << file_entries[i].taxon_id << "\n";
            }
        }
    }
    
    // Wait for transfer to complete
    cudaStreamSynchronize(load_stream);
    
    // Update status
    tier->loaded = true;
    tier->access_count++;
    tier->last_access_time = ++access_counter;
    current_gpu_usage += gpu_size;
    stats.tier_loads++;
}

void HierarchicalGPUDatabase::evict_lru_tier() {
    if (lru_order.empty()) return;
    
    // Find least recently used tier that's currently loaded
    for (auto it = lru_order.rbegin(); it != lru_order.rend(); ++it) {
        size_t tier_idx = *it;
        auto& tier = tiers[tier_idx];
        
        if (tier->loaded) {
            std::cout << "Evicting tier " << tier_idx << "\n";
            
            // Free GPU memory
            size_t freed_memory = tier->num_entries * sizeof(KmerEntry);
            tier->clear_gpu_memory();
            current_gpu_usage -= freed_memory;
            
            // Remove from LRU order
            lru_order.erase(std::next(it).base());
            stats.tier_evictions++;
            return;
        }
    }
}

void HierarchicalGPUDatabase::update_lru(size_t tier_idx) {
    // Remove from current position
    lru_order.erase(std::remove(lru_order.begin(), lru_order.end(), tier_idx),
                   lru_order.end());
    
    // Add to front (most recently used)
    lru_order.push_front(tier_idx);
}

bool HierarchicalGPUDatabase::has_memory_for_tier(size_t tier_idx) const {
    size_t needed = estimate_tier_gpu_memory(*tiers[tier_idx]);
    size_t available = (config.max_gpu_memory_gb * 1024ULL * 1024 * 1024) - current_gpu_usage;
    return needed <= available;
}

size_t HierarchicalGPUDatabase::estimate_tier_gpu_memory(const DatabaseTier& tier) const {
    return tier.num_entries * sizeof(KmerEntry);
}

// GPU kernel for hierarchical lookup
#ifdef __CUDACC__
__global__ void hierarchical_lookup_kernel(
    const uint64_t* query_hashes,
    const KmerEntry* tier_entries,
    const size_t tier_size,
    uint32_t* results,
    const size_t num_queries,
    const bool sorted_by_hash
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;
    
    // Skip if already found in previous tier
    if (results[tid] != 0) return;
    
    uint64_t target_hash = query_hashes[tid];
    
    if (sorted_by_hash) {
        // Binary search in sorted tier
        int left = 0, right = tier_size - 1;
        
        while (left <= right) {
            int mid = left + (right - left) / 2;
            uint64_t mid_hash = tier_entries[mid].hash;
            
            if (mid_hash == target_hash) {
                results[tid] = tier_entries[mid].taxon_id;
                return;
            } else if (mid_hash < target_hash) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
    } else {
        // Linear search - tiers are sorted by frequency
        for (size_t i = 0; i < tier_size; i++) {
            if (tier_entries[i].hash == target_hash) {
                results[tid] = tier_entries[i].taxon_id;
                return;
            }
        }
    }
}
#endif // __CUDACC__

void HierarchicalGPUDatabase::lookup_batch_gpu(const uint64_t* d_hashes, 
                                              uint32_t* d_taxon_ids, 
                                              size_t count,
                                              cudaStream_t stream) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Initialize all results to 0 (not found)
    cudaMemset(d_taxon_ids, 0, count * sizeof(uint32_t));
    
    // Simple approach: check all tiers for all queries
    // Since tiers are sorted by frequency, we check most frequent first
    size_t found_count = 0;
    
    for (size_t tier_idx = 0; tier_idx < tiers.size(); tier_idx++) {
        // Ensure tier is loaded
        ensure_tier_loaded(tier_idx);
        auto& tier = tiers[tier_idx];
        
        // Launch lookup kernel for all queries against this tier
        int block_size = 256;
        int num_blocks = (count + block_size - 1) / block_size;
        
#ifdef __CUDACC__
        hierarchical_lookup_kernel<<<num_blocks, block_size, 0, stream>>>(
            d_hashes, tier->d_entries, tier->num_entries,
            d_taxon_ids, count, !sorted_by_frequency  // pass sorted_by_hash flag
        );
#endif
        cudaStreamSynchronize(stream);
    }
    
    // Debug: count how many were found
    std::vector<uint32_t> h_results(count);
    cudaMemcpy(h_results.data(), d_taxon_ids, count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    for (auto taxon : h_results) {
        if (taxon > 0) found_count++;
    }
    
    if (found_count == 0 && count > 0) {
        static int debug_count = 0;
        if (debug_count++ < 5) {
            std::cout << "DEBUG: No hits found in batch of " << count << " queries\n";
            std::cout << "  Tiers checked: " << tiers.size() << "\n";
            // Sample first few hashes
            std::vector<uint64_t> h_hashes(std::min(count, size_t(5)));
            cudaMemcpy(h_hashes.data(), d_hashes, h_hashes.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost);
            for (size_t i = 0; i < h_hashes.size(); i++) {
                std::cout << "  Query hash[" << i << "]: " << h_hashes[i] << "\n";
            }
        }
    }
    
    // Update statistics
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    
    stats.total_lookups += count;
    if (stats.total_lookups > 0) {
        stats.avg_lookup_time_us = (stats.avg_lookup_time_us * (stats.total_lookups - count) + 
                                   duration.count()) / stats.total_lookups;
    } else {
        stats.avg_lookup_time_us = duration.count();
    }
}

void HierarchicalGPUDatabase::preload_frequent_tiers(size_t top_n) {
    // Sort tiers by access frequency
    std::vector<std::pair<uint64_t, size_t>> tier_frequencies;
    for (size_t i = 0; i < tiers.size(); i++) {
        tier_frequencies.push_back({tiers[i]->access_count, i});
    }
    
    std::sort(tier_frequencies.rbegin(), tier_frequencies.rend());
    
    // Load top N tiers
    std::cout << "Preloading " << std::min(top_n, tier_frequencies.size()) << " frequent tiers\n";
    
    for (size_t i = 0; i < std::min(top_n, tier_frequencies.size()); i++) {
        size_t tier_idx = tier_frequencies[i].second;
        if (!tiers[tier_idx]->loaded && has_memory_for_tier(tier_idx)) {
            load_tier(tier_idx);
            update_lru(tier_idx);
        }
    }
}

void HierarchicalGPUDatabase::print_statistics() const {
    std::cout << "\n=== Hierarchical Database Statistics ===\n";
    std::cout << "Total lookups: " << stats.total_lookups << "\n";
    std::cout << "Cache hit rate: " << (stats.get_cache_hit_rate() * 100.0) << "%\n";
    std::cout << "Tier loads: " << stats.tier_loads << "\n";
    std::cout << "Tier evictions: " << stats.tier_evictions << "\n";
    std::cout << "Average lookup time: " << stats.avg_lookup_time_us << " μs\n";
    std::cout << "GPU memory usage: " << (current_gpu_usage / 1024.0 / 1024.0) << " MB\n";
    std::cout << "Active tiers: " << lru_order.size() << "/" << tiers.size() << "\n";
}

void HierarchicalGPUDatabase::clear_cache() {
    // Clear all loaded tiers from GPU memory
    for (auto& tier : tiers) {
        if (tier->loaded) {
            tier->clear_gpu_memory();
        }
    }
    
    // Reset tracking
    current_gpu_usage = 0;
    lru_order.clear();
}

void HierarchicalGPUDatabase::load_database(const std::string& db_prefix) {
    // Clear any existing data
    clear_cache();
    tiers.clear();
    
    std::cout << "Loading hierarchical database from: " << db_prefix << "\n";
    
    // Read manifest file
    std::string manifest_path = db_prefix + "/manifest.json";
    std::ifstream manifest_file(manifest_path);
    if (!manifest_file.is_open()) {
        throw std::runtime_error("Cannot open manifest file: " + manifest_path);
    }
    
    // Parse manifest JSON (simple parsing for our structured format)
    std::string line;
    size_t total_kmers = 0;
    size_t num_tiers = 0;
    
    while (std::getline(manifest_file, line)) {
        if (line.find("\"total_kmers\":") != std::string::npos) {
            sscanf(line.c_str(), "  \"total_kmers\": %zu,", &total_kmers);
        } else if (line.find("\"num_tiers\":") != std::string::npos) {
            sscanf(line.c_str(), "  \"num_tiers\": %zu,", &num_tiers);
        } else if (line.find("\"sorted_by_frequency\":") != std::string::npos) {
            sorted_by_frequency = (line.find("true") != std::string::npos);
        }
    }
    manifest_file.close();
    
    // Re-read manifest to get tier information
    manifest_file.open(manifest_path);
    bool in_tiers_section = false;
    size_t current_tier_idx = 0;
    
    while (std::getline(manifest_file, line)) {
        if (line.find("\"tiers\":") != std::string::npos) {
            in_tiers_section = true;
            continue;
        }
        
        if (in_tiers_section && line.find("\"tier_id\":") != std::string::npos) {
            auto tier = std::make_unique<DatabaseTier>();
            
            // Read tier metadata
            size_t tier_id;
            sscanf(line.c_str(), "      \"tier_id\": %zu,", &tier_id);
            
            // Read subsequent lines for this tier
            std::getline(manifest_file, line); // filename
            char filename[256];
            sscanf(line.c_str(), "      \"filename\": \"%255[^\"]\",", filename);
            tier->filename = db_prefix + "/" + std::string(filename);
            
            std::getline(manifest_file, line); // num_entries
            sscanf(line.c_str(), "      \"num_entries\": %zu,", &tier->num_entries);
            
            std::getline(manifest_file, line); // size_bytes
            sscanf(line.c_str(), "      \"size_bytes\": %zu,", &tier->size_bytes);
            
            std::getline(manifest_file, line); // min_hash
            sscanf(line.c_str(), "      \"min_hash\": %lu,", &tier->start_hash);
            
            std::getline(manifest_file, line); // max_hash
            sscanf(line.c_str(), "      \"max_hash\": %lu", &tier->end_hash);
            
            tiers.push_back(std::move(tier));
            current_tier_idx++;
            
            if (current_tier_idx >= num_tiers) break;
        }
    }
    manifest_file.close();
    
    // Initialize stats
    stats = HierarchicalDBStats();
    stats.total_kmers = total_kmers;
    stats.num_levels = 1;
    stats.l1_buckets = tiers.size();
    stats.l1_kmers = total_kmers;
    stats.l2_shards = 0;
    stats.l2_kmers = 0;
    
    // Calculate total memory requirement
    stats.memory_bytes = 0;
    for (const auto& tier : tiers) {
        stats.memory_bytes += tier->size_bytes;
    }
    
    std::cout << "Loaded database metadata:\n";
    std::cout << "  Total k-mers: " << stats.total_kmers << "\n";
    std::cout << "  Number of tiers: " << tiers.size() << "\n";
    std::cout << "  Total size: " << (stats.memory_bytes / 1024.0 / 1024.0 / 1024.0) << " GB\n";
    std::cout << "  Sorted by: " << (sorted_by_frequency ? "frequency" : "hash") << "\n";
    
    // Preload frequent tiers if configured
    if (config.preload_frequent_tiers && config.cache_tiers > 0) {
        preload_frequent_tiers(config.cache_tiers);
    }
}

} // namespace biogpu

#endif // HIERARCHICAL_GPU_DATABASE_H