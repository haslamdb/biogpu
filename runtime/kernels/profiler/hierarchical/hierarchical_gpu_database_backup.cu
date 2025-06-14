// hierarchical_gpu_database.cu - Implementation of HierarchicalGPUDatabase methods
#include "hierarchical_gpu_database.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cuda_runtime.h>

namespace biogpu {

// Clear all cached tiers from GPU memory
void HierarchicalGPUDatabase::clear_cache() {
    for (auto& tier : tiers) {
        tier->clear_gpu_memory();
    }
    current_gpu_usage = 0;
    lru_order.clear();
}

// Load hierarchical database from files
void HierarchicalGPUDatabase::load_database(const std::string& db_prefix) {
    clear_cache();  // Clear any existing data
    
    // Load metadata file
    std::string metadata_file = db_prefix + ".meta";
    std::ifstream meta_in(metadata_file, std::ios::binary);
    if (!meta_in) {
        throw std::runtime_error("Cannot open metadata file: " + metadata_file);
    }
    
    // Read header
    size_t num_tiers;
    size_t total_entries;
    meta_in.read(reinterpret_cast<char*>(&num_tiers), sizeof(num_tiers));
    meta_in.read(reinterpret_cast<char*>(&total_entries), sizeof(total_entries));
    meta_in.read(reinterpret_cast<char*>(&config.tier_size_mb), sizeof(config.tier_size_mb));
    
    // Read tier information
    tiers.clear();
    tiers.reserve(num_tiers);
    
    stats.total_kmers = 0;
    
    for (size_t i = 0; i < num_tiers; i++) {
        auto tier = std::make_unique<DatabaseTier>();
        
        // Read tier metadata
        meta_in.read(reinterpret_cast<char*>(&tier->start_hash), sizeof(tier->start_hash));
        meta_in.read(reinterpret_cast<char*>(&tier->end_hash), sizeof(tier->end_hash));
        meta_in.read(reinterpret_cast<char*>(&tier->num_entries), sizeof(tier->num_entries));
        meta_in.read(reinterpret_cast<char*>(&tier->size_bytes), sizeof(tier->size_bytes));
        meta_in.read(reinterpret_cast<char*>(&tier->access_count), sizeof(tier->access_count));
        
        // Set file path for this tier
        tier->filename = db_prefix + ".tier" + std::to_string(i);
        tier->loaded = false;
        tier->d_entries = nullptr;
        
        stats.total_kmers += tier->num_entries;
        
        tiers.push_back(std::move(tier));
    }
    
    meta_in.close();
    
    // Initialize stats
    stats.num_levels = 1;  // Single level for now
    stats.l1_buckets = num_tiers;
    stats.l1_kmers = stats.total_kmers;
    stats.memory_bytes = 0;
    
    std::cout << "Loaded hierarchical database metadata:\n";
    std::cout << "  Tiers: " << tiers.size() << "\n";
    std::cout << "  Total k-mers: " << stats.total_kmers << "\n";
    std::cout << "  Tier size: " << config.tier_size_mb << " MB\n";
}

// Save hierarchical database to files
void HierarchicalGPUDatabase::save(const std::string& db_prefix) const {
    // Create metadata file
    std::string metadata_file = db_prefix + ".meta";
    std::ofstream meta_out(metadata_file, std::ios::binary);
    if (!meta_out) {
        throw std::runtime_error("Cannot create metadata file: " + metadata_file);
    }
    
    // Write header
    size_t num_tiers = tiers.size();
    meta_out.write(reinterpret_cast<const char*>(&num_tiers), sizeof(num_tiers));
    meta_out.write(reinterpret_cast<const char*>(&total_kmers), sizeof(total_kmers));
    meta_out.write(reinterpret_cast<const char*>(&max_tier_size), sizeof(max_tier_size));
    
    // Write tier information and save tier data
    for (const auto& tier : tiers) {
        meta_out.write(reinterpret_cast<const char*>(&tier.tier_id), sizeof(tier.tier_id));
        meta_out.write(reinterpret_cast<const char*>(&tier.num_kmers), sizeof(tier.num_kmers));
        meta_out.write(reinterpret_cast<const char*>(&tier.min_hash), sizeof(tier.min_hash));
        meta_out.write(reinterpret_cast<const char*>(&tier.max_hash), sizeof(tier.max_hash));
        meta_out.write(reinterpret_cast<const char*>(&tier.access_count), sizeof(tier.access_count));
        meta_out.write(reinterpret_cast<const char*>(&tier.last_access), sizeof(tier.last_access));
        
        // Save tier data to separate file
        std::string tier_file = db_prefix + ".tier" + std::to_string(tier.tier_id);
        std::ofstream tier_out(tier_file, std::ios::binary);
        if (!tier_out) {
            throw std::runtime_error("Cannot create tier file: " + tier_file);
        }
        
        // Write k-mer data for this tier
        if (tier.gpu_loaded && tier.d_kmers) {
            // Copy from GPU to host temporarily
            std::vector<KmerEntry> h_kmers(tier.num_kmers);
            cudaMemcpy(h_kmers.data(), tier.d_kmers, 
                      tier.num_kmers * sizeof(KmerEntry), 
                      cudaMemcpyDeviceToHost);
            
            tier_out.write(reinterpret_cast<const char*>(h_kmers.data()), 
                          tier.num_kmers * sizeof(KmerEntry));
        } else {
            // If not loaded, copy from existing file
            std::ifstream tier_in(tier.file_path, std::ios::binary);
            if (tier_in) {
                tier_out << tier_in.rdbuf();
            }
        }
        
        tier_out.close();
    }
    
    meta_out.close();
    
    std::cout << "Saved hierarchical database to: " << db_prefix << "\n";
    std::cout << "  Metadata: " << metadata_file << "\n";
    std::cout << "  " << tiers.size() << " tier files created\n";
}

// Load specific tier data into GPU memory
void HierarchicalGPUDatabase::load_tier_to_gpu(size_t tier_idx) {
    if (tier_idx >= tiers.size()) return;
    
    auto& tier = tiers[tier_idx];
    if (tier.gpu_loaded) return;  // Already loaded
    
    // Read tier data from file
    std::ifstream tier_file(tier.file_path, std::ios::binary);
    if (!tier_file) {
        throw std::runtime_error("Cannot open tier file: " + tier.file_path);
    }
    
    // Allocate host memory
    std::vector<KmerEntry> h_kmers(tier.num_kmers);
    
    // Read k-mer data
    tier_file.read(reinterpret_cast<char*>(h_kmers.data()), 
                   tier.num_kmers * sizeof(KmerEntry));
    tier_file.close();
    
    // Allocate GPU memory
    size_t tier_bytes = tier.num_kmers * sizeof(KmerEntry);
    cudaMalloc(&tier.d_kmers, tier_bytes);
    
    // Copy to GPU
    cudaMemcpy(tier.d_kmers, h_kmers.data(), tier_bytes, cudaMemcpyHostToDevice);
    
    tier.gpu_loaded = true;
    loaded_tier_count++;
    tier.last_access = std::chrono::steady_clock::now();
}

} // namespace biogpu