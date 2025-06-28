// compact_gpu_taxonomy.cu
// Implementation of compact GPU taxonomy data structures and processing

#include "compact_gpu_taxonomy.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <cstring>
#include <algorithm>
#include <queue>
#include <set>

namespace BioGPU {
namespace CompactTaxonomy {

// CUDA error checking macro
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                  << " - " << cudaGetErrorString(error) << std::endl; \
        return false; \
    } \
} while(0)

// Hash function for taxonomy lookup
__device__ inline uint32_t jenkins_hash(uint32_t key) {
    key = (key + 0x7ed55d16) + (key << 12);
    key = (key ^ 0xc761c23c) ^ (key >> 19);
    key = (key + 0x165667b1) + (key << 5);
    key = (key + 0xd3a2646c) ^ (key << 9);
    key = (key + 0xfd7046c5) + (key << 3);
    key = (key ^ 0xb55a4f09) ^ (key >> 16);
    return key;
}

// Host version of jenkins hash
__host__ inline uint32_t jenkins_hash_host(uint32_t key) {
    key = (key + 0x7ed55d16) + (key << 12);
    key = (key ^ 0xc761c23c) ^ (key >> 19);
    key = (key + 0x165667b1) + (key << 5);
    key = (key + 0xd3a2646c) ^ (key << 9);
    key = (key + 0xfd7046c5) + (key << 3);
    key = (key ^ 0xb55a4f09) ^ (key >> 16);
    return key;
}

// Device kernel for taxonomy lookup
__global__ void lookup_taxon_kernel(
    const TaxonHashEntry* entries,
    uint32_t table_size,
    uint32_t hash_mask,
    const uint32_t* query_taxons,
    uint32_t* parent_results,
    uint8_t* depth_results,
    uint8_t* rank_results,
    int num_queries
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;
    
    uint32_t taxon_id = query_taxons[tid];
    uint32_t hash = jenkins_hash(taxon_id) & hash_mask;
    
    // Linear probing
    for (int probe = 0; probe < 32; probe++) {
        uint32_t pos = (hash + probe) & hash_mask;
        const TaxonHashEntry& entry = entries[pos];
        
        if (entry.taxon_id == taxon_id) {
            parent_results[tid] = entry.parent_id;
            depth_results[tid] = entry.depth;
            rank_results[tid] = entry.rank_level;
            return;
        }
        
        if (entry.taxon_id == 0) {
            // Not found
            parent_results[tid] = 0;
            depth_results[tid] = 0;
            rank_results[tid] = 255;
            return;
        }
    }
    
    // Failed to find after max probes
    parent_results[tid] = 0;
    depth_results[tid] = 0;
    rank_results[tid] = 255;
}

// Device kernel for LCA computation
__global__ void compute_lca_kernel(
    const TaxonHashEntry* entries,
    uint32_t table_size,
    uint32_t hash_mask,
    const uint32_t* taxon1_array,
    const uint32_t* taxon2_array,
    uint32_t* lca_results,
    int num_pairs
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_pairs) return;
    
    uint32_t taxon1 = taxon1_array[tid];
    uint32_t taxon2 = taxon2_array[tid];
    
    if (taxon1 == taxon2) {
        lca_results[tid] = taxon1;
        return;
    }
    
    // Simple LCA algorithm using hash table lookups
    uint32_t current1 = taxon1;
    uint32_t current2 = taxon2;
    
    // Get depths
    uint8_t depth1 = 0, depth2 = 0;
    
    // Lookup depth for taxon1
    uint32_t hash = jenkins_hash(current1) & hash_mask;
    for (int probe = 0; probe < 32; probe++) {
        uint32_t pos = (hash + probe) & hash_mask;
        if (entries[pos].taxon_id == current1) {
            depth1 = entries[pos].depth;
            break;
        }
        if (entries[pos].taxon_id == 0) break;
    }
    
    // Lookup depth for taxon2
    hash = jenkins_hash(current2) & hash_mask;
    for (int probe = 0; probe < 32; probe++) {
        uint32_t pos = (hash + probe) & hash_mask;
        if (entries[pos].taxon_id == current2) {
            depth2 = entries[pos].depth;
            break;
        }
        if (entries[pos].taxon_id == 0) break;
    }
    
    // Bring both to same depth
    while (depth1 > depth2 && current1 != 1) {
        hash = jenkins_hash(current1) & hash_mask;
        for (int probe = 0; probe < 32; probe++) {
            uint32_t pos = (hash + probe) & hash_mask;
            if (entries[pos].taxon_id == current1) {
                current1 = entries[pos].parent_id;
                depth1--;
                break;
            }
        }
    }
    
    while (depth2 > depth1 && current2 != 1) {
        hash = jenkins_hash(current2) & hash_mask;
        for (int probe = 0; probe < 32; probe++) {
            uint32_t pos = (hash + probe) & hash_mask;
            if (entries[pos].taxon_id == current2) {
                current2 = entries[pos].parent_id;
                depth2--;
                break;
            }
        }
    }
    
    // Now walk up together
    while (current1 != current2 && current1 != 1 && current2 != 1) {
        // Move current1 up
        hash = jenkins_hash(current1) & hash_mask;
        for (int probe = 0; probe < 32; probe++) {
            uint32_t pos = (hash + probe) & hash_mask;
            if (entries[pos].taxon_id == current1) {
                current1 = entries[pos].parent_id;
                break;
            }
        }
        
        // Move current2 up
        hash = jenkins_hash(current2) & hash_mask;
        for (int probe = 0; probe < 32; probe++) {
            uint32_t pos = (hash + probe) & hash_mask;
            if (entries[pos].taxon_id == current2) {
                current2 = entries[pos].parent_id;
                break;
            }
        }
    }
    
    lca_results[tid] = (current1 == current2) ? current1 : 1;
}

bool CompactGPUTaxonomy::build_gpu_structures() {
    std::cout << "Building GPU hash table structures..." << std::endl;
    
    // Calculate hash table size (next power of 2 above 2x number of nodes)
    uint32_t table_size = 1;
    while (table_size < taxonomy_nodes.size() * 2) {
        table_size <<= 1;
    }
    uint32_t hash_mask = table_size - 1;
    
    // Create host-side hash table
    std::vector<TaxonHashEntry> host_entries(table_size);
    std::memset(host_entries.data(), 0, table_size * sizeof(TaxonHashEntry));
    
    // Build concatenated names string
    std::string all_names;
    std::unordered_map<uint32_t, uint16_t> name_offsets;
    
    for (const auto& [taxon_id, name] : taxon_names) {
        name_offsets[taxon_id] = static_cast<uint16_t>(all_names.size());
        all_names += name + '\0';
    }
    
    // Insert all nodes into hash table
    int collisions = 0;
    for (const auto& [taxon_id, node] : taxonomy_nodes) {
        uint32_t hash = jenkins_hash_host(taxon_id) & hash_mask;
        
        // Linear probing
        bool inserted = false;
        for (int probe = 0; probe < 32; probe++) {
            uint32_t pos = (hash + probe) & hash_mask;
            
            if (host_entries[pos].taxon_id == 0) {
                host_entries[pos].taxon_id = taxon_id;
                host_entries[pos].parent_id = node.parent_id;
                host_entries[pos].depth = node.depth_from_root;
                host_entries[pos].rank_level = node.rank_level;
                
                auto name_it = name_offsets.find(taxon_id);
                host_entries[pos].name_offset = (name_it != name_offsets.end()) ? 
                                                name_it->second : 0;
                inserted = true;
                break;
            }
            collisions++;
        }
        
        if (!inserted) {
            std::cerr << "Failed to insert taxon " << taxon_id << " into hash table" << std::endl;
            return false;
        }
    }
    
    std::cout << "Hash table built with " << collisions << " collisions" << std::endl;
    std::cout << "Load factor: " << (100.0 * taxonomy_nodes.size() / table_size) << "%" << std::endl;
    
    // Allocate GPU memory
    TaxonHashTable* h_hash_table = new TaxonHashTable;
    h_hash_table->table_size = table_size;
    h_hash_table->num_entries = taxonomy_nodes.size();
    h_hash_table->hash_mask = hash_mask;
    h_hash_table->names_data_size = all_names.size();
    
    // Allocate and copy hash entries
    CUDA_CHECK(cudaMalloc(&h_hash_table->entries, table_size * sizeof(TaxonHashEntry)));
    CUDA_CHECK(cudaMemcpy(h_hash_table->entries, host_entries.data(), 
                          table_size * sizeof(TaxonHashEntry), cudaMemcpyHostToDevice));
    
    // Allocate and copy names data
    if (!all_names.empty()) {
        CUDA_CHECK(cudaMalloc(&h_hash_table->names_data, all_names.size()));
        CUDA_CHECK(cudaMemcpy(h_hash_table->names_data, all_names.c_str(), 
                              all_names.size(), cudaMemcpyHostToDevice));
    }
    
    // Copy hash table structure to GPU
    CUDA_CHECK(cudaMalloc(&d_hash_table, sizeof(TaxonHashTable)));
    CUDA_CHECK(cudaMemcpy(d_hash_table, h_hash_table, sizeof(TaxonHashTable), cudaMemcpyHostToDevice));
    
    delete h_hash_table;
    
    // Build distance cache if enabled
    if (use_distance_cache) {
        if (!build_distance_cache()) {
            std::cerr << "Warning: Failed to build distance cache" << std::endl;
            // Non-fatal - continue without cache
        }
    }
    
    std::cout << "✓ GPU structures built successfully" << std::endl;
    return true;
}

bool CompactGPUTaxonomy::build_distance_cache() {
    std::cout << "Building phylogenetic distance cache..." << std::endl;
    
    // Select common taxon pairs for caching
    // For now, cache distances between major taxonomic groups
    std::vector<std::pair<uint32_t, uint32_t>> pairs_to_cache;
    
    // Find species, genera, families, etc.
    std::vector<uint32_t> species, genera, families;
    for (const auto& [taxon_id, rank] : taxon_ranks) {
        if (rank == "species" && species.size() < 1000) {
            species.push_back(taxon_id);
        } else if (rank == "genus" && genera.size() < 500) {
            genera.push_back(taxon_id);
        } else if (rank == "family" && families.size() < 200) {
            families.push_back(taxon_id);
        }
    }
    
    // Cache species-to-genus distances
    for (uint32_t sp : species) {
        for (uint32_t gen : genera) {
            pairs_to_cache.push_back({sp, gen});
            if (pairs_to_cache.size() >= 50000) break;
        }
        if (pairs_to_cache.size() >= 50000) break;
    }
    
    if (pairs_to_cache.empty()) {
        std::cout << "No pairs to cache" << std::endl;
        return true;
    }
    
    // Allocate cache structure
    PhyloDistanceCache* h_cache = new PhyloDistanceCache;
    h_cache->num_cached_pairs = pairs_to_cache.size();
    h_cache->max_pairs = pairs_to_cache.size();
    
    // Prepare data
    std::vector<uint32_t> taxon_pairs_flat;
    std::vector<uint8_t> distances;
    std::vector<uint32_t> lca_results;
    
    for (const auto& [t1, t2] : pairs_to_cache) {
        taxon_pairs_flat.push_back(t1);
        taxon_pairs_flat.push_back(t2);
        
        // Compute distance and LCA on host
        uint32_t lca = find_lca_pair(t1, t2);
        uint8_t dist = calculate_distance_sum(t1, t2, lca);
        
        distances.push_back(dist);
        lca_results.push_back(lca);
    }
    
    // Allocate GPU memory
    CUDA_CHECK(cudaMalloc(&h_cache->taxon_pairs, taxon_pairs_flat.size() * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&h_cache->distances, distances.size() * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&h_cache->lca_results, lca_results.size() * sizeof(uint32_t)));
    
    // Copy to GPU
    CUDA_CHECK(cudaMemcpy(h_cache->taxon_pairs, taxon_pairs_flat.data(), 
                          taxon_pairs_flat.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(h_cache->distances, distances.data(), 
                          distances.size() * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(h_cache->lca_results, lca_results.data(), 
                          lca_results.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Copy cache structure to GPU
    CUDA_CHECK(cudaMalloc(&d_distance_cache, sizeof(PhyloDistanceCache)));
    CUDA_CHECK(cudaMemcpy(d_distance_cache, h_cache, sizeof(PhyloDistanceCache), cudaMemcpyHostToDevice));
    
    delete h_cache;
    
    std::cout << "✓ Cached " << pairs_to_cache.size() << " taxon pair distances" << std::endl;
    return true;
}

uint32_t CompactGPUTaxonomy::find_lca_pair(uint32_t taxon1, uint32_t taxon2) {
    if (taxon1 == taxon2) return taxon1;
    
    // Get paths to root
    auto get_path = [this](uint32_t taxon) -> std::vector<uint32_t> {
        std::vector<uint32_t> path;
        uint32_t current = taxon;
        
        while (current != 1 && path.size() < 50) {
            path.push_back(current);
            auto it = taxonomy_nodes.find(current);
            if (it == taxonomy_nodes.end()) break;
            current = it->second.parent_id;
        }
        path.push_back(1);
        return path;
    };
    
    std::vector<uint32_t> path1 = get_path(taxon1);
    std::vector<uint32_t> path2 = get_path(taxon2);
    
    // Find first common ancestor
    std::set<uint32_t> ancestors1(path1.begin(), path1.end());
    
    for (uint32_t ancestor : path2) {
        if (ancestors1.find(ancestor) != ancestors1.end()) {
            return ancestor;
        }
    }
    
    return 1; // Root
}

uint8_t CompactGPUTaxonomy::calculate_distance_sum(uint32_t taxon1, uint32_t taxon2, uint32_t lca) {
    uint8_t dist1 = 0, dist2 = 0;
    
    // Distance from taxon1 to lca
    uint32_t current = taxon1;
    while (current != lca && current != 1 && dist1 < 50) {
        auto it = taxonomy_nodes.find(current);
        if (it == taxonomy_nodes.end()) break;
        current = it->second.parent_id;
        dist1++;
    }
    
    // Distance from taxon2 to lca
    current = taxon2;
    while (current != lca && current != 1 && dist2 < 50) {
        auto it = taxonomy_nodes.find(current);
        if (it == taxonomy_nodes.end()) break;
        current = it->second.parent_id;
        dist2++;
    }
    
    return dist1 + dist2;
}

void CompactGPUTaxonomy::free_gpu_memory() {
    if (d_hash_table) {
        // Get hash table from GPU to free internal pointers
        TaxonHashTable h_table;
        cudaMemcpy(&h_table, d_hash_table, sizeof(TaxonHashTable), cudaMemcpyDeviceToHost);
        
        if (h_table.entries) cudaFree(h_table.entries);
        if (h_table.names_data) cudaFree(h_table.names_data);
        
        cudaFree(d_hash_table);
        d_hash_table = nullptr;
    }
    
    if (d_distance_cache) {
        // Get cache from GPU to free internal pointers
        PhyloDistanceCache h_cache;
        cudaMemcpy(&h_cache, d_distance_cache, sizeof(PhyloDistanceCache), cudaMemcpyDeviceToHost);
        
        if (h_cache.taxon_pairs) cudaFree(h_cache.taxon_pairs);
        if (h_cache.distances) cudaFree(h_cache.distances);
        if (h_cache.lca_results) cudaFree(h_cache.lca_results);
        
        cudaFree(d_distance_cache);
        d_distance_cache = nullptr;
    }
}

// Host-side API for GPU taxonomy lookups
bool CompactGPUTaxonomy::lookup_taxons_gpu(
    const std::vector<uint32_t>& query_taxons,
    std::vector<uint32_t>& parent_results,
    std::vector<uint8_t>& depth_results,
    std::vector<uint8_t>& rank_results) {
    
    if (!d_hash_table || query_taxons.empty()) return false;
    
    // Get hash table info
    TaxonHashTable h_table;
    CUDA_CHECK(cudaMemcpy(&h_table, d_hash_table, sizeof(TaxonHashTable), cudaMemcpyDeviceToHost));
    
    // Allocate GPU memory for queries and results
    uint32_t* d_queries;
    uint32_t* d_parents;
    uint8_t* d_depths;
    uint8_t* d_ranks;
    
    size_t num_queries = query_taxons.size();
    
    CUDA_CHECK(cudaMalloc(&d_queries, num_queries * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_parents, num_queries * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_depths, num_queries * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&d_ranks, num_queries * sizeof(uint8_t)));
    
    // Copy queries to GPU
    CUDA_CHECK(cudaMemcpy(d_queries, query_taxons.data(), 
                          num_queries * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Launch kernel
    int threads = 256;
    int blocks = (num_queries + threads - 1) / threads;
    
    lookup_taxon_kernel<<<blocks, threads>>>(
        h_table.entries,
        h_table.table_size,
        h_table.hash_mask,
        d_queries,
        d_parents,
        d_depths,
        d_ranks,
        num_queries
    );
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Copy results back
    parent_results.resize(num_queries);
    depth_results.resize(num_queries);
    rank_results.resize(num_queries);
    
    CUDA_CHECK(cudaMemcpy(parent_results.data(), d_parents, 
                          num_queries * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(depth_results.data(), d_depths, 
                          num_queries * sizeof(uint8_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(rank_results.data(), d_ranks, 
                          num_queries * sizeof(uint8_t), cudaMemcpyDeviceToHost));
    
    // Cleanup
    cudaFree(d_queries);
    cudaFree(d_parents);
    cudaFree(d_depths);
    cudaFree(d_ranks);
    
    return true;
}

bool CompactGPUTaxonomy::compute_lca_batch_gpu(
    const std::vector<std::pair<uint32_t, uint32_t>>& taxon_pairs,
    std::vector<uint32_t>& lca_results) {
    
    if (!d_hash_table || taxon_pairs.empty()) return false;
    
    // Get hash table info
    TaxonHashTable h_table;
    CUDA_CHECK(cudaMemcpy(&h_table, d_hash_table, sizeof(TaxonHashTable), cudaMemcpyDeviceToHost));
    
    // Prepare data
    std::vector<uint32_t> taxon1_array, taxon2_array;
    for (const auto& [t1, t2] : taxon_pairs) {
        taxon1_array.push_back(t1);
        taxon2_array.push_back(t2);
    }
    
    // Allocate GPU memory
    uint32_t* d_taxon1;
    uint32_t* d_taxon2;
    uint32_t* d_lca_results;
    
    size_t num_pairs = taxon_pairs.size();
    
    CUDA_CHECK(cudaMalloc(&d_taxon1, num_pairs * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_taxon2, num_pairs * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_lca_results, num_pairs * sizeof(uint32_t)));
    
    // Copy to GPU
    CUDA_CHECK(cudaMemcpy(d_taxon1, taxon1_array.data(), 
                          num_pairs * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_taxon2, taxon2_array.data(), 
                          num_pairs * sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Launch kernel
    int threads = 256;
    int blocks = (num_pairs + threads - 1) / threads;
    
    compute_lca_kernel<<<blocks, threads>>>(
        h_table.entries,
        h_table.table_size,
        h_table.hash_mask,
        d_taxon1,
        d_taxon2,
        d_lca_results,
        num_pairs
    );
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Copy results back
    lca_results.resize(num_pairs);
    CUDA_CHECK(cudaMemcpy(lca_results.data(), d_lca_results, 
                          num_pairs * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    // Cleanup
    cudaFree(d_taxon1);
    cudaFree(d_taxon2);
    cudaFree(d_lca_results);
    
    return true;
}

} // namespace CompactTaxonomy
} // namespace BioGPU