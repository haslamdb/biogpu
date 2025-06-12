// gpu_kmer_database.cu - Implementation of GPU-optimized k-mer database
#include "gpu_kmer_database.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cuda_runtime.h>

namespace biogpu {

// Hash function for GPU - maps hash to bucket index
__device__ inline uint64_t hash_to_bucket(uint64_t hash, uint64_t table_size) {
    // Use the hash directly since it's already well-distributed by MurmurHash3
    return hash & (table_size - 1);  // Assumes table_size is power of 2
}

// GPU kernel for k-mer lookup using open addressing with linear probing
__global__ void kmer_lookup_kernel(
    const uint64_t* query_hashes,
    const KmerEntry* hash_table,
    const HashTableParams* params,
    uint32_t* taxon_ids,
    size_t num_queries
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;
    
    const uint64_t query_hash = query_hashes[tid];
    const uint64_t bucket = hash_to_bucket(query_hash, params->table_size);
    const uint64_t bucket_size = params->bucket_size;
    
    // Linear probing within bucket
    uint32_t result_taxon = 0;
    
    for (uint64_t i = 0; i < bucket_size; i++) {
        uint64_t idx = bucket * bucket_size + i;
        const KmerEntry& entry = hash_table[idx];
        
        if (entry.hash == query_hash) {
            result_taxon = entry.taxon_id;
            break;
        } else if (entry.hash == 0) {
            // Empty slot, k-mer not found
            break;
        }
    }
    
    taxon_ids[tid] = result_taxon;
}

// Classify reads based on minimizer hits using LCA approach
__global__ void classify_reads_kernel(
    const uint32_t* minimizer_taxons,
    const uint32_t* minimizer_offsets,
    const uint32_t* minimizer_counts,
    const uint32_t* parent_map,
    uint32_t* read_classifications,
    float* confidence_scores,
    uint32_t num_reads
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const uint32_t start = minimizer_offsets[tid];
    const uint32_t count = minimizer_counts[tid];
    
    if (count == 0) {
        read_classifications[tid] = 0;
        confidence_scores[tid] = 0.0f;
        return;
    }
    
    // Count hits per taxon
    const int MAX_TAXA = 32;  // Local array size limit
    uint32_t taxa[MAX_TAXA];
    uint32_t counts[MAX_TAXA];
    int num_taxa = 0;
    
    for (uint32_t i = 0; i < count && i < 100; i++) {  // Limit to first 100 minimizers
        uint32_t taxon = minimizer_taxons[start + i];
        if (taxon == 0) continue;
        
        // Find or add taxon
        int idx = -1;
        for (int j = 0; j < num_taxa; j++) {
            if (taxa[j] == taxon) {
                idx = j;
                break;
            }
        }
        
        if (idx >= 0) {
            counts[idx]++;
        } else if (num_taxa < MAX_TAXA) {
            taxa[num_taxa] = taxon;
            counts[num_taxa] = 1;
            num_taxa++;
        }
    }
    
    // Find taxon with most hits
    uint32_t best_taxon = 0;
    uint32_t max_hits = 0;
    uint32_t total_hits = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        total_hits += counts[i];
        if (counts[i] > max_hits) {
            max_hits = counts[i];
            best_taxon = taxa[i];
        }
    }
    
    read_classifications[tid] = best_taxon;
    confidence_scores[tid] = total_hits > 0 ? (float)max_hits / total_hits : 0.0f;
}

// Update abundance counts
__global__ void update_abundance_kernel(
    const uint32_t* classifications,
    uint32_t* taxon_counts,
    uint32_t num_reads,
    uint32_t max_taxon_id
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    uint32_t taxon = classifications[tid];
    if (taxon > 0 && taxon <= max_taxon_id) {
        atomicAdd(&taxon_counts[taxon], 1);
    }
}

// Host implementation
GPUKmerDatabase::GPUKmerDatabase() {
    params.table_size = 1ULL << 20;  // Start with 1M buckets
    params.bucket_size = 8;          // 8 entries per bucket
    params.total_entries = 0;
}

GPUKmerDatabase::~GPUKmerDatabase() {
    if (d_hash_table) cudaFree(d_hash_table);
    if (d_params) cudaFree(d_params);
    if (d_taxon_counts) cudaFree(d_taxon_counts);
    if (d_parent_map) cudaFree(d_parent_map);
}

void GPUKmerDatabase::add_kmer(uint64_t hash, uint32_t taxon_id) {
    h_hash_to_taxon[hash] = taxon_id;
    max_taxon_id = std::max(max_taxon_id, (size_t)taxon_id);
}

void GPUKmerDatabase::finalize_and_upload() {
    // Determine optimal table size (2x entries for 50% load factor)
    size_t num_entries = h_hash_to_taxon.size();
    params.total_entries = num_entries;
    
    // Find next power of 2
    size_t table_size = 1;
    while (table_size < num_entries * 2) {
        table_size <<= 1;
    }
    params.table_size = table_size;
    
    std::cout << "Building GPU hash table with " << table_size 
              << " buckets for " << num_entries << " k-mers\n";
    
    // Build hash table with linear probing
    size_t total_slots = params.table_size * params.bucket_size;
    h_entries.clear();
    h_entries.resize(total_slots, {0, 0, 0});
    
    size_t collisions = 0;
    for (const auto& [hash, taxon] : h_hash_to_taxon) {
        uint64_t bucket = hash & (params.table_size - 1);
        bool inserted = false;
        
        for (uint64_t i = 0; i < params.bucket_size; i++) {
            uint64_t idx = bucket * params.bucket_size + i;
            if (h_entries[idx].hash == 0) {
                h_entries[idx].hash = hash;
                h_entries[idx].taxon_id = taxon;
                inserted = true;
                if (i > 0) collisions++;
                break;
            }
        }
        
        if (!inserted) {
            std::cerr << "Warning: Bucket overflow, k-mer dropped\n";
        }
    }
    
    std::cout << "Hash table built with " << collisions 
              << " collisions (" << (100.0 * collisions / num_entries) << "%)\n";
    
    // Allocate GPU memory
    size_t table_bytes = total_slots * sizeof(KmerEntry);
    if (table_bytes > allocated_entries) {
        if (d_hash_table) cudaFree(d_hash_table);
        cudaMalloc(&d_hash_table, table_bytes);
        allocated_entries = table_bytes;
    }
    
    // Upload to GPU
    cudaMemcpy(d_hash_table, h_entries.data(), table_bytes, cudaMemcpyHostToDevice);
    
    if (!d_params) {
        cudaMalloc(&d_params, sizeof(HashTableParams));
    }
    cudaMemcpy(d_params, &params, sizeof(HashTableParams), cudaMemcpyHostToDevice);
    
    // Allocate taxon counts
    if (!d_taxon_counts) {
        cudaMalloc(&d_taxon_counts, (max_taxon_id + 1) * sizeof(uint32_t));
        cudaMemset(d_taxon_counts, 0, (max_taxon_id + 1) * sizeof(uint32_t));
    }
    
    std::cout << "Database uploaded to GPU (" 
              << table_bytes / (1024.0 * 1024.0) << " MB)\n";
}

void GPUKmerDatabase::lookup_batch_gpu(const uint64_t* d_hashes, 
                                      uint32_t* d_taxon_ids, 
                                      size_t count,
                                      cudaStream_t stream) {
    int block_size = 256;
    int num_blocks = (count + block_size - 1) / block_size;
    
    kmer_lookup_kernel<<<num_blocks, block_size, 0, stream>>>(
        d_hashes, d_hash_table, d_params, d_taxon_ids, count
    );
}

size_t GPUKmerDatabase::get_database_size_bytes() const {
    return params.table_size * params.bucket_size * sizeof(KmerEntry);
}

void GPUKmerDatabase::save_binary(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Write header
    file.write(reinterpret_cast<const char*>(&params), sizeof(params));
    file.write(reinterpret_cast<const char*>(&max_taxon_id), sizeof(max_taxon_id));
    
    // Write hash table
    size_t table_entries = params.table_size * params.bucket_size;
    file.write(reinterpret_cast<const char*>(h_entries.data()), 
               table_entries * sizeof(KmerEntry));
    
    std::cout << "Saved database to " << filename << "\n";
}

void GPUKmerDatabase::load_binary(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }
    
    // Read header
    file.read(reinterpret_cast<char*>(&params), sizeof(params));
    file.read(reinterpret_cast<char*>(&max_taxon_id), sizeof(max_taxon_id));
    
    // Read hash table
    size_t table_entries = params.table_size * params.bucket_size;
    h_entries.resize(table_entries);
    file.read(reinterpret_cast<char*>(h_entries.data()), 
              table_entries * sizeof(KmerEntry));
    
    std::cout << "Loaded database from " << filename << "\n";
    std::cout << "  Table size: " << params.table_size << " buckets\n";
    std::cout << "  Total k-mers: " << params.total_entries << "\n";
    
    // Upload to GPU
    finalize_and_upload();
}

} // namespace biogpu