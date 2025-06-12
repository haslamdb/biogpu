// gpu_kmer_database.h - GPU-optimized k-mer database for minimizer lookup
#ifndef GPU_KMER_DATABASE_H
#define GPU_KMER_DATABASE_H

#include <cuda_runtime.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

namespace biogpu {

// Database entry structure optimized for GPU access
struct KmerEntry {
    uint64_t hash;          // MurmurHash3 hash of k-mer
    uint32_t taxon_id;      // Taxonomic ID
    uint32_t reserved;      // Padding for alignment
};

// GPU-friendly hash table parameters
struct HashTableParams {
    uint64_t table_size;    // Size of hash table (power of 2)
    uint64_t bucket_size;   // Entries per bucket (for collision handling)
    uint64_t total_entries; // Total number of k-mers in database
};

class GPUKmerDatabase {
private:
    // Host data structures
    std::vector<KmerEntry> h_entries;
    std::unordered_map<uint64_t, uint32_t> h_hash_to_taxon;  // For building
    HashTableParams params;
    
    // Device data structures
    KmerEntry* d_hash_table = nullptr;
    HashTableParams* d_params = nullptr;
    uint32_t* d_taxon_counts = nullptr;  // For abundance tracking
    
    // Taxonomy data
    std::vector<uint32_t> parent_map;
    uint32_t* d_parent_map = nullptr;
    
    size_t allocated_entries = 0;
    size_t max_taxon_id = 0;
    
public:
    GPUKmerDatabase();
    ~GPUKmerDatabase();
    
    // Build database from k-mer file
    void build_from_file(const std::string& kmer_file, 
                        const std::string& taxonomy_file);
    
    // Add k-mers during construction
    void add_kmer(uint64_t hash, uint32_t taxon_id);
    
    // Finalize and transfer to GPU
    void finalize_and_upload();
    
    // Upload existing hash table to GPU (used when loading from binary)
    void upload_to_gpu();
    
    // GPU batch lookup - returns taxon IDs for each hash
    void lookup_batch_gpu(const uint64_t* d_hashes, 
                         uint32_t* d_taxon_ids, 
                         size_t count,
                         cudaStream_t stream = 0);
    
    // Get taxonomy information
    void load_taxonomy(const std::string& taxonomy_file);
    uint32_t get_parent(uint32_t taxon_id) const;
    
    // Statistics
    size_t get_num_kmers() const { return h_entries.size(); }
    size_t get_database_size_bytes() const;
    
    // Save/load optimized binary format
    void save_binary(const std::string& filename) const;
    void load_binary(const std::string& filename);
};

// CUDA kernels for database operations
__global__ void kmer_lookup_kernel(
    const uint64_t* query_hashes,
    const KmerEntry* hash_table,
    const HashTableParams* params,
    uint32_t* taxon_ids,
    size_t num_queries
);

__global__ void classify_reads_kernel(
    const uint32_t* minimizer_taxons,
    const uint32_t* minimizer_offsets,
    const uint32_t* minimizer_counts,
    const uint32_t* parent_map,
    uint32_t* read_classifications,
    float* confidence_scores,
    uint32_t num_reads
);

__global__ void update_abundance_kernel(
    const uint32_t* classifications,
    uint32_t* taxon_counts,
    uint32_t num_reads,
    uint32_t max_taxon_id
);

} // namespace biogpu

#endif // GPU_KMER_DATABASE_H