// gpu_kraken_classifier.h
// Header file for GPU-accelerated Kraken2-style taxonomic classifier with PAIRED-END support
// Contains only declarations - implementations in gpu_kraken_classifier.cu

#ifndef GPU_KRAKEN_CLASSIFIER_H
#define GPU_KRAKEN_CLASSIFIER_H

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>

// Forward declarations and structures
struct PairedRead {
    std::string read1;
    std::string read2;
    std::string read_id;
    bool is_paired;
    
    PairedRead(const std::string& r1, const std::string& r2 = "", 
               const std::string& id = "") 
        : read1(r1), read2(r2), read_id(id), is_paired(!r2.empty()) {}
};

struct GPUCompactHashTable {
    uint32_t* hash_cells;
    uint32_t table_size;
    uint32_t hash_mask;
    uint32_t lca_bits;
    uint32_t hash_bits;
};

struct ClassificationParams {
    int k = 35;
    int ell = 31;
    int spaces = 7;
    float confidence_threshold = 0.0f;
    bool use_spaced_seeds = true;
    int max_ambiguous_bases = 5;
    bool use_paired_end_bonus = true;
    float paired_concordance_weight = 2.0f;
    float min_pair_concordance = 0.5f;
    bool require_both_reads_classified = false;
};

struct PairedReadClassification {
    uint32_t taxon_id;
    float confidence_score;
    uint32_t read1_votes;
    uint32_t read2_votes;
    uint32_t read1_kmers;
    uint32_t read2_kmers;
    uint32_t concordant_votes;
    float pair_concordance;
    bool got_paired_bonus;
    uint32_t read1_best_taxon;
    uint32_t read2_best_taxon;
    float read1_confidence;
    float read2_confidence;
};

struct TaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    char name[64];
};

// CUDA kernel declarations (implemented in .cu file)
__global__ void classify_paired_reads_kernel(
    const char* reads_r1,
    const char* reads_r2,
    const uint32_t* read_offsets_r1,
    const uint32_t* read_offsets_r2,
    const uint32_t* read_lengths_r1,
    const uint32_t* read_lengths_r2,
    const bool* is_paired_flags,
    const GPUCompactHashTable* hash_table,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup,
    uint32_t* pair_votes_r1,
    uint32_t* pair_votes_r2,
    uint32_t* concordant_votes,
    PairedReadClassification* results,
    int num_pairs,
    ClassificationParams params
);

__global__ void compute_paired_concordance_kernel(
    const uint32_t* pair_votes_r1,
    const uint32_t* pair_votes_r2,
    const bool* is_paired_flags,
    uint32_t* concordant_votes,
    PairedReadClassification* results,
    int num_pairs,
    ClassificationParams params
);

// Device function declarations
__device__ uint32_t lookup_lca_gpu(const GPUCompactHashTable* cht, uint64_t minimizer_hash);
__device__ uint64_t extract_minimizer_with_spaced_seeds(const char* sequence, int pos, int k, int ell, int spaces);
__device__ uint32_t find_paired_classification_with_confidence(
    const uint32_t* votes_r1,
    const uint32_t* votes_r2,
    const uint32_t* concordant_votes,
    int total_kmers_r1,
    int total_kmers_r2,
    float confidence_threshold,
    float paired_weight,
    bool use_paired_bonus,
    const TaxonomyNode* taxonomy_tree,
    const uint32_t* parent_lookup
);

// Utility function declarations (implementations in .cu file)
__device__ __host__ uint64_t hash_minimizer(const char* seq, int len);
__device__ __host__ uint64_t apply_spaced_seed_mask(uint64_t hash, int spaces);
__device__ __host__ bool has_ambiguous_bases(const char* seq, int len);

// Inline implementations for cross-compilation unit usage
__device__ __host__ inline uint32_t jenkins_hash(uint64_t key) {
    uint32_t hash = (uint32_t)(key ^ (key >> 32));
    hash += (hash << 10);
    hash ^= (hash >> 6);
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

__device__ __host__ inline uint32_t compute_compact_hash(uint64_t minimizer_hash) {
    return jenkins_hash(minimizer_hash) & 0x7FFFFFFF;
}

__device__ inline uint32_t lookup_lca_gpu_impl(const GPUCompactHashTable* cht, uint64_t minimizer_hash) {
    uint32_t compact_hash = compute_compact_hash(minimizer_hash);
    uint32_t pos = compact_hash & cht->hash_mask;
    uint32_t lca_mask = (1U << cht->lca_bits) - 1;
    
    for (int probe = 0; probe < 32; probe++) {
        uint32_t cell = cht->hash_cells[pos];
        if (cell == 0) return 0;
        
        uint32_t stored_hash = cell >> cht->lca_bits;
        uint32_t expected_hash = compact_hash >> cht->lca_bits;
        
        if (stored_hash == expected_hash) {
            return cell & lca_mask;
        }
        
        pos = (pos + 1) & cht->hash_mask;
    }
    
    return 0;
}

// Main classifier class - DECLARATION ONLY
class PairedEndGPUKrakenClassifier {
private:
    // GPU data structures
    GPUCompactHashTable* d_hash_table;
    TaxonomyNode* d_taxonomy_tree;
    uint32_t* d_parent_lookup;
    uint32_t num_taxonomy_nodes;
    
    // Host taxonomy mapping
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    
    // Classification parameters
    ClassificationParams params;
    
    // GPU memory for paired-end batch processing
    char* d_reads_r1;
    char* d_reads_r2;
    uint32_t* d_read_offsets_r1;
    uint32_t* d_read_offsets_r2;
    uint32_t* d_read_lengths_r1;
    uint32_t* d_read_lengths_r2;
    bool* d_is_paired_flags;
    PairedReadClassification* d_results;
    
    // Paired-end voting arrays
    uint32_t* d_pair_votes_r1;
    uint32_t* d_pair_votes_r2;
    uint32_t* d_concordant_votes;
    
    // Configuration constants
    static const int DEFAULT_BATCH_SIZE = 10000;
    static const int MAX_READ_LENGTH = 1000;
    static const int THREADS_PER_BLOCK = 256;
    static const int MAX_TAXA_PER_PAIR = 128;
    
    bool database_loaded = false;
    
public:
    // Constructor and destructor - DECLARATIONS ONLY
    PairedEndGPUKrakenClassifier(const ClassificationParams& config = ClassificationParams());
    ~PairedEndGPUKrakenClassifier();
    
    // Public interface methods - DECLARATIONS ONLY
    bool load_database(const std::string& database_directory);
    
    std::vector<PairedReadClassification> classify_paired_reads(
        const std::vector<PairedRead>& paired_reads
    );
    
    std::vector<PairedReadClassification> classify_paired_reads_batch(
        const std::vector<PairedRead>& paired_reads,
        int batch_size = DEFAULT_BATCH_SIZE
    );
    
    std::vector<PairedReadClassification> classify_reads(
        const std::vector<std::string>& reads
    );
    
    // Configuration methods
    void set_confidence_threshold(float threshold);
    void set_paired_concordance_weight(float weight);
    void enable_paired_end_bonus(bool enable);
    ClassificationParams get_params() const;
    
    // Utility methods
    std::string get_taxon_name(uint32_t taxon_id) const;
    void print_database_stats() const;
    void print_paired_classification_stats(const std::vector<PairedReadClassification>& results) const;
    bool is_database_loaded() const;
    
private:
    // Private method declarations
    bool allocate_gpu_memory(int max_pairs);
    void free_gpu_memory();
    bool load_hash_table(const std::string& hash_table_file);
    bool load_taxonomy_tree(const std::string& taxonomy_file);
    bool load_taxonomy_tree_from_db(const std::string& db_file);
    bool read_database_header(const std::string& db_file, ClassificationParams& db_params);
    void transfer_paired_reads_to_gpu(const std::vector<PairedRead>& paired_reads);
    void retrieve_paired_results_from_gpu(std::vector<PairedReadClassification>& results, int num_pairs);
};

#endif // GPU_KRAKEN_CLASSIFIER_H