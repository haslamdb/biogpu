#ifndef GPU_KRAKEN_DATABASE_BUILDER_H
#define GPU_KRAKEN_DATABASE_BUILDER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cuda_runtime.h>

// Core type definitions
#include "gpu_kraken_types.h"
#include "gpu_kraken_classifier.h"
#include "gpu_minimizer_extraction.cuh"

// Forward declarations
struct GPUGenomeInfo;
struct GPUMinimizerHit;
struct GPUTaxonomyNode;
struct LCACandidate;


class GPUKrakenDatabaseBuilder {
private:
    // Configuration
    ClassificationParams params;
    std::string output_directory;
    MinimizerParams minimizer_params;

    // Host data
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    std::vector<LCACandidate> all_lca_candidates;

    // GPU Memory
    char* d_sequence_data;
    GPUGenomeInfo* d_genome_info;
    GPUMinimizerHit* d_minimizer_hits;
    uint32_t* d_minimizer_counts;
    uint32_t* d_global_hit_counter;
    LCACandidate* d_lca_candidates;
    
    // Memory and Batching Controls
    int MAX_SEQUENCE_BATCH = 25;
    int MAX_MINIMIZERS_PER_BATCH = 5000000;
    bool auto_scale_memory = true;
    size_t max_gpu_memory_usage_fraction = 80;

    // Statistics
    GPUBuildStats stats;
    uint64_t total_minimizers_capped = 0;
    uint64_t total_batches_capped = 0;

public:
    // Constructor and Destructor
    GPUKrakenDatabaseBuilder(const std::string& output_dir, const ClassificationParams& config = ClassificationParams());
    ~GPUKrakenDatabaseBuilder();

    // **CRITICAL** Call this right after constructor
    bool initialize_cuda_context();

    // Main build pipeline
    bool build_database_from_genomes(const std::string& genome_library_path, const std::string& taxonomy_path = "");
    
    // Configuration
    void set_batch_size(int sequences_per_batch);
    void set_minimizer_capacity(int minimizers_per_batch);
    void enable_auto_memory_scaling(bool enable = true, size_t memory_fraction = 80);

    // Public Utilities
    void print_build_statistics();

private:
    // Internal pipeline steps
    bool load_genome_files(const std::string& library_path);
    bool load_taxonomy_data(const std::string& taxonomy_path);
    bool process_genomes_gpu();
    bool save_database();

    // GPU Memory Management
    bool allocate_gpu_memory();
    void free_gpu_memory();
    void check_and_adjust_memory_safe();
    void validate_memory_settings();

    // GPU Processing
    bool process_sequence_batch(const std::vector<std::string>& sequences, const std::vector<uint32_t>& taxon_ids, int batch_offset);
    bool extract_minimizers_gpu(const char* d_sequences, const GPUGenomeInfo* d_genomes, int num_sequences, GPUMinimizerHit* d_hits, uint32_t* total_hits);
    bool compute_lca_assignments_gpu(const GPUMinimizerHit* d_hits, int num_hits, LCACandidate* d_candidates, int* num_candidates);

    // Host Helpers
    std::vector<std::string> load_sequences_from_fasta(const std::string& fasta_path);
    uint32_t extract_taxon_from_filename(const std::string& filename);
    uint32_t compute_simple_lca_host(uint32_t taxon1, uint32_t taxon2);
};

#endif // GPU_KRAKEN_DATABASE_BUILDER_H