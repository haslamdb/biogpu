// gpu_kraken_database_builder.h
// Header file for GPU-accelerated Kraken2-style database builder

#ifndef GPU_KRAKEN_DATABASE_BUILDER_H
#define GPU_KRAKEN_DATABASE_BUILDER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cuda_runtime.h>

// Include necessary headers for types
#include "gpu_kraken_types.h"  // Common type definitions
#include "gpu_kraken_classifier.h"  // For ClassificationParams
#include "gpu_minimizer_extraction.cuh"  // For MinimizerParams

// Forward declarations for types used but not defined here
struct GPUGenomeInfo;
struct GPUMinimizerHit;
struct GPUTaxonomyNode;
class ConcatenatedFnaProcessor;

// Main database builder class declaration
class GPUKrakenDatabaseBuilder {
private:
    // Configuration
    ClassificationParams params;
    std::string output_directory;
    
    // Host data
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    std::vector<LCACandidate> all_lca_candidates;
    
    // GPU memory management
    char* d_sequence_data;
    GPUGenomeInfo* d_genome_info;
    GPUMinimizerHit* d_minimizer_hits;
    uint32_t* d_minimizer_counts;
    LCACandidate* d_lca_candidates;
    GPUTaxonomyNode* d_taxonomy_nodes;
    
    // UPDATED: Higher minimizer capacity for comprehensive databases
    int MAX_SEQUENCE_BATCH = 25;           // Keep sequences reasonable
    int MAX_MINIMIZERS_PER_BATCH = 5000000; // INCREASED from 1M to 5M
    static const int MAX_READ_LENGTH = 1000;
    static const int THREADS_PER_BLOCK = 256;
    static const int MAX_TAXA_PER_PAIR = 128;
    
    // NEW: Dynamic memory scaling based on available GPU memory
    bool auto_scale_memory = true;
    size_t max_gpu_memory_usage_fraction = 80; // Use up to 80% of GPU memory
    
    // NEW: Statistics to track if we're hitting limits
    uint64_t total_minimizers_capped = 0;
    uint64_t total_batches_capped = 0;
    
    // Convert existing params to minimizer params
    MinimizerParams minimizer_params;
    
    // NEW: Kraken2-inspired parameters
    uint64_t min_clear_hash_value = 0;  // For hash-based subsampling
    double subsampling_rate = 1.0;      // Fraction of minimizers to keep
    uint64_t toggle_mask = 0xe37e28c4271b5a2dULL;  // Kraken2-style toggle mask
    
    // Statistics
    GPUBuildStats stats;
    
    // NEW: Enhanced taxonomy processor and phylogenetic data
    EnhancedNCBITaxonomyProcessor enhanced_taxonomy_processor;
    EnhancedBuildStats enhanced_stats;
    std::vector<PhylogeneticLCACandidate> phylogenetic_candidates;
    ContributingTaxaArrays contributing_taxa_arrays;
    SpeciesTrackingData species_tracking;
    
public:
    GPUKrakenDatabaseBuilder(const std::string& output_dir,
                            const ClassificationParams& config = ClassificationParams());
    ~GPUKrakenDatabaseBuilder();
    
    // FIXED: Separate CUDA initialization from constructor to avoid heap corruption
    bool initialize_cuda_context();
    
    // Main pipeline methods (keep all original interfaces)
    bool build_database_from_genomes(
        const std::string& genome_library_path,
        const std::string& taxonomy_path = ""
    );
    
    bool build_database_from_file_list(
        const std::string& file_list_path,
        const std::string& taxonomy_path = ""
    );
    
    // Step-by-step processing
    bool load_genome_files(const std::string& library_path);
    bool load_taxonomy_data(const std::string& taxonomy_path);
    bool process_genomes_gpu();
    bool save_database();
    
    // Configuration
    void set_batch_size(int sequences_per_batch);
    
    // NEW: Method to set minimizer capacity explicitly
    void set_minimizer_capacity(int minimizers_per_batch);
    
    // NEW: Auto-scale based on GPU memory
    void enable_auto_memory_scaling(bool enable = true, size_t memory_fraction = 80);
    
    // NEW: Set subsampling rate (like Kraken2's --max-db-size functionality)
    void set_subsampling_rate(double rate);
    
    void print_build_statistics();
    
    // Testing function for integration
    void test_minimizer_extraction_integration();
    
    // NEW: Enhanced database building methods
    bool build_database_from_concatenated_fna(
        const std::string& fna_file_path,
        const std::string& taxonomy_nodes_path = "",
        const std::string& taxonomy_names_path = "",
        const std::string& compact_taxonomy_path = "");
    
    // NEW: Alternative method that uses pre-parsed FNA processor
    bool build_database_with_fna_processor(
        ConcatenatedFnaProcessor& fna_processor,
        const std::string& taxonomy_nodes_path = "",
        const std::string& taxonomy_names_path = "",
        const std::string& compact_taxonomy_path = "");
    
    // NEW: Pre-build compact taxonomy
    bool prebuild_compact_taxonomy(const std::string& nodes_path, 
                                  const std::string& names_path,
                                  const std::string& output_path);
    
    // NEW: Streaming support for large concatenated FNA files
    bool build_database_from_streaming_fna(const std::string& fna_file_path,
                                          const std::string& taxonomy_path = "");
    
private:
    // GPU processing methods
    bool allocate_gpu_memory();
    void free_gpu_memory();
    void check_and_adjust_memory();
    void check_and_adjust_memory_safe();
    
    bool process_sequence_batch(
        const std::vector<std::string>& sequences,
        const std::vector<uint32_t>& taxon_ids,
        int batch_offset
    );
    
    // UPDATED: Use improved minimizer extraction
    bool extract_minimizers_gpu(
        const char* d_sequences,
        const GPUGenomeInfo* d_genomes,
        int num_sequences,
        GPUMinimizerHit* d_hits,
        uint32_t* total_hits
    );
    
    // FIXED: Proper LCA computation without over-aggressive deduplication
    bool compute_lca_assignments_gpu(
        const GPUMinimizerHit* d_hits,
        int num_hits,
        LCACandidate* d_candidates,
        int* num_candidates
    );
    
    bool build_compact_hash_table_gpu(
        const LCACandidate* d_candidates,
        int num_candidates
    );
    
    // File processing (keep all original methods)
    std::vector<std::string> load_sequences_from_fasta(const std::string& fasta_path);
    bool parse_taxonomy_files(const std::string& taxonomy_path);
    
    // Utility methods
    uint32_t extract_taxon_from_filename(const std::string& filename);
    std::vector<std::string> find_genome_files(const std::string& directory);
    
    // NEW: Enhanced processing methods
    bool process_concatenated_fna_file(const std::string& fna_file_path);
    bool enhance_candidates_with_phylogenetics();
    bool save_enhanced_database();
    
    // NEW: Helper methods for saving enhanced database
    void build_contributing_taxa_arrays();
    uint32_t add_to_contributing_taxa_arrays(const PhylogeneticLCACandidate& candidate);
    bool save_streamlined_hash_table(const std::vector<StreamlinedMinimizerMetadata>& metadata, 
                                      const std::string& filename);
    bool save_contributing_taxa_arrays(const std::string& filename);
};

#endif // GPU_KRAKEN_DATABASE_BUILDER_H