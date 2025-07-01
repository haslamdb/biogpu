// gpu_kraken_types.h
// Common type definitions for GPU Kraken database builder
// Shared across all modules to prevent duplicate definitions

#ifndef GPU_KRAKEN_TYPES_H
#define GPU_KRAKEN_TYPES_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <fstream>
#include <iostream>
#include <set>
#include <chrono>

// ===========================
// Forward Declarations
// ===========================

struct MinimizerParams;
struct GPUTaxonomyNode;
struct SpeciesTrackingData;
struct PhylogeneticLCACandidate;
struct ContributingTaxaArrays;

// ===========================
// CUDA and GPU Types
// ===========================

// GPU genome information structure
struct GPUGenomeInfo {
    uint32_t genome_id;
    uint32_t sequence_offset;        // Use this consistently
    uint32_t sequence_length;
    uint32_t minimizer_count;
    uint32_t taxon_id;
};

// GPU minimizer hit structure
struct GPUMinimizerHit {
    uint64_t minimizer_hash;  // 8 bytes
    uint32_t genome_id;       // 4 bytes
    uint32_t position;        // 4 bytes
    uint16_t strand;          // 2 bytes - bits 0-1: strand (0=forward, 1=reverse)
                              //          bits 2-3: classification (0=unique, 1=canonical, 2=redundant)
                              //          bits 4-15: reserved for future use
    uint16_t taxon_id;        // 2 bytes
    uint16_t ml_weight;       // 2 bytes - ML confidence score (1.0 scaled to uint16_t)
    uint16_t feature_flags;   // 2 bytes - encoded features
}; // Total: 24 bytes

// Strand and classification flag constants
namespace MinimizerFlags {
    // Strand flags (bits 0-1)
    constexpr uint16_t STRAND_MASK = 0x0003;
    constexpr uint16_t STRAND_FORWARD = 0x0000;
    constexpr uint16_t STRAND_REVERSE = 0x0001;
    
    // Classification flags (bits 2-3)
    constexpr uint16_t CLASSIFICATION_MASK = 0x000C;
    constexpr uint16_t CLASSIFICATION_SHIFT = 2;
    constexpr uint16_t CLASSIFICATION_UNIQUE = 0x0000;      // Unique to one species
    constexpr uint16_t CLASSIFICATION_CANONICAL = 0x0004;   // Canonical for species (most common)
    constexpr uint16_t CLASSIFICATION_REDUNDANT = 0x0008;   // Redundant within species
    
    // Helper functions
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t get_strand(uint16_t flags) {
        return flags & STRAND_MASK;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t get_classification(uint16_t flags) {
        return (flags & CLASSIFICATION_MASK) >> CLASSIFICATION_SHIFT;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t set_classification(uint16_t flags, uint16_t classification) {
        return (flags & ~CLASSIFICATION_MASK) | ((classification << CLASSIFICATION_SHIFT) & CLASSIFICATION_MASK);
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline bool is_unique(uint16_t flags) {
        return get_classification(flags) == 0;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline bool is_canonical(uint16_t flags) {
        return get_classification(flags) == 1;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline bool is_redundant(uint16_t flags) {
        return get_classification(flags) == 2;
    }
    
    // Feature flags for ml_weight and feature_flags fields
    // For feature_flags field:
    // Bits 0-2: GC content category (0-7 representing GC% ranges)
    constexpr uint16_t GC_CONTENT_MASK = 0x0007;
    constexpr uint16_t GC_CONTENT_SHIFT = 0;
    
    // Bits 3-5: Sequence complexity score (0-7)
    constexpr uint16_t COMPLEXITY_MASK = 0x0038;
    constexpr uint16_t COMPLEXITY_SHIFT = 3;
    
    // Bit 6: Position bias indicator (clustered=1, uniform=0)
    constexpr uint16_t POSITION_BIAS_MASK = 0x0040;
    constexpr uint16_t POSITION_BIAS_CLUSTERED = 0x0040;
    constexpr uint16_t POSITION_BIAS_UNIFORM = 0x0000;
    
    // Bit 7: Contamination risk flag
    constexpr uint16_t CONTAMINATION_RISK_MASK = 0x0080;
    constexpr uint16_t CONTAMINATION_RISK_FLAG = 0x0080;
    
    // Bits 8-15: Reserved for future use
    constexpr uint16_t RESERVED_MASK = 0xFF00;
    
    // Helper functions for feature flags
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t get_gc_content_category(uint16_t feature_flags) {
        return (feature_flags & GC_CONTENT_MASK) >> GC_CONTENT_SHIFT;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t set_gc_content_category(uint16_t feature_flags, uint16_t category) {
        return (feature_flags & ~GC_CONTENT_MASK) | ((category << GC_CONTENT_SHIFT) & GC_CONTENT_MASK);
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t get_complexity_score(uint16_t feature_flags) {
        return (feature_flags & COMPLEXITY_MASK) >> COMPLEXITY_SHIFT;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t set_complexity_score(uint16_t feature_flags, uint16_t score) {
        return (feature_flags & ~COMPLEXITY_MASK) | ((score << COMPLEXITY_SHIFT) & COMPLEXITY_MASK);
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline bool has_position_bias(uint16_t feature_flags) {
        return (feature_flags & POSITION_BIAS_MASK) != 0;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t set_position_bias(uint16_t feature_flags, bool clustered) {
        if (clustered) {
            return feature_flags | POSITION_BIAS_CLUSTERED;
        } else {
            return feature_flags & ~POSITION_BIAS_MASK;
        }
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline bool has_contamination_risk(uint16_t feature_flags) {
        return (feature_flags & CONTAMINATION_RISK_MASK) != 0;
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t set_contamination_risk(uint16_t feature_flags, bool risk) {
        if (risk) {
            return feature_flags | CONTAMINATION_RISK_FLAG;
        } else {
            return feature_flags & ~CONTAMINATION_RISK_MASK;
        }
    }
    
    // ML weight conversion helpers (for ml_weight field)
    // Convert float confidence (0.0-1.0) to uint16_t (0-65535)
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline uint16_t float_to_ml_weight(float confidence) {
        if (confidence <= 0.0f) return 0;
        if (confidence >= 1.0f) return 65535;
        return static_cast<uint16_t>(confidence * 65535.0f);
    }
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    inline float ml_weight_to_float(uint16_t weight) {
        return static_cast<float>(weight) / 65535.0f;
    }
}

// CUDA-compatible minimizer parameters
struct MinimizerParams {
    uint32_t k;                 // K-mer size
    uint32_t ell;               // Minimizer length
    uint32_t spaces;            // Spaced k-mer pattern
    uint64_t xor_mask;          // XOR mask for hashing
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    MinimizerParams() : k(31), ell(31), spaces(7), xor_mask(0x3c8bfbb395c60474ULL) {}
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    MinimizerParams(uint32_t k_val, uint32_t ell_val, uint32_t spaces_val, uint64_t mask) :
        k(k_val), ell(ell_val), spaces(spaces_val), xor_mask(mask) {}
};

// LCA candidate structure (host and device compatible)
struct LCACandidate {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
    
    // Constructors
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    LCACandidate() : 
        minimizer_hash(0), lca_taxon(0), genome_count(0), uniqueness_score(0.0f) {}
    
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    LCACandidate(uint64_t hash, uint32_t taxon, uint32_t count, float score) :
        minimizer_hash(hash), lca_taxon(taxon), genome_count(count), uniqueness_score(score) {}
};

// GPU launch configuration
struct LaunchConfig {
    int blocks_x = 1;
    int blocks_y = 1;
    int blocks_z = 1;
    int threads_x = 256;
    int threads_y = 1;
    int threads_z = 1;
    size_t shared_memory_bytes = 0;
    void* stream = nullptr;
    
    LaunchConfig() = default;
    LaunchConfig(int blocks, int threads) : blocks_x(blocks), threads_x(threads) {}
    LaunchConfig(int bx, int by, int tx, int ty) : blocks_x(bx), blocks_y(by), threads_x(tx), threads_y(ty) {}
};

// ===========================
// Memory Management Types
// ===========================

// Memory configuration structure
struct MemoryConfig {
    size_t max_memory_fraction = 80;        // Percentage of GPU memory to use
    size_t reserved_memory_mb = 500;        // Reserved memory in MB
    size_t minimizer_capacity = 5000000;    // Default minimizer capacity
    size_t sequence_batch_size = 25;        // Default sequence batch size
    bool enable_memory_pooling = true;      // Enable memory pooling
    bool auto_scale_enabled = true;         // Enable auto-scaling
};

// Memory statistics structure
struct MemoryStats {
    size_t total_gpu_memory = 0;
    size_t available_memory = 0;
    size_t allocated_memory = 0;
    size_t current_sequence_memory = 0;
    size_t current_minimizer_memory = 0;
    size_t current_metadata_memory = 0;
    size_t peak_usage = 0;
    double memory_efficiency = 1.0;
};

// GPU batch processing data
struct GPUBatchData {
    char* d_sequence_data = nullptr;
    struct GPUGenomeInfo* d_genome_info = nullptr;
    struct GPUMinimizerHit* d_minimizer_hits = nullptr;
    uint32_t* d_hit_counts = nullptr;
    uint32_t* d_global_counter = nullptr;
    LCACandidate* d_lca_candidates = nullptr;
    
    size_t sequence_buffer_size = 0;
    size_t max_genomes = 0;
    size_t max_minimizers = 0;
    
    bool is_allocated = false;
    
    GPUBatchData() = default;
    ~GPUBatchData() = default;
    
    // Non-copyable
    GPUBatchData(const GPUBatchData&) = delete;
    GPUBatchData& operator=(const GPUBatchData&) = delete;
    
    // Movable
    GPUBatchData(GPUBatchData&& other) noexcept;
    GPUBatchData& operator=(GPUBatchData&& other) noexcept;
};

// ===========================
// Species and Phylogenetic Types
// ===========================

// Species tracking during genome processing
struct SpeciesTrackingData {
    std::unordered_map<std::string, uint32_t> sequence_id_to_species;
    std::unordered_map<uint32_t, uint16_t> species_genome_counts;
    std::unordered_map<uint32_t, std::string> species_names;
    
    // Methods
    void add_genome(const std::string& sequence_id, uint32_t species_taxid, const std::string& species_name);
    uint32_t get_species_for_sequence(const std::string& sequence_id) const;
    uint16_t get_genome_count_for_species(uint32_t species_taxid) const;
    size_t total_species() const;
    size_t total_genomes() const;
    void clear();
    
    // Default constructor
    SpeciesTrackingData() = default;
    
    // Copy and move constructors
    SpeciesTrackingData(const SpeciesTrackingData&) = default;
    SpeciesTrackingData& operator=(const SpeciesTrackingData&) = default;
    SpeciesTrackingData(SpeciesTrackingData&&) = default;
    SpeciesTrackingData& operator=(SpeciesTrackingData&&) = default;
};

// Phylogenetic LCA candidate with enhanced metadata
struct PhylogeneticLCACandidate {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
    
    // Phylogenetic extensions
    std::vector<uint32_t> contributing_species;
    std::vector<uint16_t> genome_counts_per_species;
    uint8_t phylogenetic_spread;
    uint8_t max_phylogenetic_distance;
    
    PhylogeneticLCACandidate();
    PhylogeneticLCACandidate(const LCACandidate& basic);
    
    // Copy and move semantics
    PhylogeneticLCACandidate(const PhylogeneticLCACandidate&) = default;
    PhylogeneticLCACandidate& operator=(const PhylogeneticLCACandidate&) = default;
    PhylogeneticLCACandidate(PhylogeneticLCACandidate&&) = default;
    PhylogeneticLCACandidate& operator=(PhylogeneticLCACandidate&&) = default;
};

// Contributing taxa storage for variable-length data
struct ContributingTaxaArrays {
    std::vector<uint32_t> taxa_ids;
    std::vector<uint8_t> phylogenetic_distances;
    std::vector<uint16_t> genome_counts_per_taxon;
    
    uint32_t add_entry(uint32_t taxon_id, uint8_t distance, uint16_t genome_count);
    size_t total_entries() const;
    void clear();
    
    ContributingTaxaArrays() = default;
    ContributingTaxaArrays(const ContributingTaxaArrays&) = default;
    ContributingTaxaArrays& operator=(const ContributingTaxaArrays&) = default;
    ContributingTaxaArrays(ContributingTaxaArrays&&) = default;
    ContributingTaxaArrays& operator=(ContributingTaxaArrays&&) = default;
};

// ===========================
// Classification Parameters
// ===========================

struct ClassificationParams {
    uint32_t k = 31;                    // K-mer size
    uint32_t ell = 31;                  // Minimizer length  
    uint32_t spaces = 7;                // Spaced k-mer pattern
    double confidence_threshold = 0.0;   // Classification confidence threshold
    int num_threads = 1;                // Number of CPU threads
    bool quick_mode = false;            // Quick classification mode
    bool paired_end_processing = false; // Handle paired-end reads
    bool use_paired_end_bonus = true;   // Apply paired-end concordance bonus
    float paired_concordance_weight = 2.0f; // Weight for paired concordance
    float min_pair_concordance = 0.5f;  // Minimum concordance for bonus
    
    ClassificationParams() = default;
    ClassificationParams(uint32_t k_val, uint32_t ell_val, uint32_t spaces_val) 
        : k(k_val), ell(ell_val), spaces(spaces_val) {}
};

// ===========================
// Build Statistics
// ===========================

struct GPUBuildStats {
    uint64_t total_sequences = 0;
    uint64_t total_bases = 0;
    uint64_t total_kmers_processed = 0;
    uint64_t valid_minimizers_extracted = 0;
    uint64_t unique_minimizers = 0;
    uint64_t lca_assignments = 0;
    
    // Timing information
    double sequence_processing_time = 0.0;
    double minimizer_extraction_time = 0.0;
    double lca_computation_time = 0.0;
    double database_construction_time = 0.0;
    
    // Performance metrics
    double processing_rate_mb_per_sec = 0.0;
    double minimizer_extraction_rate = 0.0;
    
    GPUBuildStats() = default;
    
    void reset() {
        *this = GPUBuildStats();
    }
    
    void print_summary() const;
};

// Enhanced build statistics with phylogenetic data
struct EnhancedBuildStats : public GPUBuildStats {
    // Additional phylogenetic statistics
    uint64_t species_represented = 0;
    uint64_t minimizers_with_phylo_data = 0;
    uint64_t phylogenetic_lca_computations = 0;
    double phylogenetic_processing_time = 0.0;
    size_t contributing_taxa_array_size = 0;
    
    // Quality metrics
    double phylogenetic_coverage_percentage = 0.0;
    double average_taxa_per_minimizer = 0.0;
    
    void print_enhanced_stats() const;
    void calculate_derived_stats();
};

// ===========================
// Database Format Definitions
// ===========================

// Streamlined minimizer metadata (28 bytes per minimizer)
struct StreamlinedMinimizerMetadata {
    uint64_t minimizer_hash;                 // 8 bytes
    uint32_t lca_taxon;                      // 4 bytes - backward compatibility
    uint32_t total_genome_count;             // 4 bytes
    uint32_t contributing_taxa_offset;       // 4 bytes - offset into external array
    uint16_t num_contributing_taxa;          // 2 bytes
    uint16_t ml_weight;                      // 2 bytes - ML confidence score (0-65535)
    uint16_t feature_flags;                  // 2 bytes - encoded features (GC, complexity, etc.)
    uint8_t phylogenetic_spread;             // 1 byte
    uint8_t max_phylogenetic_distance;       // 1 byte
    uint16_t reserved;                       // 2 bytes - for future use, total = 32 bytes
    
    StreamlinedMinimizerMetadata() : 
        minimizer_hash(0), lca_taxon(0), total_genome_count(0),
        contributing_taxa_offset(0), num_contributing_taxa(0),
        ml_weight(0), feature_flags(0),
        phylogenetic_spread(0), max_phylogenetic_distance(0), reserved(0) {}
        
    // Helper methods for ML fields
    float get_ml_confidence() const {
        return ml_weight / 65535.0f;
    }
    
    void set_ml_confidence(float confidence) {
        ml_weight = static_cast<uint16_t>(confidence * 65535.0f);
    }
    
    // Use the MinimizerFlags namespace helpers for feature_flags
};

// ===========================
// Genome Processing Types
// ===========================

struct GenomeFileInfo {
    std::string file_path;
    uint32_t taxon_id;
    size_t file_size;
    size_t sequence_count;
    std::string species_name;
    bool is_valid;
    
    GenomeFileInfo() : taxon_id(0), file_size(0), sequence_count(0), is_valid(false) {}
    GenomeFileInfo(const std::string& path, uint32_t taxon, const std::string& species = "") 
        : file_path(path), taxon_id(taxon), species_name(species), file_size(0), sequence_count(0), is_valid(true) {}
};

struct SequenceBatch {
    std::vector<std::string> sequences;
    std::vector<uint32_t> taxon_ids;
    std::vector<std::string> sequence_ids;
    size_t total_length;
    
    SequenceBatch() : total_length(0) {}
    
    void clear() {
        sequences.clear();
        taxon_ids.clear();
        sequence_ids.clear();
        total_length = 0;
    }
    
    size_t size() const { return sequences.size(); }
    bool empty() const { return sequences.empty(); }
};

// ===========================
// Error Handling Types
// ===========================

enum class DatabaseBuildError {
    SUCCESS = 0,
    INVALID_INPUT,
    MEMORY_ALLOCATION_FAILED,
    CUDA_ERROR,
    FILE_IO_ERROR,
    TAXONOMY_LOAD_ERROR,
    INVALID_CONFIGURATION,
    PROCESSING_ERROR,
    SERIALIZATION_ERROR,
    VALIDATION_ERROR
};

struct BuildErrorInfo {
    DatabaseBuildError error_code;
    std::string error_message;
    std::string error_context;
    int line_number;
    std::string file_name;
    
    BuildErrorInfo() : error_code(DatabaseBuildError::SUCCESS), line_number(0) {}
    BuildErrorInfo(DatabaseBuildError code, const std::string& message, const std::string& context = "") 
        : error_code(code), error_message(message), error_context(context), line_number(0) {}
    
    bool is_success() const { return error_code == DatabaseBuildError::SUCCESS; }
    void print_error() const;
};

// ===========================
// Database Format Types
// ===========================

enum class DatabaseFormat {
    STANDARD_KRAKEN2,
    ENHANCED_PHYLO,
    COMPACT_BINARY,
    CUSTOM_FORMAT
};

// ===========================
// File Processing Types
// ===========================

struct FileProcessingConfig {
    bool validate_sequences = true;
    bool skip_invalid_files = true;
    size_t max_file_size_mb = 10000;
    bool enable_parallel_loading = true;
    int num_worker_threads = 4;
    std::string supported_extensions = ".fna,.fa,.fasta";
    
    // Additional fields used in implementation
    size_t max_file_count = 10000;
    bool progress_reporting = false;
    size_t progress_interval = 100;
    size_t max_file_size = 10000000000ULL;  // 10GB
    size_t max_sequence_length = 100000000;  // 100MB
};

struct FileProcessingStats {
    double processing_time = 0.0;
    size_t files_processed = 0;
    size_t total_bytes = 0;
    size_t sequences_loaded = 0;
    size_t invalid_files_skipped = 0;
    
    // Additional fields used in implementation
    size_t files_found = 0;
    size_t files_skipped = 0;
    size_t total_sequences = 0;
    size_t total_bases = 0;
    size_t processing_errors = 0;
    double average_file_size = 0.0;
    
    void reset() {
        processing_time = 0.0;
        files_processed = 0;
        total_bytes = 0;
        sequences_loaded = 0;
        invalid_files_skipped = 0;
        files_found = 0;
        files_skipped = 0;
        total_sequences = 0;
        total_bases = 0;
        processing_errors = 0;
        average_file_size = 0.0;
    }
    
    void print_summary() const {
        std::cout << "File Processing Stats:" << std::endl;
        std::cout << "  Files processed: " << files_processed << std::endl;
        std::cout << "  Sequences loaded: " << sequences_loaded << std::endl;
        std::cout << "  Total bytes: " << (total_bytes / 1024 / 1024) << " MB" << std::endl;
        std::cout << "  Processing time: " << processing_time << "s" << std::endl;
        std::cout << "  Invalid files skipped: " << invalid_files_skipped << std::endl;
    }
};

// ===========================
// Progress Monitoring Types
// ===========================

enum class BuildStage {
    INITIALIZATION,
    FILE_DISCOVERY,
    TAXONOMY_LOADING,
    MEMORY_ALLOCATION,
    SEQUENCE_PROCESSING,
    MINIMIZER_EXTRACTION,
    LCA_COMPUTATION,
    PHYLOGENETIC_ANALYSIS,
    DATABASE_SERIALIZATION,
    VALIDATION,
    CLEANUP,
    COMPLETE
};

struct ProgressInfo {
    BuildStage current_stage;
    double stage_progress_percent;
    double overall_progress_percent;
    std::string current_task;
    size_t items_processed;
    size_t total_items;
    double estimated_time_remaining_seconds;
    
    ProgressInfo() : current_stage(BuildStage::INITIALIZATION), stage_progress_percent(0.0),
                    overall_progress_percent(0.0), items_processed(0), total_items(0),
                    estimated_time_remaining_seconds(0.0) {}
    
    void update_progress(size_t processed, size_t total);
    std::string get_stage_name() const;
    void print_progress() const;
};

// ===========================
// GPU Resource Types
// ===========================

struct GPUResourceInfo {
    int device_id;
    std::string device_name;
    size_t total_memory;
    size_t free_memory;
    int major_compute_capability;
    int minor_compute_capability;
    int multiprocessor_count;
    int max_threads_per_multiprocessor;
    bool supports_required_features;
    
    GPUResourceInfo() : device_id(0), total_memory(0), free_memory(0),
                       major_compute_capability(0), minor_compute_capability(0),
                       multiprocessor_count(0), max_threads_per_multiprocessor(0),
                       supports_required_features(false) {}
    
    void print_info() const;
    bool is_sufficient_for_build() const;
};

// ===========================
// Validation Types
// ===========================

struct ValidationResult {
    bool is_valid;
    std::vector<std::string> warnings;
    std::vector<std::string> errors;
    std::string validation_context;
    
    ValidationResult() : is_valid(true) {}
    
    void add_warning(const std::string& warning) {
        warnings.push_back(warning);
    }
    
    void add_error(const std::string& error) {
        errors.push_back(error);
        is_valid = false;
    }
    
    void print_results() const;
    bool has_warnings() const { return !warnings.empty(); }
    bool has_errors() const { return !errors.empty(); }
};

// ===========================
// Constants and Limits
// ===========================

namespace Constants {
    // Biological constants
    constexpr uint32_t MIN_K_VALUE = 15;
    constexpr uint32_t MAX_K_VALUE = 31;
    constexpr uint32_t DEFAULT_K_VALUE = 31;
    constexpr uint32_t DEFAULT_ELL_VALUE = 31;
    constexpr uint32_t DEFAULT_SPACES_VALUE = 7;
    
    // File processing limits
    constexpr size_t MAX_SEQUENCE_LENGTH = 100000000;      // 100MB per sequence
    constexpr size_t MAX_FILE_SIZE = 10000000000ULL;       // 10GB per file
    constexpr size_t MAX_FILES_PER_DIRECTORY = 100000;     // 100K files max
    constexpr size_t MIN_SEQUENCE_LENGTH = 50;             // 50bp minimum
    
    // Memory limits
    constexpr size_t MAX_GPU_MEMORY_USAGE_PERCENT = 95;
    constexpr size_t DEFAULT_GPU_MEMORY_USAGE_PERCENT = 80;
    constexpr size_t MIN_GPU_MEMORY_MB = 1024;             // 1GB minimum
    
    // Performance constants
    constexpr int DEFAULT_THREADS_PER_BLOCK = 256;
    constexpr int MAX_BATCH_SIZE = 100;
    constexpr int DEFAULT_BATCH_SIZE = 25;
    
    // Database constants
    constexpr uint32_t ROOT_TAXON_ID = 1;
    constexpr uint32_t INVALID_TAXON_ID = 0;
    constexpr uint64_t INVALID_HASH = UINT64_MAX;
    
    // String constants
    const std::string DEFAULT_DATABASE_NAME = "kraken_database";
    const std::string HASH_TABLE_FILENAME = "hash_table.k2d";
    const std::string TAXONOMY_FILENAME = "taxonomy.tsv";
    const std::string CONFIG_FILENAME = "config.txt";
}

// ===========================
// Utility Functions
// ===========================

namespace TypeUtils {
    // String conversion utilities
    std::string build_stage_to_string(BuildStage stage);
    std::string error_code_to_string(DatabaseBuildError error);
    std::string format_file_size(size_t bytes);
    std::string format_duration(double seconds);
    std::string format_number_with_commas(uint64_t number);
    
    // Validation utilities
    bool is_valid_k_value(uint32_t k);
    bool is_valid_taxon_id(uint32_t taxon_id);
    bool is_valid_file_path(const std::string& path);
    
    // Hash utilities
    uint64_t combine_hashes(uint64_t hash1, uint64_t hash2);
    uint32_t simple_string_hash(const std::string& str);
    
    // Memory utilities
    size_t calculate_memory_requirement(size_t count, size_t item_size, double overhead_factor = 1.1);
    bool is_power_of_two(size_t value);
    size_t next_power_of_two(size_t value);
    
    // Implementation of commonly used functions
    inline std::string build_stage_to_string(BuildStage stage) {
        switch (stage) {
            case BuildStage::INITIALIZATION: return "Initialization";
            case BuildStage::FILE_DISCOVERY: return "File Discovery";
            case BuildStage::TAXONOMY_LOADING: return "Taxonomy Loading";
            case BuildStage::MEMORY_ALLOCATION: return "Memory Allocation";
            case BuildStage::SEQUENCE_PROCESSING: return "Sequence Processing";
            case BuildStage::MINIMIZER_EXTRACTION: return "Minimizer Extraction";
            case BuildStage::LCA_COMPUTATION: return "LCA Computation";
            case BuildStage::PHYLOGENETIC_ANALYSIS: return "Phylogenetic Analysis";
            case BuildStage::DATABASE_SERIALIZATION: return "Database Serialization";
            case BuildStage::VALIDATION: return "Validation";
            case BuildStage::CLEANUP: return "Cleanup";
            case BuildStage::COMPLETE: return "Complete";
            default: return "Unknown";
        }
    }
    
    inline std::string error_code_to_string(DatabaseBuildError error) {
        switch (error) {
            case DatabaseBuildError::SUCCESS: return "Success";
            case DatabaseBuildError::INVALID_INPUT: return "Invalid Input";
            case DatabaseBuildError::MEMORY_ALLOCATION_FAILED: return "Memory Allocation Failed";
            case DatabaseBuildError::CUDA_ERROR: return "CUDA Error";
            case DatabaseBuildError::FILE_IO_ERROR: return "File I/O Error";
            case DatabaseBuildError::TAXONOMY_LOAD_ERROR: return "Taxonomy Load Error";
            case DatabaseBuildError::INVALID_CONFIGURATION: return "Invalid Configuration";
            case DatabaseBuildError::PROCESSING_ERROR: return "Processing Error";
            case DatabaseBuildError::SERIALIZATION_ERROR: return "Serialization Error";
            case DatabaseBuildError::VALIDATION_ERROR: return "Validation Error";
            default: return "Unknown Error";
        }
    }
    
    inline std::string format_file_size(size_t bytes) {
        const char* units[] = {"B", "KB", "MB", "GB", "TB"};
        int unit = 0;
        double size = static_cast<double>(bytes);
        
        while (size >= 1024 && unit < 4) {
            size /= 1024;
            unit++;
        }
        
        char buffer[64];
        snprintf(buffer, sizeof(buffer), "%.1f %s", size, units[unit]);
        return std::string(buffer);
    }
    
    inline std::string format_duration(double seconds) {
        if (seconds < 60) {
            char buffer[32];
            snprintf(buffer, sizeof(buffer), "%.2fs", seconds);
            return std::string(buffer);
        } else if (seconds < 3600) {
            int minutes = static_cast<int>(seconds / 60);
            double remaining_seconds = seconds - minutes * 60;
            char buffer[32];
            snprintf(buffer, sizeof(buffer), "%dm %.1fs", minutes, remaining_seconds);
            return std::string(buffer);
        } else {
            int hours = static_cast<int>(seconds / 3600);
            int minutes = static_cast<int>((seconds - hours * 3600) / 60);
            char buffer[32];
            snprintf(buffer, sizeof(buffer), "%dh %dm", hours, minutes);
            return std::string(buffer);
        }
    }
    
    inline bool is_valid_k_value(uint32_t k) {
        return k >= Constants::MIN_K_VALUE && k <= Constants::MAX_K_VALUE;
    }
    
    inline bool is_valid_taxon_id(uint32_t taxon_id) {
        return taxon_id > 0 && taxon_id != Constants::INVALID_TAXON_ID;
    }
    
    inline uint64_t combine_hashes(uint64_t hash1, uint64_t hash2) {
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
    
    inline size_t calculate_memory_requirement(size_t count, size_t item_size, double overhead_factor) {
        return static_cast<size_t>(count * item_size * overhead_factor);
    }
    
    inline bool is_power_of_two(size_t value) {
        return value != 0 && (value & (value - 1)) == 0;
    }
    
    inline size_t next_power_of_two(size_t value) {
        if (value == 0) return 1;
        value--;
        value |= value >> 1;
        value |= value >> 2;
        value |= value >> 4;
        value |= value >> 8;
        value |= value >> 16;
        value |= value >> 32;
        return value + 1;
    }
}

// ===========================
// Compatibility Macros
// ===========================

// CUDA compatibility
#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#define DEVICE __device__
#define HOST __host__
#define CUDA_CALLABLE HOST_DEVICE
#else
#define HOST_DEVICE
#define DEVICE
#define HOST
#define CUDA_CALLABLE
#endif

// Compiler-specific optimizations
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

// Memory alignment
#ifdef _MSC_VER
#define ALIGN(x) __declspec(align(x))
#else
#define ALIGN(x) __attribute__((aligned(x)))
#endif

// Force inline for performance-critical functions
#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#elif defined(__GNUC__)
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

// ===========================
// Inline Method Implementations
// ===========================

// SpeciesTrackingData method implementations
inline void SpeciesTrackingData::add_genome(const std::string& sequence_id, uint32_t species_taxid, const std::string& species_name) {
    sequence_id_to_species[sequence_id] = species_taxid;
    species_genome_counts[species_taxid]++;
    if (species_names.find(species_taxid) == species_names.end()) {
        species_names[species_taxid] = species_name;
    }
}

inline uint32_t SpeciesTrackingData::get_species_for_sequence(const std::string& sequence_id) const {
    auto it = sequence_id_to_species.find(sequence_id);
    return (it != sequence_id_to_species.end()) ? it->second : 0;
}

inline uint16_t SpeciesTrackingData::get_genome_count_for_species(uint32_t species_taxid) const {
    auto it = species_genome_counts.find(species_taxid);
    return (it != species_genome_counts.end()) ? it->second : 0;
}

inline size_t SpeciesTrackingData::total_species() const {
    return species_genome_counts.size();
}

inline size_t SpeciesTrackingData::total_genomes() const {
    size_t total = 0;
    for (const auto& [species, count] : species_genome_counts) {
        total += count;
    }
    return total;
}

inline void SpeciesTrackingData::clear() {
    sequence_id_to_species.clear();
    species_genome_counts.clear();
    species_names.clear();
}

// PhylogeneticLCACandidate method implementations
inline PhylogeneticLCACandidate::PhylogeneticLCACandidate() 
    : minimizer_hash(0), lca_taxon(0), genome_count(0), uniqueness_score(0.0f),
      phylogenetic_spread(0), max_phylogenetic_distance(0) {
}

inline PhylogeneticLCACandidate::PhylogeneticLCACandidate(const LCACandidate& basic) 
    : minimizer_hash(basic.minimizer_hash), lca_taxon(basic.lca_taxon),
      genome_count(basic.genome_count), uniqueness_score(basic.uniqueness_score),
      phylogenetic_spread(0), max_phylogenetic_distance(0) {
}

// ContributingTaxaArrays method implementations
inline uint32_t ContributingTaxaArrays::add_entry(uint32_t taxon_id, uint8_t distance, uint16_t genome_count) {
    uint32_t offset = taxa_ids.size();
    taxa_ids.push_back(taxon_id);
    phylogenetic_distances.push_back(distance);
    genome_counts_per_taxon.push_back(genome_count);
    return offset;
}

inline size_t ContributingTaxaArrays::total_entries() const {
    return taxa_ids.size();
}

inline void ContributingTaxaArrays::clear() {
    taxa_ids.clear();
    phylogenetic_distances.clear();
    genome_counts_per_taxon.clear();
}

// GPUBuildStats method implementations
inline void GPUBuildStats::print_summary() const {
    std::cout << "\n=== BUILD STATISTICS ===" << std::endl;
    std::cout << "Total sequences: " << total_sequences << std::endl;
    std::cout << "Total bases: " << total_bases << std::endl;
    std::cout << "Total k-mers processed: " << total_kmers_processed << std::endl;
    std::cout << "Valid minimizers extracted: " << valid_minimizers_extracted << std::endl;
    std::cout << "Unique minimizers: " << unique_minimizers << std::endl;
    std::cout << "LCA assignments: " << lca_assignments << std::endl;
    
    std::cout << "\nTiming Information:" << std::endl;
    std::cout << "  Sequence processing: " << sequence_processing_time << "s" << std::endl;
    std::cout << "  Minimizer extraction: " << minimizer_extraction_time << "s" << std::endl;
    std::cout << "  LCA computation: " << lca_computation_time << "s" << std::endl;
    std::cout << "  Database construction: " << database_construction_time << "s" << std::endl;
    
    if (total_kmers_processed > 0) {
        double compression = (double)valid_minimizers_extracted / total_kmers_processed;
        std::cout << "Compression ratio: " << compression << " (" 
                  << (1.0/compression) << "x reduction)" << std::endl;
    }
}

// EnhancedBuildStats method implementations
inline void EnhancedBuildStats::print_enhanced_stats() const {
    std::cout << "\n=== ENHANCED BUILD STATISTICS ===" << std::endl;
    std::cout << "Species represented: " << species_represented << std::endl;
    std::cout << "Minimizers with phylogenetic data: " << minimizers_with_phylo_data << std::endl;
    std::cout << "Phylogenetic LCA computations: " << phylogenetic_lca_computations << std::endl;
    std::cout << "Contributing taxa array entries: " << contributing_taxa_array_size << std::endl;
    std::cout << "Phylogenetic processing time: " << phylogenetic_processing_time << "s" << std::endl;
    
    if (unique_minimizers > 0) {
        double phylo_coverage = (double)minimizers_with_phylo_data / unique_minimizers * 100.0;
        std::cout << "Phylogenetic coverage: " << phylo_coverage << "%" << std::endl;
    }
    
    if (contributing_taxa_array_size > 0 && minimizers_with_phylo_data > 0) {
        double avg_taxa_per_minimizer = (double)contributing_taxa_array_size / minimizers_with_phylo_data;
        std::cout << "Average taxa per minimizer: " << avg_taxa_per_minimizer << std::endl;
    }
}

inline void EnhancedBuildStats::calculate_derived_stats() {
    if (unique_minimizers > 0) {
        phylogenetic_coverage_percentage = (double)minimizers_with_phylo_data / unique_minimizers * 100.0;
    }
    
    if (minimizers_with_phylo_data > 0) {
        average_taxa_per_minimizer = (double)contributing_taxa_array_size / minimizers_with_phylo_data;
    }
}

// GPUBatchData move constructor and assignment
inline GPUBatchData::GPUBatchData(GPUBatchData&& other) noexcept 
    : d_sequence_data(other.d_sequence_data),
      d_genome_info(other.d_genome_info),
      d_minimizer_hits(other.d_minimizer_hits),
      d_hit_counts(other.d_hit_counts),
      d_global_counter(other.d_global_counter),
      d_lca_candidates(other.d_lca_candidates),
      sequence_buffer_size(other.sequence_buffer_size),
      max_genomes(other.max_genomes),
      max_minimizers(other.max_minimizers),
      is_allocated(other.is_allocated) {
    
    // Reset other object
    other.d_sequence_data = nullptr;
    other.d_genome_info = nullptr;
    other.d_minimizer_hits = nullptr;
    other.d_hit_counts = nullptr;
    other.d_global_counter = nullptr;
    other.d_lca_candidates = nullptr;
    other.sequence_buffer_size = 0;
    other.max_genomes = 0;
    other.max_minimizers = 0;
    other.is_allocated = false;
}

inline GPUBatchData& GPUBatchData::operator=(GPUBatchData&& other) noexcept {
    if (this != &other) {
        // Move data
        d_sequence_data = other.d_sequence_data;
        d_genome_info = other.d_genome_info;
        d_minimizer_hits = other.d_minimizer_hits;
        d_hit_counts = other.d_hit_counts;
        d_global_counter = other.d_global_counter;
        d_lca_candidates = other.d_lca_candidates;
        sequence_buffer_size = other.sequence_buffer_size;
        max_genomes = other.max_genomes;
        max_minimizers = other.max_minimizers;
        is_allocated = other.is_allocated;
        
        // Reset other object
        other.d_sequence_data = nullptr;
        other.d_genome_info = nullptr;
        other.d_minimizer_hits = nullptr;
        other.d_hit_counts = nullptr;
        other.d_global_counter = nullptr;
        other.d_lca_candidates = nullptr;
        other.sequence_buffer_size = 0;
        other.max_genomes = 0;
        other.max_minimizers = 0;
        other.is_allocated = false;
    }
    return *this;
}

// Utility function implementations
inline void ProgressInfo::update_progress(size_t processed, size_t total) {
    items_processed = processed;
    total_items = total;
    if (total > 0) {
        stage_progress_percent = (double)processed / total * 100.0;
    }
}

inline std::string ProgressInfo::get_stage_name() const {
    return TypeUtils::build_stage_to_string(current_stage);
}

inline void ProgressInfo::print_progress() const {
    std::cout << "[" << get_stage_name() << "] " 
              << stage_progress_percent << "% (" 
              << items_processed << "/" << total_items << ")" << std::endl;
}

inline void GPUResourceInfo::print_info() const {
    std::cout << "GPU Device " << device_id << ": " << device_name << std::endl;
    std::cout << "  Total memory: " << (total_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Free memory: " << (free_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Compute capability: " << major_compute_capability << "." << minor_compute_capability << std::endl;
    std::cout << "  Multiprocessors: " << multiprocessor_count << std::endl;
    std::cout << "  Supports required features: " << (supports_required_features ? "Yes" : "No") << std::endl;
}

inline bool GPUResourceInfo::is_sufficient_for_build() const {
    return supports_required_features && 
           (total_memory / 1024 / 1024) >= Constants::MIN_GPU_MEMORY_MB &&
           major_compute_capability >= 3;
}

inline void ValidationResult::print_results() const {
    if (is_valid) {
        std::cout << "✓ Validation passed";
    } else {
        std::cout << "✗ Validation failed";
    }
    
    if (!validation_context.empty()) {
        std::cout << " (" << validation_context << ")";
    }
    std::cout << std::endl;
    
    for (const auto& warning : warnings) {
        std::cout << "  Warning: " << warning << std::endl;
    }
    
    for (const auto& error : errors) {
        std::cout << "  Error: " << error << std::endl;
    }
}

inline void BuildErrorInfo::print_error() const {
    std::cout << "Error: " << error_message;
    if (!error_context.empty()) {
        std::cout << " (Context: " << error_context << ")";
    }
    if (line_number > 0) {
        std::cout << " at " << file_name << ":" << line_number;
    }
    std::cout << std::endl;
}

#endif // GPU_KRAKEN_TYPES_H