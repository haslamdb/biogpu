// core/gpu_database_builder_core.h
// Main orchestration class for GPU Kraken database builder
// Coordinates all specialized modules and provides public API

#ifndef GPU_DATABASE_BUILDER_CORE_H
#define GPU_DATABASE_BUILDER_CORE_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>

// Include necessary headers
#include "../processing/genome_file_processor.h"
#include "../taxonomy/taxonomy_processor.h"
#include "../output/database_serializer.h"
#include "../gpu_kraken_types.h"

// Forward declarations
class GPUMemoryManager;
struct ClassificationParams;
struct GPUBuildStats;
struct EnhancedBuildStats;

// Configuration for the entire database building process
struct DatabaseBuildConfig {
    // Core parameters
    uint32_t k_value = 31;
    uint32_t ell_value = 31;
    uint32_t spaces_value = 7;
    double subsampling_rate = 1.0;
    uint64_t toggle_mask = 0xe37e28c4271b5a2dULL;
    
    // Memory configuration
    MemoryConfig memory_config;
    
    // File processing configuration
    FileProcessingConfig file_config;
    
    // Output format selection
    DatabaseFormat output_format = DatabaseFormat::STANDARD_KRAKEN2;
    bool save_enhanced_format = false;
    bool save_compatibility_layer = true;
    
    // Build behavior
    bool enable_phylogenetic_analysis = false;
    bool validate_output = true;
    bool generate_checksums = true;
    
    // Performance settings
    bool enable_progress_reporting = true;
    int progress_report_interval = 1000;
    bool enable_detailed_timing = true;
    
    // Debugging and validation
    bool enable_memory_validation = false;
    bool enable_intermediate_saves = false;
    std::string debug_output_dir;
};

// Main GPU Kraken Database Builder coordinating class
class GPUKrakenDatabaseBuilder {
private:
    // Configuration
    DatabaseBuildConfig config_;
    std::string output_directory_;
    
    // Specialized modules
    std::unique_ptr<GPUMemoryManager> memory_manager_;
    std::unique_ptr<GenomeFileProcessor> file_processor_;
    std::unique_ptr<EnhancedNCBITaxonomyProcessor> taxonomy_processor_;
    std::unique_ptr<StandardDatabaseSerializer> standard_serializer_;
    std::unique_ptr<EnhancedDatabaseSerializer> enhanced_serializer_;
    
    // Build state
    bool initialized_;
    bool cuda_context_ready_;
    bool modules_initialized_;
    
    // Data storage
    std::vector<std::string> genome_files_;
    std::vector<uint32_t> genome_taxon_ids_;
    std::unordered_map<uint32_t, std::string> taxon_names_;
    std::unordered_map<uint32_t, uint32_t> taxon_parents_;
    std::vector<LCACandidate> all_lca_candidates_;
    std::vector<PhylogeneticLCACandidate> phylogenetic_candidates_;
    ContributingTaxaArrays contributing_taxa_arrays_;
    SpeciesTrackingData species_tracking_;
    
    // Statistics
    GPUBuildStats build_stats_;
    EnhancedBuildStats enhanced_stats_;
    
    // GPU kernels access
    MinimizerParams minimizer_params_;
    
public:
    // Constructor and destructor
    explicit GPUKrakenDatabaseBuilder(const std::string& output_dir,
                                    const DatabaseBuildConfig& config = DatabaseBuildConfig());
    ~GPUKrakenDatabaseBuilder();
    
    // Main pipeline methods (public API)
    bool build_database_from_genomes(
        const std::string& genome_library_path,
        const std::string& taxonomy_path = ""
    );
    
    bool build_database_from_file_list(
        const std::string& file_list_path,
        const std::string& taxonomy_path = ""
    );
    
    bool build_database_from_concatenated_fna(
        const std::string& fna_file_path,
        const std::string& taxonomy_nodes_path = "",
        const std::string& taxonomy_names_path = "",
        const std::string& compact_taxonomy_path = ""
    );
    
    bool build_database_from_streaming_fna(
        const std::string& fna_file_path,
        const std::string& taxonomy_path = ""
    );
    
    // Configuration methods
    bool configure_memory_settings(const MemoryConfig& memory_config);
    bool configure_file_processing(const FileProcessingConfig& file_config);
    bool configure_output_format(DatabaseFormat format, bool enhanced = false);
    bool configure_phylogenetic_analysis(bool enable);
    
    // Advanced configuration
    void set_minimizer_parameters(uint32_t k, uint32_t ell, uint32_t spaces);
    void set_subsampling_rate(double rate);
    void enable_debug_mode(const std::string& debug_dir);
    void enable_detailed_validation(bool enable);
    
    // Status and monitoring
    bool is_initialized() const { return initialized_ && modules_initialized_; }
    bool is_cuda_ready() const { return cuda_context_ready_; }
    
    // Statistics access
    const GPUBuildStats& get_build_statistics() const { return build_stats_; }
    const EnhancedBuildStats& get_enhanced_statistics() const { return enhanced_stats_; }
    
    // Progress monitoring
    void print_build_progress() const;
    void print_memory_usage() const;
    void print_module_status() const;
    
    // Validation and testing
    bool validate_configuration() const;
    bool test_gpu_functionality();
    bool validate_database_output();
    
    // Utility methods
    static DatabaseBuildConfig create_default_config();
    static DatabaseBuildConfig create_high_memory_config();
    static DatabaseBuildConfig create_streaming_config();
    static DatabaseBuildConfig create_enhanced_phylo_config();
    
private:
    // Initialization methods
    bool initialize_all_modules();
    bool initialize_cuda_context();
    bool initialize_memory_manager();
    bool initialize_file_processor();
    bool initialize_taxonomy_processor();
    bool initialize_serializers();
    
    // Core processing pipeline
    bool execute_build_pipeline();
    bool load_and_validate_inputs(const std::string& input_path, const std::string& taxonomy_path);
    bool process_genomes_with_gpu();
    bool enhance_with_phylogenetic_data();
    bool save_database_outputs();
    
    // Step-by-step processing methods
    bool load_genome_files(const std::string& library_path);
    bool load_taxonomy_data(const std::string& taxonomy_path);
    bool process_sequence_batches();
    bool process_sequence_batch(const std::vector<std::string>& sequences, 
                               const std::vector<uint32_t>& taxon_ids);
    bool compute_lca_assignments();
    bool merge_and_deduplicate_candidates();
    
    // Module coordination
    bool coordinate_memory_allocation();
    bool validate_inter_module_consistency();
    bool handle_module_errors();
    
    // Error handling and recovery
    bool attempt_error_recovery(const std::string& error_context);
    void cleanup_on_failure();
    void save_intermediate_state();
    bool load_intermediate_state();
    
    // Progress and timing
    void update_build_progress(const std::string& stage, double progress_percent);
    void record_timing_checkpoint(const std::string& stage_name);
    
    // Validation helpers
    bool validate_input_parameters() const;
    bool validate_gpu_resources() const;
    bool validate_file_inputs() const;
    bool validate_taxonomy_consistency() const;
    
    // Statistics aggregation
    void aggregate_module_statistics();
    void calculate_derived_statistics();
    void update_performance_metrics();
};

// Factory class for creating configured database builders
class DatabaseBuilderFactory {
public:
    // Standard configurations
    static std::unique_ptr<GPUKrakenDatabaseBuilder> create_standard_builder(
        const std::string& output_dir);
    
    static std::unique_ptr<GPUKrakenDatabaseBuilder> create_high_capacity_builder(
        const std::string& output_dir);
    
    static std::unique_ptr<GPUKrakenDatabaseBuilder> create_streaming_builder(
        const std::string& output_dir);
    
    static std::unique_ptr<GPUKrakenDatabaseBuilder> create_phylogenetic_builder(
        const std::string& output_dir);
    
    // Custom configuration
    static std::unique_ptr<GPUKrakenDatabaseBuilder> create_custom_builder(
        const std::string& output_dir,
        const DatabaseBuildConfig& config);
    
    // Automatic configuration based on system resources
    static std::unique_ptr<GPUKrakenDatabaseBuilder> create_auto_configured_builder(
        const std::string& output_dir,
        const std::string& input_estimate_path = "");
    
private:
    static DatabaseBuildConfig detect_optimal_configuration();
    static MemoryConfig detect_optimal_memory_config();
    static FileProcessingConfig detect_optimal_file_config();
};

// Build session manager for handling multiple builds and caching
class DatabaseBuildSession {
private:
    std::vector<std::unique_ptr<GPUKrakenDatabaseBuilder>> active_builders_;
    std::string session_cache_dir_;
    bool enable_build_caching_;
    
public:
    explicit DatabaseBuildSession(const std::string& cache_dir = "", bool enable_caching = false);
    ~DatabaseBuildSession();
    
    // Session management
    bool start_session();
    bool end_session();
    void cleanup_session();
    
    // Builder management
    size_t add_builder(std::unique_ptr<GPUKrakenDatabaseBuilder> builder);
    bool remove_builder(size_t builder_id);
    GPUKrakenDatabaseBuilder* get_builder(size_t builder_id);
    
    // Batch operations
    bool execute_all_builds();
    bool execute_builds_parallel();
    
    // Caching and reuse
    bool cache_intermediate_results(size_t builder_id);
    bool reuse_cached_results(size_t builder_id, const std::string& cache_key);
    
    // Session statistics
    void print_session_summary() const;
    size_t get_active_builder_count() const { return active_builders_.size(); }
};

// Utility functions for configuration and setup
namespace DatabaseBuilderUtils {
    // Configuration validation
    bool validate_build_config(const DatabaseBuildConfig& config);
    bool validate_system_requirements(const DatabaseBuildConfig& config);
    
    // Resource estimation
    size_t estimate_memory_requirements(const std::string& input_path, const DatabaseBuildConfig& config);
    double estimate_build_time(const std::string& input_path, const DatabaseBuildConfig& config);
    
    // Configuration optimization
    DatabaseBuildConfig optimize_config_for_input(const std::string& input_path);
    MemoryConfig optimize_memory_for_gpu();
    
    // Diagnostics
    bool run_system_diagnostics();
    bool test_gpu_performance();
    void print_system_capabilities();
    
    // Migration utilities
    bool migrate_old_database_format(const std::string& old_db_path, const std::string& new_db_path);
    bool upgrade_database_version(const std::string& db_path);
}

#endif // GPU_DATABASE_BUILDER_CORE_H
