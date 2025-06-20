// output/database_serializer.h
// Database serialization and output formatting
// Handles multiple database formats and backward compatibility

#ifndef DATABASE_SERIALIZER_H
#define DATABASE_SERIALIZER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <memory>

// Forward declarations
struct LCACandidate;
struct PhylogeneticLCACandidate;
struct ContributingTaxaArrays;
struct StreamlinedMinimizerMetadata;
struct GPUBuildStats;
struct EnhancedBuildStats;

// Database format specifications
enum class DatabaseFormat {
    STANDARD_KRAKEN2,    // Standard Kraken2 format
    ENHANCED_PHYLO,      // Enhanced format with phylogenetic data
    COMPACT_BINARY,      // Compact binary format
    STREAMING_FRIENDLY   // Optimized for streaming access
};

// Database metadata structure
struct DatabaseMetadata {
    uint64_t version = 2;
    uint32_t k_value = 31;
    uint32_t ell_value = 31;
    uint32_t spaces_value = 7;
    double subsampling_rate = 1.0;
    uint64_t min_clear_hash_value = 0;
    uint64_t toggle_mask = 0xe37e28c4271b5a2dULL;
    
    // Database statistics
    uint64_t total_sequences = 0;
    uint64_t total_bases = 0;
    uint64_t unique_minimizers = 0;
    uint64_t species_count = 0;
    uint64_t genome_count = 0;
    
    // Timing information
    double build_time = 0.0;
    double processing_rate = 0.0;
    
    // File paths and checksums
    std::string hash_table_file;
    std::string taxonomy_file;
    std::string config_file;
    std::string checksum_hash;
};

// Standard Kraken2 database serializer
class StandardDatabaseSerializer {
private:
    std::string output_directory_;
    DatabaseMetadata metadata_;
    
public:
    explicit StandardDatabaseSerializer(const std::string& output_dir);
    
    // Main serialization methods
    bool save_standard_database(
        const std::vector<LCACandidate>& candidates,
        const std::unordered_map<uint32_t, std::string>& taxon_names,
        const std::unordered_map<uint32_t, uint32_t>& taxon_parents,
        const GPUBuildStats& stats
    );
    
    // Individual component saving
    bool save_hash_table(const std::vector<LCACandidate>& candidates);
    bool save_taxonomy_file(const std::unordered_map<uint32_t, std::string>& taxon_names,
                           const std::unordered_map<uint32_t, uint32_t>& taxon_parents);
    bool save_config_file(const GPUBuildStats& stats);
    bool save_build_statistics(const GPUBuildStats& stats);
    
    // Validation and verification
    bool validate_database_consistency();
    bool generate_database_checksum();
    
    // Metadata access
    const DatabaseMetadata& get_metadata() const { return metadata_; }
    void set_metadata(const DatabaseMetadata& metadata) { metadata_ = metadata; }
    
private:
    // Helper methods
    std::string get_hash_table_path() const;
    std::string get_taxonomy_path() const;
    std::string get_config_path() const;
    std::string get_stats_path() const;
    
    bool create_output_directory();
    bool validate_output_paths();
};

// Enhanced database serializer with phylogenetic data
class EnhancedDatabaseSerializer {
private:
    std::string output_directory_;
    DatabaseMetadata metadata_;
    
public:
    explicit EnhancedDatabaseSerializer(const std::string& output_dir);
    
    // Enhanced serialization methods
    bool save_enhanced_database(
        const std::vector<PhylogeneticLCACandidate>& phylo_candidates,
        const ContributingTaxaArrays& contributing_taxa,
        const std::unordered_map<uint32_t, std::string>& taxon_names,
        const EnhancedBuildStats& stats
    );
    
    // Enhanced format components
    bool save_streamlined_hash_table(const std::vector<StreamlinedMinimizerMetadata>& metadata);
    bool save_contributing_taxa_arrays(const ContributingTaxaArrays& taxa_arrays);
    bool save_enhanced_config(const EnhancedBuildStats& stats);
    bool save_phylogenetic_summary(const std::vector<PhylogeneticLCACandidate>& candidates);
    bool save_species_mapping(const std::unordered_map<uint32_t, std::string>& species_names,
                             const std::unordered_map<uint32_t, uint16_t>& genome_counts);
    
    // Backward compatibility
    bool save_standard_compatibility_layer(const std::vector<PhylogeneticLCACandidate>& phylo_candidates);
    
    // Format validation
    bool validate_enhanced_format();
    bool check_phylogenetic_data_integrity();
    
private:
    // File path helpers
    std::string get_enhanced_hash_table_path() const;
    std::string get_contributing_taxa_path() const;
    std::string get_enhanced_config_path() const;
    std::string get_phylo_summary_path() const;
    std::string get_species_mapping_path() const;
    
    // Validation helpers
    bool validate_streamlined_metadata(const std::vector<StreamlinedMinimizerMetadata>& metadata);
    bool validate_contributing_taxa_consistency(const ContributingTaxaArrays& taxa_arrays);
};

// Compact binary database serializer for minimal disk usage
class CompactBinarySerializer {
private:
    std::string output_directory_;
    bool use_compression_;
    
public:
    explicit CompactBinarySerializer(const std::string& output_dir, bool compress = true);
    
    // Compact binary format
    bool save_compact_database(
        const std::vector<LCACandidate>& candidates,
        const std::unordered_map<uint32_t, std::string>& taxon_names,
        const DatabaseMetadata& metadata
    );
    
    // Load compact database (for verification)
    bool load_compact_database(
        std::vector<LCACandidate>& candidates,
        std::unordered_map<uint32_t, std::string>& taxon_names,
        DatabaseMetadata& metadata
    );
    
    // Compression utilities
    bool compress_database_files();
    bool decompress_database_files();
    
    // Size estimation
    size_t estimate_compressed_size(const std::vector<LCACandidate>& candidates);
    
private:
    // Binary format helpers
    bool write_binary_header(std::ofstream& out, const DatabaseMetadata& metadata);
    bool read_binary_header(std::ifstream& in, DatabaseMetadata& metadata);
    bool write_candidates_binary(std::ofstream& out, const std::vector<LCACandidate>& candidates);
    bool read_candidates_binary(std::ifstream& in, std::vector<LCACandidate>& candidates);
};

// Database format converter utility
class DatabaseFormatConverter {
private:
    std::string input_directory_;
    std::string output_directory_;
    
public:
    DatabaseFormatConverter(const std::string& input_dir, const std::string& output_dir);
    
    // Format conversion methods
    bool convert_standard_to_enhanced(const std::string& taxonomy_nodes_path = "",
                                     const std::string& taxonomy_names_path = "");
    bool convert_enhanced_to_standard();
    bool convert_to_compact_binary();
    bool convert_to_streaming_friendly();
    
    // Detection methods
    DatabaseFormat detect_database_format(const std::string& database_dir);
    bool is_valid_database_directory(const std::string& database_dir);
    
    // Migration utilities
    bool migrate_v1_to_v2();
    bool upgrade_database_format();
    
private:
    // Conversion helpers
    bool load_standard_database(std::vector<LCACandidate>& candidates,
                               std::unordered_map<uint32_t, std::string>& taxon_names,
                               std::unordered_map<uint32_t, uint32_t>& taxon_parents);
    bool load_enhanced_database(std::vector<PhylogeneticLCACandidate>& phylo_candidates,
                               ContributingTaxaArrays& contributing_taxa);
    
    bool create_phylogenetic_data_from_standard(const std::vector<LCACandidate>& standard_candidates,
                                               std::vector<PhylogeneticLCACandidate>& phylo_candidates);
};

// Database validation and verification utilities
namespace DatabaseValidation {
    // Integrity checks
    bool validate_hash_table_integrity(const std::string& hash_table_path);
    bool validate_taxonomy_consistency(const std::string& taxonomy_path);
    bool check_database_completeness(const std::string& database_dir);
    
    // Cross-validation
    bool cross_validate_hash_taxonomy(const std::string& hash_table_path, 
                                     const std::string& taxonomy_path);
    bool verify_minimizer_uniqueness(const std::string& hash_table_path);
    
    // Performance validation
    bool benchmark_database_loading(const std::string& database_dir);
    bool test_query_performance(const std::string& database_dir, 
                               const std::string& test_sequences_path);
    
    // Checksum utilities
    std::string calculate_file_checksum(const std::string& file_path);
    bool verify_database_checksums(const std::string& database_dir);
    bool generate_manifest_file(const std::string& database_dir);
}

// Database compression and optimization utilities
namespace DatabaseOptimization {
    // Size optimization
    bool optimize_hash_table_layout(const std::string& database_dir);
    bool compress_taxonomy_data(const std::string& database_dir);
    bool remove_redundant_entries(const std::string& database_dir);
    
    // Access optimization
    bool create_index_files(const std::string& database_dir);
    bool optimize_for_sequential_access(const std::string& database_dir);
    bool create_memory_mapped_version(const std::string& database_dir);
    
    // Statistics and analysis
    size_t calculate_database_size(const std::string& database_dir);
    double estimate_memory_usage(const std::string& database_dir);
    void print_database_statistics(const std::string& database_dir);
}

#endif // DATABASE_SERIALIZER_H
