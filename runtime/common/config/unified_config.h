// runtime/common/config/unified_config.h
#ifndef UNIFIED_CONFIG_H
#define UNIFIED_CONFIG_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include "../pipeline/pipeline_base.h"
#include "../pipeline/pipeline_coordinator.h"

namespace BioGPU {

// Resistance pipeline specific config
struct ResistanceConfig : public PipelineBaseConfig {
    // Database settings
    std::string card_database_path = "";
    std::string megares_database_path = "";
    std::string custom_database_path = "";
    
    // Detection parameters
    float min_identity = 0.95f;
    float min_coverage = 0.80f;
    int min_alignment_length = 50;
    int min_kmer_hits = 10;
    int batch_size = 10000;  // Keep existing batch size for compatibility
    
    // Three-stage pipeline parameters
    bool use_bloom_filter = true;
    int bloom_kmer_size = 15;      // Stage 1: Bloom filter k-mer size
    int nucleotide_kmer_size = 31; // Stage 2: Nucleotide matching k-mer size
    int protein_kmer_size = 10;    // Stage 3: Protein search k-mer size
    
    // Filtering
    bool enable_bloom_filter = true;
    size_t bloom_filter_size = 1073741824;  // 1GB
    int bloom_filter_hashes = 7;
    
    // Output options
    bool generate_html_report = true;
    bool generate_json_report = true;
    bool generate_allele_csv = true;
    bool output_sam_format = false;
    bool generate_clinical_report = true;
    
    // QRDR regions to check
    std::vector<std::string> target_genes = {"gyrA", "gyrB", "parC", "parE"};
    std::vector<std::string> target_species = {"Escherichia coli", "Klebsiella pneumoniae", 
                                               "Salmonella enterica", "Staphylococcus aureus"};
    
    // Reporting thresholds
    bool report_partial_matches = false;
    float partial_match_threshold = 0.8f;
    float high_confidence_threshold = 0.95f;
    
    // Allele frequency analysis
    bool calculate_allele_frequencies = true;
    int min_coverage_for_frequency = 10;  // Minimum reads to calculate frequency
    float min_frequency_to_report = 0.01f; // 1% minimum frequency
};

// Genes pipeline specific config
struct GenesConfig : public PipelineBaseConfig {
    // Database settings
    std::string amr_cds_path = "AMR_CDS.fa";
    std::string amr_prot_path = "AMRProt.fa";
    std::string gene_catalog_path = "ReferenceGeneCatalog.txt";
    
    // Algorithm parameters
    int minimizer_k = 15;
    int minimizer_w = 10;
    int max_mismatches = 2;
    int min_exact_matches = 5;
    
    // Translation settings
    int genetic_code_table = 11;  // Bacterial code
    bool translate_all_frames = true;
    
    // Reporting
    float min_gene_coverage = 0.95f;
    float min_gene_identity = 0.95f;
};

// Profiler pipeline specific config
struct ProfilerConfig : public PipelineBaseConfig {
    // Database settings
    std::string taxonomy_db_path = "";
    std::string reference_genomes_path = "";
    
    // Algorithm parameters
    int sketch_size = 1000;
    int kmer_size = 21;
    float min_unique_kmers = 0.1f;
    
    // Classification
    float confidence_threshold = 0.8f;
    int min_reads_for_species = 10;
    bool report_unclassified = false;
    
    // Abundance estimation
    bool estimate_abundance = true;
    std::string abundance_method = "relative";  // "relative" or "absolute"
};

// Master configuration containing all settings
struct UnifiedConfig {
    // General settings
    PipelineConfig pipeline_config;
    
    // Pipeline-specific configs
    ResistanceConfig resistance_config;
    GenesConfig genes_config;
    ProfilerConfig profiler_config;
    
    // Input/Output settings
    struct IOConfig {
        std::string input_format = "fastq";  // "fastq", "fasta"
        bool paired_end = false;
        int quality_threshold = 20;
        int min_read_length = 50;
        std::string output_format = "tsv";  // "tsv", "json", "xml"
        bool compress_output = false;
    } io_config;
    
    // Resource limits
    struct ResourceConfig {
        size_t max_memory_gb = 0;  // 0 = auto
        int max_threads = 0;  // 0 = auto
        int gpu_device_id = 0;
        size_t gpu_memory_limit_gb = 0;  // 0 = auto
    } resource_config;
    
    // Logging and debugging
    struct LogConfig {
        std::string log_level = "info";  // "debug", "info", "warning", "error"
        std::string log_file = "";
        bool log_to_console = true;
        bool save_intermediate_files = false;
        std::string temp_directory = "/tmp";
    } log_config;
};

// Configuration loader class
class ConfigLoader {
public:
    // Load configuration from file
    static std::unique_ptr<UnifiedConfig> loadFromFile(const std::string& filename);
    
    // Load configuration from JSON string
    static std::unique_ptr<UnifiedConfig> loadFromJSON(const std::string& json_str);
    
    // Save configuration to file
    static bool saveToFile(const UnifiedConfig& config, const std::string& filename);
    
    // Validate configuration
    static bool validate(const UnifiedConfig& config, std::vector<std::string>& errors);
    
    // Apply command-line overrides
    static void applyOverrides(UnifiedConfig& config, 
                               const std::map<std::string, std::string>& overrides);
    
    // Get default configuration
    static UnifiedConfig getDefault();
    
private:
    // Helper methods for parsing
    static void parseGeneralSettings(const void* json_obj, UnifiedConfig& config);
    static void parseResistanceConfig(const void* json_obj, ResistanceConfig& config);
    static void parseGenesConfig(const void* json_obj, GenesConfig& config);
    static void parseProfilerConfig(const void* json_obj, ProfilerConfig& config);
    static void parseIOConfig(const void* json_obj, UnifiedConfig::IOConfig& config);
    static void parseResourceConfig(const void* json_obj, UnifiedConfig::ResourceConfig& config);
    static void parseLogConfig(const void* json_obj, UnifiedConfig::LogConfig& config);
    
    // Validation helpers
    static bool validatePaths(const UnifiedConfig& config, std::vector<std::string>& errors);
    static bool validateParameters(const UnifiedConfig& config, std::vector<std::string>& errors);
    static bool validateResources(const UnifiedConfig& config, std::vector<std::string>& errors);
};

// Configuration manager singleton
class ConfigManager {
private:
    static ConfigManager* instance;
    std::unique_ptr<UnifiedConfig> config;
    std::string config_file_path;
    
    ConfigManager() = default;
    
public:
    static ConfigManager& getInstance() {
        if (!instance) {
            instance = new ConfigManager();
        }
        return *instance;
    }
    
    // Load configuration
    bool loadConfig(const std::string& filename) {
        config = ConfigLoader::loadFromFile(filename);
        if (config) {
            config_file_path = filename;
            return true;
        }
        return false;
    }
    
    // Get configuration
    UnifiedConfig* getConfig() { return config.get(); }
    const UnifiedConfig* getConfig() const { return config.get(); }
    
    // Update configuration
    void setConfig(std::unique_ptr<UnifiedConfig> new_config) {
        config = std::move(new_config);
    }
    
    // Reload configuration
    bool reload() {
        if (!config_file_path.empty()) {
            return loadConfig(config_file_path);
        }
        return false;
    }
    
    // Get specific pipeline configs
    const ResistanceConfig& getResistanceConfig() const { 
        return config->resistance_config; 
    }
    
    const GenesConfig& getGenesConfig() const { 
        return config->genes_config; 
    }
    
    const ProfilerConfig& getProfilerConfig() const { 
        return config->profiler_config; 
    }
};

} // namespace BioGPU

#endif // UNIFIED_CONFIG_H