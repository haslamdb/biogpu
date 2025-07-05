// runtime/common/config/unified_config.cpp
#include "unified_config.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>
#include "../../third_party/json.hpp"  // Using nlohmann/json

using json = nlohmann::json;

namespace BioGPU {

// Initialize singleton instance
ConfigManager* ConfigManager::instance = nullptr;

// Helper function to check if file exists
static bool fileExists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// Helper function to expand environment variables
static std::string expandPath(const std::string& path) {
    if (path.empty() || path[0] != '$') {
        return path;
    }
    
    size_t end = path.find_first_of("/\\");
    std::string var_name = path.substr(1, end - 1);
    const char* value = std::getenv(var_name.c_str());
    
    if (value) {
        return std::string(value) + (end != std::string::npos ? path.substr(end) : "");
    }
    return path;
}

std::unique_ptr<UnifiedConfig> ConfigLoader::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open config file: " << filename << std::endl;
        return nullptr;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    return loadFromJSON(buffer.str());
}

std::unique_ptr<UnifiedConfig> ConfigLoader::loadFromJSON(const std::string& json_str) {
    try {
        json j = json::parse(json_str);
        auto config = std::make_unique<UnifiedConfig>();
        
        // Parse each section
        parseGeneralSettings(&j, *config);
        
        if (j.contains("resistance")) {
            parseResistanceConfig(&j["resistance"], config->resistance_config);
        }
        
        if (j.contains("genes")) {
            parseGenesConfig(&j["genes"], config->genes_config);
        }
        
        if (j.contains("profiler")) {
            parseProfilerConfig(&j["profiler"], config->profiler_config);
        }
        
        if (j.contains("io")) {
            parseIOConfig(&j["io"], config->io_config);
        }
        
        if (j.contains("resources")) {
            parseResourceConfig(&j["resources"], config->resource_config);
        }
        
        if (j.contains("logging")) {
            parseLogConfig(&j["logging"], config->log_config);
        }
        
        return config;
        
    } catch (const json::exception& e) {
        std::cerr << "Error parsing JSON config: " << e.what() << std::endl;
        return nullptr;
    }
}

bool ConfigLoader::saveToFile(const UnifiedConfig& config, const std::string& filename) {
    json j;
    
    // General pipeline settings
    j["pipeline"]["batch_size"] = config.pipeline_config.batch_size;
    j["pipeline"]["max_gpu_sequences"] = config.pipeline_config.max_gpu_sequences;
    j["pipeline"]["max_gpu_bases"] = config.pipeline_config.max_gpu_bases;
    j["pipeline"]["num_gpu_streams"] = config.pipeline_config.num_gpu_streams;
    j["pipeline"]["enable_resistance"] = config.pipeline_config.enable_resistance;
    j["pipeline"]["enable_genes"] = config.pipeline_config.enable_genes;
    j["pipeline"]["enable_profiler"] = config.pipeline_config.enable_profiler;
    j["pipeline"]["output_dir"] = config.pipeline_config.output_dir;
    
    // Resistance config
    j["resistance"]["card_database_path"] = config.resistance_config.card_database_path;
    j["resistance"]["min_coverage"] = config.resistance_config.min_coverage;
    j["resistance"]["min_identity"] = config.resistance_config.min_identity;
    j["resistance"]["min_kmer_hits"] = config.resistance_config.min_kmer_hits;
    j["resistance"]["batch_size"] = config.resistance_config.batch_size;
    
    // Three-stage pipeline parameters
    j["resistance"]["use_bloom_filter"] = config.resistance_config.use_bloom_filter;
    j["resistance"]["bloom_kmer_size"] = config.resistance_config.bloom_kmer_size;
    j["resistance"]["nucleotide_kmer_size"] = config.resistance_config.nucleotide_kmer_size;
    j["resistance"]["protein_kmer_size"] = config.resistance_config.protein_kmer_size;
    j["resistance"]["enable_bloom_filter"] = config.resistance_config.enable_bloom_filter;
    j["resistance"]["bloom_filter_size"] = config.resistance_config.bloom_filter_size;
    j["resistance"]["bloom_filter_hashes"] = config.resistance_config.bloom_filter_hashes;
    
    // Output options
    j["resistance"]["generate_html_report"] = config.resistance_config.generate_html_report;
    j["resistance"]["generate_json_report"] = config.resistance_config.generate_json_report;
    j["resistance"]["generate_allele_csv"] = config.resistance_config.generate_allele_csv;
    j["resistance"]["output_sam_format"] = config.resistance_config.output_sam_format;
    j["resistance"]["generate_clinical_report"] = config.resistance_config.generate_clinical_report;
    
    // QRDR targets
    j["resistance"]["target_genes"] = config.resistance_config.target_genes;
    j["resistance"]["target_species"] = config.resistance_config.target_species;
    
    // Thresholds
    j["resistance"]["high_confidence_threshold"] = config.resistance_config.high_confidence_threshold;
    j["resistance"]["calculate_allele_frequencies"] = config.resistance_config.calculate_allele_frequencies;
    j["resistance"]["min_coverage_for_frequency"] = config.resistance_config.min_coverage_for_frequency;
    j["resistance"]["min_frequency_to_report"] = config.resistance_config.min_frequency_to_report;
    
    // Genes config
    j["genes"]["amr_cds_path"] = config.genes_config.amr_cds_path;
    j["genes"]["amr_prot_path"] = config.genes_config.amr_prot_path;
    j["genes"]["minimizer_k"] = config.genes_config.minimizer_k;
    j["genes"]["minimizer_w"] = config.genes_config.minimizer_w;
    j["genes"]["min_gene_coverage"] = config.genes_config.min_gene_coverage;
    j["genes"]["min_gene_identity"] = config.genes_config.min_gene_identity;
    
    // Profiler config
    j["profiler"]["taxonomy_db_path"] = config.profiler_config.taxonomy_db_path;
    j["profiler"]["sketch_size"] = config.profiler_config.sketch_size;
    j["profiler"]["kmer_size"] = config.profiler_config.kmer_size;
    j["profiler"]["confidence_threshold"] = config.profiler_config.confidence_threshold;
    
    // IO config
    j["io"]["input_format"] = config.io_config.input_format;
    j["io"]["paired_end"] = config.io_config.paired_end;
    j["io"]["quality_threshold"] = config.io_config.quality_threshold;
    j["io"]["output_format"] = config.io_config.output_format;
    
    // Resource config
    j["resources"]["max_memory_gb"] = config.resource_config.max_memory_gb;
    j["resources"]["max_threads"] = config.resource_config.max_threads;
    j["resources"]["gpu_device_id"] = config.resource_config.gpu_device_id;
    
    // Logging config
    j["logging"]["log_level"] = config.log_config.log_level;
    j["logging"]["log_file"] = config.log_config.log_file;
    j["logging"]["log_to_console"] = config.log_config.log_to_console;
    
    // Write to file
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot write to config file: " << filename << std::endl;
        return false;
    }
    
    file << j.dump(4);  // Pretty print with 4 spaces
    return true;
}

void ConfigLoader::parseGeneralSettings(const void* json_obj, UnifiedConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    if (j.contains("pipeline")) {
        const auto& p = j["pipeline"];
        
        if (p.contains("batch_size")) 
            config.pipeline_config.batch_size = p["batch_size"];
        if (p.contains("max_gpu_sequences")) 
            config.pipeline_config.max_gpu_sequences = p["max_gpu_sequences"];
        if (p.contains("max_gpu_bases")) 
            config.pipeline_config.max_gpu_bases = p["max_gpu_bases"];
        if (p.contains("num_gpu_streams")) 
            config.pipeline_config.num_gpu_streams = p["num_gpu_streams"];
        if (p.contains("enable_resistance")) 
            config.pipeline_config.enable_resistance = p["enable_resistance"];
        if (p.contains("enable_genes")) 
            config.pipeline_config.enable_genes = p["enable_genes"];
        if (p.contains("enable_profiler")) 
            config.pipeline_config.enable_profiler = p["enable_profiler"];
        if (p.contains("output_dir")) 
            config.pipeline_config.output_dir = expandPath(p["output_dir"]);
        if (p.contains("use_pinned_memory")) 
            config.pipeline_config.use_pinned_memory = p["use_pinned_memory"];
        if (p.contains("num_worker_threads")) 
            config.pipeline_config.num_worker_threads = p["num_worker_threads"];
    }
}

void ConfigLoader::parseResistanceConfig(const void* json_obj, ResistanceConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    // Database paths
    if (j.contains("card_database_path"))
        config.card_database_path = expandPath(j["card_database_path"]);
    if (j.contains("megares_database_path"))
        config.megares_database_path = expandPath(j["megares_database_path"]);
    if (j.contains("custom_database_path"))
        config.custom_database_path = expandPath(j["custom_database_path"]);
    
    // Detection parameters
    if (j.contains("min_coverage"))
        config.min_coverage = j["min_coverage"];
    if (j.contains("min_identity"))
        config.min_identity = j["min_identity"];
    if (j.contains("min_alignment_length"))
        config.min_alignment_length = j["min_alignment_length"];
    if (j.contains("min_kmer_hits"))
        config.min_kmer_hits = j["min_kmer_hits"];
    if (j.contains("batch_size"))
        config.batch_size = j["batch_size"];
    
    // Three-stage pipeline parameters
    if (j.contains("use_bloom_filter"))
        config.use_bloom_filter = j["use_bloom_filter"];
    if (j.contains("bloom_kmer_size"))
        config.bloom_kmer_size = j["bloom_kmer_size"];
    if (j.contains("nucleotide_kmer_size"))
        config.nucleotide_kmer_size = j["nucleotide_kmer_size"];
    if (j.contains("protein_kmer_size"))
        config.protein_kmer_size = j["protein_kmer_size"];
    
    // Bloom filter
    if (j.contains("enable_bloom_filter"))
        config.enable_bloom_filter = j["enable_bloom_filter"];
    if (j.contains("bloom_filter_size"))
        config.bloom_filter_size = j["bloom_filter_size"];
    if (j.contains("bloom_filter_hashes"))
        config.bloom_filter_hashes = j["bloom_filter_hashes"];
    
    // Output options
    if (j.contains("generate_html_report"))
        config.generate_html_report = j["generate_html_report"];
    if (j.contains("generate_json_report"))
        config.generate_json_report = j["generate_json_report"];
    if (j.contains("generate_allele_csv"))
        config.generate_allele_csv = j["generate_allele_csv"];
    if (j.contains("output_sam_format"))
        config.output_sam_format = j["output_sam_format"];
    if (j.contains("generate_clinical_report"))
        config.generate_clinical_report = j["generate_clinical_report"];
    
    // QRDR regions
    if (j.contains("target_genes")) {
        config.target_genes.clear();
        for (const auto& gene : j["target_genes"]) {
            config.target_genes.push_back(gene);
        }
    }
    if (j.contains("target_species")) {
        config.target_species.clear();
        for (const auto& species : j["target_species"]) {
            config.target_species.push_back(species);
        }
    }
    
    // Reporting thresholds
    if (j.contains("report_partial_matches"))
        config.report_partial_matches = j["report_partial_matches"];
    if (j.contains("partial_match_threshold"))
        config.partial_match_threshold = j["partial_match_threshold"];
    if (j.contains("high_confidence_threshold"))
        config.high_confidence_threshold = j["high_confidence_threshold"];
    
    // Allele frequency analysis
    if (j.contains("calculate_allele_frequencies"))
        config.calculate_allele_frequencies = j["calculate_allele_frequencies"];
    if (j.contains("min_coverage_for_frequency"))
        config.min_coverage_for_frequency = j["min_coverage_for_frequency"];
    if (j.contains("min_frequency_to_report"))
        config.min_frequency_to_report = j["min_frequency_to_report"];
    
    // Base class parameters
    if (j.contains("max_sequences_per_batch"))
        config.max_sequences_per_batch = j["max_sequences_per_batch"];
    if (j.contains("verbose"))
        config.verbose = j["verbose"];
}

void ConfigLoader::parseGenesConfig(const void* json_obj, GenesConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    // Database paths
    if (j.contains("amr_cds_path"))
        config.amr_cds_path = expandPath(j["amr_cds_path"]);
    if (j.contains("amr_prot_path"))
        config.amr_prot_path = expandPath(j["amr_prot_path"]);
    if (j.contains("gene_catalog_path"))
        config.gene_catalog_path = expandPath(j["gene_catalog_path"]);
    
    // Algorithm parameters
    if (j.contains("minimizer_k"))
        config.minimizer_k = j["minimizer_k"];
    if (j.contains("minimizer_w"))
        config.minimizer_w = j["minimizer_w"];
    if (j.contains("max_mismatches"))
        config.max_mismatches = j["max_mismatches"];
    if (j.contains("min_exact_matches"))
        config.min_exact_matches = j["min_exact_matches"];
    
    // Translation settings
    if (j.contains("genetic_code_table"))
        config.genetic_code_table = j["genetic_code_table"];
    if (j.contains("translate_all_frames"))
        config.translate_all_frames = j["translate_all_frames"];
    
    // Reporting
    if (j.contains("min_gene_coverage"))
        config.min_gene_coverage = j["min_gene_coverage"];
    if (j.contains("min_gene_identity"))
        config.min_gene_identity = j["min_gene_identity"];
}

void ConfigLoader::parseProfilerConfig(const void* json_obj, ProfilerConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    // Database paths
    if (j.contains("taxonomy_db_path"))
        config.taxonomy_db_path = expandPath(j["taxonomy_db_path"]);
    if (j.contains("reference_genomes_path"))
        config.reference_genomes_path = expandPath(j["reference_genomes_path"]);
    
    // Algorithm parameters
    if (j.contains("sketch_size"))
        config.sketch_size = j["sketch_size"];
    if (j.contains("kmer_size"))
        config.kmer_size = j["kmer_size"];
    if (j.contains("min_unique_kmers"))
        config.min_unique_kmers = j["min_unique_kmers"];
    
    // Classification
    if (j.contains("confidence_threshold"))
        config.confidence_threshold = j["confidence_threshold"];
    if (j.contains("min_reads_for_species"))
        config.min_reads_for_species = j["min_reads_for_species"];
    if (j.contains("report_unclassified"))
        config.report_unclassified = j["report_unclassified"];
    
    // Abundance
    if (j.contains("estimate_abundance"))
        config.estimate_abundance = j["estimate_abundance"];
    if (j.contains("abundance_method"))
        config.abundance_method = j["abundance_method"];
}

void ConfigLoader::parseIOConfig(const void* json_obj, UnifiedConfig::IOConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    if (j.contains("input_format"))
        config.input_format = j["input_format"];
    if (j.contains("paired_end"))
        config.paired_end = j["paired_end"];
    if (j.contains("quality_threshold"))
        config.quality_threshold = j["quality_threshold"];
    if (j.contains("min_read_length"))
        config.min_read_length = j["min_read_length"];
    if (j.contains("output_format"))
        config.output_format = j["output_format"];
    if (j.contains("compress_output"))
        config.compress_output = j["compress_output"];
}

void ConfigLoader::parseResourceConfig(const void* json_obj, UnifiedConfig::ResourceConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    if (j.contains("max_memory_gb"))
        config.max_memory_gb = j["max_memory_gb"];
    if (j.contains("max_threads"))
        config.max_threads = j["max_threads"];
    if (j.contains("gpu_device_id"))
        config.gpu_device_id = j["gpu_device_id"];
    if (j.contains("gpu_memory_limit_gb"))
        config.gpu_memory_limit_gb = j["gpu_memory_limit_gb"];
}

void ConfigLoader::parseLogConfig(const void* json_obj, UnifiedConfig::LogConfig& config) {
    const json& j = *static_cast<const json*>(json_obj);
    
    if (j.contains("log_level"))
        config.log_level = j["log_level"];
    if (j.contains("log_file"))
        config.log_file = expandPath(j["log_file"]);
    if (j.contains("log_to_console"))
        config.log_to_console = j["log_to_console"];
    if (j.contains("save_intermediate_files"))
        config.save_intermediate_files = j["save_intermediate_files"];
    if (j.contains("temp_directory"))
        config.temp_directory = expandPath(j["temp_directory"]);
}

bool ConfigLoader::validate(const UnifiedConfig& config, std::vector<std::string>& errors) {
    errors.clear();
    
    bool valid = true;
    valid &= validatePaths(config, errors);
    valid &= validateParameters(config, errors);
    valid &= validateResources(config, errors);
    
    return valid;
}

bool ConfigLoader::validatePaths(const UnifiedConfig& config, std::vector<std::string>& errors) {
    bool valid = true;
    
    // Check resistance database paths if pipeline is enabled
    if (config.pipeline_config.enable_resistance) {
        if (!config.resistance_config.card_database_path.empty() && 
            !fileExists(config.resistance_config.card_database_path)) {
            errors.push_back("CARD database not found: " + config.resistance_config.card_database_path);
            valid = false;
        }
    }
    
    // Check genes database paths if pipeline is enabled
    if (config.pipeline_config.enable_genes) {
        if (!fileExists(config.genes_config.amr_cds_path)) {
            errors.push_back("AMR CDS file not found: " + config.genes_config.amr_cds_path);
            valid = false;
        }
        if (!fileExists(config.genes_config.amr_prot_path)) {
            errors.push_back("AMR protein file not found: " + config.genes_config.amr_prot_path);
            valid = false;
        }
    }
    
    // Check profiler database paths if pipeline is enabled
    if (config.pipeline_config.enable_profiler) {
        if (!config.profiler_config.taxonomy_db_path.empty() && 
            !fileExists(config.profiler_config.taxonomy_db_path)) {
            errors.push_back("Taxonomy database not found: " + config.profiler_config.taxonomy_db_path);
            valid = false;
        }
    }
    
    return valid;
}

bool ConfigLoader::validateParameters(const UnifiedConfig& config, std::vector<std::string>& errors) {
    bool valid = true;
    
    // Validate batch sizes
    if (config.pipeline_config.batch_size == 0) {
        errors.push_back("Batch size must be greater than 0");
        valid = false;
    }
    
    // Validate resistance parameters
    if (config.resistance_config.min_coverage < 0 || config.resistance_config.min_coverage > 1) {
        errors.push_back("Resistance min_coverage must be between 0 and 1");
        valid = false;
    }
    
    // Validate genes parameters
    if (config.genes_config.minimizer_k >= config.genes_config.minimizer_w) {
        errors.push_back("Minimizer k must be less than minimizer w");
        valid = false;
    }
    
    // Validate profiler parameters
    if (config.profiler_config.confidence_threshold < 0 || 
        config.profiler_config.confidence_threshold > 1) {
        errors.push_back("Profiler confidence threshold must be between 0 and 1");
        valid = false;
    }
    
    return valid;
}

bool ConfigLoader::validateResources(const UnifiedConfig& config, std::vector<std::string>& errors) {
    bool valid = true;
    
    // Check GPU device
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    if (config.resource_config.gpu_device_id >= device_count) {
        errors.push_back("Invalid GPU device ID: " + std::to_string(config.resource_config.gpu_device_id));
        valid = false;
    }
    
    // Check thread count
    if (config.resource_config.max_threads < 0) {
        errors.push_back("Max threads cannot be negative");
        valid = false;
    }
    
    return valid;
}

void ConfigLoader::applyOverrides(UnifiedConfig& config, 
                                 const std::map<std::string, std::string>& overrides) {
    for (const auto& [key, value] : overrides) {
        // Parse dot-separated keys like "resistance.min_coverage"
        std::vector<std::string> parts;
        std::stringstream ss(key);
        std::string part;
        while (std::getline(ss, part, '.')) {
            parts.push_back(part);
        }
        
        if (parts.empty()) continue;
        
        // Apply overrides based on key structure
        if (parts[0] == "pipeline" && parts.size() == 2) {
            if (parts[1] == "batch_size") {
                config.pipeline_config.batch_size = std::stoul(value);
            } else if (parts[1] == "output_dir") {
                config.pipeline_config.output_dir = value;
            }
            // Add more pipeline overrides as needed
        } else if (parts[0] == "resistance" && parts.size() == 2) {
            if (parts[1] == "min_coverage") {
                config.resistance_config.min_coverage = std::stof(value);
            } else if (parts[1] == "min_identity") {
                config.resistance_config.min_identity = std::stof(value);
            }
            // Add more resistance overrides as needed
        }
        // Add more sections as needed
    }
}

UnifiedConfig ConfigLoader::getDefault() {
    UnifiedConfig config;
    
    // Set sensible defaults
    config.pipeline_config.batch_size = 100000;
    config.pipeline_config.max_gpu_sequences = 500000;
    config.pipeline_config.max_gpu_bases = 500000000;
    config.pipeline_config.output_dir = "output";
    
    config.resistance_config.min_coverage = 0.80f;
    config.resistance_config.min_identity = 0.95f;
    config.resistance_config.min_kmer_hits = 10;
    config.resistance_config.batch_size = 10000;
    config.resistance_config.use_bloom_filter = true;
    config.resistance_config.bloom_kmer_size = 15;
    config.resistance_config.nucleotide_kmer_size = 31;
    config.resistance_config.protein_kmer_size = 10;
    config.resistance_config.enable_bloom_filter = true;
    config.resistance_config.bloom_filter_size = 1073741824;  // 1GB
    config.resistance_config.bloom_filter_hashes = 7;
    config.resistance_config.generate_html_report = true;
    config.resistance_config.generate_json_report = true;
    config.resistance_config.generate_allele_csv = true;
    config.resistance_config.output_sam_format = false;
    config.resistance_config.generate_clinical_report = true;
    config.resistance_config.high_confidence_threshold = 0.95f;
    config.resistance_config.calculate_allele_frequencies = true;
    config.resistance_config.min_coverage_for_frequency = 10;
    config.resistance_config.min_frequency_to_report = 0.01f;
    
    config.genes_config.minimizer_k = 15;
    config.genes_config.minimizer_w = 10;
    config.genes_config.min_gene_coverage = 0.95f;
    config.genes_config.min_gene_identity = 0.95f;
    
    config.profiler_config.sketch_size = 1000;
    config.profiler_config.kmer_size = 21;
    config.profiler_config.confidence_threshold = 0.8f;
    
    return config;
}

} // namespace BioGPU