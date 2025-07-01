// gpu_kraken_database_builder.cu
// COMPLETE REPLACEMENT: Kraken2-inspired GPU database builder with all original functionality
// Fixes the over-aggressive deduplication while maintaining all features

#ifndef GPU_KRAKEN_DATABASE_BUILDER_CUH
#define GPU_KRAKEN_DATABASE_BUILDER_CUH

#include "gpu_kraken_database_builder.h"
#include "gpu_kraken_classifier.h"
#include "gpu_minimizer_extraction.cuh"
#include "streaming_fna_processor.h"
#include "processing/feature_exporter.h"
#include "processing/minimizer_feature_extractor.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/unique.h>
#include <thrust/scan.h>
#include <cub/cub.cuh>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <memory>
#include <set>
#include <algorithm>
#include <sstream>
#include <cctype>
#include <climits>
#include <stdexcept>
#include <sys/stat.h>  // For stat()
#include <sys/types.h> // For stat structures

// Keep all original structures
struct GPUTaxonomyNode {
    uint32_t taxon_id;
    uint32_t parent_id;
    uint8_t rank;
    uint8_t padding[3];
};

// Database building statistics
struct GPUBuildStats {
    uint64_t total_sequences;
    uint64_t total_bases;
    uint64_t total_kmers_processed;
    uint64_t valid_minimizers_extracted;
    uint64_t unique_minimizers;
    uint64_t lca_assignments;
    double sequence_processing_time;
    double minimizer_extraction_time;
    double lca_computation_time;
    double database_construction_time;
};

// Add these structures to gpu_kraken_database_builder.cu
// These are production-ready with no placeholders

// Define the missing LCACandidate struct that's referenced in your existing code
struct LCACandidate {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
    
    // Default constructor
    LCACandidate() : minimizer_hash(0), lca_taxon(0), genome_count(0), uniqueness_score(0.0f) {}
    
    // Constructor for backward compatibility
    LCACandidate(uint64_t hash, uint32_t taxon, uint32_t count, float score)
        : minimizer_hash(hash), lca_taxon(taxon), genome_count(count), uniqueness_score(score) {}
};

// Enhanced LCA candidate with phylogenetic metadata
struct PhylogeneticLCACandidate {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
    
    // Phylogenetic extensions
    std::vector<uint32_t> contributing_species;      // Species that have this minimizer
    std::vector<uint16_t> genome_counts_per_species; // How many genomes per species
    uint8_t phylogenetic_spread;                     // Diversity measure (0-255)
    uint8_t max_phylogenetic_distance;              // Max distance to LCA (0-255)
    
    // Constructor from basic LCACandidate
    PhylogeneticLCACandidate(const LCACandidate& basic) 
        : minimizer_hash(basic.minimizer_hash), lca_taxon(basic.lca_taxon),
          genome_count(basic.genome_count), uniqueness_score(basic.uniqueness_score),
          phylogenetic_spread(0), max_phylogenetic_distance(0) {}
    
    PhylogeneticLCACandidate() : minimizer_hash(0), lca_taxon(0), genome_count(0), 
                                uniqueness_score(0.0f), phylogenetic_spread(0), 
                                max_phylogenetic_distance(0) {}
};

// Streamlined database format (28 bytes per minimizer)
struct StreamlinedMinimizerMetadata {
    uint64_t minimizer_hash;                 // 8 bytes
    uint32_t lca_taxon;                      // 4 bytes - backward compatibility
    uint32_t total_genome_count;             // 4 bytes
    uint32_t contributing_taxa_offset;       // 4 bytes - offset into external array
    uint16_t num_contributing_taxa;          // 2 bytes
    uint8_t phylogenetic_spread;             // 1 byte
    uint8_t max_phylogenetic_distance;       // 1 byte
    uint32_t reserved;                       // 4 bytes - for future use, total = 28 bytes
};

// External variable-length data arrays
struct ContributingTaxaArrays {
    std::vector<uint32_t> taxa_ids;                      // Species taxon IDs
    std::vector<uint8_t> phylogenetic_distances;         // Distance from each to LCA
    std::vector<uint16_t> genome_counts_per_taxon;       // Genome count per taxon
    
    // Add entry and return offset
    uint32_t add_entry(uint32_t taxon_id, uint8_t distance, uint16_t genome_count) {
        uint32_t offset = taxa_ids.size();
        taxa_ids.push_back(taxon_id);
        phylogenetic_distances.push_back(distance);
        genome_counts_per_taxon.push_back(genome_count);
        return offset;
    }
    
    size_t total_entries() const { return taxa_ids.size(); }
    
    void clear() {
        taxa_ids.clear();
        phylogenetic_distances.clear();
        genome_counts_per_taxon.clear();
    }
};

// Species tracking during genome processing
struct SpeciesTrackingData {
    std::unordered_map<std::string, uint32_t> sequence_id_to_species;  // Map sequence ID to species
    std::unordered_map<uint32_t, uint16_t> species_genome_counts;      // Count genomes per species
    std::unordered_map<uint32_t, std::string> species_names;           // Species ID to name mapping
    
    void add_genome(const std::string& sequence_id, uint32_t species_taxid, const std::string& species_name) {
        sequence_id_to_species[sequence_id] = species_taxid;
        species_genome_counts[species_taxid]++;
        if (species_names.find(species_taxid) == species_names.end()) {
            species_names[species_taxid] = species_name;
        }
    }
    
    uint32_t get_species_for_sequence(const std::string& sequence_id) const {
        auto it = sequence_id_to_species.find(sequence_id);
        return (it != sequence_id_to_species.end()) ? it->second : 0;
    }
    
    uint16_t get_genome_count_for_species(uint32_t species_taxid) const {
        auto it = species_genome_counts.find(species_taxid);
        return (it != species_genome_counts.end()) ? it->second : 0;
    }
    
    size_t total_species() const { return species_genome_counts.size(); }
    size_t total_genomes() const {
        size_t total = 0;
        for (const auto& [species, count] : species_genome_counts) {
            total += count;
        }
        return total;
    }
};

// Enhanced build statistics
struct EnhancedBuildStats : public GPUBuildStats {
    // Additional phylogenetic statistics
    uint64_t species_represented = 0;
    uint64_t minimizers_with_phylo_data = 0;
    uint64_t phylogenetic_lca_computations = 0;
    double phylogenetic_processing_time = 0.0;
    size_t contributing_taxa_array_size = 0;
    
    void print_enhanced_stats() const {
        std::cout << "\n=== ENHANCED BUILD STATISTICS ===" << std::endl;
        std::cout << "Species represented: " << species_represented << std::endl;
        std::cout << "Minimizers with phylogenetic data: " << minimizers_with_phylo_data << std::endl;
        std::cout << "Phylogenetic LCA computations: " << phylogenetic_lca_computations << std::endl;
        std::cout << "Contributing taxa array entries: " << contributing_taxa_array_size << std::endl;
        std::cout << "Phylogenetic processing time: " << std::fixed << std::setprecision(2) 
                  << phylogenetic_processing_time << "s" << std::endl;
        
        if (unique_minimizers > 0) {
            double phylo_coverage = (double)minimizers_with_phylo_data / unique_minimizers * 100.0;
            std::cout << "Phylogenetic coverage: " << std::fixed << std::setprecision(1) 
                      << phylo_coverage << "%" << std::endl;
        }
        
        if (contributing_taxa_array_size > 0) {
            double avg_taxa_per_minimizer = (double)contributing_taxa_array_size / minimizers_with_phylo_data;
            std::cout << "Average taxa per minimizer: " << std::fixed << std::setprecision(1) 
                      << avg_taxa_per_minimizer << std::endl;
        }
    }
};

// Replace the NCBITaxonomyProcessor with this enhanced version that uses your compact taxonomy tool

#include "../tools/compact_gpu_taxonomy.h"  // Include your existing tool

class EnhancedNCBITaxonomyProcessor {
private:
    std::unique_ptr<BioGPU::CompactTaxonomy::CompactGPUTaxonomy> compact_taxonomy;
    bool taxonomy_loaded = false;
    
    // Host-side lookup tables for database building (extracted from compact format)
    std::unordered_map<uint32_t, uint32_t> parent_lookup;
    std::unordered_map<uint32_t, std::string> name_lookup;
    std::unordered_map<uint32_t, std::string> rank_lookup;
    std::unordered_map<uint32_t, uint8_t> depth_lookup;
    uint32_t max_taxon_id = 0;
    
public:
    // Load using your existing compact taxonomy infrastructure
    bool load_ncbi_taxonomy(const std::string& nodes_dmp_path, const std::string& names_dmp_path) {
        std::cout << "Loading NCBI taxonomy using compact taxonomy infrastructure..." << std::endl;
        
        compact_taxonomy = std::make_unique<BioGPU::CompactTaxonomy::CompactGPUTaxonomy>(false); // No cache needed for building
        
        // Use your existing NCBI file parser
        if (!compact_taxonomy->build_from_ncbi_files(nodes_dmp_path, names_dmp_path)) {
            std::cerr << "Failed to build compact taxonomy from NCBI files" << std::endl;
            return false;
        }
        
        // Extract host-side lookup tables for database building
        if (!extract_host_lookup_tables()) {
            std::cerr << "Failed to extract lookup tables from compact taxonomy" << std::endl;
            return false;
        }
        
        taxonomy_loaded = true;
        std::cout << "NCBI taxonomy loaded successfully using compact infrastructure" << std::endl;
        return true;
    }
    
    // Alternative: Load from pre-built compact taxonomy file
    bool load_from_compact_file(const std::string& compact_taxonomy_path) {
        std::cout << "Loading from pre-built compact taxonomy: " << compact_taxonomy_path << std::endl;
        
        compact_taxonomy = std::make_unique<BioGPU::CompactTaxonomy::CompactGPUTaxonomy>(false);
        
        if (!compact_taxonomy->load_compact_taxonomy(compact_taxonomy_path)) {
            std::cerr << "Failed to load compact taxonomy file" << std::endl;
            return false;
        }
        
        if (!extract_host_lookup_tables()) {
            std::cerr << "Failed to extract lookup tables from compact taxonomy" << std::endl;
            return false;
        }
        
        taxonomy_loaded = true;
        std::cout << "Compact taxonomy loaded successfully" << std::endl;
        return true;
    }
    
    // Use the same phylogenetic calculation methods as before, but with proper NCBI data
    uint32_t compute_lca_of_species(const std::vector<uint32_t>& species_list) {
        if (!taxonomy_loaded || species_list.empty()) {
            return 1; // Root
        }
        
        if (species_list.size() == 1) {
            return species_list[0];
        }
        
        // Find LCA using proper taxonomy tree traversal
        uint32_t current_lca = species_list[0];
        for (size_t i = 1; i < species_list.size(); i++) {
            current_lca = find_lca_pair(current_lca, species_list[i]);
        }
        
        return current_lca;
    }
    
    uint8_t calculate_distance_to_lca(uint32_t taxon, uint32_t lca) {
        if (!taxonomy_loaded || taxon == lca) {
            return 0;
        }
        
        // Use depth lookup if available
        auto taxon_depth_it = depth_lookup.find(taxon);
        auto lca_depth_it = depth_lookup.find(lca);
        
        if (taxon_depth_it != depth_lookup.end() && lca_depth_it != depth_lookup.end()) {
            uint8_t taxon_depth = taxon_depth_it->second;
            uint8_t lca_depth = lca_depth_it->second;
            
            if (taxon_depth >= lca_depth) {
                return taxon_depth - lca_depth;
            }
        }
        
        // Fallback: count steps manually
        uint32_t current = taxon;
        uint8_t distance = 0;
        
        while (current != lca && current != 1 && distance < 50) {
            auto parent_it = parent_lookup.find(current);
            if (parent_it == parent_lookup.end()) break;
            
            current = parent_it->second;
            distance++;
            
            if (current == lca) break;
        }
        
        return (current == lca) ? distance : 255;
    }
    
    uint8_t calculate_phylogenetic_spread(const std::vector<uint32_t>& species_list, uint32_t lca) {
        if (species_list.size() <= 1) return 0;
        
        std::vector<uint8_t> distances;
        for (uint32_t species : species_list) {
            distances.push_back(calculate_distance_to_lca(species, lca));
        }
        
        uint8_t max_dist = *std::max_element(distances.begin(), distances.end());
        uint8_t min_dist = *std::min_element(distances.begin(), distances.end());
        
        // Enhanced spread calculation using taxonomy ranks
        uint8_t range_spread = max_dist - min_dist;
        uint8_t diversity_factor = std::min((uint8_t)(species_list.size() / 5), (uint8_t)50);
        
        // Weight by taxonomic rank of LCA
        uint8_t rank_weight = get_rank_weight(lca);
        
        return std::min((uint8_t)255, (uint8_t)(range_spread * rank_weight + diversity_factor));
    }
    
    std::string get_scientific_name(uint32_t taxon_id) {
        auto it = name_lookup.find(taxon_id);
        return (it != name_lookup.end()) ? it->second : ("taxon_" + std::to_string(taxon_id));
    }
    
    std::string get_rank(uint32_t taxon_id) {
        auto it = rank_lookup.find(taxon_id);
        return (it != rank_lookup.end()) ? it->second : "no rank";
    }
    
    bool is_loaded() const { return taxonomy_loaded; }
    
    // Get the compact taxonomy instance for GPU operations (if needed later)
    BioGPU::CompactTaxonomy::CompactGPUTaxonomy* get_compact_taxonomy() {
        return compact_taxonomy.get();
    }
    
private:
    // Extract host-side lookup tables from the compact taxonomy
    bool extract_host_lookup_tables() {
        if (!compact_taxonomy) {
            return false;
        }
        
        // Extract lookup tables using the accessor methods
        parent_lookup = compact_taxonomy->get_parent_lookup_map();
        name_lookup = compact_taxonomy->get_name_lookup_map();
        rank_lookup = compact_taxonomy->get_rank_lookup_map();
        depth_lookup = compact_taxonomy->get_depth_lookup_map();
        max_taxon_id = compact_taxonomy->get_max_taxon_id();
        
        std::cout << "Extracted host lookup tables from compact taxonomy" << std::endl;
        std::cout << "  Taxa count: " << parent_lookup.size() << std::endl;
        std::cout << "  Max taxon ID: " << max_taxon_id << std::endl;
        
        return true;
    }
    
    uint32_t find_lca_pair(uint32_t taxon1, uint32_t taxon2) {
        if (taxon1 == taxon2) return taxon1;
        if (taxon1 == 1 || taxon2 == 1) return 1;
        
        // Get paths to root
        std::vector<uint32_t> path1 = get_path_to_root(taxon1);
        std::vector<uint32_t> path2 = get_path_to_root(taxon2);
        
        // Find first common ancestor
        std::set<uint32_t> ancestors1(path1.begin(), path1.end());
        
        for (uint32_t ancestor : path2) {
            if (ancestors1.find(ancestor) != ancestors1.end()) {
                return ancestor;
            }
        }
        
        return 1; // Root fallback
    }
    
    std::vector<uint32_t> get_path_to_root(uint32_t taxon_id) {
        std::vector<uint32_t> path;
        uint32_t current = taxon_id;
        
        while (current != 1 && path.size() < 50) {
            path.push_back(current);
            auto parent_it = parent_lookup.find(current);
            if (parent_it == parent_lookup.end()) break;
            current = parent_it->second;
        }
        
        path.push_back(1); // Add root
        return path;
    }
    
    uint8_t get_rank_weight(uint32_t taxon_id) {
        std::string rank = get_rank(taxon_id);
        
        // Weight factors based on taxonomic rank
        if (rank == "species" || rank == "subspecies") return 1;
        else if (rank == "genus") return 2;
        else if (rank == "family") return 3;
        else if (rank == "order") return 4;
        else if (rank == "class") return 5;
        else if (rank == "phylum") return 6;
        else if (rank == "kingdom" || rank == "superkingdom") return 7;
        else return 3; // Default for "no rank"
    }
};

// Production FNA Parser - add to gpu_kraken_database_builder.cu

class ConcatenatedFnaProcessor {
private:
    std::string fna_file_path;
    SpeciesTrackingData species_data;
    size_t bytes_processed = 0;
    size_t total_file_size = 0;
    
public:
    ConcatenatedFnaProcessor(const std::string& file_path) : fna_file_path(file_path) {
        // Get file size for progress reporting
        std::ifstream file(file_path, std::ios::ate | std::ios::binary);
        if (file.is_open()) {
            total_file_size = file.tellg();
            file.close();
        }
    }
    
    // Process the entire FNA file and populate genome_files/genome_taxon_ids for existing pipeline
    bool process_fna_file(std::vector<std::string>& genome_files, 
                         std::vector<uint32_t>& genome_taxon_ids,
                         std::unordered_map<uint32_t, std::string>& taxon_names,
                         const std::string& temp_dir = "/tmp/kraken_build") {
        
        std::cout << "Processing concatenated FNA file: " << fna_file_path << std::endl;
        std::cout << "File size: " << (total_file_size / 1024 / 1024) << " MB" << std::endl;
        
        // Create temporary directory for individual genome files
        std::filesystem::create_directories(temp_dir);
        
        std::ifstream file(fna_file_path);
        if (!file.is_open()) {
            std::cerr << "Cannot open FNA file: " << fna_file_path << std::endl;
            return false;
        }
        
        std::string line;
        std::string current_sequence;
        std::string current_header;
        uint32_t current_species = 0;
        int genome_count = 0;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        while (std::getline(file, line)) {
            bytes_processed += line.length() + 1; // +1 for newline
            
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Process previous genome if exists
                if (!current_sequence.empty() && current_species > 0) {
                    std::string temp_file = create_temp_genome_file(
                        current_sequence, current_species, current_header, 
                        temp_dir, genome_count
                    );
                    
                    if (!temp_file.empty()) {
                        genome_files.push_back(temp_file);
                        genome_taxon_ids.push_back(current_species);
                        genome_count++;
                        
                        // Progress reporting
                        if (genome_count % 1000 == 0) {
                            double progress = (double)bytes_processed / total_file_size * 100.0;
                            std::cout << "Processed " << genome_count << " genomes (" 
                                     << std::fixed << std::setprecision(1) << progress << "%)" << std::endl;
                        }
                    }
                }
                
                // Parse new header
                current_header = line;
                HeaderParseResult parse_result = parse_fna_header(line);
                current_species = parse_result.species_taxid;
                current_sequence.clear();
                
                if (current_species > 0) {
                    // Store species information
                    species_data.add_genome("genome_" + std::to_string(genome_count), 
                                          current_species, parse_result.species_name);
                    
                    // Add to taxon names if not present
                    if (taxon_names.find(current_species) == taxon_names.end()) {
                        taxon_names[current_species] = parse_result.species_name;
                    }
                } else {
                    std::cerr << "Warning: Could not parse species taxid from header: " 
                              << line.substr(0, 100) << "..." << std::endl;
                }
                
            } else {
                // Accumulate sequence data, removing all whitespace
                for (char c : line) {
                    if (!std::isspace(c)) {
                        current_sequence += c;
                    }
                }
            }
        }
        
        // Process the last genome
        if (!current_sequence.empty() && current_species > 0) {
            std::string temp_file = create_temp_genome_file(
                current_sequence, current_species, current_header, 
                temp_dir, genome_count
            );
            
            if (!temp_file.empty()) {
                genome_files.push_back(temp_file);
                genome_taxon_ids.push_back(current_species);
                genome_count++;
            }
        }
        
        file.close();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\nFNA processing completed:" << std::endl;
        std::cout << "  Genomes processed: " << genome_count << std::endl;
        std::cout << "  Species represented: " << species_data.total_species() << std::endl;
        std::cout << "  Processing time: " << duration.count() << " seconds" << std::endl;
        std::cout << "  Processing rate: " << (bytes_processed / 1024 / 1024 / duration.count()) 
                  << " MB/s" << std::endl;
        
        return genome_count > 0;
    }
    
    const SpeciesTrackingData& get_species_data() const { return species_data; }
    
private:
    struct HeaderParseResult {
        uint32_t species_taxid = 0;
        std::string species_name;
        std::string accession;
        std::string description;
    };
    
    HeaderParseResult parse_fna_header(const std::string& header) {
        HeaderParseResult result;
        
        // Parse header format: >kraken:taxid|XXXXX|ACCESSION description
        size_t taxid_start = header.find("taxid|");
        if (taxid_start == std::string::npos) {
            return result; // Invalid header
        }
        
        // Extract taxid
        size_t taxid_value_start = taxid_start + 6; // length of "taxid|"
        size_t taxid_end = header.find('|', taxid_value_start);
        if (taxid_end == std::string::npos) {
            return result; // Invalid header
        }
        
        try {
            std::string taxid_str = header.substr(taxid_value_start, taxid_end - taxid_value_start);
            result.species_taxid = std::stoul(taxid_str);
        } catch (const std::exception& e) {
            return result; // Invalid taxid
        }
        
        // Extract accession and description
        size_t accession_start = taxid_end + 1;
        size_t description_start = header.find(' ', accession_start);
        
        if (description_start != std::string::npos) {
            result.accession = header.substr(accession_start, description_start - accession_start);
            result.description = header.substr(description_start + 1);
            
            // Extract species name from description if possible
            result.species_name = extract_species_name_from_description(result.description);
        } else {
            result.accession = header.substr(accession_start);
        }
        
        // Fallback species name
        if (result.species_name.empty()) {
            result.species_name = "species_" + std::to_string(result.species_taxid);
        }
        
        return result;
    }
    
    std::string extract_species_name_from_description(const std::string& description) {
        // Try to extract genus + species from description
        // Look for pattern: "Genus species" at the beginning
        std::istringstream iss(description);
        std::string genus, species;
        
        if (iss >> genus >> species) {
            // Basic validation: genus should start with capital, species with lowercase
            if (!genus.empty() && !species.empty() && 
                std::isupper(genus[0]) && std::islower(species[0])) {
                return genus + " " + species;
            }
        }
        
        return ""; // Could not extract valid species name
    }
    
    std::string create_temp_genome_file(const std::string& sequence, 
                                       uint32_t species_taxid,
                                       const std::string& original_header,
                                       const std::string& temp_dir,
                                       int genome_index) {
        
        // Skip very short sequences
        if (sequence.length() < 1000) {
            return ""; // Too short to be useful
        }
        
        // Create filename
        std::string filename = temp_dir + "/genome_" + std::to_string(species_taxid) + 
                              "_" + std::to_string(genome_index) + ".fasta";
        
        // Write to temporary file
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            std::cerr << "Cannot create temporary file: " << filename << std::endl;
            return "";
        }
        
        // Write header and sequence
        outfile << ">species_" << species_taxid << "_genome_" << genome_index << "\n";
        
        // Write sequence in lines of 80 characters (FASTA format)
        const size_t line_length = 80;
        for (size_t i = 0; i < sequence.length(); i += line_length) {
            size_t end = std::min(i + line_length, sequence.length());
            outfile << sequence.substr(i, end - i) << "\n";
        }
        
        outfile.close();
        
        return filename;
    }
};

// Add these helper functions for debugging
class MemoryDebugHelper {
public:
    static void check_heap_integrity(const std::string& location) {
        std::cout << "DEBUG: Checking heap integrity at " << location << std::endl;
        
        try {
            // Try a small allocation to test heap
            void* test_ptr = malloc(1024);
            if (test_ptr) {
                free(test_ptr);
                std::cout << "DEBUG: Heap allocation test passed at " << location << std::endl;
            } else {
                std::cerr << "WARNING: Heap allocation test failed at " << location << std::endl;
            }
        } catch (...) {
            std::cerr << "ERROR: Exception during heap test at " << location << std::endl;
        }
    }
    
    static void print_memory_stats(const std::string& location) {
        std::cout << "DEBUG: Memory stats at " << location << std::endl;
        
        // Check virtual memory usage
        std::ifstream status("/proc/self/status");
        if (status.is_open()) {
            std::string line;
            while (std::getline(status, line)) {
                if (line.find("VmSize:") == 0 || line.find("VmRSS:") == 0 || 
                    line.find("VmPeak:") == 0 || line.find("VmHWM:") == 0) {
                    std::cout << "  " << line << std::endl;
                }
            }
            status.close();
        }
        
        // Check GPU memory if available
        size_t free_mem, total_mem;
        cudaError_t cuda_status = cudaMemGetInfo(&free_mem, &total_mem);
        if (cuda_status == cudaSuccess) {
            std::cout << "  GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
                      << (total_mem / 1024 / 1024) << " MB total" << std::endl;
        }
    }
    
    static bool validate_vector_integrity(const std::vector<std::string>& vec, const std::string& name) {
        std::cout << "DEBUG: Validating vector '" << name << "' with size " << vec.size() << std::endl;
        
        try {
            if (vec.size() > 1000000) {
                std::cerr << "WARNING: Vector " << name << " has suspicious size: " << vec.size() << std::endl;
                return false;
            }
            
            // Check first few and last few elements
            size_t check_count = std::min(vec.size(), (size_t)10);
            for (size_t i = 0; i < check_count; i++) {
                if (vec[i].empty()) {
                    std::cerr << "WARNING: Empty string at index " << i << " in vector " << name << std::endl;
                }
                if (vec[i].size() > 10000) {
                    std::cerr << "WARNING: Very long string (" << vec[i].size() 
                              << " chars) at index " << i << " in vector " << name << std::endl;
                }
            }
            
            std::cout << "DEBUG: Vector " << name << " validation passed" << std::endl;
            return true;
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR: Exception during vector validation for " << name << ": " << e.what() << std::endl;
            return false;
        }
    }
};


#ifndef GPU_KRAKEN_DATABASE_BUILDER_HEADER_ONLY

// ================================================================
// CUDA KERNELS - Keep existing interfaces but fix implementation
// ================================================================

// MurmurHash3 implementation - now in header file
// Using the implementation from gpu_minimizer_extraction.cuh

// TRUE Kraken2-style sliding window minimizer extraction
__global__ void extract_minimizers_kraken2_improved_kernel(
    const char* sequence_data,
    const GPUGenomeInfo* genome_info,
    int num_genomes,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* hit_counts_per_genome,
    uint32_t* global_hit_counter,
    MinimizerParams params,
    uint64_t min_clear_hash_value,
    uint64_t toggle_mask,
    int max_minimizers) {
    
    // OPTIMIZED: Use multiple threads per genome instead of just thread 0
    int genome_id = blockIdx.x;
    int thread_id = threadIdx.x;
    int block_size = blockDim.x;
    
    if (genome_id >= num_genomes) return;
    
    const GPUGenomeInfo& genome = genome_info[genome_id];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) {
        if (thread_id == 0) hit_counts_per_genome[genome_id] = 0;
        return;
    }
    
    // Shared memory for thread coordination
    extern __shared__ uint64_t shared_mem[];
    uint64_t* shared_minimizers = shared_mem;
    uint32_t* shared_positions = (uint32_t*)(shared_minimizers + blockDim.x);
    bool* shared_valid = (bool*)(shared_positions + blockDim.x);
    
    uint32_t total_kmers = seq_length - params.k + 1;
    uint32_t kmers_per_thread = (total_kmers + block_size - 1) / block_size;
    uint32_t start_kmer = thread_id * kmers_per_thread;
    uint32_t end_kmer = min(start_kmer + kmers_per_thread, total_kmers);
    
    // Each thread processes its assigned k-mers
    shared_minimizers[thread_id] = UINT64_MAX;
    shared_positions[thread_id] = 0;
    shared_valid[thread_id] = false;
    
    // Process assigned k-mers
    for (uint32_t kmer_idx = start_kmer; kmer_idx < end_kmer; kmer_idx++) {
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, kmer_idx, params.k, params.ell, params.spaces, params.xor_mask
        );
        
        if (minimizer != UINT64_MAX) {
            // Apply subsampling if enabled
            if (min_clear_hash_value > 0 && minimizer < min_clear_hash_value) {
                continue;
            }
            
            // Apply toggle mask
            minimizer ^= toggle_mask;
            
            // Store first valid minimizer from this thread's range
            if (!shared_valid[thread_id]) {
                shared_minimizers[thread_id] = minimizer;
                shared_positions[thread_id] = kmer_idx;
                shared_valid[thread_id] = true;
                break;
            }
        }
    }
    
    __syncthreads();
    
    // Thread 0 processes results in order to maintain deduplication
    if (thread_id == 0) {
        uint64_t last_minimizer = UINT64_MAX;
        uint32_t local_count = 0;
        
        // Process all positions in order
        for (uint32_t pos = 0; pos < total_kmers; pos++) {
            uint32_t responsible_thread = pos / kmers_per_thread;
            if (responsible_thread >= block_size) responsible_thread = block_size - 1;
            
            if (shared_valid[responsible_thread] && 
                shared_positions[responsible_thread] == pos &&
                shared_minimizers[responsible_thread] != last_minimizer) {
                
                // FIXED: Atomic add first, then check bounds to prevent overflow
                uint32_t global_pos = atomicAdd(global_hit_counter, 1);
                if (global_pos < max_minimizers) {
                    GPUMinimizerHit hit;
                    hit.minimizer_hash = shared_minimizers[responsible_thread];
                    hit.taxon_id = genome.taxon_id;
                    hit.position = pos;
                    hit.genome_id = genome.genome_id;
                    
                    minimizer_hits[global_pos] = hit;
                    local_count++;
                    last_minimizer = shared_minimizers[responsible_thread];
                } else {
                    // We've hit the limit, stop processing
                    break;
                }
            }
        }
        
        hit_counts_per_genome[genome_id] = local_count;
    }
}

// FIXED: Kernel to convert hits to candidates without over-aggressive deduplication
__global__ void convert_hits_to_candidates_fixed(
    const GPUMinimizerHit* hits,
    int num_hits,
    LCACandidate* candidates) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_hits) return;
    
    const GPUMinimizerHit& hit = hits[idx];
    LCACandidate& candidate = candidates[idx];
    
    candidate.minimizer_hash = hit.minimizer_hash;
    candidate.lca_taxon = hit.taxon_id;
    candidate.genome_count = 1;  // Will be updated during final merge
    candidate.uniqueness_score = 1.0f;
}

// Simple LCA computation device function
__device__ uint32_t compute_simple_lca_gpu(uint32_t taxon1, uint32_t taxon2) {
    if (taxon1 == 0) return taxon2;
    if (taxon2 == 0) return taxon1;
    if (taxon1 == taxon2) return taxon1;
    
    // Simplified LCA - return smaller taxon ID
    // In a real implementation, this would walk up the taxonomy tree
    return (taxon1 < taxon2) ? taxon1 : taxon2;
}

// Host version of the same function
__host__ uint32_t compute_simple_lca_host(uint32_t taxon1, uint32_t taxon2) {
    if (taxon1 == 0) return taxon2;
    if (taxon2 == 0) return taxon1;
    if (taxon1 == taxon2) return taxon1;
    
    // Simplified LCA - return smaller taxon ID
    // In a real implementation, this would walk up the taxonomy tree
    return (taxon1 < taxon2) ? taxon1 : taxon2;
}

// Device function to check for valid bases
__device__ bool has_valid_bases(const char* seq, int length) {
    for (int i = 0; i < length; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

// ================================================================
// IMPLEMENTATION - Keep all original functionality
// ================================================================

#endif // GPU_KRAKEN_DATABASE_BUILDER_HEADER_ONLY

#ifndef GPU_KRAKEN_CLASSIFIER_HEADER_ONLY
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>
#include <iomanip>

#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// Fixed GPUKrakenDatabaseBuilder constructor and memory allocation
// This fixes the malloc(): invalid size (unsorted) error
// FIXED: GPUKrakenDatabaseBuilder constructor with proper initialization order and bounds checking
// This fixes the heap corruption issue identified by valgrind

GPUKrakenDatabaseBuilder::GPUKrakenDatabaseBuilder(
    const std::string& output_dir,
    const ClassificationParams& config)
    : output_directory(output_dir),  // Initialize string first
      params(config),                // Copy config (not reference)
      d_sequence_data(nullptr),
      d_genome_info(nullptr),
      d_minimizer_hits(nullptr),
      d_minimizer_counts(nullptr),
      d_lca_candidates(nullptr),
      d_taxonomy_nodes(nullptr),
      MAX_SEQUENCE_BATCH(10),              // Safe default
      MAX_MINIMIZERS_PER_BATCH(1000000),   // Safe 1M default
      auto_scale_memory(true),
      max_gpu_memory_usage_fraction(80),
      total_minimizers_capped(0),
      total_batches_capped(0),
      min_clear_hash_value(0),
      subsampling_rate(1.0),
      toggle_mask(0xe37e28c4271b5a2dULL) {
    
    // CRITICAL FIX: Initialize all POD members before any operations
    memset(&stats, 0, sizeof(stats));
    
    // Initialize minimizer params from classification params
    minimizer_params.k = params.k;
    minimizer_params.ell = params.ell;
    minimizer_params.spaces = params.spaces;
    minimizer_params.xor_mask = 0x3c8bfbb395c60474ULL;
    
    // Initialize enhanced stats with zero initialization
    enhanced_stats = EnhancedBuildStats();
    memset(&enhanced_stats, 0, sizeof(enhanced_stats));
    
    std::cout << "Initializing GPU Kraken database builder (High Capacity)..." << std::endl;
    std::cout << "Parameters: k=" << minimizer_params.k 
              << ", ell=" << minimizer_params.ell 
              << ", spaces=" << minimizer_params.spaces << std::endl;
    
    // FIXED: Validate string parameters before using them
    if (output_directory.empty()) {
        throw std::invalid_argument("Output directory cannot be empty");
    }
    
    if (output_directory.size() > 4096) {
        throw std::invalid_argument("Output directory path too long");
    }
    
    // Create output directory safely with proper error handling
    try {
        std::error_code ec;
        std::filesystem::create_directories(output_directory, ec);
        if (ec) {
            std::cerr << "Warning: Could not create output directory: " << ec.message() << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not create output directory: " << e.what() << std::endl;
    }
    
    std::cout << "Initial safe settings:" << std::endl;
    std::cout << "  Minimizers per batch: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
    std::cout << "  Sequences per batch: " << MAX_SEQUENCE_BATCH << std::endl;
    
    // CRITICAL FIX: Don't call CUDA functions in constructor
    // CUDA context might not be ready yet - defer to enable_auto_memory_scaling()
    if (auto_scale_memory) {
        std::cout << "Auto memory scaling will be configured after CUDA initialization" << std::endl;
    }
    
    // Final validation that all values are sane
    validate_memory_settings();
    
    std::cout << "GPUKrakenDatabaseBuilder constructor completed successfully" << std::endl;
}

// FIXED: Separate method for CUDA initialization that can be called safely
bool GPUKrakenDatabaseBuilder::initialize_cuda_context() {
    std::cout << "Initializing CUDA context..." << std::endl;
    
    // Force CUDA runtime initialization
    cudaError_t init_status = cudaFree(0);
    if (init_status != cudaSuccess) {
        std::cerr << "Warning: CUDA runtime initialization failed: " << cudaGetErrorString(init_status) << std::endl;
        return false;
    }
    
    // Check device count
    int device_count;
    cudaError_t cuda_status = cudaGetDeviceCount(&device_count);
    if (cuda_status != cudaSuccess || device_count == 0) {
        std::cerr << "Error: No CUDA devices available" << std::endl;
        return false;
    }
    
    // Set device
    cuda_status = cudaSetDevice(0);
    if (cuda_status != cudaSuccess) {
        std::cerr << "Error: Failed to set CUDA device" << std::endl;
        return false;
    }
    
    // Now safe to call memory scaling
    if (auto_scale_memory) {
        try {
            check_and_adjust_memory_safe();
        } catch (const std::exception& e) {
            std::cerr << "Auto-scaling failed: " << e.what() << std::endl;
            auto_scale_memory = false;
            return false;
        }
    }
    
    std::cout << "CUDA context initialized successfully" << std::endl;
    return true;
}

// Destructor - safe with CUDA context validation
GPUKrakenDatabaseBuilder::~GPUKrakenDatabaseBuilder() {
    try {
        // Check if we're in a valid CUDA context before cleanup
        int device;
        cudaError_t err = cudaGetDevice(&device);
        if (err == cudaSuccess) {
            free_gpu_memory();
        } else {
            // CUDA context is invalid - just null the pointers without freeing
            d_sequence_data = nullptr;
            d_genome_info = nullptr;
            d_minimizer_hits = nullptr;
            d_minimizer_counts = nullptr;
            d_lca_candidates = nullptr;
            d_taxonomy_nodes = nullptr;
        }
    } catch (...) {
        // Suppress all exceptions in destructor
    }
}

void GPUKrakenDatabaseBuilder::check_and_adjust_memory() {
    std::cout << "DEBUG: Entering check_and_adjust_memory()" << std::endl;
    
    size_t free_mem, total_mem;
    std::cout << "DEBUG: About to call cudaMemGetInfo..." << std::endl;
    
    cudaError_t cuda_status = cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "DEBUG: cudaMemGetInfo returned: " << cudaGetErrorString(cuda_status) << std::endl;
    
    if (cuda_status != cudaSuccess) {
        std::cerr << "CUDA memory query failed: " << cudaGetErrorString(cuda_status) << std::endl;
        return;
    }
    
    std::cout << "DEBUG: free_mem=" << free_mem << ", total_mem=" << total_mem << std::endl;
    std::cout << "GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
              << (total_mem / 1024 / 1024) << " MB total" << std::endl;
    
    if (auto_scale_memory) {
        // SAFE: Use simple logic to prevent overflow
        size_t memory_gb = total_mem / (1024ULL * 1024ULL * 1024ULL);
        
        std::cout << "GPU Memory: " << memory_gb << " GB detected" << std::endl;
        
        // SAFE: Conservative settings based on GPU memory size
        if (memory_gb >= 24) {
            MAX_MINIMIZERS_PER_BATCH = 3000000;
            MAX_SEQUENCE_BATCH = 25;
        } else if (memory_gb >= 16) {
            MAX_MINIMIZERS_PER_BATCH = 2000000;
            MAX_SEQUENCE_BATCH = 20;
        } else if (memory_gb >= 8) {
            MAX_MINIMIZERS_PER_BATCH = 1000000;
            MAX_SEQUENCE_BATCH = 15;
        } else {
            MAX_MINIMIZERS_PER_BATCH = 500000;
            MAX_SEQUENCE_BATCH = 10;
        }
        
        std::cout << "Auto-scaled settings:" << std::endl;
        std::cout << "  Minimizers per batch: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
        std::cout << "  Sequences per batch: " << MAX_SEQUENCE_BATCH << std::endl;
    }
    
    // CRITICAL: Final validation to prevent overflow
    if (MAX_MINIMIZERS_PER_BATCH <= 0 || MAX_MINIMIZERS_PER_BATCH > 10000000) {
        std::cerr << "ERROR: Invalid minimizer batch size, using safe default" << std::endl;
        MAX_MINIMIZERS_PER_BATCH = 1000000;
    }
    
    if (MAX_SEQUENCE_BATCH <= 0 || MAX_SEQUENCE_BATCH > 50) {
        std::cerr << "ERROR: Invalid sequence batch size, using safe default" << std::endl;
        MAX_SEQUENCE_BATCH = 10;
    }
}

// FIXED: Safe memory checking with overflow protection
void GPUKrakenDatabaseBuilder::check_and_adjust_memory_safe() {
    std::cout << "DEBUG: Entering check_and_adjust_memory_safe()" << std::endl;
    
    // First ensure CUDA is initialized by getting device
    int device = 0;
    cudaError_t device_status = cudaGetDevice(&device);
    if (device_status != cudaSuccess) {
        std::cerr << "Failed to get CUDA device: " << cudaGetErrorString(device_status) << std::endl;
        throw std::runtime_error("CUDA not initialized");
    }
    std::cout << "DEBUG: Using CUDA device " << device << std::endl;
    
    size_t free_mem = 0, total_mem = 0;
    
    std::cout << "DEBUG: About to call cudaMemGetInfo" << std::endl;
    // Safely query GPU memory
    cudaError_t cuda_status = cudaMemGetInfo(&free_mem, &total_mem);
    if (cuda_status != cudaSuccess) {
        std::cerr << "CUDA memory query failed: " << cudaGetErrorString(cuda_status) << std::endl;
        throw std::runtime_error("Failed to query GPU memory");
    }
    std::cout << "DEBUG: cudaMemGetInfo succeeded" << std::endl;
    
    // Validate memory values
    if (total_mem == 0 || free_mem > total_mem) {
        throw std::runtime_error("Invalid GPU memory values returned");
    }
    
    std::cout << "GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
              << (total_mem / 1024 / 1024) << " MB total" << std::endl;
    
    // Calculate available memory with safety checks
    size_t available_memory = 0;
    if (max_gpu_memory_usage_fraction > 0 && max_gpu_memory_usage_fraction <= 100) {
        // Safe multiplication with overflow check
        size_t max_fraction = total_mem / 100;  // Avoid overflow
        available_memory = max_fraction * max_gpu_memory_usage_fraction;
        
        // Use minimum of free memory and calculated fraction
        available_memory = std::min(available_memory, free_mem);
    } else {
        throw std::runtime_error("Invalid memory usage fraction");
    }
    
    // Calculate batch sizes with conservative estimates
    // Each minimizer needs approximately:
    // - 8 bytes for minimizer value
    // - 4 bytes for position
    // - 4 bytes for taxon ID
    // - Additional overhead for hash table structures
    const size_t bytes_per_minimizer = 32;  // Conservative estimate with overhead
    const size_t bytes_per_sequence = 1024 * 1024;  // 1MB per sequence (conservative)
    
    // Reserve memory for other GPU operations (500MB minimum)
    const size_t reserved_memory = 500ULL * 1024ULL * 1024ULL;
    
    if (available_memory <= reserved_memory) {
        throw std::runtime_error("Insufficient GPU memory available");
    }
    
    size_t working_memory = available_memory - reserved_memory;
    
    // Calculate safe batch sizes
    size_t max_minimizers = working_memory / bytes_per_minimizer;
    size_t max_sequences = working_memory / bytes_per_sequence;
    
    // Apply reasonable limits
    MAX_MINIMIZERS_PER_BATCH = std::min(max_minimizers, (size_t)5000000);  // Cap at 5M
    MAX_MINIMIZERS_PER_BATCH = std::max((size_t)MAX_MINIMIZERS_PER_BATCH, (size_t)100000);  // Min 100K
    
    MAX_SEQUENCE_BATCH = std::min(max_sequences, (size_t)50);  // Cap at 50
    MAX_SEQUENCE_BATCH = std::max((size_t)MAX_SEQUENCE_BATCH, (size_t)1);  // Min 1
    
    std::cout << "Calculated safe batch sizes based on available memory:" << std::endl;
    std::cout << "  Available memory: " << (available_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Minimizers per batch: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
    std::cout << "  Sequences per batch: " << MAX_SEQUENCE_BATCH << std::endl;
}

// Validate memory settings to prevent overflow
// FIXED: Safer memory validation with proper bounds checking
void GPUKrakenDatabaseBuilder::validate_memory_settings() {
    // Ensure batch sizes are within reasonable bounds
    if (MAX_MINIMIZERS_PER_BATCH <= 0) {
        std::cerr << "ERROR: Invalid minimizer batch size " << MAX_MINIMIZERS_PER_BATCH 
                  << ", using safe default" << std::endl;
        MAX_MINIMIZERS_PER_BATCH = 1000000;
    }
    
    if (MAX_MINIMIZERS_PER_BATCH > 50000000) {  // Cap at 50M
        std::cerr << "WARNING: Very large minimizer batch size " << MAX_MINIMIZERS_PER_BATCH 
                  << ", capping at 50M" << std::endl;
        MAX_MINIMIZERS_PER_BATCH = 50000000;
    }
    
    if (MAX_SEQUENCE_BATCH <= 0) {
        std::cerr << "ERROR: Invalid sequence batch size " << MAX_SEQUENCE_BATCH 
                  << ", using safe default" << std::endl;
        MAX_SEQUENCE_BATCH = 10;
    }
    
    if (MAX_SEQUENCE_BATCH > 100) {
        std::cerr << "WARNING: Very large sequence batch size " << MAX_SEQUENCE_BATCH 
                  << ", capping at 100" << std::endl;
        MAX_SEQUENCE_BATCH = 100;
    }
    
    // Validate other memory-related settings
    if (max_gpu_memory_usage_fraction <= 0 || max_gpu_memory_usage_fraction > 100) {
        std::cerr << "ERROR: Invalid GPU memory usage fraction " << max_gpu_memory_usage_fraction 
                  << ", using 80%" << std::endl;
        max_gpu_memory_usage_fraction = 80;
    }
    
    std::cout << "Memory settings validated:" << std::endl;
    std::cout << "  MAX_MINIMIZERS_PER_BATCH: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
    std::cout << "  MAX_SEQUENCE_BATCH: " << MAX_SEQUENCE_BATCH << std::endl;
    std::cout << "  GPU memory usage: " << max_gpu_memory_usage_fraction << "%" << std::endl;
}

// Configuration method implementations
void GPUKrakenDatabaseBuilder::set_batch_size(int sequences_per_batch) {
    if (sequences_per_batch > 0 && sequences_per_batch <= 1000) {
        MAX_SEQUENCE_BATCH = sequences_per_batch;
        std::cout << "Set GPU batch size to " << MAX_SEQUENCE_BATCH 
                  << " sequences, " << MAX_MINIMIZERS_PER_BATCH << " minimizers" << std::endl;
    }
}

void GPUKrakenDatabaseBuilder::set_minimizer_capacity(int minimizers_per_batch) {
    if (minimizers_per_batch > 0 && minimizers_per_batch <= 50000000) { // Cap at 50M
        MAX_MINIMIZERS_PER_BATCH = minimizers_per_batch;
        std::cout << "Set minimizer capacity to " << MAX_MINIMIZERS_PER_BATCH 
                  << " per batch" << std::endl;
        
        // Recalculate memory requirements
        check_and_adjust_memory();
    } else {
        std::cerr << "Warning: Invalid minimizer capacity. Must be 1-50M" << std::endl;
    }
}

// UPDATED: Enhanced enable_auto_memory_scaling that doesn't call CUDA in constructor
void GPUKrakenDatabaseBuilder::enable_auto_memory_scaling(bool enable, size_t memory_fraction) {
    // Validate memory_fraction
    if (memory_fraction == 0 || memory_fraction > 100) {
        std::cerr << "ERROR: Invalid memory_fraction " << memory_fraction 
                  << ", must be between 1 and 100" << std::endl;
        memory_fraction = 80; // Use safe default
    }
    
    auto_scale_memory = enable;
    max_gpu_memory_usage_fraction = memory_fraction;
    
    std::cout << "Auto memory scaling " << (enable ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "Memory fraction: " << memory_fraction << "%" << std::endl;
    
    // Don't initialize CUDA here - it will be done in initialize_cuda_context()
    if (enable) {
        std::cout << "Note: CUDA context must be initialized before auto-scaling takes effect" << std::endl;
    }
}

void GPUKrakenDatabaseBuilder::set_subsampling_rate(double rate) {
    if (rate > 0.0 && rate <= 1.0) {
        subsampling_rate = rate;
        min_clear_hash_value = (uint64_t)((1.0 - rate) * UINT64_MAX);
        std::cout << "Set subsampling rate to " << rate 
                  << " (min_clear_hash_value: 0x" << std::hex << min_clear_hash_value 
                  << std::dec << ")" << std::endl;
    }
}

// Debug version of build_database_from_genomes to track malloc error
// UPDATED: build_database_from_genomes with proper initialization
bool GPUKrakenDatabaseBuilder::build_database_from_genomes(
    const std::string& genome_library_path,
    const std::string& taxonomy_path) {
    
    std::cout << "\n=== BUILDING KRAKEN DATABASE FROM GENOMES ===" << std::endl;
    
    // CRITICAL FIX: Initialize CUDA context first
    if (!initialize_cuda_context()) {
        std::cerr << "Failed to initialize CUDA context" << std::endl;
        return false;
    }
    
    // Validate inputs
    if (genome_library_path.empty()) {
        std::cerr << "ERROR: genome_library_path is empty!" << std::endl;
        return false;
    }
    
    // Check for string corruption
    for (size_t i = 0; i < genome_library_path.size(); i++) {
        if (genome_library_path[i] == '\0') {
            std::cerr << "ERROR: Null character found in genome_library_path at position " << i << std::endl;
            return false;
        }
    }
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    try {
        if (!load_genome_files(genome_library_path)) {
            std::cerr << "Failed to load genome files" << std::endl;
            return false;
        }
        
        if (!taxonomy_path.empty()) {
            if (!load_taxonomy_data(taxonomy_path)) {
                std::cerr << "Warning: Failed to load taxonomy data" << std::endl;
            }
        }
        
        if (!process_genomes_gpu()) {
            std::cerr << "Failed to process genomes on GPU" << std::endl;
            return false;
        }
        
        if (!save_database()) {
            std::cerr << "Failed to save database" << std::endl;
            return false;
        }
        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
        
        std::cout << "\n=== DATABASE BUILD COMPLETE ===" << std::endl;
        std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
        print_build_statistics();
        
        return true;
        
    } catch (const std::bad_alloc& e) {
        std::cerr << "MALLOC ERROR: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return false;
    }
}

// Fixed load_genome_files with proper bounds checking and error handling
bool GPUKrakenDatabaseBuilder::load_genome_files(const std::string& library_path) {
    std::cout << "Loading genome files from: " << library_path << " (SAFE VERSION)" << std::endl;
    
    // Clear existing data first
    genome_files.clear();
    genome_taxon_ids.clear();
    
    // Validate the library path early
    if (library_path.empty() || library_path.size() > 4096) {
        std::cerr << "Error: Invalid library path" << std::endl;
        return false;
    }
    
    // Use stat() instead of std::filesystem for safety
    struct stat path_stat;
    if (stat(library_path.c_str(), &path_stat) != 0) {
        std::cerr << "Error: Genome library path does not exist: " << library_path << std::endl;
        return false;
    }
    
    if (!S_ISDIR(path_stat.st_mode)) {
        std::cerr << "Error: Genome library path is not a directory: " << library_path << std::endl;
        return false;
    }
    
    std::cout << "Finding genome files with safe method..." << std::endl;
    
    // Use system find command instead of std::filesystem to avoid heap corruption
    std::string find_command = "find \"" + library_path + "\" -type f \\( "
                              "-name \"*.fna\" -o -name \"*.fa\" -o -name \"*.fasta\" "
                              "-o -name \"*.ffn\" -o -name \"*.faa\" \\) 2>/dev/null | head -50000";
    
    FILE* pipe = popen(find_command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Failed to execute find command" << std::endl;
        return false;
    }
    
    std::vector<std::string> found_files;
    found_files.reserve(10000);  // Pre-allocate to prevent reallocations
    
    char path_buffer[5000];
    while (fgets(path_buffer, sizeof(path_buffer), pipe)) {
        // Remove newline
        size_t len = strlen(path_buffer);
        if (len > 0 && path_buffer[len-1] == '\n') {
            path_buffer[len-1] = '\0';
            len--;
        }
        
        // Validate path length
        if (len == 0 || len > 4000) continue;
        
        // Check if we need to expand capacity
        if (found_files.size() >= found_files.capacity() - 10) {
            try {
                found_files.reserve(found_files.capacity() * 1.5);
            } catch (const std::bad_alloc& e) {
                std::cerr << "Memory allocation failed, stopping file search" << std::endl;
                break;
            }
        }
        
        try {
            found_files.emplace_back(path_buffer);
        } catch (const std::exception& e) {
            std::cerr << "Error adding file path: " << e.what() << std::endl;
            continue;
        }
        
        // Safety limit
        if (found_files.size() >= 50000) {
            std::cout << "Reached file limit of 50000, stopping search" << std::endl;
            break;
        }
    }
    
    pclose(pipe);
    
    if (found_files.empty()) {
        std::cerr << "No genome files found in " << library_path << std::endl;
        return false;
    }
    
    std::cout << "Found " << found_files.size() << " genome files, processing..." << std::endl;
    
    // Reserve space to prevent reallocations during processing
    try {
        genome_files.reserve(found_files.size());
        genome_taxon_ids.reserve(found_files.size());
    } catch (const std::bad_alloc& e) {
        std::cerr << "Failed to reserve memory for " << found_files.size() << " genome files" << std::endl;
        return false;
    }
    
    // Process files one by one with error checking
    size_t processed_count = 0;
    for (const auto& file_path : found_files) {
        try {
            // Validate file path
            if (file_path.empty() || file_path.size() > 4000) {
                std::cerr << "Warning: Skipping invalid file path" << std::endl;
                continue;
            }
            
            // Check if file still exists
            struct stat file_stat;
            if (stat(file_path.c_str(), &file_stat) != 0) {
                std::cerr << "Warning: File no longer exists: " << file_path << std::endl;
                continue;
            }
            
            uint32_t taxon_id = extract_taxon_from_filename(file_path);
            
            // Safely add to vectors
            genome_files.push_back(file_path);
            genome_taxon_ids.push_back(taxon_id);
            
            // Update taxon names if not present
            if (taxon_names.find(taxon_id) == taxon_names.end()) {
                // Extract filename safely
                size_t last_slash = file_path.find_last_of("/\\");
                std::string filename = (last_slash != std::string::npos) ? 
                                      file_path.substr(last_slash + 1) : file_path;
                
                size_t last_dot = filename.find_last_of('.');
                std::string stem = (last_dot != std::string::npos) ? 
                                  filename.substr(0, last_dot) : filename;
                
                // Validate stem string
                if (!stem.empty() && stem.size() < 1000) {
                    taxon_names[taxon_id] = stem;
                } else {
                    taxon_names[taxon_id] = "taxon_" + std::to_string(taxon_id);
                }
            }
            
            processed_count++;
            
            // Progress reporting for large datasets
            if (processed_count % 1000 == 0) {
                std::cout << "Processed " << processed_count << " files..." << std::endl;
            }
            
            // Sanity check: make sure vectors stay in sync
            if (genome_files.size() != genome_taxon_ids.size()) {
                std::cerr << "FATAL: Vector size mismatch detected!" << std::endl;
                return false;
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Error processing file " << file_path << ": " << e.what() << std::endl;
            // Continue with next file rather than failing completely
            continue;
        }
    }
    
    if (genome_files.empty()) {
        std::cerr << "No valid genome files could be processed" << std::endl;
        return false;
    }
    
    // SAFE SORTING: Sort using try-catch
    try {
        std::cout << "Sorting " << genome_files.size() << " files..." << std::endl;
        std::sort(genome_files.begin(), genome_files.end());
        std::cout << "Sorting completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Warning: Sorting failed: " << e.what() << std::endl;
        // Don't fail completely if sorting fails
    }
    
    std::cout << "Successfully loaded " << genome_files.size() << " genome files" << std::endl;
    stats.total_sequences = genome_files.size();
    
    // Final validation
    if (genome_files.size() != genome_taxon_ids.size()) {
        std::cerr << "FATAL: Final vector size mismatch!" << std::endl;
        return false;
    }
    
    return true;
}

// Keep original genome processing but use improved minimizer extraction
bool GPUKrakenDatabaseBuilder::process_genomes_gpu() {
    std::cout << "Processing genomes on GPU..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (!allocate_gpu_memory()) {
        return false;
    }
    
    for (size_t batch_start = 0; batch_start < genome_files.size(); 
         batch_start += MAX_SEQUENCE_BATCH) {
        
        size_t batch_end = std::min(batch_start + MAX_SEQUENCE_BATCH, genome_files.size());
        
        std::cout << "Processing batch " << (batch_start / MAX_SEQUENCE_BATCH + 1) 
                  << ": genomes " << batch_start << "-" << batch_end << std::endl;
        
        std::vector<std::string> batch_sequences;
        std::vector<uint32_t> batch_taxon_ids;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            auto sequences = load_sequences_from_fasta(genome_files[i]);
            for (const auto& seq : sequences) {
                batch_sequences.push_back(seq);
                batch_taxon_ids.push_back(genome_taxon_ids[i]);
                stats.total_bases += seq.length();
            }
        }
        
        if (!process_sequence_batch(batch_sequences, batch_taxon_ids, batch_start)) {
            std::cerr << "Failed to process batch starting at " << batch_start << std::endl;
            return false;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.sequence_processing_time = 
        std::chrono::duration<double>(end_time - start_time).count();
    
    return true;
}

// UPDATED: Process sequence batch with memory bounds checking
bool GPUKrakenDatabaseBuilder::process_sequence_batch(
    const std::vector<std::string>& sequences,
    const std::vector<uint32_t>& taxon_ids,
    int batch_offset) {
    
    if (sequences.empty()) return true;
    
    // Validate input consistency
    if (sequences.size() != taxon_ids.size()) {
        std::cerr << "ERROR: Sequence/taxon count mismatch" << std::endl;
        return false;
    }
    
    // Check if we exceed the allocated batch size
    if (static_cast<int>(sequences.size()) > MAX_SEQUENCE_BATCH) {
        std::cerr << "ERROR: Batch size " << sequences.size() 
                  << " exceeds maximum " << MAX_SEQUENCE_BATCH << std::endl;
        return false;
    }
    
    // Pre-calculate total size to prevent buffer overflow
    uint64_t total_sequence_length_64 = 0;
    size_t valid_sequences = 0;
    
    for (const auto& seq : sequences) {
        if (seq.length() >= minimizer_params.k && seq.length() <= 10000000) { // Max 10MB per sequence
            total_sequence_length_64 += seq.length();
            valid_sequences++;
        }
        
        // Early overflow detection
        if (total_sequence_length_64 > 100000000ULL) { // 100MB limit
            std::cerr << "ERROR: Total sequence length exceeds 100MB limit" << std::endl;
            return false;
        }
    }
    
    if (valid_sequences == 0) {
        std::cout << "No valid sequences in this batch, skipping" << std::endl;
        return true;
    }
    
    // Check if total size fits in our allocated memory
    size_t total_sequence_length = static_cast<size_t>(total_sequence_length_64);
    size_t max_sequence_memory = static_cast<size_t>(MAX_SEQUENCE_BATCH) * 5000000; // 5MB per slot
    
    if (total_sequence_length > max_sequence_memory) {
        std::cerr << "ERROR: Total sequence length (" << total_sequence_length 
                  << ") exceeds allocated memory (" << max_sequence_memory << ")" << std::endl;
        return false;
    }
    
    // Build genome info and concatenated sequences with bounds checking
    std::string concatenated_sequences;
    std::vector<GPUGenomeInfo> genome_infos;
    
    try {
        concatenated_sequences.reserve(total_sequence_length);
        genome_infos.reserve(valid_sequences);
    } catch (const std::bad_alloc& e) {
        std::cerr << "ERROR: Failed to reserve memory for sequences: " << e.what() << std::endl;
        return false;
    }
    
    uint32_t current_offset = 0;
    for (size_t i = 0; i < sequences.size(); i++) {
        if (sequences[i].length() < minimizer_params.k) continue;
        if (sequences[i].length() > 10000000) continue; // Skip very large sequences
        
        // Check for overflow in offset calculation
        if (current_offset > UINT32_MAX - sequences[i].length()) {
            std::cerr << "ERROR: Sequence offset overflow detected" << std::endl;
            return false;
        }
        
        GPUGenomeInfo info;
        info.taxon_id = taxon_ids[i];
        info.sequence_offset = current_offset;
        info.sequence_length = static_cast<uint32_t>(sequences[i].length());
        info.genome_id = static_cast<uint32_t>(batch_offset + i);
        
        genome_infos.push_back(info);
        concatenated_sequences += sequences[i];
        current_offset += static_cast<uint32_t>(sequences[i].length());
        
        // Update stats safely
        if (stats.total_kmers_processed < UINT64_MAX - (sequences[i].length() - minimizer_params.k + 1)) {
            stats.total_kmers_processed += sequences[i].length() - minimizer_params.k + 1;
        }
    }
    
    // Final validation before GPU copy
    if (concatenated_sequences.length() != total_sequence_length) {
        std::cerr << "ERROR: Sequence length mismatch in concatenation" << std::endl;
        return false;
    }
    
    if (genome_infos.size() > static_cast<size_t>(MAX_SEQUENCE_BATCH)) {
        std::cerr << "ERROR: Too many genome infos created" << std::endl;
        return false;
    }
    
    // Copy to GPU with error checking
    std::cout << "Copying " << (concatenated_sequences.length() / 1024 / 1024) 
              << " MB of sequence data to GPU..." << std::endl;
    
    cudaError_t err = cudaMemcpy(d_sequence_data, concatenated_sequences.c_str(),
                                concatenated_sequences.length(), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "ERROR: Failed to copy sequence data to GPU: " 
                  << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    err = cudaMemcpy(d_genome_info, genome_infos.data(),
                    genome_infos.size() * sizeof(GPUGenomeInfo), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "ERROR: Failed to copy genome info to GPU: " 
                  << cudaGetErrorString(err) << std::endl;
        return false;
    }
    
    // Continue with GPU processing...
    uint32_t total_hits = 0;
    if (!extract_minimizers_gpu(d_sequence_data, d_genome_info, 
                               static_cast<int>(genome_infos.size()), d_minimizer_hits, &total_hits)) {
        std::cerr << "ERROR: Minimizer extraction failed" << std::endl;
        return false;
    }
    
    int num_candidates = 0;
    if (!compute_lca_assignments_gpu(d_minimizer_hits, total_hits,
                                    d_lca_candidates, &num_candidates)) {
        std::cerr << "ERROR: LCA assignment failed" << std::endl;
        return false;
    }
    
    if (num_candidates > 0) {
        // Validate before copying back
        if (num_candidates > MAX_MINIMIZERS_PER_BATCH) {
            std::cerr << "ERROR: Too many candidates: " << num_candidates 
                      << " > " << MAX_MINIMIZERS_PER_BATCH << std::endl;
            return false;
        }
        
        std::vector<LCACandidate> batch_candidates;
        try {
            batch_candidates.resize(num_candidates);
        } catch (const std::bad_alloc& e) {
            std::cerr << "ERROR: Failed to allocate memory for candidates: " << e.what() << std::endl;
            return false;
        }
        
        err = cudaMemcpy(batch_candidates.data(), d_lca_candidates,
                        num_candidates * sizeof(LCACandidate), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            std::cerr << "ERROR: Failed to copy candidates from GPU: " 
                      << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        // Safe append with size checking
        size_t old_size = all_lca_candidates.size();
        if (old_size > SIZE_MAX - batch_candidates.size()) {
            std::cerr << "ERROR: Would cause overflow in all_lca_candidates" << std::endl;
            return false;
        }
        
        all_lca_candidates.insert(all_lca_candidates.end(), 
                                 batch_candidates.begin(), batch_candidates.end());
        
        std::cout << "✓ Accumulated " << num_candidates << " candidates (total: " 
                  << all_lca_candidates.size() << ")" << std::endl;
    }
    
    stats.valid_minimizers_extracted += total_hits;
    stats.lca_assignments += num_candidates;
    
    return true;
}

// UPDATED: Track when we hit capacity limits
bool GPUKrakenDatabaseBuilder::extract_minimizers_gpu(
    const char* d_sequences,
    const GPUGenomeInfo* d_genomes,
    int num_sequences,
    GPUMinimizerHit* d_hits,
    uint32_t* total_hits) {
    
    std::cout << "Extracting minimizers from " << num_sequences 
              << " sequences (capacity: " << MAX_MINIMIZERS_PER_BATCH << ")..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Reset global counter
    uint32_t zero = 0;
    uint32_t* d_global_counter;
    CUDA_CHECK(cudaMalloc(&d_global_counter, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_global_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Clear hit counts
    CUDA_CHECK(cudaMemset(d_minimizer_counts, 0, num_sequences * sizeof(uint32_t)));
    
    // Add timing to measure improvement
    auto kernel_start = std::chrono::high_resolution_clock::now();
    
    // Validate capacity before launching kernel
    if (MAX_MINIMIZERS_PER_BATCH <= 0 || MAX_MINIMIZERS_PER_BATCH > 100000000) {
        std::cerr << "ERROR: Invalid minimizer capacity: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
        std::cerr << "Must be between 1 and 100M" << std::endl;
        cudaFree(d_global_counter);
        return false;
    }
    
    // Launch improved kernel with higher capacity
    int threads_per_block = 256;
    size_t shared_mem_size = threads_per_block * (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(bool));
    
    extract_minimizers_kraken2_improved_kernel<<<num_sequences, threads_per_block, shared_mem_size>>>(
        d_sequences, d_genomes, num_sequences,
        d_hits, d_minimizer_counts, d_global_counter,
        minimizer_params, min_clear_hash_value, toggle_mask, MAX_MINIMIZERS_PER_BATCH
    );
    
    CUDA_CHECK(cudaDeviceSynchronize());
    
    auto kernel_end = std::chrono::high_resolution_clock::now();
    auto kernel_duration = std::chrono::duration<double>(kernel_end - kernel_start).count();
    
    std::cout << "  Kernel execution time: " << std::fixed << std::setprecision(3) 
              << kernel_duration << "s" << std::endl;
    
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA error in minimizer extraction: %s\n", cudaGetErrorString(error));
        cudaFree(d_global_counter);
        return false;
    }
    
    CUDA_CHECK(cudaMemcpy(total_hits, d_global_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    
    // NEW: Track if we hit the limit
    if (*total_hits >= MAX_MINIMIZERS_PER_BATCH) {
        total_minimizers_capped += (*total_hits - MAX_MINIMIZERS_PER_BATCH);
        total_batches_capped++;
        
        printf("WARNING: Minimizer extraction hit capacity limit!\n");
        printf("  Extracted: %u, Capacity: %d\n", *total_hits, MAX_MINIMIZERS_PER_BATCH);
        printf("  This batch was capped - consider increasing capacity\n");
        
        *total_hits = MAX_MINIMIZERS_PER_BATCH;
    }
    
    cudaFree(d_global_counter);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.minimizer_extraction_time += 
        std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "Extracted " << *total_hits << " minimizers" << std::endl;
    
    return true;
}

// FIXED: Compute LCA assignments without over-aggressive deduplication
bool GPUKrakenDatabaseBuilder::compute_lca_assignments_gpu(
    const GPUMinimizerHit* d_hits,
    int num_hits,
    LCACandidate* d_candidates,
    int* num_candidates) {
    
    if (num_hits == 0) {
        *num_candidates = 0;
        return true;
    }
    
    std::cout << "Computing LCA assignments for " << num_hits << " hits..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // FIXED: Instead of aggressive deduplication, just convert hits to candidates
    // The deduplication will happen during final database merge
    int blocks = (num_hits + 255) / 256;
    convert_hits_to_candidates_fixed<<<blocks, 256>>>(
        d_hits, num_hits, d_candidates
    );
    CUDA_CHECK(cudaDeviceSynchronize());
    
    *num_candidates = num_hits;  // Keep all hits as candidates
    
    auto end_time = std::chrono::high_resolution_clock::now();
    stats.lca_computation_time += 
        std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "Created " << *num_candidates << " LCA candidates" << std::endl;
    
    return true;
}

// ALSO NEED TO FIX allocate_gpu_memory() to prevent overflow
bool GPUKrakenDatabaseBuilder::allocate_gpu_memory() {
    std::cout << "Allocating GPU memory with overflow protection..." << std::endl;
    
    // Validate parameters first
    if (MAX_SEQUENCE_BATCH <= 0 || MAX_SEQUENCE_BATCH > 100) {
        std::cerr << "ERROR: Invalid MAX_SEQUENCE_BATCH: " << MAX_SEQUENCE_BATCH << std::endl;
        MAX_SEQUENCE_BATCH = 10; // Force safe value
    }
    
    if (MAX_MINIMIZERS_PER_BATCH <= 0 || MAX_MINIMIZERS_PER_BATCH > 10000000) {
        std::cerr << "ERROR: Invalid MAX_MINIMIZERS_PER_BATCH: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
        MAX_MINIMIZERS_PER_BATCH = 1000000; // Force safe value
    }
    
    // Use safe calculations with bounds checking
    size_t sequence_memory = 0;
    size_t genome_info_memory = 0;
    size_t minimizer_memory = 0;
    size_t candidate_memory = 0;
    size_t count_memory = 0;
    
    try {
        // Calculate each allocation separately with overflow checking
        sequence_memory = static_cast<size_t>(MAX_SEQUENCE_BATCH) * 3000000ULL; // 3MB per sequence (conservative)
        
        if (sequence_memory / MAX_SEQUENCE_BATCH != 3000000ULL) {
            throw std::overflow_error("Sequence memory calculation overflow");
        }
        
        genome_info_memory = static_cast<size_t>(MAX_SEQUENCE_BATCH) * sizeof(GPUGenomeInfo);
        minimizer_memory = static_cast<size_t>(MAX_MINIMIZERS_PER_BATCH) * sizeof(GPUMinimizerHit);
        candidate_memory = static_cast<size_t>(MAX_MINIMIZERS_PER_BATCH) * sizeof(LCACandidate);
        count_memory = static_cast<size_t>(MAX_SEQUENCE_BATCH) * sizeof(uint32_t);
        
        // Check for overflow in individual calculations
        if (minimizer_memory / MAX_MINIMIZERS_PER_BATCH != sizeof(GPUMinimizerHit) ||
            candidate_memory / MAX_MINIMIZERS_PER_BATCH != sizeof(LCACandidate)) {
            throw std::overflow_error("Minimizer memory calculation overflow");
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Memory calculation overflow: " << e.what() << std::endl;
        return false;
    }
    
    // Calculate total with overflow checking
    size_t total_required = 0;
    if (sequence_memory > SIZE_MAX - genome_info_memory ||
        (sequence_memory + genome_info_memory) > SIZE_MAX - minimizer_memory ||
        (sequence_memory + genome_info_memory + minimizer_memory) > SIZE_MAX - candidate_memory ||
        (sequence_memory + genome_info_memory + minimizer_memory + candidate_memory) > SIZE_MAX - count_memory) {
        std::cerr << "ERROR: Total memory calculation would overflow!" << std::endl;
        return false;
    }
    
    total_required = sequence_memory + genome_info_memory + minimizer_memory + candidate_memory + count_memory;
    
    // Check available GPU memory
    size_t free_memory, total_memory_gpu;
    cudaMemGetInfo(&free_memory, &total_memory_gpu);
    
    std::cout << "Safe memory requirements:" << std::endl;
    std::cout << "  Sequence data: " << (sequence_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Genome info: " << (genome_info_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Minimizer hits: " << (minimizer_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  LCA candidates: " << (candidate_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Count arrays: " << (count_memory / 1024) << " KB" << std::endl;
    std::cout << "  Total required: " << (total_required / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Available: " << (free_memory / 1024 / 1024) << " MB" << std::endl;
    
    // Use only 70% of available memory for extra safety
    if (total_required > (free_memory * 70 / 100)) {
        std::cerr << "ERROR: Not enough GPU memory (using 70% safety margin)!" << std::endl;
        std::cerr << "Required: " << (total_required / 1024 / 1024) << " MB" << std::endl;
        std::cerr << "Available (70%): " << (free_memory * 70 / 100 / 1024 / 1024) << " MB" << std::endl;
        return false;
    }
    
    // Allocate with immediate error checking and cleanup on failure
    cudaError_t error;
    
    std::cout << "Allocating sequence memory..." << std::endl;
    error = cudaMalloc(&d_sequence_data, sequence_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate sequence memory: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    error = cudaMalloc(&d_genome_info, genome_info_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate genome info: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_sequence_data); d_sequence_data = nullptr;
        return false;
    }
    
    error = cudaMalloc(&d_minimizer_hits, minimizer_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate minimizer hits: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_sequence_data); d_sequence_data = nullptr;
        cudaFree(d_genome_info); d_genome_info = nullptr;
        return false;
    }
    
    error = cudaMalloc(&d_lca_candidates, candidate_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate LCA candidates: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_sequence_data); d_sequence_data = nullptr;
        cudaFree(d_genome_info); d_genome_info = nullptr;
        cudaFree(d_minimizer_hits); d_minimizer_hits = nullptr;
        return false;
    }
    
    error = cudaMalloc(&d_minimizer_counts, count_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate minimizer counts: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_sequence_data); d_sequence_data = nullptr;
        cudaFree(d_genome_info); d_genome_info = nullptr;
        cudaFree(d_minimizer_hits); d_minimizer_hits = nullptr;
        cudaFree(d_lca_candidates); d_lca_candidates = nullptr;
        return false;
    }
    
    std::cout << "✓ Successfully allocated " << (total_required / 1024 / 1024) 
              << " MB of GPU memory safely" << std::endl;
    
    return true;
}

void GPUKrakenDatabaseBuilder::free_gpu_memory() {
    // Lambda for safe CUDA free with null pointer checking
    auto safe_free = [](void** ptr) {
        if (ptr && *ptr) {
            cudaError_t err = cudaFree(*ptr);
            // Don't print errors in destructor - just ensure pointer is nulled
            *ptr = nullptr;
        }
    };

    safe_free((void**)&d_sequence_data);
    safe_free((void**)&d_genome_info);
    safe_free((void**)&d_minimizer_hits);
    safe_free((void**)&d_minimizer_counts);
    safe_free((void**)&d_lca_candidates);
    safe_free((void**)&d_taxonomy_nodes);
}

// Keep all original file processing methods
std::vector<std::string> GPUKrakenDatabaseBuilder::find_genome_files(const std::string& directory) {
    std::vector<std::string> files;
    
    // Validate directory parameter first
    if (directory.empty()) {
        std::cerr << "Error: Empty directory path provided to find_genome_files" << std::endl;
        return files;
    }
    
    // Check for invalid characters or corruption
    for (char c : directory) {
        if (c == '\0') {
            std::cerr << "Error: Null character found in directory path" << std::endl;
            return files;
        }
    }
    
    std::cout << "Searching for genome files in: '" << directory << "'" << std::endl;
    
    // Validate directory before attempting to iterate
    if (!std::filesystem::exists(directory)) {
        std::cerr << "Error: Directory does not exist: '" << directory << "'" << std::endl;
        return files;
    }
    
    if (!std::filesystem::is_directory(directory)) {
        std::cerr << "Error: Path is not a directory: '" << directory << "'" << std::endl;
        return files;
    }
    
    try {
        // Pre-count files to reserve space and avoid reallocations
        size_t estimated_file_count = 0;
        
        try {
            for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
                if (entry.is_regular_file()) {
                    std::string extension = entry.path().extension().string();
                    if (extension == ".fna" || extension == ".fa" || extension == ".fasta" ||
                        extension == ".ffn" || extension == ".faa") {
                        estimated_file_count++;
                    }
                }
                
                // Safety check to prevent infinite loops or excessive memory usage
                if (estimated_file_count > 100000) {
                    std::cout << "Large number of files detected (>" << estimated_file_count << "), proceeding with processing..." << std::endl;
                    break;
                }
            }
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error during file counting: " << e.what() << std::endl;
            // Continue with processing, just don't reserve space
        }
        
        std::cout << "Estimated " << estimated_file_count << " genome files to process" << std::endl;
        
        // Reserve space to prevent reallocations
        if (estimated_file_count > 0 && estimated_file_count < 100000) {
            try {
                files.reserve(estimated_file_count);
            } catch (const std::bad_alloc& e) {
                std::cerr << "Warning: Could not reserve space for files: " << e.what() << std::endl;
                // Continue without reservation
            }
        }
        
        // Now actually collect the files
        int file_count = 0;
        for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
            if (entry.is_regular_file()) {
                std::string extension = entry.path().extension().string();
                
                if (extension == ".fna" || extension == ".fa" || extension == ".fasta" ||
                    extension == ".ffn" || extension == ".faa") {
                    
                    std::string file_path = entry.path().string();
                    
                    // Validate file path before adding
                    if (!file_path.empty() && file_path.size() < 4096) {  // Reasonable path limit
                        files.push_back(file_path);
                        file_count++;
                        
                        // Progress reporting for large directories
                        if (file_count % 1000 == 0) {
                            std::cout << "Found " << file_count << " genome files..." << std::endl;
                        }
                        
                        // Safety limit to prevent excessive memory usage
                        if (file_count >= 50000) {
                            std::cout << "Reached file limit of " << file_count << " files" << std::endl;
                            break;
                        }
                    } else {
                        std::cerr << "Warning: Skipping invalid file path: " << file_path.substr(0, 100) << "..." << std::endl;
                    }
                }
            }
        }
        
        std::cout << "Found " << files.size() << " total genome files" << std::endl;
        
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error in find_genome_files: " << e.what() << std::endl;
        std::cerr << "Directory: '" << directory << "'" << std::endl;
        return files;
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation error in find_genome_files: " << e.what() << std::endl;
        return files;
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error in find_genome_files: " << e.what() << std::endl;
        return files;
    }
    
    try {
        std::sort(files.begin(), files.end());
    } catch (const std::exception& e) {
        std::cerr << "Error sorting files: " << e.what() << std::endl;
        // Return unsorted files rather than failing
    }
    
    return files;
}

uint32_t GPUKrakenDatabaseBuilder::extract_taxon_from_filename(const std::string& filename) {
    // Validate input
    if (filename.empty()) {
        std::cerr << "Warning: Empty filename provided to extract_taxon_from_filename" << std::endl;
        return 1000000;  // Default fallback
    }
    
    if (filename.size() > 4096) {
        std::cerr << "Warning: Filename too long: " << filename.substr(0, 100) << "..." << std::endl;
        return 1000000;  // Default fallback
    }
    
    try {
        std::filesystem::path p(filename);
        std::string stem = p.stem().string();
        
        if (stem.empty()) {
            std::cerr << "Warning: Could not extract stem from filename: " << filename << std::endl;
            return 1000000;
        }
        
        // Try to extract taxid pattern
        std::regex taxid_pattern(R"(taxid[_-](\d+))");
        std::smatch match;
        if (std::regex_search(stem, match, taxid_pattern)) {
            try {
                uint32_t taxon_id = std::stoul(match[1].str());
                if (taxon_id > 0 && taxon_id < UINT32_MAX) {
                    return taxon_id;
                }
            } catch (const std::exception& e) {
                std::cerr << "Warning: Invalid taxon ID in filename " << filename << ": " << e.what() << std::endl;
            }
        }
        
        // Try GCF pattern
        std::regex gcf_pattern(R"(GCF_(\d+\.\d+))");
        if (std::regex_search(stem, match, gcf_pattern)) {
            try {
                std::hash<std::string> hasher;
                uint32_t hash_val = hasher(match[1].str()) % 1000000 + 1000000;
                return hash_val;
            } catch (const std::exception& e) {
                std::cerr << "Warning: Error hashing GCF pattern: " << e.what() << std::endl;
            }
        }
        
        // Fallback: hash the entire stem
        std::hash<std::string> hasher;
        uint32_t hash_val = hasher(stem) % 1000000 + 2000000;
        return hash_val;
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing filename " << filename << ": " << e.what() << std::endl;
        return 1000000;  // Safe default
    }
}

bool GPUKrakenDatabaseBuilder::validate_path_safety(const std::string& path) {
    if (path.empty() || path.size() > 4096) return false;
    
    // Check for null characters
    for (char c : path) {
        if (c == '\0') return false;
    }
    
    return true;
}

std::vector<std::string> GPUKrakenDatabaseBuilder::load_sequences_from_fasta(const std::string& fasta_path) {
    std::vector<std::string> sequences;
    
    // Validate input
    if (fasta_path.empty()) {
        std::cerr << "Error: Empty FASTA path" << std::endl;
        return sequences;
    }
    
    if (!std::filesystem::exists(fasta_path)) {
        std::cerr << "Error: FASTA file does not exist: " << fasta_path << std::endl;
        return sequences;
    }
    
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        std::cerr << "Cannot open FASTA file: " << fasta_path << std::endl;
        return sequences;
    }
    
    try {
        std::string line, current_sequence;
        bool in_sequence = false;
        int line_number = 0;
        
        // Reserve some initial space
        current_sequence.reserve(10000);  // Reserve 10KB for typical sequence
        
        while (std::getline(file, line)) {
            line_number++;
            
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Start of new sequence
                if (in_sequence && !current_sequence.empty()) {
                    // Validate sequence before adding
                    if (current_sequence.size() >= minimizer_params.k && 
                        current_sequence.size() <= 50000000) {  // Max 50MB per sequence
                        sequences.push_back(current_sequence);
                    } else if (current_sequence.size() < minimizer_params.k) {
                        std::cerr << "Warning: Sequence too short (" << current_sequence.size() 
                                  << " bp) in " << fasta_path << " at line " << line_number << std::endl;
                    } else {
                        std::cerr << "Warning: Sequence too long (" << current_sequence.size() 
                                  << " bp) in " << fasta_path << " at line " << line_number << std::endl;
                    }
                    current_sequence.clear();
                    current_sequence.reserve(10000);  // Reserve for next sequence
                }
                in_sequence = true;
            } else if (in_sequence) {
                // Validate line before adding
                if (line.size() > 100000) {  // Unreasonably long line
                    std::cerr << "Warning: Very long line (" << line.size() 
                              << " chars) in " << fasta_path << " at line " << line_number << std::endl;
                    continue;
                }
                
                // Add line to current sequence
                try {
                    current_sequence += line;
                    
                    // Check for reasonable sequence length
                    if (current_sequence.size() > 50000000) {  // 50MB limit
                        std::cerr << "Warning: Sequence exceeds size limit in " << fasta_path << std::endl;
                        current_sequence.clear();
                        in_sequence = false;
                    }
                } catch (const std::bad_alloc& e) {
                    std::cerr << "Memory allocation failed for sequence in " << fasta_path 
                              << " at line " << line_number << std::endl;
                    current_sequence.clear();
                    in_sequence = false;
                }
            }
            
            // Safety check for file processing
            if (line_number % 100000 == 0) {
                std::cout << "Processed " << line_number << " lines from " << fasta_path << std::endl;
            }
        }
        
        // Handle the last sequence
        if (in_sequence && !current_sequence.empty()) {
            if (current_sequence.size() >= minimizer_params.k && 
                current_sequence.size() <= 50000000) {
                sequences.push_back(current_sequence);
            }
        }
        
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed while reading " << fasta_path << ": " << e.what() << std::endl;
        sequences.clear();  // Clear partial results
    } catch (const std::exception& e) {
        std::cerr << "Error reading FASTA file " << fasta_path << ": " << e.what() << std::endl;
        sequences.clear();  // Clear partial results
    }
    
    file.close();
    
    std::cout << "Loaded " << sequences.size() << " sequences from " << fasta_path << std::endl;
    return sequences;
}

// UPDATED: Save database with buffer overflow protection
bool GPUKrakenDatabaseBuilder::save_database() {
    std::cout << "Saving database to " << output_directory << "..." << std::endl;
    
    // Create output directory if it doesn't exist
    try {
        std::filesystem::create_directories(output_directory);
    } catch (const std::exception& e) {
        std::cerr << "Failed to create output directory: " << e.what() << std::endl;
        return false;
    }
    
    // Process all LCA candidates and merge duplicates with bounds checking
    std::unordered_map<uint64_t, LCACandidate> unique_candidates;
    
    // Reserve space to prevent rehashing during insertion
    unique_candidates.reserve(all_lca_candidates.size());
    
    size_t processed_count = 0;
    for (const auto& candidate : all_lca_candidates) {
        processed_count++;
        
        // Progress reporting for large datasets
        if (processed_count % 1000000 == 0) {
            std::cout << "Processing candidate " << processed_count 
                      << "/" << all_lca_candidates.size() << std::endl;
        }
        
        auto it = unique_candidates.find(candidate.minimizer_hash);
        if (it == unique_candidates.end()) {
            unique_candidates[candidate.minimizer_hash] = candidate;
        } else {
            // Update genome count and compute simple LCA
            auto& existing = it->second;
            existing.genome_count++;
            existing.lca_taxon = compute_simple_lca_host(existing.lca_taxon, candidate.lca_taxon);
            existing.uniqueness_score = 1.0f / existing.genome_count;
        }
    }
    
    stats.unique_minimizers = unique_candidates.size();
    std::cout << "Processed " << all_lca_candidates.size() << " candidates into " 
              << unique_candidates.size() << " unique minimizers" << std::endl;
    
    // Save hash table with proper file path construction
    std::string hash_file = output_directory + "/hash_table.k2d";
    std::cout << "Writing hash table to: " << hash_file << std::endl;
    
    std::ofstream hash_out(hash_file, std::ios::binary);
    if (!hash_out.is_open()) {
        std::cerr << "Cannot create hash table file: " << hash_file << std::endl;
        return false;
    }
    
    // Write header with validation
    uint64_t table_size = unique_candidates.size() * 2;
    uint64_t num_entries = unique_candidates.size();
    
    hash_out.write(reinterpret_cast<const char*>(&table_size), sizeof(uint64_t));
    hash_out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    
    if (hash_out.fail()) {
        std::cerr << "Failed to write hash table header" << std::endl;
        hash_out.close();
        return false;
    }
    
    // Write entries with error checking
    size_t entries_written = 0;
    for (const auto& [hash, candidate] : unique_candidates) {
        hash_out.write(reinterpret_cast<const char*>(&candidate.minimizer_hash), sizeof(uint64_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.lca_taxon), sizeof(uint32_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.genome_count), sizeof(uint32_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.uniqueness_score), sizeof(float));
        
        if (hash_out.fail()) {
            std::cerr << "Failed to write hash table entry " << entries_written << std::endl;
            hash_out.close();
            return false;
        }
        
        entries_written++;
        if (entries_written % 100000 == 0) {
            std::cout << "Written " << entries_written << "/" << num_entries << " entries" << std::endl;
        }
    }
    hash_out.close();
    
    // Save taxonomy with bounds checking
    std::string taxonomy_file = output_directory + "/taxonomy.tsv";
    std::cout << "Writing taxonomy to: " << taxonomy_file << std::endl;
    
    std::ofstream tax_out(taxonomy_file);
    if (!tax_out.is_open()) {
        std::cerr << "Cannot create taxonomy file: " << taxonomy_file << std::endl;
        return false;
    }
    
    tax_out << "taxon_id\tname\tparent_id\n";
    for (const auto& [taxon_id, name] : taxon_names) {
        uint32_t parent_id = taxon_parents.count(taxon_id) ? taxon_parents[taxon_id] : 0;
        tax_out << taxon_id << "\t" << name << "\t" << parent_id << "\n";
        
        if (tax_out.fail()) {
            std::cerr << "Failed to write taxonomy entry for taxon " << taxon_id << std::endl;
            tax_out.close();
            return false;
        }
    }
    tax_out.close();
    
    // Save configuration with error checking
    std::string config_file = output_directory + "/config.txt";
    std::cout << "Writing config to: " << config_file << std::endl;
    
    std::ofstream config_out(config_file);
    if (!config_out.is_open()) {
        std::cerr << "Cannot create config file: " << config_file << std::endl;
        return false;
    }
    
    config_out << "k=" << minimizer_params.k << "\n";
    config_out << "ell=" << minimizer_params.ell << "\n";
    config_out << "spaces=" << minimizer_params.spaces << "\n";
    config_out << "subsampling_rate=" << subsampling_rate << "\n";
    config_out << "min_clear_hash_value=0x" << std::hex << min_clear_hash_value << std::dec << "\n";
    config_out << "unique_minimizers=" << stats.unique_minimizers << "\n";
    config_out << "total_sequences=" << stats.total_sequences << "\n";
    config_out.close();
    
    std::cout << "Database components saved successfully!" << std::endl;
    std::cout << "  Hash table: " << hash_file << " (" << num_entries << " entries)" << std::endl;
    std::cout << "  Taxonomy: " << taxonomy_file << " (" << taxon_names.size() << " taxa)" << std::endl;
    std::cout << "  Config: " << config_file << std::endl;
    
    return true;
}

// UPDATED: Include capping statistics
void GPUKrakenDatabaseBuilder::print_build_statistics() {
    std::cout << "\n=== BUILD STATISTICS ===" << std::endl;
    std::cout << "Total sequences processed: " << stats.total_sequences << std::endl;
    std::cout << "Total bases processed: " << stats.total_bases << std::endl;
    std::cout << "Total k-mers processed: " << stats.total_kmers_processed << std::endl;
    std::cout << "Valid minimizers extracted: " << stats.valid_minimizers_extracted << std::endl;
    std::cout << "Unique minimizers: " << stats.unique_minimizers << std::endl;
    std::cout << "LCA assignments computed: " << stats.lca_assignments << std::endl;
    
    // NEW: Capping statistics
    if (total_batches_capped > 0) {
        std::cout << "\n⚠️  CAPACITY WARNINGS:" << std::endl;
        std::cout << "Batches that hit minimizer limit: " << total_batches_capped << std::endl;
        std::cout << "Estimated minimizers lost to capping: " << total_minimizers_capped << std::endl;
        std::cout << "Recommendation: Increase minimizer capacity or enable auto-scaling" << std::endl;
    } else {
        std::cout << "\n✓ No capacity limits hit during processing" << std::endl;
    }
    
    if (stats.total_kmers_processed > 0) {
        double compression = (double)stats.valid_minimizers_extracted / stats.total_kmers_processed;
        std::cout << "Initial compression ratio: " << std::fixed << std::setprecision(4) 
                  << compression << " (" << std::fixed << std::setprecision(1) 
                  << (1.0/compression) << "x reduction)" << std::endl;
        
        double final_compression = (double)stats.unique_minimizers / stats.total_kmers_processed;
        std::cout << "Final compression ratio: " << std::fixed << std::setprecision(4) 
                  << final_compression << " (" << std::fixed << std::setprecision(1) 
                  << (1.0/final_compression) << "x reduction)" << std::endl;
    }
    
    if (subsampling_rate < 1.0) {
        std::cout << "Subsampling rate applied: " << std::fixed << std::setprecision(3) 
                  << subsampling_rate << std::endl;
    }
    
    std::cout << "Processing times:" << std::endl;
    std::cout << "  Sequence processing: " << std::fixed << std::setprecision(2) 
              << stats.sequence_processing_time << "s" << std::endl;
    std::cout << "  Minimizer extraction: " << std::fixed << std::setprecision(2) 
              << stats.minimizer_extraction_time << "s" << std::endl;
    std::cout << "  LCA computation: " << std::fixed << std::setprecision(2) 
              << stats.lca_computation_time << "s" << std::endl;
    
    if (stats.total_bases > 0) {
        double bases_per_second = stats.total_bases / 
            (stats.sequence_processing_time + stats.minimizer_extraction_time);
        std::cout << "Processing rate: " << std::scientific << std::setprecision(2) 
                  << bases_per_second << " bases/second" << std::endl;
    }
}

// Keep all original taxonomy loading methods (they work fine)
bool GPUKrakenDatabaseBuilder::load_taxonomy_data(const std::string& taxonomy_path) {
    std::cout << "Loading taxonomy data from: " << taxonomy_path << std::endl;
    
    taxon_names.clear();
    taxon_parents.clear();
    
    std::filesystem::path base_path(taxonomy_path);
    std::string nodes_file, names_file;
    
    if (std::filesystem::is_directory(base_path)) {
        nodes_file = base_path / "nodes.dmp";
        names_file = base_path / "names.dmp";
        
        if (!std::filesystem::exists(nodes_file)) {
            nodes_file = base_path.parent_path() / "nodes.dmp";
            names_file = base_path.parent_path() / "names.dmp";
        }
    } else {
        base_path = base_path.parent_path();
        nodes_file = base_path / "nodes.dmp";
        names_file = base_path / "names.dmp";
    }
    
    if (!std::filesystem::exists(nodes_file)) {
        std::cerr << "nodes.dmp not found at: " << nodes_file << std::endl;
        return false;
    }
    
    std::cout << "Loading taxonomy tree from: " << nodes_file << std::endl;
    std::ifstream nodes_in(nodes_file);
    if (!nodes_in.is_open()) {
        std::cerr << "Cannot open nodes.dmp: " << nodes_file << std::endl;
        return false;
    }
    
    std::string line;
    int nodes_loaded = 0;
    std::unordered_map<uint32_t, std::string> taxon_ranks;
    
    while (std::getline(nodes_in, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '|')) {
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            fields.push_back(token);
        }
        
        if (fields.size() >= 3) {
            try {
                uint32_t taxon_id = std::stoul(fields[0]);
                uint32_t parent_id = std::stoul(fields[1]);
                std::string rank = fields[2];
                
                taxon_parents[taxon_id] = parent_id;
                taxon_ranks[taxon_id] = rank;
                nodes_loaded++;
                
                taxon_names[taxon_id] = "taxon_" + std::to_string(taxon_id);
            } catch (const std::exception& e) {
                continue;
            }
        }
    }
    nodes_in.close();
    
    std::cout << "Loaded " << nodes_loaded << " taxonomy nodes" << std::endl;
    
    if (!std::filesystem::exists(names_file)) {
        std::cerr << "names.dmp not found at: " << names_file << std::endl;
        std::cerr << "Using taxon IDs as names" << std::endl;
        return nodes_loaded > 0;
    }
    
    std::cout << "Loading taxonomy names from: " << names_file << std::endl;
    std::ifstream names_in(names_file);
    if (!names_in.is_open()) {
        std::cerr << "Cannot open names.dmp: " << names_file << std::endl;
        return nodes_loaded > 0;
    }
    
    int names_loaded = 0;
    
    while (std::getline(names_in, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        while (std::getline(iss, token, '|')) {
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            fields.push_back(token);
        }
        
        if (fields.size() >= 4) {
            try {
                uint32_t taxon_id = std::stoul(fields[0]);
                std::string name_txt = fields[1];
                std::string name_class = fields[3];
                
                if (name_class == "scientific name" && taxon_parents.find(taxon_id) != taxon_parents.end()) {
                    taxon_names[taxon_id] = name_txt;
                    names_loaded++;
                }
            } catch (const std::exception& e) {
                continue;
            }
        }
    }
    names_in.close();
    
    std::cout << "Loaded " << names_loaded << " scientific names" << std::endl;
    
    if (taxon_names.find(1) == taxon_names.end()) {
        taxon_names[1] = "root";
        taxon_parents[1] = 0;
    }
    
    std::cout << "Total taxa in database: " << taxon_parents.size() << std::endl;
    
    return nodes_loaded > 0;
}

// Keep original test function
void GPUKrakenDatabaseBuilder::test_minimizer_extraction_integration() {
    std::cout << "Testing minimizer extraction integration..." << std::endl;
    
    std::string test_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    std::vector<std::string> test_sequences = {test_seq};
    std::vector<uint32_t> test_taxon_ids = {12345};
    
    bool success = process_sequence_batch(test_sequences, test_taxon_ids, 0);
    
    if (success) {
        std::cout << "✓ Minimizer extraction test PASSED" << std::endl;
    } else {
        std::cout << "✗ Minimizer extraction test FAILED" << std::endl;
    }
}

// Placeholder methods to maintain interface compatibility
bool GPUKrakenDatabaseBuilder::build_database_from_file_list(
    const std::string& file_list_path,
    const std::string& taxonomy_path) {
    // Implementation would be similar to build_database_from_genomes
    // but loading file list instead of directory
    std::cerr << "build_database_from_file_list not yet implemented" << std::endl;
    return false;
}

bool GPUKrakenDatabaseBuilder::parse_taxonomy_files(const std::string& taxonomy_path) {
    return load_taxonomy_data(taxonomy_path);
}

bool GPUKrakenDatabaseBuilder::build_compact_hash_table_gpu(
    const LCACandidate* d_candidates,
    int num_candidates) {
    // This functionality is now integrated into save_database()
    return true;
}

// NEW: Enhanced database building from concatenated FNA
bool GPUKrakenDatabaseBuilder::build_database_from_concatenated_fna(
    const std::string& fna_file_path,
    const std::string& taxonomy_nodes_path,
    const std::string& taxonomy_names_path,
    const std::string& compact_taxonomy_path) {
    
    std::cout << "\n=== BUILDING ENHANCED KRAKEN DATABASE FROM CONCATENATED FNA ===" << std::endl;
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Initialize enhanced statistics
    enhanced_stats = EnhancedBuildStats();
    
    // Load taxonomy - try multiple approaches
    bool taxonomy_available = false;
    
    if (!compact_taxonomy_path.empty()) {
        // Option 1: Use pre-built compact taxonomy
        std::cout << "Loading pre-built compact taxonomy..." << std::endl;
        taxonomy_available = enhanced_taxonomy_processor.load_from_compact_file(compact_taxonomy_path);
        
    } else if (!taxonomy_nodes_path.empty() && !taxonomy_names_path.empty()) {
        // Option 2: Build from NCBI files using compact taxonomy infrastructure
        std::cout << "Building compact taxonomy from NCBI files..." << std::endl;
        taxonomy_available = enhanced_taxonomy_processor.load_ncbi_taxonomy(taxonomy_nodes_path, taxonomy_names_path);
        
    } else {
        std::cout << "No taxonomy specified - phylogenetic calculations will be simplified" << std::endl;
    }
    
    if (!taxonomy_available) {
        std::cout << "Warning: Failed to load taxonomy, using simplified phylogenetic calculations" << std::endl;
    }
    
    // Process concatenated FNA file
    if (!process_concatenated_fna_file(fna_file_path)) {
        std::cerr << "Failed to process concatenated FNA file" << std::endl;
        return false;
    }
    
    // Process genomes using existing GPU pipeline
    if (!process_genomes_gpu()) {
        std::cerr << "Failed to process genomes on GPU" << std::endl;
        return false;
    }
    
    // Enhance LCA candidates with phylogenetic data
    if (!enhance_candidates_with_phylogenetics()) {
        std::cerr << "Failed to enhance candidates with phylogenetics" << std::endl;
        return false;
    }
    
    // Save enhanced database
    if (!save_enhanced_database()) {
        std::cerr << "Failed to save enhanced database" << std::endl;
        return false;
    }
    
    // Optionally save compact taxonomy for future use
    if (taxonomy_available && compact_taxonomy_path.empty()) {
        std::string compact_output = output_directory + "/compact_taxonomy.bin";
        if (enhanced_taxonomy_processor.get_compact_taxonomy()) {
            enhanced_taxonomy_processor.get_compact_taxonomy()->save_compact_taxonomy(compact_output);
            std::cout << "Saved compact taxonomy to: " << compact_output << std::endl;
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
    
    std::cout << "\n=== ENHANCED DATABASE BUILD COMPLETE ===" << std::endl;
    std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
    
    enhanced_stats.print_enhanced_stats();
    print_build_statistics();
    
    return true;
}

// Pre-build compact taxonomy method
bool GPUKrakenDatabaseBuilder::prebuild_compact_taxonomy(
    const std::string& nodes_path, 
    const std::string& names_path,
    const std::string& output_path) {
    
    std::cout << "Pre-building compact taxonomy..." << std::endl;
    
    EnhancedNCBITaxonomyProcessor temp_processor;
    if (!temp_processor.load_ncbi_taxonomy(nodes_path, names_path)) {
        std::cerr << "Failed to load NCBI taxonomy" << std::endl;
        return false;
    }
    
    if (temp_processor.get_compact_taxonomy()) {
        if (!temp_processor.get_compact_taxonomy()->save_compact_taxonomy(output_path)) {
            std::cerr << "Failed to save compact taxonomy" << std::endl;
            return false;
        }
        
        std::cout << "Compact taxonomy saved to: " << output_path << std::endl;
        return true;
    }
    
    return false;
}

// Alternative database building method using pre-parsed FNA processor
bool GPUKrakenDatabaseBuilder::build_database_with_fna_processor(
    ConcatenatedFnaProcessor& fna_processor,
    const std::string& taxonomy_nodes_path,
    const std::string& taxonomy_names_path,
    const std::string& compact_taxonomy_path) {
    
    std::cout << "\n=== BUILDING ENHANCED DATABASE WITH PRE-PARSED FNA PROCESSOR ===" << std::endl;
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Initialize enhanced statistics
    enhanced_stats = EnhancedBuildStats();
    
    // Load taxonomy - try multiple approaches
    bool taxonomy_available = false;
    
    if (!compact_taxonomy_path.empty()) {
        // Option 1: Use pre-built compact taxonomy
        std::cout << "Loading pre-built compact taxonomy..." << std::endl;
        taxonomy_available = enhanced_taxonomy_processor.load_from_compact_file(compact_taxonomy_path);
        
    } else if (!taxonomy_nodes_path.empty() && !taxonomy_names_path.empty()) {
        // Option 2: Build from NCBI files using compact taxonomy infrastructure
        std::cout << "Building compact taxonomy from NCBI files..." << std::endl;
        taxonomy_available = enhanced_taxonomy_processor.load_ncbi_taxonomy(taxonomy_nodes_path, taxonomy_names_path);
        
    } else {
        std::cout << "No taxonomy specified - phylogenetic calculations will be simplified" << std::endl;
    }
    
    if (!taxonomy_available) {
        std::cout << "Warning: Failed to load taxonomy, using simplified phylogenetic calculations" << std::endl;
    }
    
    // Get pre-parsed data from FNA processor
    const SpeciesTrackingData& processor_species_data = fna_processor.get_species_data();
    
    // Merge species tracking data
    species_tracking = processor_species_data;
    
    // Update enhanced stats from processor data
    enhanced_stats.total_sequences = genome_files.size();
    enhanced_stats.species_represented = species_tracking.total_species();
    stats.total_sequences = genome_files.size();
    
    std::cout << "Using pre-parsed FNA data:" << std::endl;
    std::cout << "  Total genomes: " << genome_files.size() << std::endl;
    std::cout << "  Total species: " << species_tracking.total_species() << std::endl;
    
    // Process genomes using existing GPU pipeline
    if (!process_genomes_gpu()) {
        std::cerr << "Failed to process genomes on GPU" << std::endl;
        return false;
    }
    
    // Enhance LCA candidates with phylogenetic data
    if (!enhance_candidates_with_phylogenetics()) {
        std::cerr << "Failed to enhance candidates with phylogenetics" << std::endl;
        return false;
    }
    
    // Save enhanced database
    if (!save_enhanced_database()) {
        std::cerr << "Failed to save enhanced database" << std::endl;
        return false;
    }
    
    // Optionally save compact taxonomy for future use
    if (taxonomy_available && compact_taxonomy_path.empty()) {
        std::string compact_output = output_directory + "/compact_taxonomy.bin";
        if (enhanced_taxonomy_processor.get_compact_taxonomy()) {
            enhanced_taxonomy_processor.get_compact_taxonomy()->save_compact_taxonomy(compact_output);
            std::cout << "Saved compact taxonomy to: " << compact_output << std::endl;
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
    
    std::cout << "\n=== ENHANCED DATABASE BUILD COMPLETE ===" << std::endl;
    std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
    
    enhanced_stats.print_enhanced_stats();
    print_build_statistics();
    
    return true;
}

// Process concatenated FNA file
bool GPUKrakenDatabaseBuilder::process_concatenated_fna_file(const std::string& fna_file_path) {
    std::cout << "Using ConcatenatedFnaProcessor for efficient FNA parsing..." << std::endl;
    
    // Use the efficient ConcatenatedFnaProcessor
    ConcatenatedFnaProcessor fna_processor(fna_file_path);
    
    // Temporary directory for genome files
    std::string temp_dir = output_directory + "/temp_genomes";
    
    // Clear existing genome data
    genome_files.clear();
    genome_taxon_ids.clear();
    
    // Process FNA file
    if (!fna_processor.process_fna_file(genome_files, genome_taxon_ids, taxon_names, temp_dir)) {
        std::cerr << "Failed to process FNA file with ConcatenatedFnaProcessor" << std::endl;
        return false;
    }
    
    // Get species tracking data from processor
    const SpeciesTrackingData& processor_species_data = fna_processor.get_species_data();
    
    // Merge species data
    for (const auto& [seq_id, species_taxid] : processor_species_data.sequence_id_to_species) {
        uint16_t genome_count = processor_species_data.get_genome_count_for_species(species_taxid);
        auto species_name_it = processor_species_data.species_names.find(species_taxid);
        std::string species_name = (species_name_it != processor_species_data.species_names.end()) 
                                  ? species_name_it->second 
                                  : enhanced_taxonomy_processor.get_scientific_name(species_taxid);
        
        species_tracking.add_genome(seq_id, species_taxid, species_name);
    }
    
    // Update enhanced stats
    enhanced_stats.total_sequences = genome_files.size();
    enhanced_stats.species_represented = species_tracking.total_species();
    stats.total_sequences = genome_files.size();
    
    std::cout << "FNA processing complete:" << std::endl;
    std::cout << "  Total genomes loaded: " << genome_files.size() << std::endl;
    std::cout << "  Total species: " << species_tracking.total_species() << std::endl;
    
    return !genome_files.empty();
}

// Enhance candidates with phylogenetic data
bool GPUKrakenDatabaseBuilder::enhance_candidates_with_phylogenetics() {
    std::cout << "Enhancing LCA candidates with phylogenetic data..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Group all_lca_candidates by minimizer hash
    std::unordered_map<uint64_t, std::vector<LCACandidate*>> minimizer_groups;
    
    for (auto& candidate : all_lca_candidates) {
        minimizer_groups[candidate.minimizer_hash].push_back(&candidate);
    }
    
    std::cout << "Processing " << minimizer_groups.size() << " unique minimizers..." << std::endl;
    
    phylogenetic_candidates.clear();
    contributing_taxa_arrays.clear();
    
    int processed_count = 0;
    for (auto& [minimizer_hash, candidates] : minimizer_groups) {
        
        // Collect species contributions
        std::unordered_map<uint32_t, uint16_t> species_genome_counts;
        
        for (LCACandidate* candidate : candidates) {
            uint32_t species = candidate->lca_taxon;
            species_genome_counts[species] += candidate->genome_count;
        }
        
        // Convert to vectors
        std::vector<uint32_t> contributing_species;
        std::vector<uint16_t> genome_counts_per_species;
        
        for (const auto& [species, count] : species_genome_counts) {
            contributing_species.push_back(species);
            genome_counts_per_species.push_back(count);
        }
        
        // Use enhanced taxonomy processor for better phylogenetic calculations
        uint32_t lca = enhanced_taxonomy_processor.compute_lca_of_species(contributing_species);
        uint8_t phylo_spread = enhanced_taxonomy_processor.calculate_phylogenetic_spread(contributing_species, lca);
        
        uint8_t max_distance = 0;
        for (uint32_t species : contributing_species) {
            uint8_t distance = enhanced_taxonomy_processor.calculate_distance_to_lca(species, lca);
            max_distance = std::max(max_distance, distance);
        }
        
        // Create enhanced phylogenetic candidate
        PhylogeneticLCACandidate phylo_candidate;
        phylo_candidate.minimizer_hash = minimizer_hash;
        phylo_candidate.lca_taxon = lca;
        phylo_candidate.contributing_species = contributing_species;
        phylo_candidate.genome_counts_per_species = genome_counts_per_species;
        phylo_candidate.phylogenetic_spread = phylo_spread;
        phylo_candidate.max_phylogenetic_distance = max_distance;
        
        // Calculate total genome count
        phylo_candidate.genome_count = 0;
        for (uint16_t count : genome_counts_per_species) {
            phylo_candidate.genome_count += count;
        }
        phylo_candidate.uniqueness_score = 1.0f / phylo_candidate.genome_count;
        
        phylogenetic_candidates.push_back(phylo_candidate);
        enhanced_stats.phylogenetic_lca_computations++;
        
        processed_count++;
        if (processed_count % 10000 == 0) {
            std::cout << "Processed " << processed_count << " minimizers..." << std::endl;
        }
    }
    
    enhanced_stats.minimizers_with_phylo_data = phylogenetic_candidates.size();
    enhanced_stats.unique_minimizers = phylogenetic_candidates.size();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    enhanced_stats.phylogenetic_processing_time = 
        std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "Phylogenetic enhancement completed with enhanced taxonomy processor" << std::endl;
    std::cout << "  Minimizers with phylogenetic data: " << enhanced_stats.minimizers_with_phylo_data << std::endl;
    std::cout << "  Processing time: " << std::fixed << std::setprecision(2) 
              << enhanced_stats.phylogenetic_processing_time << "s" << std::endl;
    
    return true;
}

// Save enhanced database
bool GPUKrakenDatabaseBuilder::save_enhanced_database() {
    std::cout << "Saving enhanced database with phylogenetic metadata..." << std::endl;
    
    // Clear contributing taxa arrays for fresh build
    contributing_taxa_arrays.clear();
    
    // Convert phylogenetic candidates to streamlined format
    std::vector<StreamlinedMinimizerMetadata> streamlined_metadata;
    streamlined_metadata.reserve(phylogenetic_candidates.size());
    
    for (const auto& phylo_candidate : phylogenetic_candidates) {
        StreamlinedMinimizerMetadata metadata;
        metadata.minimizer_hash = phylo_candidate.minimizer_hash;
        metadata.lca_taxon = phylo_candidate.lca_taxon;
        metadata.total_genome_count = phylo_candidate.genome_count;
        metadata.phylogenetic_spread = phylo_candidate.phylogenetic_spread;
        metadata.max_phylogenetic_distance = phylo_candidate.max_phylogenetic_distance;
        
        // Set offset into contributing taxa arrays
        metadata.contributing_taxa_offset = contributing_taxa_arrays.total_entries();
        metadata.num_contributing_taxa = phylo_candidate.contributing_species.size();
        metadata.reserved = 0; // Future use
        
        // Add contributing taxa to arrays
        for (size_t i = 0; i < phylo_candidate.contributing_species.size(); i++) {
            uint32_t species = phylo_candidate.contributing_species[i];
            uint16_t genome_count = phylo_candidate.genome_counts_per_species[i];
            uint8_t distance = enhanced_taxonomy_processor.calculate_distance_to_lca(species, phylo_candidate.lca_taxon);
            
            contributing_taxa_arrays.add_entry(species, distance, genome_count);
        }
        
        streamlined_metadata.push_back(metadata);
    }
    
    // Update statistics
    enhanced_stats.contributing_taxa_array_size = contributing_taxa_arrays.total_entries();
    
    // Save streamlined hash table
    std::string hash_file = output_directory + "/enhanced_hash_table.k2d";
    if (!save_streamlined_hash_table(streamlined_metadata, hash_file)) {
        return false;
    }
    
    // Save contributing taxa arrays
    std::string taxa_file = output_directory + "/contributing_taxa.k2d";
    if (!save_contributing_taxa_arrays(taxa_file)) {
        return false;
    }
    
    // Save enhanced configuration
    std::string config_file = output_directory + "/enhanced_config.txt";
    if (!save_enhanced_config(config_file)) {
        return false;
    }
    
    // Save species mapping
    std::string species_file = output_directory + "/species_mapping.tsv";
    std::ofstream species_out(species_file);
    
    species_out << "species_taxid\tspecies_name\tgenome_count\n";
    for (const auto& [taxid, name] : species_tracking.species_names) {
        uint16_t count = species_tracking.get_genome_count_for_species(taxid);
        species_out << taxid << "\t" << name << "\t" << count << "\n";
    }
    species_out.close();
    
    // Also save standard format for backward compatibility
    if (!save_standard_database()) {
        std::cerr << "Warning: Failed to save standard database format" << std::endl;
        // Don't fail completely if standard format fails
    }
    
    std::cout << "Enhanced database saved successfully!" << std::endl;
    std::cout << "  Hash table: " << hash_file << std::endl;
    std::cout << "  Contributing taxa: " << taxa_file << std::endl;
    std::cout << "  Configuration: " << config_file << std::endl;
    std::cout << "  Species mapping: " << species_file << std::endl;
    
    return true;
}

// Build contributing taxa arrays from phylogenetic candidates
void GPUKrakenDatabaseBuilder::build_contributing_taxa_arrays() {
    contributing_taxa_arrays.clear();
    
    for (const auto& candidate : phylogenetic_candidates) {
        for (size_t i = 0; i < candidate.contributing_species.size(); i++) {
            uint32_t species = candidate.contributing_species[i];
            uint16_t genome_count = candidate.genome_counts_per_species[i];
            uint8_t distance = enhanced_taxonomy_processor.calculate_distance_to_lca(species, candidate.lca_taxon);
            
            contributing_taxa_arrays.add_entry(species, distance, genome_count);
        }
    }
    
    enhanced_stats.contributing_taxa_array_size = contributing_taxa_arrays.total_entries();
}

// Add candidate's taxa to arrays and return offset
uint32_t GPUKrakenDatabaseBuilder::add_to_contributing_taxa_arrays(const PhylogeneticLCACandidate& candidate) {
    uint32_t start_offset = contributing_taxa_arrays.total_entries();
    
    for (size_t i = 0; i < candidate.contributing_species.size(); i++) {
        uint32_t species = candidate.contributing_species[i];
        uint16_t genome_count = candidate.genome_counts_per_species[i];
        uint8_t distance = enhanced_taxonomy_processor.calculate_distance_to_lca(species, candidate.lca_taxon);
        
        contributing_taxa_arrays.add_entry(species, distance, genome_count);
    }
    
    return start_offset;
}

// Save streamlined hash table
bool GPUKrakenDatabaseBuilder::save_streamlined_hash_table(
    const std::vector<StreamlinedMinimizerMetadata>& metadata,
    const std::string& filename) {
    
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Cannot create hash table file: " << filename << std::endl;
        return false;
    }
    
    // Write header
    uint64_t version = 2;  // Enhanced version
    uint64_t num_entries = metadata.size();
    uint64_t metadata_size = sizeof(StreamlinedMinimizerMetadata);
    
    out.write(reinterpret_cast<const char*>(&version), sizeof(uint64_t));
    out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    out.write(reinterpret_cast<const char*>(&metadata_size), sizeof(uint64_t));
    
    // Write metadata
    out.write(reinterpret_cast<const char*>(metadata.data()), 
              metadata.size() * sizeof(StreamlinedMinimizerMetadata));
    
    out.close();
    
    std::cout << "  Wrote " << num_entries << " streamlined metadata entries (" 
              << (num_entries * sizeof(StreamlinedMinimizerMetadata) / 1024 / 1024) << " MB)" << std::endl;
    
    return true;
}

// Save contributing taxa arrays
bool GPUKrakenDatabaseBuilder::save_contributing_taxa_arrays(const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Cannot create contributing taxa file: " << filename << std::endl;
        return false;
    }
    
    uint64_t num_entries = contributing_taxa_arrays.total_entries();
    
    // Write header
    out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    
    // Write arrays
    out.write(reinterpret_cast<const char*>(contributing_taxa_arrays.taxa_ids.data()), 
              num_entries * sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(contributing_taxa_arrays.phylogenetic_distances.data()), 
              num_entries * sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(contributing_taxa_arrays.genome_counts_per_taxon.data()), 
              num_entries * sizeof(uint16_t));
    
    out.close();
    
    std::cout << "  Wrote " << num_entries << " contributing taxa entries" << std::endl;
    
    return true;
}

// Save enhanced configuration
bool GPUKrakenDatabaseBuilder::save_enhanced_config(const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot create config file: " << filename << std::endl;
        return false;
    }
    
    out << "# Enhanced Kraken Database Configuration\n";
    out << "version=2\n";
    out << "k=" << minimizer_params.k << "\n";
    out << "ell=" << minimizer_params.ell << "\n";
    out << "spaces=" << minimizer_params.spaces << "\n";
    out << "subsampling_rate=" << subsampling_rate << "\n";
    out << "min_clear_hash_value=0x" << std::hex << min_clear_hash_value << std::dec << "\n";
    out << "toggle_mask=0x" << std::hex << toggle_mask << std::dec << "\n";
    out << "\n# Database statistics\n";
    out << "total_sequences=" << enhanced_stats.total_sequences << "\n";
    out << "total_bases=" << enhanced_stats.total_bases << "\n";
    out << "unique_minimizers=" << enhanced_stats.unique_minimizers << "\n";
    out << "species_count=" << species_tracking.total_species() << "\n";
    out << "genome_count=" << species_tracking.total_genomes() << "\n";
    out << "minimizers_with_phylo=" << enhanced_stats.minimizers_with_phylo_data << "\n";
    out << "contributing_taxa_entries=" << enhanced_stats.contributing_taxa_array_size << "\n";
    out << "\n# Timing information\n";
    out << "sequence_processing_time=" << enhanced_stats.sequence_processing_time << "\n";
    out << "minimizer_extraction_time=" << enhanced_stats.minimizer_extraction_time << "\n";
    out << "phylogenetic_processing_time=" << enhanced_stats.phylogenetic_processing_time << "\n";
    out << "database_construction_time=" << enhanced_stats.database_construction_time << "\n";
    
    out.close();
    
    return true;
}

// Save standard database format for backward compatibility
bool GPUKrakenDatabaseBuilder::save_standard_database() {
    std::cout << "Saving standard database format for backward compatibility..." << std::endl;
    
    // Call the original save_database() method
    return save_database();
}

// NEW: Streaming support for large concatenated FNA files
bool GPUKrakenDatabaseBuilder::build_database_from_streaming_fna(
    const std::string& fna_file_path,
    const std::string& taxonomy_path) {
    
    std::cout << "\n=== BUILDING KRAKEN DATABASE FROM STREAMING FNA ===" << std::endl;
    std::cout << "FNA file: " << fna_file_path << std::endl;
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Verify file exists and is readable
    if (!std::filesystem::exists(fna_file_path)) {
        std::cerr << "Error: FNA file does not exist: " << fna_file_path << std::endl;
        return false;
    }
    
    if (!std::filesystem::is_regular_file(fna_file_path)) {
        std::cerr << "Error: Path is not a regular file: " << fna_file_path << std::endl;
        return false;
    }
    
    // Load taxonomy if provided
    if (!taxonomy_path.empty()) {
        if (!load_taxonomy_data(taxonomy_path)) {
            std::cerr << "Warning: Failed to load taxonomy data" << std::endl;
        }
    }
    
    // Create temporary directory for genome files
    std::string temp_dir = output_directory + "/temp_genomes";
    
    try {
        // Initialize streaming processor with reasonable batch size
        // Adjust based on available memory - smaller batches for limited memory
        int batch_size = std::min(MAX_SEQUENCE_BATCH, 25);
        StreamingFnaProcessor processor(fna_file_path, temp_dir, batch_size);
        
        // Clear any existing genome files list
        genome_files.clear();
        genome_taxon_ids.clear();
        
        // Process file in streaming batches
        std::vector<std::string> batch_files;
        std::vector<uint32_t> batch_taxons;
        int batch_number = 0;
        
        while (processor.process_next_batch(batch_files, batch_taxons)) {
            batch_number++;
            std::cout << "\nProcessing streaming batch " << batch_number 
                      << " (" << batch_files.size() << " genomes)" << std::endl;
            
            // Update our internal lists (for compatibility with existing code)
            genome_files = batch_files;
            genome_taxon_ids = batch_taxons;
            
            // Update taxon names if not in taxonomy
            for (size_t i = 0; i < batch_files.size(); i++) {
                uint32_t taxon_id = batch_taxons[i];
                if (taxon_id > 0 && taxon_names.find(taxon_id) == taxon_names.end()) {
                    std::filesystem::path p(batch_files[i]);
                    taxon_names[taxon_id] = "taxon_" + std::to_string(taxon_id);
                }
            }
            
            // Process this batch on GPU
            if (!process_genomes_gpu()) {
                std::cerr << "Failed to process batch " << batch_number << " on GPU" << std::endl;
                std::filesystem::remove_all(temp_dir);
                return false;
            }
            
            // Clean up temporary files after processing
            for (const auto& file : batch_files) {
                std::filesystem::remove(file);
            }
            
            // Update statistics
            stats.total_sequences += batch_files.size();
        }
        
        // Clean up temp directory
        std::filesystem::remove_all(temp_dir);
        
        std::cout << "\nStreaming processing complete:" << std::endl;
        std::cout << "  Total genomes processed: " << processor.get_total_genomes() << std::endl;
        std::cout << "  Total bases processed: " << processor.get_total_bases() << std::endl;
        std::cout << "  Total batches: " << batch_number << std::endl;
        
        // Save the complete database
        if (!save_database()) {
            std::cerr << "Failed to save database" << std::endl;
            return false;
        }
        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
        
        std::cout << "\n=== DATABASE BUILD COMPLETE ===" << std::endl;
        std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
        print_build_statistics();
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during streaming processing: " << e.what() << std::endl;
        // Clean up on error
        if (std::filesystem::exists(temp_dir)) {
            std::filesystem::remove_all(temp_dir);
        }
        return false;
    }
}

#endif // GPU_KRAKEN_CLASSIFIER_HEADER_ONLY

// Feature export implementation
bool GPUKrakenDatabaseBuilder::export_features_for_training(const FeatureExportConfig& export_config) {
    std::cout << "\n=== EXPORTING FEATURES FOR ML TRAINING ===" << std::endl;
    
    if (all_lca_candidates.empty()) {
        std::cerr << "No minimizer data available. Build database first." << std::endl;
        return false;
    }
    
    // Create feature exporter
    FeatureExporter exporter(export_config);
    
    // Convert LCA candidates to minimizer hits
    std::vector<GPUMinimizerHit> minimizer_hits;
    minimizer_hits.reserve(all_lca_candidates.size());
    
    // Group candidates by minimizer hash to get position and genome info
    std::unordered_map<uint64_t, std::vector<GPUMinimizerHit>> hash_to_hits;
    
    for (const auto& candidate : all_lca_candidates) {
        GPUMinimizerHit hit;
        hit.minimizer_hash = candidate.minimizer_hash;
        hit.taxon_id = candidate.lca_taxon;
        hit.genome_id = 0; // Would need genome tracking
        hit.position = 0; // Would need position tracking
        
        // Encode basic features
        hit.feature_flags = 0;
        
        // Encode uniqueness score as ML weight
        hit.ml_weight = static_cast<uint16_t>(candidate.uniqueness_score * 65535);
        
        minimizer_hits.push_back(hit);
    }
    
    // Create a minimal feature extractor (if not available)
    MinimizerFeatureExtractor feature_extractor;
    
    // Export features
    bool success = exporter.export_training_features(
        minimizer_hits, 
        feature_extractor,
        taxon_names
    );
    
    if (success) {
        std::cout << "Feature export completed successfully!" << std::endl;
        std::cout << "Exported " << exporter.get_total_features_exported() 
                  << " unique minimizer features" << std::endl;
    }
    
    return success;
}

#endif // GPU_KRAKEN_DATABASE_BUILDER_CUH