// gpu_kraken_database_builder.cu
// COMPLETE REPLACEMENT: Kraken2-inspired GPU database builder with all original functionality
// Fixes the over-aggressive deduplication while maintaining all features

#ifndef GPU_KRAKEN_DATABASE_BUILDER_CUH
#define GPU_KRAKEN_DATABASE_BUILDER_CUH

#include "gpu_kraken_classifier.cu"
#include "gpu_minimizer_extraction.cuh"
#include "streaming_fna_processor.h"
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
    
    // Configuration - declarations only (implementations moved outside class)
    void set_batch_size(int sequences_per_batch);
    void set_minimizer_capacity(int minimizers_per_batch);
    void enable_auto_memory_scaling(bool enable = true, size_t memory_fraction = 80);
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
    bool save_enhanced_config(const std::string& filename);
    bool save_standard_database();
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

// Constructor - keep original but add new features
GPUKrakenDatabaseBuilder::GPUKrakenDatabaseBuilder(
    const std::string& output_dir,
    const ClassificationParams& config)
    : output_directory(output_dir), params(config),
      d_sequence_data(nullptr), d_genome_info(nullptr),
      d_minimizer_hits(nullptr), d_minimizer_counts(nullptr),
      d_lca_candidates(nullptr), d_taxonomy_nodes(nullptr) {
    
    // Convert classification params to minimizer params
    minimizer_params.k = params.k;
    minimizer_params.ell = params.ell;
    minimizer_params.spaces = params.spaces;
    minimizer_params.xor_mask = 0x3c8bfbb395c60474ULL;  // Kraken2-style XOR mask
    
    std::cout << "Initializing GPU Kraken database builder (High Capacity)..." << std::endl;
    std::cout << "Parameters: k=" << minimizer_params.k 
              << ", ell=" << minimizer_params.ell 
              << ", spaces=" << minimizer_params.spaces << std::endl;
    
    // Create output directory
    std::filesystem::create_directories(output_directory);
    
    // Enable auto-scaling by default for better memory utilization
    enable_auto_memory_scaling(true, 80);
    
    // Initialize statistics
    memset(&stats, 0, sizeof(stats));
    total_minimizers_capped = 0;
    total_batches_capped = 0;
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

// UPDATED: Smarter memory management
void GPUKrakenDatabaseBuilder::check_and_adjust_memory() {
    size_t free_mem, total_mem;
    cudaError_t cuda_status = cudaMemGetInfo(&free_mem, &total_mem);
    if (cuda_status != cudaSuccess) {
        std::cerr << "CUDA memory query failed: " << cudaGetErrorString(cuda_status) << std::endl;
        return;
    }
    
    std::cout << "GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
              << (total_mem / 1024 / 1024) << " MB total" << std::endl;
    
    if (auto_scale_memory) {
        // Calculate optimal batch sizes based on available memory
        size_t usable_memory = (total_mem * max_gpu_memory_usage_fraction) / 100;
        
        // Memory breakdown:
        size_t sequence_memory_per_genome = 10000000; // ~10MB avg per genome
        size_t minimizer_hit_size = sizeof(GPUMinimizerHit); // 24 bytes
        size_t lca_candidate_size = sizeof(LCACandidate);     // 20 bytes
        size_t genome_info_size = sizeof(GPUGenomeInfo);      // 16 bytes
        
        // Reserve memory for other allocations (20% buffer)
        size_t available_for_minimizers = usable_memory * 0.6; // 60% for minimizers
        size_t available_for_sequences = usable_memory * 0.2;  // 20% for sequences
        
        // Calculate optimal minimizer capacity
        size_t memory_per_minimizer = minimizer_hit_size + lca_candidate_size;
        int optimal_minimizer_capacity = available_for_minimizers / memory_per_minimizer;
        
        // Calculate optimal sequence batch size
        int optimal_sequence_batch = available_for_sequences / sequence_memory_per_genome;
        
        // Apply reasonable limits
        optimal_minimizer_capacity = std::min(optimal_minimizer_capacity, (int)50000000); // Max 50M
        optimal_minimizer_capacity = std::max(optimal_minimizer_capacity, (int)1000000);  // Min 1M
        
        optimal_sequence_batch = std::min(optimal_sequence_batch, 100);  // Max 100 genomes
        optimal_sequence_batch = std::max(optimal_sequence_batch, 5);    // Min 5 genomes
        
        MAX_MINIMIZERS_PER_BATCH = optimal_minimizer_capacity;
        MAX_SEQUENCE_BATCH = optimal_sequence_batch;
        
        std::cout << "Auto-scaled batch sizes:" << std::endl;
        std::cout << "  Minimizers per batch: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
        std::cout << "  Sequences per batch: " << MAX_SEQUENCE_BATCH << std::endl;
        std::cout << "  Estimated memory usage: " << (usable_memory / 1024 / 1024) << " MB" << std::endl;
    } else {
        // Manual validation of current settings
        size_t sequence_memory = MAX_SEQUENCE_BATCH * 10000000; 
        size_t minimizer_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit);
        size_t candidate_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(LCACandidate);
        size_t genome_info_memory = MAX_SEQUENCE_BATCH * sizeof(GPUGenomeInfo);
        size_t total_needed = sequence_memory + minimizer_memory + candidate_memory + genome_info_memory;
        
        std::cout << "Memory requirements (manual settings):" << std::endl;
        std::cout << "  Sequence data: " << (sequence_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "  Minimizer hits: " << (minimizer_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "  LCA candidates: " << (candidate_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "  Total required: " << (total_needed / 1024 / 1024) << " MB" << std::endl;
        
        size_t safe_memory = free_mem * 0.8;
        
        if (total_needed > safe_memory) {
            std::cout << "WARNING: Current settings may exceed available GPU memory!" << std::endl;
            std::cout << "Consider enabling auto-scaling or reducing batch sizes." << std::endl;
            
            // Auto-adjust if memory requirements are too high
            double scale = (double)safe_memory / total_needed;
            MAX_MINIMIZERS_PER_BATCH = (int)(MAX_MINIMIZERS_PER_BATCH * scale);
            MAX_SEQUENCE_BATCH = (int)(MAX_SEQUENCE_BATCH * scale);
            
            std::cout << "Auto-adjusted to fit memory:" << std::endl;
            std::cout << "  Minimizers: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
            std::cout << "  Sequences: " << MAX_SEQUENCE_BATCH << std::endl;
        }
    }
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

void GPUKrakenDatabaseBuilder::enable_auto_memory_scaling(bool enable, size_t memory_fraction) {
    auto_scale_memory = enable;
    max_gpu_memory_usage_fraction = memory_fraction;
    if (enable) {
        std::cout << "Enabled auto memory scaling (using " << memory_fraction 
                  << "% of GPU memory)" << std::endl;
        check_and_adjust_memory();
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

// Keep original main pipeline method
bool GPUKrakenDatabaseBuilder::build_database_from_genomes(
    const std::string& genome_library_path,
    const std::string& taxonomy_path) {
    
    std::cout << "\n=== BUILDING KRAKEN DATABASE FROM GENOMES (Kraken2-inspired) ===" << std::endl;
    std::cout << "Input genome_library_path: '" << genome_library_path << "'" << std::endl;
    auto total_start = std::chrono::high_resolution_clock::now();
    
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
}

// Keep original file loading logic
bool GPUKrakenDatabaseBuilder::load_genome_files(const std::string& library_path) {
    std::cout << "Loading genome files from: " << library_path << std::endl;
    
    // Validate the library path
    if (library_path.empty()) {
        std::cerr << "Error: Genome library path is empty!" << std::endl;
        return false;
    }
    
    if (!std::filesystem::exists(library_path)) {
        std::cerr << "Error: Genome library path does not exist: " << library_path << std::endl;
        return false;
    }
    
    if (!std::filesystem::is_directory(library_path)) {
        std::cerr << "Error: Genome library path is not a directory: " << library_path << std::endl;
        return false;
    }
    
    genome_files = find_genome_files(library_path);
    
    if (genome_files.empty()) {
        std::cerr << "No genome files found in " << library_path << std::endl;
        return false;
    }
    
    genome_taxon_ids.reserve(genome_files.size());
    for (const auto& file : genome_files) {
        uint32_t taxon_id = extract_taxon_from_filename(file);
        genome_taxon_ids.push_back(taxon_id);
        
        if (taxon_names.find(taxon_id) == taxon_names.end()) {
            std::filesystem::path p(file);
            taxon_names[taxon_id] = p.stem().string();
        }
    }
    
    std::cout << "Loaded " << genome_files.size() << " genome files" << std::endl;
    stats.total_sequences = genome_files.size();
    
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
    
    // Pre-calculate total size to prevent buffer overflow
    size_t total_sequence_length = 0;
    for (const auto& seq : sequences) {
        total_sequence_length += seq.length();
        
        // Skip sequences that are too short or too long
        if (seq.length() < minimizer_params.k) {
            continue; // Skip, but don't fail
        }
        if (seq.length() > 50000000) { // 50MB limit per sequence
            std::cerr << "Warning: Skipping very large sequence (" 
                      << seq.length() << " bp)" << std::endl;
            continue;
        }
    }
    
    // Check if total size exceeds our allocated memory
    size_t max_sequence_memory = size_t(MAX_SEQUENCE_BATCH) * 10000000;
    if (total_sequence_length > max_sequence_memory) {
        std::cerr << "ERROR: Total sequence length (" << total_sequence_length 
                  << ") exceeds allocated memory (" << max_sequence_memory << ")" << std::endl;
        return false;
    }
    
    // Build genome info and concatenated sequences with bounds checking
    std::string concatenated_sequences;
    std::vector<GPUGenomeInfo> genome_infos;
    
    concatenated_sequences.reserve(total_sequence_length);
    genome_infos.reserve(sequences.size());
    
    uint32_t current_offset = 0;
    for (size_t i = 0; i < sequences.size(); i++) {
        if (sequences[i].length() < minimizer_params.k) continue;
        if (sequences[i].length() > 50000000) continue;
        
        GPUGenomeInfo info;
        info.taxon_id = taxon_ids[i];
        info.sequence_offset = current_offset;
        info.sequence_length = sequences[i].length();
        info.genome_id = batch_offset + i;
        
        genome_infos.push_back(info);
        concatenated_sequences += sequences[i];
        current_offset += sequences[i].length();
        
        // Update stats safely
        stats.total_kmers_processed += sequences[i].length() - minimizer_params.k + 1;
    }
    
    if (genome_infos.empty()) {
        std::cout << "No valid sequences in this batch, skipping" << std::endl;
        return true;
    }
    
    // Final size validation before GPU copy
    if (concatenated_sequences.length() > max_sequence_memory) {
        std::cerr << "ERROR: Final concatenated size exceeds memory limit" << std::endl;
        return false;
    }
    
    // Copy to GPU with error checking
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
                               genome_infos.size(), d_minimizer_hits, &total_hits)) {
        return false;
    }
    
    int num_candidates = 0;
    if (!compute_lca_assignments_gpu(d_minimizer_hits, total_hits,
                                    d_lca_candidates, &num_candidates)) {
        return false;
    }
    
    if (num_candidates > 0) {
        std::vector<LCACandidate> batch_candidates(num_candidates);
        err = cudaMemcpy(batch_candidates.data(), d_lca_candidates,
                        num_candidates * sizeof(LCACandidate), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            std::cerr << "ERROR: Failed to copy candidates from GPU: " 
                      << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        all_lca_candidates.insert(all_lca_candidates.end(), 
                                 batch_candidates.begin(), batch_candidates.end());
        
        std::cout << "Accumulated " << num_candidates << " candidates (total: " 
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

// Keep all original utility functions
bool GPUKrakenDatabaseBuilder::allocate_gpu_memory() {
    std::cout << "Allocating GPU memory..." << std::endl;
    
    size_t free_memory, total_memory_gpu;
    cudaMemGetInfo(&free_memory, &total_memory_gpu);
    std::cout << "GPU Memory: " << (free_memory / 1024 / 1024) << " MB free / " 
              << (total_memory_gpu / 1024 / 1024) << " MB total" << std::endl;
    
    // Validate parameters before allocation
    if (MAX_SEQUENCE_BATCH <= 0 || MAX_SEQUENCE_BATCH > 1000) {
        std::cerr << "ERROR: Invalid sequence batch size: " << MAX_SEQUENCE_BATCH << std::endl;
        return false;
    }
    
    if (MAX_MINIMIZERS_PER_BATCH <= 0 || MAX_MINIMIZERS_PER_BATCH > 100000000) {
        std::cerr << "ERROR: Invalid minimizer capacity: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
        return false;
    }
    
    size_t sequence_memory = size_t(MAX_SEQUENCE_BATCH) * 10000000;
    size_t genome_info_memory = MAX_SEQUENCE_BATCH * sizeof(GPUGenomeInfo);
    size_t minimizer_memory = size_t(MAX_MINIMIZERS_PER_BATCH) * sizeof(GPUMinimizerHit);
    size_t candidate_memory = size_t(MAX_MINIMIZERS_PER_BATCH) * sizeof(LCACandidate);
    
    size_t total_required = sequence_memory + genome_info_memory + 
                           minimizer_memory + candidate_memory;
    
    std::cout << "Memory requirements:" << std::endl;
    std::cout << "  Sequence data: " << (sequence_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Minimizer hits: " << (minimizer_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  LCA candidates: " << (candidate_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  Total required: " << (total_required / 1024 / 1024) << " MB" << std::endl;
    
    if (total_required > free_memory * 0.9) {
        std::cerr << "ERROR: Not enough GPU memory!" << std::endl;
        std::cerr << "Required: " << (total_required / 1024 / 1024) << " MB" << std::endl;
        std::cerr << "Available: " << (free_memory * 0.9 / 1024 / 1024) << " MB" << std::endl;
        return false;
    }
    
    CUDA_CHECK(cudaMalloc(&d_sequence_data, sequence_memory));
    CUDA_CHECK(cudaMalloc(&d_genome_info, genome_info_memory));
    CUDA_CHECK(cudaMalloc(&d_minimizer_hits, minimizer_memory));
    CUDA_CHECK(cudaMalloc(&d_minimizer_counts, MAX_SEQUENCE_BATCH * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_lca_candidates, candidate_memory));
    
    std::cout << "Successfully allocated " << (total_required / 1024 / 1024) 
              << " MB of GPU memory" << std::endl;
    
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
    
    // Validate directory before attempting to iterate
    if (directory.empty() || !std::filesystem::exists(directory) || !std::filesystem::is_directory(directory)) {
        std::cerr << "Error: Invalid directory for find_genome_files: '" << directory << "'" << std::endl;
        return files;
    }
    
    try {
        for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
            if (entry.is_regular_file()) {
                std::string extension = entry.path().extension().string();
                
                if (extension == ".fna" || extension == ".fa" || extension == ".fasta" ||
                    extension == ".ffn" || extension == ".faa") {
                    files.push_back(entry.path().string());
                }
            }
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error in find_genome_files: " << e.what() << std::endl;
        std::cerr << "Directory: '" << directory << "'" << std::endl;
        return files;
    }
    
    std::sort(files.begin(), files.end());
    return files;
}

uint32_t GPUKrakenDatabaseBuilder::extract_taxon_from_filename(const std::string& filename) {
    std::filesystem::path p(filename);
    std::string stem = p.stem().string();
    
    std::regex taxid_pattern(R"(taxid[_-](\d+))");
    std::smatch match;
    if (std::regex_search(stem, match, taxid_pattern)) {
        return std::stoul(match[1].str());
    }
    
    std::regex gcf_pattern(R"(GCF_(\d+\.\d+))");
    if (std::regex_search(stem, match, gcf_pattern)) {
        std::hash<std::string> hasher;
        return hasher(match[1].str()) % 1000000 + 1000000;
    }
    
    std::hash<std::string> hasher;
    return hasher(stem) % 1000000 + 2000000;
}

std::vector<std::string> GPUKrakenDatabaseBuilder::load_sequences_from_fasta(const std::string& fasta_path) {
    std::vector<std::string> sequences;
    std::ifstream file(fasta_path);
    
    if (!file.is_open()) {
        std::cerr << "Cannot open FASTA file: " << fasta_path << std::endl;
        return sequences;
    }
    
    std::string line, current_sequence;
    bool in_sequence = false;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (in_sequence && !current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
            in_sequence = true;
        } else if (in_sequence) {
            current_sequence += line;
        }
    }
    
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }
    
    file.close();
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

#endif // GPU_KRAKEN_DATABASE_BUILDER_CUH