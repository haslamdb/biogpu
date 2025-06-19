// gpu_kraken_types.h
// Common type definitions for GPU Kraken classifier and database builder

#ifndef GPU_KRAKEN_TYPES_H
#define GPU_KRAKEN_TYPES_H

#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <set>
#include <iostream>
#include <iomanip>

// LCA candidate structure for database building
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

// GPU build statistics
struct GPUBuildStats {
    uint64_t total_sequences = 0;
    uint64_t total_bases = 0;
    uint64_t total_minimizers = 0;
    uint64_t unique_minimizers = 0;
    uint64_t total_genomes = 0;
    uint64_t cuda_kernel_time_ms = 0;
    uint64_t host_processing_time_ms = 0;
    double database_construction_time = 0.0;
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

// Enhanced NCBI taxonomy processor (placeholder for now)
struct EnhancedNCBITaxonomyProcessor {
    // This would be implemented to use the compact taxonomy tool
};

#endif // GPU_KRAKEN_TYPES_H