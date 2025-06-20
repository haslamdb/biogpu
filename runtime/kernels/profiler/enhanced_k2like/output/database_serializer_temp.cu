// output/database_serializer.cu
// Implementation of database serialization and output formatting
// Extracted from monolithic database builder

#include "database_serializer.h"
#include "../gpu_kraken_types.h"
#include "../taxonomy/taxonomy_processor.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <chrono>

// Simple host-side LCA computation (simplified version)
inline uint32_t compute_simple_lca_host(uint32_t taxon1, uint32_t taxon2) {
    // For now, return the minimum taxon as a simple heuristic
    // In a real implementation, this would traverse the taxonomy tree
    if (taxon1 == 0) return taxon2;
    if (taxon2 == 0) return taxon1;
    return std::min(taxon1, taxon2);
}

// ===========================
// StandardDatabaseSerializer Implementation
// ===========================

StandardDatabaseSerializer::StandardDatabaseSerializer(const std::string& output_dir) 
    : output_directory_(output_dir) {
    
    // Initialize metadata with defaults
    metadata_ = DatabaseMetadata();
    
    // Create output directory if it doesn't exist
    create_output_directory();
}

bool StandardDatabaseSerializer::save_standard_database(
    const std::vector<LCACandidate>& candidates,
    const std::unordered_map<uint32_t, std::string>& taxon_names,
    const std::unordered_map<uint32_t, uint32_t>& taxon_parents,
    const GPUBuildStats& stats) {
    
    std::cout << "\n=== SAVING STANDARD KRAKEN DATABASE ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Update metadata from stats
    metadata_.total_sequences = stats.total_sequences;
    metadata_.total_bases = stats.total_bases;
    metadata_.unique_minimizers = stats.unique_minimizers;
    
    // Process candidates and merge duplicates
    std::unordered_map<uint64_t, LCACandidate> unique_candidates;
    unique_candidates.reserve(candidates.size());
    
    std::cout << "Processing " << candidates.size() << " candidates..." << std::endl;
    
    for (const auto& candidate : candidates) {
        auto it = unique_candidates.find(candidate.minimizer_hash);
        if (it == unique_candidates.end()) {
            unique_candidates[candidate.minimizer_hash] = candidate;
        } else {
            // Merge duplicates
            auto& existing = it->second;
            existing.genome_count++;
            existing.lca_taxon = compute_simple_lca_host(existing.lca_taxon, candidate.lca_taxon);
            existing.uniqueness_score = 1.0f / existing.genome_count;
        }
    }
    
    // Convert back to vector for saving
    std::vector<LCACandidate> merged_candidates;
    merged_candidates.reserve(unique_candidates.size());
    for (const auto& [hash, candidate] : unique_candidates) {
        merged_candidates.push_back(candidate);
    }
    
    std::cout << "Merged to " << merged_candidates.size() << " unique minimizers" << std::endl;
    
    // Save individual components
    if (!save_hash_table(merged_candidates)) {
        std::cerr << "Failed to save hash table" << std::endl;
        return false;
    }
    
    if (!save_taxonomy_file(taxon_names, taxon_parents)) {
        std::cerr << "Failed to save taxonomy file" << std::endl;
        return false;
    }
    
    if (!save_config_file(stats)) {
        std::cerr << "Failed to save config file" << std::endl;
        return false;
    }
    
    if (!save_build_statistics(stats)) {
        std::cerr << "Failed to save build statistics" << std::endl;
        return false;
    }
    
    // Generate checksums and validate
    if (!generate_database_checksum()) {
        std::cerr << "Warning: Failed to generate database checksum" << std::endl;
    }
    
    if (!validate_database_consistency()) {
        std::cerr << "Warning: Database consistency validation failed" << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "\n=== DATABASE SAVE COMPLETE ===" << std::endl;
    std::cout << "Save time: " << std::fixed << std::setprecision(2) << duration << " seconds" << std::endl;
    std::cout << "Hash table: " << get_hash_table_path() << std::endl;
    std::cout << "Taxonomy: " << get_taxonomy_path() << std::endl;
    std::cout << "Config: " << get_config_path() << std::endl;
    
    return true;
}

bool StandardDatabaseSerializer::save_hash_table(const std::vector<LCACandidate>& candidates) {
    std::string hash_file = get_hash_table_path();
    std::cout << "Writing hash table to: " << hash_file << std::endl;
    
    std::ofstream hash_out(hash_file, std::ios::binary);
    if (!hash_out.is_open()) {
        std::cerr << "Cannot create hash table file: " << hash_file << std::endl;
        return false;
    }
    
    // Write header
    uint64_t table_size = candidates.size() * 2; // Kraken2 format
    uint64_t num_entries = candidates.size();
    
    hash_out.write(reinterpret_cast<const char*>(&table_size), sizeof(uint64_t));
    hash_out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    
    if (hash_out.fail()) {
        std::cerr << "Failed to write hash table header" << std::endl;
        return false;
    }
    
    // Write entries with progress reporting
    size_t entries_written = 0;
    for (const auto& candidate : candidates) {
        hash_out.write(reinterpret_cast<const char*>(&candidate.minimizer_hash), sizeof(uint64_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.lca_taxon), sizeof(uint32_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.genome_count), sizeof(uint32_t));
        hash_out.write(reinterpret_cast<const char*>(&candidate.uniqueness_score), sizeof(float));
        
        if (hash_out.fail()) {
            std::cerr << "Failed to write hash table entry " << entries_written << std::endl;
            return false;
        }
        
        entries_written++;
        if (entries_written % 100000 == 0) {
            std::cout << "  Written " << entries_written << "/" << num_entries << " entries" << std::endl;
        }
    }
    
    hash_out.close();
    
    // Update metadata
    metadata_.hash_table_file = hash_file;
    metadata_.unique_minimizers = num_entries;
    
    std::cout << "✓ Hash table saved: " << num_entries << " entries (" 
              << (num_entries * sizeof(LCACandidate) / 1024 / 1024) << " MB)" << std::endl;
    
    return true;
}

bool StandardDatabaseSerializer::save_taxonomy_file(
    const std::unordered_map<uint32_t, std::string>& taxon_names,
    const std::unordered_map<uint32_t, uint32_t>& taxon_parents) {
    
    std::string taxonomy_file = get_taxonomy_path();
    std::cout << "Writing taxonomy to: " << taxonomy_file << std::endl;
    
    std::ofstream tax_out(taxonomy_file);
    if (!tax_out.is_open()) {
        std::cerr << "Cannot create taxonomy file: " << taxonomy_file << std::endl;
        return false;
    }
    
    // Write header
    tax_out << "taxon_id\tname\tparent_id\n";
    
    // Write entries
    size_t entries_written = 0;
    for (const auto& [taxon_id, name] : taxon_names) {
        uint32_t parent_id = 0;
        auto parent_it = taxon_parents.find(taxon_id);
        if (parent_it != taxon_parents.end()) {
            parent_id = parent_it->second;
        }
        
        tax_out << taxon_id << "\t" << name << "\t" << parent_id << "\n";
        
        if (tax_out.fail()) {
            std::cerr << "Failed to write taxonomy entry for taxon " << taxon_id << std::endl;
            return false;
        }
        
        entries_written++;
    }
    
    tax_out.close();
    
    // Update metadata
    metadata_.taxonomy_file = taxonomy_file;
    
    std::cout << "✓ Taxonomy saved: " << entries_written << " taxa" << std::endl;
    
    return true;
}

bool StandardDatabaseSerializer::save_config_file(const GPUBuildStats& stats) {
    std::string config_file = get_config_path();
    std::cout << "Writing config to: " << config_file << std::endl;
    
    std::ofstream config_out(config_file);
    if (!config_out.is_open()) {
        std::cerr << "Cannot create config file: " << config_file << std::endl;
        return false;
    }
    
    // Write configuration parameters
    config_out << "# Kraken Database Configuration\n";
    config_out << "version=1\n";
    config_out << "k=" << metadata_.k_value << "\n";
    config_out << "ell=" << metadata_.ell_value << "\n";
    config_out << "spaces=" << metadata_.spaces_value << "\n";
    config_out << "subsampling_rate=" << metadata_.subsampling_rate << "\n";
    config_out << "min_clear_hash_value=0x" << std::hex << metadata_.min_clear_hash_value << std::dec << "\n";
    config_out << "toggle_mask=0x" << std::hex << metadata_.toggle_mask << std::dec << "\n";
    
    config_out << "\n# Database Statistics\n";
    config_out << "unique_minimizers=" << metadata_.unique_minimizers << "\n";
    config_out << "total_sequences=" << metadata_.total_sequences << "\n";
    config_out << "total_bases=" << metadata_.total_bases << "\n";
    
    config_out << "\n# Build Information\n";
    config_out << "build_time=" << std::fixed << std::setprecision(2) << metadata_.build_time << "\n";
    config_out << "processing_rate=" << std::scientific << std::setprecision(2) << metadata_.processing_rate << "\n";
    
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    config_out << "build_date=" << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S") << "\n";
    
    config_out.close();
    
    // Update metadata
    metadata_.config_file = config_file;
    
    std::cout << "✓ Config file saved" << std::endl;
    
    return true;
}

bool StandardDatabaseSerializer::save_build_statistics(const GPUBuildStats& stats) {
    std::string stats_file = get_stats_path();
    std::cout << "Writing build statistics to: " << stats_file << std::endl;
    
    std::ofstream stats_out(stats_file);
    if (!stats_out.is_open()) {
        std::cerr << "Cannot create statistics file: " << stats_file << std::endl;
        return false;
    }
    
    stats_out << "# Build Statistics\n";
    stats_out << "total_sequences=" << stats.total_sequences << "\n";
    stats_out << "total_bases=" << stats.total_bases << "\n";
    stats_out << "total_kmers_processed=" << stats.total_kmers_processed << "\n";
    stats_out << "valid_minimizers_extracted=" << stats.valid_minimizers_extracted << "\n";
    stats_out << "unique_minimizers=" << stats.unique_minimizers << "\n";
    stats_out << "lca_assignments=" << stats.lca_assignments << "\n";
    
    stats_out << "\n# Timing Information\n";
    stats_out << "sequence_processing_time=" << std::fixed << std::setprecision(3) << stats.sequence_processing_time << "\n";
    stats_out << "minimizer_extraction_time=" << std::fixed << std::setprecision(3) << stats.minimizer_extraction_time << "\n";
    stats_out << "lca_computation_time=" << std::fixed << std::setprecision(3) << stats.lca_computation_time << "\n";
    stats_out << "database_construction_time=" << std::fixed << std::setprecision(3) << stats.database_construction_time << "\n";
    
    // Calculate derived statistics
    if (stats.total_kmers_processed > 0) {
        double compression_ratio = (double)stats.valid_minimizers_extracted / stats.total_kmers_processed;
        stats_out << "\n# Compression Statistics\n";
        stats_out << "initial_compression_ratio=" << std::fixed << std::setprecision(6) << compression_ratio << "\n";
        
        if (stats.unique_minimizers > 0) {
            double final_compression = (double)stats.unique_minimizers / stats.total_kmers_processed;
            stats_out << "final_compression_ratio=" << std::fixed << std::setprecision(6) << final_compression << "\n";
        }
    }
    
    if (stats.total_bases > 0) {
        double total_time = stats.sequence_processing_time + stats.minimizer_extraction_time + stats.lca_computation_time;
        if (total_time > 0) {
            double processing_rate = stats.total_bases / total_time;
            stats_out << "\n# Performance Statistics\n";
            stats_out << "bases_per_second=" << std::scientific << std::setprecision(2) << processing_rate << "\n";
        }
    }
    
    stats_out.close();
    
    std::cout << "✓ Build statistics saved" << std::endl;
    
    return true;
}

// Helper methods
std::string StandardDatabaseSerializer::get_hash_table_path() const {
    return output_directory_ + "/hash_table.k2d";
}

std::string StandardDatabaseSerializer::get_taxonomy_path() const {
    return output_directory_ + "/taxonomy.tsv";
}

std::string StandardDatabaseSerializer::get_config_path() const {
    return output_directory_ + "/config.txt";
}

std::string StandardDatabaseSerializer::get_stats_path() const {
    return output_directory_ + "/build_stats.txt";
}

bool StandardDatabaseSerializer::create_output_directory() {
    try {
        std::filesystem::create_directories(output_directory_);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create output directory: " << e.what() << std::endl;
        return false;
    }
}

bool StandardDatabaseSerializer::validate_database_consistency() {
    // Check that all required files exist
    std::vector<std::string> required_files = {
        get_hash_table_path(),
        get_taxonomy_path(),
        get_config_path()
    };
    
    for (const auto& file : required_files) {
        if (!std::filesystem::exists(file)) {
            std::cerr << "Missing required database file: " << file << std::endl;
            return false;
        }
    }
    
    // Additional validation can be added here
    return true;
}

bool StandardDatabaseSerializer::generate_database_checksum() {
    // This would generate checksums for all database files
    // Implementation depends on specific checksum requirements
    return true;
}

// ===========================
// EnhancedDatabaseSerializer Implementation
// ===========================

EnhancedDatabaseSerializer::EnhancedDatabaseSerializer(const std::string& output_dir) 
    : output_directory_(output_dir) {
    
    metadata_ = DatabaseMetadata();
    metadata_.version = 2; // Enhanced version
    
    try {
        std::filesystem::create_directories(output_directory_);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not create output directory: " << e.what() << std::endl;
    }
}

bool EnhancedDatabaseSerializer::save_enhanced_database(
    const std::vector<PhylogeneticLCACandidate>& phylo_candidates,
    const ContributingTaxaArrays& contributing_taxa,
    const std::unordered_map<uint32_t, std::string>& taxon_names,
    const EnhancedBuildStats& stats) {
    
    std::cout << "\n=== SAVING ENHANCED DATABASE WITH PHYLOGENETIC DATA ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Update metadata
    metadata_.total_sequences = stats.total_sequences;
    metadata_.total_bases = stats.total_bases;
    metadata_.unique_minimizers = stats.unique_minimizers;
    metadata_.species_count = stats.species_represented;
    
    // Convert phylogenetic candidates to streamlined format
    std::vector<StreamlinedMinimizerMetadata> streamlined_metadata;
    streamlined_metadata.reserve(phylo_candidates.size());
    
    for (const auto& phylo_candidate : phylo_candidates) {
        StreamlinedMinimizerMetadata metadata;
        metadata.minimizer_hash = phylo_candidate.minimizer_hash;
        metadata.lca_taxon = phylo_candidate.lca_taxon;
        metadata.total_genome_count = phylo_candidate.genome_count;
        metadata.phylogenetic_spread = phylo_candidate.phylogenetic_spread;
        metadata.max_phylogenetic_distance = phylo_candidate.max_phylogenetic_distance;
        metadata.num_contributing_taxa = phylo_candidate.contributing_species.size();
        metadata.contributing_taxa_offset = 0; // Will be set when saving arrays
        metadata.reserved = 0;
        
        streamlined_metadata.push_back(metadata);
    }
    
    // Save components
    if (!save_streamlined_hash_table(streamlined_metadata)) {
        std::cerr << "Failed to save streamlined hash table" << std::endl;
        return false;
    }
    
    if (!save_contributing_taxa_arrays(contributing_taxa)) {
        std::cerr << "Failed to save contributing taxa arrays" << std::endl;
        return false;
    }
    
    if (!save_enhanced_config(stats)) {
        std::cerr << "Failed to save enhanced config" << std::endl;
        return false;
    }
    
    if (!save_phylogenetic_summary(phylo_candidates)) {
        std::cerr << "Failed to save phylogenetic summary" << std::endl;
        return false;
    }
    
    // Extract species information for mapping
    std::unordered_map<uint32_t, uint16_t> genome_counts;
    for (const auto& candidate : phylo_candidates) {
        for (size_t i = 0; i < candidate.contributing_species.size(); i++) {
            uint32_t species = candidate.contributing_species[i];
            genome_counts[species] += candidate.genome_counts_per_species[i];
        }
    }
    
    if (!save_species_mapping(taxon_names, genome_counts)) {
        std::cerr << "Failed to save species mapping" << std::endl;
        return false;
    }
    
    // Save backward compatibility layer
    if (!save_standard_compatibility_layer(phylo_candidates)) {
        std::cerr << "Warning: Failed to save standard compatibility layer" << std::endl;
    }
    
    // Validation
    if (!validate_enhanced_format()) {
        std::cerr << "Warning: Enhanced format validation failed" << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "\n=== ENHANCED DATABASE SAVE COMPLETE ===" << std::endl;
    std::cout << "Save time: " << std::fixed << std::setprecision(2) << duration << " seconds" << std::endl;
    std::cout << "Enhanced hash table: " << get_enhanced_hash_table_path() << std::endl;
    std::cout << "Contributing taxa: " << get_contributing_taxa_path() << std::endl;
    std::cout << "Enhanced config: " << get_enhanced_config_path() << std::endl;
    std::cout << "Species mapping: " << get_species_mapping_path() << std::endl;
    
    return true;
}

bool EnhancedDatabaseSerializer::save_streamlined_hash_table(const std::vector<StreamlinedMinimizerMetadata>& metadata) {
    std::string filename = get_enhanced_hash_table_path();
    std::cout << "Writing streamlined hash table to: " << filename << std::endl;
    
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
    
    if (out.fail()) {
        std::cerr << "Failed to write enhanced hash table header" << std::endl;
        return false;
    }
    
    // Write metadata entries
    out.write(reinterpret_cast<const char*>(metadata.data()), 
              metadata.size() * sizeof(StreamlinedMinimizerMetadata));
    
    if (out.fail()) {
        std::cerr << "Failed to write enhanced hash table data" << std::endl;
        return false;
    }
    
    out.close();
    
    std::cout << "✓ Streamlined hash table saved: " << num_entries << " entries (" 
              << (num_entries * sizeof(StreamlinedMinimizerMetadata) / 1024 / 1024) << " MB)" << std::endl;
    
    return true;
}

bool EnhancedDatabaseSerializer::save_contributing_taxa_arrays(const ContributingTaxaArrays& taxa_arrays) {
    std::string filename = get_contributing_taxa_path();
    std::cout << "Writing contributing taxa arrays to: " << filename << std::endl;
    
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Cannot create contributing taxa file: " << filename << std::endl;
        return false;
    }
    
    uint64_t num_entries = taxa_arrays.total_entries();
    
    // Write header
    out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    
    if (out.fail()) {
        std::cerr << "Failed to write contributing taxa header" << std::endl;
        return false;
    }
    
    // Write arrays
    out.write(reinterpret_cast<const char*>(taxa_arrays.taxa_ids.data()), 
              num_entries * sizeof(uint32_t));
    out.write(reinterpret_cast<const char*>(taxa_arrays.phylogenetic_distances.data()), 
              num_entries * sizeof(uint8_t));
    out.write(reinterpret_cast<const char*>(taxa_arrays.genome_counts_per_taxon.data()), 
              num_entries * sizeof(uint16_t));
    
    if (out.fail()) {
        std::cerr << "Failed to write contributing taxa arrays" << std::endl;
        return false;
    }
    
    out.close();
    
    std::cout << "✓ Contributing taxa arrays saved: " << num_entries << " entries" << std::endl;
    
    return true;
}

bool EnhancedDatabaseSerializer::save_enhanced_config(const EnhancedBuildStats& stats) {
    std::string filename = get_enhanced_config_path();
    std::cout << "Writing enhanced config to: " << filename << std::endl;
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot create enhanced config file: " << filename << std::endl;
        return false;
    }
    
    out << "# Enhanced Kraken Database Configuration\n";
    out << "version=2\n";
    out << "k=" << metadata_.k_value << "\n";
    out << "ell=" << metadata_.ell_value << "\n";
    out << "spaces=" << metadata_.spaces_value << "\n";
    out << "subsampling_rate=" << metadata_.subsampling_rate << "\n";
    out << "min_clear_hash_value=0x" << std::hex << metadata_.min_clear_hash_value << std::dec << "\n";
    out << "toggle_mask=0x" << std::hex << metadata_.toggle_mask << std::dec << "\n";
    
    out << "\n# Database Statistics\n";
    out << "total_sequences=" << stats.total_sequences << "\n";
    out << "total_bases=" << stats.total_bases << "\n";
    out << "unique_minimizers=" << stats.unique_minimizers << "\n";
    out << "species_count=" << stats.species_represented << "\n";
    out << "minimizers_with_phylo=" << stats.minimizers_with_phylo_data << "\n";
    out << "contributing_taxa_entries=" << stats.contributing_taxa_array_size << "\n";
    
    out << "\n# Timing Information\n";
    out << "sequence_processing_time=" << stats.sequence_processing_time << "\n";
    out << "minimizer_extraction_time=" << stats.minimizer_extraction_time << "\n";
    out << "phylogenetic_processing_time=" << stats.phylogenetic_processing_time << "\n";
    out << "database_construction_time=" << stats.database_construction_time << "\n";
    
    out.close();
    
    std::cout << "✓ Enhanced config saved" << std::endl;
    
    return true;
}

bool EnhancedDatabaseSerializer::save_species_mapping(
    const std::unordered_map<uint32_t, std::string>& species_names,
    const std::unordered_map<uint32_t, uint16_t>& genome_counts) {
    
    std::string filename = get_species_mapping_path();
    std::cout << "Writing species mapping to: " << filename << std::endl;
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot create species mapping file: " << filename << std::endl;
        return false;
    }
    
    out << "species_taxid\tspecies_name\tgenome_count\n";
    
    for (const auto& [taxid, name] : species_names) {
        uint16_t count = 0;
        auto count_it = genome_counts.find(taxid);
        if (count_it != genome_counts.end()) {
            count = count_it->second;
        }
        
        out << taxid << "\t" << name << "\t" << count << "\n";
    }
    
    out.close();
    
    std::cout << "✓ Species mapping saved: " << species_names.size() << " species" << std::endl;
    
    return true;
}

bool EnhancedDatabaseSerializer::save_standard_compatibility_layer(const std::vector<PhylogeneticLCACandidate>& phylo_candidates) {
    // Convert phylogenetic candidates to standard format for backward compatibility
    std::vector<LCACandidate> standard_candidates;
    standard_candidates.reserve(phylo_candidates.size());
    
    for (const auto& phylo : phylo_candidates) {
        LCACandidate standard;
        standard.minimizer_hash = phylo.minimizer_hash;
        standard.lca_taxon = phylo.lca_taxon;
        standard.genome_count = phylo.genome_count;
        standard.uniqueness_score = phylo.uniqueness_score;
        
        standard_candidates.push_back(standard);
    }
    
    // Create standard serializer for compatibility layer
    StandardDatabaseSerializer compat_serializer(output_directory_);
    
    // Save only the hash table in standard format
    return compat_serializer.save_hash_table(standard_candidates);
}

// Helper methods for enhanced serializer
std::string EnhancedDatabaseSerializer::get_enhanced_hash_table_path() const {
    return output_directory_ + "/enhanced_hash_table.k2d";
}

std::string EnhancedDatabaseSerializer::get_contributing_taxa_path() const {
    return output_directory_ + "/contributing_taxa.k2d";
}

std::string EnhancedDatabaseSerializer::get_enhanced_config_path() const {
    return output_directory_ + "/enhanced_config.txt";
}

std::string EnhancedDatabaseSerializer::get_species_mapping_path() const {
    return output_directory_ + "/species_mapping.tsv";
}

bool EnhancedDatabaseSerializer::validate_enhanced_format() {
    // Validate that all enhanced format files exist and are consistent
    std::vector<std::string> required_files = {
        get_enhanced_hash_table_path(),
        get_contributing_taxa_path(),
        get_enhanced_config_path(),
        get_species_mapping_path()
    };
    
    for (const auto& file : required_files) {
        if (!std::filesystem::exists(file)) {
            std::cerr << "Missing enhanced database file: " << file << std::endl;
            return false;
        }
    }
    
    return true;
}

// Additional placeholder implementations for other classes would go here...
// This provides the core functionality needed for the refactoring
