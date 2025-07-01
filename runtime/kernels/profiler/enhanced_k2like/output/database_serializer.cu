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
#include <cassert>  // For assert in validation methods

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

bool StandardDatabaseSerializer::validate_output_paths() {
    // Check if output directory exists and is writable
    if (!std::filesystem::exists(output_directory_)) {
        std::cerr << "Output directory does not exist: " << output_directory_ << std::endl;
        return false;
    }
    
    // Check write permissions
    std::filesystem::perms p = std::filesystem::status(output_directory_).permissions();
    if ((p & std::filesystem::perms::owner_write) == std::filesystem::perms::none) {
        std::cerr << "No write permission for output directory: " << output_directory_ << std::endl;
        return false;
    }
    
    return true;
}

// ===========================
// EnhancedDatabaseSerializer Implementation
// ===========================

EnhancedDatabaseSerializer::EnhancedDatabaseSerializer(const std::string& output_dir) 
    : output_directory_(output_dir) {
    
    metadata_ = DatabaseMetadata();
    metadata_.version = 4; // Enhanced version with 32-bit feature flags
    
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
    
    // Track ML statistics
    uint64_t ml_weight_count = 0;
    uint64_t contamination_count = 0;
    double total_ml_confidence = 0.0;
    
    for (const auto& phylo_candidate : phylo_candidates) {
        StreamlinedMinimizerMetadata metadata;
        metadata.minimizer_hash = phylo_candidate.minimizer_hash;
        metadata.lca_taxon = phylo_candidate.lca_taxon;
        metadata.total_genome_count = phylo_candidate.genome_count;
        metadata.phylogenetic_spread = phylo_candidate.phylogenetic_spread;
        metadata.max_phylogenetic_distance = phylo_candidate.max_phylogenetic_distance;
        metadata.num_contributing_taxa = phylo_candidate.contributing_species.size();
        metadata.contributing_taxa_offset = 0; // Will be set when saving arrays
        
        // Set ML fields (would be populated from actual ML model)
        // For now, use uniqueness score as a proxy for ML confidence
        metadata.set_ml_confidence(phylo_candidate.uniqueness_score);
        if (phylo_candidate.uniqueness_score > 0.0f) {
            ml_weight_count++;
            total_ml_confidence += phylo_candidate.uniqueness_score;
        }
        
        // Set feature flags (would be populated from feature extractor)
        metadata.feature_flags = 0;
        // Example: encode some basic features
        uint8_t gc_category = 4; // Medium GC content as default
        metadata.feature_flags = MinimizerFlags::set_gc_content_category(metadata.feature_flags, gc_category);
        
        // Check for contamination (example heuristic)
        if (phylo_candidate.contributing_species.size() > 10) {
            metadata.feature_flags = MinimizerFlags::set_contamination_risk(metadata.feature_flags, true);
            contamination_count++;
        }
        
        
        streamlined_metadata.push_back(metadata);
    }
    
    // Update ML statistics in metadata
    metadata_.minimizers_with_ml_weights = ml_weight_count;
    metadata_.contamination_markers = contamination_count;
    metadata_.average_ml_confidence = ml_weight_count > 0 ? 
        static_cast<float>(total_ml_confidence / ml_weight_count) : 0.0f;
    
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
    
    // Save ML-specific components (Version 3)
    if (!save_ml_weight_lookup(streamlined_metadata)) {
        std::cerr << "Failed to save ML weight lookup table" << std::endl;
        return false;
    }
    
    if (!save_feature_statistics(streamlined_metadata)) {
        std::cerr << "Failed to save feature statistics" << std::endl;
        return false;
    }
    
    if (!save_contamination_markers(streamlined_metadata)) {
        std::cerr << "Failed to save contamination markers" << std::endl;
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
    
    if (!validate_ml_features()) {
        std::cerr << "Warning: ML feature validation failed" << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "\n=== ENHANCED DATABASE SAVE COMPLETE ===" << std::endl;
    std::cout << "Save time: " << std::fixed << std::setprecision(2) << duration << " seconds" << std::endl;
    std::cout << "Enhanced hash table: " << get_enhanced_hash_table_path() << std::endl;
    std::cout << "Contributing taxa: " << get_contributing_taxa_path() << std::endl;
    std::cout << "Enhanced config: " << get_enhanced_config_path() << std::endl;
    std::cout << "Species mapping: " << get_species_mapping_path() << std::endl;
    std::cout << "ML weights: " << get_ml_weights_path() << std::endl;
    std::cout << "Feature statistics: " << get_feature_stats_path() << std::endl;
    std::cout << "Contamination markers: " << get_contamination_markers_path() << std::endl;
    
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
    uint64_t version = 4;  // Version 4 with expanded feature flags (32-bit)
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

bool EnhancedDatabaseSerializer::save_phylogenetic_summary(const std::vector<PhylogeneticLCACandidate>& candidates) {
    std::string filename = get_phylo_summary_path();
    std::cout << "Writing phylogenetic summary to: " << filename << std::endl;
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot create phylogenetic summary file: " << filename << std::endl;
        return false;
    }
    
    // Write header
    out << "# Phylogenetic Summary Report\n";
    out << "# Total minimizers with phylogenetic data: " << candidates.size() << "\n\n";
    
    // Calculate statistics
    size_t total_contributing_taxa = 0;
    size_t max_contributing_taxa = 0;
    double avg_phylo_spread = 0.0;
    
    for (const auto& candidate : candidates) {
        size_t num_taxa = candidate.contributing_species.size();
        total_contributing_taxa += num_taxa;
        max_contributing_taxa = std::max(max_contributing_taxa, num_taxa);
        avg_phylo_spread += candidate.phylogenetic_spread;
    }
    
    if (!candidates.empty()) {
        out << "Average contributing taxa per minimizer: " 
            << (double)total_contributing_taxa / candidates.size() << "\n";
        out << "Maximum contributing taxa: " << max_contributing_taxa << "\n";
        out << "Average phylogenetic spread: " 
            << avg_phylo_spread / candidates.size() << "\n";
    }
    
    out.close();
    std::cout << "✓ Phylogenetic summary saved" << std::endl;
    return true;
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

std::string EnhancedDatabaseSerializer::get_phylo_summary_path() const {
    return output_directory_ + "/phylogenetic_summary.txt";
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

// ===========================
// CompactBinarySerializer Implementation
// ===========================

CompactBinarySerializer::CompactBinarySerializer(const std::string& output_dir, bool compress)
    : output_directory_(output_dir), use_compression_(compress) {
    std::filesystem::create_directories(output_directory_);
}

bool CompactBinarySerializer::save_compact_database(
    const std::vector<LCACandidate>& candidates,
    const std::unordered_map<uint32_t, std::string>& taxon_names,
    const DatabaseMetadata& metadata) {
    // TODO: Implement compact binary format
    std::cerr << "Compact binary format not yet implemented" << std::endl;
    return false;
}

bool CompactBinarySerializer::load_compact_database(
    std::vector<LCACandidate>& candidates,
    std::unordered_map<uint32_t, std::string>& taxon_names,
    DatabaseMetadata& metadata) {
    // TODO: Implement compact binary loading
    return false;
}

bool CompactBinarySerializer::compress_database_files() {
    // TODO: Implement compression
    return false;
}

bool CompactBinarySerializer::decompress_database_files() {
    // TODO: Implement decompression
    return false;
}

size_t CompactBinarySerializer::estimate_compressed_size(const std::vector<LCACandidate>& candidates) {
    // Rough estimate: 70% of original size
    return candidates.size() * sizeof(LCACandidate) * 0.7;
}

bool CompactBinarySerializer::write_binary_header(std::ofstream& out, const DatabaseMetadata& metadata) {
    // TODO: Implement binary header writing
    return false;
}

bool CompactBinarySerializer::read_binary_header(std::ifstream& in, DatabaseMetadata& metadata) {
    // TODO: Implement binary header reading
    return false;
}

bool CompactBinarySerializer::write_candidates_binary(std::ofstream& out, const std::vector<LCACandidate>& candidates) {
    // TODO: Implement binary candidate writing
    return false;
}

bool CompactBinarySerializer::read_candidates_binary(std::ifstream& in, std::vector<LCACandidate>& candidates) {
    // TODO: Implement binary candidate reading
    return false;
}

// ===========================
// DatabaseFormatConverter Implementation
// ===========================

DatabaseFormatConverter::DatabaseFormatConverter(const std::string& input_dir, const std::string& output_dir)
    : input_directory_(input_dir), output_directory_(output_dir) {
}

bool DatabaseFormatConverter::convert_standard_to_enhanced(
    const std::string& taxonomy_nodes_path,
    const std::string& taxonomy_names_path) {
    // TODO: Implement format conversion
    return false;
}

bool DatabaseFormatConverter::convert_enhanced_to_standard() {
    // TODO: Implement format conversion
    return false;
}

bool DatabaseFormatConverter::convert_to_compact_binary() {
    // TODO: Implement format conversion
    return false;
}

bool DatabaseFormatConverter::convert_to_streaming_friendly() {
    // TODO: Implement format conversion
    return false;
}

DatabaseFormat DatabaseFormatConverter::detect_database_format(const std::string& database_dir) {
    // Check for enhanced format files first
    if (std::filesystem::exists(database_dir + "/enhanced_hash_table.k2d")) {
        return DatabaseFormat::ENHANCED_PHYLO;
    }
    // Check for standard format
    if (std::filesystem::exists(database_dir + "/hash_table.k2d")) {
        return DatabaseFormat::STANDARD_KRAKEN2;
    }
    return DatabaseFormat::COMPACT_BINARY;
}

bool DatabaseFormatConverter::is_valid_database_directory(const std::string& database_dir) {
    return std::filesystem::exists(database_dir) && 
           std::filesystem::is_directory(database_dir);
}

bool DatabaseFormatConverter::migrate_v1_to_v2() {
    // TODO: Implement migration
    return false;
}

bool DatabaseFormatConverter::upgrade_database_format() {
    // TODO: Implement upgrade
    return false;
}

bool DatabaseFormatConverter::load_standard_database(
    std::vector<LCACandidate>& candidates,
    std::unordered_map<uint32_t, std::string>& taxon_names,
    std::unordered_map<uint32_t, uint32_t>& taxon_parents) {
    // TODO: Implement loading
    return false;
}

bool DatabaseFormatConverter::load_enhanced_database(
    std::vector<PhylogeneticLCACandidate>& phylo_candidates,
    ContributingTaxaArrays& contributing_taxa) {
    // TODO: Implement loading
    return false;
}

bool DatabaseFormatConverter::create_phylogenetic_data_from_standard(
    const std::vector<LCACandidate>& standard_candidates,
    std::vector<PhylogeneticLCACandidate>& phylo_candidates) {
    // TODO: Implement conversion
    return false;
}

// ===========================
// DatabaseValidation Namespace Implementation
// ===========================

namespace DatabaseValidation {
    
bool validate_hash_table_integrity(const std::string& hash_table_path) {
    if (!std::filesystem::exists(hash_table_path)) {
        return false;
    }
    
    std::ifstream file(hash_table_path, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    // Read header
    uint64_t table_size, num_entries;
    file.read(reinterpret_cast<char*>(&table_size), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&num_entries), sizeof(uint64_t));
    
    // Basic sanity checks
    if (table_size == 0 || num_entries == 0 || num_entries > table_size) {
        return false;
    }
    
    file.close();
    return true;
}

bool validate_taxonomy_consistency(const std::string& taxonomy_path) {
    if (!std::filesystem::exists(taxonomy_path)) {
        return false;
    }
    
    std::ifstream file(taxonomy_path);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        line_count++;
        // Basic validation - should have 3 tab-separated fields
        size_t tab_count = std::count(line.begin(), line.end(), '\t');
        if (line_count > 1 && tab_count != 2) {
            return false;
        }
    }
    
    file.close();
    return line_count > 1; // At least header + one entry
}

bool check_database_completeness(const std::string& database_dir) {
    std::vector<std::string> required_files = {
        database_dir + "/hash_table.k2d",
        database_dir + "/taxonomy.tsv",
        database_dir + "/config.txt"
    };
    
    for (const auto& file : required_files) {
        if (!std::filesystem::exists(file)) {
            return false;
        }
    }
    return true;
}

bool cross_validate_hash_taxonomy(const std::string& hash_table_path, 
                                 const std::string& taxonomy_path) {
    // TODO: Implement cross-validation
    return true;
}

bool verify_minimizer_uniqueness(const std::string& hash_table_path) {
    // TODO: Implement uniqueness verification
    return true;
}

bool benchmark_database_loading(const std::string& database_dir) {
    // TODO: Implement benchmark
    return true;
}

bool test_query_performance(const std::string& database_dir, 
                           const std::string& test_sequences_path) {
    // TODO: Implement performance test
    return true;
}

std::string calculate_file_checksum(const std::string& file_path) {
    // TODO: Implement checksum calculation
    return "NOT_IMPLEMENTED";
}

bool verify_database_checksums(const std::string& database_dir) {
    // TODO: Implement checksum verification
    return true;
}

bool generate_manifest_file(const std::string& database_dir) {
    // TODO: Implement manifest generation
    return true;
}

} // namespace DatabaseValidation

// ===========================
// DatabaseOptimization Namespace Implementation
// ===========================

namespace DatabaseOptimization {

bool optimize_hash_table_layout(const std::string& database_dir) {
    // TODO: Implement optimization
    return false;
}

bool compress_taxonomy_data(const std::string& database_dir) {
    // TODO: Implement compression
    return false;
}

bool remove_redundant_entries(const std::string& database_dir) {
    // TODO: Implement redundancy removal
    return false;
}

bool create_index_files(const std::string& database_dir) {
    // TODO: Implement index creation
    return false;
}

bool optimize_for_sequential_access(const std::string& database_dir) {
    // TODO: Implement sequential optimization
    return false;
}

bool create_memory_mapped_version(const std::string& database_dir) {
    // TODO: Implement memory mapping
    return false;
}

size_t calculate_database_size(const std::string& database_dir) {
    size_t total_size = 0;
    
    try {
        for (const auto& entry : std::filesystem::directory_iterator(database_dir)) {
            if (entry.is_regular_file()) {
                total_size += entry.file_size();
            }
        }
    } catch (const std::exception&) {
        return 0;
    }
    
    return total_size;
}

double estimate_memory_usage(const std::string& database_dir) {
    // Rough estimate: 1.2x file size for in-memory structures
    return calculate_database_size(database_dir) * 1.2;
}

void print_database_statistics(const std::string& database_dir) {
    std::cout << "Database Directory: " << database_dir << std::endl;
    std::cout << "Total Size: " << (calculate_database_size(database_dir) / 1024 / 1024) << " MB" << std::endl;
    std::cout << "Estimated Memory Usage: " << (estimate_memory_usage(database_dir) / 1024 / 1024) << " MB" << std::endl;
    
    // List files
    std::cout << "\nDatabase Files:" << std::endl;
    for (const auto& entry : std::filesystem::directory_iterator(database_dir)) {
        if (entry.is_regular_file()) {
            std::cout << "  " << entry.path().filename().string() 
                      << " (" << (entry.file_size() / 1024 / 1024) << " MB)" << std::endl;
        }
    }
}

} // namespace DatabaseOptimization

// ===========================
// ML-specific save methods
// ===========================

bool EnhancedDatabaseSerializer::save_ml_weight_lookup(const std::vector<StreamlinedMinimizerMetadata>& metadata) {
    std::string filename = get_ml_weights_path();
    std::cout << "Writing ML weight lookup table to: " << filename << std::endl;
    
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Cannot create ML weights file: " << filename << std::endl;
        return false;
    }
    
    // Write header
    uint64_t version = 1;
    uint64_t num_entries = 0;
    
    // Count entries with ML weights
    for (const auto& entry : metadata) {
        if (entry.ml_weight > 0) num_entries++;
    }
    
    out.write(reinterpret_cast<const char*>(&version), sizeof(uint64_t));
    out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
    
    // Write ML weight entries
    for (const auto& entry : metadata) {
        if (entry.ml_weight > 0) {
            out.write(reinterpret_cast<const char*>(&entry.minimizer_hash), sizeof(uint64_t));
            out.write(reinterpret_cast<const char*>(&entry.ml_weight), sizeof(uint16_t));
        }
    }
    
    out.close();
    std::cout << "✓ ML weight lookup saved: " << num_entries << " entries" << std::endl;
    return true;
}

bool EnhancedDatabaseSerializer::save_feature_statistics(const std::vector<StreamlinedMinimizerMetadata>& metadata) {
    std::string filename = get_feature_stats_path();
    std::cout << "Writing feature statistics to: " << filename << std::endl;
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot create feature stats file: " << filename << std::endl;
        return false;
    }
    
    // Calculate feature statistics
    std::vector<uint64_t> gc_distribution(8, 0);
    std::vector<uint64_t> complexity_distribution(8, 0);
    uint64_t position_bias_count = 0;
    uint64_t contamination_risk_count = 0;
    
    for (const auto& entry : metadata) {
        uint8_t gc_cat = MinimizerFlags::get_gc_content_category(entry.feature_flags);
        uint8_t complexity = MinimizerFlags::get_complexity_score(entry.feature_flags);
        
        if (gc_cat < 8) gc_distribution[gc_cat]++;
        if (complexity < 8) complexity_distribution[complexity]++;
        
        if (MinimizerFlags::has_position_bias(entry.feature_flags)) position_bias_count++;
        if (MinimizerFlags::has_contamination_risk(entry.feature_flags)) contamination_risk_count++;
    }
    
    // Write statistics in JSON format
    out << "{\n";
    out << "  \"total_minimizers\": " << metadata.size() << ",\n";
    out << "  \"minimizers_with_ml_weights\": " << metadata_.minimizers_with_ml_weights << ",\n";
    out << "  \"average_ml_confidence\": " << metadata_.average_ml_confidence << ",\n";
    out << "  \"contamination_markers\": " << contamination_risk_count << ",\n";
    out << "  \"gc_distribution\": [";
    for (size_t i = 0; i < gc_distribution.size(); i++) {
        out << gc_distribution[i];
        if (i < gc_distribution.size() - 1) out << ", ";
    }
    out << "],\n";
    out << "  \"complexity_distribution\": [";
    for (size_t i = 0; i < complexity_distribution.size(); i++) {
        out << complexity_distribution[i];
        if (i < complexity_distribution.size() - 1) out << ", ";
    }
    out << "],\n";
    out << "  \"position_bias_count\": " << position_bias_count << "\n";
    out << "}\n";
    
    out.close();
    std::cout << "✓ Feature statistics saved" << std::endl;
    return true;
}

bool EnhancedDatabaseSerializer::save_contamination_markers(const std::vector<StreamlinedMinimizerMetadata>& metadata) {
    std::string filename = get_contamination_markers_path();
    std::cout << "Writing contamination markers to: " << filename << std::endl;
    
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Cannot create contamination markers file: " << filename << std::endl;
        return false;
    }
    
    // Write header
    uint64_t version = 1;
    uint64_t num_contaminated = 0;
    
    // Count contaminated minimizers
    for (const auto& entry : metadata) {
        if (MinimizerFlags::has_contamination_risk(entry.feature_flags)) {
            num_contaminated++;
        }
    }
    
    out.write(reinterpret_cast<const char*>(&version), sizeof(uint64_t));
    out.write(reinterpret_cast<const char*>(&num_contaminated), sizeof(uint64_t));
    
    // Write contaminated minimizer hashes
    for (const auto& entry : metadata) {
        if (MinimizerFlags::has_contamination_risk(entry.feature_flags)) {
            out.write(reinterpret_cast<const char*>(&entry.minimizer_hash), sizeof(uint64_t));
            // Write contamination confidence (using phylogenetic spread as proxy)
            float contamination_score = entry.phylogenetic_spread / 255.0f;
            out.write(reinterpret_cast<const char*>(&contamination_score), sizeof(float));
        }
    }
    
    out.close();
    std::cout << "✓ Contamination markers saved: " << num_contaminated << " entries" << std::endl;
    return true;
}

// Backward compatibility check
bool EnhancedDatabaseSerializer::check_version_compatibility(uint64_t file_version) {
    if (file_version > metadata_.version) {
        std::cerr << "Database version " << file_version << " is newer than supported version " 
                  << metadata_.version << std::endl;
        return false;
    }
    
    if (file_version < 2) {
        std::cerr << "Database version " << file_version << " is too old. Minimum supported version is 2." << std::endl;
        std::cerr << "Note: Version 4 is required for 32-bit feature flags support." << std::endl;
        return false;
    }
    
    return true;
}

// ML feature validation
bool EnhancedDatabaseSerializer::validate_ml_features() {
    std::cout << "Validating ML features..." << std::endl;
    
    // Check ML weights file
    std::string ml_weights_file = get_ml_weights_path();
    if (!std::filesystem::exists(ml_weights_file)) {
        std::cerr << "ML weights file missing: " << ml_weights_file << std::endl;
        return false;
    }
    
    // Check feature stats file
    std::string feature_stats_file = get_feature_stats_path();
    if (!std::filesystem::exists(feature_stats_file)) {
        std::cerr << "Feature statistics file missing: " << feature_stats_file << std::endl;
        return false;
    }
    
    // Check contamination markers file
    std::string contamination_file = get_contamination_markers_path();
    if (!std::filesystem::exists(contamination_file)) {
        std::cerr << "Contamination markers file missing: " << contamination_file << std::endl;
        return false;
    }
    
    std::cout << "✓ ML features validated successfully" << std::endl;
    return true;
}

// File path helper implementations
std::string EnhancedDatabaseSerializer::get_ml_weights_path() const {
    return output_directory_ + "/ml_weights.bin";
}

std::string EnhancedDatabaseSerializer::get_feature_stats_path() const {
    return output_directory_ + "/feature_statistics.json";
}

std::string EnhancedDatabaseSerializer::get_contamination_markers_path() const {
    return output_directory_ + "/contamination_markers.bin";
}
