// processing/genome_file_processor.h
// Genome file discovery, validation, and sequence loading
// Handles multiple input formats and provides clean error handling

#ifndef GENOME_FILE_PROCESSOR_H
#define GENOME_FILE_PROCESSOR_H

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <fstream>
#include "../gpu_kraken_types.h"

// Note: FileProcessingConfig and FileProcessingStats are defined in gpu_kraken_types.h

// Individual genome file processor
class GenomeFileProcessor {
private:
    FileProcessingConfig config_;
    FileProcessingStats stats_;
    
public:
    explicit GenomeFileProcessor(const FileProcessingConfig& config = FileProcessingConfig());
    
    // File discovery methods
    std::vector<std::string> find_genome_files(const std::string& directory);
    std::vector<std::string> load_file_list(const std::string& file_list_path);
    bool validate_genome_file(const std::string& file_path);
    
    // Sequence loading methods
    std::vector<std::string> load_sequences_from_fasta(const std::string& fasta_path);
    bool count_sequences_in_file(const std::string& file_path, size_t& sequence_count, size_t& total_bases);
    
    // Taxon extraction and metadata
    uint32_t extract_taxon_from_filename(const std::string& filename);
    std::string extract_species_name_from_file(const std::string& file_path);
    
    // Batch processing
    bool process_genome_files_batch(
        const std::vector<std::string>& file_paths,
        std::vector<std::string>& all_sequences,
        std::vector<uint32_t>& sequence_taxon_ids
    );
    
    // Statistics and monitoring
    const FileProcessingStats& get_statistics() const { return stats_; }
    void print_processing_summary() const;
    void reset_statistics();
    
private:
    // Internal validation methods
    bool validate_path_safety(const std::string& path);
    bool is_valid_genome_file_extension(const std::string& extension);
    bool validate_fasta_format(const std::string& file_path);
    bool validate_sequence_content(const std::string& sequence);
    
    // Internal processing helpers
    bool process_single_fasta_file(
        const std::string& file_path,
        std::vector<std::string>& sequences,
        std::vector<uint32_t>& taxon_ids,
        uint32_t default_taxon_id
    );
    
    void update_processing_stats(const std::string& file_path, size_t sequences_added, size_t bases_added);
};

// Concatenated FNA file processor
class ConcatenatedFnaProcessor {
private:
    std::string fna_file_path_;
    SpeciesTrackingData species_data_;
    size_t bytes_processed_;
    size_t total_file_size_;
    FileProcessingConfig config_;
    
public:
    explicit ConcatenatedFnaProcessor(const std::string& file_path, 
                                    const FileProcessingConfig& config = FileProcessingConfig());
    
    // Main processing method
    bool process_fna_file(
        std::vector<std::string>& genome_files, 
        std::vector<uint32_t>& genome_taxon_ids,
        std::unordered_map<uint32_t, std::string>& taxon_names,
        const std::string& temp_dir = "/tmp/kraken_build"
    );
    
    // Access to species data
    const SpeciesTrackingData& get_species_data() const { return species_data_; }
    
    // Progress monitoring
    double get_progress_percentage() const;
    size_t get_bytes_processed() const { return bytes_processed_; }
    size_t get_total_file_size() const { return total_file_size_; }
    
private:
    // Header parsing structures and methods
    struct HeaderParseResult {
        uint32_t species_taxid = 0;
        std::string species_name;
        std::string accession;
        std::string description;
    };
    
    HeaderParseResult parse_fna_header(const std::string& header);
    std::string extract_species_name_from_description(const std::string& description);
    
    // Temporary file creation
    std::string create_temp_genome_file(
        const std::string& sequence, 
        uint32_t species_taxid,
        const std::string& original_header,
        const std::string& temp_dir,
        int genome_index
    );
    
    // Validation helpers
    bool validate_sequence_content(const std::string& sequence);
    bool is_reasonable_sequence_length(size_t length);
};

// Streaming FNA processor for very large files
class StreamingFnaProcessor {
private:
    std::string fna_file_path_;
    std::string temp_directory_;
    size_t batch_size_;
    size_t current_genome_count_;
    size_t total_bases_processed_;
    bool processing_active_;
    
    // Enhanced members for better streaming
    std::ifstream file_stream_;
    std::string current_buffer_;
    bool end_of_file_reached_;
    
    // Chunking configuration
    static constexpr size_t BUFFER_SIZE = 64 * 1024 * 1024;  // 64MB read buffer
    static constexpr size_t MAX_SEQUENCE_SIZE = 50 * 1024 * 1024;  // 50MB max per genome
    
    // State management for incomplete sequences
    struct IncompleteSequence {
        std::string header;
        std::string sequence;
        uint32_t taxon_id;
        
        void clear() {
            header.clear();
            sequence.clear();
            taxon_id = 0;
        }
        
        bool empty() const {
            return header.empty() && sequence.empty();
        }
    };
    
    IncompleteSequence incomplete_state_;
    std::vector<char> read_buffer_;
    
    // Statistics
    size_t total_genomes_read_;
    size_t sequences_too_large_;
    size_t sequences_with_invalid_taxon_;
    
public:
    StreamingFnaProcessor(const std::string& fna_path, 
                         const std::string& temp_dir,
                         size_t batch_size = 25);
    ~StreamingFnaProcessor();
    
    // Streaming processing methods
    bool process_next_batch(std::vector<std::string>& batch_files, 
                           std::vector<uint32_t>& batch_taxons);
    bool has_more_data() const;
    void reset_processing();
    
    // Statistics
    size_t get_total_genomes() const { return current_genome_count_; }
    size_t get_total_bases() const { return total_bases_processed_; }
    size_t get_sequences_too_large() const { return sequences_too_large_; }
    size_t get_sequences_with_invalid_taxon() const { return sequences_with_invalid_taxon_; }
    
private:
    // Enhanced internal methods
    bool process_batch_from_buffer(std::vector<std::string>& batch_files, 
                                  std::vector<uint32_t>& batch_taxons);
    bool read_next_chunk_to_buffer();
    void cleanup_temp_files();
    uint32_t parse_taxon_from_header(const std::string& header);
    bool write_genome_to_temp_file(const std::string& header,
                                  const std::string& sequence,
                                  uint32_t taxon_id,
                                  std::string& temp_file_path);
};

// Utility functions for file processing
namespace FileProcessingUtils {
    // File validation utilities
    bool is_fasta_file(const std::string& file_path);
    bool is_compressed_file(const std::string& file_path);
    size_t estimate_uncompressed_size(const std::string& compressed_file);
    
    // Sequence validation utilities
    bool validate_dna_sequence(const std::string& sequence);
    bool has_valid_dna_characters(const std::string& sequence);
    double calculate_gc_content(const std::string& sequence);
    
    // File system utilities
    bool create_safe_directory(const std::string& directory_path);
    bool cleanup_directory(const std::string& directory_path);
    std::string generate_temp_filename(const std::string& base_dir, const std::string& prefix, int index);
    
    // Progress and monitoring utilities
    void print_file_processing_progress(size_t current, size_t total, const std::string& current_file);
    std::string format_file_size(size_t bytes);
    std::string format_processing_rate(size_t bytes, double seconds);
}

#endif // GENOME_FILE_PROCESSOR_H
