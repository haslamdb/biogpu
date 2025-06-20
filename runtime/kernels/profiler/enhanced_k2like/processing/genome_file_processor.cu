// processing/genome_file_processor.cu
// Implementation of genome file processing functionality
// Extracted from monolithic database builder for better modularity

#include "genome_file_processor.h"
#include "../gpu_kraken_types.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <algorithm>
#include <regex>
#include <sstream>
#include <cctype>
#include <sys/stat.h>
#include <cstdio>

// ===========================
// GenomeFileProcessor Implementation
// ===========================

GenomeFileProcessor::GenomeFileProcessor(const FileProcessingConfig& config) 
    : config_(config) {
    reset_statistics();
}

std::vector<std::string> GenomeFileProcessor::find_genome_files(const std::string& directory) {
    std::vector<std::string> files;
    
    std::cout << "Finding genome files in: " << directory << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Validate directory parameter first
    if (!validate_path_safety(directory)) {
        std::cerr << "Error: Invalid directory path provided" << std::endl;
        return files;
    }
    
    // Use system find command for safety and performance
    std::string find_command = "find \"" + directory + "\" -type f \\( "
                              "-name \"*.fna\" -o -name \"*.fa\" -o -name \"*.fasta\" "
                              "-o -name \"*.ffn\" -o -name \"*.faa\" \\) 2>/dev/null | head -" + 
                              std::to_string(config_.max_file_count);
    
    FILE* pipe = popen(find_command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Failed to execute find command" << std::endl;
        return files;
    }
    
    files.reserve(std::min(config_.max_file_count, size_t(10000)));
    
    char path_buffer[5000];
    while (fgets(path_buffer, sizeof(path_buffer), pipe)) {
        // Remove newline
        size_t len = strlen(path_buffer);
        if (len > 0 && path_buffer[len-1] == '\n') {
            path_buffer[len-1] = '\0';
            len--;
        }
        
        // Validate path length and content
        if (len == 0 || len > 4000) continue;
        
        try {
            std::string file_path(path_buffer);
            if (validate_genome_file(file_path)) {
                files.emplace_back(file_path);
                stats_.files_found++;
                
                // Progress reporting
                if (config_.progress_reporting && stats_.files_found % config_.progress_interval == 0) {
                    std::cout << "Found " << stats_.files_found << " genome files..." << std::endl;
                }
            }
            
            // Safety limit
            if (files.size() >= config_.max_file_count) {
                std::cout << "Reached file limit of " << config_.max_file_count << std::endl;
                break;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error processing file path: " << e.what() << std::endl;
            continue;
        }
    }
    
    pclose(pipe);
    
    // Sort files for consistent processing order
    try {
        std::sort(files.begin(), files.end());
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to sort files: " << e.what() << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time).count();
    stats_.processing_time += duration;
    
    std::cout << "Found " << files.size() << " genome files in " << duration << " seconds" << std::endl;
    
    return files;
}

bool GenomeFileProcessor::validate_genome_file(const std::string& file_path) {
    // Basic path validation
    if (!validate_path_safety(file_path)) {
        return false;
    }
    
    // Check if file exists and is readable
    struct stat file_stat;
    if (stat(file_path.c_str(), &file_stat) != 0) {
        return false;
    }
    
    // Check file size
    if (file_stat.st_size > config_.max_file_size) {
        if (config_.progress_reporting) {
            std::cerr << "Warning: File too large, skipping: " << file_path << std::endl;
        }
        return false;
    }
    
    // Check file extension
    std::filesystem::path p(file_path);
    if (!is_valid_genome_file_extension(p.extension().string())) {
        return false;
    }
    
    return true;
}

std::vector<std::string> GenomeFileProcessor::load_sequences_from_fasta(const std::string& fasta_path) {
    std::vector<std::string> sequences;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Validate input
    if (!validate_genome_file(fasta_path)) {
        std::cerr << "Invalid FASTA file: " << fasta_path << std::endl;
        return sequences;
    }
    
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        std::cerr << "Cannot open FASTA file: " << fasta_path << std::endl;
        stats_.processing_errors++;
        return sequences;
    }
    
    try {
        std::string line, current_sequence;
        bool in_sequence = false;
        int line_number = 0;
        size_t file_sequences = 0;
        size_t file_bases = 0;
        
        current_sequence.reserve(50000); // Reserve space for typical sequence
        
        while (std::getline(file, line)) {
            line_number++;
            
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Process previous sequence
                if (in_sequence && !current_sequence.empty()) {
                    if (current_sequence.size() >= 50 && // Minimum reasonable size
                        current_sequence.size() <= config_.max_sequence_length) {
                        
                        if (!config_.validate_sequences || validate_sequence_content(current_sequence)) {
                            sequences.push_back(current_sequence);
                            file_sequences++;
                            file_bases += current_sequence.size();
                        }
                    }
                    current_sequence.clear();
                    current_sequence.reserve(50000);
                }
                in_sequence = true;
            } else if (in_sequence) {
                // Validate line length
                if (line.size() > 100000) {
                    std::cerr << "Warning: Very long line in " << fasta_path 
                              << " at line " << line_number << std::endl;
                    continue;
                }
                
                // Add line to sequence, removing whitespace
                for (char c : line) {
                    if (!std::isspace(c)) {
                        current_sequence += c;
                    }
                }
                
                // Check sequence length limit
                if (current_sequence.size() > config_.max_sequence_length) {
                    std::cerr << "Warning: Sequence exceeds size limit in " << fasta_path << std::endl;
                    current_sequence.clear();
                    in_sequence = false;
                }
            }
            
            // Progress for very large files
            if (line_number % 100000 == 0 && config_.progress_reporting) {
                std::cout << "Processed " << line_number << " lines from " << fasta_path << std::endl;
            }
        }
        
        // Handle the last sequence
        if (in_sequence && !current_sequence.empty()) {
            if (current_sequence.size() >= 50 && 
                current_sequence.size() <= config_.max_sequence_length) {
                
                if (!config_.validate_sequences || validate_sequence_content(current_sequence)) {
                    sequences.push_back(current_sequence);
                    file_sequences++;
                    file_bases += current_sequence.size();
                }
            }
        }
        
        // Update statistics
        stats_.files_processed++;
        stats_.total_sequences += file_sequences;
        stats_.total_bases += file_bases;
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end_time - start_time).count();
        stats_.processing_time += duration;
        
        if (config_.progress_reporting) {
            std::cout << "Loaded " << sequences.size() << " sequences (" 
                      << (file_bases / 1024 / 1024) << " MB) from " << fasta_path << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error reading FASTA file " << fasta_path << ": " << e.what() << std::endl;
        stats_.processing_errors++;
        sequences.clear();
    }
    
    file.close();
    return sequences;
}

uint32_t GenomeFileProcessor::extract_taxon_from_filename(const std::string& filename) {
    if (!validate_path_safety(filename)) {
        return 1000000; // Default fallback
    }
    
    try {
        std::filesystem::path p(filename);
        std::string stem = p.stem().string();
        
        if (stem.empty()) {
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
            } catch (const std::exception&) {
                // Fall through to other patterns
            }
        }
        
        // Try GCF pattern
        std::regex gcf_pattern(R"(GCF_(\d+\.\d+))");
        if (std::regex_search(stem, match, gcf_pattern)) {
            try {
                std::hash<std::string> hasher;
                uint32_t hash_val = hasher(match[1].str()) % 1000000 + 1000000;
                return hash_val;
            } catch (const std::exception&) {
                // Fall through to default
            }
        }
        
        // Fallback: hash the stem
        std::hash<std::string> hasher;
        uint32_t hash_val = hasher(stem) % 1000000 + 2000000;
        return hash_val;
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing filename " << filename << ": " << e.what() << std::endl;
        return 1000000;
    }
}

bool GenomeFileProcessor::process_genome_files_batch(
    const std::vector<std::string>& file_paths,
    std::vector<std::string>& all_sequences,
    std::vector<uint32_t>& sequence_taxon_ids) {
    
    std::cout << "Processing batch of " << file_paths.size() << " genome files..." << std::endl;
    
    for (const auto& file_path : file_paths) {
        uint32_t taxon_id = extract_taxon_from_filename(file_path);
        std::vector<std::string> file_sequences = load_sequences_from_fasta(file_path);
        
        // Add sequences with their taxon IDs
        for (const auto& sequence : file_sequences) {
            all_sequences.push_back(sequence);
            sequence_taxon_ids.push_back(taxon_id);
        }
        
        if (file_sequences.empty()) {
            stats_.files_skipped++;
        }
    }
    
    return true;
}

void GenomeFileProcessor::print_processing_summary() const {
    std::cout << "\n=== FILE PROCESSING SUMMARY ===" << std::endl;
    std::cout << "Files found: " << stats_.files_found << std::endl;
    std::cout << "Files processed: " << stats_.files_processed << std::endl;
    std::cout << "Files skipped: " << stats_.files_skipped << std::endl;
    std::cout << "Total sequences: " << stats_.total_sequences << std::endl;
    std::cout << "Total bases: " << (stats_.total_bases / 1024 / 1024) << " MB" << std::endl;
    std::cout << "Processing errors: " << stats_.processing_errors << std::endl;
    std::cout << "Processing time: " << std::fixed << std::setprecision(2) 
              << stats_.processing_time << " seconds" << std::endl;
    
    if (stats_.processing_time > 0) {
        double rate = stats_.total_bases / stats_.processing_time / 1024 / 1024;
        std::cout << "Processing rate: " << std::fixed << std::setprecision(2) 
                  << rate << " MB/s" << std::endl;
    }
}

// Private helper methods
bool GenomeFileProcessor::validate_path_safety(const std::string& path) {
    if (path.empty() || path.size() > 4096) return false;
    
    // Check for null characters
    for (char c : path) {
        if (c == '\0') return false;
    }
    
    return true;
}

bool GenomeFileProcessor::is_valid_genome_file_extension(const std::string& extension) {
    return (extension == ".fna" || extension == ".fa" || extension == ".fasta" ||
            extension == ".ffn" || extension == ".faa");
}

bool GenomeFileProcessor::validate_sequence_content(const std::string& sequence) {
    if (sequence.empty()) return false;
    
    // Check for valid DNA characters
    for (char c : sequence) {
        char upper_c = std::toupper(c);
        if (upper_c != 'A' && upper_c != 'C' && upper_c != 'G' && upper_c != 'T' && upper_c != 'N') {
            return false;
        }
    }
    
    return true;
}

void GenomeFileProcessor::reset_statistics() {
    stats_ = FileProcessingStats();
}

// ===========================
// ConcatenatedFnaProcessor Implementation  
// ===========================

ConcatenatedFnaProcessor::ConcatenatedFnaProcessor(const std::string& file_path, 
                                                 const FileProcessingConfig& config)
    : fna_file_path_(file_path), bytes_processed_(0), total_file_size_(0), config_(config) {
    
    // Get file size for progress reporting
    std::ifstream file(file_path, std::ios::ate | std::ios::binary);
    if (file.is_open()) {
        total_file_size_ = file.tellg();
        file.close();
    }
}

bool ConcatenatedFnaProcessor::process_fna_file(
    std::vector<std::string>& genome_files, 
    std::vector<uint32_t>& genome_taxon_ids,
    std::unordered_map<uint32_t, std::string>& taxon_names,
    const std::string& temp_dir) {
    
    std::cout << "Processing concatenated FNA file: " << fna_file_path_ << std::endl;
    std::cout << "File size: " << (total_file_size_ / 1024 / 1024) << " MB" << std::endl;
    
    // Create temporary directory
    try {
        std::filesystem::create_directories(temp_dir);
    } catch (const std::exception& e) {
        std::cerr << "Failed to create temp directory: " << e.what() << std::endl;
        return false;
    }
    
    std::ifstream file(fna_file_path_);
    if (!file.is_open()) {
        std::cerr << "Cannot open FNA file: " << fna_file_path_ << std::endl;
        return false;
    }
    
    std::string line;
    std::string current_sequence;
    std::string current_header;
    uint32_t current_species = 0;
    int genome_count = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    while (std::getline(file, line)) {
        bytes_processed_ += line.length() + 1;
        
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Process previous genome
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
                    if (config_.progress_reporting && genome_count % config_.progress_interval == 0) {
                        double progress = (double)bytes_processed_ / total_file_size_ * 100.0;
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
            current_sequence.reserve(1000000); // Reserve 1MB
            
            if (current_species > 0) {
                // Store species information
                species_data_.add_genome("genome_" + std::to_string(genome_count), 
                                       current_species, parse_result.species_name);
                
                // Add to taxon names
                if (taxon_names.find(current_species) == taxon_names.end()) {
                    taxon_names[current_species] = parse_result.species_name;
                }
            }
            
        } else {
            // Accumulate sequence data
            for (char c : line) {
                if (!std::isspace(c)) {
                    current_sequence += c;
                }
            }
            
            // Check size limit
            if (current_sequence.size() > config_.max_sequence_length) {
                std::cerr << "Warning: Sequence too large, truncating" << std::endl;
                current_sequence.clear();
                current_species = 0;
            }
        }
    }
    
    // Process last genome
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
    std::cout << "  Species represented: " << species_data_.total_species() << std::endl;
    std::cout << "  Processing time: " << duration.count() << " seconds" << std::endl;
    
    if (duration.count() > 0) {
        std::cout << "  Processing rate: " << (bytes_processed_ / 1024 / 1024 / duration.count()) 
                  << " MB/s" << std::endl;
    }
    
    return genome_count > 0;
}

// Header parsing implementation
ConcatenatedFnaProcessor::HeaderParseResult ConcatenatedFnaProcessor::parse_fna_header(const std::string& header) {
    HeaderParseResult result;
    
    // Parse header format: >kraken:taxid|XXXXX|ACCESSION description
    size_t taxid_start = header.find("taxid|");
    if (taxid_start == std::string::npos) {
        return result;
    }
    
    size_t taxid_value_start = taxid_start + 6;
    size_t taxid_end = header.find('|', taxid_value_start);
    if (taxid_end == std::string::npos) {
        return result;
    }
    
    try {
        std::string taxid_str = header.substr(taxid_value_start, taxid_end - taxid_value_start);
        result.species_taxid = std::stoul(taxid_str);
    } catch (const std::exception&) {
        return result;
    }
    
    // Extract accession and description
    size_t accession_start = taxid_end + 1;
    size_t description_start = header.find(' ', accession_start);
    
    if (description_start != std::string::npos) {
        result.accession = header.substr(accession_start, description_start - accession_start);
        result.description = header.substr(description_start + 1);
        result.species_name = extract_species_name_from_description(result.description);
    } else {
        result.accession = header.substr(accession_start);
    }
    
    if (result.species_name.empty()) {
        result.species_name = "species_" + std::to_string(result.species_taxid);
    }
    
    return result;
}

std::string ConcatenatedFnaProcessor::extract_species_name_from_description(const std::string& description) {
    std::istringstream iss(description);
    std::string genus, species;
    
    if (iss >> genus >> species) {
        if (!genus.empty() && !species.empty() && 
            std::isupper(genus[0]) && std::islower(species[0])) {
            return genus + " " + species;
        }
    }
    
    return "";
}

std::string ConcatenatedFnaProcessor::create_temp_genome_file(
    const std::string& sequence, 
    uint32_t species_taxid,
    const std::string& original_header,
    const std::string& temp_dir,
    int genome_index) {
    
    // Skip very short sequences
    if (sequence.length() < 1000) {
        return "";
    }
    
    std::string filename = temp_dir + "/genome_" + std::to_string(species_taxid) + 
                          "_" + std::to_string(genome_index) + ".fasta";
    
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        return "";
    }
    
    // Write header and sequence
    outfile << ">species_" << species_taxid << "_genome_" << genome_index << "\n";
    
    // Write sequence in 80-character lines
    const size_t line_length = 80;
    for (size_t i = 0; i < sequence.length(); i += line_length) {
        size_t end = std::min(i + line_length, sequence.length());
        outfile << sequence.substr(i, end - i) << "\n";
    }
    
    outfile.close();
    return filename;
}

double ConcatenatedFnaProcessor::get_progress_percentage() const {
    if (total_file_size_ == 0) return 0.0;
    return (double)bytes_processed_ / total_file_size_ * 100.0;
}

// ===========================
// Utility Functions Implementation
// ===========================

namespace FileProcessingUtils {
    
    bool is_fasta_file(const std::string& file_path) {
        std::filesystem::path p(file_path);
        std::string ext = p.extension().string();
        return (ext == ".fna" || ext == ".fa" || ext == ".fasta" || 
                ext == ".ffn" || ext == ".faa");
    }
    
    bool validate_dna_sequence(const std::string& sequence) {
        for (char c : sequence) {
            char upper_c = std::toupper(c);
            if (upper_c != 'A' && upper_c != 'C' && upper_c != 'G' && 
                upper_c != 'T' && upper_c != 'N') {
                return false;
            }
        }
        return true;
    }
    
    bool create_safe_directory(const std::string& directory_path) {
        try {
            std::filesystem::create_directories(directory_path);
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }
    
    std::string format_file_size(size_t bytes) {
        const char* units[] = {"B", "KB", "MB", "GB", "TB"};
        int unit = 0;
        double size = bytes;
        
        while (size >= 1024 && unit < 4) {
            size /= 1024;
            unit++;
        }
        
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << size << " " << units[unit];
        return oss.str();
    }
    
    void print_file_processing_progress(size_t current, size_t total, const std::string& current_file) {
        if (total > 0) {
            double percent = (double)current / total * 100.0;
            std::cout << "Processing file " << current << "/" << total 
                      << " (" << std::fixed << std::setprecision(1) << percent << "%): " 
                      << std::filesystem::path(current_file).filename().string() << std::endl;
        }
    }
}
