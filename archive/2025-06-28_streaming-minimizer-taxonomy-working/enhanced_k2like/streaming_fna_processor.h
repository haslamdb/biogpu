// streaming_fna_processor.h
// Streaming processor for large concatenated FNA files
#ifndef STREAMING_FNA_PROCESSOR_H
#define STREAMING_FNA_PROCESSOR_H

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <sstream>

class StreamingFnaProcessor {
private:
    std::ifstream fna_file;
    std::string temp_dir;
    int genomes_per_batch;
    size_t total_genomes_processed;
    size_t total_bases_processed;
    
public:
    StreamingFnaProcessor(const std::string& fna_path, const std::string& temp_directory, int batch_size = 25) 
        : temp_dir(temp_directory), genomes_per_batch(batch_size), 
          total_genomes_processed(0), total_bases_processed(0) {
        fna_file.open(fna_path);
        if (!fna_file.is_open()) {
            throw std::runtime_error("Failed to open FNA file: " + fna_path);
        }
        std::filesystem::create_directories(temp_dir);
        std::cout << "Initialized streaming processor for: " << fna_path << std::endl;
        std::cout << "Temporary directory: " << temp_dir << std::endl;
        std::cout << "Genomes per batch: " << genomes_per_batch << std::endl;
    }
    
    ~StreamingFnaProcessor() {
        if (fna_file.is_open()) {
            fna_file.close();
        }
    }
    
    // Process next batch of genomes, writing temporary files
    bool process_next_batch(std::vector<std::string>& genome_files, 
                           std::vector<uint32_t>& taxon_ids) {
        genome_files.clear();
        taxon_ids.clear();
        
        std::string line, current_sequence, current_header;
        uint32_t current_taxon = 0;
        int genomes_in_batch = 0;
        
        // Continue from where we left off
        while (std::getline(fna_file, line) && genomes_in_batch < genomes_per_batch) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Save previous genome if exists
                if (!current_sequence.empty()) {
                    std::string temp_file = save_temp_genome(current_header, current_sequence, 
                                                            current_taxon, total_genomes_processed);
                    genome_files.push_back(temp_file);
                    taxon_ids.push_back(current_taxon);
                    genomes_in_batch++;
                    total_genomes_processed++;
                    total_bases_processed += current_sequence.length();
                }
                
                // Parse new header
                current_header = line;
                current_taxon = extract_taxon_from_header(line);
                current_sequence.clear();
            } else {
                // Append sequence data (remove whitespace)
                for (char c : line) {
                    if (!std::isspace(c)) {
                        current_sequence += c;
                    }
                }
            }
        }
        
        // Handle last genome in file
        if (!current_sequence.empty() && genomes_in_batch < genomes_per_batch) {
            std::string temp_file = save_temp_genome(current_header, current_sequence, 
                                                    current_taxon, total_genomes_processed);
            genome_files.push_back(temp_file);
            taxon_ids.push_back(current_taxon);
            total_genomes_processed++;
            total_bases_processed += current_sequence.length();
        }
        
        if (!genome_files.empty()) {
            std::cout << "Loaded batch with " << genome_files.size() << " genomes ("
                      << "total processed: " << total_genomes_processed << ")" << std::endl;
        }
        
        return !genome_files.empty();
    }
    
    // Get processing statistics
    size_t get_total_genomes() const { return total_genomes_processed; }
    size_t get_total_bases() const { return total_bases_processed; }
    
private:
    std::string save_temp_genome(const std::string& header, const std::string& sequence, 
                                uint32_t taxon_id, size_t genome_index) {
        std::string filename = temp_dir + "/genome_" + std::to_string(taxon_id) + 
                              "_" + std::to_string(genome_index) + ".fna";
        
        std::ofstream out(filename);
        if (!out.is_open()) {
            throw std::runtime_error("Failed to create temp file: " + filename);
        }
        
        // Write header
        out << header << "\n";
        
        // Write sequence in 80-character lines for standard FASTA format
        for (size_t i = 0; i < sequence.length(); i += 80) {
            out << sequence.substr(i, 80) << "\n";
        }
        out.close();
        
        return filename;
    }
    
    uint32_t extract_taxon_from_header(const std::string& header) {
        // Try multiple formats:
        // Format 1: >kraken:taxid|1234|GCF_000001.1 Organism name
        // Format 2: >taxon_1234|organism_name
        // Format 3: >1234 organism_name
        
        // Try format 1: kraken:taxid|1234|
        size_t taxid_pos = header.find("taxid|");
        if (taxid_pos != std::string::npos) {
            size_t start = taxid_pos + 6;
            size_t end = header.find('|', start);
            if (end != std::string::npos) {
                try {
                    return std::stoul(header.substr(start, end - start));
                } catch (...) {}
            }
        }
        
        // Try format 2: >taxon_1234|
        if (header.find(">taxon_") == 0) {
            size_t start = 7;  // length of ">taxon_"
            size_t end = header.find('|', start);
            if (end == std::string::npos) {
                end = header.find(' ', start);
            }
            if (end != std::string::npos) {
                try {
                    return std::stoul(header.substr(start, end - start));
                } catch (...) {}
            }
        }
        
        // Try format 3: >1234 (taxon ID at start)
        if (header.length() > 1 && std::isdigit(header[1])) {
            size_t start = 1;
            size_t end = start;
            while (end < header.length() && std::isdigit(header[end])) {
                end++;
            }
            if (end > start) {
                try {
                    return std::stoul(header.substr(start, end - start));
                } catch (...) {}
            }
        }
        
        // If no taxon ID found, return 0 (unknown)
        std::cerr << "Warning: Could not extract taxon ID from header: " << header << std::endl;
        return 0;
    }
};

#endif // STREAMING_FNA_PROCESSOR_H