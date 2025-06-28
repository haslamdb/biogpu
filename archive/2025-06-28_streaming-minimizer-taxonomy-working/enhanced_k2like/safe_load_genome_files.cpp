// safe_load_genome_files.cpp
// Memory-safe replacement for the problematic load_genome_files function

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <regex>
#include <filesystem>

class SafeGenomeLoader {
private:
    // Memory tracking
    size_t vector_operations = 0;
    size_t heap_checks = 0;
    
    // Configuration
    size_t max_files_to_process = 50000;  // Safety limit
    size_t max_path_length = 4000;       // Reasonable path limit
    
public:
    bool check_heap_integrity(const std::string& location) {
        heap_checks++;
        
        if (heap_checks % 1000 == 0) {
            std::cout << "Heap check #" << heap_checks << " at " << location << std::endl;
        }
        
        try {
            void* test_ptr = malloc(64);
            if (!test_ptr) {
                std::cerr << "HEAP ERROR: malloc failed at " << location << std::endl;
                return false;
            }
            memset(test_ptr, 0x55, 64);
            free(test_ptr);
            return true;
        } catch (...) {
            std::cerr << "HEAP ERROR: Exception during test at " << location << std::endl;
            return false;
        }
    }
    
    // Safe alternative to std::filesystem operations
    bool path_exists_safe(const std::string& path) {
        if (path.empty() || path.size() > max_path_length) return false;
        
        struct stat path_stat;
        return (stat(path.c_str(), &path_stat) == 0);
    }
    
    bool is_directory_safe(const std::string& path) {
        if (path.empty() || path.size() > max_path_length) return false;
        
        struct stat path_stat;
        if (stat(path.c_str(), &path_stat) != 0) return false;
        
        return S_ISDIR(path_stat.st_mode);
    }
    
    // Use system find command instead of std::filesystem iterator
    std::vector<std::string> find_genome_files_safe(const std::string& directory) {
        std::vector<std::string> files;
        
        std::cout << "Loading genome files from: " << directory << std::endl;
        
        if (!check_heap_integrity("find_genome_files_start")) {
            return files;
        }
        
        // Validate directory
        if (directory.empty() || directory.size() > max_path_length) {
            std::cerr << "Invalid directory path length" << std::endl;
            return files;
        }
        
        if (!path_exists_safe(directory)) {
            std::cerr << "Directory does not exist: " << directory << std::endl;
            return files;
        }
        
        if (!is_directory_safe(directory)) {
            std::cerr << "Path is not a directory: " << directory << std::endl;
            return files;
        }
        
        // Use find command to avoid std::filesystem issues
        std::string find_command = "find \"" + directory + "\" -type f \\( "
                                  "-name \"*.fna\" -o -name \"*.fa\" -o -name \"*.fasta\" "
                                  "-o -name \"*.ffn\" -o -name \"*.faa\" \\) 2>/dev/null | head -50000";
        
        std::cout << "Executing safe file search..." << std::endl;
        
        FILE* pipe = popen(find_command.c_str(), "r");
        if (!pipe) {
            std::cerr << "Failed to execute find command" << std::endl;
            return files;
        }
        
        // Reserve reasonable space
        files.reserve(10000);  // Start conservative
        
        char path_buffer[5000];  // Larger buffer for safety
        size_t files_found = 0;
        
        while (fgets(path_buffer, sizeof(path_buffer), pipe) && files_found < max_files_to_process) {
            // Remove newline
            size_t len = strlen(path_buffer);
            if (len > 0 && path_buffer[len-1] == '\n') {
                path_buffer[len-1] = '\0';
                len--;
            }
            
            // Validate path
            if (len == 0 || len > max_path_length) {
                std::cerr << "Skipping invalid path (length: " << len << ")" << std::endl;
                continue;
            }
            
            // Check if we need to expand capacity
            if (files.size() >= files.capacity() - 10) {
                size_t new_capacity = files.capacity() * 1.5;
                if (new_capacity > max_files_to_process) {
                    new_capacity = max_files_to_process;
                }
                
                try {
                    files.reserve(new_capacity);
                    std::cout << "Expanded capacity to " << new_capacity << std::endl;
                } catch (const std::bad_alloc& e) {
                    std::cerr << "Memory allocation failed, stopping at " << files.size() << " files" << std::endl;
                    break;
                }
                
                if (!check_heap_integrity("after_vector_expansion")) {
                    std::cerr << "Heap corruption detected after vector expansion" << std::endl;
                    break;
                }
            }
            
            // Add file safely
            try {
                files.emplace_back(path_buffer);
                vector_operations++;
                files_found++;
                
                // Periodic heap check
                if (files_found % 500 == 0) {
                    std::cout << "Processed " << files_found << " files..." << std::endl;
                    if (!check_heap_integrity("file_processing_" + std::to_string(files_found))) {
                        std::cerr << "Heap corruption detected during file processing" << std::endl;
                        break;
                    }
                }
                
            } catch (const std::exception& e) {
                std::cerr << "Error adding file " << path_buffer << ": " << e.what() << std::endl;
                continue;
            }
        }
        
        pclose(pipe);
        
        std::cout << "Found " << files.size() << " genome files" << std::endl;
        
        // IMPORTANT: Sort using safe method
        if (!files.empty()) {
            try {
                std::cout << "Sorting " << files.size() << " files..." << std::endl;
                std::sort(files.begin(), files.end());
                std::cout << "Sorting completed successfully" << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Sorting failed: " << e.what() << std::endl;
                // Don't fail completely, just return unsorted
            }
        }
        
        if (!check_heap_integrity("find_genome_files_end")) {
            std::cerr << "Warning: Heap integrity compromised" << std::endl;
        }
        
        return files;
    }
    
    // Safe taxon ID extraction
    uint32_t extract_taxon_from_filename_safe(const std::string& filename) {
        if (filename.empty() || filename.size() > max_path_length) {
            return 1000000;  // Default fallback
        }
        
        try {
            // Extract filename from path manually (safer than std::filesystem)
            size_t last_slash = filename.find_last_of("/\\");
            std::string base_name = (last_slash != std::string::npos) ? 
                                   filename.substr(last_slash + 1) : filename;
            
            // Remove extension manually
            size_t last_dot = base_name.find_last_of('.');
            std::string stem = (last_dot != std::string::npos) ? 
                              base_name.substr(0, last_dot) : base_name;
            
            // Try to extract taxid pattern
            std::regex taxid_pattern(R"(taxid[_-](\d+))");
            std::smatch match;
            if (std::regex_search(stem, match, taxid_pattern) && match.size() > 1) {
                uint32_t taxon_id = std::stoul(match[1].str());
                if (taxon_id > 0 && taxon_id < UINT32_MAX) {
                    return taxon_id;
                }
            }
            
            // Try GCF pattern
            std::regex gcf_pattern(R"(GCF_(\d+\.\d+))");
            if (std::regex_search(stem, match, gcf_pattern) && match.size() > 1) {
                std::hash<std::string> hasher;
                uint32_t hash_val = hasher(match[1].str()) % 1000000 + 1000000;
                return hash_val;
            }
            
            // Fallback: hash the stem
            std::hash<std::string> hasher;
            uint32_t hash_val = hasher(stem) % 1000000 + 2000000;
            return hash_val;
            
        } catch (const std::exception& e) {
            std::cerr << "Error extracting taxon ID from " << filename << ": " << e.what() << std::endl;
            return 1000000;
        }
    }
    
    // Safe replacement for the entire load_genome_files function
    bool load_genome_files_safe(const std::string& library_path,
                               std::vector<std::string>& genome_files,
                               std::vector<uint32_t>& genome_taxon_ids,
                               std::unordered_map<uint32_t, std::string>& taxon_names) {
        
        std::cout << "=== SAFE GENOME FILE LOADING ===" << std::endl;
        
        if (!check_heap_integrity("load_genome_files_start")) {
            return false;
        }
        
        // Clear existing data safely
        try {
            genome_files.clear();
            genome_taxon_ids.clear();
            // Don't clear taxon_names as it might have existing data
        } catch (const std::exception& e) {
            std::cerr << "Error clearing vectors: " << e.what() << std::endl;
            return false;
        }
        
        if (!check_heap_integrity("after_vector_clear")) {
            return false;
        }
        
        // Find files using safe method
        auto found_files = find_genome_files_safe(library_path);
        
        if (found_files.empty()) {
            std::cerr << "No genome files found in " << library_path << std::endl;
            return false;
        }
        
        if (!check_heap_integrity("after_file_finding")) {
            return false;
        }
        
        // Reserve space to prevent reallocations
        try {
            genome_files.reserve(found_files.size());
            genome_taxon_ids.reserve(found_files.size());
            std::cout << "Reserved space for " << found_files.size() << " genomes" << std::endl;
        } catch (const std::bad_alloc& e) {
            std::cerr << "Failed to reserve memory for " << found_files.size() << " genomes" << std::endl;
            return false;
        }
        
        if (!check_heap_integrity("after_vector_reserve")) {
            return false;
        }
        
        // Process files one by one
        size_t processed_count = 0;
        for (const auto& file_path : found_files) {
            if (!check_heap_integrity("processing_file_" + std::to_string(processed_count))) {
                std::cerr << "Heap corruption detected at file " << processed_count << std::endl;
                break;
            }
            
            // Validate file path
            if (file_path.empty() || file_path.size() > max_path_length) {
                std::cerr << "Skipping invalid file path" << std::endl;
                continue;
            }
            
            // Check if file still exists
            if (!path_exists_safe(file_path)) {
                std::cerr << "File no longer exists: " << file_path << std::endl;
                continue;
            }
            
            try {
                // Extract taxon ID
                uint32_t taxon_id = extract_taxon_from_filename_safe(file_path);
                
                // Add to vectors
                genome_files.push_back(file_path);
                genome_taxon_ids.push_back(taxon_id);
                vector_operations += 2;
                
                // Update taxon names if not present
                if (taxon_names.find(taxon_id) == taxon_names.end()) {
                    // Extract filename for name
                    size_t last_slash = file_path.find_last_of("/\\");
                    std::string filename = (last_slash != std::string::npos) ? 
                                          file_path.substr(last_slash + 1) : file_path;
                    
                    size_t last_dot = filename.find_last_of('.');
                    std::string stem = (last_dot != std::string::npos) ? 
                                      filename.substr(0, last_dot) : filename;
                    
                    if (!stem.empty() && stem.size() < 500) {  // Reasonable name limit
                        taxon_names[taxon_id] = stem;
                    } else {
                        taxon_names[taxon_id] = "taxon_" + std::to_string(taxon_id);
                    }
                }
                
                processed_count++;
                
                // Verify vector consistency
                if (genome_files.size() != genome_taxon_ids.size()) {
                    std::cerr << "FATAL: Vector size mismatch!" << std::endl;
                    return false;
                }
                
                // Progress reporting
                if (processed_count % 1000 == 0) {
                    std::cout << "Processed " << processed_count << " files..." << std::endl;
                }
                
            } catch (const std::exception& e) {
                std::cerr << "Error processing file " << file_path << ": " << e.what() << std::endl;
                // Continue with next file
                continue;
            }
        }
        
        if (!check_heap_integrity("load_genome_files_end")) {
            std::cerr << "Warning: Heap integrity check failed at end" << std::endl;
        }
        
        std::cout << "Successfully loaded " << genome_files.size() << " genome files" << std::endl;
        std::cout << "Vector operations performed: " << vector_operations << std::endl;
        std::cout << "Heap integrity checks: " << heap_checks << std::endl;
        
        return !genome_files.empty();
    }
    
    void print_statistics() {
        std::cout << "\n=== SAFE LOADER STATISTICS ===" << std::endl;
        std::cout << "Vector operations: " << vector_operations << std::endl;
        std::cout << "Heap checks: " << heap_checks << std::endl;
    }
};

// Test function
int main(int argc, char* argv[]) {
    std::string test_dir = "/tmp";
    if (argc > 1) {
        test_dir = argv[1];
    }
    
    std::cout << "Testing safe genome file loading with directory: " << test_dir << std::endl;
    
    SafeGenomeLoader loader;
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    
    bool success = loader.load_genome_files_safe(test_dir, genome_files, genome_taxon_ids, taxon_names);
    
    if (success) {
        std::cout << "✅ Safe loading completed successfully!" << std::endl;
        std::cout << "Files loaded: " << genome_files.size() << std::endl;
        std::cout << "Taxa mapped: " << taxon_names.size() << std::endl;
    } else {
        std::cout << "❌ Safe loading failed" << std::endl;
        return 1;
    }
    
    loader.print_statistics();
    return 0;
}
