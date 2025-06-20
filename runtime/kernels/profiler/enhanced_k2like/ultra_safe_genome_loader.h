// ultra_safe_genome_loader.h
// EMERGENCY FIX: Pure C-style implementation to avoid ALL STL heap issues

#ifndef ULTRA_SAFE_GENOME_LOADER_H
#define ULTRA_SAFE_GENOME_LOADER_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>

// C-style structures to replace STL containers
struct FileEntry {
    char* path;
    uint32_t taxon_id;
    FileEntry* next;
};

struct FileList {
    FileEntry* head;
    FileEntry* tail;
    size_t count;
    size_t max_path_length;
};

class UltraSafeGenomeLoader {
private:
    FileList* file_list;
    char* temp_buffer;
    size_t buffer_size;
    
    // Memory tracking
    size_t allocations;
    size_t deallocations;
    
public:
    UltraSafeGenomeLoader() : file_list(nullptr), temp_buffer(nullptr), buffer_size(8192), allocations(0), deallocations(0) {
        printf("UltraSafeGenomeLoader: Initializing with pure C memory management\n");
        
        // Allocate file list structure
        file_list = (FileList*)calloc(1, sizeof(FileList));
        if (!file_list) {
            fprintf(stderr, "FATAL: Cannot allocate file list structure\n");
            exit(1);
        }
        allocations++;
        
        file_list->head = nullptr;
        file_list->tail = nullptr;
        file_list->count = 0;
        file_list->max_path_length = 4096;
        
        // Allocate temporary buffer
        temp_buffer = (char*)malloc(buffer_size);
        if (!temp_buffer) {
            fprintf(stderr, "FATAL: Cannot allocate temporary buffer\n");
            exit(1);
        }
        allocations++;
        
        printf("UltraSafeGenomeLoader: Initialized successfully\n");
    }
    
    ~UltraSafeGenomeLoader() {
        printf("UltraSafeGenomeLoader: Cleaning up\n");
        cleanup();
        printf("UltraSafeGenomeLoader: Allocations=%zu, Deallocations=%zu\n", allocations, deallocations);
    }
    
    void cleanup() {
        if (file_list) {
            clear_file_list();
            free(file_list);
            file_list = nullptr;
            deallocations++;
        }
        
        if (temp_buffer) {
            free(temp_buffer);
            temp_buffer = nullptr;
            deallocations++;
        }
    }
    
    void clear_file_list() {
        if (!file_list) return;
        
        FileEntry* current = file_list->head;
        while (current) {
            FileEntry* next = current->next;
            
            if (current->path) {
                free(current->path);
                deallocations++;
            }
            
            free(current);
            deallocations++;
            
            current = next;
        }
        
        file_list->head = nullptr;
        file_list->tail = nullptr;
        file_list->count = 0;
    }
    
    bool validate_path(const char* path) {
        if (!path) return false;
        
        size_t len = strlen(path);
        if (len == 0 || len > file_list->max_path_length) {
            return false;
        }
        
        // Check for null bytes in middle of string
        for (size_t i = 0; i < len; i++) {
            if (path[i] == '\0') return false;
        }
        
        return true;
    }
    
    bool add_file_entry(const char* path, uint32_t taxon_id) {
        if (!validate_path(path)) {
            fprintf(stderr, "Invalid path rejected: %s\n", path ? path : "NULL");
            return false;
        }
        
        // Allocate new entry
        FileEntry* entry = (FileEntry*)calloc(1, sizeof(FileEntry));
        if (!entry) {
            fprintf(stderr, "Failed to allocate file entry\n");
            return false;
        }
        allocations++;
        
        // Copy path
        size_t path_len = strlen(path);
        entry->path = (char*)malloc(path_len + 1);
        if (!entry->path) {
            fprintf(stderr, "Failed to allocate path string\n");
            free(entry);
            deallocations++;
            return false;
        }
        allocations++;
        
        strcpy(entry->path, path);
        entry->taxon_id = taxon_id;
        entry->next = nullptr;
        
        // Add to list
        if (!file_list->head) {
            file_list->head = entry;
            file_list->tail = entry;
        } else {
            file_list->tail->next = entry;
            file_list->tail = entry;
        }
        
        file_list->count++;
        return true;
    }
    
    uint32_t extract_taxon_from_path(const char* path) {
        if (!path) return 1000000;
        
        // Simple hash-based extraction
        uint32_t hash = 5381;
        const char* p = path;
        
        while (*p) {
            hash = ((hash << 5) + hash) + (unsigned char)*p;
            p++;
        }
        
        // Return hash in safe range
        return (hash % 900000) + 100000;
    }
    
    bool find_genome_files_ultra_safe(const char* directory) {
        printf("UltraSafeGenomeLoader: Scanning directory %s\n", directory);
        
        if (!directory || strlen(directory) == 0) {
            fprintf(stderr, "Invalid directory path\n");
            return false;
        }
        
        // Check if directory exists using stat
        struct stat dir_stat;
        if (stat(directory, &dir_stat) != 0) {
            fprintf(stderr, "Directory does not exist: %s\n", directory);
            return false;
        }
        
        if (!S_ISDIR(dir_stat.st_mode)) {
            fprintf(stderr, "Path is not a directory: %s\n", directory);
            return false;
        }
        
        // Use find command to get file list
        char find_cmd[8192];
        int ret = snprintf(find_cmd, sizeof(find_cmd),
            "find \"%s\" -type f \\( -name \"*.fna\" -o -name \"*.fa\" -o -name \"*.fasta\" \\) 2>/dev/null | head -10000",
            directory);
        
        if (ret >= sizeof(find_cmd) || ret < 0) {
            fprintf(stderr, "Find command too long\n");
            return false;
        }
        
        printf("Executing: %s\n", find_cmd);
        
        FILE* pipe = popen(find_cmd, "r");
        if (!pipe) {
            fprintf(stderr, "Failed to execute find command\n");
            return false;
        }
        
        size_t files_processed = 0;
        while (fgets(temp_buffer, buffer_size, pipe) && files_processed < 10000) {
            // Remove newline
            size_t len = strlen(temp_buffer);
            if (len > 0 && temp_buffer[len-1] == '\n') {
                temp_buffer[len-1] = '\0';
                len--;
            }
            
            if (len == 0) continue;
            
            // Validate path
            if (!validate_path(temp_buffer)) {
                fprintf(stderr, "Skipping invalid path: %s\n", temp_buffer);
                continue;
            }
            
            // Extract taxon ID
            uint32_t taxon_id = extract_taxon_from_path(temp_buffer);
            
            // Add to list
            if (add_file_entry(temp_buffer, taxon_id)) {
                files_processed++;
                
                if (files_processed % 100 == 0) {
                    printf("Processed %zu files...\n", files_processed);
                }
            } else {
                fprintf(stderr, "Failed to add file: %s\n", temp_buffer);
            }
        }
        
        pclose(pipe);
        
        printf("UltraSafeGenomeLoader: Found %zu genome files\n", file_list->count);
        return file_list->count > 0;
    }
    
    // Convert to STL containers ONLY when needed
    bool export_to_vectors(std::vector<std::string>& genome_files, 
                          std::vector<uint32_t>& genome_taxon_ids) {
        printf("UltraSafeGenomeLoader: Converting to STL vectors\n");
        
        // Clear output vectors
        genome_files.clear();
        genome_taxon_ids.clear();
        
        if (!file_list || file_list->count == 0) {
            printf("No files to export\n");
            return false;
        }
        
        try {
            // Reserve space
            genome_files.reserve(file_list->count);
            genome_taxon_ids.reserve(file_list->count);
            
            // Convert entries
            FileEntry* current = file_list->head;
            size_t exported = 0;
            
            while (current && exported < file_list->count) {
                if (current->path && validate_path(current->path)) {
                    genome_files.emplace_back(current->path);
                    genome_taxon_ids.push_back(current->taxon_id);
                    exported++;
                }
                current = current->next;
            }
            
            printf("Exported %zu files to STL vectors\n", exported);
            return exported > 0;
            
        } catch (const std::exception& e) {
            fprintf(stderr, "STL export failed: %s\n", e.what());
            return false;
        }
    }
    
    size_t get_file_count() const {
        return file_list ? file_list->count : 0;
    }
    
    void print_statistics() {
        printf("\n=== ULTRA SAFE LOADER STATISTICS ===\n");
        printf("Files found: %zu\n", get_file_count());
        printf("Memory allocations: %zu\n", allocations);
        printf("Memory deallocations: %zu\n", deallocations);
        printf("Memory balance: %s\n", (allocations == deallocations) ? "CLEAN" : "LEAK DETECTED");
    }
};

// EMERGENCY REPLACEMENT FOR load_genome_files
bool load_genome_files_emergency_fix(const std::string& library_path,
                                    std::vector<std::string>& genome_files,
                                    std::vector<uint32_t>& genome_taxon_ids,
                                    std::unordered_map<uint32_t, std::string>& taxon_names) {
    
    printf("\n=== EMERGENCY ULTRA-SAFE GENOME LOADING ===\n");
    
    // Use C-style loader to avoid STL heap issues
    UltraSafeGenomeLoader loader;
    
    // Find files using pure C memory management
    if (!loader.find_genome_files_ultra_safe(library_path.c_str())) {
        fprintf(stderr, "Failed to find genome files\n");
        return false;
    }
    
    loader.print_statistics();
    
    // Export to STL containers only at the end
    bool success = loader.export_to_vectors(genome_files, genome_taxon_ids);
    
    if (success) {
        // Update taxon names
        for (size_t i = 0; i < genome_files.size(); i++) {
            uint32_t taxon_id = genome_taxon_ids[i];
            
            if (taxon_names.find(taxon_id) == taxon_names.end()) {
                // Extract simple name
                std::string path = genome_files[i];
                size_t last_slash = path.find_last_of("/\\");
                std::string filename = (last_slash != std::string::npos) ? 
                                      path.substr(last_slash + 1) : path;
                
                size_t last_dot = filename.find_last_of('.');
                std::string stem = (last_dot != std::string::npos) ? 
                                  filename.substr(0, last_dot) : filename;
                
                taxon_names[taxon_id] = stem.empty() ? ("taxon_" + std::to_string(taxon_id)) : stem;
            }
        }
        
        printf("Successfully loaded %zu genome files\n", genome_files.size());
    }
    
    return success;
}

#endif // ULTRA_SAFE_GENOME_LOADER_H
