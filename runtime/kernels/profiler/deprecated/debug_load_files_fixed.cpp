// debug_load_files_fixed.cpp
// FIXED: Added missing headers and improved memory safety checks

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sys/stat.h>
#include <cstdlib>
#include <cstring>
#include <algorithm>  // ADDED: For std::sort
#include <exception>
#include <stdexcept>

// Memory debugging utilities
class SafeMemoryTracker {
private:
    static size_t allocation_count;
    static size_t deallocation_count;
    
public:
    static void track_allocation() { allocation_count++; }
    static void track_deallocation() { deallocation_count++; }
    
    static bool check_heap_integrity(const std::string& location) {
        std::cout << "HEAP CHECK: " << location 
                  << " (allocs: " << allocation_count 
                  << ", deallocs: " << deallocation_count << ")" << std::endl;
        
        // Try a small allocation/deallocation
        try {
            void* test_ptr = malloc(16);
            if (!test_ptr) {
                std::cerr << "ERROR: malloc(16) failed at " << location << std::endl;
                return false;
            }
            memset(test_ptr, 0xAA, 16);
            free(test_ptr);
            track_allocation();
            track_deallocation();
            return true;
        } catch (...) {
            std::cerr << "ERROR: Exception during heap test at " << location << std::endl;
            return false;
        }
    }
    
    static void print_stats() {
        std::cout << "Memory stats: " << allocation_count << " allocations, " 
                  << deallocation_count << " deallocations" << std::endl;
        if (allocation_count != deallocation_count) {
            std::cout << "WARNING: Memory leak detected!" << std::endl;
        }
    }
};

size_t SafeMemoryTracker::allocation_count = 0;
size_t SafeMemoryTracker::deallocation_count = 0;

// Safe filesystem operations
class SafeFilesystem {
public:
    // Safe check if path exists using stat() instead of std::filesystem
    static bool exists_safe(const std::string& path) {
        if (path.empty() || path.size() > 4096) return false;
        
        struct stat path_stat;
        return (stat(path.c_str(), &path_stat) == 0);
    }
    
    // Safe check if path is directory
    static bool is_directory_safe(const std::string& path) {
        if (path.empty() || path.size() > 4096) return false;
        
        struct stat path_stat;
        if (stat(path.c_str(), &path_stat) != 0) return false;
        
        return S_ISDIR(path_stat.st_mode);
    }
    
    // Safe alternative to std::filesystem::recursive_directory_iterator
    static std::vector<std::string> find_files_safe(const std::string& directory, 
                                                    const std::vector<std::string>& extensions) {
        std::vector<std::string> files;
        
        if (!exists_safe(directory) || !is_directory_safe(directory)) {
            std::cerr << "Invalid directory: " << directory << std::endl;
            return files;
        }
        
        std::cout << "Using safe file finding for: " << directory << std::endl;
        
        // Use system command instead of std::filesystem to avoid heap corruption
        std::string find_command = "find \"" + directory + "\" -type f \\( ";
        for (size_t i = 0; i < extensions.size(); i++) {
            if (i > 0) find_command += " -o ";
            find_command += "-name \"*" + extensions[i] + "\"";
        }
        find_command += " \\) 2>/dev/null";
        
        std::cout << "Find command: " << find_command << std::endl;
        
        FILE* pipe = popen(find_command.c_str(), "r");
        if (!pipe) {
            std::cerr << "Failed to run find command" << std::endl;
            return files;
        }
        
        char path_buffer[4096];
        while (fgets(path_buffer, sizeof(path_buffer), pipe)) {
            // Remove newline
            size_t len = strlen(path_buffer);
            if (len > 0 && path_buffer[len-1] == '\n') {
                path_buffer[len-1] = '\0';
            }
            
            // Validate path
            if (len > 0 && len < 4000) {
                try {
                    files.emplace_back(path_buffer);
                    SafeMemoryTracker::track_allocation();
                } catch (const std::exception& e) {
                    std::cerr << "Error adding file: " << e.what() << std::endl;
                    break;
                }
            }
        }
        
        pclose(pipe);
        return files;
    }
};

// Test safe string operations
bool test_string_operations() {
    std::cout << "\n=== TESTING STRING OPERATIONS ===" << std::endl;
    SafeMemoryTracker::check_heap_integrity("string_test_start");
    
    try {
        // Test basic string operations
        std::string test_str = "test_genome_file.fna";
        SafeMemoryTracker::track_allocation();
        
        std::cout << "Original string: " << test_str << std::endl;
        
        // Test extension extraction
        size_t dot_pos = test_str.find_last_of('.');
        if (dot_pos != std::string::npos) {
            std::string ext = test_str.substr(dot_pos);
            std::cout << "Extension: " << ext << std::endl;
            SafeMemoryTracker::track_allocation();
        }
        
        // Test vector operations
        std::vector<std::string> test_vec;
        test_vec.reserve(100);  // Pre-allocate to avoid reallocations
        
        for (int i = 0; i < 10; i++) {
            test_vec.push_back("file_" + std::to_string(i) + ".fna");
            SafeMemoryTracker::track_allocation();
        }
        
        std::cout << "Created vector with " << test_vec.size() << " strings" << std::endl;
        
        // Test sorting
        std::sort(test_vec.begin(), test_vec.end());
        std::cout << "Sorted successfully" << std::endl;
        
        SafeMemoryTracker::check_heap_integrity("string_test_end");
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "String test failed: " << e.what() << std::endl;
        return false;
    }
}

// Test vector operations with capacity management
bool test_vector_operations() {
    std::cout << "\n=== TESTING VECTOR OPERATIONS ===" << std::endl;
    SafeMemoryTracker::check_heap_integrity("vector_test_start");
    
    try {
        std::vector<std::string> large_vector;
        
        // Conservative capacity planning
        const size_t expected_size = 1000;
        large_vector.reserve(expected_size);
        SafeMemoryTracker::track_allocation();
        
        std::cout << "Reserved capacity: " << large_vector.capacity() << std::endl;
        
        // Add strings without triggering reallocations
        for (size_t i = 0; i < expected_size; i++) {
            large_vector.emplace_back("/path/to/genome_" + std::to_string(i) + ".fna");
            
            if (i % 100 == 0) {
                SafeMemoryTracker::check_heap_integrity("vector_growth_" + std::to_string(i));
            }
        }
        
        std::cout << "Added " << large_vector.size() << " strings" << std::endl;
        std::cout << "Final capacity: " << large_vector.capacity() << std::endl;
        
        // Test sorting large vector
        std::sort(large_vector.begin(), large_vector.end());
        std::cout << "Sorted " << large_vector.size() << " strings successfully" << std::endl;
        
        SafeMemoryTracker::check_heap_integrity("vector_test_end");
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Vector test failed: " << e.what() << std::endl;
        SafeMemoryTracker::check_heap_integrity("vector_test_error");
        return false;
    }
}

// Test the exact operations from load_genome_files
bool test_genome_file_operations(const std::string& test_dir) {
    std::cout << "\n=== TESTING GENOME FILE OPERATIONS ===" << std::endl;
    std::cout << "Test directory: " << test_dir << std::endl;
    
    SafeMemoryTracker::check_heap_integrity("genome_test_start");
    
    if (!SafeFilesystem::exists_safe(test_dir)) {
        std::cerr << "Test directory does not exist: " << test_dir << std::endl;
        return false;
    }
    
    if (!SafeFilesystem::is_directory_safe(test_dir)) {
        std::cerr << "Path is not a directory: " << test_dir << std::endl;
        return false;
    }
    
    try {
        // Replicate the operations from load_genome_files
        std::vector<std::string> genome_files;
        std::vector<uint32_t> genome_taxon_ids;
        
        // Use safe file finding
        std::vector<std::string> extensions = {".fna", ".fa", ".fasta", ".ffn", ".faa"};
        auto found_files = SafeFilesystem::find_files_safe(test_dir, extensions);
        
        std::cout << "Found " << found_files.size() << " files" << std::endl;
        SafeMemoryTracker::check_heap_integrity("after_file_finding");
        
        if (found_files.empty()) {
            std::cout << "No genome files found - this is normal for test" << std::endl;
            return true;
        }
        
        // Reserve space to prevent reallocations
        genome_files.reserve(found_files.size());
        genome_taxon_ids.reserve(found_files.size());
        SafeMemoryTracker::track_allocation();
        SafeMemoryTracker::track_allocation();
        
        // Process files one by one
        for (size_t i = 0; i < found_files.size() && i < 100; i++) {  // Limit for testing
            const auto& file = found_files[i];
            
            // Validate file path
            if (file.empty() || file.size() > 4000) {
                std::cerr << "Invalid file path at index " << i << std::endl;
                continue;
            }
            
            // Extract taxon ID (simplified version)
            uint32_t taxon_id = 1000000 + i;  // Simple sequential ID for testing
            
            // Add to vectors
            genome_files.push_back(file);
            genome_taxon_ids.push_back(taxon_id);
            SafeMemoryTracker::track_allocation();
            SafeMemoryTracker::track_allocation();
            
            // Progress check
            if (i % 10 == 0) {
                SafeMemoryTracker::check_heap_integrity("processing_file_" + std::to_string(i));
            }
            
            // Verify vector consistency
            if (genome_files.size() != genome_taxon_ids.size()) {
                std::cerr << "FATAL: Vector size mismatch at index " << i << std::endl;
                return false;
            }
        }
        
        std::cout << "Successfully processed " << genome_files.size() << " files" << std::endl;
        
        // Test sorting (this often triggers heap corruption)
        try {
            std::sort(genome_files.begin(), genome_files.end());
            std::cout << "Sorting successful" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Sorting failed: " << e.what() << std::endl;
            return false;
        }
        
        SafeMemoryTracker::check_heap_integrity("genome_test_end");
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Genome file test failed: " << e.what() << std::endl;
        SafeMemoryTracker::check_heap_integrity("genome_test_error");
        return false;
    }
}

// Test std::filesystem operations specifically
bool test_filesystem_operations(const std::string& test_dir) {
    std::cout << "\n=== TESTING FILESYSTEM OPERATIONS ===" << std::endl;
    
    SafeMemoryTracker::check_heap_integrity("filesystem_test_start");
    
    try {
        // Test basic filesystem operations that might cause issues
        
        // 1. Test exists()
        std::cout << "Testing std::filesystem::exists()..." << std::endl;
        bool exists = std::filesystem::exists(test_dir);
        std::cout << "Directory exists: " << exists << std::endl;
        SafeMemoryTracker::check_heap_integrity("after_exists_check");
        
        if (!exists) {
            std::cout << "Directory doesn't exist, skipping filesystem tests" << std::endl;
            return true;
        }
        
        // 2. Test is_directory()
        std::cout << "Testing std::filesystem::is_directory()..." << std::endl;
        bool is_dir = std::filesystem::is_directory(test_dir);
        std::cout << "Is directory: " << is_dir << std::endl;
        SafeMemoryTracker::check_heap_integrity("after_is_directory_check");
        
        if (!is_dir) {
            std::cout << "Path is not a directory, skipping iterator tests" << std::endl;
            return true;
        }
        
        // 3. Test directory_iterator (non-recursive first)
        std::cout << "Testing std::filesystem::directory_iterator..." << std::endl;
        int file_count = 0;
        for (const auto& entry : std::filesystem::directory_iterator(test_dir)) {
            file_count++;
            if (file_count <= 5) {  // Only print first few
                std::cout << "  Entry: " << entry.path().string() << std::endl;
            }
            if (file_count >= 100) break;  // Safety limit
        }
        std::cout << "Found " << file_count << " entries with directory_iterator" << std::endl;
        SafeMemoryTracker::check_heap_integrity("after_directory_iterator");
        
        // 4. Test recursive_directory_iterator (the suspected culprit)
        std::cout << "Testing std::filesystem::recursive_directory_iterator..." << std::endl;
        int recursive_count = 0;
        try {
            for (const auto& entry : std::filesystem::recursive_directory_iterator(test_dir)) {
                recursive_count++;
                if (recursive_count <= 5) {  // Only print first few
                    std::cout << "  Recursive entry: " << entry.path().string() << std::endl;
                }
                if (recursive_count >= 100) {
                    std::cout << "  Stopping at 100 entries for safety" << std::endl;
                    break;
                }
            }
            std::cout << "Found " << recursive_count << " entries with recursive_directory_iterator" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "recursive_directory_iterator failed: " << e.what() << std::endl;
            return false;
        }
        SafeMemoryTracker::check_heap_integrity("after_recursive_iterator");
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Filesystem test failed: " << e.what() << std::endl;
        SafeMemoryTracker::check_heap_integrity("filesystem_test_error");
        return false;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "=== MEMORY SAFETY DEBUG TEST ===" << std::endl;
    
    // Enable heap checking
    SafeMemoryTracker::check_heap_integrity("program_start");
    
    // Determine test directory
    std::string test_dir = "/tmp";  // Default safe directory
    if (argc > 1) {
        test_dir = argv[1];
    }
    
    std::cout << "Using test directory: " << test_dir << std::endl;
    
    bool all_passed = true;
    
    // Run individual tests
    if (!test_string_operations()) {
        std::cerr << "❌ String operations test failed" << std::endl;
        all_passed = false;
    } else {
        std::cout << "✅ String operations test passed" << std::endl;
    }
    
    if (!test_vector_operations()) {
        std::cerr << "❌ Vector operations test failed" << std::endl;
        all_passed = false;
    } else {
        std::cout << "✅ Vector operations test passed" << std::endl;
    }
    
    if (!test_filesystem_operations(test_dir)) {
        std::cerr << "❌ Filesystem operations test failed" << std::endl;
        all_passed = false;
    } else {
        std::cout << "✅ Filesystem operations test passed" << std::endl;
    }
    
    if (!test_genome_file_operations(test_dir)) {
        std::cerr << "❌ Genome file operations test failed" << std::endl;
        all_passed = false;
    } else {
        std::cout << "✅ Genome file operations test passed" << std::endl;
    }
    
    SafeMemoryTracker::check_heap_integrity("program_end");
    SafeMemoryTracker::print_stats();
    
    if (all_passed) {
        std::cout << "\n🎉 All tests passed! The issue may be in the larger codebase." << std::endl;
        return 0;
    } else {
        std::cout << "\n💥 Some tests failed - found the problematic operations." << std::endl;
        return 1;
    }
}
