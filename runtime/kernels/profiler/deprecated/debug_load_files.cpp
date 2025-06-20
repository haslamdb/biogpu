#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <cstdlib>

// Test 1: Basic filesystem operations
bool test_filesystem_access(const std::string& path) {
    std::cout << "Testing filesystem access to: " << path << std::endl;
    
    try {
        // Check if path exists
        if (!std::filesystem::exists(path)) {
            std::cout << "  Path does not exist" << std::endl;
            return false;
        }
        
        std::cout << "  Path exists" << std::endl;
        
        // Check if it's a directory
        if (!std::filesystem::is_directory(path)) {
            std::cout << "  Path is not a directory" << std::endl;
            return false;
        }
        
        std::cout << "  Path is a directory" << std::endl;
        
        // Try to count files (this is where crashes often occur)
        int file_count = 0;
        std::cout << "  Attempting directory iteration..." << std::endl;
        
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            file_count++;
            if (file_count > 100) break; // Safety limit
        }
        
        std::cout << "  Found " << file_count << " entries" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cout << "  EXCEPTION: " << e.what() << std::endl;
        return false;
    }
}

// Test 2: Vector operations that might cause heap corruption
bool test_vector_operations() {
    std::cout << "Testing vector operations..." << std::endl;
    
    try {
        std::vector<std::string> test_files;
        
        // Reserve space (this sometimes causes issues)
        std::cout << "  Reserving vector space..." << std::endl;
        test_files.reserve(1000);
        
        // Add some strings
        std::cout << "  Adding test strings..." << std::endl;
        for (int i = 0; i < 100; i++) {
            test_files.push_back("test_file_" + std::to_string(i) + ".fasta");
        }
        
        std::cout << "  Vector size: " << test_files.size() << std::endl;
        std::cout << "  Vector capacity: " << test_files.capacity() << std::endl;
        
        // Clear and test again
        test_files.clear();
        test_files.shrink_to_fit();
        
        std::cout << "  Vector operations successful" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cout << "  EXCEPTION in vector operations: " << e.what() << std::endl;
        return false;
    }
}

// Test 3: Combined filesystem + vector operations (mimics load_genome_files)
bool test_combined_operations(const std::string& path) {
    std::cout << "Testing combined filesystem + vector operations..." << std::endl;
    
    try {
        std::vector<std::string> found_files;
        found_files.reserve(1000);
        
        int count = 0;
        for (const auto& entry : std::filesystem::recursive_directory_iterator(path)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().string();
                std::string ext = entry.path().extension().string();
                
                if (ext == ".fna" || ext == ".fa" || ext == ".fasta") {
                    found_files.push_back(filename);
                    count++;
                    
                    if (count % 100 == 0) {
                        std::cout << "  Processed " << count << " files..." << std::endl;
                    }
                    
                    // Safety limit to prevent infinite processing
                    if (count > 10000) {
                        std::cout << "  Reached safety limit of 10000 files" << std::endl;
                        break;
                    }
                }
            }
        }
        
        std::cout << "  Found " << found_files.size() << " genome files" << std::endl;
        
        // Test sorting (another operation that can cause issues)
        std::cout << "  Sorting files..." << std::endl;
        std::sort(found_files.begin(), found_files.end());
        
        std::cout << "  Combined operations successful" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cout << "  EXCEPTION in combined operations: " << e.what() << std::endl;
        return false;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "=== FILESYSTEM CORRUPTION DEBUG TEST ===" << std::endl;
    
    // Set up malloc debugging
    setenv("MALLOC_CHECK_", "2", 1);
    
    std::string test_path = "data/type_strain_genomes_processed";
    if (argc > 1) {
        test_path = argv[1];
    }
    
    std::cout << "Testing path: " << test_path << std::endl << std::endl;
    
    // Run tests in sequence
    bool test1 = test_filesystem_access(test_path);
    std::cout << "Test 1 (filesystem): " << (test1 ? "PASSED" : "FAILED") << std::endl << std::endl;
    
    bool test2 = test_vector_operations();
    std::cout << "Test 2 (vectors): " << (test2 ? "PASSED" : "FAILED") << std::endl << std::endl;
    
    if (test1 && test2) {
        bool test3 = test_combined_operations(test_path);
        std::cout << "Test 3 (combined): " << (test3 ? "PASSED" : "FAILED") << std::endl;
        
        if (test3) {
            std::cout << std::endl << "✓ ALL TESTS PASSED - The issue is elsewhere" << std::endl;
        } else {
            std::cout << std::endl << "❌ CRASH IN COMBINED OPERATIONS" << std::endl;
        }
    } else {
        std::cout << std::endl << "❌ BASIC OPERATIONS FAILED" << std::endl;
    }
    
    return 0;
}
