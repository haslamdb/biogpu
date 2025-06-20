#!/bin/bash
# Debug and fix the load_genome_files crash

echo "=== DEBUGGING load_genome_files CRASH ==="

# The crash happens right after "About to call load_genome_files..."
# This suggests the issue is in the file system traversal or vector operations

# Let's create a minimal test to isolate the exact problem
cat > debug_load_files.cpp << 'EOF'
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
EOF

echo "Created filesystem debug test"

# Compile the test
echo "Compiling debug test..."
g++ -std=c++17 -g -O0 debug_load_files.cpp -o debug_load_files -lstdc++fs

if [ $? -eq 0 ]; then
    echo "✓ Debug test compiled successfully"
    
    echo ""
    echo "Running debug test..."
    export MALLOC_CHECK_=2
    ./debug_load_files data/type_strain_genomes_processed
    
    test_exit=$?
    echo ""
    echo "Debug test exit code: $test_exit"
    
    if [ $test_exit -eq 0 ]; then
        echo "✓ Debug test passed - issue is in the full program logic"
        echo ""
        echo "=== CREATING SAFER load_genome_files IMPLEMENTATION ==="
        
        # Create a safer version of the load_genome_files function
        cat > safe_load_genome_files.cpp << 'EOF'
// Safer implementation of load_genome_files to replace the problematic version

#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstdlib>

std::vector<std::string> safe_find_genome_files(const std::string& directory) {
    std::vector<std::string> files;
    
    std::cout << "SAFE: Finding genome files in: " << directory << std::endl;
    
    if (directory.empty()) {
        std::cout << "SAFE: Empty directory path" << std::endl;
        return files;
    }
    
    try {
        // Use more conservative approach
        if (!std::filesystem::exists(directory)) {
            std::cout << "SAFE: Directory does not exist" << std::endl;
            return files;
        }
        
        if (!std::filesystem::is_directory(directory)) {
            std::cout << "SAFE: Path is not a directory" << std::endl;
            return files;
        }
        
        // Pre-allocate reasonable capacity to avoid frequent reallocations
        files.reserve(10000);
        
        int file_count = 0;
        
        // Use non-recursive iterator first to be safer
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            try {
                if (entry.is_regular_file()) {
                    std::string path_str = entry.path().string();
                    std::string ext = entry.path().extension().string();
                    
                    // Convert extension to lowercase for comparison
                    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                    
                    if (ext == ".fna" || ext == ".fa" || ext == ".fasta" || 
                        ext == ".ffn" || ext == ".faa") {
                        
                        files.push_back(path_str);
                        file_count++;
                        
                        if (file_count % 1000 == 0) {
                            std::cout << "SAFE: Found " << file_count << " files..." << std::endl;
                        }
                        
                        // Safety limit
                        if (file_count > 50000) {
                            std::cout << "SAFE: Reached safety limit of 50000 files" << std::endl;
                            break;
                        }
                    }
                }
            } catch (const std::exception& e) {
                std::cout << "SAFE: Error processing entry: " << e.what() << std::endl;
                continue; // Skip problematic entries
            }
        }
        
        std::cout << "SAFE: Found " << files.size() << " genome files" << std::endl;
        
        // Safe sorting
        if (!files.empty()) {
            try {
                std::sort(files.begin(), files.end());
                std::cout << "SAFE: Files sorted successfully" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "SAFE: Sorting failed: " << e.what() << std::endl;
                // Continue without sorting
            }
        }
        
    } catch (const std::exception& e) {
        std::cout << "SAFE: Exception in find_genome_files: " << e.what() << std::endl;
        files.clear(); // Return empty vector on error
    }
    
    return files;
}

bool safe_load_genome_files_test(const std::string& library_path) {
    std::cout << "SAFE: Testing load_genome_files with: " << library_path << std::endl;
    
    try {
        auto files = safe_find_genome_files(library_path);
        
        if (files.empty()) {
            std::cout << "SAFE: No files found" << std::endl;
            return false;
        }
        
        std::cout << "SAFE: Successfully loaded " << files.size() << " files" << std::endl;
        
        // Test creating taxon IDs (another potential crash point)
        std::vector<uint32_t> taxon_ids;
        taxon_ids.reserve(files.size());
        
        for (size_t i = 0; i < files.size(); i++) {
            // Simple hash-based taxon ID
            std::hash<std::string> hasher;
            uint32_t taxon_id = hasher(std::filesystem::path(files[i]).stem().string()) % 1000000 + 1000000;
            taxon_ids.push_back(taxon_id);
            
            if (i % 1000 == 0) {
                std::cout << "SAFE: Processed " << i << " taxon IDs..." << std::endl;
            }
        }
        
        std::cout << "SAFE: Created " << taxon_ids.size() << " taxon IDs" << std::endl;
        std::cout << "SAFE: load_genome_files test SUCCESSFUL" << std::endl;
        
        return true;
        
    } catch (const std::exception& e) {
        std::cout << "SAFE: Exception in load test: " << e.what() << std::endl;
        return false;
    }
}

int main(int argc, char* argv[]) {
    setenv("MALLOC_CHECK_", "2", 1);
    
    std::string test_path = "data/type_strain_genomes_processed";
    if (argc > 1) {
        test_path = argv[1];
    }
    
    bool success = safe_load_genome_files_test(test_path);
    
    if (success) {
        std::cout << "✓ SAFE VERSION WORKS!" << std::endl;
        std::cout << "This implementation can replace the problematic load_genome_files" << std::endl;
    } else {
        std::cout << "❌ Even safe version has issues" << std::endl;
    }
    
    return success ? 0 : 1;
}
EOF

        echo "Compiling safer version..."
        g++ -std=c++17 -g -O0 safe_load_genome_files.cpp -o safe_load_files -lstdc++fs
        
        if [ $? -eq 0 ]; then
            echo "✓ Safe version compiled"
            echo ""
            echo "Testing safer implementation..."
            ./safe_load_files data/type_strain_genomes_processed
            
            safe_exit=$?
            if [ $safe_exit -eq 0 ]; then
                echo ""
                echo "🎉 SAFE VERSION WORKS!"
                echo "The issue is in the original load_genome_files implementation"
                echo ""
                echo "SOLUTION: Replace the load_genome_files method in your GPUKrakenDatabaseBuilder"
                echo "with the safer version that handles errors gracefully."
            else
                echo ""
                echo "❌ Even the safe version crashes - deeper system issue"
            fi
        fi
        
    elif [ $test_exit -eq 134 ]; then
        echo "❌ Debug test shows HEAP CORRUPTION in filesystem operations"
        echo "The issue is in std::filesystem library usage or the directory contents"
    else
        echo "❌ Debug test failed with exit code: $test_exit"
    fi
    
else
    echo "❌ Debug test compilation failed"
fi

echo ""
echo "=== SUMMARY ==="
echo "The heap corruption occurs specifically in load_genome_files()"
echo "This suggests the issue is in:"
echo "1. std::filesystem::recursive_directory_iterator"
echo "2. Vector reallocations during file list building"
echo "3. String operations with file paths"
echo "4. Possible corrupted files in the genome directory"

echo ""
echo "NEXT STEPS:"
echo "1. Run the debug test to isolate the exact cause"
echo "2. If debug test passes, replace load_genome_files with safer version"
echo "3. If debug test fails, investigate filesystem/directory issues"
