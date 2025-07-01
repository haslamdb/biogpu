// test_basic_database_build_debug.cu
// Debug version to test file finding

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

int main() {
    std::string test_genomes_dir = "/home/david/Documents/Code/biogpu/data/type_strain_reference_genomes";
    
    std::cout << "=== Debug: Testing file finding ===" << std::endl;
    std::cout << "Directory: " << test_genomes_dir << std::endl;
    
    // Check if directory exists
    if (!std::filesystem::exists(test_genomes_dir)) {
        std::cerr << "Directory does not exist!" << std::endl;
        return 1;
    }
    
    if (!std::filesystem::is_directory(test_genomes_dir)) {
        std::cerr << "Path is not a directory!" << std::endl;
        return 1;
    }
    
    // List all files
    std::cout << "\nAll files in directory:" << std::endl;
    int count = 0;
    try {
        for (const auto& entry : std::filesystem::directory_iterator(test_genomes_dir)) {
            if (entry.is_regular_file()) {
                std::cout << "  " << entry.path().filename().string() 
                          << " (ext: " << entry.path().extension().string() << ")" << std::endl;
                count++;
                if (count >= 5) {
                    std::cout << "  ... (showing first 5 files)" << std::endl;
                    break;
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error iterating directory: " << e.what() << std::endl;
        return 1;
    }
    
    // Test path for dangerous characters
    std::cout << "\nTesting path safety:" << std::endl;
    const std::string dangerous_chars = ";|&`$<>\\";
    if (test_genomes_dir.find_first_of(dangerous_chars) != std::string::npos) {
        std::cerr << "Path contains dangerous characters!" << std::endl;
    } else {
        std::cout << "Path is safe" << std::endl;
    }
    
    if (test_genomes_dir.find("..") != std::string::npos) {
        std::cerr << "Path contains directory traversal!" << std::endl;
    } else {
        std::cout << "No directory traversal detected" << std::endl;
    }
    
    if (test_genomes_dir.length() > 4096) {
        std::cerr << "Path is too long!" << std::endl;
    } else {
        std::cout << "Path length is OK: " << test_genomes_dir.length() << " characters" << std::endl;
    }
    
    return 0;
}