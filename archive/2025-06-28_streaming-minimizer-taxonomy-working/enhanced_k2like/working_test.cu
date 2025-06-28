#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

// Simulated functions for testing
bool load_genome_files(const std::string& path) {
    std::cout << "Simulating genome loading from: " << path << std::endl;
    return std::filesystem::exists(path);
}

int main(int argc, char* argv[]) {
    std::cout << "Working test version" << std::endl;
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " build --genome-dir <dir> --output <dir>" << std::endl;
        return 1;
    }
    
    std::string command(argv[1]);
    std::string genome_dir, output_dir;
    
    // Parse arguments
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--genome-dir" && i + 1 < argc) {
            genome_dir = argv[++i];
        } else if (std::string(argv[i]) == "--output" && i + 1 < argc) {
            output_dir = argv[++i];
        }
    }
    
    if (command == "build") {
        std::cout << "Build command:" << std::endl;
        std::cout << "  Genome dir: " << genome_dir << std::endl;
        std::cout << "  Output dir: " << output_dir << std::endl;
        
        if (!genome_dir.empty()) {
            bool success = load_genome_files(genome_dir);
            std::cout << "Genome loading: " << (success ? "SUCCESS" : "FAILED") << std::endl;
        }
        
        std::cout << "Build simulation completed!" << std::endl;
    }
    
    return 0;
}
