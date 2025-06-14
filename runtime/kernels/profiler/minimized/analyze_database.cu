#include <iostream>
#include <fstream>
#include <iomanip>

// Database analyzer to determine if streaming is needed
// Analyzes minimizer database size and provides recommendations

class DatabaseAnalyzer {
public:
    struct DatabaseStats {
        uint64_t total_file_size_mb;
        uint32_t num_organisms;
        uint64_t num_minimizer_hashes;
        uint64_t total_minimizer_entries;
        uint64_t estimated_gpu_memory_mb;
        bool needs_streaming;
        uint32_t recommended_chunks;
        uint64_t chunk_size_mb;
    };
    
    static DatabaseStats analyze_database(const std::string& db_path) {
        std::cout << "Analyzing database: " << db_path << std::endl;
        
        DatabaseStats stats = {};
        
        // Get file size
        std::ifstream file(db_path, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open database file");
        }
        
        stats.total_file_size_mb = file.tellg() / (1024 * 1024);
        file.seekg(0, std::ios::beg);
        
        // Read header
        uint32_t magic, version, k_size, m_size;
        uint64_t num_minimizer_hashes;
        
        file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        file.read(reinterpret_cast<char*>(&version), sizeof(version));
        file.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        file.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        file.read(reinterpret_cast<char*>(&stats.num_organisms), sizeof(stats.num_organisms));
        file.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        stats.num_minimizer_hashes = num_minimizer_hashes;
        
        if (magic != 0x4D494E49) {
            throw std::runtime_error("Invalid database format");
        }
        
        // Skip organism metadata to count minimizer entries
        skip_organism_metadata(file, stats.num_organisms);
        
        // Count total minimizer entries
        stats.total_minimizer_entries = count_minimizer_entries(file);
        
        // Estimate GPU memory requirements
        size_t gpu_memory_per_entry = sizeof(uint64_t) + sizeof(uint32_t) + sizeof(float);  // 16 bytes
        stats.estimated_gpu_memory_mb = (stats.total_minimizer_entries * gpu_memory_per_entry) / (1024 * 1024);
        
        // Add overhead for organism arrays and working memory
        stats.estimated_gpu_memory_mb += 500;  // 500MB overhead
        
        // Determine if streaming is needed (leave 16GB for reads and working memory)
        size_t available_gpu_memory_mb = 32 * 1024;  // Assume 32GB available for database
        stats.needs_streaming = stats.estimated_gpu_memory_mb > available_gpu_memory_mb;
        
        if (stats.needs_streaming) {
            stats.recommended_chunks = (stats.estimated_gpu_memory_mb / 2048) + 1;  // 2GB chunks
            stats.chunk_size_mb = stats.estimated_gpu_memory_mb / stats.recommended_chunks;
        } else {
            stats.recommended_chunks = 1;
            stats.chunk_size_mb = stats.estimated_gpu_memory_mb;
        }
        
        file.close();
        return stats;
    }
    
    static void print_analysis(const DatabaseStats& stats) {
        std::cout << "\n=== DATABASE ANALYSIS REPORT ===" << std::endl;
        std::cout << "File size: " << stats.total_file_size_mb << " MB" << std::endl;
        std::cout << "Organisms: " << stats.num_organisms << std::endl;
        std::cout << "Unique minimizer hashes: " << stats.num_minimizer_hashes << std::endl;
        std::cout << "Total minimizer entries: " << stats.total_minimizer_entries << std::endl;
        std::cout << "Estimated GPU memory needed: " << stats.estimated_gpu_memory_mb << " MB" << std::endl;
        
        std::cout << "\n=== RECOMMENDATION ===" << std::endl;
        if (stats.needs_streaming) {
            std::cout << "⚠️  STREAMING REQUIRED" << std::endl;
            std::cout << "Database is too large for single GPU load" << std::endl;
            std::cout << "Recommended configuration:" << std::endl;
            std::cout << "- Number of chunks: " << stats.recommended_chunks << std::endl;
            std::cout << "- Chunk size: " << stats.chunk_size_mb << " MB" << std::endl;
            std::cout << "\nUse: ./streaming_gpu_profiler" << std::endl;
        } else {
            std::cout << "✅ DIRECT GPU LOADING POSSIBLE" << std::endl;
            std::cout << "Database fits comfortably in GPU memory" << std::endl;
            std::cout << "Available GPU memory buffer: " 
                      << (32 * 1024 - stats.estimated_gpu_memory_mb) << " MB" << std::endl;
            std::cout << "\nUse: ./gpu_community_profiler" << std::endl;
        }
        
        std::cout << "\n=== PERFORMANCE ESTIMATES ===" << std::endl;
        if (stats.needs_streaming) {
            std::cout << "Expected processing time: " << (stats.recommended_chunks * 2) << "-" 
                      << (stats.recommended_chunks * 5) << " seconds per 1M reads" << std::endl;
            std::cout << "Memory efficiency: High (streaming)" << std::endl;
        } else {
            std::cout << "Expected processing time: 1-3 seconds per 1M reads" << std::endl;
            std::cout << "Memory efficiency: Maximum (direct GPU)" << std::endl;
        }
        
        // Organism density analysis
        float entries_per_organism = (float)stats.total_minimizer_entries / stats.num_organisms;
        std::cout << "Average minimizers per organism: " << std::fixed << std::setprecision(0) 
                  << entries_per_organism << std::endl;
        
        if (entries_per_organism > 100000) {
            std::cout << "Note: High minimizer density - consider increasing uniqueness filtering" << std::endl;
        } else if (entries_per_organism < 10000) {
            std::cout << "Note: Low minimizer density - good for fast processing" << std::endl;
        }
    }
    
private:
    static void skip_organism_metadata(std::ifstream& file, uint32_t num_organisms) {
        for (uint32_t i = 0; i < num_organisms; i++) {
            // Skip fixed-size fields
            file.seekg(sizeof(uint32_t) + sizeof(uint32_t) + sizeof(float) + 
                      sizeof(uint64_t) + sizeof(uint32_t), std::ios::cur);
            
            // Skip variable-length name
            uint16_t name_length;
            file.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            file.seekg(name_length, std::ios::cur);
            
            // Skip variable-length taxonomy path
            uint16_t taxonomy_length;
            file.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            file.seekg(taxonomy_length, std::ios::cur);
        }
    }
    
    static uint64_t count_minimizer_entries(std::ifstream& file) {
        uint64_t total_entries = 0;
        uint64_t hash;
        uint32_t num_entries_for_hash;
        
        std::cout << "Counting minimizer entries..." << std::endl;
        
        while (file.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
            file.read(reinterpret_cast<char*>(&num_entries_for_hash), sizeof(num_entries_for_hash));
            
            // Skip the actual entries
            file.seekg(num_entries_for_hash * (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(uint8_t)), 
                      std::ios::cur);
            
            total_entries += num_entries_for_hash;
            
            if (total_entries % 1000000 == 0) {
                std::cout << "\rCounted " << total_entries << " entries..." << std::flush;
            }
        }
        
        std::cout << "\nFinished counting: " << total_entries << " total entries" << std::endl;
        return total_entries;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <minimizer_database>" << std::endl;
        std::cerr << "\nAnalyzes database size and recommends profiling strategy" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    
    try {
        auto stats = DatabaseAnalyzer::analyze_database(database_path);
        DatabaseAnalyzer::print_analysis(stats);
        
        // Generate configuration file for streaming profiler
        if (stats.needs_streaming) {
            std::ofstream config("streaming_config.txt");
            config << "# Streaming GPU Profiler Configuration\n";
            config << "# Generated by database analyzer\n\n";
            config << "max_gpu_memory_gb=" << (stats.chunk_size_mb / 1024) << "\n";
            config << "chunk_size_mb=" << stats.chunk_size_mb << "\n";
            config << "recommended_chunks=" << stats.recommended_chunks << "\n";
            config << "total_database_size_mb=" << stats.estimated_gpu_memory_mb << "\n";
            config.close();
            
            std::cout << "\nGenerated streaming_config.txt with optimal parameters" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}