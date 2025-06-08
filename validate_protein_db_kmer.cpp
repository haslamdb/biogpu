// Add this function to translated_search_revised.cu or as a separate utility

#include <fstream>
#include <json/json.h>

extern "C" {
    bool validate_protein_db_kmer_size(const char* db_path, int expected_kmer_size) {
        std::string metadata_path = std::string(db_path) + "/metadata.json";
        std::ifstream metadata_file(metadata_path);
        
        if (!metadata_file.is_open()) {
            std::cerr << "Warning: Could not open metadata file: " << metadata_path << std::endl;
            return false;
        }
        
        Json::Value root;
        Json::Reader reader;
        
        if (!reader.parse(metadata_file, root)) {
            std::cerr << "Warning: Could not parse metadata JSON\n";
            metadata_file.close();
            return false;
        }
        
        metadata_file.close();
        
        if (!root.isMember("kmer_length")) {
            std::cerr << "Warning: Metadata does not contain kmer_length field\n";
            return false;
        }
        
        int db_kmer_size = root["kmer_length"].asInt();
        
        if (db_kmer_size != expected_kmer_size) {
            std::cerr << "ERROR: Protein database was built with k-mer size " << db_kmer_size 
                      << " but pipeline expects k-mer size " << expected_kmer_size << std::endl;
            return false;
        }
        
        return true;
    }
    
    // Alternative simpler version without JSON dependency
    bool validate_protein_db_kmer_size_simple(const char* db_path, int expected_kmer_size) {
        std::string metadata_path = std::string(db_path) + "/metadata.json";
        std::ifstream metadata_file(metadata_path);
        
        if (!metadata_file.is_open()) {
            return false;
        }
        
        std::string line;
        while (std::getline(metadata_file, line)) {
            // Look for "kmer_length": X pattern
            size_t pos = line.find("\"kmer_length\":");
            if (pos != std::string::npos) {
                // Extract the number after the colon
                size_t num_start = line.find_first_of("0123456789", pos);
                if (num_start != std::string::npos) {
                    int db_kmer_size = std::stoi(line.substr(num_start));
                    metadata_file.close();
                    
                    if (db_kmer_size != expected_kmer_size) {
                        std::cerr << "ERROR: Protein database k-mer size mismatch! "
                                  << "Database: " << db_kmer_size 
                                  << ", Expected: " << expected_kmer_size << std::endl;
                        return false;
                    }
                    return true;
                }
            }
        }
        
        metadata_file.close();
        return false;
    }
}