// extract_fq_genes_fasta.cpp
// Extracts nucleotide sequences from FQ genes JSON files and saves as FASTA
// This allows us to use the same bloom filter builder for both pipelines

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

void process_json_file(const std::string& json_path, std::ofstream& output) {
    std::ifstream file(json_path);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open " << json_path << std::endl;
        return;
    }
    
    try {
        json data;
        file >> data;
        
        // Extract species and gene from the data
        for (auto& [key, seq_data] : data.items()) {
            if (seq_data.contains("nucleotide") && seq_data.contains("species") && seq_data.contains("gene")) {
                std::string species = seq_data["species"];
                std::string gene = seq_data["gene"];
                std::string sequence = seq_data["nucleotide"];
                std::string accession = seq_data.value("accession", key);
                
                // Write FASTA header
                output << ">" << accession << " " << species << " " << gene << std::endl;
                
                // Write sequence (80 chars per line)
                for (size_t i = 0; i < sequence.length(); i += 80) {
                    output << sequence.substr(i, 80) << std::endl;
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON file " << json_path << ": " << e.what() << std::endl;
    }
    
    file.close();
}

void process_directory(const std::string& dir_path, std::ofstream& output) {
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << dir_path << std::endl;
        return;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        
        // Skip . and ..
        if (filename == "." || filename == "..") continue;
        
        std::string full_path = dir_path + "/" + filename;
        struct stat st;
        if (stat(full_path.c_str(), &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                // Process subdirectory
                process_directory(full_path, output);
            } else if (filename.find(".json") != std::string::npos) {
                std::cout << "Processing: " << full_path << std::endl;
                process_json_file(full_path, output);
            }
        }
    }
    
    closedir(dir);
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fq_genes_directory> <output.fasta>" << std::endl;
        std::cerr << "Example: " << argv[0] << " ../../data/fq_genes fq_resistance_genes.fasta" << std::endl;
        return 1;
    }
    
    std::string input_dir = argv[1];
    std::string output_file = argv[2];
    
    std::ofstream output(output_file);
    if (!output.is_open()) {
        std::cerr << "Error: Cannot create output file " << output_file << std::endl;
        return 1;
    }
    
    std::cout << "Extracting sequences from " << input_dir << " to " << output_file << std::endl;
    process_directory(input_dir, output);
    
    output.close();
    std::cout << "Done!" << std::endl;
    
    return 0;
}