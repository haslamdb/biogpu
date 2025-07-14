// build_kmer_index.cpp
// Builds a k-mer index from nucleotide sequences (FASTA files or directories)
// Creates kmer_index.bin file compatible with existing pipelines

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <dirent.h>
#include <sys/stat.h>

// K-mer configuration
const int DEFAULT_KMER_LENGTH = 15;

// Base encoding
const uint8_t BASE_A = 0;
const uint8_t BASE_C = 1;
const uint8_t BASE_G = 2;
const uint8_t BASE_T = 3;
const uint8_t BASE_INVALID = 255;

struct KmerIndexHeader {
    uint32_t kmer_length;
    uint32_t num_kmers;
    uint32_t num_sequences;
    uint32_t has_reverse_complement;
};

struct KmerEntry {
    uint64_t kmer;
    uint32_t sequence_id;
    uint32_t position;
};

// Encode a nucleotide base
inline uint8_t encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return BASE_A;
        case 'C': case 'c': return BASE_C;
        case 'G': case 'g': return BASE_G;
        case 'T': case 't': case 'U': case 'u': return BASE_T;
        default: return BASE_INVALID;
    }
}

// Encode a k-mer to 64-bit integer
uint64_t encode_kmer(const std::string& sequence, size_t start, int kmer_length) {
    uint64_t encoded = 0;
    
    for (int i = 0; i < kmer_length; i++) {
        uint8_t base = encode_base(sequence[start + i]);
        if (base == BASE_INVALID) {
            return UINT64_MAX; // Invalid k-mer
        }
        encoded = (encoded << 2) | base;
    }
    
    return encoded;
}

// Compute reverse complement of encoded k-mer
uint64_t reverse_complement_kmer(uint64_t kmer, int kmer_length) {
    uint64_t rc = 0;
    
    for (int i = 0; i < kmer_length; i++) {
        uint8_t base = kmer & 3;
        uint8_t complement;
        
        switch (base) {
            case BASE_A: complement = BASE_T; break;
            case BASE_T: complement = BASE_A; break;
            case BASE_C: complement = BASE_G; break;
            case BASE_G: complement = BASE_C; break;
            default: return UINT64_MAX;
        }
        
        rc = (rc << 2) | complement;
        kmer >>= 2;
    }
    
    return rc;
}

class KmerIndexBuilder {
private:
    int kmer_length;
    bool include_reverse_complement;
    std::vector<std::string> sequences;
    std::vector<std::string> sequence_names;
    std::unordered_map<uint64_t, std::vector<KmerEntry>> kmer_map;
    
public:
    KmerIndexBuilder(int k = DEFAULT_KMER_LENGTH, bool include_rc = true) 
        : kmer_length(k), include_reverse_complement(include_rc) {
        if (k > 31) {
            throw std::invalid_argument("K-mer length must be <= 31 for 64-bit encoding");
        }
    }
    
    // Add sequences from FASTA file
    void add_fasta_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return;
        }
        
        std::string line, current_seq, current_name;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Save previous sequence if exists
                if (!current_seq.empty()) {
                    sequences.push_back(current_seq);
                    sequence_names.push_back(current_name);
                    current_seq.clear();
                }
                current_name = line.substr(1);
            } else {
                // Append to current sequence
                current_seq += line;
            }
        }
        
        // Save last sequence
        if (!current_seq.empty()) {
            sequences.push_back(current_seq);
            sequence_names.push_back(current_name);
        }
        
        file.close();
        std::cout << "Loaded " << sequences.size() << " sequences from " << filename << std::endl;
    }
    
    // Process directory of FASTA files
    void add_fasta_directory(const std::string& dir_path) {
        DIR* dir = opendir(dir_path.c_str());
        if (!dir) {
            std::cerr << "Error: Cannot open directory " << dir_path << std::endl;
            return;
        }
        
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string filename = entry->d_name;
            
            if (filename == "." || filename == "..") continue;
            
            std::string full_path = dir_path + "/" + filename;
            struct stat st;
            if (stat(full_path.c_str(), &st) == 0) {
                if (S_ISDIR(st.st_mode)) {
                    add_fasta_directory(full_path);
                } else if (filename.find(".fa") != std::string::npos || 
                           filename.find(".fasta") != std::string::npos ||
                           filename.find(".fna") != std::string::npos) {
                    add_fasta_file(full_path);
                }
            }
        }
        
        closedir(dir);
    }
    
    // Build k-mer index from loaded sequences
    void build_index() {
        kmer_map.clear();
        
        for (uint32_t seq_id = 0; seq_id < sequences.size(); seq_id++) {
            const std::string& seq = sequences[seq_id];
            
            if (seq.length() < (size_t)kmer_length) continue;
            
            for (size_t pos = 0; pos <= seq.length() - kmer_length; pos++) {
                // Forward k-mer
                uint64_t kmer = encode_kmer(seq, pos, kmer_length);
                if (kmer != UINT64_MAX) {
                    kmer_map[kmer].push_back({kmer, seq_id, (uint32_t)pos});
                    
                    // Reverse complement
                    if (include_reverse_complement) {
                        uint64_t rc = reverse_complement_kmer(kmer, kmer_length);
                        if (rc != UINT64_MAX && rc != kmer) { // Don't add palindromes twice
                            kmer_map[rc].push_back({rc, seq_id, (uint32_t)pos});
                        }
                    }
                }
            }
        }
        
        std::cout << "Built index with " << kmer_map.size() << " unique k-mers" << std::endl;
    }
    
    // Save k-mer index to binary file
    void save_index(const std::string& output_file) {
        std::ofstream file(output_file, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create output file " << output_file << std::endl;
            return;
        }
        
        // Prepare sorted k-mer list
        std::vector<uint64_t> sorted_kmers;
        sorted_kmers.reserve(kmer_map.size());
        for (const auto& pair : kmer_map) {
            sorted_kmers.push_back(pair.first);
        }
        std::sort(sorted_kmers.begin(), sorted_kmers.end());
        
        // Write header
        KmerIndexHeader header;
        header.kmer_length = kmer_length;
        header.num_kmers = sorted_kmers.size();
        header.num_sequences = sequences.size();
        header.has_reverse_complement = include_reverse_complement ? 1 : 0;
        
        file.write(reinterpret_cast<const char*>(&header), sizeof(header));
        
        // Write sorted k-mers
        file.write(reinterpret_cast<const char*>(sorted_kmers.data()), 
                   sorted_kmers.size() * sizeof(uint64_t));
        
        // Write k-mer positions (for each k-mer, write count and then positions)
        for (const uint64_t& kmer : sorted_kmers) {
            const auto& entries = kmer_map[kmer];
            uint32_t count = entries.size();
            file.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
            
            for (const auto& entry : entries) {
                file.write(reinterpret_cast<const char*>(&entry.sequence_id), sizeof(uint32_t));
                file.write(reinterpret_cast<const char*>(&entry.position), sizeof(uint32_t));
            }
        }
        
        file.close();
        
        std::cout << "Saved k-mer index to " << output_file << std::endl;
        std::cout << "Index contains:" << std::endl;
        std::cout << "  - " << header.num_kmers << " unique k-mers" << std::endl;
        std::cout << "  - " << header.num_sequences << " sequences" << std::endl;
        std::cout << "  - K-mer length: " << header.kmer_length << std::endl;
        std::cout << "  - Includes RC: " << (header.has_reverse_complement ? "Yes" : "No") << std::endl;
    }
    
    // Save sequence metadata
    void save_metadata(const std::string& output_dir) {
        std::string metadata_file = output_dir + "/sequences.txt";
        std::ofstream file(metadata_file);
        if (!file.is_open()) {
            std::cerr << "Warning: Cannot create metadata file" << std::endl;
            return;
        }
        
        for (size_t i = 0; i < sequence_names.size(); i++) {
            file << i << "\t" << sequence_names[i] << "\t" << sequences[i].length() << std::endl;
        }
        
        file.close();
        std::cout << "Saved sequence metadata to " << metadata_file << std::endl;
    }
};

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <input> <output_dir> [options]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Arguments:" << std::endl;
    std::cerr << "  input:       FASTA file or directory containing FASTA files" << std::endl;
    std::cerr << "  output_dir:  Directory to save kmer_index.bin and metadata" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -k <int>     K-mer length (default: 15)" << std::endl;
    std::cerr << "  --no-rc      Don't include reverse complement k-mers" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Examples:" << std::endl;
    std::cerr << "  " << program_name << " ../../data/AMR_CDS.fa amr_index/" << std::endl;
    std::cerr << "  " << program_name << " ../../data/fq_genes/ fq_index/ -k 15" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string input_path = argv[1];
    std::string output_dir = argv[2];
    
    // Parse options
    int kmer_length = DEFAULT_KMER_LENGTH;
    bool include_rc = true;
    
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-k" && i + 1 < argc) {
            kmer_length = std::stoi(argv[++i]);
            if (kmer_length > 31) {
                std::cerr << "Error: K-mer length must be <= 31" << std::endl;
                return 1;
            }
        } else if (arg == "--no-rc") {
            include_rc = false;
        }
    }
    
    // Create output directory if it doesn't exist
    struct stat st;
    if (stat(output_dir.c_str(), &st) != 0) {
        std::string mkdir_cmd = "mkdir -p " + output_dir;
        system(mkdir_cmd.c_str());
    }
    
    // Build index
    KmerIndexBuilder builder(kmer_length, include_rc);
    
    // Check if input is file or directory
    if (stat(input_path.c_str(), &st) != 0) {
        std::cerr << "Error: Cannot access " << input_path << std::endl;
        return 1;
    }
    
    if (S_ISDIR(st.st_mode)) {
        std::cout << "Processing directory: " << input_path << std::endl;
        builder.add_fasta_directory(input_path);
    } else {
        std::cout << "Processing file: " << input_path << std::endl;
        builder.add_fasta_file(input_path);
    }
    
    // Build and save index
    builder.build_index();
    
    std::string index_file = output_dir + "/kmer_index.bin";
    builder.save_index(index_file);
    builder.save_metadata(output_dir);
    
    std::cout << "\nK-mer index built successfully!" << std::endl;
    std::cout << "You can now use this index to build a bloom filter:" << std::endl;
    std::cout << "  ./build_bloom_filter_simple " << index_file << " bloom_filter.bin" << std::endl;
    
    return 0;
}