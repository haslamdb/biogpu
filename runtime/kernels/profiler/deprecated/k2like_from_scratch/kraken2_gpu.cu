// kraken2_gpu.cu
// Clean implementation of Kraken2-style database building and classification
// Fixed minimizer extraction to match Kraken2 behavior

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <iomanip>

// Kraken2 parameters - matching default values
struct Kraken2Config {
    int k = 35;              // k-mer length
    int l = 31;              // minimizer length
    int s = 7;               // spaced seed parameter
    float confidence = 0.0f; // confidence threshold
    
    // Memory limits
    int max_genomes_per_batch = 50;
    int max_minimizers_per_batch = 1000000;  // 1M minimizers max
};

// Structures
struct GenomeInfo {
    uint32_t taxon_id;
    uint32_t sequence_offset;
    uint32_t sequence_length;
    uint32_t genome_id;
};

struct MinimizerHit {
    uint64_t minimizer_hash;
    uint32_t taxon_id;
    uint32_t position;
    uint32_t genome_id;
};

struct HashTableEntry {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
};

struct ClassificationResult {
    uint32_t taxon_id;
    float confidence;
    uint32_t hit_count;
    uint32_t total_kmers;
};

// Device functions for minimizer extraction
__device__ uint64_t encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1; 
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // Invalid
    }
}

__device__ uint64_t hash_sequence(const char* seq, int pos, int len) {
    uint64_t hash = 0;
    for (int i = 0; i < len; i++) {
        uint64_t base = encode_base(seq[pos + i]);
        if (base == 4) return UINT64_MAX;  // Invalid sequence
        hash = (hash << 2) | base;
    }
    return hash;
}

__device__ uint64_t reverse_complement(uint64_t hash, int len) {
    uint64_t rc = 0;
    for (int i = 0; i < len; i++) {
        uint64_t base = (hash >> (2 * i)) & 3;
        uint64_t rc_base = 3 - base;  // A<->T, C<->G
        rc = (rc << 2) | rc_base;
    }
    return rc;
}

__device__ uint64_t canonical_hash(uint64_t hash, int len) {
    uint64_t rc = reverse_complement(hash, len);
    return (hash < rc) ? hash : rc;
}

__device__ bool is_valid_sequence(const char* seq, int pos, int len) {
    for (int i = 0; i < len; i++) {
        char c = seq[pos + i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

// Kraken2-style minimizer extraction
__device__ uint64_t extract_minimizer(const char* sequence, int kmer_pos, 
                                      int k, int l, int s) {
    if (!is_valid_sequence(sequence, kmer_pos, k)) {
        return UINT64_MAX;
    }
    
    uint64_t min_hash = UINT64_MAX;
    
    // Sliding window within k-mer to find minimizer
    for (int i = 0; i <= k - l; i++) {
        uint64_t lmer_hash = hash_sequence(sequence, kmer_pos + i, l);
        if (lmer_hash == UINT64_MAX) continue;
        
        uint64_t canonical = canonical_hash(lmer_hash, l);
        
        // Apply spaced seed if enabled
        if (s > 0) {
            uint64_t spaced = 0;
            int out_pos = 0;
            for (int j = 0; j < l; j++) {
                if (j % (s + 1) == 0) {
                    uint64_t base = (canonical >> (2 * (l - 1 - j))) & 3;
                    spaced |= (base << (2 * out_pos));
                    out_pos++;
                }
            }
            canonical = spaced;
        }
        
        if (canonical < min_hash) {
            min_hash = canonical;
        }
    }
    
    return min_hash;
}

// Database building kernel
__global__ void extract_minimizers_kernel(
    const char* sequence_data,
    const GenomeInfo* genome_info,
    int num_genomes,
    MinimizerHit* minimizer_hits,
    uint32_t* global_counter,
    Kraken2Config config,
    int max_hits) {
    
    int genome_idx = blockIdx.x;
    if (genome_idx >= num_genomes || threadIdx.x != 0) return;
    
    const GenomeInfo& genome = genome_info[genome_idx];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_len = genome.sequence_length;
    
    if (seq_len < config.k) return;
    
    uint64_t last_minimizer = UINT64_MAX;
    uint32_t total_kmers = seq_len - config.k + 1;
    
    // Process each k-mer
    for (uint32_t pos = 0; pos < total_kmers; pos++) {
        uint64_t minimizer = extract_minimizer(sequence, pos, config.k, config.l, config.s);
        
        if (minimizer != UINT64_MAX) {
            // Only store if different from last minimizer (Kraken2 compression)
            if (minimizer != last_minimizer) {
                uint32_t hit_idx = atomicAdd(global_counter, 1);
                
                if (hit_idx < max_hits) {
                    MinimizerHit hit;
                    hit.minimizer_hash = minimizer;
                    hit.taxon_id = genome.taxon_id;
                    hit.position = pos;
                    hit.genome_id = genome.genome_id;
                    minimizer_hits[hit_idx] = hit;
                } else {
                    break;  // Hit limit
                }
                
                last_minimizer = minimizer;
            }
        }
    }
}

// Classification kernel
__global__ void classify_reads_kernel(
    const char* reads_data,
    const uint32_t* read_offsets,
    const uint32_t* read_lengths,
    int num_reads,
    const HashTableEntry* hash_table,
    uint32_t table_size,
    ClassificationResult* results,
    Kraken2Config config) {
    
    int read_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (read_idx >= num_reads) return;
    
    const char* read = reads_data + read_offsets[read_idx];
    uint32_t read_len = read_lengths[read_idx];
    
    ClassificationResult& result = results[read_idx];
    result.taxon_id = 0;
    result.confidence = 0.0f;
    result.hit_count = 0;
    result.total_kmers = 0;
    
    if (read_len < config.k) return;
    
    // Count hits for each taxon
    uint32_t taxon_votes[64] = {0};  // Support up to 64 taxa per read
    uint32_t taxon_ids[64] = {0};
    int num_taxa = 0;
    
    uint32_t total_kmers = read_len - config.k + 1;
    uint64_t last_minimizer = UINT64_MAX;
    
    for (uint32_t pos = 0; pos < total_kmers; pos++) {
        uint64_t minimizer = extract_minimizer(read, pos, config.k, config.l, config.s);
        
        if (minimizer != UINT64_MAX && minimizer != last_minimizer) {
            result.total_kmers++;
            
            // Look up in hash table
            uint32_t hash_pos = minimizer % table_size;
            
            // Linear probing
            for (int probe = 0; probe < 32; probe++) {
                uint32_t idx = (hash_pos + probe) % table_size;
                const HashTableEntry& entry = hash_table[idx];
                
                if (entry.minimizer_hash == 0) break;  // Empty slot
                
                if (entry.minimizer_hash == minimizer) {
                    // Found match - add vote
                    uint32_t taxon = entry.lca_taxon;
                    
                    // Find or add taxon in local array
                    bool found = false;
                    for (int i = 0; i < num_taxa; i++) {
                        if (taxon_ids[i] == taxon) {
                            taxon_votes[i]++;
                            found = true;
                            break;
                        }
                    }
                    
                    if (!found && num_taxa < 64) {
                        taxon_ids[num_taxa] = taxon;
                        taxon_votes[num_taxa] = 1;
                        num_taxa++;
                    }
                    
                    result.hit_count++;
                    break;
                }
            }
            
            last_minimizer = minimizer;
        }
    }
    
    // Find best taxon
    uint32_t best_votes = 0;
    uint32_t best_taxon = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (taxon_votes[i] > best_votes) {
            best_votes = taxon_votes[i];
            best_taxon = taxon_ids[i];
        }
    }
    
    if (best_votes > 0) {
        result.taxon_id = best_taxon;
        result.confidence = (float)best_votes / result.total_kmers;
        
        // Apply confidence threshold
        if (result.confidence < config.confidence) {
            result.taxon_id = 0;  // Unclassified
        }
    }
}

// Host classes
class Kraken2GPUDatabase {
private:
    std::vector<HashTableEntry> hash_table;
    std::unordered_map<uint32_t, std::string> taxon_names;
    uint32_t table_size;
    Kraken2Config config;
    
public:
    Kraken2GPUDatabase(const Kraken2Config& cfg = Kraken2Config()) : config(cfg) {}
    
    bool build_from_genomes(const std::string& genome_dir, const std::string& output_dir) {
        std::cout << "\n=== Building Kraken2-style Database ===" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Find genome files
        std::vector<std::string> genome_files = find_genome_files(genome_dir);
        if (genome_files.empty()) {
            std::cerr << "No genome files found in " << genome_dir << std::endl;
            return false;
        }
        
        std::cout << "Found " << genome_files.size() << " genome files" << std::endl;
        
        // Process genomes in batches
        std::vector<MinimizerHit> all_hits;
        
        for (size_t batch_start = 0; batch_start < genome_files.size(); 
             batch_start += config.max_genomes_per_batch) {
            
            size_t batch_end = std::min(batch_start + config.max_genomes_per_batch, 
                                       genome_files.size());
            
            std::cout << "Processing batch " << (batch_start / config.max_genomes_per_batch + 1) 
                      << ": genomes " << batch_start << "-" << (batch_end-1) << std::endl;
            
            auto batch_hits = process_genome_batch(genome_files, batch_start, batch_end);
            all_hits.insert(all_hits.end(), batch_hits.begin(), batch_hits.end());
            
            std::cout << "  Total hits so far: " << all_hits.size() << std::endl;
        }
        
        // Build hash table
        build_hash_table(all_hits);
        
        // Save database
        std::filesystem::create_directories(output_dir);
        save_database(output_dir);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "✓ Database built in " << duration.count() << " seconds" << std::endl;
        std::cout << "  Total minimizers: " << all_hits.size() << std::endl;
        std::cout << "  Hash table size: " << table_size << std::endl;
        std::cout << "  Unique minimizers: " << hash_table.size() << std::endl;
        
        return true;
    }
    
    bool load_database(const std::string& db_dir) {
        std::string hash_file = db_dir + "/hash_table.k2d";
        std::string tax_file = db_dir + "/taxonomy.tsv";
        
        if (!load_hash_table(hash_file) || !load_taxonomy(tax_file)) {
            return false;
        }
        
        std::cout << "Database loaded: " << hash_table.size() << " entries" << std::endl;
        return true;
    }
    
    std::vector<ClassificationResult> classify_reads(const std::vector<std::string>& reads) {
        if (hash_table.empty()) {
            std::cerr << "Database not loaded!" << std::endl;
            return {};
        }
        
        std::cout << "Classifying " << reads.size() << " reads..." << std::endl;
        
        // Prepare data for GPU
        std::string concatenated_reads;
        std::vector<uint32_t> read_offsets, read_lengths;
        
        uint32_t offset = 0;
        for (const auto& read : reads) {
            read_offsets.push_back(offset);
            read_lengths.push_back(read.length());
            concatenated_reads += read;
            offset += read.length();
        }
        
        // Allocate GPU memory
        char* d_reads;
        uint32_t* d_offsets;
        uint32_t* d_lengths;
        HashTableEntry* d_hash_table;
        ClassificationResult* d_results;
        
        cudaMalloc(&d_reads, concatenated_reads.length());
        cudaMalloc(&d_offsets, reads.size() * sizeof(uint32_t));
        cudaMalloc(&d_lengths, reads.size() * sizeof(uint32_t));
        cudaMalloc(&d_hash_table, hash_table.size() * sizeof(HashTableEntry));
        cudaMalloc(&d_results, reads.size() * sizeof(ClassificationResult));
        
        // Copy data to GPU
        cudaMemcpy(d_reads, concatenated_reads.c_str(), concatenated_reads.length(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_offsets, read_offsets.data(), read_offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_lengths, read_lengths.data(), read_lengths.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_hash_table, hash_table.data(), hash_table.size() * sizeof(HashTableEntry), cudaMemcpyHostToDevice);
        
        // Launch classification kernel
        int num_blocks = (reads.size() + 255) / 256;
        classify_reads_kernel<<<num_blocks, 256>>>(
            d_reads, d_offsets, d_lengths, reads.size(),
            d_hash_table, table_size, d_results, config
        );
        
        cudaDeviceSynchronize();
        
        // Get results
        std::vector<ClassificationResult> results(reads.size());
        cudaMemcpy(results.data(), d_results, reads.size() * sizeof(ClassificationResult), cudaMemcpyDeviceToHost);
        
        // Cleanup
        cudaFree(d_reads);
        cudaFree(d_offsets);
        cudaFree(d_lengths);
        cudaFree(d_hash_table);
        cudaFree(d_results);
        
        return results;
    }
    
private:
    std::vector<std::string> find_genome_files(const std::string& dir) {
        std::vector<std::string> files;
        
        for (const auto& entry : std::filesystem::recursive_directory_iterator(dir)) {
            if (entry.is_regular_file()) {
                std::string ext = entry.path().extension().string();
                if (ext == ".fna" || ext == ".fa" || ext == ".fasta") {
                    files.push_back(entry.path().string());
                }
            }
        }
        
        std::sort(files.begin(), files.end());
        return files;
    }
    
    std::vector<MinimizerHit> process_genome_batch(const std::vector<std::string>& genome_files,
                                                   size_t start, size_t end) {
        // Load sequences
        std::string concatenated_seqs;
        std::vector<GenomeInfo> genome_infos;
        
        uint32_t offset = 0;
        for (size_t i = start; i < end; i++) {
            auto sequences = load_fasta(genome_files[i]);
            uint32_t taxon_id = extract_taxon_id(genome_files[i]);
            
            if (taxon_names.find(taxon_id) == taxon_names.end()) {
                taxon_names[taxon_id] = std::filesystem::path(genome_files[i]).stem().string();
            }
            
            for (const auto& seq : sequences) {
                if (seq.length() >= config.k) {
                    GenomeInfo info;
                    info.taxon_id = taxon_id;
                    info.sequence_offset = offset;
                    info.sequence_length = seq.length();
                    info.genome_id = i;
                    
                    genome_infos.push_back(info);
                    concatenated_seqs += seq;
                    offset += seq.length();
                }
            }
        }
        
        if (genome_infos.empty()) {
            return {};
        }
        
        // Allocate GPU memory
        char* d_sequences;
        GenomeInfo* d_genome_info;
        MinimizerHit* d_hits;
        uint32_t* d_counter;
        
        cudaMalloc(&d_sequences, concatenated_seqs.length());
        cudaMalloc(&d_genome_info, genome_infos.size() * sizeof(GenomeInfo));
        cudaMalloc(&d_hits, config.max_minimizers_per_batch * sizeof(MinimizerHit));
        cudaMalloc(&d_counter, sizeof(uint32_t));
        
        // Copy data
        cudaMemcpy(d_sequences, concatenated_seqs.c_str(), concatenated_seqs.length(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_genome_info, genome_infos.data(), genome_infos.size() * sizeof(GenomeInfo), cudaMemcpyHostToDevice);
        
        uint32_t zero = 0;
        cudaMemcpy(d_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Launch kernel
        extract_minimizers_kernel<<<genome_infos.size(), 1>>>(
            d_sequences, d_genome_info, genome_infos.size(),
            d_hits, d_counter, config, config.max_minimizers_per_batch
        );
        
        cudaDeviceSynchronize();
        
        // Get results
        uint32_t num_hits;
        cudaMemcpy(&num_hits, d_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::vector<MinimizerHit> hits(num_hits);
        if (num_hits > 0) {
            cudaMemcpy(hits.data(), d_hits, num_hits * sizeof(MinimizerHit), cudaMemcpyDeviceToHost);
        }
        
        // Cleanup
        cudaFree(d_sequences);
        cudaFree(d_genome_info);
        cudaFree(d_hits);
        cudaFree(d_counter);
        
        std::cout << "  Extracted " << num_hits << " minimizers from " 
                  << genome_infos.size() << " sequences" << std::endl;
        
        return hits;
    }
    
    void build_hash_table(const std::vector<MinimizerHit>& hits) {
        std::cout << "Building hash table from " << hits.size() << " minimizers..." << std::endl;
        
        // Group by minimizer hash
        std::unordered_map<uint64_t, std::vector<uint32_t>> hash_to_taxa;
        
        for (const auto& hit : hits) {
            hash_to_taxa[hit.minimizer_hash].push_back(hit.taxon_id);
        }
        
        // Create hash table entries
        hash_table.clear();
        hash_table.reserve(hash_to_taxa.size());
        
        for (const auto& [hash, taxa] : hash_to_taxa) {
            HashTableEntry entry;
            entry.minimizer_hash = hash;
            entry.lca_taxon = compute_lca(taxa);
            entry.genome_count = taxa.size();
            entry.uniqueness_score = 1.0f / taxa.size();
            hash_table.push_back(entry);
        }
        
        // Set table size for hash table lookup
        table_size = next_power_of_2(hash_table.size() * 2);  // 50% load factor
        
        std::cout << "Hash table built: " << hash_table.size() << " unique minimizers" << std::endl;
    }
    
    uint32_t compute_lca(const std::vector<uint32_t>& taxa) {
        // Simplified LCA - just return first taxon
        // In a real implementation, walk up taxonomy tree
        return taxa.empty() ? 0 : taxa[0];
    }
    
    uint32_t next_power_of_2(uint32_t n) {
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;
        return n;
    }
    
    std::vector<std::string> load_fasta(const std::string& filename) {
        std::vector<std::string> sequences;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return sequences;
        }
        
        std::string line, current_seq;
        bool in_seq = false;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (in_seq && !current_seq.empty()) {
                    sequences.push_back(current_seq);
                    current_seq.clear();
                }
                in_seq = true;
            } else if (in_seq) {
                current_seq += line;
            }
        }
        
        if (!current_seq.empty()) {
            sequences.push_back(current_seq);
        }
        
        return sequences;
    }
    
    uint32_t extract_taxon_id(const std::string& filename) {
        std::hash<std::string> hasher;
        std::string basename = std::filesystem::path(filename).stem().string();
        return (hasher(basename) % 1000000) + 100000;
    }
    
    void save_database(const std::string& output_dir) {
        // Save hash table
        std::string hash_file = output_dir + "/hash_table.k2d";
        std::ofstream out(hash_file, std::ios::binary);
        
        uint64_t size = table_size;
        uint64_t entries = hash_table.size();
        
        out.write(reinterpret_cast<const char*>(&size), sizeof(uint64_t));
        out.write(reinterpret_cast<const char*>(&entries), sizeof(uint64_t));
        
        for (const auto& entry : hash_table) {
            out.write(reinterpret_cast<const char*>(&entry), sizeof(HashTableEntry));
        }
        out.close();
        
        // Save taxonomy
        std::string tax_file = output_dir + "/taxonomy.tsv";
        std::ofstream tax_out(tax_file);
        
        tax_out << "taxon_id\tname\n";
        for (const auto& [id, name] : taxon_names) {
            tax_out << id << "\t" << name << "\n";
        }
        tax_out.close();
    }
    
    bool load_hash_table(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }
        
        uint64_t size, entries;
        in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&entries), sizeof(uint64_t));
        
        table_size = size;
        hash_table.resize(entries);
        
        for (auto& entry : hash_table) {
            in.read(reinterpret_cast<char*>(&entry), sizeof(HashTableEntry));
        }
        
        return true;
    }
    
    bool load_taxonomy(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }
        
        std::string line;
        std::getline(in, line);  // Skip header
        
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            uint32_t id;
            std::string name;
            
            if (iss >> id) {
                iss.ignore(1);  // Skip tab
                std::getline(iss, name);
                taxon_names[id] = name;
            }
        }
        
        return true;
    }
};

// Simple command-line interface
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage:\n";
        std::cout << "  Build: " << argv[0] << " build <genome_dir> <output_dir>\n";
        std::cout << "  Classify: " << argv[0] << " classify <database_dir> <reads_file>\n";
        return 1;
    }
    
    std::string command = argv[1];
    
    Kraken2Config config;
    Kraken2GPUDatabase db(config);
    
    if (command == "build" && argc >= 4) {
        std::string genome_dir = argv[2];
        std::string output_dir = argv[3];
        
        if (db.build_from_genomes(genome_dir, output_dir)) {
            std::cout << "✓ Database build successful\n";
            return 0;
        } else {
            std::cout << "❌ Database build failed\n";
            return 1;
        }
        
    } else if (command == "classify" && argc >= 4) {
        std::string db_dir = argv[2];
        std::string reads_file = argv[3];
        
        if (!db.load_database(db_dir)) {
            std::cout << "❌ Failed to load database\n";
            return 1;
        }
        
        // Load reads from FASTQ
        std::vector<std::string> reads;
        std::ifstream file(reads_file);
        std::string line;
        int line_count = 0;
        
        while (std::getline(file, line) && reads.size() < 10000) {  // Limit for testing
            line_count++;
            if (line_count % 4 == 2) {  // Sequence line in FASTQ
                if (line.length() >= config.k) {
                    reads.push_back(line);
                }
            }
        }
        
        if (reads.empty()) {
            std::cout << "No valid reads found\n";
            return 1;
        }
        
        auto results = db.classify_reads(reads);
        
        // Print results
        int classified = 0;
        for (size_t i = 0; i < results.size(); i++) {
            const auto& result = results[i];
            char status = (result.taxon_id > 0) ? 'C' : 'U';
            
            std::cout << status << "\tread_" << i << "\t" << result.taxon_id 
                      << "\t" << reads[i].length() << "\t" 
                      << std::fixed << std::setprecision(3) << result.confidence 
                      << "\n";
            
            if (result.taxon_id > 0) classified++;
        }
        
        std::cout << "\nClassified: " << classified << "/" << results.size() 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * classified / results.size()) << "%)\n";
        
        return 0;
        
    } else {
        std::cout << "Invalid command or arguments\n";
        return 1;
    }
}

// Compile with:
// nvcc -std=c++17 -O3 -arch=sm_70 kraken2_gpu.cu -o kraken2_gpu
//
// Usage:
// ./kraken2_gpu build ./genomes ./database
// ./kraken2_gpu classify ./database reads.fastq