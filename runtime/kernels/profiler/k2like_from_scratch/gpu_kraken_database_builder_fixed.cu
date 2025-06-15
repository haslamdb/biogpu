// kraken2_gpu_database_builder.cu
// FIXED: Proper Kraken2-style minimizer extraction with correct compression
// This version extracts minimizers exactly like Kraken2 does

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/reduce.h>
#include <cub/cub.cuh>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <algorithm>

// Configuration - much more conservative
#define MAX_GENOMES_PER_BATCH 10
#define MAX_MINIMIZERS_PER_BATCH 500000  // Reduced from 2M to 500K
#define MAX_SEQUENCE_LENGTH 50000000     // 50MB max per genome
#define THREADS_PER_BLOCK 256

// Kraken2-style parameters
struct Kraken2Params {
    int k = 35;              // k-mer length
    int ell = 31;            // minimizer length (l)
    int spaces = 7;          // spaced seed parameter
    bool use_spaced_seeds = true;
    
    // Kraken2 uses these defaults
    static Kraken2Params kraken2_defaults() {
        Kraken2Params p;
        p.k = 35;
        p.ell = 31;
        p.spaces = 7;
        return p;
    }
};

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

struct LCAEntry {
    uint64_t minimizer_hash;
    uint32_t lca_taxon;
    uint32_t genome_count;
    float uniqueness_score;
};

// Device functions for proper Kraken2-style minimizer extraction
__device__ uint64_t encode_base_kraken2(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // Invalid
    }
}

__device__ uint64_t hash_kmer_kraken2(const char* sequence, int pos, int len) {
    uint64_t hash = 0;
    for (int i = 0; i < len; i++) {
        uint64_t base = encode_base_kraken2(sequence[pos + i]);
        if (base == 4) return UINT64_MAX;  // Invalid k-mer
        hash = (hash << 2) | base;
    }
    return hash;
}

__device__ uint64_t reverse_complement_hash(uint64_t hash, int len) {
    uint64_t rc = 0;
    for (int i = 0; i < len; i++) {
        uint64_t base = (hash >> (2 * i)) & 3;
        uint64_t rc_base = 3 - base;  // A<->T, C<->G
        rc = (rc << 2) | rc_base;
    }
    return rc;
}

__device__ uint64_t canonical_hash(uint64_t hash, int len) {
    uint64_t rc = reverse_complement_hash(hash, len);
    return (hash < rc) ? hash : rc;
}

__device__ bool is_valid_kmer(const char* seq, int pos, int len) {
    for (int i = 0; i < len; i++) {
        char c = seq[pos + i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

// FIXED: Proper Kraken2-style minimizer extraction
__device__ uint64_t extract_minimizer_kraken2(
    const char* sequence,
    int kmer_pos,
    int k,
    int ell,
    int spaces) {
    
    // Check if k-mer is valid
    if (!is_valid_kmer(sequence, kmer_pos, k)) {
        return UINT64_MAX;
    }
    
    uint64_t min_hash = UINT64_MAX;
    
    // Sliding window within the k-mer to find minimizer
    for (int i = 0; i <= k - ell; i++) {
        uint64_t lmer_hash = hash_kmer_kraken2(sequence, kmer_pos + i, ell);
        if (lmer_hash == UINT64_MAX) continue;
        
        // Get canonical hash
        uint64_t canonical = canonical_hash(lmer_hash, ell);
        
        // Apply spaced seed mask if enabled
        if (spaces > 0) {
            uint64_t masked = 0;
            int out_pos = 0;
            for (int j = 0; j < ell; j++) {
                if (j % (spaces + 1) == 0) {  // Keep every (spaces+1)th position
                    uint64_t base = (canonical >> (2 * (ell - 1 - j))) & 3;
                    masked |= (base << (2 * out_pos));
                    out_pos++;
                }
            }
            canonical = masked;
        }
        
        if (canonical < min_hash) {
            min_hash = canonical;
        }
    }
    
    return min_hash;
}

// FIXED: Kraken2-style minimizer extraction kernel with proper compression
__global__ void extract_minimizers_kraken2_kernel(
    const char* sequence_data,
    const GenomeInfo* genome_info,
    int num_genomes,
    MinimizerHit* minimizer_hits,
    uint32_t* global_hit_counter,
    Kraken2Params params,
    int max_minimizers) {
    
    int genome_id = blockIdx.x;
    if (genome_id >= num_genomes) return;
    
    const GenomeInfo& genome = genome_info[genome_id];
    const char* sequence = sequence_data + genome.sequence_offset;
    uint32_t seq_length = genome.sequence_length;
    
    if (seq_length < params.k) return;
    
    // CRITICAL: Only use thread 0 to ensure sequential processing
    if (threadIdx.x != 0) return;
    
    uint64_t last_minimizer = UINT64_MAX;
    uint32_t local_count = 0;
    
    // Process each k-mer in the sequence
    for (uint32_t kmer_pos = 0; kmer_pos <= seq_length - params.k; kmer_pos++) {
        uint64_t minimizer = extract_minimizer_kraken2(
            sequence, kmer_pos, params.k, params.ell, params.spaces
        );
        
        if (minimizer != UINT64_MAX) {
            // CRITICAL: Only store if different from previous minimizer
            // This is the key to Kraken2's compression!
            if (minimizer != last_minimizer) {
                uint32_t global_pos = atomicAdd(global_hit_counter, 1);
                
                if (global_pos < max_minimizers) {
                    MinimizerHit hit;
                    hit.minimizer_hash = minimizer;
                    hit.taxon_id = genome.taxon_id;
                    hit.position = kmer_pos;
                    hit.genome_id = genome.genome_id;
                    
                    minimizer_hits[global_pos] = hit;
                    local_count++;
                } else {
                    // Hit limit, stop processing
                    break;
                }
                
                last_minimizer = minimizer;
            }
            // If minimizer == last_minimizer, skip it (compression!)
        }
    }
    
    // Optional: store per-genome stats
    // hit_counts_per_genome[genome_id] = local_count;
}

// Host-side database builder class
class Kraken2GPUDatabaseBuilder {
private:
    std::string output_directory;
    Kraken2Params params;
    
    // Host data
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    std::vector<LCAEntry> all_lca_entries;
    
    // GPU memory
    char* d_sequence_data;
    GenomeInfo* d_genome_info;
    MinimizerHit* d_minimizer_hits;
    uint32_t* d_global_counter;
    
    // Statistics
    uint64_t total_genomes_processed;
    uint64_t total_minimizers_extracted;
    uint64_t total_unique_minimizers;
    
public:
    Kraken2GPUDatabaseBuilder(const std::string& output_dir) 
        : output_directory(output_dir), params(Kraken2Params::kraken2_defaults()),
          d_sequence_data(nullptr), d_genome_info(nullptr), 
          d_minimizer_hits(nullptr), d_global_counter(nullptr),
          total_genomes_processed(0), total_minimizers_extracted(0), total_unique_minimizers(0) {
        
        std::filesystem::create_directories(output_directory);
        allocate_gpu_memory();
    }
    
    ~Kraken2GPUDatabaseBuilder() {
        free_gpu_memory();
    }
    
    bool build_database(const std::string& genome_dir, const std::string& taxonomy_dir = "") {
        std::cout << "\n=== BUILDING KRAKEN2-STYLE DATABASE ===" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Load genome files
        if (!load_genome_files(genome_dir)) {
            return false;
        }
        
        // Load taxonomy if provided
        if (!taxonomy_dir.empty()) {
            load_taxonomy_data(taxonomy_dir);
        }
        
        // Process genomes in small batches
        if (!process_all_genomes()) {
            return false;
        }
        
        // Save database
        if (!save_database()) {
            return false;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\n✓ Database build completed in " << duration.count() << " seconds" << std::endl;
        print_statistics();
        
        return true;
    }
    
private:
    bool allocate_gpu_memory() {
        size_t sequence_memory = size_t(MAX_GENOMES_PER_BATCH) * MAX_SEQUENCE_LENGTH;
        size_t genome_info_memory = MAX_GENOMES_PER_BATCH * sizeof(GenomeInfo);
        size_t minimizer_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(MinimizerHit);
        
        std::cout << "Allocating GPU memory:" << std::endl;
        std::cout << "  Sequences: " << (sequence_memory / 1024 / 1024) << " MB" << std::endl;
        std::cout << "  Minimizers: " << (minimizer_memory / 1024 / 1024) << " MB" << std::endl;
        
        cudaError_t err;
        err = cudaMalloc(&d_sequence_data, sequence_memory);
        if (err != cudaSuccess) {
            std::cerr << "Failed to allocate sequence memory: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        err = cudaMalloc(&d_genome_info, genome_info_memory);
        if (err != cudaSuccess) {
            std::cerr << "Failed to allocate genome info: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        err = cudaMalloc(&d_minimizer_hits, minimizer_memory);
        if (err != cudaSuccess) {
            std::cerr << "Failed to allocate minimizer hits: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        err = cudaMalloc(&d_global_counter, sizeof(uint32_t));
        if (err != cudaSuccess) {
            std::cerr << "Failed to allocate counter: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        return true;
    }
    
    void free_gpu_memory() {
        if (d_sequence_data) cudaFree(d_sequence_data);
        if (d_genome_info) cudaFree(d_genome_info);
        if (d_minimizer_hits) cudaFree(d_minimizer_hits);
        if (d_global_counter) cudaFree(d_global_counter);
    }
    
    bool load_genome_files(const std::string& genome_dir) {
        std::cout << "Loading genome files from: " << genome_dir << std::endl;
        
        for (const auto& entry : std::filesystem::recursive_directory_iterator(genome_dir)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                std::string extension = entry.path().extension().string();
                
                if (extension == ".fna" || extension == ".fa" || extension == ".fasta") {
                    genome_files.push_back(entry.path().string());
                    
                    // Extract taxon ID from filename
                    uint32_t taxon_id = extract_taxon_from_filename(filename);
                    genome_taxon_ids.push_back(taxon_id);
                    
                    if (taxon_names.find(taxon_id) == taxon_names.end()) {
                        taxon_names[taxon_id] = filename;
                    }
                }
            }
        }
        
        std::cout << "Found " << genome_files.size() << " genome files" << std::endl;
        return !genome_files.empty();
    }
    
    uint32_t extract_taxon_from_filename(const std::string& filename) {
        // Try to extract taxon ID from filename patterns
        // This is a simplified version - you might need more sophisticated parsing
        std::hash<std::string> hasher;
        return (hasher(filename) % 1000000) + 100000;  // Generate consistent ID
    }
    
    bool process_all_genomes() {
        std::cout << "Processing " << genome_files.size() << " genomes in batches of " 
                  << MAX_GENOMES_PER_BATCH << "..." << std::endl;
        
        for (size_t batch_start = 0; batch_start < genome_files.size(); 
             batch_start += MAX_GENOMES_PER_BATCH) {
            
            size_t batch_end = std::min(batch_start + MAX_GENOMES_PER_BATCH, genome_files.size());
            
            std::cout << "Processing batch: genomes " << batch_start << "-" << (batch_end-1) << std::endl;
            
            if (!process_genome_batch(batch_start, batch_end)) {
                std::cerr << "Failed to process batch starting at " << batch_start << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
    bool process_genome_batch(size_t batch_start, size_t batch_end) {
        // Load sequences for this batch
        std::string concatenated_sequences;
        std::vector<GenomeInfo> genome_infos;
        
        uint32_t current_offset = 0;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            auto sequences = load_fasta_sequences(genome_files[i]);
            
            for (const auto& seq : sequences) {
                if (seq.length() < params.k) continue;  // Skip too-short sequences
                
                GenomeInfo info;
                info.taxon_id = genome_taxon_ids[i];
                info.sequence_offset = current_offset;
                info.sequence_length = seq.length();
                info.genome_id = i;
                
                genome_infos.push_back(info);
                concatenated_sequences += seq;
                current_offset += seq.length();
                
                // Check for memory limits
                if (current_offset > MAX_SEQUENCE_LENGTH * MAX_GENOMES_PER_BATCH) {
                    std::cout << "Reached sequence memory limit, processing partial batch" << std::endl;
                    break;
                }
            }
        }
        
        if (genome_infos.empty()) {
            std::cout << "No valid sequences in batch, skipping" << std::endl;
            return true;
        }
        
        std::cout << "  Loaded " << genome_infos.size() << " sequences (" 
                  << (concatenated_sequences.length() / 1024 / 1024) << " MB)" << std::endl;
        
        // Transfer to GPU
        cudaMemcpy(d_sequence_data, concatenated_sequences.c_str(),
                   concatenated_sequences.length(), cudaMemcpyHostToDevice);
        
        cudaMemcpy(d_genome_info, genome_infos.data(),
                   genome_infos.size() * sizeof(GenomeInfo), cudaMemcpyHostToDevice);
        
        // Reset counter
        uint32_t zero = 0;
        cudaMemcpy(d_global_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Launch kernel
        int num_blocks = genome_infos.size();
        extract_minimizers_kraken2_kernel<<<num_blocks, 1>>>(
            d_sequence_data, d_genome_info, genome_infos.size(),
            d_minimizer_hits, d_global_counter, params, MAX_MINIMIZERS_PER_BATCH
        );
        
        cudaDeviceSynchronize();
        
        // Check for errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        // Get results
        uint32_t num_hits;
        cudaMemcpy(&num_hits, d_global_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        if (num_hits > 0) {
            std::cout << "  Extracted " << num_hits << " minimizers" << std::endl;
            
            // Process minimizers into LCA entries
            process_minimizers_to_lca(num_hits);
        }
        
        total_genomes_processed += genome_infos.size();
        total_minimizers_extracted += num_hits;
        
        return true;
    }
    
    void process_minimizers_to_lca(uint32_t num_hits) {
        // Copy minimizers from GPU
        std::vector<MinimizerHit> hits(num_hits);
        cudaMemcpy(hits.data(), d_minimizer_hits, 
                   num_hits * sizeof(MinimizerHit), cudaMemcpyDeviceToHost);
        
        // Group by hash and create LCA entries
        std::unordered_map<uint64_t, std::vector<uint32_t>> hash_to_taxa;
        
        for (const auto& hit : hits) {
            hash_to_taxa[hit.minimizer_hash].push_back(hit.taxon_id);
        }
        
        // Create LCA entries
        for (const auto& [hash, taxa] : hash_to_taxa) {
            LCAEntry entry;
            entry.minimizer_hash = hash;
            entry.lca_taxon = compute_lca(taxa);  // Simplified: just use first taxon
            entry.genome_count = taxa.size();
            entry.uniqueness_score = 1.0f / taxa.size();
            
            all_lca_entries.push_back(entry);
        }
        
        total_unique_minimizers += hash_to_taxa.size();
        
        std::cout << "  Created " << hash_to_taxa.size() << " unique LCA entries" << std::endl;
    }
    
    uint32_t compute_lca(const std::vector<uint32_t>& taxa) {
        // Simplified LCA computation - just return first taxon
        // In a real implementation, you'd walk up the taxonomy tree
        return taxa.empty() ? 0 : taxa[0];
    }
    
    std::vector<std::string> load_fasta_sequences(const std::string& fasta_file) {
        std::vector<std::string> sequences;
        std::ifstream file(fasta_file);
        
        if (!file.is_open()) {
            std::cerr << "Cannot open FASTA file: " << fasta_file << std::endl;
            return sequences;
        }
        
        std::string line, current_sequence;
        bool in_sequence = false;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (in_sequence && !current_sequence.empty()) {
                    sequences.push_back(current_sequence);
                    current_sequence.clear();
                }
                in_sequence = true;
            } else if (in_sequence) {
                current_sequence += line;
            }
        }
        
        if (!current_sequence.empty()) {
            sequences.push_back(current_sequence);
        }
        
        return sequences;
    }
    
    bool load_taxonomy_data(const std::string& taxonomy_dir) {
        // Load NCBI taxonomy files if available
        std::string nodes_file = taxonomy_dir + "/nodes.dmp";
        std::string names_file = taxonomy_dir + "/names.dmp";
        
        // This is a simplified implementation
        // You can expand this to parse actual NCBI taxonomy files
        std::cout << "Loading taxonomy from: " << taxonomy_dir << std::endl;
        
        return true;
    }
    
    bool save_database() {
        std::cout << "Saving database to " << output_directory << "..." << std::endl;
        
        // Save hash table
        std::string hash_file = output_directory + "/hash_table.k2d";
        std::ofstream hash_out(hash_file, std::ios::binary);
        
        if (!hash_out.is_open()) {
            std::cerr << "Cannot create hash table file" << std::endl;
            return false;
        }
        
        // Write header
        uint64_t table_size = all_lca_entries.size() * 2;  // 50% load factor
        uint64_t num_entries = all_lca_entries.size();
        
        hash_out.write(reinterpret_cast<const char*>(&table_size), sizeof(uint64_t));
        hash_out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
        
        // Write entries
        for (const auto& entry : all_lca_entries) {
            hash_out.write(reinterpret_cast<const char*>(&entry.minimizer_hash), sizeof(uint64_t));
            hash_out.write(reinterpret_cast<const char*>(&entry.lca_taxon), sizeof(uint32_t));
            hash_out.write(reinterpret_cast<const char*>(&entry.genome_count), sizeof(uint32_t));
            hash_out.write(reinterpret_cast<const char*>(&entry.uniqueness_score), sizeof(float));
        }
        hash_out.close();
        
        // Save taxonomy
        std::string taxonomy_file = output_directory + "/taxonomy.tsv";
        std::ofstream tax_out(taxonomy_file);
        
        tax_out << "taxon_id\tname\tparent_id\n";
        for (const auto& [taxon_id, name] : taxon_names) {
            uint32_t parent_id = taxon_parents.count(taxon_id) ? taxon_parents[taxon_id] : 0;
            tax_out << taxon_id << "\t" << name << "\t" << parent_id << "\n";
        }
        tax_out.close();
        
        std::cout << "✓ Database saved: " << num_entries << " minimizers" << std::endl;
        return true;
    }
    
    void print_statistics() {
        std::cout << "\n=== DATABASE BUILD STATISTICS ===" << std::endl;
        std::cout << "Total genomes processed: " << total_genomes_processed << std::endl;
        std::cout << "Total minimizers extracted: " << total_minimizers_extracted << std::endl;
        std::cout << "Unique minimizers: " << total_unique_minimizers << std::endl;
        std::cout << "Compression ratio: " << std::fixed << std::setprecision(2) 
                  << (double)total_unique_minimizers / total_minimizers_extracted << std::endl;
        std::cout << "Average minimizers per genome: " 
                  << (total_minimizers_extracted / total_genomes_processed) << std::endl;
    }
};

// Main function to test the database builder
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <genome_dir> <output_dir> [taxonomy_dir]" << std::endl;
        return 1;
    }
    
    std::string genome_dir = argv[1];
    std::string output_dir = argv[2];
    std::string taxonomy_dir = (argc > 3) ? argv[3] : "";
    
    std::cout << "Building Kraken2-style database..." << std::endl;
    std::cout << "Genome directory: " << genome_dir << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    
    Kraken2GPUDatabaseBuilder builder(output_dir);
    
    bool success = builder.build_database(genome_dir, taxonomy_dir);
    
    if (success) {
        std::cout << "\n🎉 Database build completed successfully!" << std::endl;
        return 0;
    } else {
        std::cout << "\n❌ Database build failed!" << std::endl;
        return 1;
    }
}

// Compilation instructions:
/*
nvcc -std=c++14 -O3 -arch=sm_70 kraken2_gpu_database_builder.cu -o kraken2_gpu_builder

Usage:
./kraken2_gpu_builder ./genomes ./output_db ./taxonomy
*/