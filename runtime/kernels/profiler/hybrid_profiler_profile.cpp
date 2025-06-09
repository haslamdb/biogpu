// hybrid_kraken_pipeline.cpp
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <cuda_runtime.h>
#include <omp.h>
#include "minimizer_common.h"

// Forward declarations
class MinimizerExtractor {
public:
    MinimizerExtractor(int k_size = 31, int window_size = 15);
    ~MinimizerExtractor();
    std::vector<std::vector<Minimizer>> extract_minimizers(const std::vector<std::string>& sequences);
private:
    void allocate_device_memory(size_t num_reads, size_t total_sequence_length);
    int k;
    int m;
    size_t allocated_reads;
    void* d_sequences;
    void* d_sequence_offsets;
    void* d_sequence_lengths;
    void* d_minimizers;
    void* d_minimizer_counts;
};

class FastqReader {
public:
    FastqReader(const std::string& filename);
    ~FastqReader();
    bool read_batch(ReadBatch& batch, size_t batch_size);
private:
    std::string filename;
    void* gz_file;
    bool is_gzipped;
    std::ifstream text_file;
    bool getline(std::string& line);
};

// Kraken2 database structures
struct KrakenDBHeader {
    uint64_t version;
    uint64_t k_size;
    uint64_t minimizer_len;
    uint64_t taxonomy_nodes;
    uint64_t database_size;
};

// Memory-mapped Kraken2 database
class Kraken2Database {
private:
    // Memory mapped files
    void* taxonomy_mmap = nullptr;
    void* hashtable_mmap = nullptr;
    size_t taxonomy_size;
    size_t hashtable_size;
    
    // Parsed structures
    KrakenDBHeader header;
    uint32_t* parent_map;           // Taxonomy tree
    std::unordered_map<uint64_t, uint32_t> minimizer_to_taxon;  // For smaller DBs
    
    // For large databases, use direct memory access
    struct HashTableEntry {
        uint64_t minimizer;
        uint32_t taxon_id;
    } __attribute__((packed));
    
    HashTableEntry* hash_table = nullptr;
    size_t hash_table_entries;
    
public:
    Kraken2Database(const std::string& db_path) {
        load_taxonomy(db_path + "/taxo.k2d");
        load_hashtable(db_path + "/hash.k2d");
    }
    
    ~Kraken2Database() {
        if (taxonomy_mmap) munmap(taxonomy_mmap, taxonomy_size);
        if (hashtable_mmap) munmap(hashtable_mmap, hashtable_size);
    }
    
    void load_taxonomy(const std::string& filename) {
        int fd = open(filename.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("Cannot open taxonomy file");
        
        struct stat sb;
        fstat(fd, &sb);
        taxonomy_size = sb.st_size;
        
        taxonomy_mmap = mmap(nullptr, taxonomy_size, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        
        if (taxonomy_mmap == MAP_FAILED) {
            throw std::runtime_error("Cannot mmap taxonomy file");
        }
        
        // Parse taxonomy structure
        parent_map = (uint32_t*)((char*)taxonomy_mmap + sizeof(KrakenDBHeader));
    }
    
    void load_hashtable(const std::string& filename) {
        int fd = open(filename.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("Cannot open hash table file");
        
        struct stat sb;
        fstat(fd, &sb);
        hashtable_size = sb.st_size;
        
        hashtable_mmap = mmap(nullptr, hashtable_size, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        
        if (hashtable_mmap == MAP_FAILED) {
            throw std::runtime_error("Cannot mmap hash table file");
        }
        
        // Setup direct access to hash table
        hash_table = (HashTableEntry*)((char*)hashtable_mmap + sizeof(KrakenDBHeader));
        hash_table_entries = (hashtable_size - sizeof(KrakenDBHeader)) / sizeof(HashTableEntry);
        
        std::cout << "Loaded " << hash_table_entries << " hash table entries\n";
    }
    
    // Lookup minimizer using binary search (assuming sorted hash table)
    uint32_t lookup_minimizer(uint64_t minimizer) const {
        // Binary search in sorted array
        size_t left = 0;
        size_t right = hash_table_entries;
        
        while (left < right) {
            size_t mid = left + (right - left) / 2;
            if (hash_table[mid].minimizer == minimizer) {
                return hash_table[mid].taxon_id;
            } else if (hash_table[mid].minimizer < minimizer) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        return 0;  // Not found
    }
    
    // Batch lookup for better cache efficiency
    void lookup_batch(const uint64_t* minimizers, uint32_t* taxon_ids, size_t count) const {
        #pragma omp parallel for
        for (size_t i = 0; i < count; i++) {
            taxon_ids[i] = lookup_minimizer(minimizers[i]);
        }
    }
    
    // Find lowest common ancestor
    uint32_t find_lca(uint32_t taxon1, uint32_t taxon2) const {
        std::unordered_set<uint32_t> ancestors;
        
        // Trace taxon1's ancestry
        uint32_t current = taxon1;
        while (current != 0 && current != 1) {  // 1 is root
            ancestors.insert(current);
            current = parent_map[current];
        }
        
        // Find first common ancestor with taxon2
        current = taxon2;
        while (current != 0 && current != 1) {
            if (ancestors.count(current)) {
                return current;
            }
            current = parent_map[current];
        }
        
        return 1;  // Root
    }
};

// Result structure for classification
struct ClassificationResult {
    uint32_t taxon_id;
    float confidence;
    uint32_t hits;
};

// GPU kernel for processing classification results
__global__ void process_classification_kernel(
    const uint32_t* minimizer_taxons,      // Taxon IDs from CPU lookup
    const uint32_t* minimizer_offsets,     // Offset for each read
    const uint32_t* minimizer_counts,      // Number of minimizers per read
    ClassificationResult* results,         // Output classifications
    const uint32_t num_reads
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const uint32_t start = minimizer_offsets[tid];
    const uint32_t count = minimizer_counts[tid];
    
    if (count == 0) {
        results[tid] = {0, 0.0f, 0};
        return;
    }
    
    // Count hits per taxon
    const int MAX_TAXA = 32;  // Local storage limit
    uint32_t taxon_hits[MAX_TAXA];
    uint32_t unique_taxa[MAX_TAXA];
    int num_unique = 0;
    
    for (int i = 0; i < MAX_TAXA; i++) {
        taxon_hits[i] = 0;
        unique_taxa[i] = 0;
    }
    
    // Count occurrences
    for (uint32_t i = 0; i < count; i++) {
        uint32_t taxon = minimizer_taxons[start + i];
        if (taxon == 0) continue;  // Skip unclassified
        
        // Find or add taxon
        bool found = false;
        for (int j = 0; j < num_unique; j++) {
            if (unique_taxa[j] == taxon) {
                taxon_hits[j]++;
                found = true;
                break;
            }
        }
        
        if (!found && num_unique < MAX_TAXA) {
            unique_taxa[num_unique] = taxon;
            taxon_hits[num_unique] = 1;
            num_unique++;
        }
    }
    
    // Find taxon with most hits
    uint32_t best_taxon = 0;
    uint32_t max_hits = 0;
    
    for (int i = 0; i < num_unique; i++) {
        if (taxon_hits[i] > max_hits) {
            max_hits = taxon_hits[i];
            best_taxon = unique_taxa[i];
        }
    }
    
    // Calculate confidence
    float confidence = (count > 0) ? (float)max_hits / count : 0.0f;
    
    results[tid] = {best_taxon, confidence, max_hits};
}

// GPU kernel for abundance counting
__global__ void update_abundance_kernel(
    const ClassificationResult* classifications,
    const uint32_t num_reads,
    uint32_t* taxon_counts,        // Global counts per taxon
    const uint32_t max_taxon_id
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    uint32_t taxon = classifications[tid].taxon_id;
    if (taxon > 0 && taxon < max_taxon_id) {
        atomicAdd(&taxon_counts[taxon], 1);
    }
}

// Main hybrid pipeline
class HybridKrakenPipeline {
private:
    std::unique_ptr<Kraken2Database> database;
    std::unique_ptr<MinimizerExtractor> extractor;
    
    // GPU memory for results
    uint32_t* d_minimizer_taxons = nullptr;
    uint32_t* d_minimizer_offsets = nullptr;
    uint32_t* d_minimizer_counts = nullptr;
    ClassificationResult* d_classifications = nullptr;
    uint32_t* d_taxon_counts = nullptr;
    
    // Pinned host memory for fast transfers
    uint64_t* h_minimizer_hashes = nullptr;
    uint32_t* h_minimizer_taxons = nullptr;
    
    size_t allocated_reads = 0;
    size_t allocated_minimizers = 0;
    
public:
    HybridKrakenPipeline(const std::string& db_path, int k = 31, int m = 15) {
        database = std::make_unique<Kraken2Database>(db_path);
        extractor = std::make_unique<MinimizerExtractor>(k, m);
        
        // Allocate initial GPU memory
        allocate_gpu_memory(10000, 500000);  // 10k reads, 500k minimizers
        
        // Allocate taxon counting array (adjust size as needed)
        cudaMalloc(&d_taxon_counts, 1000000 * sizeof(uint32_t));
        cudaMemset(d_taxon_counts, 0, 1000000 * sizeof(uint32_t));
    }
    
    ~HybridKrakenPipeline() {
        if (d_minimizer_taxons) cudaFree(d_minimizer_taxons);
        if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
        if (d_minimizer_counts) cudaFree(d_minimizer_counts);
        if (d_classifications) cudaFree(d_classifications);
        if (d_taxon_counts) cudaFree(d_taxon_counts);
        if (h_minimizer_hashes) cudaFreeHost(h_minimizer_hashes);
        if (h_minimizer_taxons) cudaFreeHost(h_minimizer_taxons);
    }
    
    void allocate_gpu_memory(size_t num_reads, size_t num_minimizers) {
        if (num_reads > allocated_reads) {
            if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
            if (d_minimizer_counts) cudaFree(d_minimizer_counts);
            if (d_classifications) cudaFree(d_classifications);
            
            cudaMalloc(&d_minimizer_offsets, (num_reads + 1) * sizeof(uint32_t));
            cudaMalloc(&d_minimizer_counts, num_reads * sizeof(uint32_t));
            cudaMalloc(&d_classifications, num_reads * sizeof(ClassificationResult));
            
            allocated_reads = num_reads;
        }
        
        if (num_minimizers > allocated_minimizers) {
            if (d_minimizer_taxons) cudaFree(d_minimizer_taxons);
            if (h_minimizer_hashes) cudaFreeHost(h_minimizer_hashes);
            if (h_minimizer_taxons) cudaFreeHost(h_minimizer_taxons);
            
            cudaMalloc(&d_minimizer_taxons, num_minimizers * sizeof(uint32_t));
            cudaMallocHost(&h_minimizer_hashes, num_minimizers * sizeof(uint64_t));
            cudaMallocHost(&h_minimizer_taxons, num_minimizers * sizeof(uint32_t));
            
            allocated_minimizers = num_minimizers;
        }
    }
    
    std::vector<ClassificationResult> classify_batch(const std::vector<std::string>& sequences) {
        size_t num_reads = sequences.size();
        
        // Step 1: Extract minimizers on GPU
        auto minimizer_results = extractor->extract_minimizers(sequences);
        
        // Prepare minimizer data
        std::vector<uint64_t> all_hashes;
        std::vector<uint32_t> offsets(num_reads + 1);
        std::vector<uint32_t> counts(num_reads);
        
        uint32_t current_offset = 0;
        for (size_t i = 0; i < num_reads; i++) {
            offsets[i] = current_offset;
            counts[i] = minimizer_results[i].size();
            
            for (const auto& m : minimizer_results[i]) {
                all_hashes.push_back(m.hash);
            }
            
            current_offset += counts[i];
        }
        offsets[num_reads] = current_offset;
        
        // Ensure we have enough memory
        allocate_gpu_memory(num_reads, all_hashes.size());
        
        // Copy minimizer hashes to pinned memory
        std::copy(all_hashes.begin(), all_hashes.end(), h_minimizer_hashes);
        
        // Step 2: CPU database lookup
        auto lookup_start = std::chrono::high_resolution_clock::now();
        database->lookup_batch(h_minimizer_hashes, h_minimizer_taxons, all_hashes.size());
        auto lookup_end = std::chrono::high_resolution_clock::now();
        
        auto lookup_time = std::chrono::duration_cast<std::chrono::milliseconds>(lookup_end - lookup_start);
        std::cout << "Database lookup: " << lookup_time.count() << " ms for " 
                  << all_hashes.size() << " minimizers\n";
        
        // Step 3: Transfer results to GPU and classify
        cudaMemcpy(d_minimizer_taxons, h_minimizer_taxons, 
                   all_hashes.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_offsets, offsets.data(), 
                   offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_counts, counts.data(), 
                   counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Launch classification kernel
        int block_size = 256;
        int num_blocks = (num_reads + block_size - 1) / block_size;
        
        process_classification_kernel<<<num_blocks, block_size>>>(
            d_minimizer_taxons, d_minimizer_offsets, d_minimizer_counts,
            d_classifications, num_reads
        );
        
        // Update abundance counts
        update_abundance_kernel<<<num_blocks, block_size>>>(
            d_classifications, num_reads, d_taxon_counts, 1000000
        );
        
        cudaDeviceSynchronize();
        
        // Copy results back
        std::vector<ClassificationResult> results(num_reads);
        cudaMemcpy(results.data(), d_classifications, 
                   num_reads * sizeof(ClassificationResult), cudaMemcpyDeviceToHost);
        
        return results;
    }
    
    std::vector<std::pair<uint32_t, uint32_t>> get_abundance_profile() {
        std::vector<uint32_t> taxon_counts(1000000);
        cudaMemcpy(taxon_counts.data(), d_taxon_counts, 
                   1000000 * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::vector<std::pair<uint32_t, uint32_t>> profile;
        for (size_t i = 0; i < taxon_counts.size(); i++) {
            if (taxon_counts[i] > 0) {
                profile.push_back({i, taxon_counts[i]});
            }
        }
        
        // Sort by abundance
        std::sort(profile.begin(), profile.end(), 
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        return profile;
    }
};

// Example usage with streaming
int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <kraken_db_path> <fastq_file>\n";
        return 1;
    }
    
    std::string db_path = argv[1];
    std::string fastq_file = argv[2];
    
    try {
        // Initialize pipeline
        HybridKrakenPipeline pipeline(db_path);
        
        // Process file in batches
        FastqReader reader(fastq_file);
        ReadBatch batch;
        size_t total_reads = 0;
        size_t classified_reads = 0;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        while (reader.read_batch(batch, 10000)) {
            auto results = pipeline.classify_batch(batch.sequences);
            
            // Count classified reads
            for (const auto& r : results) {
                if (r.taxon_id > 0) classified_reads++;
            }
            
            total_reads += batch.size();
            
            if (total_reads % 100000 == 0) {
                std::cout << "Processed " << total_reads << " reads, "
                          << classified_reads << " classified ("
                          << (100.0 * classified_reads / total_reads) << "%)\n";
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Print results
        std::cout << "\n=== Classification Summary ===\n";
        std::cout << "Total reads: " << total_reads << "\n";
        std::cout << "Classified reads: " << classified_reads << " ("
                  << (100.0 * classified_reads / total_reads) << "%)\n";
        std::cout << "Processing time: " << duration.count() << " seconds\n";
        std::cout << "Reads per second: " << total_reads / duration.count() << "\n";
        
        // Print top taxa
        auto profile = pipeline.get_abundance_profile();
        std::cout << "\n=== Top 20 Taxa ===\n";
        for (size_t i = 0; i < std::min(size_t(20), profile.size()); i++) {
            std::cout << "Taxon " << profile[i].first 
                      << ": " << profile[i].second << " reads\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}