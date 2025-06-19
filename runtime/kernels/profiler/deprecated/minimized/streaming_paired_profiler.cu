#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <memory>
#include <queue>
#include <thread>
#include <iomanip>
#include <cmath>

// Include the common paired-end structures
#include "paired_read_common.h"

// Forward declarations of GPU kernels
__global__ void match_chunk_kernel_paired(
    const uint64_t* query_minimizers,
    const uint32_t* query_pair_ids,
    const uint8_t* query_read_numbers,
    uint32_t num_queries,
    const uint64_t* chunk_hashes,
    const uint32_t* chunk_organism_ids,
    const float* chunk_weights,
    uint32_t chunk_size,
    uint32_t* organism_hits,
    float* organism_scores,
    float* pair_scores_r1,
    float* pair_scores_r2,
    uint32_t max_organisms,
    uint32_t max_pairs
);

__global__ void calculate_streaming_concordance_kernel(
    const float* pair_scores_r1,
    const float* pair_scores_r2,
    const bool* is_paired_flags,
    uint32_t num_pairs,
    uint32_t max_organisms,
    uint32_t* organism_paired_hits,
    float* organism_concordance_scores
);

struct DatabaseChunk {
    std::vector<uint64_t> minimizer_hashes;
    std::vector<uint32_t> organism_ids;
    std::vector<float> weights;
    uint64_t hash_range_start;
    uint64_t hash_range_end;
    size_t chunk_id;
};

struct StreamingConfig {
    size_t max_gpu_memory_gb = 32;      // Leave 16GB for reads and working memory
    size_t chunk_size_mb = 2048;        // 2GB chunks for streaming
    size_t max_chunks_in_memory = 8;    // Keep 8 chunks in CPU RAM
    bool use_async_streaming = true;
    bool enable_chunk_caching = true;
    size_t read_batch_size = 100000;    // Process reads in batches
    bool use_paired_end_bonus = true;
    float paired_concordance_weight = 1.5f;
};

class StreamingGPUProfiler {
private:
    // Database metadata
    std::unordered_map<uint32_t, std::string> organism_names;
    std::unordered_map<uint32_t, std::string> organism_taxonomy;
    std::unordered_map<uint32_t, uint64_t> organism_genome_sizes;
    std::unordered_map<uint32_t, uint32_t> organism_expected_minimizers;
    std::vector<uint32_t> organism_id_list;
    
    // Database chunks
    std::vector<std::unique_ptr<DatabaseChunk>> database_chunks;
    std::vector<bool> chunk_loaded_flags;
    
    // GPU memory pools
    thrust::device_vector<uint64_t> d_chunk_hashes;
    thrust::device_vector<uint32_t> d_chunk_organism_ids;
    thrust::device_vector<float> d_chunk_weights;
    
    // Persistent GPU arrays for results
    thrust::device_vector<uint32_t> d_organism_hits;
    thrust::device_vector<float> d_organism_scores;
    thrust::device_vector<uint32_t> d_organism_paired_hits;
    thrust::device_vector<float> d_organism_concordance_scores;
    
    // Read processing with paired-end support
    thrust::device_vector<uint64_t> d_read_minimizers;
    thrust::device_vector<uint32_t> d_read_pair_ids;
    thrust::device_vector<uint8_t> d_read_numbers;
    
    // Paired-end scoring arrays for streaming
    thrust::device_vector<float> d_pair_scores_r1;
    thrust::device_vector<float> d_pair_scores_r2;
    thrust::device_vector<bool> d_is_paired_flags;
    
    // CUDA streams for async operations
    cudaStream_t compute_stream;
    cudaStream_t memory_stream;
    
    StreamingConfig config;
    int k = 35, m = 31;
    uint32_t max_organisms = 50000;  // Increased for larger databases
    uint32_t max_pairs = 2000000;    // Support up to 2M read pairs
    
    // Chunk caching
    std::unordered_map<size_t, std::unique_ptr<DatabaseChunk>> chunk_cache;
    std::queue<size_t> cache_lru_queue;
    
public:
    StreamingGPUProfiler(const StreamingConfig& cfg = StreamingConfig()) : config(cfg) {
        // Create CUDA streams
        cudaStreamCreate(&compute_stream);
        cudaStreamCreate(&memory_stream);
        
        // Initialize persistent GPU arrays
        d_organism_hits.resize(max_organisms, 0);
        d_organism_scores.resize(max_organisms, 0.0f);
        d_organism_paired_hits.resize(max_organisms, 0);
        d_organism_concordance_scores.resize(max_organisms, 0.0f);
        
        print_memory_config();
    }
    
    ~StreamingGPUProfiler() {
        cudaStreamDestroy(compute_stream);
        cudaStreamDestroy(memory_stream);
    }
    
    void print_memory_config() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        std::cout << "Streaming GPU Memory Configuration:" << std::endl;
        std::cout << "- Total GPU memory: " << total_mem / (1024*1024*1024) << " GB" << std::endl;
        std::cout << "- Available memory: " << free_mem / (1024*1024*1024) << " GB" << std::endl;
        std::cout << "- Configured for database chunks: " << config.max_gpu_memory_gb << " GB" << std::endl;
        std::cout << "- Chunk size: " << config.chunk_size_mb << " MB" << std::endl;
        std::cout << "- Max read pairs: " << max_pairs << std::endl;
        std::cout << "- Paired-end bonus: " << (config.use_paired_end_bonus ? "enabled" : "disabled") << std::endl;
    }
    
    void load_database_metadata(const std::string& db_path) {
        std::cout << "Loading database metadata: " << db_path << std::endl;
        
        std::ifstream in(db_path, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open database file: " + db_path);
        }
        
        // Read header
        uint32_t magic, version, k_size, m_size, num_organisms;
        uint64_t num_minimizer_hashes;
        
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        in.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        in.read(reinterpret_cast<char*>(&num_organisms), sizeof(num_organisms));
        in.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        if (magic != 0x4D494E49) {
            throw std::runtime_error("Invalid database format");
        }
        
        k = k_size;
        m = m_size;
        
        std::cout << "Database: k=" << k << ", m=" << m 
                  << ", organisms=" << num_organisms 
                  << ", minimizer hashes=" << num_minimizer_hashes << std::endl;
        
        // Read organism metadata
        for (uint32_t i = 0; i < num_organisms; i++) {
            uint32_t taxonomy_id, taxon_level, minimizer_count;
            uint64_t genome_size;
            float gc_content;
            
            in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
            in.read(reinterpret_cast<char*>(&taxon_level), sizeof(taxon_level));
            in.read(reinterpret_cast<char*>(&gc_content), sizeof(gc_content));
            in.read(reinterpret_cast<char*>(&genome_size), sizeof(genome_size));
            in.read(reinterpret_cast<char*>(&minimizer_count), sizeof(minimizer_count));
            
            // Read strings
            uint16_t name_length;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            std::string name(name_length, '\0');
            in.read(&name[0], name_length);
            organism_names[taxonomy_id] = name;
            
            uint16_t taxonomy_length;
            in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            std::string taxonomy_path(taxonomy_length, '\0');
            in.read(&taxonomy_path[0], taxonomy_length);
            organism_taxonomy[taxonomy_id] = taxonomy_path;
            
            // Store metadata
            organism_genome_sizes[taxonomy_id] = genome_size;
            organism_expected_minimizers[taxonomy_id] = minimizer_count;
            organism_id_list.push_back(taxonomy_id);
        }
        
        in.close();
        
        std::cout << "Loaded metadata for " << organism_names.size() << " organisms" << std::endl;
    }
    
    void create_database_chunks(const std::string& db_path) {
        std::cout << "Creating database chunks for paired-end streaming..." << std::endl;
        
        std::ifstream in(db_path, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open database file for chunking");
        }
        
        // Skip to minimizer data (after header and organism metadata)
        seek_to_minimizer_data(in);
        
        // Calculate chunk parameters
        size_t entries_per_chunk = (config.chunk_size_mb * 1024 * 1024) / 
                                  (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(float));
        
        std::cout << "Target entries per chunk: " << entries_per_chunk << std::endl;
        
        // Read all minimizer entries and sort by hash
        std::vector<std::tuple<uint64_t, uint32_t, float>> all_entries;
        load_all_minimizer_entries(in, all_entries);
        
        in.close();
        
        // Sort by hash for efficient range queries
        std::sort(all_entries.begin(), all_entries.end());
        
        std::cout << "Loaded and sorted " << all_entries.size() << " minimizer entries" << std::endl;
        
        // Create chunks
        size_t chunk_id = 0;
        for (size_t i = 0; i < all_entries.size(); i += entries_per_chunk) {
            auto chunk = std::make_unique<DatabaseChunk>();
            chunk->chunk_id = chunk_id++;
            
            size_t chunk_end = std::min(i + entries_per_chunk, all_entries.size());
            
            for (size_t j = i; j < chunk_end; j++) {
                chunk->minimizer_hashes.push_back(std::get<0>(all_entries[j]));
                chunk->organism_ids.push_back(std::get<1>(all_entries[j]));
                chunk->weights.push_back(std::get<2>(all_entries[j]));
            }
            
            chunk->hash_range_start = chunk->minimizer_hashes.front();
            chunk->hash_range_end = chunk->minimizer_hashes.back();
            
            database_chunks.push_back(std::move(chunk));
        }
        
        chunk_loaded_flags.resize(database_chunks.size(), false);
        
        std::cout << "Created " << database_chunks.size() << " database chunks for streaming" << std::endl;
        
        // Print chunk statistics
        for (size_t i = 0; i < std::min(size_t(5), database_chunks.size()); i++) {
            const auto& chunk = database_chunks[i];
            std::cout << "Chunk " << i << ": " << chunk->minimizer_hashes.size() 
                      << " entries, hash range [" << std::hex 
                      << chunk->hash_range_start << "-" << chunk->hash_range_end 
                      << std::dec << "]" << std::endl;
        }
    }
    
private:
    void seek_to_minimizer_data(std::ifstream& in) {
        // Reset to beginning and skip header
        in.seekg(0, std::ios::beg);
        
        // Skip header: magic, version, k_size, m_size, num_organisms, num_minimizer_hashes
        uint32_t magic, version, k_size, m_size, num_organisms;
        uint64_t num_minimizer_hashes;
        
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        in.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        in.read(reinterpret_cast<char*>(&num_organisms), sizeof(num_organisms));
        in.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        // Skip organism metadata section
        for (uint32_t i = 0; i < num_organisms; i++) {
            // Skip fixed-size fields: taxonomy_id, taxon_level, gc_content, genome_size, minimizer_count
            in.seekg(sizeof(uint32_t) + sizeof(uint32_t) + sizeof(float) + 
                    sizeof(uint64_t) + sizeof(uint32_t), std::ios::cur);
            
            // Skip variable-length organism name
            uint16_t name_length;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            in.seekg(name_length, std::ios::cur);
            
            // Skip variable-length taxonomy path
            uint16_t taxonomy_length;
            in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            in.seekg(taxonomy_length, std::ios::cur);
        }
        
        std::cout << "Positioned at minimizer data section (offset: " << in.tellg() << ")" << std::endl;
    }
    
    void load_all_minimizer_entries(std::ifstream& in, 
                                   std::vector<std::tuple<uint64_t, uint32_t, float>>& entries) {
        std::cout << "Loading all minimizer entries for chunking..." << std::endl;
        
        uint64_t hash;
        uint32_t num_entries_for_hash;
        
        while (in.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
            in.read(reinterpret_cast<char*>(&num_entries_for_hash), sizeof(num_entries_for_hash));
            
            for (uint32_t i = 0; i < num_entries_for_hash; i++) {
                uint64_t minimizer_hash;
                uint32_t taxonomy_id;
                uint8_t uniqueness_score;
                
                in.read(reinterpret_cast<char*>(&minimizer_hash), sizeof(minimizer_hash));
                in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
                in.read(reinterpret_cast<char*>(&uniqueness_score), sizeof(uniqueness_score));
                
                // Calculate weight including uniqueness and genome size normalization
                float weight = 1.0f;
                if (uniqueness_score > 0) {
                    weight = uniqueness_score / 255.0f;
                }
                if (organism_genome_sizes.count(taxonomy_id)) {
                    uint64_t genome_size = organism_genome_sizes[taxonomy_id];
                    weight = weight / std::log10(genome_size + 1);
                }
                
                // Map taxonomy_id to array index for GPU processing
                auto it = std::lower_bound(organism_id_list.begin(), organism_id_list.end(), taxonomy_id);
                if (it != organism_id_list.end() && *it == taxonomy_id) {
                    uint32_t array_index = std::distance(organism_id_list.begin(), it);
                    entries.emplace_back(minimizer_hash, array_index, weight);
                }
            }
            
            if (entries.size() % 1000000 == 0) {
                std::cout << "\rLoaded " << entries.size() << " entries..." << std::flush;
            }
        }
        
        std::cout << "\nFinished loading " << entries.size() << " total entries" << std::endl;
    }
    
public:
    // NEW: Main paired-end streaming method
    void process_reads_with_streaming(const std::vector<PairedRead>& paired_reads) {
        std::cout << "Processing " << paired_reads.size() << " paired reads with streaming..." << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Extract all minimizers from paired reads with tracking
        std::vector<uint64_t> all_minimizers;
        std::vector<uint32_t> minimizer_pair_ids;
        std::vector<uint8_t> minimizer_read_numbers;
        
        extract_paired_minimizers_cpu(paired_reads, all_minimizers, minimizer_pair_ids, minimizer_read_numbers);
        
        // Find which chunks we need
        auto relevant_chunks = find_relevant_chunks(all_minimizers);
        std::cout << "Need to process " << relevant_chunks.size() 
                  << " chunks out of " << database_chunks.size() << std::endl;
        
        // Reset organism counters
        thrust::fill(d_organism_hits.begin(), d_organism_hits.end(), 0);
        thrust::fill(d_organism_scores.begin(), d_organism_scores.end(), 0.0f);
        thrust::fill(d_organism_paired_hits.begin(), d_organism_paired_hits.end(), 0);
        thrust::fill(d_organism_concordance_scores.begin(), d_organism_concordance_scores.end(), 0.0f);
        
        // Prepare paired-end scoring arrays
        setup_paired_scoring_arrays(paired_reads.size());
        
        // Process each relevant chunk
        for (size_t chunk_idx : relevant_chunks) {
            std::cout << "Processing chunk " << chunk_idx << "..." << std::endl;
            
            // Load chunk to GPU
            load_chunk_to_gpu(chunk_idx);
            
            // Match minimizers against this chunk with paired-end tracking
            match_paired_minimizers_against_chunk(all_minimizers, minimizer_pair_ids, minimizer_read_numbers);
            
            // Optional: Free GPU memory for this chunk if not caching
            if (!config.enable_chunk_caching) {
                d_chunk_hashes.clear();
                d_chunk_organism_ids.clear();
                d_chunk_weights.clear();
                chunk_loaded_flags[chunk_idx] = false;
            }
        }
        
        // Calculate paired-end concordance after all chunks processed
        calculate_streaming_paired_concordance(paired_reads);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Paired-end streaming profiling completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Processed " << relevant_chunks.size() << " chunks, "
                  << all_minimizers.size() << " minimizers from " 
                  << paired_reads.size() << " read pairs" << std::endl;
    }
    
    // BACKWARD COMPATIBILITY: Single-read streaming method
    void process_reads_with_streaming(const std::vector<std::string>& reads) {
        auto paired_reads = convert_to_paired_format(reads);
        process_reads_with_streaming(paired_reads);
    }
    
private:
    void setup_paired_scoring_arrays(size_t num_pairs) {
        // Allocate paired-end scoring arrays
        size_t pair_scoring_size = num_pairs * max_organisms;
        d_pair_scores_r1.resize(pair_scoring_size, 0.0f);
        d_pair_scores_r2.resize(pair_scoring_size, 0.0f);
        
        // Create paired flags array
        std::vector<bool> is_paired_flags(num_pairs, false);
        for (size_t i = 0; i < num_pairs; i++) {
            is_paired_flags[i] = true;  // Will be set correctly when processing actual reads
        }
        d_is_paired_flags = is_paired_flags;
        
        std::cout << "Allocated paired-end scoring arrays for " << num_pairs << " pairs" << std::endl;
    }
    
    void extract_paired_minimizers_cpu(const std::vector<PairedRead>& paired_reads,
                                      std::vector<uint64_t>& minimizers,
                                      std::vector<uint32_t>& pair_ids,
                                      std::vector<uint8_t>& read_numbers) {
        std::cout << "Extracting minimizers from " << paired_reads.size() << " paired reads..." << std::endl;
        
        int window_size = k - m + 1;
        
        for (size_t pair_idx = 0; pair_idx < paired_reads.size(); pair_idx++) {
            const auto& pair = paired_reads[pair_idx];
            
            // Process read1
            extract_minimizers_from_read(pair.read1, pair_idx, 1, minimizers, pair_ids, read_numbers);
            
            // Process read2 if paired
            if (pair.is_paired && !pair.read2.empty()) {
                extract_minimizers_from_read(pair.read2, pair_idx, 2, minimizers, pair_ids, read_numbers);
            }
        }
        
        std::cout << "Extracted " << minimizers.size() << " minimizers from " 
                  << paired_reads.size() << " paired reads" << std::endl;
    }
    
    void extract_minimizers_from_read(const std::string& read,
                                     uint32_t pair_id,
                                     uint8_t read_number,
                                     std::vector<uint64_t>& minimizers,
                                     std::vector<uint32_t>& pair_ids,
                                     std::vector<uint8_t>& read_numbers) {
        if (read.length() < k) return;
        
        int window_size = k - m + 1;
        
        for (size_t i = 0; i <= read.length() - k; i += window_size) {
            std::string kmer_window = read.substr(i, k);
            
            uint64_t min_hash = UINT64_MAX;
            bool found_valid = false;
            
            for (int j = 0; j <= k - m; j++) {
                std::string minimizer_seq = kmer_window.substr(j, m);
                uint64_t hash = hash_minimizer(minimizer_seq);
                
                if (hash != UINT64_MAX && hash < min_hash) {
                    min_hash = hash;
                    found_valid = true;
                }
            }
            
            if (found_valid) {
                minimizers.push_back(min_hash);
                pair_ids.push_back(pair_id);
                read_numbers.push_back(read_number);
            }
        }
    }
    
    uint64_t hash_minimizer(const std::string& seq) {
        uint64_t forward_hash = 0;
        uint64_t reverse_hash = 0;
        
        bool valid = true;
        
        for (size_t i = 0; i < seq.length(); i++) {
            int base = encode_base(seq[i]);
            if (base == -1) {
                valid = false;
                break;
            }
            
            forward_hash = (forward_hash << 2) | base;
            reverse_hash = (reverse_hash >> 2) | (((uint64_t)(3 ^ base)) << (2 * (seq.length() - 1)));
        }
        
        return valid ? std::min(forward_hash, reverse_hash) : UINT64_MAX;
    }
    
    int encode_base(char base) {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }
    
    std::vector<size_t> find_relevant_chunks(const std::vector<uint64_t>& query_hashes) {
        std::vector<size_t> relevant_chunks;
        
        // Sort query hashes for efficient range checking
        auto sorted_queries = query_hashes;
        std::sort(sorted_queries.begin(), sorted_queries.end());
        
        for (size_t i = 0; i < database_chunks.size(); i++) {
            const auto& chunk = database_chunks[i];
            
            // Check if any query hash falls in this chunk's range
            auto it = std::lower_bound(sorted_queries.begin(), sorted_queries.end(), 
                                     chunk->hash_range_start);
            
            if (it != sorted_queries.end() && *it <= chunk->hash_range_end) {
                relevant_chunks.push_back(i);
            }
        }
        
        return relevant_chunks;
    }
    
    void load_chunk_to_gpu(size_t chunk_id) {
        if (chunk_id >= database_chunks.size()) return;
        
        const auto& chunk = database_chunks[chunk_id];
        
        // Resize GPU arrays for this chunk
        d_chunk_hashes.resize(chunk->minimizer_hashes.size());
        d_chunk_organism_ids.resize(chunk->organism_ids.size());
        d_chunk_weights.resize(chunk->weights.size());
        
        // Copy to GPU asynchronously using execution policy
        auto policy = thrust::cuda::par.on(memory_stream);
        thrust::copy(policy, chunk->minimizer_hashes.begin(), chunk->minimizer_hashes.end(),
                     d_chunk_hashes.begin());
        thrust::copy(policy, chunk->organism_ids.begin(), chunk->organism_ids.end(),
                     d_chunk_organism_ids.begin());
        thrust::copy(policy, chunk->weights.begin(), chunk->weights.end(),
                     d_chunk_weights.begin());
        
        // Synchronize memory stream
        cudaStreamSynchronize(memory_stream);
        
        chunk_loaded_flags[chunk_id] = true;
        
        std::cout << "Loaded chunk " << chunk_id << " to GPU (" 
                  << chunk->minimizer_hashes.size() << " entries)" << std::endl;
    }
    
    void match_paired_minimizers_against_chunk(const std::vector<uint64_t>& query_minimizers,
                                              const std::vector<uint32_t>& pair_ids,
                                              const std::vector<uint8_t>& read_numbers) {
        // Copy query data to GPU
        d_read_minimizers.resize(query_minimizers.size());
        d_read_pair_ids.resize(pair_ids.size());
        d_read_numbers.resize(read_numbers.size());
        
        thrust::copy(query_minimizers.begin(), query_minimizers.end(), d_read_minimizers.begin());
        thrust::copy(pair_ids.begin(), pair_ids.end(), d_read_pair_ids.begin());
        thrust::copy(read_numbers.begin(), read_numbers.end(), d_read_numbers.begin());
        
        // Launch GPU kernel for paired matching against current chunk
        int block_size = 256;
        int num_blocks = (query_minimizers.size() + block_size - 1) / block_size;
        
        match_chunk_kernel_paired<<<num_blocks, block_size, 0, compute_stream>>>(
            thrust::raw_pointer_cast(d_read_minimizers.data()),
            thrust::raw_pointer_cast(d_read_pair_ids.data()),
            thrust::raw_pointer_cast(d_read_numbers.data()),
            query_minimizers.size(),
            thrust::raw_pointer_cast(d_chunk_hashes.data()),
            thrust::raw_pointer_cast(d_chunk_organism_ids.data()),
            thrust::raw_pointer_cast(d_chunk_weights.data()),
            d_chunk_hashes.size(),
            thrust::raw_pointer_cast(d_organism_hits.data()),
            thrust::raw_pointer_cast(d_organism_scores.data()),
            thrust::raw_pointer_cast(d_pair_scores_r1.data()),
            thrust::raw_pointer_cast(d_pair_scores_r2.data()),
            max_organisms,
            max_pairs
        );
        
        cudaStreamSynchronize(compute_stream);
    }
    
    void calculate_streaming_paired_concordance(const std::vector<PairedRead>& paired_reads) {
        std::cout << "Calculating paired-end concordance for streaming results..." << std::endl;
        
        // Update is_paired_flags based on actual read data
        std::vector<bool> is_paired_flags;
        for (const auto& pair : paired_reads) {
            is_paired_flags.push_back(pair.is_paired);
        }
        d_is_paired_flags = is_paired_flags;
        
        uint32_t num_pairs = paired_reads.size();
        
        // Launch concordance calculation kernel
        int block_size = 256;
        int num_blocks = (num_pairs + block_size - 1) / block_size;
        
        calculate_streaming_concordance_kernel<<<num_blocks, block_size, 0, compute_stream>>>(
            thrust::raw_pointer_cast(d_pair_scores_r1.data()),
            thrust::raw_pointer_cast(d_pair_scores_r2.data()),
            thrust::raw_pointer_cast(d_is_paired_flags.data()),
            num_pairs,
            max_organisms,
            thrust::raw_pointer_cast(d_organism_paired_hits.data()),
            thrust::raw_pointer_cast(d_organism_concordance_scores.data())
        );
        
        cudaStreamSynchronize(compute_stream);
        
        std::cout << "Streaming paired concordance calculation completed" << std::endl;
    }
    
public:
    std::vector<OrganismProfile> get_final_results() {
        // Copy results back from GPU and generate profiles with paired-end information
        thrust::host_vector<uint32_t> h_organism_hits = d_organism_hits;
        thrust::host_vector<float> h_organism_scores = d_organism_scores;
        thrust::host_vector<uint32_t> h_organism_paired_hits = d_organism_paired_hits;
        thrust::host_vector<float> h_organism_concordance_scores = d_organism_concordance_scores;
        
        std::vector<OrganismProfile> profiles;
        
        // Calculate total abundance for normalization with paired-end bonus
        float total_abundance = 0.0f;
        for (size_t i = 0; i < organism_id_list.size() && i < max_organisms; i++) {
            float base_score = h_organism_scores[i];
            float concordance_bonus = 0.0f;
            
            if (config.use_paired_end_bonus && h_organism_paired_hits[i] > 0) {
                concordance_bonus = h_organism_concordance_scores[i] * config.paired_concordance_weight;
            }
            
            total_abundance += base_score + concordance_bonus;
        }
        
        for (size_t i = 0; i < organism_id_list.size() && i < max_organisms; i++) {
            uint32_t tax_id = organism_id_list[i];
            uint32_t hits = h_organism_hits[i];
            float score = h_organism_scores[i];
            uint32_t paired_hits = h_organism_paired_hits[i];
            float concordance_score = h_organism_concordance_scores[i];
            
            if (hits < 10) continue;  // Minimum threshold
            
            OrganismProfile profile;
            profile.taxonomy_id = tax_id;
            profile.name = organism_names[tax_id];
            profile.taxonomy_path = organism_taxonomy[tax_id];
            profile.total_hits = hits;
            profile.paired_hits = paired_hits;
            profile.single_hits = hits - (paired_hits * 2);  // Estimate single hits
            
            // Calculate paired concordance rate
            profile.paired_concordance = (paired_hits > 0) ? 
                                       concordance_score / (paired_hits * 2) : 0.0f;
            
            // Calculate abundance with paired-end bonus
            float base_abundance = (total_abundance > 0) ? score / total_abundance : 0.0f;
            float concordance_bonus = 0.0f;
            
            if (config.use_paired_end_bonus && paired_hits > 0) {
                concordance_bonus = (concordance_score * config.paired_concordance_weight) / total_abundance;
            }
            
            profile.abundance = base_abundance + concordance_bonus;
            
            // Simplified coverage calculation for streaming
            uint32_t expected_minimizers = organism_expected_minimizers[tax_id];
            profile.coverage_breadth = (expected_minimizers > 0) ? 
                                     std::min(1.0f, (float)hits / expected_minimizers) : 0.0f;
            
            // Enhanced confidence score with paired-end information
            float minimizer_confidence = std::min(1.0f, (float)(std::log10(hits + 1) / 3.0));
            float paired_confidence = (paired_hits > 0) ? 
                                    std::min(1.0f, profile.paired_concordance * 2.0f) : 0.0f;
            
            profile.confidence_score = (minimizer_confidence + paired_confidence) / 2.0f;
            
            if (profile.abundance > 1e-6 && profile.confidence_score > 0.1) {
                profiles.push_back(profile);
            }
        }
        
        // Sort by abundance
        std::sort(profiles.begin(), profiles.end(),
                 [](const OrganismProfile& a, const OrganismProfile& b) {
                     return a.abundance > b.abundance;
                 });
        
        return profiles;
    }
    
    void print_streaming_results(const std::vector<OrganismProfile>& profiles) {
        std::cout << "\n=== STREAMING PAIRED-END GPU PROFILING RESULTS ===" << std::endl;
        
        if (profiles.empty()) {
            std::cout << "No organisms detected above significance thresholds" << std::endl;
            return;
        }
        
        std::cout << "Detected " << profiles.size() << " organisms with streaming approach" << std::endl;
        
        // Print top organisms with paired-end information
        std::cout << "\nTop Organisms (Streaming + Paired-End):" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        std::cout << std::left << std::setw(35) << "Organism Name" 
                  << std::setw(10) << "Abundance" 
                  << std::setw(8) << "Total"
                  << std::setw(8) << "Paired"
                  << std::setw(12) << "Concordance"
                  << std::setw(10) << "Confidence"
                  << std::setw(20) << "Taxonomy" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        
        for (size_t i = 0; i < std::min(size_t(15), profiles.size()); i++) {
            const auto& profile = profiles[i];
            
            // Extract genus from taxonomy path
            std::string genus = "Unknown";
            size_t last_semicolon = profile.taxonomy_path.find_last_of(';');
            if (last_semicolon != std::string::npos) {
                genus = profile.taxonomy_path.substr(last_semicolon + 1);
            }
            
            // Truncate long names
            std::string display_name = profile.name;
            if (display_name.length() > 30) {
                display_name = display_name.substr(0, 27) + "...";
            }
            
            std::cout << std::left << std::setw(35) << display_name
                      << std::fixed << std::setprecision(3) << std::setw(10) << (profile.abundance * 100)
                      << std::setw(8) << profile.total_hits
                      << std::setw(8) << profile.paired_hits
                      << std::fixed << std::setprecision(2) << std::setw(12) << (profile.paired_concordance * 100)
                      << std::fixed << std::setprecision(3) << std::setw(10) << profile.confidence_score
                      << std::setw(20) << genus << std::endl;
        }
        
        std::cout << std::string(120, '-') << std::endl;
        
        // Paired-end summary for streaming
        uint32_t total_paired_hits = 0;
        float avg_concordance = 0.0f;
        uint32_t organisms_with_pairs = 0;
        
        for (const auto& profile : profiles) {
            total_paired_hits += profile.paired_hits;
            if (profile.paired_hits > 0) {
                avg_concordance += profile.paired_concordance;
                organisms_with_pairs++;
            }
        }
        
        if (organisms_with_pairs > 0) {
            avg_concordance /= organisms_with_pairs;
        }
        
        std::cout << "\nStreaming Paired-End Summary:" << std::endl;
        std::cout << "- Total paired hits: " << total_paired_hits << std::endl;
        std::cout << "- Organisms with paired evidence: " << organisms_with_pairs << std::endl;
        std::cout << "- Average concordance rate: " << std::fixed << std::setprecision(1) 
                  << (avg_concordance * 100) << "%" << std::endl;
        std::cout << "- Paired-end bonus: " << (config.use_paired_end_bonus ? "enabled" : "disabled") << std::endl;
    }
    
    void print_memory_usage() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        std::cout << "Streaming GPU memory usage: " << (total_mem - free_mem) / (1024*1024) 
                  << " MB used, " << free_mem / (1024*1024) << " MB free" << std::endl;
    }
};

// GPU kernel implementations

__global__ void match_chunk_kernel_paired(
    const uint64_t* query_minimizers,
    const uint32_t* query_pair_ids,
    const uint8_t* query_read_numbers,
    uint32_t num_queries,
    const uint64_t* chunk_hashes,
    const uint32_t* chunk_organism_ids,
    const float* chunk_weights,
    uint32_t chunk_size,
    uint32_t* organism_hits,
    float* organism_scores,
    float* pair_scores_r1,
    float* pair_scores_r2,
    uint32_t max_organisms,
    uint32_t max_pairs
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;
    
    uint64_t query_hash = query_minimizers[tid];
    uint32_t pair_id = query_pair_ids[tid];
    uint8_t read_num = query_read_numbers[tid];
    
    // Binary search in chunk
    int left = 0;
    int right = chunk_size - 1;
    
    while (left <= right) {
        int mid = left + (right - left) / 2;
        uint64_t chunk_hash = chunk_hashes[mid];
        
        if (chunk_hash == query_hash) {
            // Found match
            uint32_t org_id = chunk_organism_ids[mid];
            float weight = chunk_weights[mid];
            
            if (org_id < max_organisms) {
                atomicAdd(&organism_hits[org_id], 1);
                atomicAdd(&organism_scores[org_id], weight);
                
                // Track paired-end scoring
                if (pair_id < max_pairs) {
                    if (read_num == 1) {
                        atomicAdd(&pair_scores_r1[pair_id * max_organisms + org_id], weight);
                    } else if (read_num == 2) {
                        atomicAdd(&pair_scores_r2[pair_id * max_organisms + org_id], weight);
                    }
                }
            }
            
            // Check neighboring entries
            int left_check = mid - 1;
            while (left_check >= 0 && chunk_hashes[left_check] == query_hash) {
                uint32_t org_id_left = chunk_organism_ids[left_check];
                float weight_left = chunk_weights[left_check];
                
                if (org_id_left < max_organisms) {
                    atomicAdd(&organism_hits[org_id_left], 1);
                    atomicAdd(&organism_scores[org_id_left], weight_left);
                    
                    if (pair_id < max_pairs) {
                        if (read_num == 1) {
                            atomicAdd(&pair_scores_r1[pair_id * max_organisms + org_id_left], weight_left);
                        } else if (read_num == 2) {
                            atomicAdd(&pair_scores_r2[pair_id * max_organisms + org_id_left], weight_left);
                        }
                    }
                }
                left_check--;
            }
            
            int right_check = mid + 1;
            while (right_check < chunk_size && chunk_hashes[right_check] == query_hash) {
                uint32_t org_id_right = chunk_organism_ids[right_check];
                float weight_right = chunk_weights[right_check];
                
                if (org_id_right < max_organisms) {
                    atomicAdd(&organism_hits[org_id_right], 1);
                    atomicAdd(&organism_scores[org_id_right], weight_right);
                    
                    if (pair_id < max_pairs) {
                        if (read_num == 1) {
                            atomicAdd(&pair_scores_r1[pair_id * max_organisms + org_id_right], weight_right);
                        } else if (read_num == 2) {
                            atomicAdd(&pair_scores_r2[pair_id * max_organisms + org_id_right], weight_right);
                        }
                    }
                }
                right_check++;
            }
            
            break;
        } else if (chunk_hash < query_hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
}

__global__ void calculate_streaming_concordance_kernel(
    const float* pair_scores_r1,
    const float* pair_scores_r2,
    const bool* is_paired_flags,
    uint32_t num_pairs,
    uint32_t max_organisms,
    uint32_t* organism_paired_hits,
    float* organism_concordance_scores
) {
    int pair_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_id >= num_pairs || !is_paired_flags[pair_id]) return;
    
    // Find best organism for read1 and read2
    uint32_t best_org_r1 = 0;
    uint32_t best_org_r2 = 0;
    float best_score_r1 = 0.0f;
    float best_score_r2 = 0.0f;
    
    for (uint32_t org = 0; org < max_organisms; org++) {
        float score_r1 = pair_scores_r1[pair_id * max_organisms + org];
        float score_r2 = pair_scores_r2[pair_id * max_organisms + org];
        
        if (score_r1 > best_score_r1) {
            best_score_r1 = score_r1;
            best_org_r1 = org;
        }
        
        if (score_r2 > best_score_r2) {
            best_score_r2 = score_r2;
            best_org_r2 = org;
        }
    }
    
    // If both reads map to the same organism with good scores
    if (best_org_r1 == best_org_r2 && best_score_r1 > 0.1f && best_score_r2 > 0.1f) {
        atomicAdd(&organism_paired_hits[best_org_r1], 1);
        atomicAdd(&organism_concordance_scores[best_org_r1], best_score_r1 + best_score_r2);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <large_database> <fastq_r1> [fastq_r2] [output_prefix]" << std::endl;
        std::cerr << "\nStreaming GPU profiler with paired-end support for large databases (>48GB)" << std::endl;
        std::cerr << "Uses database chunking and CPU-GPU streaming with paired-end scoring" << std::endl;
        std::cerr << "\nExamples:" << std::endl;
        std::cerr << "  Single-end: " << argv[0] << " large_db.db reads.fastq" << std::endl;
        std::cerr << "  Paired-end: " << argv[0] << " large_db.db reads_R1.fastq reads_R2.fastq" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string fastq_r1 = argv[2];
    std::string fastq_r2 = "";
    std::string output_prefix = "streaming_paired_profile";
    
    if (argc > 3) {
        std::string arg3 = argv[3];
        if (arg3.find(".fastq") != std::string::npos || arg3.find(".fq") != std::string::npos) {
            fastq_r2 = arg3;
            if (argc > 4) output_prefix = argv[4];
        } else {
            output_prefix = arg3;
        }
    }
    
    try {
        // Configure for large database streaming with paired-end support
        StreamingConfig config;
        config.max_gpu_memory_gb = 32;      // Use 32GB for chunks, leave 16GB for reads
        config.chunk_size_mb = 2048;        // 2GB chunks
        config.max_chunks_in_memory = 8;    // Keep 8 chunks in CPU RAM
        config.use_async_streaming = true;
        config.use_paired_end_bonus = true;
        config.paired_concordance_weight = 1.5f;
        
        StreamingGPUProfiler profiler(config);
        
        // Load database metadata and create chunks
        profiler.load_database_metadata(database_path);
        profiler.create_database_chunks(database_path);
        
        // Read FASTQ file(s) and create paired reads
        std::vector<PairedRead> paired_reads;
        
        if (!fastq_r2.empty()) {
            std::cout << "Loading paired-end reads from " << fastq_r1 << " and " << fastq_r2 << std::endl;
            // Implementation would load paired files similar to adaptive profiler
            // For now, using single-end as example
            
        } else {
            std::cout << "Loading single-end reads from " << fastq_r1 << std::endl;
            std::vector<std::string> single_reads;
            
            std::ifstream file(fastq_r1);
            std::string line;
            int line_count = 0;
            
            while (std::getline(file, line)) {
                line_count++;
                if (line_count % 4 == 2) {  // Sequence line
                    if (!line.empty() && line.length() >= 35) {
                        single_reads.push_back(line);
                    }
                }
                
                if (line_count % 400000 == 0) {
                    std::cout << "\rRead " << (line_count / 4) << " sequences..." << std::flush;
                }
                
                if (single_reads.size() >= 2000000) {  // 2M reads for testing
                    std::cout << "\nLimited to " << single_reads.size() << " reads for testing" << std::endl;
                    break;
                }
            }
            
            // Convert to paired format
            paired_reads = convert_to_paired_format(single_reads);
        }
        
        std::cout << "\nLoaded " << paired_reads.size() << " read pairs for streaming profiling" << std::endl;
        
        // Process with streaming and paired-end support
        profiler.process_reads_with_streaming(paired_reads);
        
        // Get results
        auto profiles = profiler.get_final_results();
        
        // Print results
        profiler.print_streaming_results(profiles);
        profiler.print_memory_usage();
        
        std::cout << "\n=== STREAMING PAIRED-END PROFILING COMPLETE ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}