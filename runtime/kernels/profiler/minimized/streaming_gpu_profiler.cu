#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
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

// Include the missing OrganismProfile struct definition
struct OrganismProfile {
    uint32_t taxonomy_id;
    std::string name;
    std::string taxonomy_path;
    float abundance;
    float coverage_breadth;
    uint32_t unique_minimizers;
    uint32_t total_hits;
    float confidence_score;
};

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
    thrust::device_vector<uint32_t> d_organism_unique_counts;
    
    // Read processing
    thrust::device_vector<uint64_t> d_read_minimizers;
    thrust::device_vector<uint32_t> d_read_batch_offsets;
    
    // CUDA streams for async operations
    cudaStream_t compute_stream;
    cudaStream_t memory_stream;
    
    StreamingConfig config;
    int k = 35, m = 31;
    uint32_t max_organisms = 50000;  // Increased for larger databases
    
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
        d_organism_unique_counts.resize(max_organisms, 0);
        
        print_memory_config();
    }
    
    ~StreamingGPUProfiler() {
        cudaStreamDestroy(compute_stream);
        cudaStreamDestroy(memory_stream);
    }
    
    void print_memory_config() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        std::cout << "GPU Memory Configuration:" << std::endl;
        std::cout << "- Total GPU memory: " << total_mem / (1024*1024*1024) << " GB" << std::endl;
        std::cout << "- Available memory: " << free_mem / (1024*1024*1024) << " GB" << std::endl;
        std::cout << "- Configured for database chunks: " << config.max_gpu_memory_gb << " GB" << std::endl;
        std::cout << "- Chunk size: " << config.chunk_size_mb << " MB" << std::endl;
        std::cout << "- CPU chunk cache: " << config.max_chunks_in_memory << " chunks" << std::endl;
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
        std::cout << "Creating database chunks for streaming..." << std::endl;
        
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
        
        std::cout << "Created " << database_chunks.size() << " database chunks" << std::endl;
        
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
        
        // Now we're at the start of the minimizer index data
        std::cout << "Positioned at minimizer data section (offset: " << in.tellg() << ")" << std::endl;
    }
    
    void load_all_minimizer_entries(std::ifstream& in, 
                                   std::vector<std::tuple<uint64_t, uint32_t, float>>& entries) {
        std::cout << "Loading all minimizer entries for chunking..." << std::endl;
        
        // Read minimizer index in the format from our minimizer database
        uint64_t hash;
        uint32_t num_entries_for_hash;
        
        while (in.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
            in.read(reinterpret_cast<char*>(&num_entries_for_hash), sizeof(num_entries_for_hash));
            
            for (uint32_t i = 0; i < num_entries_for_hash; i++) {
                // Read minimizer entry in packed format: hash, taxonomy_id, uniqueness_score
                uint64_t minimizer_hash;
                uint32_t taxonomy_id;
                uint8_t uniqueness_score;
                
                in.read(reinterpret_cast<char*>(&minimizer_hash), sizeof(minimizer_hash));
                in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
                in.read(reinterpret_cast<char*>(&uniqueness_score), sizeof(uniqueness_score));
                
                // Verify hash consistency
                if (minimizer_hash != hash) {
                    std::cerr << "Warning: Hash mismatch in database at position " << entries.size() << std::endl;
                }
                
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
                } else {
                    std::cerr << "Warning: Unknown taxonomy_id " << taxonomy_id << " in database" << std::endl;
                }
            }
            
            if (entries.size() % 1000000 == 0) {
                std::cout << "\rLoaded " << entries.size() << " entries..." << std::flush;
            }
        }
        
        std::cout << "\nFinished loading " << entries.size() << " total entries" << std::endl;
    }
    
public:
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
        
        // Copy to GPU asynchronously
        thrust::copy_async(chunk->minimizer_hashes.begin(), chunk->minimizer_hashes.end(),
                          d_chunk_hashes.begin(), memory_stream);
        thrust::copy_async(chunk->organism_ids.begin(), chunk->organism_ids.end(),
                          d_chunk_organism_ids.begin(), memory_stream);
        thrust::copy_async(chunk->weights.begin(), chunk->weights.end(),
                          d_chunk_weights.begin(), memory_stream);
        
        // Synchronize memory stream
        cudaStreamSynchronize(memory_stream);
        
        chunk_loaded_flags[chunk_id] = true;
        
        std::cout << "Loaded chunk " << chunk_id << " to GPU (" 
                  << chunk->minimizer_hashes.size() << " entries)" << std::endl;
    }
    
    void process_reads_with_streaming(const std::vector<std::string>& reads) {
        std::cout << "Processing " << reads.size() << " reads with streaming..." << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Extract all minimizers from reads first
        std::vector<uint64_t> all_minimizers;
        extract_all_minimizers_cpu(reads, all_minimizers);
        
        // Find which chunks we need
        auto relevant_chunks = find_relevant_chunks(all_minimizers);
        std::cout << "Need to process " << relevant_chunks.size() 
                  << " chunks out of " << database_chunks.size() << std::endl;
        
        // Reset organism counters
        thrust::fill(d_organism_hits.begin(), d_organism_hits.end(), 0);
        thrust::fill(d_organism_scores.begin(), d_organism_scores.end(), 0.0f);
        
        // Process each relevant chunk
        for (size_t chunk_idx : relevant_chunks) {
            std::cout << "Processing chunk " << chunk_idx << "..." << std::endl;
            
            // Load chunk to GPU
            load_chunk_to_gpu(chunk_idx);
            
            // Match minimizers against this chunk
            match_minimizers_against_chunk(all_minimizers);
            
            // Optional: Free GPU memory for this chunk if not caching
            if (!config.enable_chunk_caching) {
                d_chunk_hashes.clear();
                d_chunk_organism_ids.clear();
                d_chunk_weights.clear();
                chunk_loaded_flags[chunk_idx] = false;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Streaming profiling completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Processed " << relevant_chunks.size() << " chunks, "
                  << all_minimizers.size() << " minimizers" << std::endl;
    }
    
private:
    void extract_all_minimizers_cpu(const std::vector<std::string>& reads,
                                   std::vector<uint64_t>& minimizers) {
        std::cout << "Extracting minimizers from " << reads.size() << " reads..." << std::endl;
        
        int window_size = k - m + 1;
        
        for (const auto& read : reads) {
            if (read.length() < k) continue;
            
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
                }
            }
        }
        
        std::cout << "Extracted " << minimizers.size() << " minimizers" << std::endl;
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
    
    void match_minimizers_against_chunk(const std::vector<uint64_t>& query_minimizers) {
        // Copy query minimizers to GPU
        d_read_minimizers.resize(query_minimizers.size());
        thrust::copy(query_minimizers.begin(), query_minimizers.end(), d_read_minimizers.begin());
        
        // Launch GPU kernel for matching against current chunk
        int block_size = 256;
        int num_blocks = (query_minimizers.size() + block_size - 1) / block_size;
        
        match_chunk_kernel<<<num_blocks, block_size, 0, compute_stream>>>(
            thrust::raw_pointer_cast(d_read_minimizers.data()),
            query_minimizers.size(),
            thrust::raw_pointer_cast(d_chunk_hashes.data()),
            thrust::raw_pointer_cast(d_chunk_organism_ids.data()),
            thrust::raw_pointer_cast(d_chunk_weights.data()),
            d_chunk_hashes.size(),
            thrust::raw_pointer_cast(d_organism_hits.data()),
            thrust::raw_pointer_cast(d_organism_scores.data())
        );
        
        cudaStreamSynchronize(compute_stream);
    }
    
public:
    std::vector<OrganismProfile> get_final_results() {
        // Copy results back from GPU and generate profiles
        thrust::host_vector<uint32_t> h_organism_hits = d_organism_hits;
        thrust::host_vector<float> h_organism_scores = d_organism_scores;
        
        std::vector<OrganismProfile> profiles;
        
        // Calculate total abundance for normalization
        float total_abundance = 0.0f;
        for (size_t i = 0; i < organism_id_list.size() && i < max_organisms; i++) {
            total_abundance += h_organism_scores[i];
        }
        
        for (size_t i = 0; i < organism_id_list.size() && i < max_organisms; i++) {
            uint32_t tax_id = organism_id_list[i];
            uint32_t hits = h_organism_hits[i];
            float score = h_organism_scores[i];
            
            if (hits < 10) continue;  // Minimum threshold
            
            OrganismProfile profile;
            profile.taxonomy_id = tax_id;
            profile.name = organism_names[tax_id];
            profile.taxonomy_path = organism_taxonomy[tax_id];
            profile.total_hits = hits;
            profile.abundance = (total_abundance > 0) ? score / total_abundance : 0.0f;
            
            // Simplified coverage calculation for streaming
            uint32_t expected_minimizers = organism_expected_minimizers[tax_id];
            profile.coverage_breadth = (expected_minimizers > 0) ? 
                                     std::min(1.0f, (float)hits / expected_minimizers) : 0.0f;
            
            profile.confidence_score = std::min(1.0f, std::log10(hits + 1) / 3.0f);
            
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
    
    void print_memory_usage() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        std::cout << "GPU memory usage: " << (total_mem - free_mem) / (1024*1024) 
                  << " MB used, " << free_mem / (1024*1024) << " MB free" << std::endl;
    }
};

// GPU kernel for chunk-based matching
__global__ void match_chunk_kernel(
    const uint64_t* query_minimizers,
    uint32_t num_queries,
    const uint64_t* chunk_hashes,
    const uint32_t* chunk_organism_ids,
    const float* chunk_weights,
    uint32_t chunk_size,
    uint32_t* organism_hits,
    float* organism_scores
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_queries) return;
    
    uint64_t query_hash = query_minimizers[tid];
    
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
            
            atomicAdd(&organism_hits[org_id], 1);
            atomicAdd(&organism_scores[org_id], weight);
            
            // Check neighboring entries
            int left_check = mid - 1;
            while (left_check >= 0 && chunk_hashes[left_check] == query_hash) {
                atomicAdd(&organism_hits[chunk_organism_ids[left_check]], 1);
                atomicAdd(&organism_scores[chunk_organism_ids[left_check]], chunk_weights[left_check]);
                left_check--;
            }
            
            int right_check = mid + 1;
            while (right_check < chunk_size && chunk_hashes[right_check] == query_hash) {
                atomicAdd(&organism_hits[chunk_organism_ids[right_check]], 1);
                atomicAdd(&organism_scores[chunk_organism_ids[right_check]], chunk_weights[right_check]);
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

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <large_database> <fastq_file> [output_prefix]" << std::endl;
        std::cerr << "\nStreaming GPU profiler for large databases (>48GB)" << std::endl;
        std::cerr << "Uses database chunking and CPU-GPU streaming" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string fastq_file = argv[2];
    std::string output_prefix = (argc > 3) ? argv[3] : "streaming_profile";
    
    try {
        // Configure for large database streaming
        StreamingConfig config;
        config.max_gpu_memory_gb = 32;      // Use 32GB for chunks, leave 16GB for reads
        config.chunk_size_mb = 2048;        // 2GB chunks
        config.max_chunks_in_memory = 8;    // Keep 8 chunks in CPU RAM
        config.use_async_streaming = true;
        
        StreamingGPUProfiler profiler(config);
        
        // Load database metadata and create chunks
        profiler.load_database_metadata(database_path);
        profiler.create_database_chunks(database_path);
        
        // Read FASTQ file
        std::vector<std::string> reads;
        std::ifstream file(fastq_file);
        std::string line;
        int line_count = 0;
        
        std::cout << "Reading FASTQ file: " << fastq_file << std::endl;
        
        while (std::getline(file, line)) {
            line_count++;
            if (line_count % 4 == 2) {  // Sequence line
                if (!line.empty() && line.length() >= 35) {
                    reads.push_back(line);
                }
            }
            
            if (line_count % 400000 == 0) {
                std::cout << "\rRead " << (line_count / 4) << " sequences..." << std::flush;
            }
            
            // Process all reads for production
            if (reads.size() >= 2000000) {  // 2M reads for testing
                std::cout << "\nLimited to " << reads.size() << " reads for testing" << std::endl;
                break;
            }
        }
        
        std::cout << "\nLoaded " << reads.size() << " reads for streaming profiling" << std::endl;
        
        // Process with streaming
        profiler.process_reads_with_streaming(reads);
        
        // Get results
        auto profiles = profiler.get_final_results();
        
        // Print results
        std::cout << "\n=== STREAMING GPU PROFILING RESULTS ===" << std::endl;
        std::cout << "Detected " << profiles.size() << " organisms" << std::endl;
        
        for (size_t i = 0; i < std::min(size_t(20), profiles.size()); i++) {
            const auto& profile = profiles[i];
            std::cout << std::fixed << std::setprecision(4)
                      << i+1 << ". " << profile.name << ": " 
                      << (profile.abundance * 100) << "% ("
                      << profile.total_hits << " hits)" << std::endl;
        }
        
        profiler.print_memory_usage();
        
        std::cout << "\n=== STREAMING PROFILING COMPLETE ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}