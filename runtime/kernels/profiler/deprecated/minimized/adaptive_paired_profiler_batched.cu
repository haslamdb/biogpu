#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <memory>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libgen.h>
#include <unistd.h>
#include <limits.h>
#include <zlib.h>
#include <filesystem>
#include <tuple>
#include <algorithm>

// Include common structures
#include "paired_read_common.h"
#include "sample_csv_parser.h"

struct DatabaseStats {
    uint64_t estimated_gpu_memory_mb;
    uint32_t num_organisms;
    uint64_t total_minimizer_entries;
    bool needs_streaming;
    uint32_t recommended_chunks;
};

// FastqRecord structure for batch processing
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Forward declare the GPU profiler class (we'll implement this inline)
class DirectGPUProfiler {
private:
    // Database metadata
    std::unordered_map<uint32_t, std::string> organism_names;
    std::unordered_map<uint32_t, std::string> organism_taxonomy;
    std::unordered_map<uint32_t, uint64_t> organism_genome_sizes;
    std::unordered_map<uint32_t, uint32_t> organism_expected_minimizers;
    std::vector<uint32_t> organism_id_list;
    
    // GPU database storage
    thrust::device_vector<uint64_t> d_database_hashes;
    thrust::device_vector<uint32_t> d_database_organism_ids;
    thrust::device_vector<float> d_database_weights;
    
    // GPU results storage - now accumulated across batches
    thrust::device_vector<uint32_t> d_organism_hits;
    thrust::device_vector<float> d_organism_scores;
    thrust::device_vector<uint32_t> d_organism_paired_hits;
    thrust::device_vector<float> d_organism_concordance_scores;
    
    int k = 35, m = 31;
    uint32_t max_organisms = 10000;
    
public:
    DirectGPUProfiler() {
        // Initialize GPU arrays
        d_organism_hits.resize(max_organisms, 0);
        d_organism_scores.resize(max_organisms, 0.0f);
        d_organism_paired_hits.resize(max_organisms, 0);
        d_organism_concordance_scores.resize(max_organisms, 0.0f);
    }
    
    void load_minimizer_database(const std::string& db_path) {
        std::cout << "Loading minimizer database directly: " << db_path << std::endl;
        auto load_start = std::chrono::high_resolution_clock::now();
        
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
        
        std::cout << "Database parameters: k=" << k << ", m=" << m 
                  << ", organisms=" << num_organisms << std::endl;
        
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
        
        // Read minimizer index and prepare for GPU - CORRECT APPROACH
        std::vector<std::tuple<uint64_t, uint32_t, float>> all_entries;
        
        std::cout << "Loading " << num_minimizer_hashes << " minimizer hash entries..." << std::endl;
        
        for (uint64_t i = 0; i < num_minimizer_hashes; i++) {
            uint64_t hash;
            uint32_t num_entries;
            
            in.read(reinterpret_cast<char*>(&hash), sizeof(hash));
            in.read(reinterpret_cast<char*>(&num_entries), sizeof(num_entries));
            
            for (uint32_t j = 0; j < num_entries; j++) {
                uint64_t minimizer_hash;
                uint32_t taxonomy_id;
                uint8_t uniqueness_score;
                
                in.read(reinterpret_cast<char*>(&minimizer_hash), sizeof(minimizer_hash));
                in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
                in.read(reinterpret_cast<char*>(&uniqueness_score), sizeof(uniqueness_score));
                
                // Calculate weight
                float weight = uniqueness_score / 255.0f;
                
                // Map taxonomy_id to array index
                auto it = std::find(organism_id_list.begin(), organism_id_list.end(), taxonomy_id);
                if (it != organism_id_list.end()) {
                    uint32_t array_index = std::distance(organism_id_list.begin(), it);
                    all_entries.emplace_back(minimizer_hash, array_index, weight);
                }
            }
            
            // Progress indicator for large databases
            if (i % 100000 == 0) {
                std::cout << "\rProcessed " << i << "/" << num_minimizer_hashes 
                         << " hash entries (" << all_entries.size() << " total minimizers)..." << std::flush;
            }
        }
        
        in.close();
        std::cout << "\rCompleted loading " << all_entries.size() << " minimizer entries" << std::endl;
        
        // Sort by hash for binary search
        std::cout << "Sorting minimizers for efficient lookup..." << std::endl;
        std::sort(all_entries.begin(), all_entries.end());
        
        // Transfer to GPU
        std::vector<uint64_t> hashes;
        std::vector<uint32_t> org_ids;
        std::vector<float> weights;
        
        hashes.reserve(all_entries.size());
        org_ids.reserve(all_entries.size());
        weights.reserve(all_entries.size());
        
        for (const auto& entry : all_entries) {
            hashes.push_back(std::get<0>(entry));
            org_ids.push_back(std::get<1>(entry));
            weights.push_back(std::get<2>(entry));
        }
        
        d_database_hashes = hashes;
        d_database_organism_ids = org_ids;
        d_database_weights = weights;
        
        auto load_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(load_end - load_start);
        
        size_t gpu_memory_mb = (all_entries.size() * (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(float))) / (1024 * 1024);
        std::cout << "Database loaded in " << duration.count() << " seconds" << std::endl;
        std::cout << "GPU memory used: ~" << gpu_memory_mb << " MB" << std::endl;
        std::cout << "Total minimizer entries: " << all_entries.size() << std::endl;
        
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "GPU memory: " << (total_mem - free_mem) / (1024*1024) << " MB used" << std::endl;
    }
    
    // Process a batch of paired reads and accumulate results
    void process_paired_batch(const std::vector<PairedRead>& batch,
                            size_t global_pair_offset) {
        std::cout << "Processing batch of " << batch.size() << " read pairs..." << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Extract minimizers from this batch
        std::vector<uint64_t> batch_minimizers;
        std::vector<uint32_t> batch_pair_ids;
        std::vector<uint8_t> batch_read_numbers;
        
        extract_paired_minimizers(batch, batch_minimizers, batch_pair_ids, batch_read_numbers);
        
        // Match against database and accumulate results
        match_minimizers_gpu_incremental(batch_minimizers, batch_pair_ids, 
                                       batch_read_numbers, batch.size(),
                                       global_pair_offset);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Batch processed in " << duration.count() << " ms" << std::endl;
    }
    
    // Reset accumulators for new sample
    void reset_accumulators() {
        thrust::fill(d_organism_hits.begin(), d_organism_hits.end(), 0);
        thrust::fill(d_organism_scores.begin(), d_organism_scores.end(), 0.0f);
        thrust::fill(d_organism_paired_hits.begin(), d_organism_paired_hits.end(), 0);
        thrust::fill(d_organism_concordance_scores.begin(), d_organism_concordance_scores.end(), 0.0f);
    }
    
    // Generate final profiles after processing all batches
    std::vector<OrganismProfile> generate_final_profiles() {
        return calculate_organism_profiles();
    }
    
private:
    void extract_paired_minimizers(const std::vector<PairedRead>& paired_reads,
                                  std::vector<uint64_t>& minimizers,
                                  std::vector<uint32_t>& pair_ids,
                                  std::vector<uint8_t>& read_numbers) {
        int window_size = k - m + 1;
        
        for (size_t pair_idx = 0; pair_idx < paired_reads.size(); pair_idx++) {
            const auto& pair = paired_reads[pair_idx];
            
            // Process read1
            extract_read_minimizers(pair.read1, pair_idx, 1, minimizers, pair_ids, read_numbers);
            
            // Process read2 if paired
            if (pair.is_paired && !pair.read2.empty()) {
                extract_read_minimizers(pair.read2, pair_idx, 2, minimizers, pair_ids, read_numbers);
            }
        }
    }
    
    void extract_read_minimizers(const std::string& read,
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
    
    void match_minimizers_gpu_incremental(const std::vector<uint64_t>& minimizers,
                                        const std::vector<uint32_t>& pair_ids,
                                        const std::vector<uint8_t>& read_numbers,
                                        size_t num_pairs,
                                        size_t global_pair_offset) {
        // Get current state from GPU
        thrust::host_vector<uint32_t> h_organism_hits = d_organism_hits;
        thrust::host_vector<float> h_organism_scores = d_organism_scores;
        thrust::host_vector<uint32_t> h_paired_hits = d_organism_paired_hits;
        thrust::host_vector<float> h_concordance_scores = d_organism_concordance_scores;
        
        // Get database from GPU for matching
        thrust::host_vector<uint64_t> h_db_hashes = d_database_hashes;
        thrust::host_vector<uint32_t> h_db_org_ids = d_database_organism_ids;
        thrust::host_vector<float> h_db_weights = d_database_weights;
        
        // Track paired concordance for this batch
        std::vector<std::vector<float>> pair_scores_r1(num_pairs, std::vector<float>(max_organisms, 0.0f));
        std::vector<std::vector<float>> pair_scores_r2(num_pairs, std::vector<float>(max_organisms, 0.0f));
        
        // Match minimizers
        for (size_t i = 0; i < minimizers.size(); i++) {
            uint64_t query_hash = minimizers[i];
            uint32_t pair_id = pair_ids[i];
            uint8_t read_num = read_numbers[i];
            
            // Binary search in database
            auto it = std::lower_bound(h_db_hashes.begin(), h_db_hashes.end(), query_hash);
            
            // Process all matches for this hash
            while (it != h_db_hashes.end() && *it == query_hash) {
                size_t idx = it - h_db_hashes.begin();
                uint32_t org_id = h_db_org_ids[idx];
                float weight = h_db_weights[idx];
                
                if (org_id < max_organisms) {
                    h_organism_hits[org_id]++;
                    h_organism_scores[org_id] += weight;
                    
                    // Track paired-end scoring for this batch
                    if (pair_id < num_pairs) {
                        if (read_num == 1) {
                            pair_scores_r1[pair_id][org_id] += weight;
                        } else if (read_num == 2) {
                            pair_scores_r2[pair_id][org_id] += weight;
                        }
                    }
                }
                
                ++it;
            }
        }
        
        // Calculate paired concordance for this batch
        for (size_t pair_id = 0; pair_id < num_pairs; pair_id++) {
            // Find best organism for each read
            uint32_t best_org_r1 = 0, best_org_r2 = 0;
            float best_score_r1 = 0.0f, best_score_r2 = 0.0f;
            
            for (uint32_t org = 0; org < max_organisms; org++) {
                if (pair_scores_r1[pair_id][org] > best_score_r1) {
                    best_score_r1 = pair_scores_r1[pair_id][org];
                    best_org_r1 = org;
                }
                if (pair_scores_r2[pair_id][org] > best_score_r2) {
                    best_score_r2 = pair_scores_r2[pair_id][org];
                    best_org_r2 = org;
                }
            }
            
            // If both reads map to same organism with good scores
            if (best_org_r1 == best_org_r2 && best_score_r1 > 0.1f && best_score_r2 > 0.1f) {
                h_paired_hits[best_org_r1]++;
                h_concordance_scores[best_org_r1] += best_score_r1 + best_score_r2;
            }
        }
        
        // Copy accumulated results back to GPU
        d_organism_hits = h_organism_hits;
        d_organism_scores = h_organism_scores;
        d_organism_paired_hits = h_paired_hits;
        d_organism_concordance_scores = h_concordance_scores;
    }
    
    std::vector<OrganismProfile> calculate_organism_profiles() {
        // Copy results from GPU
        thrust::host_vector<uint32_t> h_organism_hits = d_organism_hits;
        thrust::host_vector<float> h_organism_scores = d_organism_scores;
        thrust::host_vector<uint32_t> h_organism_paired_hits = d_organism_paired_hits;
        thrust::host_vector<float> h_organism_concordance_scores = d_organism_concordance_scores;
        
        std::vector<OrganismProfile> profiles;
        
        // Calculate total abundance
        float total_abundance = 0.0f;
        for (size_t i = 0; i < organism_id_list.size() && i < max_organisms; i++) {
            total_abundance += h_organism_scores[i];
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
            profile.single_hits = hits - (paired_hits * 2);
            
            profile.abundance = (total_abundance > 0) ? score / total_abundance : 0.0f;
            
            uint32_t expected_minimizers = organism_expected_minimizers[tax_id];
            profile.coverage_breadth = (expected_minimizers > 0) ? 
                                     std::min(1.0f, (float)hits / expected_minimizers) : 0.0f;
            
            profile.paired_concordance = (paired_hits > 0) ? 
                                       concordance_score / (paired_hits * 2) : 0.0f;
            
            // Simple confidence score
            profile.confidence_score = std::min(1.0f, (float)(std::log10(hits + 1) / 3.0));
            
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
};

class AdaptivePairedEndProfiler {
private:
    // Database metadata
    std::unordered_map<uint32_t, std::string> organism_names;
    std::unordered_map<uint32_t, std::string> organism_taxonomy;
    std::unordered_map<uint32_t, uint64_t> organism_genome_sizes;
    std::unordered_map<uint32_t, uint32_t> organism_expected_minimizers;
    std::vector<uint32_t> organism_id_list;
    
    // Analysis results
    DatabaseStats db_stats;
    bool use_streaming_mode = false;
    
    // Direct GPU profiler
    std::unique_ptr<DirectGPUProfiler> gpu_profiler;
    
    int k = 35, m = 31;
    
    // Batch processing parameters
    static const int DEFAULT_BATCH_SIZE = 50000;  // Number of read pairs per batch
    int batch_size = DEFAULT_BATCH_SIZE;
    
public:
    AdaptivePairedEndProfiler() {
        std::cout << "Adaptive Paired-End GPU Profiler (Batched Version)" << std::endl;
        std::cout << "Batch processing with direct GPU profiling" << std::endl;
    }
    
    void set_batch_size(int size) {
        batch_size = size;
        std::cout << "Batch size set to: " << batch_size << " read pairs" << std::endl;
    }
    
    void load_and_analyze_database(const std::string& db_path) {
        std::cout << "Analyzing database: " << db_path << std::endl;
        
        // Analyze database size
        db_stats = analyze_database_requirements(db_path);
        
        // Decide on processing strategy
        use_streaming_mode = db_stats.needs_streaming;
        
        std::cout << "\n=== DATABASE ANALYSIS RESULTS ===" << std::endl;
        std::cout << "Estimated GPU memory needed: " << db_stats.estimated_gpu_memory_mb << " MB" << std::endl;
        std::cout << "Number of organisms: " << db_stats.num_organisms << std::endl;
        std::cout << "Total minimizer entries: " << db_stats.total_minimizer_entries << std::endl;
        
        if (use_streaming_mode) {
            std::cout << "🔄 STREAMING MODE REQUIRED" << std::endl;
            std::cout << "Database too large for direct GPU loading" << std::endl;
            throw std::runtime_error("Streaming mode not implemented in this version");
        } else {
            std::cout << "⚡ DIRECT GPU MODE SELECTED" << std::endl;
            std::cout << "Database fits in GPU memory - maximum performance" << std::endl;
            initialize_direct_mode(db_path);
        }
        
        // Load organism metadata
        load_organism_metadata(db_path);
    }
    
private:
    DatabaseStats analyze_database_requirements(const std::string& db_path) {
        DatabaseStats stats = {};
        
        std::ifstream file(db_path, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open database file for analysis");
        }
        
        // Get file size
        uint64_t file_size_mb = file.tellg() / (1024 * 1024);
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
        
        if (magic != 0x4D494E49 && magic != 0x494E494D) {
            throw std::runtime_error("Invalid database format");
        }
        
        stats.total_minimizer_entries = num_minimizer_hashes;
        
        // Estimate memory requirements
        stats.estimated_gpu_memory_mb = file_size_mb * 2;  // Conservative estimate
        
        // Decision: use streaming if estimated memory > 30GB
        size_t available_memory_mb = 30 * 1024;
        stats.needs_streaming = stats.estimated_gpu_memory_mb > available_memory_mb;
        
        file.close();
        return stats;
    }
    
    void initialize_direct_mode(const std::string& db_path) {
        std::cout << "Initializing direct GPU mode..." << std::endl;
        
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "GPU memory available: " << free_mem / (1024*1024*1024) << " GB" << std::endl;
        
        gpu_profiler = std::make_unique<DirectGPUProfiler>();
        gpu_profiler->load_minimizer_database(db_path);
        
        std::cout << "Direct GPU mode ready" << std::endl;
    }
    
    void load_organism_metadata(const std::string& db_path) {
        std::cout << "Loading organism metadata..." << std::endl;
        
        std::ifstream in(db_path, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open database for metadata loading");
        }
        
        // Skip header
        uint32_t magic, version, k_size, m_size, num_organisms;
        uint64_t num_minimizer_hashes;
        
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        in.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        in.read(reinterpret_cast<char*>(&num_organisms), sizeof(num_organisms));
        in.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        k = k_size;
        m = m_size;
        
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
            uint16_t name_length, taxonomy_length;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            std::string name(name_length, '\0');
            in.read(&name[0], name_length);
            
            in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            std::string taxonomy_path(taxonomy_length, '\0');
            in.read(&taxonomy_path[0], taxonomy_length);
            
            // Store metadata
            organism_names[taxonomy_id] = name;
            organism_taxonomy[taxonomy_id] = taxonomy_path;
            organism_genome_sizes[taxonomy_id] = genome_size;
            organism_expected_minimizers[taxonomy_id] = minimizer_count;
            organism_id_list.push_back(taxonomy_id);
        }
        
        in.close();
        std::cout << "Loaded metadata for " << organism_names.size() << " organisms" << std::endl;
    }
    
public:
    // Batch reading function similar to resistance pipeline
    bool readBatch(gzFile gz_r1, gzFile gz_r2, 
                  std::vector<FastqRecord>& batch_r1,
                  std::vector<FastqRecord>& batch_r2, 
                  int max_size) {
        char buffer[1024];
        
        for (int i = 0; i < max_size; i++) {
            FastqRecord rec1, rec2;
            
            // Read R1
            if (gzgets(gz_r1, buffer, 1024) == NULL) return i > 0;
            rec1.header = std::string(buffer);
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
            rec1.sequence = std::string(buffer);
            rec1.sequence.pop_back(); // Remove newline
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false; // +
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
            rec1.quality = std::string(buffer);
            rec1.quality.pop_back();
            
            // Read R2
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.header = std::string(buffer);
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.sequence = std::string(buffer);
            rec2.sequence.pop_back();
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false; // +
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.quality = std::string(buffer);
            rec2.quality.pop_back();
            
            batch_r1.push_back(rec1);
            batch_r2.push_back(rec2);
        }
        
        return true;
    }
    
    // NEW: Process paired FASTQ files in batches
    std::vector<OrganismProfile> profile_paired_fastq_batched(const std::string& r1_path, 
                                                             const std::string& r2_path = "") {
        std::cout << "\n=== STARTING BATCHED PAIRED-END COMMUNITY PROFILING ===" << std::endl;
        std::cout << "Processing files: " << r1_path << std::endl;
        if (!r2_path.empty()) {
            std::cout << "                  " << r2_path << std::endl;
        }
        std::cout << "Batch size: " << batch_size << " read pairs" << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Open FASTQ files
        gzFile gz_r1 = gzopen(r1_path.c_str(), "r");
        gzFile gz_r2 = r2_path.empty() ? nullptr : gzopen(r2_path.c_str(), "r");
        
        if (!gz_r1 || (!r2_path.empty() && !gz_r2)) {
            throw std::runtime_error("Failed to open input files");
        }
        
        // Reset GPU profiler accumulators
        gpu_profiler->reset_accumulators();
        
        size_t total_pairs_processed = 0;
        int batch_num = 0;
        
        // Process in batches
        while (true) {
            std::vector<FastqRecord> batch_r1, batch_r2;
            
            bool has_data;
            if (gz_r2) {
                has_data = readBatch(gz_r1, gz_r2, batch_r1, batch_r2, batch_size);
            } else {
                // Handle single file (interleaved) - not implemented yet
                throw std::runtime_error("Interleaved FASTQ not yet supported");
            }
            
            if (!has_data || batch_r1.empty()) break;
            
            // Convert FastqRecords to PairedReads
            std::vector<PairedRead> paired_batch;
            for (size_t i = 0; i < batch_r1.size(); i++) {
                std::string read_id = extract_read_id(batch_r1[i].header);
                PairedRead pair(batch_r1[i].sequence, 
                              i < batch_r2.size() ? batch_r2[i].sequence : "",
                              read_id);
                paired_batch.push_back(pair);
            }
            
            // Process this batch
            batch_num++;
            std::cout << "\nProcessing batch " << batch_num << " (" 
                     << paired_batch.size() << " pairs)..." << std::endl;
            
            gpu_profiler->process_paired_batch(paired_batch, total_pairs_processed);
            
            total_pairs_processed += paired_batch.size();
            
            // Progress update
            if (total_pairs_processed % 100000 == 0) {
                std::cout << "Processed " << total_pairs_processed << " read pairs..." << std::endl;
            }
        }
        
        // Close files
        gzclose(gz_r1);
        if (gz_r2) gzclose(gz_r2);
        
        // Generate final profiles
        auto profiles = gpu_profiler->generate_final_profiles();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\nBatched profiling completed:" << std::endl;
        std::cout << "- Total read pairs: " << total_pairs_processed << std::endl;
        std::cout << "- Total batches: " << batch_num << std::endl;
        std::cout << "- Total time: " << duration.count() << " seconds" << std::endl;
        std::cout << "- Performance: " << (total_pairs_processed / duration.count()) 
                 << " pairs/second" << std::endl;
        
        return profiles;
    }
    
    // DEPRECATED: Old method that loads all reads at once
    std::vector<PairedRead> load_paired_fastq(const std::string& r1_path, 
                                             const std::string& r2_path = "") {
        std::cerr << "WARNING: load_paired_fastq is deprecated and may cause memory issues!" << std::endl;
        std::cerr << "Use profile_paired_fastq_batched() instead." << std::endl;
        throw std::runtime_error("Method deprecated - use batched processing");
    }
    
    std::string extract_read_id(const std::string& header) {
        if (header.empty() || header[0] != '@') {
            return "";
        }
        size_t space_pos = header.find(' ');
        if (space_pos != std::string::npos) {
            return header.substr(1, space_pos - 1);
        } else {
            size_t newline_pos = header.find('\n');
            if (newline_pos != std::string::npos) {
                return header.substr(1, newline_pos - 1);
            }
            return header.substr(1);
        }
    }
    
    void print_results(const std::vector<OrganismProfile>& profiles) {
        std::cout << "\n=== ADAPTIVE PAIRED-END PROFILING RESULTS ===" << std::endl;
        std::cout << "Processing mode: BATCHED GPU (memory efficient)" << std::endl;
        
        if (profiles.empty()) {
            std::cout << "No organisms detected above significance thresholds" << std::endl;
            return;
        }
        
        std::cout << "\nTop organisms detected:" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        std::cout << std::left << std::setw(40) << "Organism Name" 
                  << std::setw(10) << "Abundance" 
                  << std::setw(8) << "Total"
                  << std::setw(8) << "Paired"
                  << std::setw(12) << "Concordance"
                  << std::setw(10) << "Confidence" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        
        for (size_t i = 0; i < std::min(size_t(15), profiles.size()); i++) {
            const auto& profile = profiles[i];
            
            std::string display_name = profile.name;
            if (display_name.length() > 35) {
                display_name = display_name.substr(0, 32) + "...";
            }
            
            std::cout << std::left << std::setw(40) << display_name
                      << std::fixed << std::setprecision(3) << std::setw(10) << (profile.abundance * 100)
                      << std::setw(8) << profile.total_hits
                      << std::setw(8) << profile.paired_hits
                      << std::fixed << std::setprecision(2) << std::setw(12) << (profile.paired_concordance * 100)
                      << std::fixed << std::setprecision(3) << std::setw(10) << profile.confidence_score
                      << std::endl;
        }
    }
    
    void write_results(const std::vector<OrganismProfile>& profiles, const std::string& output_prefix) {
        std::string output_file = output_prefix + "_abundance.tsv";
        std::ofstream out(output_file);
        
        out << "taxonomy_id\torganism_name\ttaxonomy_path\trelative_abundance\t"
            << "coverage_breadth\ttotal_hits\tpaired_hits\tsingle_hits\t"
            << "paired_concordance\tconfidence_score\n";
        
        for (const auto& profile : profiles) {
            out << profile.taxonomy_id << "\t"
                << profile.name << "\t"
                << profile.taxonomy_path << "\t"
                << std::scientific << profile.abundance << "\t"
                << std::fixed << std::setprecision(6) << profile.coverage_breadth << "\t"
                << profile.total_hits << "\t"
                << profile.paired_hits << "\t"
                << profile.single_hits << "\t"
                << std::fixed << std::setprecision(4) << profile.paired_concordance << "\t"
                << std::fixed << std::setprecision(4) << profile.confidence_score << "\n";
        }
        
        out.close();
        std::cout << "Results written to: " << output_file << std::endl;
    }
};

// Enhanced BatchAdaptivePairedEndProfiler with batch processing capability
class BatchAdaptivePairedEndProfiler : public AdaptivePairedEndProfiler {
private:
    bool database_loaded = false;
    std::string current_database_path;
    
public:
    BatchAdaptivePairedEndProfiler() {
        std::cout << "Batch Adaptive Paired-End GPU Profiler" << std::endl;
        std::cout << "Optimized for processing multiple samples with shared database" << std::endl;
    }
    
    // Load database once - keep in GPU memory for all subsequent samples
    void load_database_once(const std::string& db_path) {
        if (database_loaded && current_database_path == db_path) {
            std::cout << "Database already loaded: " << db_path << std::endl;
            return;
        }
        
        std::cout << "\n=== LOADING DATABASE FOR BATCH PROCESSING ===" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Load and analyze database using parent class method
        load_and_analyze_database(db_path);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        database_loaded = true;
        current_database_path = db_path;
        
        std::cout << "Database loaded and resident in GPU memory (" << duration.count() << "s)" << std::endl;
        std::cout << "Ready for rapid batch processing!" << std::endl;
    }
    
    // Process a single sample (database already loaded)
    std::vector<OrganismProfile> process_single_sample(
        const std::string& sample_name,
        const std::string& r1_path, 
        const std::string& r2_path = ""
    ) {
        if (!database_loaded) {
            throw std::runtime_error("Database not loaded. Call load_database_once() first.");
        }
        
        std::cout << "\n=== PROCESSING SAMPLE: " << sample_name << " ===" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Use batched processing
        auto profiles = profile_paired_fastq_batched(r1_path, r2_path);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Sample processed in " << duration.count() << " ms" << std::endl;
        
        return profiles;
    }
    
    // Process multiple samples in batch mode
    void process_sample_batch(
        const std::string& database_path,
        const std::vector<std::tuple<std::string, std::string, std::string>>& samples,
        const std::string& output_dir = "batch_results"
    ) {
        std::cout << "\n=== BATCH PROCESSING " << samples.size() << " SAMPLES ===" << std::endl;
        
        // Load database once
        load_database_once(database_path);
        
        // Create output directory
        std::filesystem::create_directories(output_dir);
        
        auto batch_start = std::chrono::high_resolution_clock::now();
        
        // Process each sample
        for (size_t i = 0; i < samples.size(); i++) {
            const auto& [sample_name, r1_path, r2_path] = samples[i];
            
            std::cout << "\n--- Sample " << (i + 1) << "/" << samples.size() << " ---" << std::endl;
            
            try {
                // Process this sample
                auto profiles = process_single_sample(sample_name, r1_path, r2_path);
                
                // Write results for this sample
                std::string output_prefix = output_dir + "/" + sample_name + "/" + sample_name;
                std::filesystem::create_directories(output_dir + "/" + sample_name);
                write_results(profiles, output_prefix);
                
                // Print brief summary
                if (!profiles.empty()) {
                    std::cout << "Top organism: " << profiles[0].name 
                             << " (" << std::fixed << std::setprecision(2) 
                             << profiles[0].abundance * 100 << "%)" << std::endl;
                }
                
            } catch (const std::exception& e) {
                std::cerr << "Error processing sample " << sample_name << ": " << e.what() << std::endl;
                continue;
            }
        }
        
        auto batch_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(batch_end - batch_start);
        
        std::cout << "\n=== BATCH PROCESSING COMPLETE ===" << std::endl;
        std::cout << samples.size() << " files processed in " << total_duration.count() << " seconds" << std::endl;
        std::cout << "Average time per sample: " << (total_duration.count() / (float)samples.size()) << " seconds" << std::endl;
        std::cout << "Results written to: " << output_dir << "/" << std::endl;
        
        // Generate batch summary report
        generate_batch_summary(samples, output_dir);
    }
    
    // Process multiple samples using CSV format
    void process_sample_batch_csv(
        const std::string& database_path,
        const std::string& csv_path,
        const std::string& output_dir = "batch_results",
        bool stop_on_error = false
    ) {
        std::cout << "\n=== BATCH PROCESSING FROM CSV ===" << std::endl;
        
        // Create batch processor
        BioGPU::BatchProcessor batch_processor(output_dir, true);
        
        // Load samples from CSV
        if (!batch_processor.loadSamples(csv_path)) {
            throw std::runtime_error("Failed to load samples from CSV: " + csv_path);
        }
        
        std::cout << "Loaded " << batch_processor.getParser().getSampleCount() << " samples from CSV" << std::endl;
        batch_processor.getParser().printSummary();
        
        // Load database once
        load_database_once(database_path);
        
        // Process batch with lambda function
        auto process_func = [&](const BioGPU::SampleInfo& sample, const std::string& output_prefix) -> int {
            try {
                // Process this sample
                auto profiles = process_single_sample(sample.sample_name, 
                                                    sample.read1_path, 
                                                    sample.read2_path);
                
                // Write results
                write_results(profiles, output_prefix);
                
                // Write summary
                std::string summary_file = output_prefix + "_summary.tsv";
                std::ofstream summary(summary_file);
                summary << "sample_name\ttotal_organisms_detected\ttop_organism\ttop_abundance\n";
                
                if (!profiles.empty()) {
                    summary << sample.sample_name << "\t" 
                           << profiles.size() << "\t"
                           << profiles[0].name << "\t"
                           << std::scientific << profiles[0].abundance << "\n";
                }
                summary.close();
                
                return 0;  // Success
            } catch (const std::exception& e) {
                std::cerr << "Error processing sample " << sample.sample_name 
                         << ": " << e.what() << std::endl;
                return stop_on_error ? -1 : 1;  // Error
            }
        };
        
        // Process all samples
        bool success = batch_processor.processBatch(process_func);
        
        if (success) {
            std::cout << "\n✓ Batch processing completed successfully" << std::endl;
        } else {
            std::cout << "\n✗ Batch processing completed with errors" << std::endl;
        }
    }
    
private:
    void generate_batch_summary(const std::vector<std::tuple<std::string, std::string, std::string>>& samples,
                               const std::string& output_dir) {
        std::string summary_file = output_dir + "/batch_summary.tsv";
        std::ofstream out(summary_file);
        
        out << "sample_name\ttotal_organisms_detected\ttop_organism\ttop_abundance\n";
        
        for (const auto& [sample_name, r1_path, r2_path] : samples) {
            std::string result_file = output_dir + "/" + sample_name + "/" + sample_name + "_abundance.tsv";
            std::ifstream in(result_file);
            
            if (in.is_open()) {
                std::string line;
                std::getline(in, line);  // Skip header
                
                int count = 0;
                std::string top_organism = "";
                float top_abundance = 0.0f;
                
                while (std::getline(in, line)) {
                    count++;
                    if (count == 1) {
                        std::istringstream iss(line);
                        std::string tax_id, name, tax_path;
                        float abundance;
                        
                        std::getline(iss, tax_id, '\t');
                        std::getline(iss, name, '\t');
                        std::getline(iss, tax_path, '\t');
                        iss >> abundance;
                        
                        top_organism = name;
                        top_abundance = abundance;
                    }
                }
                
                out << sample_name << "\t" << count << "\t" 
                    << top_organism << "\t" << std::scientific << top_abundance << "\n";
                
                in.close();
            }
        }
        
        out.close();
        std::cout << "Batch summary written to: " << summary_file << std::endl;
    }
};

// Main function for testing
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <database> <mode> [options]" << std::endl;
        std::cerr << "Modes:" << std::endl;
        std::cerr << "  single <r1.fastq.gz> <r2.fastq.gz> <output_prefix>" << std::endl;
        std::cerr << "  batch <csv_file> <output_dir>" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string mode = argv[2];
    
    try {
        if (mode == "single") {
            if (argc < 6) {
                std::cerr << "Single mode requires: <r1.fastq.gz> <r2.fastq.gz> <output_prefix>" << std::endl;
                return 1;
            }
            
            std::string r1_path = argv[3];
            std::string r2_path = argv[4];
            std::string output_prefix = argv[5];
            
            AdaptivePairedEndProfiler profiler;
            profiler.load_and_analyze_database(database_path);
            
            auto profiles = profiler.profile_paired_fastq_batched(r1_path, r2_path);
            
            profiler.print_results(profiles);
            profiler.write_results(profiles, output_prefix);
            
        } else if (mode == "batch") {
            if (argc < 5) {
                std::cerr << "Batch mode requires: <csv_file> <output_dir>" << std::endl;
                return 1;
            }
            
            std::string csv_path = argv[3];
            std::string output_dir = argv[4];
            
            BatchAdaptivePairedEndProfiler batch_profiler;
            batch_profiler.process_sample_batch_csv(database_path, csv_path, output_dir);
            
        } else {
            std::cerr << "Unknown mode: " << mode << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}