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
    
    // GPU results storage
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
        
        // Read minimizer index and prepare for GPU
        std::vector<std::tuple<uint64_t, uint32_t, float>> all_entries;
        
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
        }
        
        in.close();
        
        // Sort by hash for binary search
        std::sort(all_entries.begin(), all_entries.end());
        
        // Transfer to GPU
        std::vector<uint64_t> hashes;
        std::vector<uint32_t> org_ids;
        std::vector<float> weights;
        
        for (const auto& entry : all_entries) {
            hashes.push_back(std::get<0>(entry));
            org_ids.push_back(std::get<1>(entry));
            weights.push_back(std::get<2>(entry));
        }
        
        d_database_hashes = hashes;
        d_database_organism_ids = org_ids;
        d_database_weights = weights;
        
        std::cout << "Loaded " << all_entries.size() << " minimizer entries to GPU" << std::endl;
        
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "GPU memory: " << (total_mem - free_mem) / (1024*1024) << " MB used" << std::endl;
    }
    
    std::vector<OrganismProfile> profile_paired_reads_directly(const std::vector<PairedRead>& paired_reads) {
        std::cout << "Processing " << paired_reads.size() << " paired reads directly on GPU..." << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Extract minimizers from paired reads
        std::vector<uint64_t> all_minimizers;
        std::vector<uint32_t> minimizer_pair_ids;
        std::vector<uint8_t> minimizer_read_numbers;
        
        extract_paired_minimizers(paired_reads, all_minimizers, minimizer_pair_ids, minimizer_read_numbers);
        
        // Match against database
        match_minimizers_gpu(all_minimizers, minimizer_pair_ids, minimizer_read_numbers, paired_reads.size());
        
        // Generate profiles
        auto profiles = calculate_organism_profiles();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Direct GPU profiling completed in " << duration.count() << " ms" << std::endl;
        
        return profiles;
    }
    
private:
    void extract_paired_minimizers(const std::vector<PairedRead>& paired_reads,
                                  std::vector<uint64_t>& minimizers,
                                  std::vector<uint32_t>& pair_ids,
                                  std::vector<uint8_t>& read_numbers) {
        std::cout << "Extracting minimizers from paired reads..." << std::endl;
        
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
        
        std::cout << "Extracted " << minimizers.size() << " minimizers" << std::endl;
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
    
    void match_minimizers_gpu(const std::vector<uint64_t>& minimizers,
                             const std::vector<uint32_t>& pair_ids,
                             const std::vector<uint8_t>& read_numbers,
                             size_t num_pairs) {
        std::cout << "Matching " << minimizers.size() << " minimizers against database..." << std::endl;
        
        // Reset counters
        thrust::fill(d_organism_hits.begin(), d_organism_hits.end(), 0);
        thrust::fill(d_organism_scores.begin(), d_organism_scores.end(), 0.0f);
        thrust::fill(d_organism_paired_hits.begin(), d_organism_paired_hits.end(), 0);
        thrust::fill(d_organism_concordance_scores.begin(), d_organism_concordance_scores.end(), 0.0f);
        
        // Simple CPU matching for now (can be moved to GPU later)
        thrust::host_vector<uint64_t> h_db_hashes = d_database_hashes;
        thrust::host_vector<uint32_t> h_db_org_ids = d_database_organism_ids;
        thrust::host_vector<float> h_db_weights = d_database_weights;
        
        thrust::host_vector<uint32_t> h_organism_hits(max_organisms, 0);
        thrust::host_vector<float> h_organism_scores(max_organisms, 0.0f);
        
        // Track paired concordance
        std::vector<std::vector<float>> pair_scores_r1(num_pairs, std::vector<float>(max_organisms, 0.0f));
        std::vector<std::vector<float>> pair_scores_r2(num_pairs, std::vector<float>(max_organisms, 0.0f));
        
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
                    
                    // Track paired-end scoring
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
        
        // Calculate paired concordance
        thrust::host_vector<uint32_t> h_paired_hits(max_organisms, 0);
        thrust::host_vector<float> h_concordance_scores(max_organisms, 0.0f);
        
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
        
        // Copy results back to GPU
        d_organism_hits = h_organism_hits;
        d_organism_scores = h_organism_scores;
        d_organism_paired_hits = h_paired_hits;
        d_organism_concordance_scores = h_concordance_scores;
        
        std::cout << "Matching completed" << std::endl;
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
    
public:
    AdaptivePairedEndProfiler() {
        std::cout << "Adaptive Paired-End GPU Profiler" << std::endl;
        std::cout << "Direct processing without temporary files" << std::endl;
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
    // FIXED: Direct FASTQ loading with gzip support
    std::vector<PairedRead> load_paired_fastq(const std::string& r1_path, 
                                             const std::string& r2_path = "") {
        std::cout << "Loading paired-end reads directly..." << std::endl;
        
        if (r2_path.empty()) {
            return load_interleaved_fastq(r1_path);
        } else {
            return load_separate_fastq_files(r1_path, r2_path);
        }
    }
    
private:
    std::vector<PairedRead> load_separate_fastq_files(const std::string& r1_path, 
                                                     const std::string& r2_path) {
        std::cout << "Loading paired-end reads from " << r1_path << " and " << r2_path << std::endl;
        std::vector<PairedRead> paired_reads;
        
        // Check if files are gzipped
        bool r1_gzipped = (r1_path.substr(r1_path.length() - 3) == ".gz");
        bool r2_gzipped = (r2_path.substr(r2_path.length() - 3) == ".gz");
        
        if (r1_gzipped || r2_gzipped) {
            return load_gzipped_paired_files(r1_path, r2_path);
        }
        
        std::ifstream r1_file(r1_path);
        std::ifstream r2_file(r2_path);
        
        if (!r1_file.is_open() || !r2_file.is_open()) {
            throw std::runtime_error("Cannot open paired-end FASTQ files");
        }
        
        std::string r1_line, r2_line;
        std::string r1_header, r1_seq, r2_header, r2_seq;
        int line_count = 0;
        
        while (std::getline(r1_file, r1_line) && std::getline(r2_file, r2_line)) {
            int line_type = (line_count % 4);
            
            if (line_type == 0) {
                // Header lines
                r1_header = r1_line;
                r2_header = r2_line;
            } else if (line_type == 1) {
                // Sequence lines
                r1_seq = r1_line;
                r2_seq = r2_line;
                
                if (!r1_seq.empty() && !r2_seq.empty() && 
                    r1_seq.length() >= k && r2_seq.length() >= k) {
                    
                    std::string read_id = extract_read_id(r1_header);
                    PairedRead pair(r1_seq, r2_seq, read_id);
                    paired_reads.push_back(pair);
                }
            }
            // Skip quality lines (line_type 2 and 3)
            
            line_count++;
            
            if (paired_reads.size() % 50000 == 0) {
                std::cout << "\rLoaded " << paired_reads.size() << " paired reads..." << std::flush;
            }
            
            // Limit for testing
            if (paired_reads.size() >= 500000) {
                std::cout << "\nLimited to " << paired_reads.size() << " pairs for testing" << std::endl;
                break;
            }
        }
        
        std::cout << "\nLoaded " << paired_reads.size() << " paired-end reads" << std::endl;
        return paired_reads;
    }
    
    std::vector<PairedRead> load_gzipped_paired_files(const std::string& r1_path, 
                                                     const std::string& r2_path) {
        std::cout << "Loading gzipped paired-end files..." << std::endl;
        std::vector<PairedRead> paired_reads;
        
        gzFile r1_file = gzopen(r1_path.c_str(), "rb");
        gzFile r2_file = gzopen(r2_path.c_str(), "rb");
        
        if (!r1_file || !r2_file) {
            throw std::runtime_error("Cannot open gzipped FASTQ files");
        }
        
        char r1_buffer[1024], r2_buffer[1024];
        std::string r1_header, r1_seq, r2_header, r2_seq;
        int line_count = 0;
        
        while (gzgets(r1_file, r1_buffer, sizeof(r1_buffer)) && 
               gzgets(r2_file, r2_buffer, sizeof(r2_buffer))) {
            
            // Remove newlines
            std::string r1_line(r1_buffer);
            std::string r2_line(r2_buffer);
            if (!r1_line.empty() && r1_line.back() == '\n') r1_line.pop_back();
            if (!r2_line.empty() && r2_line.back() == '\n') r2_line.pop_back();
            
            int line_type = (line_count % 4);
            
            if (line_type == 0) {
                // Header lines
                r1_header = r1_line;
                r2_header = r2_line;
            } else if (line_type == 1) {
                // Sequence lines
                r1_seq = r1_line;
                r2_seq = r2_line;
                
                if (!r1_seq.empty() && !r2_seq.empty() && 
                    r1_seq.length() >= k && r2_seq.length() >= k) {
                    
                    std::string read_id = extract_read_id(r1_header);
                    PairedRead pair(r1_seq, r2_seq, read_id);
                    paired_reads.push_back(pair);
                }
            }
            
            line_count++;
            
            if (paired_reads.size() % 50000 == 0) {
                std::cout << "\rLoaded " << paired_reads.size() << " paired reads..." << std::flush;
            }
            
            // Limit for testing
            if (paired_reads.size() >= 500000) {
                std::cout << "\nLimited to " << paired_reads.size() << " pairs for testing" << std::endl;
                break;
            }
        }
        
        gzclose(r1_file);
        gzclose(r2_file);
        
        std::cout << "\nLoaded " << paired_reads.size() << " paired-end reads from gzipped files" << std::endl;
        return paired_reads;
    }
    
    std::vector<PairedRead> load_interleaved_fastq(const std::string& fastq_path) {
        std::cout << "Loading interleaved paired-end reads from " << fastq_path << std::endl;
        return std::vector<PairedRead>();  // Placeholder
    }
    
    std::string extract_read_id(const std::string& header) {
        if (header.empty() || header[0] != '@') {
            return "";
        }
        size_t space_pos = header.find(' ');
        if (space_pos != std::string::npos) {
            return header.substr(1, space_pos - 1);
        } else {
            return header.substr(1);
        }
    }
    
public:
    std::vector<OrganismProfile> profile_paired_end_community(const std::vector<PairedRead>& paired_reads) {
        std::cout << "\n=== STARTING PAIRED-END COMMUNITY PROFILING ===" << std::endl;
        std::cout << "Processing " << paired_reads.size() << " read pairs" << std::endl;
        std::cout << "Mode: DIRECT GPU (no temporary files)" << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Process directly with GPU profiler
        auto profiles = gpu_profiler->profile_paired_reads_directly(paired_reads);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Paired-end profiling completed in " << duration.count() << " ms" << std::endl;
        
        return profiles;
    }
    
    void print_results(const std::vector<OrganismProfile>& profiles) {
        std::cout << "\n=== ADAPTIVE PAIRED-END PROFILING RESULTS ===" << std::endl;
        std::cout << "Processing mode: DIRECT GPU (no temp files)" << std::endl;
        
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
// Keeps database in GPU memory for processing multiple samples
class BatchAdaptivePairedEndProfiler : public AdaptivePairedEndProfiler {
private:
    std::unique_ptr<DirectGPUProfiler> batch_gpu_profiler;
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
        
        // Load paired-end reads (only step that changes per sample)
        auto paired_reads = load_paired_fastq(r1_path, r2_path);
        std::cout << "Loaded " << paired_reads.size() << " paired reads" << std::endl;
        
        // Profile using pre-loaded database (this should be very fast)
        auto profiles = profile_paired_end_community(paired_reads);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Sample processed in " << duration.count() << " ms (database already loaded)" << std::endl;
        
        return profiles;
    }
    
    // Process multiple samples in batch mode (text format - for backward compatibility)
    void process_sample_batch(
        const std::string& database_path,
        const std::vector<std::tuple<std::string, std::string, std::string>>& samples, // name, r1, r2
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
                // Process this sample (fast - database already loaded)
                auto profiles = process_single_sample(sample_name, r1_path, r2_path);
                
                // Write results for this sample
                std::string output_prefix = output_dir + "/" + sample_name;
                write_sample_results(profiles, output_prefix, sample_name);
                
                // Print brief summary
                print_sample_summary(profiles, sample_name);
                
            } catch (const std::exception& e) {
                std::cerr << "Error processing sample " << sample_name << ": " << e.what() << std::endl;
                continue;
            }
        }
        
        auto batch_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(batch_end - batch_start);
        
        std::cout << "\n=== BATCH PROCESSING COMPLETE ===" << std::endl;
        std::cout << "Processed " << samples.size() << " samples in " << total_duration.count() << " seconds" << std::endl;
        std::cout << "Average time per sample: " << (total_duration.count() / (float)samples.size()) << " seconds" << std::endl;
        std::cout << "Results written to: " << output_dir << "/" << std::endl;
        
        // Generate batch summary report
        generate_batch_summary(samples, output_dir);
    }
    
    // Process multiple samples using CSV format (consistent with resistance pipeline)
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
                
                // Write results for this sample
                write_sample_results(profiles, output_prefix, sample.sample_name);
                
                // Print brief summary
                print_sample_summary(profiles, sample.sample_name);
                
                return 0; // Success
            } catch (const std::exception& e) {
                std::cerr << "Error processing sample " << sample.sample_name 
                         << ": " << e.what() << std::endl;
                return 1; // Failure
            }
        };
        
        // Run batch processing
        int failed_count = batch_processor.processBatch(process_func, stop_on_error);
        
        if (failed_count > 0) {
            std::cout << "\nWARNING: " << failed_count << " samples failed processing" << std::endl;
        }
        
        // Generate overall batch summary
        generate_batch_summary_csv(batch_processor.getParser(), output_dir);
    }
    
    // Load samples from a file list
    static std::vector<std::tuple<std::string, std::string, std::string>> 
    load_sample_list(const std::string& sample_list_file) {
        std::vector<std::tuple<std::string, std::string, std::string>> samples;
        
        std::ifstream file(sample_list_file);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open sample list file: " + sample_list_file);
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;  // Skip empty lines and comments
            
            std::istringstream iss(line);
            std::string sample_name, r1_path, r2_path;
            
            if (iss >> sample_name >> r1_path >> r2_path) {
                samples.emplace_back(sample_name, r1_path, r2_path);
            } else if (iss.clear(), iss.str(line), iss >> sample_name >> r1_path) {
                // Single-end sample
                samples.emplace_back(sample_name, r1_path, "");
            }
        }
        
        std::cout << "Loaded " << samples.size() << " samples from " << sample_list_file << std::endl;
        return samples;
    }
    
private:
    void write_sample_results(const std::vector<OrganismProfile>& profiles, 
                             const std::string& output_prefix,
                             const std::string& sample_name) {
        
        // Write abundance table
        std::string abundance_file = output_prefix + "_abundance.tsv";
        std::ofstream out(abundance_file);
        
        out << "sample_name\ttaxonomy_id\torganism_name\ttaxonomy_path\t"
            << "relative_abundance\tcoverage_breadth\ttotal_hits\tpaired_hits\t"
            << "single_hits\tpaired_concordance\tconfidence_score\n";
        
        for (const auto& profile : profiles) {
            out << sample_name << "\t"
                << profile.taxonomy_id << "\t"
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
    }
    
    void print_sample_summary(const std::vector<OrganismProfile>& profiles,
                             const std::string& sample_name) {
        if (profiles.empty()) {
            std::cout << "No organisms detected in " << sample_name << std::endl;
            return;
        }
        
        std::cout << "Top organisms in " << sample_name << ":" << std::endl;
        for (size_t i = 0; i < std::min(size_t(5), profiles.size()); i++) {
            const auto& profile = profiles[i];
            std::cout << "  " << (i+1) << ". " << profile.name 
                      << " (" << std::fixed << std::setprecision(1) << (profile.abundance * 100) << "%)" << std::endl;
        }
    }
    
    void generate_batch_summary(
        const std::vector<std::tuple<std::string, std::string, std::string>>& samples,
        const std::string& output_dir
    ) {
        std::cout << "Generating batch summary report..." << std::endl;
        
        std::string summary_file = output_dir + "/batch_summary.tsv";
        std::ofstream summary(summary_file);
        
        summary << "sample_name\ttotal_organisms_detected\ttop_organism\ttop_abundance\n";
        
        // Read each sample's results and extract summary info
        for (const auto& [sample_name, r1_path, r2_path] : samples) {
            std::string abundance_file = output_dir + "/" + sample_name + "_abundance.tsv";
            
            std::ifstream results(abundance_file);
            if (!results.is_open()) continue;
            
            std::string line;
            std::getline(results, line);  // Skip header
            
            int organism_count = 0;
            std::string top_organism = "None";
            float top_abundance = 0.0f;
            
            while (std::getline(results, line)) {
                if (organism_count == 0) {
                    // First line is the most abundant
                    std::istringstream iss(line);
                    std::string field;
                    for (int i = 0; i < 4; i++) std::getline(iss, field, '\t');  // Skip to abundance
                    
                    if (std::getline(iss, field, '\t')) {
                        top_abundance = std::stof(field);
                        
                        // Get organism name (3rd field)
                        std::istringstream iss2(line);
                        for (int i = 0; i < 2; i++) std::getline(iss2, field, '\t');
                        std::getline(iss2, top_organism, '\t');
                    }
                }
                organism_count++;
            }
            
            summary << sample_name << "\t" << organism_count << "\t" 
                    << top_organism << "\t" << std::scientific << top_abundance << "\n";
        }
        
        summary.close();
        std::cout << "Batch summary written to: " << summary_file << std::endl;
    }
    
    void generate_batch_summary_csv(
        const BioGPU::SampleCSVParser& parser,
        const std::string& output_dir
    ) {
        std::cout << "Generating batch summary report..." << std::endl;
        
        std::string summary_file = output_dir + "/batch_summary.tsv";
        std::ofstream summary(summary_file);
        
        summary << "sample_name\ttotal_organisms_detected\ttop_organism\ttop_abundance\n";
        
        // Read each sample's results and extract summary info
        for (size_t i = 0; i < parser.getSampleCount(); i++) {
            const auto* sample = parser.getSample(i);
            if (!sample) continue;
            
            std::string abundance_file = output_dir + "/" + sample->sample_name + "/" + 
                                       sample->sample_name + "_abundance.tsv";
            
            std::ifstream results(abundance_file);
            if (!results.is_open()) continue;
            
            std::string line;
            std::getline(results, line);  // Skip header
            
            int organism_count = 0;
            std::string top_organism = "None";
            float top_abundance = 0.0f;
            
            while (std::getline(results, line)) {
                if (organism_count == 0) {
                    // First line is the most abundant
                    std::istringstream iss(line);
                    std::string field;
                    for (int i = 0; i < 4; i++) std::getline(iss, field, '\t');  // Skip to abundance
                    
                    if (std::getline(iss, field, '\t')) {
                        top_abundance = std::stof(field);
                        
                        // Get organism name (3rd field)
                        std::istringstream iss2(line);
                        for (int i = 0; i < 2; i++) std::getline(iss2, field, '\t');
                        std::getline(iss2, top_organism, '\t');
                    }
                }
                organism_count++;
            }
            
            summary << sample->sample_name << "\t" << organism_count << "\t" 
                    << top_organism << "\t" << std::scientific << top_abundance << "\n";
        }
        
        summary.close();
        std::cout << "Batch summary written to: " << summary_file << std::endl;
    }
};

// Enhanced main function for batch processing
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <database> <mode> [options]" << std::endl;
        std::cerr << "\nModes:" << std::endl;
        std::cerr << "  single <R1.fastq> [R2.fastq] [output_prefix]  - Process single sample" << std::endl;
        std::cerr << "  batch <sample_list> [output_dir]              - Process batch of samples" << std::endl;
        std::cerr << "\nBatch file formats:" << std::endl;
        std::cerr << "  CSV format (recommended - consistent with resistance pipeline):" << std::endl;
        std::cerr << "    SampleName,FilePath,R1 file,R2 file" << std::endl;
        std::cerr << "    Sample1,~/data/,sample1_R1.fq.gz,sample1_R2.fq.gz" << std::endl;
        std::cerr << "  Text format (tab-separated):" << std::endl;
        std::cerr << "    sample_name    R1_file    R2_file" << std::endl;
        std::cerr << "\nExamples:" << std::endl;
        std::cerr << "  # Single sample" << std::endl;
        std::cerr << "  " << argv[0] << " microbes.db single reads_R1.fq.gz reads_R2.fq.gz" << std::endl;
        std::cerr << "  # Batch processing (CSV)" << std::endl;
        std::cerr << "  " << argv[0] << " microbes.db batch samples.csv batch_results/" << std::endl;
        std::cerr << "  # Batch processing (text)" << std::endl;
        std::cerr << "  " << argv[0] << " microbes.db batch samples.txt batch_results/" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string mode = argv[2];
    
    try {
        BatchAdaptivePairedEndProfiler profiler;
        
        if (mode == "single") {
            // Single sample mode (backwards compatible)
            if (argc < 4) {
                std::cerr << "Single mode requires: <database> single <R1.fastq> [R2.fastq] [output_prefix]" << std::endl;
                return 1;
            }
            
            std::string r1_path = argv[3];
            std::string r2_path = "";
            std::string output_prefix = "single_profile";
            
            // Parse optional arguments
            if (argc > 4) {
                std::string arg4 = argv[4];
                if (arg4.find(".fastq") != std::string::npos || arg4.find(".fq") != std::string::npos) {
                    r2_path = arg4;
                    if (argc > 5) output_prefix = argv[5];
                } else {
                    output_prefix = arg4;
                }
            }
            
            // Load database once
            profiler.load_database_once(database_path);
            
            // Process single sample
            auto profiles = profiler.process_single_sample("single_sample", r1_path, r2_path);
            
            // Write results (reuse existing methods)
            profiler.print_results(profiles);
            profiler.write_results(profiles, output_prefix);
            
        } else if (mode == "batch") {
            // Batch processing mode
            if (argc < 4) {
                std::cerr << "Batch mode requires: <database> batch <sample_list> [output_dir]" << std::endl;
                std::cerr << "Sample list can be a .csv or .txt file" << std::endl;
                return 1;
            }
            
            std::string sample_list_file = argv[3];
            std::string output_dir = (argc > 4) ? argv[4] : "batch_results";
            
            // Check file extension to determine format
            bool is_csv = false;
            if (sample_list_file.size() >= 4) {
                std::string ext = sample_list_file.substr(sample_list_file.size() - 4);
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                is_csv = (ext == ".csv");
            }
            
            if (is_csv) {
                // Use CSV parser (consistent with resistance pipeline)
                profiler.process_sample_batch_csv(database_path, sample_list_file, output_dir);
            } else {
                // Use legacy text format parser
                auto samples = BatchAdaptivePairedEndProfiler::load_sample_list(sample_list_file);
                profiler.process_sample_batch(database_path, samples, output_dir);
            }
            
        } else {
            std::cerr << "Unknown mode: " << mode << std::endl;
            std::cerr << "Supported modes: single, batch" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}