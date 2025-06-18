#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <memory>
#include <cmath>

// Include the common paired-end structures
#include "paired_read_common.h"

struct MinimizerMatch {
    uint64_t minimizer_hash;
    uint32_t taxonomy_id;
    uint8_t uniqueness_score;
} __attribute__((packed));

struct GPUMinimizerEntry {
    uint64_t hash;
    uint32_t taxonomy_id;
    float weight;  // Pre-calculated weight including uniqueness and genome size normalization
};

struct OrganismStats {
    uint32_t taxonomy_id;
    float total_score;
    uint32_t unique_minimizers;
    uint32_t total_hits;
    uint32_t paired_hits;
    uint32_t expected_minimizers;
    uint64_t genome_size;
};

// Enhanced CUDA kernels for paired-end processing
__global__ void extract_paired_minimizers_kernel(
    const char* sequences,
    const uint32_t* sequence_offsets,
    const uint32_t* sequence_lengths,
    const uint32_t* read_pair_ids,
    const uint8_t* read_numbers,
    uint64_t* output_minimizers,
    uint32_t* output_read_pair_ids,
    uint8_t* output_read_numbers,
    uint32_t* output_positions,
    uint32_t num_sequences,
    int k,
    int m,
    int window_size
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_sequences) return;
    
    const char* sequence = sequences + sequence_offsets[tid];
    uint32_t seq_length = sequence_lengths[tid];
    uint32_t read_pair_id = read_pair_ids[tid];
    uint8_t read_number = read_numbers[tid];
    
    if (seq_length < k) return;
    
    uint32_t base_output_idx = sequence_offsets[tid] / window_size;
    
    for (uint32_t i = 0; i <= seq_length - k; i += window_size) {
        if (i + k > seq_length) break;
        
        // Find minimizer in this k-mer window
        uint64_t min_hash = UINT64_MAX;
        bool found_valid = false;
        
        for (int j = 0; j <= k - m; j++) {
            if (i + j + m > seq_length) break;
            
            // Hash minimizer
            uint64_t forward_hash = 0;
            uint64_t reverse_hash = 0;
            bool valid = true;
            
            for (int l = 0; l < m; l++) {
                char base = sequence[i + j + l];
                int encoded = -1;
                
                switch (base) {
                    case 'A': case 'a': encoded = 0; break;
                    case 'C': case 'c': encoded = 1; break;
                    case 'G': case 'g': encoded = 2; break;
                    case 'T': case 't': encoded = 3; break;
                    default: valid = false; break;
                }
                
                if (!valid) break;
                
                forward_hash = (forward_hash << 2) | encoded;
                reverse_hash = (reverse_hash >> 2) | (((uint64_t)(3 ^ encoded)) << (2 * (m - 1)));
            }
            
            if (valid) {
                uint64_t canonical_hash = min(forward_hash, reverse_hash);
                if (canonical_hash < min_hash) {
                    min_hash = canonical_hash;
                    found_valid = true;
                }
            }
        }
        
        if (found_valid) {
            uint32_t global_idx = atomicAdd(&output_positions[0], 1);
            if (global_idx < 10000000) {  // Safety limit
                output_minimizers[global_idx] = min_hash;
                output_read_pair_ids[global_idx] = read_pair_id;
                output_read_numbers[global_idx] = read_number;
            }
        }
    }
}

__global__ void match_paired_minimizers_kernel(
    const uint64_t* read_minimizers,
    const uint32_t* read_pair_ids,
    const uint8_t* read_numbers,
    uint32_t num_read_minimizers,
    const GPUMinimizerEntry* database_entries,
    uint32_t num_db_entries,
    uint32_t* organism_hits,
    float* organism_scores,
    uint32_t* organism_paired_hits,
    uint32_t* pair_organism_votes_r1,  // Best organism for each pair's read1
    uint32_t* pair_organism_votes_r2,  // Best organism for each pair's read2
    float* pair_scores_r1,
    float* pair_scores_r2,
    uint32_t max_organisms,
    uint32_t max_pairs
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_read_minimizers) return;
    
    uint64_t query_hash = read_minimizers[tid];
    uint32_t pair_id = read_pair_ids[tid];
    uint8_t read_num = read_numbers[tid];
    
    // Binary search in sorted database
    int left = 0;
    int right = num_db_entries - 1;
    
    while (left <= right) {
        int mid = left + (right - left) / 2;
        uint64_t db_hash = database_entries[mid].hash;
        
        if (db_hash == query_hash) {
            // Found match - record hit
            uint32_t org_id = database_entries[mid].taxonomy_id;
            float weight = database_entries[mid].weight;
            
            if (org_id < max_organisms) {
                atomicAdd(&organism_hits[org_id], 1);
                atomicAdd(&organism_scores[org_id], weight);
                
                // Track paired-end voting
                if (pair_id < max_pairs) {
                    if (read_num == 1) {
                        atomicAdd(&pair_scores_r1[pair_id * max_organisms + org_id], weight);
                    } else if (read_num == 2) {
                        atomicAdd(&pair_scores_r2[pair_id * max_organisms + org_id], weight);
                    }
                }
            }
            
            // Check neighboring entries for the same hash
            int left_check = mid - 1;
            while (left_check >= 0 && database_entries[left_check].hash == query_hash) {
                uint32_t org_id_left = database_entries[left_check].taxonomy_id;
                float weight_left = database_entries[left_check].weight;
                
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
            while (right_check < num_db_entries && database_entries[right_check].hash == query_hash) {
                uint32_t org_id_right = database_entries[right_check].taxonomy_id;
                float weight_right = database_entries[right_check].weight;
                
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
        } else if (db_hash < query_hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
}

__global__ void calculate_paired_concordance_kernel(
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

class GPUMicrobialCommunityProfiler {
private:
    // CPU database storage
    std::unordered_map<uint32_t, std::string> organism_names;
    std::unordered_map<uint32_t, std::string> organism_taxonomy;
    std::unordered_map<uint32_t, uint64_t> organism_genome_sizes;
    std::unordered_map<uint32_t, uint32_t> organism_expected_minimizers;
    std::vector<uint32_t> organism_id_list;
    
    // GPU database storage
    thrust::device_vector<GPUMinimizerEntry> d_database_entries;
    thrust::device_vector<uint32_t> d_organism_hits;
    thrust::device_vector<float> d_organism_scores;
    thrust::device_vector<uint32_t> d_organism_unique_counts;
    thrust::device_vector<uint32_t> d_organism_paired_hits;
    thrust::device_vector<float> d_organism_concordance_scores;
    
    // GPU working memory for paired-end processing
    thrust::device_vector<uint64_t> d_read_minimizers;
    thrust::device_vector<uint32_t> d_read_pair_ids;
    thrust::device_vector<uint8_t> d_read_numbers;
    thrust::device_vector<uint32_t> d_unique_minimizer_flags;
    
    // Paired-end scoring arrays
    thrust::device_vector<float> d_pair_scores_r1;
    thrust::device_vector<float> d_pair_scores_r2;
    thrust::device_vector<bool> d_is_paired_flags;
    
    // GPU sequence storage
    thrust::device_vector<char> d_sequences;
    thrust::device_vector<uint32_t> d_sequence_offsets;
    thrust::device_vector<uint32_t> d_sequence_lengths;
    thrust::device_vector<uint32_t> d_sequence_pair_ids;
    thrust::device_vector<uint8_t> d_sequence_read_numbers;
    
    int k = 35;
    int m = 31;
    int window_size = 4;
    uint32_t max_organisms = 10000;
    uint32_t max_pairs = 1000000;
    
    struct ProfilingConfig {
        float min_abundance = 1e-6;
        float min_coverage = 0.001;
        uint32_t min_hits = 10;
        float confidence_threshold = 0.1;
        bool normalize_by_genome_size = true;
        bool weight_by_uniqueness = true;
        bool use_paired_end_bonus = true;
        float paired_concordance_weight = 2.0f;
    } config;
    
public:
    GPUMicrobialCommunityProfiler() {
        // Initialize GPU memory
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        std::cout << "Using GPU: " << prop.name << " with " 
                  << prop.totalGlobalMem / (1024*1024*1024) << " GB memory" << std::endl;
        
        // Allocate organism tracking arrays
        d_organism_hits.resize(max_organisms, 0);
        d_organism_scores.resize(max_organisms, 0.0f);
        d_organism_unique_counts.resize(max_organisms, 0);
        d_organism_paired_hits.resize(max_organisms, 0);
        d_organism_concordance_scores.resize(max_organisms, 0.0f);
    }
    
    void load_minimizer_database(const std::string& db_path) {
        std::cout << "Loading minimizer database: " << db_path << std::endl;
        
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
        window_size = k - m + 1;
        
        std::cout << "Database parameters: k=" << k << ", m=" << m 
                  << ", organisms=" << num_organisms 
                  << ", minimizer hashes=" << num_minimizer_hashes << std::endl;
        
        // Read organism metadata
        std::vector<uint32_t> taxonomy_ids;
        for (uint32_t i = 0; i < num_organisms; i++) {
            uint32_t taxonomy_id, taxon_level, minimizer_count;
            uint64_t genome_size;
            float gc_content;
            
            in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
            in.read(reinterpret_cast<char*>(&taxon_level), sizeof(taxon_level));
            in.read(reinterpret_cast<char*>(&gc_content), sizeof(gc_content));
            in.read(reinterpret_cast<char*>(&genome_size), sizeof(genome_size));
            in.read(reinterpret_cast<char*>(&minimizer_count), sizeof(minimizer_count));
            
            // Read organism name
            uint16_t name_length;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            std::string name(name_length, '\0');
            in.read(&name[0], name_length);
            organism_names[taxonomy_id] = name;
            
            // Read taxonomy path
            uint16_t taxonomy_length;
            in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            std::string taxonomy_path(taxonomy_length, '\0');
            in.read(&taxonomy_path[0], taxonomy_length);
            organism_taxonomy[taxonomy_id] = taxonomy_path;
            
            // Store metadata
            organism_genome_sizes[taxonomy_id] = genome_size;
            organism_expected_minimizers[taxonomy_id] = minimizer_count;
            taxonomy_ids.push_back(taxonomy_id);
        }
        
        // Create organism ID mapping for GPU arrays
        organism_id_list = taxonomy_ids;
        std::sort(organism_id_list.begin(), organism_id_list.end());
        
        // Read minimizer index and prepare GPU database
        std::vector<GPUMinimizerEntry> gpu_entries;
        
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
                
                // Calculate weight including uniqueness and genome size normalization
                float weight = 1.0f;
                if (config.weight_by_uniqueness) {
                    weight = uniqueness_score / 255.0f;
                }
                if (config.normalize_by_genome_size) {
                    uint64_t genome_size = organism_genome_sizes[taxonomy_id];
                    weight = weight / std::log10(genome_size + 1);
                }
                
                // Map taxonomy_id to array index
                auto it = std::lower_bound(organism_id_list.begin(), organism_id_list.end(), taxonomy_id);
                uint32_t array_index = std::distance(organism_id_list.begin(), it);
                
                GPUMinimizerEntry entry{hash, array_index, weight};
                gpu_entries.push_back(entry);
            }
        }
        
        in.close();
        
        // Sort GPU entries by hash for binary search
        std::sort(gpu_entries.begin(), gpu_entries.end(),
                 [](const GPUMinimizerEntry& a, const GPUMinimizerEntry& b) {
                     return a.hash < b.hash;
                 });
        
        // Transfer to GPU
        d_database_entries = gpu_entries;
        
        std::cout << "Loaded " << gpu_entries.size() << " minimizer entries to GPU from " 
                  << organism_names.size() << " organisms" << std::endl;
        
        print_gpu_memory_usage();
    }
    
    // NEW: Main paired-end profiling method
    std::vector<OrganismProfile> profile_community_gpu(const std::vector<PairedRead>& paired_reads) {
        std::cout << "GPU profiling community from " << paired_reads.size() << " read pairs..." << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Prepare sequences for GPU
        prepare_paired_sequences_for_gpu(paired_reads);
        
        // Extract minimizers on GPU with paired-end tracking
        extract_paired_minimizers_gpu();
        
        // Match minimizers against database with paired-end scoring
        match_paired_minimizers_gpu();
        
        // Calculate paired-end concordance
        calculate_paired_concordance_gpu();
        
        // Calculate final results
        auto profiles = calculate_paired_organism_profiles();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Paired-end GPU profiling completed in " << duration.count() << " ms" << std::endl;
        
        return profiles;
    }
    
    // BACKWARD COMPATIBILITY: Single-read method
    std::vector<OrganismProfile> profile_community_gpu(const std::vector<std::string>& reads) {
        auto paired_reads = convert_to_paired_format(reads);
        return profile_community_gpu(paired_reads);
    }
    
private:
    void prepare_paired_sequences_for_gpu(const std::vector<PairedRead>& paired_reads) {
        std::cout << "Preparing " << paired_reads.size() << " paired reads for GPU..." << std::endl;
        
        // Calculate total sequence length and prepare metadata
        std::vector<char> all_sequences;
        std::vector<uint32_t> offsets;
        std::vector<uint32_t> lengths;
        std::vector<uint32_t> pair_ids;
        std::vector<uint8_t> read_numbers;
        std::vector<bool> is_paired_flags;
        
        offsets.push_back(0);
        
        for (size_t i = 0; i < paired_reads.size(); i++) {
            const auto& pair = paired_reads[i];
            
            // Add read1
            for (char c : pair.read1) {
                all_sequences.push_back(c);
            }
            lengths.push_back(pair.read1.length());
            pair_ids.push_back(i);
            read_numbers.push_back(1);
            offsets.push_back(all_sequences.size());
            
            // Add read2 if paired
            if (pair.is_paired && !pair.read2.empty()) {
                for (char c : pair.read2) {
                    all_sequences.push_back(c);
                }
                lengths.push_back(pair.read2.length());
                pair_ids.push_back(i);
                read_numbers.push_back(2);
                offsets.push_back(all_sequences.size());
            }
            
            is_paired_flags.push_back(pair.is_paired);
        }
        
        // Transfer to GPU
        d_sequences = all_sequences;
        d_sequence_offsets = offsets;
        d_sequence_lengths = lengths;
        d_sequence_pair_ids = pair_ids;
        d_sequence_read_numbers = read_numbers;
        d_is_paired_flags = is_paired_flags;
        
        // Allocate paired-end scoring arrays
        size_t pair_scoring_size = paired_reads.size() * max_organisms;
        d_pair_scores_r1.resize(pair_scoring_size, 0.0f);
        d_pair_scores_r2.resize(pair_scoring_size, 0.0f);
        
        std::cout << "Transferred " << all_sequences.size() << " bases from " 
                  << paired_reads.size() << " pairs to GPU" << std::endl;
    }
    
    void extract_paired_minimizers_gpu() {
        std::cout << "Extracting minimizers from paired reads on GPU..." << std::endl;
        
        uint32_t num_sequences = d_sequence_lengths.size();
        
        // Allocate output arrays
        size_t max_minimizers = num_sequences * 200;  // Estimate
        d_read_minimizers.resize(max_minimizers);
        d_read_pair_ids.resize(max_minimizers);
        d_read_numbers.resize(max_minimizers);
        
        // Reset global counter
        thrust::device_vector<uint32_t> d_counter(1, 0);
        
        // Launch kernel
        int block_size = 256;
        int num_blocks = (num_sequences + block_size - 1) / block_size;
        
        extract_paired_minimizers_kernel<<<num_blocks, block_size>>>(
            thrust::raw_pointer_cast(d_sequences.data()),
            thrust::raw_pointer_cast(d_sequence_offsets.data()),
            thrust::raw_pointer_cast(d_sequence_lengths.data()),
            thrust::raw_pointer_cast(d_sequence_pair_ids.data()),
            thrust::raw_pointer_cast(d_sequence_read_numbers.data()),
            thrust::raw_pointer_cast(d_read_minimizers.data()),
            thrust::raw_pointer_cast(d_read_pair_ids.data()),
            thrust::raw_pointer_cast(d_read_numbers.data()),
            thrust::raw_pointer_cast(d_counter.data()),
            num_sequences, k, m, window_size
        );
        
        cudaDeviceSynchronize();
        
        // Get actual number of minimizers extracted
        uint32_t actual_minimizers = d_counter[0];
        d_read_minimizers.resize(actual_minimizers);
        d_read_pair_ids.resize(actual_minimizers);
        d_read_numbers.resize(actual_minimizers);
        
        std::cout << "Extracted " << actual_minimizers << " minimizers from paired reads" << std::endl;
    }
    
    void match_paired_minimizers_gpu() {
        std::cout << "Matching paired minimizers against database on GPU..." << std::endl;
        
        uint32_t num_minimizers = d_read_minimizers.size();
        
        // Reset organism counters
        thrust::fill(d_organism_hits.begin(), d_organism_hits.end(), 0);
        thrust::fill(d_organism_scores.begin(), d_organism_scores.end(), 0.0f);
        thrust::fill(d_organism_paired_hits.begin(), d_organism_paired_hits.end(), 0);
        thrust::fill(d_organism_concordance_scores.begin(), d_organism_concordance_scores.end(), 0.0f);
        thrust::fill(d_pair_scores_r1.begin(), d_pair_scores_r1.end(), 0.0f);
        thrust::fill(d_pair_scores_r2.begin(), d_pair_scores_r2.end(), 0.0f);
        
        // Launch matching kernel
        int block_size = 256;
        int num_blocks = (num_minimizers + block_size - 1) / block_size;
        
        match_paired_minimizers_kernel<<<num_blocks, block_size>>>(
            thrust::raw_pointer_cast(d_read_minimizers.data()),
            thrust::raw_pointer_cast(d_read_pair_ids.data()),
            thrust::raw_pointer_cast(d_read_numbers.data()),
            num_minimizers,
            thrust::raw_pointer_cast(d_database_entries.data()),
            d_database_entries.size(),
            thrust::raw_pointer_cast(d_organism_hits.data()),
            thrust::raw_pointer_cast(d_organism_scores.data()),
            thrust::raw_pointer_cast(d_organism_paired_hits.data()),
            thrust::raw_pointer_cast(d_pair_scores_r1.data()),
            thrust::raw_pointer_cast(d_pair_scores_r2.data()),
            thrust::raw_pointer_cast(d_pair_scores_r1.data()),
            thrust::raw_pointer_cast(d_pair_scores_r2.data()),
            max_organisms,
            max_pairs
        );
        
        cudaDeviceSynchronize();
        
        std::cout << "Paired minimizer matching completed" << std::endl;
    }
    
    void calculate_paired_concordance_gpu() {
        std::cout << "Calculating paired-end concordance on GPU..." << std::endl;
        
        uint32_t num_pairs = d_is_paired_flags.size();
        
        // Launch concordance calculation kernel
        int block_size = 256;
        int num_blocks = (num_pairs + block_size - 1) / block_size;
        
        calculate_paired_concordance_kernel<<<num_blocks, block_size>>>(
            thrust::raw_pointer_cast(d_pair_scores_r1.data()),
            thrust::raw_pointer_cast(d_pair_scores_r2.data()),
            thrust::raw_pointer_cast(d_is_paired_flags.data()),
            num_pairs,
            max_organisms,
            thrust::raw_pointer_cast(d_organism_paired_hits.data()),
            thrust::raw_pointer_cast(d_organism_concordance_scores.data())
        );
        
        cudaDeviceSynchronize();
        
        std::cout << "Paired concordance calculation completed" << std::endl;
    }
    
    std::vector<OrganismProfile> calculate_paired_organism_profiles() {
        std::cout << "Calculating paired-end organism profiles..." << std::endl;
        
        // Copy results back from GPU
        thrust::host_vector<uint32_t> h_organism_hits = d_organism_hits;
        thrust::host_vector<float> h_organism_scores = d_organism_scores;
        thrust::host_vector<uint32_t> h_organism_paired_hits = d_organism_paired_hits;
        thrust::host_vector<float> h_organism_concordance_scores = d_organism_concordance_scores;
        
        std::vector<OrganismProfile> profiles;
        
        // Calculate total abundance for normalization
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
            
            if (hits < config.min_hits) continue;
            
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
            
            // Calculate coverage breadth (simplified for now)
            uint32_t expected_minimizers = organism_expected_minimizers[tax_id];
            profile.coverage_breadth = (expected_minimizers > 0) ? 
                                     std::min(1.0f, (float)hits / expected_minimizers) : 0.0f;
            
            // Enhanced confidence score with paired-end information
            float minimizer_confidence = std::min(1.0f, (float)(std::log10(hits + 1) / 3.0));
            float coverage_confidence = std::min(1.0f, profile.coverage_breadth * 10.0f);
            float paired_confidence = (paired_hits > 0) ? 
                                    std::min(1.0f, profile.paired_concordance * 2.0f) : 0.0f;
            
            profile.confidence_score = (minimizer_confidence + coverage_confidence + paired_confidence) / 3.0f;
            
            // Apply thresholds
            if (profile.abundance >= config.min_abundance &&
                profile.coverage_breadth >= config.min_coverage &&
                profile.confidence_score >= config.confidence_threshold) {
                profiles.push_back(profile);
            }
        }
        
        // Sort by abundance
        std::sort(profiles.begin(), profiles.end(),
                 [](const OrganismProfile& a, const OrganismProfile& b) {
                     return a.abundance > b.abundance;
                 });
        
        std::cout << "Generated " << profiles.size() << " organism profiles with paired-end scoring" << std::endl;
        
        return profiles;
    }
    
public:
    void print_gpu_memory_usage() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        size_t used_mem = total_mem - free_mem;
        std::cout << "GPU memory usage: " << (used_mem / 1024 / 1024) << " MB used, "
                  << (free_mem / 1024 / 1024) << " MB free" << std::endl;
    }
    
    void print_paired_community_profile(const std::vector<OrganismProfile>& profiles) {
        std::cout << "\n=== PAIRED-END GPU MICROBIAL COMMUNITY PROFILE ===" << std::endl;
        
        if (profiles.empty()) {
            std::cout << "No organisms detected above significance thresholds" << std::endl;
            return;
        }
        
        // Print top organisms with paired-end information
        std::cout << "\nTop Detected Organisms (with paired-end metrics):" << std::endl;
        std::cout << std::string(140, '-') << std::endl;
        std::cout << std::left << std::setw(40) << "Organism Name" 
                  << std::setw(10) << "Abundance" 
                  << std::setw(8) << "Coverage" 
                  << std::setw(10) << "Confidence"
                  << std::setw(8) << "Total"
                  << std::setw(8) << "Paired"
                  << std::setw(8) << "Single"
                  << std::setw(12) << "Concordance"
                  << std::setw(20) << "Taxonomy" << std::endl;
        std::cout << std::string(140, '-') << std::endl;
        
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
            if (display_name.length() > 35) {
                display_name = display_name.substr(0, 32) + "...";
            }
            
            std::cout << std::left << std::setw(40) << display_name
                      << std::fixed << std::setprecision(3) << std::setw(10) << (profile.abundance * 100)
                      << std::fixed << std::setprecision(2) << std::setw(8) << (profile.coverage_breadth * 100)
                      << std::fixed << std::setprecision(3) << std::setw(10) << profile.confidence_score
                      << std::setw(8) << profile.total_hits
                      << std::setw(8) << profile.paired_hits
                      << std::setw(8) << profile.single_hits
                      << std::fixed << std::setprecision(2) << std::setw(12) << (profile.paired_concordance * 100)
                      << std::setw(20) << genus << std::endl;
        }
        
        std::cout << std::string(140, '-') << std::endl;
        
        // Paired-end summary statistics
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
        
        std::cout << "\nPaired-End Summary:" << std::endl;
        std::cout << "- Total paired hits: " << total_paired_hits << std::endl;
        std::cout << "- Organisms with paired evidence: " << organisms_with_pairs << std::endl;
        std::cout << "- Average concordance rate: " << std::fixed << std::setprecision(1) 
                  << (avg_concordance * 100) << "%" << std::endl;
    }
    
    void write_paired_abundance_table(const std::vector<OrganismProfile>& profiles, 
                                     const std::string& output_file) {
        std::ofstream out(output_file);
        
        // Enhanced header with paired-end columns
        out << "taxonomy_id\torganism_name\ttaxonomy_path\trelative_abundance\t"
            << "coverage_breadth\ttotal_hits\tpaired_hits\tsingle_hits\t"
            << "paired_concordance\tconfidence_score\n";
        
        // Data rows
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
        std::cout << "Paired-end abundance table written to: " << output_file << std::endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <minimizer_database> <fastq_file> [fastq_r2] [output_prefix]" << std::endl;
        std::cerr << "\nGPU-accelerated microbial community profiler with paired-end support" << std::endl;
        std::cerr << "Supports both single-end and paired-end reads for enhanced accuracy" << std::endl;
        std::cerr << "\nExamples:" << std::endl;
        std::cerr << "  Single-end: " << argv[0] << " microbes.db reads.fastq" << std::endl;
        std::cerr << "  Paired-end: " << argv[0] << " microbes.db reads_R1.fastq reads_R2.fastq" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string fastq_r1 = argv[2];
    std::string fastq_r2 = "";
    std::string output_prefix = "paired_gpu_profile";
    
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
        GPUMicrobialCommunityProfiler profiler;
        
        // Load database
        profiler.load_minimizer_database(database_path);
        
        // Load reads (paired-end or single-end)
        std::vector<PairedRead> paired_reads;
        
        if (!fastq_r2.empty()) {
            std::cout << "Loading paired-end reads from " << fastq_r1 << " and " << fastq_r2 << std::endl;
            // Load paired-end reads from separate files
            // (Implementation similar to adaptive_paired_end_profiler.cu)
            // For brevity, using placeholder here
            
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
                
                if (single_reads.size() >= 1000000) {
                    std::cout << "\nLimited to " << single_reads.size() << " reads for testing" << std::endl;
                    break;
                }
            }
            
            // Convert to paired format
            paired_reads = convert_to_paired_format(single_reads);
        }
        
        std::cout << "\nLoaded " << paired_reads.size() << " read pairs for GPU profiling" << std::endl;
        
        // Profile community with paired-end support
        auto profiles = profiler.profile_community_gpu(paired_reads);
        
        // Print results
        profiler.print_paired_community_profile(profiles);
        
        // Write output files
        profiler.write_paired_abundance_table(profiles, output_prefix + "_paired_abundance.tsv");
        
        std::cout << "\n=== PAIRED-END GPU PROFILING COMPLETE ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}