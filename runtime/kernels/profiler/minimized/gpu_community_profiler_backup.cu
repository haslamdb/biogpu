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

// GPU-accelerated microbial community profiler
// Uses CUDA kernels for high-throughput minimizer matching and scoring

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
    uint32_t expected_minimizers;
    uint64_t genome_size;
};

// CUDA kernels for GPU processing
__global__ void extract_minimizers_kernel(
    const char* sequences,
    const uint32_t* sequence_offsets,
    const uint32_t* sequence_lengths,
    uint64_t* output_minimizers,
    uint32_t* output_read_ids,
    uint32_t* output_positions,
    uint32_t num_reads,
    int k,
    int m,
    int window_size
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const char* sequence = sequences + sequence_offsets[tid];
    uint32_t seq_length = sequence_lengths[tid];
    
    if (seq_length < k) return;
    
    uint32_t output_idx = sequence_offsets[tid] / window_size;  // Approximate output position
    
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
            uint32_t global_idx = atomicAdd(&output_positions[0], 1);  // Global counter
            if (global_idx < 10000000) {  // Safety limit
                output_minimizers[global_idx] = min_hash;
                output_read_ids[global_idx] = tid;
            }
        }
    }
}

__global__ void match_minimizers_kernel(
    const uint64_t* read_minimizers,
    const uint32_t* read_ids,
    uint32_t num_read_minimizers,
    const GPUMinimizerEntry* database_entries,
    uint32_t num_db_entries,
    uint32_t* organism_hits,
    float* organism_scores,
    uint32_t* unique_minimizer_flags,
    uint32_t max_organisms
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_read_minimizers) return;
    
    uint64_t query_hash = read_minimizers[tid];
    
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
                
                // Mark this minimizer as seen for this organism
                unique_minimizer_flags[tid * max_organisms + org_id] = 1;
            }
            
            // Check neighboring entries for the same hash
            int left_check = mid - 1;
            while (left_check >= 0 && database_entries[left_check].hash == query_hash) {
                uint32_t org_id_left = database_entries[left_check].taxonomy_id;
                float weight_left = database_entries[left_check].weight;
                
                if (org_id_left < max_organisms) {
                    atomicAdd(&organism_hits[org_id_left], 1);
                    atomicAdd(&organism_scores[org_id_left], weight_left);
                    unique_minimizer_flags[tid * max_organisms + org_id_left] = 1;
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
                    unique_minimizer_flags[tid * max_organisms + org_id_right] = 1;
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

__global__ void calculate_unique_minimizers_kernel(
    const uint32_t* unique_minimizer_flags,
    uint32_t num_minimizers,
    uint32_t max_organisms,
    uint32_t* organism_unique_counts
) {
    int org_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (org_id >= max_organisms) return;
    
    uint32_t unique_count = 0;
    for (uint32_t i = 0; i < num_minimizers; i++) {
        if (unique_minimizer_flags[i * max_organisms + org_id] == 1) {
            unique_count++;
        }
    }
    
    organism_unique_counts[org_id] = unique_count;
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
    
    // GPU working memory
    thrust::device_vector<uint64_t> d_read_minimizers;
    thrust::device_vector<uint32_t> d_read_ids;
    thrust::device_vector<uint32_t> d_unique_minimizer_flags;
    
    // GPU sequence storage
    thrust::device_vector<char> d_sequences;
    thrust::device_vector<uint32_t> d_sequence_offsets;
    thrust::device_vector<uint32_t> d_sequence_lengths;
    
    int k = 35;
    int m = 31;
    int window_size = 4;
    uint32_t max_organisms = 10000;
    
    struct ProfilingConfig {
        float min_abundance = 1e-6;
        float min_coverage = 0.001;
        uint32_t min_hits = 10;
        float confidence_threshold = 0.1;
        bool normalize_by_genome_size = true;
        bool weight_by_uniqueness = true;
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
    
    void print_gpu_memory_usage() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        size_t used_mem = total_mem - free_mem;
        std::cout << "GPU memory usage: " << (used_mem / 1024 / 1024) << " MB used, "
                  << (free_mem / 1024 / 1024) << " MB free" << std::endl;
    }
    
    std::vector<OrganismProfile> profile_community_gpu(const std::vector<std::string>& reads) {
        std::cout << "GPU profiling community composition from " << reads.size() << " reads..." << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Prepare sequences for GPU
        prepare_sequences_for_gpu(reads);
        
        // Extract minimizers on GPU
        extract_minimizers_gpu();
        
        // Match minimizers against database on GPU
        match_minimizers_gpu();
        
        // Calculate final results
        auto profiles = calculate_organism_profiles();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "GPU profiling completed in " << duration.count() << " ms" << std::endl;
        
        return profiles;
    }
    
private:
    void prepare_sequences_for_gpu(const std::vector<std::string>& reads) {
        std::cout << "Preparing " << reads.size() << " sequences for GPU..." << std::endl;
        
        // Calculate total sequence length and prepare offsets
        std::vector<char> all_sequences;
        std::vector<uint32_t> offsets;
        std::vector<uint32_t> lengths;
        
        offsets.push_back(0);
        for (const auto& read : reads) {
            for (char c : read) {
                all_sequences.push_back(c);
            }
            lengths.push_back(read.length());
            offsets.push_back(all_sequences.size());
        }
        
        // Transfer to GPU
        d_sequences = all_sequences;
        d_sequence_offsets = offsets;
        d_sequence_lengths = lengths;
        
        std::cout << "Transferred " << all_sequences.size() << " bases to GPU" << std::endl;
    }
    
    void extract_minimizers_gpu() {
        std::cout << "Extracting minimizers on GPU..." << std::endl;
        
        uint32_t num_reads = d_sequence_lengths.size();
        
        // Allocate output arrays
        size_t max_minimizers = num_reads * 200;  // Estimate
        d_read_minimizers.resize(max_minimizers);
        d_read_ids.resize(max_minimizers);
        
        // Reset global counter
        thrust::device_vector<uint32_t> d_counter(1, 0);
        
        // Launch kernel
        int block_size = 256;
        int num_blocks = (num_reads + block_size - 1) / block_size;
        
        extract_minimizers_kernel<<<num_blocks, block_size>>>(
            thrust::raw_pointer_cast(d_sequences.data()),
            thrust::raw_pointer_cast(d_sequence_offsets.data()),
            thrust::raw_pointer_cast(d_sequence_lengths.data()),
            thrust::raw_pointer_cast(d_read_minimizers.data()),
            thrust::raw_pointer_cast(d_read_ids.data()),
            thrust::raw_pointer_cast(d_counter.data()),
            num_reads, k, m, window_size
        );
        
        cudaDeviceSynchronize();
        
        // Get actual number of minimizers extracted
        uint32_t actual_minimizers = d_counter[0];
        d_read_minimizers.resize(actual_minimizers);
        d_read_ids.resize(actual_minimizers);
        
        std::cout << "Extracted " << actual_minimizers << " minimizers" << std::endl;
    }
    
    void match_minimizers_gpu() {
        std::cout << "Matching minimizers against database on GPU..." << std::endl;
        
        uint32_t num_minimizers = d_read_minimizers.size();
        
        // Reset organism counters
        thrust::fill(d_organism_hits.begin(), d_organism_hits.end(), 0);
        thrust::fill(d_organism_scores.begin(), d_organism_scores.end(), 0.0f);
        thrust::fill(d_organism_unique_counts.begin(), d_organism_unique_counts.end(), 0);
        
        // Allocate unique minimizer tracking
        d_unique_minimizer_flags.resize(num_minimizers * max_organisms, 0);
        
        // Launch matching kernel
        int block_size = 256;
        int num_blocks = (num_minimizers + block_size - 1) / block_size;
        
        match_minimizers_kernel<<<num_blocks, block_size>>>(
            thrust::raw_pointer_cast(d_read_minimizers.data()),
            thrust::raw_pointer_cast(d_read_ids.data()),
            num_minimizers,
            thrust::raw_pointer_cast(d_database_entries.data()),
            d_database_entries.size(),
            thrust::raw_pointer_cast(d_organism_hits.data()),
            thrust::raw_pointer_cast(d_organism_scores.data()),
            thrust::raw_pointer_cast(d_unique_minimizer_flags.data()),
            max_organisms
        );
        
        cudaDeviceSynchronize();
        
        // Calculate unique minimizer counts
        num_blocks = (max_organisms + block_size - 1) / block_size;
        
        calculate_unique_minimizers_kernel<<<num_blocks, block_size>>>(
            thrust::raw_pointer_cast(d_unique_minimizer_flags.data()),
            num_minimizers,
            max_organisms,
            thrust::raw_pointer_cast(d_organism_unique_counts.data())
        );
        
        cudaDeviceSynchronize();
        
        std::cout << "Minimizer matching completed" << std::endl;
    }
    
    std::vector<OrganismProfile> calculate_organism_profiles() {
        std::cout << "Calculating organism profiles..." << std::endl;
        
        // Copy results back from GPU
        thrust::host_vector<uint32_t> h_organism_hits = d_organism_hits;
        thrust::host_vector<float> h_organism_scores = d_organism_scores;
        thrust::host_vector<uint32_t> h_organism_unique_counts = d_organism_unique_counts;
        
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
            uint32_t unique_minimizers = h_organism_unique_counts[i];
            
            if (hits < config.min_hits) continue;
            
            OrganismProfile profile;
            profile.taxonomy_id = tax_id;
            profile.name = organism_names[tax_id];
            profile.taxonomy_path = organism_taxonomy[tax_id];
            profile.total_hits = hits;
            profile.unique_minimizers = unique_minimizers;
            
            // Calculate normalized abundance
            profile.abundance = (total_abundance > 0) ? score / total_abundance : 0.0f;
            
            // Calculate coverage breadth
            uint32_t expected_minimizers = organism_expected_minimizers[tax_id];
            profile.coverage_breadth = (expected_minimizers > 0) ? 
                                     (float)unique_minimizers / expected_minimizers : 0.0f;
            
            // Calculate confidence score
            float minimizer_confidence = std::min(1.0f, (float)(std::log10(unique_minimizers + 1) / 3.0));
            float coverage_confidence = std::min(1.0f, profile.coverage_breadth * 10.0f);
            profile.confidence_score = (minimizer_confidence + coverage_confidence) / 2.0f;
            
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
        
        std::cout << "Generated " << profiles.size() << " organism profiles above thresholds" << std::endl;
        
        return profiles;
    }
    
public:
    void print_community_profile(const std::vector<OrganismProfile>& profiles) {
        std::cout << "\n=== GPU MICROBIAL COMMUNITY PROFILE ===" << std::endl;
        
        if (profiles.empty()) {
            std::cout << "No organisms detected above significance thresholds" << std::endl;
            return;
        }
        
        // Calculate diversity metrics
        float shannon_diversity = 0.0f;
        float simpson_diversity = 0.0f;
        
        for (const auto& profile : profiles) {
            if (profile.abundance > 0) {
                shannon_diversity -= profile.abundance * std::log(profile.abundance);
                simpson_diversity += profile.abundance * profile.abundance;
            }
        }
        simpson_diversity = 1.0f - simpson_diversity;
        
        std::cout << "\nCommunity Diversity:" << std::endl;
        std::cout << "- Species richness: " << profiles.size() << std::endl;
        std::cout << "- Shannon diversity: " << std::fixed << std::setprecision(3) << shannon_diversity << std::endl;
        std::cout << "- Simpson diversity: " << std::fixed << std::setprecision(3) << simpson_diversity << std::endl;
        
        // Print top organisms
        std::cout << "\nTop Detected Organisms:" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        std::cout << std::left << std::setw(50) << "Organism Name" 
                  << std::setw(12) << "Abundance" 
                  << std::setw(12) << "Coverage" 
                  << std::setw(12) << "Confidence"
                  << std::setw(10) << "Hits"
                  << std::setw(22) << "Taxonomy" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        
        for (size_t i = 0; i < std::min(size_t(20), profiles.size()); i++) {
            const auto& profile = profiles[i];
            
            // Extract genus from taxonomy path
            std::string genus = "Unknown";
            size_t last_semicolon = profile.taxonomy_path.find_last_of(';');
            if (last_semicolon != std::string::npos) {
                genus = profile.taxonomy_path.substr(last_semicolon + 1);
            }
            
            // Truncate long names
            std::string display_name = profile.name;
            if (display_name.length() > 45) {
                display_name = display_name.substr(0, 42) + "...";
            }
            
            std::cout << std::left << std::setw(50) << display_name
                      << std::fixed << std::setprecision(4) << std::setw(12) << (profile.abundance * 100)
                      << std::fixed << std::setprecision(3) << std::setw(12) << (profile.coverage_breadth * 100)
                      << std::fixed << std::setprecision(3) << std::setw(12) << profile.confidence_score
                      << std::setw(10) << profile.total_hits
                      << std::setw(22) << genus << std::endl;
        }
        
        std::cout << std::string(120, '-') << std::endl;
        
        // Performance summary
        std::cout << "\nGPU Performance Summary:" << std::endl;
        size_t total_minimizers = 0;
        size_t total_hits = 0;
        for (const auto& profile : profiles) {
            total_hits += profile.total_hits;
        }
        
        std::cout << "- Total minimizer hits: " << total_hits << std::endl;
        std::cout << "- Organisms above thresholds: " << profiles.size() << std::endl;
    }
    
    void write_abundance_table(const std::vector<OrganismProfile>& profiles, 
                              const std::string& output_file) {
        std::ofstream out(output_file);
        
        // Header
        out << "taxonomy_id\torganism_name\ttaxonomy_path\trelative_abundance\t"
            << "coverage_breadth\tunique_minimizers\ttotal_hits\tconfidence_score\n";
        
        // Data rows
        for (const auto& profile : profiles) {
            out << profile.taxonomy_id << "\t"
                << profile.name << "\t"
                << profile.taxonomy_path << "\t"
                << std::scientific << profile.abundance << "\t"
                << std::fixed << std::setprecision(6) << profile.coverage_breadth << "\t"
                << profile.unique_minimizers << "\t"
                << profile.total_hits << "\t"
                << std::fixed << std::setprecision(4) << profile.confidence_score << "\n";
        }
        
        out.close();
        std::cout << "Abundance table written to: " << output_file << std::endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <minimizer_database> <fastq_file> [output_prefix]" << std::endl;
        std::cerr << "\nGPU-accelerated microbial community profiler" << std::endl;
        std::cerr << "High-performance CUDA-based minimizer matching" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string fastq_file = argv[2];
    std::string output_prefix = (argc > 3) ? argv[3] : "gpu_community_profile";
    
    try {
        GPUMicrobialCommunityProfiler profiler;
        
        // Load database
        profiler.load_minimizer_database(database_path);
        
        // Read FASTQ file
        std::vector<std::string> reads;
        std::ifstream file(fastq_file);
        std::string line;
        int line_count = 0;
        
        std::cout << "Reading FASTQ file: " << fastq_file << std::endl;
        
        while (std::getline(file, line)) {
            line_count++;
            if (line_count % 4 == 2) {  // Sequence line
                if (!line.empty() && line.length() >= 35) {  // Minimum length check
                    reads.push_back(line);
                }
            }
            
            // Progress indicator
            if (line_count % 400000 == 0) {
                std::cout << "\rRead " << (line_count / 4) << " sequences..." << std::flush;
            }
            
            // Process all reads for production use
            // Remove this limit for full-scale analysis
            if (reads.size() >= 1000000) {  // 1M reads for testing
                std::cout << "\nLimited to " << reads.size() << " reads for testing" << std::endl;
                break;
            }
        }
        
        std::cout << "\nLoaded " << reads.size() << " valid reads for GPU profiling" << std::endl;
        
        // Profile community on GPU
        auto profiles = profiler.profile_community_gpu(reads);
        
        // Print results
        profiler.print_community_profile(profiles);
        
        // Write output files
        profiler.write_abundance_table(profiles, output_prefix + "_abundance_table.tsv");
        
        std::cout << "\n=== GPU PROFILING COMPLETE ===" << std::endl;
        std::cout << "Results written to: " << output_prefix << "_abundance_table.tsv" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}