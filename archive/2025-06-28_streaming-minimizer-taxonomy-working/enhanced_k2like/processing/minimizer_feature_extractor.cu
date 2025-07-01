// processing/minimizer_feature_extractor.cu
// Two-pass minimizer feature extraction for enhanced ML-ready metadata
// Calculates uniqueness, conservation, and contamination features

#include "minimizer_feature_extractor.h"
#include "../gpu_kraken_types.h"
#include "../gpu/gpu_database_kernels.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/transform.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>

// ===========================
// Device Functions
// ===========================

__device__ float calculate_shannon_entropy(uint32_t* counts, int num_taxa) {
    if (num_taxa <= 1) return 0.0f;
    
    uint32_t total = 0;
    for (int i = 0; i < num_taxa; i++) {
        total += counts[i];
    }
    
    if (total == 0) return 0.0f;
    
    float entropy = 0.0f;
    for (int i = 0; i < num_taxa; i++) {
        if (counts[i] > 0) {
            float p = (float)counts[i] / total;
            entropy -= p * log2f(p);
        }
    }
    
    return entropy;
}

__device__ uint16_t encode_ml_weight(float weight) {
    // Encode float weight (0.0-1.0) to uint16_t (0-65535)
    weight = fmaxf(0.0f, fminf(1.0f, weight)); // Clamp to [0,1]
    return (uint16_t)(weight * 65535.0f);
}

__device__ uint8_t calculate_conservation_category(float conservation_score) {
    // Map conservation score to 3-bit category (0-7)
    if (conservation_score < 0.1f) return 0;      // Highly specific
    else if (conservation_score < 0.2f) return 1; // Species-specific
    else if (conservation_score < 0.3f) return 2; // Genus-specific
    else if (conservation_score < 0.5f) return 3; // Family-specific
    else if (conservation_score < 0.7f) return 4; // Order-specific
    else if (conservation_score < 0.85f) return 5; // Class-specific
    else if (conservation_score < 0.95f) return 6; // Phylum-specific
    else return 7;                                // Universal
}

// ===========================
// CUDA Kernels
// ===========================

__global__ void compute_minimizer_features_kernel(
    GPUMinimizerHit* minimizer_hits,
    const MinimizerFeatureStats* statistics,
    const uint64_t* unique_minimizers,
    int num_hits,
    int num_unique_minimizers,
    float total_genomes) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    GPUMinimizerHit& hit = minimizer_hits[tid];
    uint64_t minimizer_hash = hit.minimizer_hash;
    
    // Binary search for minimizer in unique list
    int left = 0, right = num_unique_minimizers - 1;
    int stats_idx = -1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        if (unique_minimizers[mid] == minimizer_hash) {
            stats_idx = mid;
            break;
        } else if (unique_minimizers[mid] < minimizer_hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    if (stats_idx >= 0) {
        const MinimizerFeatureStats& stats = statistics[stats_idx];
        
        // Calculate uniqueness score (inverse frequency)
        float frequency = (float)stats.total_occurrences / total_genomes;
        float uniqueness = 1.0f - fminf(frequency, 1.0f);
        hit.ml_weight = encode_ml_weight(uniqueness);
        
        // Calculate conservation score based on taxonomic distribution
        float conservation = (float)stats.unique_taxa / stats.total_occurrences;
        uint8_t conservation_cat = calculate_conservation_category(conservation);
        
        // Update feature flags
        uint16_t new_flags = hit.feature_flags;
        
        // Clear conservation bits (8-10) and set new value
        new_flags &= ~(0x7 << 8);
        new_flags |= (conservation_cat << 8);
        
        // Set uniqueness indicator (bit 11)
        if (uniqueness > 0.9f) {
            new_flags |= (1 << 11);
        }
        
        // Set rare minimizer flag (bit 12)
        if (stats.total_occurrences < 5) {
            new_flags |= (1 << 12);
        }
        
        hit.feature_flags = new_flags;
    }
}

__global__ void mark_position_bias_kernel(
    GPUMinimizerHit* minimizer_hits,
    int num_hits,
    int clustering_window) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    // Check for position clustering
    int nearby_count = 0;
    uint32_t current_pos = minimizer_hits[tid].position;
    uint16_t current_genome = minimizer_hits[tid].genome_id;
    
    // Look at nearby minimizers (simplified - in production, use sorted position index)
    for (int i = max(0, tid - 50); i < min(num_hits, tid + 50); i++) {
        if (i == tid) continue;
        
        if (minimizer_hits[i].genome_id == current_genome) {
            uint32_t pos_diff = abs((int)minimizer_hits[i].position - (int)current_pos);
            if (pos_diff < clustering_window) {
                nearby_count++;
            }
        }
    }
    
    // Set position bias bit if clustered
    if (nearby_count > 3) {
        minimizer_hits[tid].feature_flags |= (1 << 6);
    }
}

// ===========================
// MinimizerFeatureExtractor Implementation
// ===========================

MinimizerFeatureExtractor::MinimizerFeatureExtractor(size_t max_minimizers, size_t max_genomes)
    : max_minimizers_(max_minimizers), max_genomes_(max_genomes), 
      d_statistics_(nullptr), d_unique_minimizers_(nullptr),
      num_unique_minimizers_(0), total_genomes_processed_(0) {
    
    // Pre-allocate statistics storage
    cudaMalloc(&d_statistics_, max_minimizers * sizeof(MinimizerFeatureStats));
    cudaMalloc(&d_unique_minimizers_, max_minimizers * sizeof(uint64_t));
    
    // Initialize contamination patterns
    initialize_contamination_patterns();
}

MinimizerFeatureExtractor::~MinimizerFeatureExtractor() {
    if (d_statistics_) cudaFree(d_statistics_);
    if (d_unique_minimizers_) cudaFree(d_unique_minimizers_);
}

bool MinimizerFeatureExtractor::process_first_pass(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids) {
    
    std::cout << "First pass: Collecting minimizer statistics..." << std::endl;
    
    // Collect statistics on host first (can be optimized with GPU reduction)
    std::vector<GPUMinimizerHit> h_hits(num_hits);
    cudaMemcpy(h_hits.data(), d_minimizer_hits, num_hits * sizeof(GPUMinimizerHit), 
               cudaMemcpyDeviceToHost);
    
    // Build statistics map
    for (const auto& hit : h_hits) {
        auto& stats = minimizer_stats_map_[hit.minimizer_hash];
        stats.total_occurrences++;
        stats.taxon_occurrences[hit.taxon_id]++;
        stats.unique_taxa = stats.taxon_occurrences.size();
        stats.position_sum += hit.position;
        
        // Track genome occurrence
        uint32_t genome_taxon = (hit.genome_id < genome_taxon_ids.size()) ? 
                                genome_taxon_ids[hit.genome_id] : hit.taxon_id;
        stats.genome_taxa.insert(genome_taxon);
    }
    
    total_genomes_processed_ += genome_taxon_ids.size();
    
    // Calculate basic features during first pass
    calculate_basic_features(d_minimizer_hits, num_hits);
    
    return true;
}

bool MinimizerFeatureExtractor::process_second_pass(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits) {
    
    std::cout << "Second pass: Computing advanced features..." << std::endl;
    
    // Convert statistics map to arrays for GPU
    std::vector<uint64_t> unique_minimizers;
    std::vector<MinimizerFeatureStats> statistics;
    
    for (const auto& [hash, stats] : minimizer_stats_map_) {
        unique_minimizers.push_back(hash);
        statistics.push_back(stats);
    }
    
    num_unique_minimizers_ = unique_minimizers.size();
    
    // Sort minimizers for binary search on GPU
    std::sort(unique_minimizers.begin(), unique_minimizers.end());
    
    // Copy to GPU
    cudaMemcpy(d_unique_minimizers_, unique_minimizers.data(), 
               num_unique_minimizers_ * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_statistics_, statistics.data(), 
               num_unique_minimizers_ * sizeof(MinimizerFeatureStats), cudaMemcpyHostToDevice);
    
    // Launch feature computation kernel
    int block_size = 256;
    int grid_size = (num_hits + block_size - 1) / block_size;
    
    compute_minimizer_features_kernel<<<grid_size, block_size>>>(
        d_minimizer_hits, d_statistics_, d_unique_minimizers_,
        num_hits, num_unique_minimizers_, (float)total_genomes_processed_
    );
    
    cudaDeviceSynchronize();
    
    // Mark position bias
    mark_position_bias_kernel<<<grid_size, block_size>>>(
        d_minimizer_hits, num_hits, 1000  // 1kb clustering window
    );
    
    cudaDeviceSynchronize();
    
    // Mark contamination
    mark_contamination_features(d_minimizer_hits, num_hits);
    
    return true;
}

float MinimizerFeatureExtractor::calculate_uniqueness_score(
    uint64_t minimizer_hash) const {
    
    auto it = minimizer_stats_map_.find(minimizer_hash);
    if (it == minimizer_stats_map_.end()) {
        return 1.0f; // Never seen = unique
    }
    
    const auto& stats = it->second;
    float frequency = (float)stats.total_occurrences / total_genomes_processed_;
    return 1.0f - std::min(frequency, 1.0f);
}

float MinimizerFeatureExtractor::calculate_conservation_score(
    uint64_t minimizer_hash) const {
    
    auto it = minimizer_stats_map_.find(minimizer_hash);
    if (it == minimizer_stats_map_.end()) {
        return 0.0f; // Never seen = not conserved
    }
    
    const auto& stats = it->second;
    
    // Conservation based on taxonomic spread vs occurrences
    float taxonomic_diversity = (float)stats.unique_taxa / stats.total_occurrences;
    
    // Also consider phylogenetic spread if available
    if (!stats.genome_taxa.empty()) {
        // Simple metric: ratio of unique taxa to total occurrences
        float genome_diversity = (float)stats.genome_taxa.size() / stats.total_occurrences;
        taxonomic_diversity = (taxonomic_diversity + genome_diversity) / 2.0f;
    }
    
    return taxonomic_diversity;
}

bool MinimizerFeatureExtractor::identify_contamination_markers(
    const std::vector<GPUMinimizerHit>& minimizer_hits,
    std::vector<bool>& contamination_flags) const {
    
    contamination_flags.resize(minimizer_hits.size(), false);
    
    for (size_t i = 0; i < minimizer_hits.size(); i++) {
        uint64_t hash = minimizer_hits[i].minimizer_hash;
        
        // Check against known contamination patterns
        if (human_minimizer_hashes_.find(hash) != human_minimizer_hashes_.end() ||
            adapter_minimizer_hashes_.find(hash) != adapter_minimizer_hashes_.end()) {
            contamination_flags[i] = true;
        }
    }
    
    return true;
}

FeatureExtractionStats MinimizerFeatureExtractor::get_statistics() const {
    FeatureExtractionStats stats;
    stats.total_unique_minimizers = minimizer_stats_map_.size();
    stats.total_genomes_processed = total_genomes_processed_;
    
    // Calculate distribution statistics
    for (const auto& [hash, mstats] : minimizer_stats_map_) {
        if (mstats.total_occurrences == 1) {
            stats.singleton_minimizers++;
        } else if (mstats.total_occurrences < 10) {
            stats.rare_minimizers++;
        } else if (mstats.total_occurrences > total_genomes_processed_ * 0.5) {
            stats.common_minimizers++;
        }
        
        stats.avg_taxonomic_spread += mstats.unique_taxa;
    }
    
    if (!minimizer_stats_map_.empty()) {
        stats.avg_taxonomic_spread /= minimizer_stats_map_.size();
    }
    
    return stats;
}

void MinimizerFeatureExtractor::export_feature_statistics(const std::string& output_file) const {
    std::ofstream out(output_file);
    if (!out.is_open()) return;
    
    out << "minimizer_hash\ttotal_occurrences\tunique_taxa\tuniqueness_score\tconservation_score\n";
    
    for (const auto& [hash, stats] : minimizer_stats_map_) {
        float uniqueness = calculate_uniqueness_score(hash);
        float conservation = calculate_conservation_score(hash);
        
        out << hash << "\t" 
            << stats.total_occurrences << "\t"
            << stats.unique_taxa << "\t"
            << uniqueness << "\t"
            << conservation << "\n";
    }
    
    out.close();
}

// Private methods

void MinimizerFeatureExtractor::calculate_basic_features(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits) {
    
    // This would calculate GC content and complexity features
    // For now, using placeholder implementation
    
    int block_size = 256;
    int grid_size = (num_hits + block_size - 1) / block_size;
    
    // Launch kernel to calculate sequence-based features
    // (GC content and complexity are already set in extraction phase)
}

void MinimizerFeatureExtractor::mark_contamination_features(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits) {
    
    // For efficiency, this would use a GPU-based hash table lookup
    // Simplified version: copy to host, check, copy back
    
    std::vector<GPUMinimizerHit> h_hits(num_hits);
    cudaMemcpy(h_hits.data(), d_minimizer_hits, num_hits * sizeof(GPUMinimizerHit),
               cudaMemcpyDeviceToHost);
    
    for (auto& hit : h_hits) {
        if (human_minimizer_hashes_.find(hit.minimizer_hash) != human_minimizer_hashes_.end() ||
            adapter_minimizer_hashes_.find(hit.minimizer_hash) != adapter_minimizer_hashes_.end()) {
            hit.feature_flags |= (1 << 7); // Set contamination bit
        }
    }
    
    cudaMemcpy(d_minimizer_hits, h_hits.data(), num_hits * sizeof(GPUMinimizerHit),
               cudaMemcpyHostToDevice);
}

void MinimizerFeatureExtractor::initialize_contamination_patterns() {
    // In production, these would be loaded from files
    // Example placeholders for common adapter sequences
    
    // Common Illumina adapters (simplified - would compute all k-mers)
    adapter_minimizer_hashes_.insert(0x1234567890ABCDEF); // Placeholder
    adapter_minimizer_hashes_.insert(0xFEDCBA0987654321); // Placeholder
    
    // Human-associated minimizers (would be loaded from database)
    human_minimizer_hashes_.insert(0xAAAAAAAAAAAAAAAA); // Placeholder
    human_minimizer_hashes_.insert(0xBBBBBBBBBBBBBBBB); // Placeholder
}

void MinimizerFeatureExtractor::load_contamination_database(const std::string& contamination_db_path) {
    std::ifstream db_file(contamination_db_path);
    if (!db_file.is_open()) {
        std::cerr << "Warning: Could not open contamination database: " << contamination_db_path << std::endl;
        return;
    }
    
    std::string line;
    while (std::getline(db_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        // Expected format: TYPE<TAB>HASH
        size_t tab_pos = line.find('\t');
        if (tab_pos != std::string::npos) {
            std::string type = line.substr(0, tab_pos);
            uint64_t hash = std::stoull(line.substr(tab_pos + 1), nullptr, 16);
            
            if (type == "HUMAN") {
                human_minimizer_hashes_.insert(hash);
            } else if (type == "ADAPTER") {
                adapter_minimizer_hashes_.insert(hash);
            }
        }
    }
    
    db_file.close();
    std::cout << "Loaded contamination database: " 
              << human_minimizer_hashes_.size() << " human markers, "
              << adapter_minimizer_hashes_.size() << " adapter markers" << std::endl;
}

void MinimizerFeatureExtractor::set_contamination_patterns(
    const std::unordered_set<uint64_t>& human_hashes,
    const std::unordered_set<uint64_t>& adapter_hashes) {
    human_minimizer_hashes_ = human_hashes;
    adapter_minimizer_hashes_ = adapter_hashes;
}

// ===========================
// Integration Helper Functions
// ===========================

namespace FeatureExtractionUtils {
    
    bool run_two_pass_feature_extraction(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits,
        const std::vector<uint32_t>& genome_taxon_ids,
        size_t max_minimizers) {
        
        MinimizerFeatureExtractor extractor(max_minimizers, genome_taxon_ids.size());
        
        // First pass
        if (!extractor.process_first_pass(d_minimizer_hits, num_hits, genome_taxon_ids)) {
            return false;
        }
        
        // Second pass
        if (!extractor.process_second_pass(d_minimizer_hits, num_hits)) {
            return false;
        }
        
        // Export statistics
        auto stats = extractor.get_statistics();
        std::cout << "Feature extraction complete:" << std::endl;
        std::cout << "  Unique minimizers: " << stats.total_unique_minimizers << std::endl;
        std::cout << "  Singleton minimizers: " << stats.singleton_minimizers << std::endl;
        std::cout << "  Average taxonomic spread: " << stats.avg_taxonomic_spread << std::endl;
        
        return true;
    }
    
    void print_feature_distribution(const std::vector<GPUMinimizerHit>& minimizer_hits) {
        std::unordered_map<int, int> gc_distribution;
        std::unordered_map<int, int> complexity_distribution;
        int contamination_count = 0;
        int position_bias_count = 0;
        
        for (const auto& hit : minimizer_hits) {
            // Extract GC category (bits 0-2)
            int gc_cat = hit.feature_flags & 0x7;
            gc_distribution[gc_cat]++;
            
            // Extract complexity (bits 3-5)
            int complexity = (hit.feature_flags >> 3) & 0x7;
            complexity_distribution[complexity]++;
            
            // Check contamination (bit 7)
            if (hit.feature_flags & (1 << 7)) {
                contamination_count++;
            }
            
            // Check position bias (bit 6)
            if (hit.feature_flags & (1 << 6)) {
                position_bias_count++;
            }
        }
        
        std::cout << "\nFeature Distribution:" << std::endl;
        std::cout << "GC Content Categories:" << std::endl;
        for (const auto& [cat, count] : gc_distribution) {
            std::cout << "  Category " << cat << ": " << count << std::endl;
        }
        
        std::cout << "Complexity Categories:" << std::endl;
        for (const auto& [cat, count] : complexity_distribution) {
            std::cout << "  Category " << cat << ": " << count << std::endl;
        }
        
        std::cout << "Contamination markers: " << contamination_count << std::endl;
        std::cout << "Position-biased minimizers: " << position_bias_count << std::endl;
    }
}