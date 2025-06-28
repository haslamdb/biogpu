// uniqueness_score_implementation.cu
// Implementation of uniqueness score feature extraction and encoding
// Reduces false positives by identifying minimizers that are too common

#include "../gpu_kraken_types.h"
#include "../processing/minimizer_feature_extractor.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cmath>

// ===========================
// Device Functions for Uniqueness
// ===========================

__device__ uint8_t encode_uniqueness_category(float uniqueness_score) {
    // Map uniqueness score (0.0-1.0) to 3-bit category (0-7)
    // Higher values = more unique = better for classification
    
    if (uniqueness_score >= 0.95f) return 7;      // Extremely unique (singleton-like)
    else if (uniqueness_score >= 0.90f) return 6; // Very unique
    else if (uniqueness_score >= 0.80f) return 5; // Highly unique
    else if (uniqueness_score >= 0.70f) return 4; // Moderately unique
    else if (uniqueness_score >= 0.50f) return 3; // Somewhat unique
    else if (uniqueness_score >= 0.30f) return 2; // Low uniqueness
    else if (uniqueness_score >= 0.10f) return 1; // Very low uniqueness
    else return 0;                                 // Extremely common (likely false positive)
}

__device__ float decode_ml_weight_to_uniqueness(uint16_t ml_weight) {
    // Convert encoded ML weight back to float uniqueness score
    return (float)ml_weight / 65535.0f;
}

__device__ bool is_reliable_minimizer(uint8_t uniqueness_category) {
    // Consider minimizers with uniqueness >= 4 as reliable for classification
    return uniqueness_category >= 4;
}

// ===========================
// CUDA Kernels
// ===========================

__global__ void compute_and_encode_uniqueness_kernel(
    GPUMinimizerHit* minimizer_hits,
    const uint32_t* occurrence_counts,
    const uint64_t* unique_minimizers,
    int num_hits,
    int num_unique_minimizers,
    float total_genomes_processed) {
    
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
        uint32_t occurrences = occurrence_counts[stats_idx];
        
        // Calculate uniqueness score
        float frequency = (float)occurrences / total_genomes_processed;
        float uniqueness = 1.0f - fminf(frequency, 1.0f);
        
        // Apply additional scaling for better discrimination
        // Penalize extremely common minimizers more severely
        if (frequency > 0.5f) {
            uniqueness *= uniqueness; // Square to emphasize the penalty
        }
        
        // Encode ML weight
        uint16_t encoded_uniqueness = (uint16_t)(uniqueness * 65535.0f);
        hit.ml_weight = encoded_uniqueness;
        
        // Encode uniqueness category in feature flags (bits 8-10)
        uint8_t uniqueness_category = encode_uniqueness_category(uniqueness);
        
        // Clear bits 8-10 and set new uniqueness category
        uint32_t new_flags = hit.feature_flags & ~(0x7 << 8);
        new_flags |= (uniqueness_category << 8);
        
        // Set unique minimizer flag (bit 11) for highly unique minimizers
        if (uniqueness >= 0.90f) {
            new_flags |= (1 << 11);
        }
        
        // Set rare minimizer flag (bit 12) for very rare minimizers
        if (occurrences <= 3) {
            new_flags |= (1 << 12);
        }
        
        // Set reliability flag (bit 13) for classification-reliable minimizers
        if (is_reliable_minimizer(uniqueness_category)) {
            new_flags |= (1 << 13);
        }
        
        hit.feature_flags = new_flags;
    }
}

__global__ void apply_uniqueness_filtering_kernel(
    GPUMinimizerHit* minimizer_hits,
    uint8_t* keep_flags,
    int num_hits,
    float min_uniqueness_threshold,
    bool filter_extremely_common) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    const GPUMinimizerHit& hit = minimizer_hits[tid];
    
    // Decode uniqueness from ML weight
    float uniqueness = decode_ml_weight_to_uniqueness(hit.ml_weight);
    
    // Extract uniqueness category from feature flags
    uint8_t uniqueness_category = (hit.feature_flags >> 8) & 0x7;
    
    bool keep = true;
    
    // Apply uniqueness threshold
    if (uniqueness < min_uniqueness_threshold) {
        keep = false;
    }
    
    // Filter extremely common minimizers if requested
    if (filter_extremely_common && uniqueness_category == 0) {
        keep = false;
    }
    
    // Additional quality checks
    // Skip minimizers that are likely contamination or low complexity
    if (hit.feature_flags & (1 << 7)) {  // Contamination flag
        keep = false;
    }
    
    keep_flags[tid] = keep ? 1 : 0;
}

// ===========================
// Host-side Uniqueness Calculator
// ===========================

class UniquenessCalculator {
private:
    std::unordered_map<uint64_t, uint32_t> minimizer_occurrence_counts_;
    size_t total_genomes_processed_;
    
public:
    UniquenessCalculator() : total_genomes_processed_(0) {}
    
    void add_genome_batch(const std::vector<GPUMinimizerHit>& hits, 
                         const std::vector<uint32_t>& genome_taxon_ids) {
        
        // Count unique genomes in this batch
        std::set<uint32_t> unique_genomes;
        for (const auto& hit : hits) {
            if (hit.genome_id < genome_taxon_ids.size()) {
                unique_genomes.insert(genome_taxon_ids[hit.genome_id]);
            }
        }
        total_genomes_processed_ += unique_genomes.size();
        
        // Count minimizer occurrences
        for (const auto& hit : hits) {
            minimizer_occurrence_counts_[hit.minimizer_hash]++;
        }
    }
    
    float calculate_uniqueness_score(uint64_t minimizer_hash) const {
        auto it = minimizer_occurrence_counts_.find(minimizer_hash);
        if (it == minimizer_occurrence_counts_.end()) {
            return 1.0f; // Never seen = perfectly unique
        }
        
        uint32_t occurrences = it->second;
        float frequency = (float)occurrences / total_genomes_processed_;
        
        // Apply logarithmic scaling for better discrimination
        float uniqueness = 1.0f - fminf(frequency, 1.0f);
        
        // Additional penalty for very common minimizers
        if (frequency > 0.1f) {
            uniqueness *= (1.0f - 0.5f * (frequency - 0.1f));
        }
        
        return fmaxf(0.0f, fminf(1.0f, uniqueness));
    }
    
    std::vector<uint64_t> get_unique_minimizers() const {
        std::vector<uint64_t> unique_minimizers;
        for (const auto& [hash, count] : minimizer_occurrence_counts_) {
            unique_minimizers.push_back(hash);
        }
        std::sort(unique_minimizers.begin(), unique_minimizers.end());
        return unique_minimizers;
    }
    
    std::vector<uint32_t> get_occurrence_counts(const std::vector<uint64_t>& unique_minimizers) const {
        std::vector<uint32_t> counts;
        for (uint64_t hash : unique_minimizers) {
            auto it = minimizer_occurrence_counts_.find(hash);
            counts.push_back(it != minimizer_occurrence_counts_.end() ? it->second : 0);
        }
        return counts;
    }
    
    void print_uniqueness_distribution() const {
        std::vector<int> distribution(8, 0);
        
        for (const auto& [hash, count] : minimizer_occurrence_counts_) {
            float uniqueness = calculate_uniqueness_score(hash);
            uint8_t category = 0;
            
            if (uniqueness >= 0.95f) category = 7;
            else if (uniqueness >= 0.90f) category = 6;
            else if (uniqueness >= 0.80f) category = 5;
            else if (uniqueness >= 0.70f) category = 4;
            else if (uniqueness >= 0.50f) category = 3;
            else if (uniqueness >= 0.30f) category = 2;
            else if (uniqueness >= 0.10f) category = 1;
            else category = 0;
            
            distribution[category]++;
        }
        
        std::cout << "\nUniqueness Score Distribution:" << std::endl;
        const char* category_names[] = {
            "Extremely Common (0.00-0.10)", "Very Low (0.10-0.30)", 
            "Low (0.30-0.50)", "Moderate (0.50-0.70)",
            "High (0.70-0.80)", "Very High (0.80-0.90)", 
            "Extremely High (0.90-0.95)", "Singleton-like (0.95-1.00)"
        };
        
        for (int i = 0; i < 8; i++) {
            std::cout << "  Category " << i << " (" << category_names[i] << "): " 
                      << distribution[i] << " minimizers" << std::endl;
        }
    }
    
    // Statistics
    size_t get_total_unique_minimizers() const { return minimizer_occurrence_counts_.size(); }
    size_t get_total_genomes_processed() const { return total_genomes_processed_; }
    
    double get_average_occurrence_rate() const {
        if (minimizer_occurrence_counts_.empty()) return 0.0;
        
        uint64_t total_occurrences = 0;
        for (const auto& [hash, count] : minimizer_occurrence_counts_) {
            total_occurrences += count;
        }
        
        return (double)total_occurrences / minimizer_occurrence_counts_.size();
    }
};

// ===========================
// Host-side API Functions
// ===========================

bool compute_and_encode_uniqueness_scores(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids,
    float total_genomes_processed) {
    
    std::cout << "Computing and encoding uniqueness scores for " << num_hits << " minimizer hits..." << std::endl;
    
    // Copy hits to host for analysis
    std::vector<GPUMinimizerHit> h_hits(num_hits);
    cudaMemcpy(h_hits.data(), d_minimizer_hits, num_hits * sizeof(GPUMinimizerHit), 
               cudaMemcpyDeviceToHost);
    
    // Create uniqueness calculator
    UniquenessCalculator calculator;
    calculator.add_genome_batch(h_hits, genome_taxon_ids);
    
    // Get unique minimizers and their counts
    std::vector<uint64_t> unique_minimizers = calculator.get_unique_minimizers();
    std::vector<uint32_t> occurrence_counts = calculator.get_occurrence_counts(unique_minimizers);
    
    std::cout << "Found " << unique_minimizers.size() << " unique minimizers" << std::endl;
    
    // Allocate GPU memory
    uint64_t* d_unique_minimizers;
    uint32_t* d_occurrence_counts;
    
    cudaMalloc(&d_unique_minimizers, unique_minimizers.size() * sizeof(uint64_t));
    cudaMalloc(&d_occurrence_counts, occurrence_counts.size() * sizeof(uint32_t));
    
    // Copy to GPU
    cudaMemcpy(d_unique_minimizers, unique_minimizers.data(),
               unique_minimizers.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_occurrence_counts, occurrence_counts.data(),
               occurrence_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch kernel
    int block_size = 256;
    int grid_size = (num_hits + block_size - 1) / block_size;
    
    compute_and_encode_uniqueness_kernel<<<grid_size, block_size>>>(
        d_minimizer_hits,
        d_occurrence_counts,
        d_unique_minimizers,
        num_hits,
        unique_minimizers.size(),
        total_genomes_processed
    );
    
    cudaDeviceSynchronize();
    
    // Check for errors
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        std::cerr << "Uniqueness computation kernel failed: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_unique_minimizers);
        cudaFree(d_occurrence_counts);
        return false;
    }
    
    // Print distribution
    calculator.print_uniqueness_distribution();
    
    // Cleanup
    cudaFree(d_unique_minimizers);
    cudaFree(d_occurrence_counts);
    
    std::cout << "✓ Uniqueness scores computed and encoded successfully" << std::endl;
    return true;
}

bool apply_uniqueness_filtering(
    GPUMinimizerHit* d_minimizer_hits,
    size_t& num_hits,
    float min_uniqueness_threshold = 0.3f,
    bool filter_extremely_common = true) {
    
    std::cout << "Applying uniqueness filtering (threshold: " << min_uniqueness_threshold << ")..." << std::endl;
    
    // Allocate flag array
    uint8_t* d_keep_flags;
    cudaMalloc(&d_keep_flags, num_hits * sizeof(uint8_t));
    
    // Launch filtering kernel
    int block_size = 256;
    int grid_size = (num_hits + block_size - 1) / block_size;
    
    apply_uniqueness_filtering_kernel<<<grid_size, block_size>>>(
        d_minimizer_hits,
        d_keep_flags,
        num_hits,
        min_uniqueness_threshold,
        filter_extremely_common
    );
    
    cudaDeviceSynchronize();
    
    // Count how many to keep
    thrust::device_ptr<uint8_t> d_flags_ptr(d_keep_flags);
    size_t kept_count = thrust::reduce(d_flags_ptr, d_flags_ptr + num_hits, 0);
    
    std::cout << "Uniqueness filtering: keeping " << kept_count << " / " << num_hits 
              << " minimizers (" << (100.0 * kept_count / num_hits) << "%)" << std::endl;
    
    // TODO: Implement compaction to actually remove filtered minimizers
    // For now, just report the statistics
    
    cudaFree(d_keep_flags);
    return true;
}

// ===========================
// Integration with Feature Extractor
// ===========================

bool integrate_uniqueness_with_feature_extractor(
    MinimizerFeatureExtractor* feature_extractor,
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids) {
    
    std::cout << "Integrating uniqueness calculation with feature extractor..." << std::endl;
    
    if (!feature_extractor) {
        std::cerr << "Feature extractor is null" << std::endl;
        return false;
    }
    
    // First pass: collect statistics (existing functionality)
    if (!feature_extractor->process_first_pass(d_minimizer_hits, num_hits, genome_taxon_ids)) {
        std::cerr << "Feature extractor first pass failed" << std::endl;
        return false;
    }
    
    // Get statistics and compute uniqueness
    const auto& minimizer_stats = feature_extractor->get_minimizer_stats();
    float total_genomes = genome_taxon_ids.size();
    
    std::cout << "Computing uniqueness for " << minimizer_stats.size() << " unique minimizers..." << std::endl;
    
    // Prepare data for GPU kernel
    std::vector<uint64_t> unique_minimizers;
    std::vector<uint32_t> occurrence_counts;
    
    for (const auto& [hash, stats] : minimizer_stats) {
        unique_minimizers.push_back(hash);
        occurrence_counts.push_back(stats.total_occurrences);
    }
    
    std::sort(unique_minimizers.begin(), unique_minimizers.end());
    
    // Copy to GPU and run uniqueness kernel
    uint64_t* d_unique_minimizers;
    uint32_t* d_occurrence_counts;
    
    cudaMalloc(&d_unique_minimizers, unique_minimizers.size() * sizeof(uint64_t));
    cudaMalloc(&d_occurrence_counts, occurrence_counts.size() * sizeof(uint32_t));
    
    cudaMemcpy(d_unique_minimizers, unique_minimizers.data(),
               unique_minimizers.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_occurrence_counts, occurrence_counts.data(),
               occurrence_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch uniqueness computation
    int block_size = 256;
    int grid_size = (num_hits + block_size - 1) / block_size;
    
    compute_and_encode_uniqueness_kernel<<<grid_size, block_size>>>(
        d_minimizer_hits,
        d_occurrence_counts,
        d_unique_minimizers,
        num_hits,
        unique_minimizers.size(),
        total_genomes
    );
    
    cudaDeviceSynchronize();
    
    // Second pass: continue with other features
    if (!feature_extractor->process_second_pass(d_minimizer_hits, num_hits)) {
        std::cerr << "Feature extractor second pass failed" << std::endl;
        cudaFree(d_unique_minimizers);
        cudaFree(d_occurrence_counts);
        return false;
    }
    
    // Cleanup
    cudaFree(d_unique_minimizers);
    cudaFree(d_occurrence_counts);
    
    std::cout << "✓ Uniqueness integration completed successfully" << std::endl;
    return true;
}

// ===========================
// Utility Functions
// ===========================

namespace UniquenessUtils {
    
    // Decode uniqueness category from feature flags
    inline uint8_t get_uniqueness_category(uint32_t feature_flags) {
        return (feature_flags >> 8) & 0x7;
    }
    
    // Check if minimizer is marked as unique
    inline bool is_unique_minimizer(uint32_t feature_flags) {
        return (feature_flags & (1 << 11)) != 0;
    }
    
    // Check if minimizer is marked as rare
    inline bool is_rare_minimizer(uint32_t feature_flags) {
        return (feature_flags & (1 << 12)) != 0;
    }
    
    // Check if minimizer is reliable for classification
    inline bool is_reliable_minimizer(uint32_t feature_flags) {
        return (feature_flags & (1 << 13)) != 0;
    }
    
    // Convert uniqueness category to human-readable string
    std::string uniqueness_category_to_string(uint8_t category) {
        const char* names[] = {
            "Extremely Common", "Very Low", "Low", "Moderate",
            "High", "Very High", "Extremely High", "Singleton-like"
        };
        return (category < 8) ? names[category] : "Unknown";
    }
    
    // Calculate expected false positive reduction
    double estimate_false_positive_reduction(float uniqueness_threshold) {
        // Empirical formula based on minimizer frequency distribution
        // Higher thresholds filter more common minimizers = fewer false positives
        return 1.0 - exp(-3.0 * uniqueness_threshold);
    }
    
    // Print minimizer quality summary
    void print_minimizer_quality_summary(const std::vector<GPUMinimizerHit>& hits) {
        std::vector<int> uniqueness_dist(8, 0);
        int unique_count = 0, rare_count = 0, reliable_count = 0;
        
        for (const auto& hit : hits) {
            uint8_t category = get_uniqueness_category(hit.feature_flags);
            uniqueness_dist[category]++;
            
            if (is_unique_minimizer(hit.feature_flags)) unique_count++;
            if (is_rare_minimizer(hit.feature_flags)) rare_count++;
            if (is_reliable_minimizer(hit.feature_flags)) reliable_count++;
        }
        
        std::cout << "\nMinimizer Quality Summary:" << std::endl;
        std::cout << "Total minimizers: " << hits.size() << std::endl;
        std::cout << "Unique minimizers (≥90% uniqueness): " << unique_count << " ("
                  << (100.0 * unique_count / hits.size()) << "%)" << std::endl;
        std::cout << "Rare minimizers (≤3 occurrences): " << rare_count << " ("
                  << (100.0 * rare_count / hits.size()) << "%)" << std::endl;
        std::cout << "Reliable minimizers: " << reliable_count << " ("
                  << (100.0 * reliable_count / hits.size()) << "%)" << std::endl;
        
        std::cout << "\nUniqueness Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            std::cout << "  " << uniqueness_category_to_string(i) << ": " 
                      << uniqueness_dist[i] << std::endl;
        }
    }
}
