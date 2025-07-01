// working_uniqueness_implementation.cu
// Working implementation that integrates with existing codebase
// Avoids namespace conflicts and provides the actual GPU kernels

#include "namespace_conflict_resolution.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <iomanip>

// ===========================
// Device Functions (GPU)
// ===========================

__device__ uint8_t encode_uniqueness_category_gpu(float uniqueness_score) {
    if (uniqueness_score >= 0.95f) return 7;      // Singleton-like
    else if (uniqueness_score >= 0.90f) return 6; // Extremely high
    else if (uniqueness_score >= 0.80f) return 5; // Very high
    else if (uniqueness_score >= 0.70f) return 4; // High
    else if (uniqueness_score >= 0.50f) return 3; // Moderate
    else if (uniqueness_score >= 0.30f) return 2; // Low
    else if (uniqueness_score >= 0.10f) return 1; // Very low
    else return 0;                                 // Extremely common
}

// ===========================
// GPU Kernels
// ===========================

__global__ void compute_uniqueness_kernel(
    GPUMinimizerHit* minimizer_hits,
    const uint64_t* unique_minimizer_hashes,
    const uint32_t* occurrence_counts,
    int num_hits,
    int num_unique_minimizers,
    float total_genomes) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    GPUMinimizerHit& hit = minimizer_hits[tid];
    uint64_t minimizer_hash = hit.minimizer_hash;
    
    // Binary search to find this minimizer's index
    int left = 0, right = num_unique_minimizers - 1;
    int found_idx = -1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        if (unique_minimizer_hashes[mid] == minimizer_hash) {
            found_idx = mid;
            break;
        } else if (unique_minimizer_hashes[mid] < minimizer_hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    if (found_idx >= 0) {
        uint32_t occurrences = occurrence_counts[found_idx];
        
        // Calculate uniqueness score
        float frequency = (float)occurrences / total_genomes;
        float uniqueness = 1.0f - fminf(frequency, 1.0f);
        
        // Apply logarithmic scaling for better discrimination
        if (frequency > 0.1f) {
            uniqueness *= (1.0f - 0.3f * (frequency - 0.1f));
        }
        
        // Encode into ML weight
        hit.ml_weight = (uint16_t)(fmaxf(0.0f, fminf(1.0f, uniqueness)) * 65535.0f);
        
        // Encode uniqueness category
        uint8_t category = encode_uniqueness_category_gpu(uniqueness);
        uint32_t new_flags = hit.feature_flags & ~(0x7 << 8);  // Clear bits 8-10
        new_flags |= (category << 8);  // Set uniqueness category
        
        // Set quality flags
        if (uniqueness >= 0.90f) {
            new_flags |= (1 << 11);  // Unique flag
        }
        if (occurrences <= 3) {
            new_flags |= (1 << 12);  // Rare flag
        }
        if (category >= 4) {  // High uniqueness or better
            new_flags |= (1 << 13);  // Reliable flag
        }
        
        hit.feature_flags = new_flags;
    }
}

__global__ void count_minimizer_occurrences_kernel(
    const GPUMinimizerHit* minimizer_hits,
    const uint64_t* unique_minimizer_hashes,
    uint32_t* occurrence_counts,
    int num_hits,
    int num_unique_minimizers) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    uint64_t minimizer_hash = minimizer_hits[tid].minimizer_hash;
    
    // Binary search to find index
    int left = 0, right = num_unique_minimizers - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        if (unique_minimizer_hashes[mid] == minimizer_hash) {
            atomicAdd(&occurrence_counts[mid], 1);
            break;
        } else if (unique_minimizer_hashes[mid] < minimizer_hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
}

// ===========================
// Host Implementation
// ===========================

class UniquenessProcessor {
private:
    std::unordered_map<uint64_t, uint32_t> occurrence_map_;
    size_t total_genomes_;
    
public:
    UniquenessProcessor() : total_genomes_(0) {}
    
    void process_minimizer_batch(const std::vector<GPUMinimizerHit>& hits,
                                const std::vector<uint32_t>& genome_taxon_ids) {
        // Count unique genomes in this batch
        std::set<uint32_t> unique_genomes;
        for (const auto& hit : hits) {
            if (hit.genome_id < genome_taxon_ids.size()) {
                unique_genomes.insert(genome_taxon_ids[hit.genome_id]);
            }
        }
        total_genomes_ += unique_genomes.size();
        
        // Count occurrences
        for (const auto& hit : hits) {
            occurrence_map_[hit.minimizer_hash]++;
        }
    }
    
    bool compute_gpu_uniqueness(GPUMinimizerHit* d_minimizer_hits, size_t num_hits) {
        std::cout << "Computing uniqueness scores for " << num_hits << " minimizers..." << std::endl;
        
        // Create sorted arrays for GPU
        std::vector<uint64_t> unique_hashes;
        std::vector<uint32_t> occurrence_counts;
        
        for (const auto& [hash, count] : occurrence_map_) {
            unique_hashes.push_back(hash);
            occurrence_counts.push_back(count);
        }
        
        // Sort by hash for binary search
        std::vector<size_t> indices(unique_hashes.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                  [&](size_t a, size_t b) { return unique_hashes[a] < unique_hashes[b]; });
        
        std::vector<uint64_t> sorted_hashes;
        std::vector<uint32_t> sorted_counts;
        for (size_t idx : indices) {
            sorted_hashes.push_back(unique_hashes[idx]);
            sorted_counts.push_back(occurrence_counts[idx]);
        }
        
        std::cout << "Found " << sorted_hashes.size() << " unique minimizers across " 
                  << total_genomes_ << " genomes" << std::endl;
        
        // Allocate GPU memory
        uint64_t* d_hashes;
        uint32_t* d_counts;
        
        cudaMalloc(&d_hashes, sorted_hashes.size() * sizeof(uint64_t));
        cudaMalloc(&d_counts, sorted_counts.size() * sizeof(uint32_t));
        
        // Copy to GPU
        cudaMemcpy(d_hashes, sorted_hashes.data(), 
                   sorted_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_counts, sorted_counts.data(), 
                   sorted_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Launch kernel
        int block_size = 256;
        int grid_size = (num_hits + block_size - 1) / block_size;
        
        compute_uniqueness_kernel<<<grid_size, block_size>>>(
            d_minimizer_hits,
            d_hashes,
            d_counts,
            num_hits,
            sorted_hashes.size(),
            (float)total_genomes_
        );
        
        cudaDeviceSynchronize();
        
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "Uniqueness kernel failed: " << cudaGetErrorString(error) << std::endl;
            cudaFree(d_hashes);
            cudaFree(d_counts);
            return false;
        }
        
        // Print distribution
        print_uniqueness_distribution(sorted_counts);
        
        cudaFree(d_hashes);
        cudaFree(d_counts);
        
        std::cout << "✓ Uniqueness computation completed successfully" << std::endl;
        return true;
    }
    
    void print_uniqueness_distribution(const std::vector<uint32_t>& counts) {
        std::vector<int> distribution(8, 0);
        
        for (uint32_t count : counts) {
            float frequency = (float)count / total_genomes_;
            float uniqueness = 1.0f - std::min(frequency, 1.0f);
            
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
        
        std::cout << "\nUniqueness Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            if (distribution[i] > 0) {
                std::cout << "  " << MinimizerFlags::uniqueness_category_name(i) 
                          << ": " << distribution[i] << " minimizers" << std::endl;
            }
        }
    }
    
    size_t get_unique_minimizer_count() const { return occurrence_map_.size(); }
    size_t get_total_genomes() const { return total_genomes_; }
};

// ===========================
// Public API Functions
// ===========================

// Global processor instance
static UniquenessProcessor g_uniqueness_processor;

bool compute_and_encode_uniqueness_scores(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids,
    float total_genomes_processed) {
    
    if (num_hits == 0) return true;
    
    std::cout << "Starting uniqueness score computation..." << std::endl;
    
    // Copy hits to host for analysis
    std::vector<GPUMinimizerHit> h_hits(num_hits);
    cudaMemcpy(h_hits.data(), d_minimizer_hits, num_hits * sizeof(GPUMinimizerHit), 
               cudaMemcpyDeviceToHost);
    
    // Process with uniqueness processor
    g_uniqueness_processor.process_minimizer_batch(h_hits, genome_taxon_ids);
    
    // Compute uniqueness on GPU
    bool success = g_uniqueness_processor.compute_gpu_uniqueness(d_minimizer_hits, num_hits);
    
    if (success) {
        std::cout << "Uniqueness scoring completed for " << num_hits << " minimizers" << std::endl;
        std::cout << "  Unique minimizers: " << g_uniqueness_processor.get_unique_minimizer_count() << std::endl;
        std::cout << "  Total genomes: " << g_uniqueness_processor.get_total_genomes() << std::endl;
    }
    
    return success;
}

bool integrate_uniqueness_with_feature_extractor(
    MinimizerFeatureExtractor* feature_extractor,
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids) {
    
    std::cout << "Integrating uniqueness with feature extractor..." << std::endl;
    
    // If no feature extractor, use standalone computation
    if (!feature_extractor) {
        return compute_and_encode_uniqueness_scores(
            d_minimizer_hits, num_hits, genome_taxon_ids, genome_taxon_ids.size()
        );
    }
    
    // Run feature extractor first pass
    if (!feature_extractor->process_first_pass(d_minimizer_hits, num_hits, genome_taxon_ids)) {
        std::cerr << "Feature extractor first pass failed" << std::endl;
        return false;
    }
    
    // Compute uniqueness scores
    if (!compute_and_encode_uniqueness_scores(
            d_minimizer_hits, num_hits, genome_taxon_ids, genome_taxon_ids.size())) {
        std::cerr << "Uniqueness computation failed" << std::endl;
        return false;
    }
    
    // Run feature extractor second pass
    if (!feature_extractor->process_second_pass(d_minimizer_hits, num_hits)) {
        std::cerr << "Feature extractor second pass failed" << std::endl;
        return false;
    }
    
    std::cout << "✓ Uniqueness integration completed successfully" << std::endl;
    return true;
}

// ===========================
// Statistics Collection
// ===========================

struct UniquenessStats {
    size_t total_minimizers = 0;
    size_t unique_minimizers = 0;       // ≥90% uniqueness
    size_t rare_minimizers = 0;         // ≤3 occurrences
    size_t reliable_minimizers = 0;     // Suitable for classification
    double average_uniqueness = 0.0;
    std::vector<int> category_distribution;
    
    UniquenessStats() : category_distribution(8, 0) {}
    
    void print() const {
        std::cout << "\n=== UNIQUENESS STATISTICS ===" << std::endl;
        std::cout << "Total minimizers: " << total_minimizers << std::endl;
        std::cout << "Unique minimizers (≥90%): " << unique_minimizers << " ("
                  << (100.0 * unique_minimizers / total_minimizers) << "%)" << std::endl;
        std::cout << "Rare minimizers (≤3 occurrences): " << rare_minimizers << " ("
                  << (100.0 * rare_minimizers / total_minimizers) << "%)" << std::endl;
        std::cout << "Reliable minimizers: " << reliable_minimizers << " ("
                  << (100.0 * reliable_minimizers / total_minimizers) << "%)" << std::endl;
        std::cout << "Average uniqueness score: " << std::fixed << std::setprecision(3) 
                  << average_uniqueness << std::endl;
        
        std::cout << "\nCategory Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            if (category_distribution[i] > 0) {
                std::cout << "  " << MinimizerFlags::uniqueness_category_name(i) 
                          << ": " << category_distribution[i] << std::endl;
            }
        }
    }
};

UniquenessStats collect_uniqueness_statistics(const std::vector<GPUMinimizerHit>& hits) {
    UniquenessStats stats;
    stats.total_minimizers = hits.size();
    
    double uniqueness_sum = 0.0;
    
    for (const auto& hit : hits) {
        uint8_t category = MinimizerFlags::get_uniqueness_category_safe(hit.feature_flags);
        stats.category_distribution[category]++;
        
        if (MinimizerFlags::is_unique_minimizer_safe(hit.feature_flags)) {
            stats.unique_minimizers++;
        }
        if (MinimizerFlags::is_rare_minimizer_safe(hit.feature_flags)) {
            stats.rare_minimizers++;
        }
        if (MinimizerFlags::is_reliable_minimizer_safe(hit.feature_flags)) {
            stats.reliable_minimizers++;
        }
        
        float uniqueness = MinimizerFlags::ml_weight_to_uniqueness(hit.ml_weight);
        uniqueness_sum += uniqueness;
    }
    
    if (stats.total_minimizers > 0) {
        stats.average_uniqueness = uniqueness_sum / stats.total_minimizers;
    }
    
    return stats;
}

// ===========================
// Testing Function
// ===========================

bool test_uniqueness_implementation() {
    std::cout << "Testing uniqueness implementation..." << std::endl;
    
    // Create test data
    std::vector<GPUMinimizerHit> test_hits;
    GPUMinimizerHit hit1;
    hit1.minimizer_hash = 0x1111111111111111;
    hit1.genome_id = 0;
    hit1.position = 10;
    hit1.strand = 0;  // forward strand
    hit1.taxon_id = 1000;
    hit1.ml_weight = 0;
    hit1.feature_flags = 0;
    test_hits.push_back(hit1);  // Will be common
    
    GPUMinimizerHit hit2;
    hit2.minimizer_hash = 0x2222222222222222;
    hit2.genome_id = 1;
    hit2.position = 20;
    hit2.strand = 0;  // forward strand
    hit2.taxon_id = 1001;
    hit2.ml_weight = 0;
    hit2.feature_flags = 0;
    test_hits.push_back(hit2);  // Will be rare
    
    GPUMinimizerHit hit3;
    hit3.minimizer_hash = 0x3333333333333333;
    hit3.genome_id = 2;
    hit3.position = 30;
    hit3.strand = 0;  // forward strand
    hit3.taxon_id = 1002;
    hit3.ml_weight = 0;
    hit3.feature_flags = 0;
    test_hits.push_back(hit3);  // Will be unique
    
    std::vector<uint32_t> genome_ids = {1000, 1001, 1002};
    
    // Allocate GPU memory
    GPUMinimizerHit* d_hits;
    cudaMalloc(&d_hits, test_hits.size() * sizeof(GPUMinimizerHit));
    cudaMemcpy(d_hits, test_hits.data(), test_hits.size() * sizeof(GPUMinimizerHit), 
               cudaMemcpyHostToDevice);
    
    // Test computation
    bool success = compute_and_encode_uniqueness_scores(d_hits, test_hits.size(), genome_ids, 3.0f);
    
    if (success) {
        // Copy back and check results
        cudaMemcpy(test_hits.data(), d_hits, test_hits.size() * sizeof(GPUMinimizerHit),
                   cudaMemcpyDeviceToHost);
        
        std::cout << "Test results:" << std::endl;
        for (size_t i = 0; i < test_hits.size(); i++) {
            uint8_t category = MinimizerFlags::get_uniqueness_category_safe(test_hits[i].feature_flags);
            bool unique = MinimizerFlags::is_unique_minimizer_safe(test_hits[i].feature_flags);
            bool reliable = MinimizerFlags::is_reliable_minimizer_safe(test_hits[i].feature_flags);
            
            std::cout << "  Hit " << i << ": category=" << (int)category 
                      << ", unique=" << unique << ", reliable=" << reliable << std::endl;
        }
        
        UniquenessStats stats = collect_uniqueness_statistics(test_hits);
        stats.print();
    }
    
    cudaFree(d_hits);
    return success;
}