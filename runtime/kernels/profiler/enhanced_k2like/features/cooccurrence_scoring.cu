#include "../gpu_kraken_types.h"
#include "../features/namespace_conflict_resolution.h"
#include "../features/enhanced_minimizer_flags.h"
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <algorithm>

// Configuration constants
constexpr int COOCCURRENCE_WINDOW_BP = 1000;
constexpr int MAX_COOCCURRENCES_TO_TRACK = 100;  // Per minimizer
constexpr float EXPECTED_RANDOM_COOCCURRENCE_RATE = 0.001f;  // Based on genome size/minimizer density

// Structure for efficient co-occurrence tracking
struct MinimizerCooccurrenceEntry {
    uint64_t minimizer_hash;
    uint16_t genome_id;
    uint32_t position;
    uint32_t hit_index;  // Original index in hits array
};

// Structure for tracking pairwise co-occurrences
struct CooccurrencePair {
    uint64_t hash1;
    uint64_t hash2;
    uint32_t count;
    float distance_sum;  // For weighted scoring
};

// Device function to calculate position-weighted co-occurrence score
__device__ float calculate_distance_weight(int distance, int window_size) {
    // Closer minimizers get higher weight (exponential decay)
    float normalized_distance = (float)distance / window_size;
    return expf(-2.0f * normalized_distance);  // e^(-2x) decay
}

// Kernel to sort minimizer hits by genome and position for efficient window search
__global__ void prepare_sorted_indices_kernel(
    const GPUMinimizerHit* minimizer_hits,
    MinimizerCooccurrenceEntry* sorted_entries,
    int num_hits) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    const GPUMinimizerHit& hit = minimizer_hits[tid];
    MinimizerCooccurrenceEntry& entry = sorted_entries[tid];
    
    entry.minimizer_hash = hit.minimizer_hash;
    entry.genome_id = hit.genome_id;
    entry.position = hit.position;
    entry.hit_index = tid;
}

// Kernel to compute co-occurrence matrix using sorted data
__global__ void build_cooccurrence_matrix_kernel(
    const MinimizerCooccurrenceEntry* sorted_entries,
    const int* genome_boundaries,  // Start index for each genome in sorted array
    CooccurrencePair* cooccurrence_pairs,
    uint32_t* pair_count,
    int num_entries,
    int num_genomes,
    int window_size_bp,
    int max_pairs) {
    
    int genome_id = blockIdx.x;
    if (genome_id >= num_genomes) return;
    
    int thread_id = threadIdx.x;
    int threads_per_block = blockDim.x;
    
    // Get boundaries for this genome
    int genome_start = genome_boundaries[genome_id];
    int genome_end = (genome_id + 1 < num_genomes) ? 
                     genome_boundaries[genome_id + 1] : num_entries;
    
    int genome_size = genome_end - genome_start;
    if (genome_size < 2) return;  // Need at least 2 minimizers
    
    // Each thread processes a subset of positions
    for (int i = genome_start + thread_id; i < genome_end; i += threads_per_block) {
        const MinimizerCooccurrenceEntry& center = sorted_entries[i];
        
        // Binary search for window boundaries
        // Find start of window
        int window_start = i;
        int left = genome_start;
        int right = i - 1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            int distance = center.position - sorted_entries[mid].position;
            if (distance <= window_size_bp) {
                window_start = mid;
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        
        // Find end of window
        int window_end = i;
        left = i + 1;
        right = genome_end - 1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            int distance = sorted_entries[mid].position - center.position;
            if (distance <= window_size_bp) {
                window_end = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        
        // Process all pairs within window
        for (int j = window_start; j <= window_end; j++) {
            if (i == j) continue;
            
            const MinimizerCooccurrenceEntry& neighbor = sorted_entries[j];
            int distance = abs((int)neighbor.position - (int)center.position);
            
            if (distance <= window_size_bp) {
                // Ensure consistent ordering (smaller hash first)
                uint64_t hash1 = min(center.minimizer_hash, neighbor.minimizer_hash);
                uint64_t hash2 = max(center.minimizer_hash, neighbor.minimizer_hash);
                
                // Calculate position-based weight
                float weight = calculate_distance_weight(distance, window_size_bp);
                
                // Atomically add to global pair list
                uint32_t pair_idx = atomicAdd(pair_count, 1);
                if (pair_idx < max_pairs) {
                    CooccurrencePair& pair = cooccurrence_pairs[pair_idx];
                    pair.hash1 = hash1;
                    pair.hash2 = hash2;
                    pair.count = 1;
                    pair.distance_sum = weight;
                }
            }
        }
    }
}

// Kernel to aggregate co-occurrence pairs and compute final scores
__global__ void aggregate_cooccurrence_scores_kernel(
    const CooccurrencePair* all_pairs,
    const uint64_t* unique_minimizers,
    const uint32_t* minimizer_occurrence_counts,
    float* cooccurrence_scores,
    int num_pairs,
    int num_unique_minimizers,
    float total_genome_length,
    int window_size_bp) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_unique_minimizers) return;
    
    uint64_t target_hash = unique_minimizers[tid];
    uint32_t target_occurrences = minimizer_occurrence_counts[tid];
    
    if (target_occurrences == 0) {
        cooccurrence_scores[tid] = 0.0f;
        return;
    }
    
    // Count weighted co-occurrences for this minimizer
    float total_cooccurrence_weight = 0.0f;
    int unique_partners = 0;
    
    // Track unique partner minimizers (simplified - in production use hash set)
    uint64_t recent_partners[32];  // Track last 32 unique partners
    int partner_count = 0;
    
    for (int i = 0; i < num_pairs; i++) {
        const CooccurrencePair& pair = all_pairs[i];
        
        bool involves_target = (pair.hash1 == target_hash || pair.hash2 == target_hash);
        if (!involves_target) continue;
        
        uint64_t partner_hash = (pair.hash1 == target_hash) ? pair.hash2 : pair.hash1;
        
        // Check if this is a new partner
        bool is_new_partner = true;
        for (int j = 0; j < partner_count && j < 32; j++) {
            if (recent_partners[j] == partner_hash) {
                is_new_partner = false;
                break;
            }
        }
        
        if (is_new_partner && partner_count < 32) {
            recent_partners[partner_count++] = partner_hash;
            unique_partners++;
        }
        
        total_cooccurrence_weight += pair.distance_sum;
    }
    
    // Calculate expected random co-occurrences
    float minimizer_density = (float)num_unique_minimizers / total_genome_length;
    float expected_cooccurrences_per_window = minimizer_density * window_size_bp * target_occurrences;
    float expected_random_weight = expected_cooccurrences_per_window * 0.5f;  // Average weight
    
    // Compute enrichment score
    float enrichment = 1.0f;
    if (expected_random_weight > 0.0f) {
        enrichment = total_cooccurrence_weight / expected_random_weight;
    }
    
    // Consider partner diversity (minimizers with diverse partners are more reliable)
    float diversity_factor = 1.0f - expf(-(float)unique_partners / 10.0f);
    
    // Final score combines enrichment and diversity
    float final_score = fminf(1.0f, (enrichment * diversity_factor) / 10.0f);
    
    cooccurrence_scores[tid] = final_score;
}

// Kernel to encode co-occurrence scores into feature flags
__global__ void encode_cooccurrence_in_feature_flags_kernel(
    GPUMinimizerHit* minimizer_hits,
    const uint64_t* unique_minimizers,
    const float* cooccurrence_scores,
    int num_hits,
    int num_unique_minimizers) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    GPUMinimizerHit& hit = minimizer_hits[tid];
    
    // Binary search for this minimizer's score
    int idx = thrust::lower_bound(thrust::seq,
                                  unique_minimizers,
                                  unique_minimizers + num_unique_minimizers,
                                  hit.minimizer_hash) - unique_minimizers;
    
    float score = 0.0f;
    if (idx < num_unique_minimizers && unique_minimizers[idx] == hit.minimizer_hash) {
        score = cooccurrence_scores[idx];
    }
    
    // Encode score into 3-bit category (bits 13-15)
    uint8_t category;
    if (score < 0.05f) category = 0;      // No significant co-occurrence
    else if (score < 0.15f) category = 1; // Very low
    else if (score < 0.30f) category = 2; // Low
    else if (score < 0.50f) category = 3; // Moderate
    else if (score < 0.70f) category = 4; // High
    else if (score < 0.85f) category = 5; // Very high
    else if (score < 0.95f) category = 6; // Extremely high
    else category = 7;                     // Perfect co-occurrence pattern
    
    // Update feature flags (bits 13-15)
    hit.feature_flags = EnhancedMinimizerFlags::set_cooccurrence_category(hit.feature_flags, category);
}

// Host-side wrapper class
class CooccurrenceScorer {
private:
    size_t max_pairs_;
    float total_genome_length_;
    
    // Device memory
    MinimizerCooccurrenceEntry* d_sorted_entries_;
    int* d_genome_boundaries_;
    CooccurrencePair* d_cooccurrence_pairs_;
    uint32_t* d_pair_count_;
    
public:
    CooccurrenceScorer(size_t estimated_hits = 10000000) 
        : max_pairs_(estimated_hits * 10),  // Estimate ~10 co-occurrences per minimizer
          total_genome_length_(0),
          d_sorted_entries_(nullptr),
          d_genome_boundaries_(nullptr),
          d_cooccurrence_pairs_(nullptr),
          d_pair_count_(nullptr) {}
    
    ~CooccurrenceScorer() {
        cleanup();
    }
    
    void cleanup() {
        if (d_sorted_entries_) cudaFree(d_sorted_entries_);
        if (d_genome_boundaries_) cudaFree(d_genome_boundaries_);
        if (d_cooccurrence_pairs_) cudaFree(d_cooccurrence_pairs_);
        if (d_pair_count_) cudaFree(d_pair_count_);
        
        d_sorted_entries_ = nullptr;
        d_genome_boundaries_ = nullptr;
        d_cooccurrence_pairs_ = nullptr;
        d_pair_count_ = nullptr;
    }
    
    bool compute_cooccurrence_scores(
        GPUMinimizerHit* d_minimizer_hits,
        size_t num_hits,
        const std::vector<uint64_t>& unique_minimizers,
        const std::vector<uint32_t>& occurrence_counts,
        const std::vector<GPUGenomeInfo>& genome_info,
        int window_size_bp = COOCCURRENCE_WINDOW_BP) {
        
        std::cout << "Computing co-occurrence scores for " << num_hits << " minimizers..." << std::endl;
        
        // Allocate memory for sorted entries
        cudaMalloc(&d_sorted_entries_, num_hits * sizeof(MinimizerCooccurrenceEntry));
        
        // Prepare sorted entries
        int block_size = 256;
        int grid_size = (num_hits + block_size - 1) / block_size;
        
        prepare_sorted_indices_kernel<<<grid_size, block_size>>>(
            d_minimizer_hits, d_sorted_entries_, num_hits
        );
        
        // Sort by genome ID and position using Thrust
        thrust::device_ptr<MinimizerCooccurrenceEntry> d_entries_ptr(d_sorted_entries_);
        thrust::sort(d_entries_ptr, d_entries_ptr + num_hits,
            [] __device__ (const MinimizerCooccurrenceEntry& a, const MinimizerCooccurrenceEntry& b) {
                if (a.genome_id != b.genome_id) return a.genome_id < b.genome_id;
                return a.position < b.position;
            });
        
        // Find genome boundaries
        std::vector<int> h_genome_boundaries;
        h_genome_boundaries.push_back(0);
        
        std::vector<MinimizerCooccurrenceEntry> h_sorted_entries(num_hits);
        cudaMemcpy(h_sorted_entries.data(), d_sorted_entries_, 
                   num_hits * sizeof(MinimizerCooccurrenceEntry), cudaMemcpyDeviceToHost);
        
        uint16_t current_genome = h_sorted_entries[0].genome_id;
        for (size_t i = 1; i < num_hits; i++) {
            if (h_sorted_entries[i].genome_id != current_genome) {
                h_genome_boundaries.push_back(i);
                current_genome = h_sorted_entries[i].genome_id;
            }
        }
        
        // Calculate total genome length
        total_genome_length_ = 0;
        for (const auto& genome : genome_info) {
            total_genome_length_ += genome.sequence_length;
        }
        
        // Copy genome boundaries to device
        cudaMalloc(&d_genome_boundaries_, h_genome_boundaries.size() * sizeof(int));
        cudaMemcpy(d_genome_boundaries_, h_genome_boundaries.data(),
                   h_genome_boundaries.size() * sizeof(int), cudaMemcpyHostToDevice);
        
        // Allocate memory for co-occurrence pairs
        cudaMalloc(&d_cooccurrence_pairs_, max_pairs_ * sizeof(CooccurrencePair));
        cudaMalloc(&d_pair_count_, sizeof(uint32_t));
        cudaMemset(d_pair_count_, 0, sizeof(uint32_t));
        
        // Build co-occurrence matrix
        int num_genomes = h_genome_boundaries.size();
        dim3 grid(num_genomes);
        dim3 block(256);
        
        build_cooccurrence_matrix_kernel<<<grid, block>>>(
            d_sorted_entries_,
            d_genome_boundaries_,
            d_cooccurrence_pairs_,
            d_pair_count_,
            num_hits,
            num_genomes,
            window_size_bp,
            max_pairs_
        );
        
        // Get actual pair count
        uint32_t h_pair_count;
        cudaMemcpy(&h_pair_count, d_pair_count_, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        std::cout << "Found " << h_pair_count << " co-occurrence pairs" << std::endl;
        
        // Sort and aggregate pairs
        if (h_pair_count > 0) {
            thrust::device_ptr<CooccurrencePair> d_pairs_ptr(d_cooccurrence_pairs_);
            thrust::sort(d_pairs_ptr, d_pairs_ptr + std::min(h_pair_count, (uint32_t)max_pairs_),
                [] __device__ (const CooccurrencePair& a, const CooccurrencePair& b) {
                    if (a.hash1 != b.hash1) return a.hash1 < b.hash1;
                    return a.hash2 < b.hash2;
                });
            
            // Reduce by key to aggregate duplicate pairs
            // (This is simplified - in production use thrust::reduce_by_key)
        }
        
        // Upload unique minimizers and compute scores
        uint64_t* d_unique_minimizers;
        uint32_t* d_occurrence_counts;
        float* d_cooccurrence_scores;
        
        cudaMalloc(&d_unique_minimizers, unique_minimizers.size() * sizeof(uint64_t));
        cudaMalloc(&d_occurrence_counts, occurrence_counts.size() * sizeof(uint32_t));
        cudaMalloc(&d_cooccurrence_scores, unique_minimizers.size() * sizeof(float));
        
        cudaMemcpy(d_unique_minimizers, unique_minimizers.data(),
                   unique_minimizers.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_occurrence_counts, occurrence_counts.data(),
                   occurrence_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Aggregate scores
        grid_size = (unique_minimizers.size() + block_size - 1) / block_size;
        aggregate_cooccurrence_scores_kernel<<<grid_size, block_size>>>(
            d_cooccurrence_pairs_,
            d_unique_minimizers,
            d_occurrence_counts,
            d_cooccurrence_scores,
            std::min(h_pair_count, (uint32_t)max_pairs_),
            unique_minimizers.size(),
            total_genome_length_,
            window_size_bp
        );
        
        // Encode scores in feature flags
        grid_size = (num_hits + block_size - 1) / block_size;
        encode_cooccurrence_in_feature_flags_kernel<<<grid_size, block_size>>>(
            d_minimizer_hits,
            d_unique_minimizers,
            d_cooccurrence_scores,
            num_hits,
            unique_minimizers.size()
        );
        
        cudaDeviceSynchronize();
        
        // Cleanup
        cudaFree(d_unique_minimizers);
        cudaFree(d_occurrence_counts);
        cudaFree(d_cooccurrence_scores);
        cleanup();
        
        std::cout << "✓ Co-occurrence scoring completed" << std::endl;
        return true;
    }
};

// Integration function
bool compute_and_encode_cooccurrence_scores(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint64_t>& unique_minimizers,
    const std::vector<uint32_t>& occurrence_counts,
    const std::vector<GPUGenomeInfo>& genome_info,
    int window_size_bp = COOCCURRENCE_WINDOW_BP) {
    
    CooccurrenceScorer scorer(num_hits);
    return scorer.compute_cooccurrence_scores(
        d_minimizer_hits,
        num_hits,
        unique_minimizers,
        occurrence_counts,
        genome_info,
        window_size_bp
    );
}
