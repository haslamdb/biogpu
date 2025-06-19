// gpu_phylogenetic_calculator.cuh
// GPU-optimized phylogenetic distance calculations

#ifndef GPU_PHYLOGENETIC_CALCULATOR_CUH
#define GPU_PHYLOGENETIC_CALCULATOR_CUH

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// GPU-friendly phylogenetic data structure
struct GPUPhylogeneticData {
    uint32_t* parent_lookup;    // [taxon_id] -> parent_taxon_id
    uint8_t* depth_lookup;      // [taxon_id] -> depth from root
    uint8_t* rank_lookup;       // [taxon_id] -> taxonomic rank level
    uint32_t max_taxon_id;      // Maximum valid taxon ID
};

// Device function to find LCA of two taxa
__device__ inline uint32_t find_lca_gpu_phylo(
    uint32_t taxon1, 
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    uint32_t max_taxon_id) {
    
    if (taxon1 == taxon2) return taxon1;
    if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return 1;  // Root
    if (taxon1 == 0 || taxon2 == 0) return 1;  // Root
    
    // Get depths
    uint8_t depth1 = depth_lookup[taxon1];
    uint8_t depth2 = depth_lookup[taxon2];
    
    // Bring both taxa to the same depth
    uint32_t current1 = taxon1;
    uint32_t current2 = taxon2;
    
    // Move deeper taxon up to match shallower one
    while (depth1 > depth2 && current1 != 0 && current1 <= max_taxon_id) {
        current1 = parent_lookup[current1];
        depth1--;
    }
    
    while (depth2 > depth1 && current2 != 0 && current2 <= max_taxon_id) {
        current2 = parent_lookup[current2];
        depth2--;
    }
    
    // Now find common ancestor
    for (int steps = 0; steps < 50; steps++) {
        if (current1 == current2) return current1;
        if (current1 == 0 || current2 == 0) break;
        if (current1 > max_taxon_id || current2 > max_taxon_id) break;
        
        current1 = parent_lookup[current1];
        current2 = parent_lookup[current2];
    }
    
    return 1;  // Root fallback
}

// Device function to calculate normalized phylogenetic distance
__device__ inline float calculate_phylogenetic_distance_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    uint32_t max_taxon_id) {
    
    if (taxon1 == taxon2) return 0.0f;
    if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return 1.0f;
    
    // Find LCA
    uint32_t lca = find_lca_gpu_phylo(taxon1, taxon2, parent_lookup, depth_lookup, max_taxon_id);
    
    // Get depths
    uint8_t depth1 = depth_lookup[taxon1];
    uint8_t depth2 = depth_lookup[taxon2];
    uint8_t lca_depth = (lca <= max_taxon_id) ? depth_lookup[lca] : 0;
    
    // Calculate total distance (sum of distances to LCA)
    int total_distance = (depth1 - lca_depth) + (depth2 - lca_depth);
    
    // Normalize by maximum reasonable distance (20 taxonomic levels)
    return fminf(1.0f, (float)total_distance / 20.0f);
}

// Device function to check if two taxa are closely related
__device__ inline bool are_taxa_closely_related_gpu(
    uint32_t taxon1,
    uint32_t taxon2,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    const uint8_t* rank_lookup,
    uint32_t max_taxon_id,
    uint8_t max_rank_distance = 3) {  // Allow up to 3 rank levels difference
    
    if (taxon1 == taxon2) return true;
    if (taxon1 > max_taxon_id || taxon2 > max_taxon_id) return false;
    
    // Find LCA
    uint32_t lca = find_lca_gpu_phylo(taxon1, taxon2, parent_lookup, depth_lookup, max_taxon_id);
    
    // Check if LCA is at a reasonable rank (not too distant)
    if (lca <= max_taxon_id) {
        uint8_t lca_rank = rank_lookup[lca];
        
        // If LCA is at genus level or below, consider closely related
        // Rank encoding: 7=genus, 8=species, etc.
        if (lca_rank >= 7) return true;
        
        // If LCA is at family level and both taxa are within max_rank_distance, still related
        if (lca_rank >= 6) {  // Family level
            uint8_t depth1 = depth_lookup[taxon1];
            uint8_t depth2 = depth_lookup[taxon2];
            uint8_t lca_depth = depth_lookup[lca];
            
            int max_distance = max((depth1 - lca_depth), (depth2 - lca_depth));
            return max_distance <= max_rank_distance;
        }
    }
    
    return false;
}

// Device function for weighted phylogenetic consistency calculation
__device__ inline float calculate_weighted_phylogenetic_consistency_gpu(
    const uint32_t* hit_taxa,
    const uint32_t* hit_votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    uint32_t max_taxon_id,
    float distance_decay_rate = 2.0f) {
    
    if (num_taxa <= 1) return 1.0f;
    
    float total_weighted_score = 0.0f;
    uint32_t total_votes = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (hit_votes[i] > 0) {
            // Calculate phylogenetic distance to primary taxon
            float distance = calculate_phylogenetic_distance_gpu(
                primary_taxon, hit_taxa[i], parent_lookup, depth_lookup, max_taxon_id
            );
            
            // Calculate exponential decay weight
            float weight = expf(-distance * distance_decay_rate);
            
            total_weighted_score += weight * hit_votes[i];
            total_votes += hit_votes[i];
        }
    }
    
    return total_votes > 0 ? total_weighted_score / total_votes : 0.0f;
}

// Device function for rank-based weighting
__device__ inline float calculate_rank_weighted_consistency_gpu(
    const uint32_t* hit_taxa,
    const uint32_t* hit_votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    const uint8_t* rank_lookup,
    uint32_t max_taxon_id) {
    
    if (num_taxa <= 1) return 1.0f;
    
    float total_weighted_score = 0.0f;
    uint32_t total_votes = 0;
    
    for (int i = 0; i < num_taxa; i++) {
        if (hit_votes[i] > 0) {
            // Find LCA and its rank
            uint32_t lca = find_lca_gpu_phylo(
                primary_taxon, hit_taxa[i], parent_lookup, depth_lookup, max_taxon_id
            );
            
            float weight = 1.0f;
            if (lca <= max_taxon_id) {
                uint8_t lca_rank = rank_lookup[lca];
                
                // Weight based on LCA rank level
                switch (lca_rank) {
                    case 8:  // species
                    case 9:  // subspecies
                        weight = 1.0f; break;
                    case 7:  // genus
                        weight = 0.9f; break;
                    case 6:  // family
                        weight = 0.7f; break;
                    case 5:  // order
                        weight = 0.5f; break;
                    case 4:  // class
                        weight = 0.3f; break;
                    case 3:  // phylum
                        weight = 0.1f; break;
                    default:
                        weight = 0.05f; break;
                }
            } else {
                weight = 0.01f;  // Very distant
            }
            
            total_weighted_score += weight * hit_votes[i];
            total_votes += hit_votes[i];
        }
    }
    
    return total_votes > 0 ? total_weighted_score / total_votes : 0.0f;
}

// Device function for hybrid phylogenetic consistency
__device__ inline float calculate_hybrid_phylogenetic_consistency_gpu(
    const uint32_t* hit_taxa,
    const uint32_t* hit_votes,
    int num_taxa,
    uint32_t primary_taxon,
    const uint32_t* parent_lookup,
    const uint8_t* depth_lookup,
    const uint8_t* rank_lookup,
    uint32_t max_taxon_id,
    float distance_weight = 0.7f,
    float rank_weight = 0.3f) {
    
    if (num_taxa <= 1) return 1.0f;
    
    // Calculate distance-based consistency
    float distance_consistency = calculate_weighted_phylogenetic_consistency_gpu(
        hit_taxa, hit_votes, num_taxa, primary_taxon,
        parent_lookup, depth_lookup, max_taxon_id
    );
    
    // Calculate rank-based consistency
    float rank_consistency = calculate_rank_weighted_consistency_gpu(
        hit_taxa, hit_votes, num_taxa, primary_taxon,
        parent_lookup, depth_lookup, rank_lookup, max_taxon_id
    );
    
    // Combine with weights
    return distance_weight * distance_consistency + rank_weight * rank_consistency;
}

// Host function to copy phylogenetic data to GPU
__host__ inline bool copy_phylogenetic_data_to_gpu(
    const std::vector<uint32_t>& parent_lookup,
    const std::vector<uint8_t>& depth_lookup,
    const std::vector<uint8_t>& rank_lookup,
    uint32_t max_taxon_id,
    GPUPhylogeneticData* d_phylo_data) {
    
    // Allocate GPU memory
    size_t parent_size = (max_taxon_id + 1) * sizeof(uint32_t);
    size_t depth_size = (max_taxon_id + 1) * sizeof(uint8_t);
    size_t rank_size = (max_taxon_id + 1) * sizeof(uint8_t);
    
    cudaError_t error;
    
    error = cudaMalloc(&d_phylo_data->parent_lookup, parent_size);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate parent lookup: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    error = cudaMalloc(&d_phylo_data->depth_lookup, depth_size);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate depth lookup: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    error = cudaMalloc(&d_phylo_data->rank_lookup, rank_size);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate rank lookup: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    // Copy data to GPU
    error = cudaMemcpy(d_phylo_data->parent_lookup, parent_lookup.data(), parent_size, cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        std::cerr << "Failed to copy parent lookup: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    error = cudaMemcpy(d_phylo_data->depth_lookup, depth_lookup.data(), depth_size, cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        std::cerr << "Failed to copy depth lookup: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    error = cudaMemcpy(d_phylo_data->rank_lookup, rank_lookup.data(), rank_size, cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        std::cerr << "Failed to copy rank lookup: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    d_phylo_data->max_taxon_id = max_taxon_id;
    
    std::cout << "Copied phylogenetic data to GPU: " << max_taxon_id + 1 << " taxa, "
              << (parent_size + depth_size + rank_size) / 1024 / 1024 << " MB" << std::endl;
    
    return true;
}

// Host function to free GPU phylogenetic data
__host__ inline void free_phylogenetic_data_gpu(GPUPhylogeneticData* d_phylo_data) {
    if (d_phylo_data->parent_lookup) {
        cudaFree(d_phylo_data->parent_lookup);
        d_phylo_data->parent_lookup = nullptr;
    }
    if (d_phylo_data->depth_lookup) {
        cudaFree(d_phylo_data->depth_lookup);
        d_phylo_data->depth_lookup = nullptr;
    }
    if (d_phylo_data->rank_lookup) {
        cudaFree(d_phylo_data->rank_lookup);
        d_phylo_data->rank_lookup = nullptr;
    }
    d_phylo_data->max_taxon_id = 0;
}

#endif // GPU_PHYLOGENETIC_CALCULATOR_CUH