// enhanced_mutation_detection_unified.cu
// Updated mutation detection that uses the global FQ resistance mapper

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdint.h>
#include <cstdio>
#include "global_fq_resistance_mapper.h"

// Copy GPU structure definition from header
struct FQResistanceMutationGPU_Device {
    uint32_t species_gene_id;
    uint16_t position;
    char wildtype_aa;
    char mutant_aas[8];
    uint8_t num_mutants;
    float resistance_score;
};

// Device-side resistance checking
__device__ bool check_fq_resistance_gpu(
    const FQResistanceMutationGPU_Device* mutations,
    uint32_t num_mutations,
    uint32_t species_gene_id,
    uint16_t position,
    char observed_aa,
    char* expected_wildtype
) {
    // Binary search for the mutation position
    int left = 0;
    int right = num_mutations - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        const FQResistanceMutationGPU_Device& mut = mutations[mid];
        
        // Compare by species_gene_id first, then position
        if (mut.species_gene_id < species_gene_id) {
            left = mid + 1;
        } else if (mut.species_gene_id > species_gene_id) {
            right = mid - 1;
        } else {
            // Same species_gene_id, compare position
            if (mut.position < position) {
                left = mid + 1;
            } else if (mut.position > position) {
                right = mid - 1;
            } else {
                // Found the position
                *expected_wildtype = mut.wildtype_aa;
                
                // Check if observed AA is in the resistant list
                for (int i = 0; i < mut.num_mutants; i++) {
                    if (mut.mutant_aas[i] == observed_aa) {
                        return true;
                    }
                }
                return false;
            }
        }
    }
    
    // Position not found in resistance database
    *expected_wildtype = 'X';
    return false;
}

// Enhanced protein match structure
struct EnhancedProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint32_t species_gene_id;  // Combined ID for FQ mapper
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    
    // FQ resistance specific fields
    uint8_t num_fq_mutations;
    uint8_t num_total_mutations;
    uint16_t fq_positions[10];
    char wildtype_aas[10];
    char observed_aas[10];
    bool is_fq_resistant[10];
    
    // Alignment quality
    bool covers_qrdr;
    char query_sequence[100];
    char ref_sequence[100];
};

// Main mutation detection kernel using global mapper
__global__ void detect_fq_mutations_kernel(
    const EnhancedProteinMatch* protein_matches,
    const uint32_t* match_counts,
    const FQResistanceMutationGPU_Device* fq_mutations,
    const uint32_t num_fq_mutations,
    const int num_reads,
    EnhancedProteinMatch* output_matches,
    uint32_t* output_counts
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    uint32_t num_matches = match_counts[tid];
    if (num_matches == 0) {
        output_counts[tid] = 0;
        return;
    }
    
    uint32_t output_count = 0;
    
    for (uint32_t i = 0; i < num_matches && output_count < 32; i++) {
        const EnhancedProteinMatch& match = protein_matches[tid * 32 + i];
        
        // Copy basic match info
        EnhancedProteinMatch enhanced = match;
        enhanced.num_fq_mutations = 0;
        enhanced.num_total_mutations = 0;
        
        // Check each mutation position against FQ database
        for (int j = 0; j < match.num_total_mutations && j < 10; j++) {
            uint16_t global_position = match.ref_start + match.fq_positions[j] + 1; // Convert to 1-based
            char observed_aa = match.observed_aas[j];
            char expected_wildtype = 'X';
            
            bool is_fq = check_fq_resistance_gpu(
                fq_mutations, num_fq_mutations,
                match.species_gene_id, global_position,
                observed_aa, &expected_wildtype
            );
            
            if (is_fq || expected_wildtype != 'X') {
                // This is a known FQ resistance position
                enhanced.fq_positions[enhanced.num_fq_mutations] = global_position;
                enhanced.wildtype_aas[enhanced.num_fq_mutations] = expected_wildtype;
                enhanced.observed_aas[enhanced.num_fq_mutations] = observed_aa;
                enhanced.is_fq_resistant[enhanced.num_fq_mutations] = is_fq;
                enhanced.num_fq_mutations++;
            }
        }
        
        // Check QRDR coverage based on gene and position
        enhanced.covers_qrdr = false;
        if (match.species_gene_id != UINT32_MAX) {
            // Simple QRDR check - you can make this more sophisticated
            uint16_t start = match.ref_start + 1;
            uint16_t end = start + match.match_length;
            
            // Check common QRDR positions (this should be loaded from database too)
            if ((start <= 83 && end >= 83) || (start <= 87 && end >= 87)) {
                enhanced.covers_qrdr = true;
            }
        }
        
        output_matches[tid * 32 + output_count] = enhanced;
        output_count++;
    }
    
    output_counts[tid] = output_count;
}

// Host wrapper function
extern "C" {
    int detect_fq_mutations_unified(
        const void* d_protein_matches,
        const uint32_t* d_match_counts,
        int num_reads,
        void* d_output_matches,
        uint32_t* d_output_counts
    ) {
        // Get GPU mutations from global mapper
        void* gpu_mutations = get_fq_gpu_mutations();
        uint32_t num_mutations = get_fq_gpu_mutation_count();
        
        if (!gpu_mutations || num_mutations == 0) {
            printf("ERROR: FQ mutation database not loaded on GPU\n");
            return -1;
        }
        
        // Launch kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        detect_fq_mutations_kernel<<<grid_size, block_size>>>(
            static_cast<const EnhancedProteinMatch*>(d_protein_matches),
            d_match_counts,
            static_cast<const FQResistanceMutationGPU_Device*>(gpu_mutations),
            num_mutations,
            num_reads,
            static_cast<EnhancedProteinMatch*>(d_output_matches),
            d_output_counts
        );
        
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            printf("ERROR: CUDA kernel failed: %s\n", cudaGetErrorString(error));
            return -1;
        }
        
        cudaDeviceSynchronize();
        return 0;
    }
    
    // Host-side function to check FQ mutations
    void check_fq_mutations_host(
        const char* species,
        const char* gene,
        const EnhancedProteinMatch* match,
        bool* fq_results,
        int* num_fq_mutations
    ) {
        *num_fq_mutations = 0;
        
        for (int i = 0; i < match->num_total_mutations && i < 10; i++) {
            uint16_t position = match->ref_start + match->fq_positions[i] + 1; // 1-based
            char observed = match->observed_aas[i];
            
            bool is_fq = is_fq_resistance_mutation(species, gene, position, 
                                                  match->wildtype_aas[i], observed);
            
            if (is_fq) {
                fq_results[*num_fq_mutations] = true;
                (*num_fq_mutations)++;
            }
        }
    }
}