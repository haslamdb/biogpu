// arg_assignment_gpu.cu
#ifndef ARG_ASSIGNMENT_GPU_CU
#define ARG_ASSIGNMENT_GPU_CU

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>
#include <cub/cub.cuh>

namespace cg = cooperative_groups;

// ARG assignment result structure
struct ARGAssignment {
    uint32_t read_id;
    uint32_t arg_id;
    uint32_t protein_id;
    int8_t frame;
    uint16_t alignment_start;
    uint16_t alignment_end;
    uint16_t alignment_length;
    float identity;
    float coverage;  // Fraction of ARG covered
    float confidence_score;
    bool is_complete_gene;  // Full-length match
    uint8_t num_mutations;
    char drug_class[32];
    char gene_name[64];
};

// Coverage statistics per ARG
struct ARGCoverage {
    uint32_t arg_id;
    uint32_t* position_counts;  // Coverage depth at each position
    uint32_t total_reads;
    uint32_t unique_reads;
    float mean_coverage;
    float coverage_uniformity;  // Evenness of coverage
    uint16_t covered_length;
    uint16_t gene_length;
};

// Enhanced protein match processing for ARG assignment
__global__ void process_arg_assignments_kernel(
    const ProteinMatch* protein_matches,
    const uint32_t* match_counts,
    const uint32_t* protein_to_arg_map,
    const ARGEntry* arg_database,
    const int num_reads,
    ARGAssignment* assignments,
    uint32_t* assignment_counts,
    ARGCoverage* coverage_stats,
    const float min_identity = 0.90f,
    const float min_coverage = 0.80f
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    uint32_t num_matches = match_counts[tid];
    if (num_matches == 0) {
        assignment_counts[tid] = 0;
        return;
    }
    
    uint32_t assignments_for_read = 0;
    
    // Process each protein match for this read
    for (uint32_t i = 0; i < num_matches && i < 32; i++) {
        const ProteinMatch& match = protein_matches[tid * 32 + i];
        
        // Check if this protein is an ARG
        uint32_t arg_id = protein_to_arg_map[match.protein_id];
        if (arg_id == UINT32_MAX) continue;  // Not an ARG
        
        const ARGEntry& arg = arg_database[arg_id];
        
        // Apply stringent filtering for ARG assignment
        if (match.identity < min_identity) continue;
        
        // Calculate coverage of the ARG
        float coverage = (float)match.match_length / arg.length;
        if (coverage < min_coverage) continue;
        
        // Create assignment
        ARGAssignment& assignment = assignments[tid * 10 + assignments_for_read];
        assignment.read_id = tid;
        assignment.arg_id = arg_id;
        assignment.protein_id = match.protein_id;
        assignment.frame = match.frame;
        assignment.alignment_start = match.ref_start;
        assignment.alignment_end = match.ref_start + match.match_length;
        assignment.alignment_length = match.match_length;
        assignment.identity = match.identity;
        assignment.coverage = coverage;
        assignment.num_mutations = match.num_mutations;
        
        // Calculate confidence score based on multiple factors
        float length_score = min(1.0f, (float)match.match_length / 50.0f);
        float identity_score = (match.identity - 0.8f) / 0.2f;  // Scale 0.8-1.0 to 0-1
        float coverage_score = coverage;
        assignment.confidence_score = (length_score + identity_score + coverage_score) / 3.0f;
        
        // Check if this is a complete gene match
        assignment.is_complete_gene = (coverage >= 0.95f && match.identity >= 0.95f);
        
        // Copy metadata (would be done differently in practice)
        for (int j = 0; j < 32; j++) {
            assignment.drug_class[j] = arg.drug_class[j];
        }
        for (int j = 0; j < 64; j++) {
            assignment.gene_name[j] = arg.gene_name[j];
        }
        
        // Update coverage statistics (atomic operations)
        if (coverage_stats != nullptr) {
            atomicAdd(&coverage_stats[arg_id].total_reads, 1);
            
            // Update position-specific coverage
            for (uint16_t pos = assignment.alignment_start; pos < assignment.alignment_end; pos++) {
                atomicAdd(&coverage_stats[arg_id].position_counts[pos], 1);
            }
        }
        
        assignments_for_read++;
        if (assignments_for_read >= 10) break;  // Max assignments per read
    }
    
    assignment_counts[tid] = assignments_for_read;
}

// Calculate coverage statistics
__global__ void calculate_coverage_stats_kernel(
    ARGCoverage* coverage_stats,
    const ARGEntry* arg_database,
    const int num_args
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_args) return;
    
    ARGCoverage& stats = coverage_stats[tid];
    const ARGEntry& arg = arg_database[tid];
    
    // Calculate coverage metrics
    uint32_t total_coverage = 0;
    uint32_t covered_positions = 0;
    float sum_squared_diff = 0.0f;
    
    for (uint16_t i = 0; i < arg.length; i++) {
        uint32_t depth = stats.position_counts[i];
        if (depth > 0) {
            covered_positions++;
            total_coverage += depth;
        }
    }
    
    stats.covered_length = covered_positions;
    stats.gene_length = arg.length;
    
    if (covered_positions > 0) {
        stats.mean_coverage = (float)total_coverage / arg.length;
        
        // Calculate uniformity (coefficient of variation)
        for (uint16_t i = 0; i < arg.length; i++) {
            float diff = stats.position_counts[i] - stats.mean_coverage;
            sum_squared_diff += diff * diff;
        }
        
        float variance = sum_squared_diff / arg.length;
        float std_dev = sqrtf(variance);
        stats.coverage_uniformity = 1.0f - (std_dev / (stats.mean_coverage + 1e-6f));
    } else {
        stats.mean_coverage = 0.0f;
        stats.coverage_uniformity = 0.0f;
    }
}

#endif // ARG_ASSIGNMENT_GPU_CU
