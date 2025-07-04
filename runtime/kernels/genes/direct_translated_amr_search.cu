// direct_translated_amr_search.cu
#ifndef DIRECT_TRANSLATED_AMR_SEARCH_CU
#define DIRECT_TRANSLATED_AMR_SEARCH_CU

#include <cuda_runtime.h>
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

// Genetic code lookup table (stored in constant memory for fast access)
__constant__ char CODON_TABLE[64] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
};

// Nucleotide to index mapping
__device__ __forceinline__ int nucToIndex(char nuc) {
    switch(nuc) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1; // N or invalid
    }
}

// Fast 6-frame translation kernel
__global__ void translateReadsKernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    char* translated_reads,  // Output: 6 frames per read
    int* frame_offsets,      // Where each frame starts
    int* frame_lengths,      // Length of each translated frame
    const int num_reads,
    const int max_read_length
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int reads_per_thread = (num_reads + gridDim.x * blockDim.x - 1) / (gridDim.x * blockDim.x);
    
    // Process multiple reads per thread for better utilization
    for (int r = 0; r < reads_per_thread; r++) {
        const int read_idx = tid * reads_per_thread + r;
        if (read_idx >= num_reads) return;
        
        const int read_start = read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        const char* read = &reads[read_start];
        
        // Calculate output position for this read's 6 frames
        const int output_base = read_idx * max_read_length * 6;
        
        // Translate all 6 frames
        for (int frame = 0; frame < 6; frame++) {
            const bool is_reverse = (frame >= 3);
            const int frame_shift = frame % 3;
            const int frame_output_start = output_base + frame * max_read_length;
            
            frame_offsets[read_idx * 6 + frame] = frame_output_start;
            
            int aa_pos = 0;
            
            // Forward frames (0, 1, 2) or reverse frames (3, 4, 5)
            if (!is_reverse) {
                // Forward translation
                for (int pos = frame_shift; pos + 2 < read_len; pos += 3) {
                    int codon_idx = 0;
                    bool valid_codon = true;
                    
                    // Build codon index
                    for (int i = 0; i < 3; i++) {
                        int nuc_idx = nucToIndex(read[pos + i]);
                        if (nuc_idx < 0) {
                            valid_codon = false;
                            break;
                        }
                        codon_idx = (codon_idx << 2) | nuc_idx;
                    }
                    
                    if (valid_codon) {
                        translated_reads[frame_output_start + aa_pos] = CODON_TABLE[codon_idx];
                    } else {
                        translated_reads[frame_output_start + aa_pos] = 'X'; // Unknown amino acid
                    }
                    aa_pos++;
                }
            } else {
                // Reverse complement translation
                for (int pos = read_len - 1 - frame_shift; pos >= 2; pos -= 3) {
                    int codon_idx = 0;
                    bool valid_codon = true;
                    
                    // Build codon index from reverse complement
                    for (int i = 0; i < 3; i++) {
                        char nuc = read[pos - i];
                        // Complement
                        switch(nuc) {
                            case 'A': case 'a': nuc = 'T'; break;
                            case 'T': case 't': nuc = 'A'; break;
                            case 'C': case 'c': nuc = 'G'; break;
                            case 'G': case 'g': nuc = 'C'; break;
                        }
                        
                        int nuc_idx = nucToIndex(nuc);
                        if (nuc_idx < 0) {
                            valid_codon = false;
                            break;
                        }
                        codon_idx = (codon_idx << 2) | nuc_idx;
                    }
                    
                    if (valid_codon) {
                        translated_reads[frame_output_start + aa_pos] = CODON_TABLE[codon_idx];
                    } else {
                        translated_reads[frame_output_start + aa_pos] = 'X';
                    }
                    aa_pos++;
                }
            }
            
            frame_lengths[read_idx * 6 + frame] = aa_pos;
        }
    }
}

// Direct protein search kernel using sliding window
__global__ void searchAMRProteinsKernel(
    const char* translated_reads,     // 6-frame translated reads
    const int* frame_offsets,
    const int* frame_lengths,
    const char* amr_proteins,         // AMR protein database
    const int* protein_offsets,
    const int* protein_lengths,
    const float* identity_thresholds, // Per-protein thresholds
    AMRMatch* matches,                // Output matches
    int* match_counts,                // Matches per read
    const int num_reads,
    const int num_proteins,
    const int min_match_length = 20,  // Minimum amino acid match
    const float default_identity = 0.85f
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int read_idx = tid / 6;  // Each thread handles one frame
    const int frame = tid % 6;
    
    if (read_idx >= num_reads) return;
    
    // Get this frame's translated sequence
    const int frame_start = frame_offsets[read_idx * 6 + frame];
    const int frame_len = frame_lengths[read_idx * 6 + frame];
    const char* frame_seq = &translated_reads[frame_start];
    
    // Shared memory for protein sequences (cache frequently accessed proteins)
    extern __shared__ char shared_mem[];
    char* cached_proteins = shared_mem;
    const int cache_size = 4096; // Adjust based on shared memory availability
    
    int matches_found = 0;
    const int max_matches_per_read = 10;
    
    // Search against each AMR protein
    for (int prot_idx = 0; prot_idx < num_proteins; prot_idx++) {
        const int prot_start = protein_offsets[prot_idx];
        const int prot_len = protein_lengths[prot_idx];
        const char* protein = &amr_proteins[prot_start];
        const float threshold = identity_thresholds ? identity_thresholds[prot_idx] : default_identity;
        
        // Skip if protein is longer than frame
        if (prot_len > frame_len) continue;
        
        // Sliding window search
        for (int start_pos = 0; start_pos <= frame_len - min_match_length; start_pos++) {
            // Quick check: look for exact matches in critical positions
            // This is a heuristic - you might want to check specific conserved residues
            bool potential_match = true;
            
            // Check first few amino acids for quick rejection
            for (int i = 0; i < min(5, prot_len); i++) {
                if (frame_seq[start_pos + i] != protein[i] && 
                    frame_seq[start_pos + i] != 'X' && 
                    protein[i] != 'X') {
                    potential_match = false;
                    break;
                }
            }
            
            if (!potential_match) continue;
            
            // Perform detailed alignment
            int matches = 0;
            int mismatches = 0;
            int max_extension = min(prot_len, frame_len - start_pos);
            
            for (int i = 0; i < max_extension; i++) {
                if (frame_seq[start_pos + i] == protein[i]) {
                    matches++;
                } else if (frame_seq[start_pos + i] != 'X' && protein[i] != 'X') {
                    mismatches++;
                    // Early termination if too many mismatches
                    if (mismatches > (1.0f - threshold) * max_extension) break;
                }
            }
            
            int aligned_length = matches + mismatches;
            if (aligned_length < min_match_length) continue;
            
            float identity = (float)matches / aligned_length;
            
            if (identity >= threshold) {
                // Found a match!
                if (matches_found < max_matches_per_read) {
                    int match_idx = read_idx * max_matches_per_read + matches_found;
                    AMRMatch& match = matches[match_idx];
                    
                    match.read_id = read_idx;
                    match.protein_id = prot_idx;
                    match.frame = frame;
                    match.read_start = start_pos;
                    match.protein_start = 0;
                    match.match_length = aligned_length;
                    match.identity = identity;
                    match.num_matches = matches;
                    match.num_mismatches = mismatches;
                    match.is_complete_protein = (aligned_length >= prot_len * 0.95f);
                    
                    matches_found++;
                }
                
                // Skip ahead to avoid overlapping matches
                start_pos += aligned_length / 2;
            }
        }
    }
    
    // Store match count for this read
    if (frame == 0) {
        match_counts[read_idx] = 0;
    }
    __syncthreads();
    
    // Accumulate matches from all frames
    atomicAdd(&match_counts[read_idx], matches_found);
}

// Optimized version using protein profiles for faster search
__global__ void searchAMRProfilesKernel(
    const char* translated_reads,
    const int* frame_offsets,
    const int* frame_lengths,
    const ProteinProfile* amr_profiles,  // Pre-computed profiles with conserved positions
    const int num_reads,
    const int num_profiles,
    AMRProfileMatch* matches,
    int* match_counts
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int read_idx = tid / 6;
    const int frame = tid % 6;
    
    if (read_idx >= num_reads) return;
    
    const int frame_start = frame_offsets[read_idx * 6 + frame];
    const int frame_len = frame_lengths[read_idx * 6 + frame];
    const char* frame_seq = &translated_reads[frame_start];
    
    // Use shared memory for protein profiles
    extern __shared__ ProteinProfile shared_profiles[];
    
    // Cooperative loading of profiles
    auto block = cg::this_thread_block();
    for (int i = threadIdx.x; i < min(num_profiles, 32); i += blockDim.x) {
        shared_profiles[i] = amr_profiles[i];
    }
    block.sync();
    
    int matches_found = 0;
    
    // Search using profile matching (faster for large databases)
    for (int prof_idx = 0; prof_idx < num_profiles; prof_idx++) {
        const ProteinProfile& profile = (prof_idx < 32) ? 
            shared_profiles[prof_idx] : amr_profiles[prof_idx];
        
        // Check conserved positions first
        bool conserved_match = true;
        for (int i = 0; i < profile.num_conserved && i < frame_len; i++) {
            int pos = profile.conserved_positions[i];
            if (pos < frame_len && 
                frame_seq[pos] != profile.conserved_residues[i]) {
                conserved_match = false;
                break;
            }
        }
        
        if (!conserved_match) continue;
        
        // Detailed scoring if conserved positions match
        float score = scoreProteinMatch(frame_seq, frame_len, profile);
        
        if (score >= profile.score_threshold) {
            // Record match
            if (matches_found < 10) {
                int match_idx = read_idx * 10 + matches_found;
                matches[match_idx].read_id = read_idx;
                matches[match_idx].profile_id = prof_idx;
                matches[match_idx].frame = frame;
                matches[match_idx].score = score;
                matches_found++;
            }
        }
    }
    
    if (frame == 0) match_counts[read_idx] = 0;
    __syncthreads();
    atomicAdd(&match_counts[read_idx], matches_found);
}

// Helper structures
struct AMRMatch {
    uint32_t read_id;
    uint32_t protein_id;
    int8_t frame;
    uint16_t read_start;
    uint16_t protein_start;
    uint16_t match_length;
    float identity;
    uint16_t num_matches;
    uint16_t num_mismatches;
    bool is_complete_protein;
};

struct ProteinProfile {
    uint32_t protein_id;
    uint16_t length;
    uint8_t num_conserved;
    uint16_t conserved_positions[20];  // Key positions that must match
    char conserved_residues[20];       // Expected residues at those positions
    float score_threshold;
    char drug_class[32];
};

struct AMRProfileMatch {
    uint32_t read_id;
    uint32_t profile_id;
    int8_t frame;
    float score;
};

// Device function for profile scoring
__device__ float scoreProteinMatch(
    const char* sequence,
    int seq_length,
    const ProteinProfile& profile
) {
    // Implement position-specific scoring matrix (PSSM) or similar
    // This is a simplified version
    float score = 0.0f;
    int positions_checked = 0;
    
    for (int i = 0; i < profile.num_conserved; i++) {
        if (profile.conserved_positions[i] < seq_length) {
            if (sequence[profile.conserved_positions[i]] == profile.conserved_residues[i]) {
                score += 2.0f;  // Conserved position match
            } else {
                score -= 1.0f;  // Conserved position mismatch
            }
            positions_checked++;
        }
    }
    
    return positions_checked > 0 ? score / positions_checked : 0.0f;
}

#endif // DIRECT_TRANSLATED_AMR_SEARCH_CU
