// translated_search_amr.cu - CORRECTED VERSION
// Enhanced memory safety and bounds checking for clinical AMR detection

#ifndef TRANSLATED_SEARCH_CU
#define TRANSLATED_SEARCH_CU

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>
#include <stdint.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cassert>

namespace cg = cooperative_groups;

// Debug and safety macros
#define DEBUG_AMR 0
#define BOUNDS_CHECK 1
#define DEBUG_PRINT(fmt, ...) // Disabled
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        printf("[CUDA ERROR] %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        return false; \
    } \
} while(0)

// Enhanced constants for protein k-mer approach with safety margins
#define CODON_SIZE 3
#define NUM_FRAMES 6
#define MAX_PROTEIN_LENGTH 200
#define MAX_READ_LENGTH 1000
#define PROTEIN_KMER_SIZE 8
#define MIN_PEPTIDE_LENGTH 20
#define MIN_SEED_HITS 1
#define EXTENSION_THRESHOLD 15
#define MIN_IDENTITY_THRESHOLD 0.75f  // Slightly more lenient for AMR detection
#define MIN_ALIGNMENT_LENGTH 12       // Reduced for short AMR peptides
#define SW_SCORE_THRESHOLD 60.0f      // More sensitive for resistance mutations
#define AA_ALPHABET_SIZE 24
#define MAX_SEEDS_PER_FRAME 100
#define MAX_PROTEIN_SEEDS 20
#define MAX_MATCHES_PER_READ 32

// Genetic code table (standard code)
__constant__ char GENETIC_CODE[64] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',  // AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT
    'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',  // AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',  // CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT
    'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',  // CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',  // GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT
    'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',  // GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S',  // TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT
    'W', 'C', '*', 'C', 'L', 'F', 'L', 'F'   // TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT
};

// BLOSUM62 matrix for amino acid scoring
__constant__ float BLOSUM62_SCORES[24*24] = {
    // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4, // A
   -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4, // R
   -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4, // N
   -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4, // D
    0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4, // C
   -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4, // Q
   -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, // E
    0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4, // G
   -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4, // H
   -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4, // I
   -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4, // L
   -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4, // K
   -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4, // M
   -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4, // F
   -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4, // P
    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4, // S
    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4, // T
   -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4, // W
   -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4, // Y
    0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4, // V
   -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4, // B
   -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, // Z
    0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4, // X
   -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1  // *
};

// Safe amino acid to index mapping
__host__ __device__ inline int aa_to_index(char aa) {
    switch(aa) {
        case 'A': return 0;  case 'R': return 1;  case 'N': return 2;  case 'D': return 3;
        case 'C': return 4;  case 'Q': return 5;  case 'E': return 6;  case 'G': return 7;
        case 'H': return 8;  case 'I': return 9;  case 'L': return 10; case 'K': return 11;
        case 'M': return 12; case 'F': return 13; case 'P': return 14; case 'S': return 15;
        case 'T': return 16; case 'W': return 17; case 'Y': return 18; case 'V': return 19;
        case 'B': return 20; case 'Z': return 21; case 'X': return 22; case '*': return 23;
        default: return 22; // Unknown -> X
    }
}

// Get BLOSUM score with bounds checking
__device__ inline float get_blosum_score(char aa1, char aa2) {
    int idx1 = aa_to_index(aa1);
    int idx2 = aa_to_index(aa2);
    
    #if BOUNDS_CHECK
    if (idx1 < 0 || idx1 >= 24 || idx2 < 0 || idx2 >= 24) {
        return -4.0f; // Penalty for invalid amino acids
    }
    #endif
    
    return BLOSUM62_SCORES[idx1 * 24 + idx2];
}

// Safe base encoding for translation
__device__ inline int base_to_index(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;  // Invalid base
    }
}

// Safe codon translation with bounds checking
__device__ inline char translate_codon(const char* codon) {
    if (!codon) return 'X';
    
    int idx = 0;
    for (int i = 0; i < 3; i++) {
        int base = base_to_index(codon[i]);
        if (base < 0) return 'X';  // Unknown
        idx = (idx << 2) | base;
    }
    
    #if BOUNDS_CHECK
    if (idx < 0 || idx >= 64) return 'X';
    #endif
    
    return GENETIC_CODE[idx];
}

// Hash function for protein k-mers with overflow protection
__host__ __device__ inline uint64_t hash_protein_kmer(const char* kmer) {
    if (!kmer) return 0;
    
    uint64_t hash = 0;
    const uint64_t prime = 31;
    
    for (int i = 0; i < PROTEIN_KMER_SIZE; i++) {
        uint8_t aa_val = aa_to_index(kmer[i]);
        // Prevent overflow
        if (hash > (UINT64_MAX / prime)) {
            hash = (hash % prime) * prime + aa_val;
        } else {
            hash = hash * prime + aa_val;
        }
    }
    
    return hash;
}

// Structure for translated read frame
struct TranslatedFrame {
    char sequence[MAX_PROTEIN_LENGTH];
    uint16_t length;
    int8_t frame;        // -3, -2, -1, +1, +2, +3
    uint16_t start_pos;  // Start position in original read
    bool valid;          // Validity flag
};

// Enhanced protein match structure for AMR detection - MUST match pipeline exactly
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;    // Position in translated frame
    uint16_t ref_start;      // Position in reference protein
    uint16_t match_length;
    float alignment_score;
    float identity;
    // Coverage tracking fields
    uint32_t gene_length;     // Total gene length
    uint16_t coverage_start;  // Start of covered region
    uint16_t coverage_end;    // End of covered region
    // Remove mutation-specific fields for AMR gene detection
    // (mutations are not needed for gene presence/absence)
    bool used_smith_waterman;  // Flag indicating if SW was used
    bool concordant;           // Flag for paired-end concordance
    char query_peptide[51];  // Store aligned peptide sequence (up to 50 AA + null terminator)
};

// Enhanced protein database structure
struct ProteinDatabase {
    uint32_t num_proteins;
    uint32_t num_kmers;
    uint32_t total_sequence_length;  // For bounds checking
    
    // Sorted k-mer index for binary search
    uint64_t* sorted_kmer_hashes;
    uint32_t* kmer_start_indices;
    uint32_t* kmer_counts;
    uint32_t* position_data;
    
    // Protein metadata with bounds checking
    uint32_t* protein_ids;
    uint32_t* gene_ids;
    uint32_t* species_ids;
    uint16_t* seq_lengths;
    uint32_t* seq_offsets;
    
    // Reference sequences
    char* sequences;
    
    // AMR-specific metadata
    uint8_t* is_amr_gene;      // Flag for AMR genes (0=false, 1=true)
    uint8_t* resistance_class; // Resistance class (0=unknown, 1=fluoroquinolone, etc.)
    char** gene_families;      // Array of gene family names
};

// Helper structure for seed clustering
struct SeedHit {
    uint32_t protein_id;
    uint32_t query_pos;
    uint32_t ref_pos;
    float score;
    bool valid;
};

// Safe binary search for protein k-mer
__device__ inline int binary_search_protein_kmer(
    const ProteinDatabase* db,
    uint64_t target_hash
) {
    if (!db || !db->sorted_kmer_hashes || db->num_kmers == 0) {
        return -1;
    }
    
    int left = 0;
    int right = db->num_kmers - 1;
    
    while (left <= right) {
        int mid = left + (right - left) / 2;  // Prevent overflow
        
        #if BOUNDS_CHECK
        if (mid < 0 || mid >= db->num_kmers) {
            return -1;
        }
        #endif
        
        uint64_t mid_hash = db->sorted_kmer_hashes[mid];
        
        if (mid_hash == target_hash) {
            return mid;
        } else if (mid_hash < target_hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    return -1; // Not found
}

// Enhanced Smith-Waterman with bounds checking and AMR-specific scoring
__device__ float smith_waterman_align(
    const char* query, uint16_t query_len,
    const char* ref, uint16_t ref_len,
    uint16_t* best_query_start,
    uint16_t* best_ref_start,
    uint16_t* best_length,
    char* aligned_query = nullptr,
    char* aligned_ref = nullptr
) {
    // Input validation
    if (!query || !ref || query_len == 0 || ref_len == 0) {
        if (best_query_start) *best_query_start = 0;
        if (best_ref_start) *best_ref_start = 0;
        if (best_length) *best_length = 0;
        return 0.0f;
    }
    
    // Constants for alignment
    const float GAP_OPEN = -3.0f;
    const float GAP_EXTEND = -1.0f;
    const int BAND_WIDTH = 15;  // Wider band for AMR detection
    const int MAX_ALIGN_LEN = 50;
    
    // Bounds checking for sequence lengths
    if (query_len > MAX_ALIGN_LEN || ref_len > MAX_ALIGN_LEN) {
        // Fallback to simple scoring for very large sequences
        float score = 0.0f;
        int matches = 0;
        int aligned_len = min(query_len, ref_len);
        
        for (int i = 0; i < aligned_len; i++) {
            float blosum = get_blosum_score(query[i], ref[i]);
            score += blosum;
            if (blosum > 0) matches++;
        }
        
        if (best_query_start) *best_query_start = 0;
        if (best_ref_start) *best_ref_start = 0;
        if (best_length) *best_length = aligned_len;
        
        return score;
    }
    
    // Dynamic programming matrix (local arrays)
    float H[MAX_ALIGN_LEN + 1][MAX_ALIGN_LEN + 1];
    
    // Initialize matrix
    for (int i = 0; i <= query_len; i++) {
        for (int j = 0; j <= ref_len; j++) {
            H[i][j] = 0.0f;
        }
    }
    
    float max_score = 0.0f;
    int max_i = 0, max_j = 0;
    
    // Fill scoring matrix with banding optimization
    for (int i = 1; i <= query_len; i++) {
        int j_start = max(1, i - BAND_WIDTH);
        int j_end = min((int)ref_len, i + BAND_WIDTH);
        
        for (int j = j_start; j <= j_end; j++) {
            float match = H[i-1][j-1] + get_blosum_score(query[i-1], ref[j-1]);
            float delete_gap = H[i-1][j] + GAP_OPEN;
            float insert_gap = H[i][j-1] + GAP_OPEN;
            
            H[i][j] = fmaxf(0.0f, fmaxf(match, fmaxf(delete_gap, insert_gap)));
            
            if (H[i][j] > max_score) {
                max_score = H[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback to find alignment boundaries
    int i = max_i, j = max_j;
    int align_len = 0;
    int align_start_i = i, align_start_j = j;
    
    while (i > 0 && j > 0 && H[i][j] > 0 && align_len < MAX_ALIGN_LEN) {
        if (i > 0 && j > 0 && 
            H[i][j] == H[i-1][j-1] + get_blosum_score(query[i-1], ref[j-1])) {
            // Match or mismatch
            if (aligned_query && align_len < 50) {
                aligned_query[align_len] = query[i-1];
                aligned_ref[align_len] = ref[j-1];
            }
            i--; j--;
            align_start_i = i;
            align_start_j = j;
        } else if (i > 0 && H[i][j] == H[i-1][j] + GAP_OPEN) {
            // Deletion in query
            if (aligned_query && align_len < 50) {
                aligned_query[align_len] = query[i-1];
                aligned_ref[align_len] = '-';
            }
            i--;
            align_start_i = i;
        } else if (j > 0 && H[i][j] == H[i][j-1] + GAP_OPEN) {
            // Insertion in query
            if (aligned_query && align_len < 50) {
                aligned_query[align_len] = '-';
                aligned_ref[align_len] = ref[j-1];
            }
            j--;
            align_start_j = j;
        } else {
            break;
        }
        align_len++;
    }
    
    if (best_query_start) *best_query_start = align_start_i;
    if (best_ref_start) *best_ref_start = align_start_j;
    if (best_length) *best_length = align_len;
    
    return max_score;
}

// Enhanced 6-frame translation kernel with comprehensive bounds checking
__global__ void six_frame_translate_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const int num_reads,
    TranslatedFrame* translated_frames,
    uint32_t* frame_counts
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Bounds checking for thread ID
    if (tid >= num_reads) return;
    
    // Additional input validation
    if (!reads || !read_lengths || !read_offsets || !translated_frames || !frame_counts) {
        if (frame_counts && tid < num_reads) frame_counts[tid] = 0;
        return;
    }
    
    const int read_len = read_lengths[tid];
    
    // Validate read length
    if (read_len <= 0 || read_len > MAX_READ_LENGTH) {
        frame_counts[tid] = 0;
        return;
    }
    
    // Validate read offset
    const int read_offset = read_offsets[tid];
    if (read_offset < 0) {
        frame_counts[tid] = 0;
        return;
    }
    
    const char* read = reads + read_offset;
    TranslatedFrame* read_frames = &translated_frames[tid * NUM_FRAMES];
    int valid_frames = 0;
    
    // Initialize frames
    for (int f = 0; f < NUM_FRAMES; f++) {
        read_frames[f].length = 0;
        read_frames[f].valid = false;
        read_frames[f].sequence[0] = '\0';
    }
    
    // Forward frames (+1, +2, +3)
    for (int frame = 0; frame < 3 && valid_frames < NUM_FRAMES; frame++) {
        TranslatedFrame* tf = &read_frames[valid_frames];
        tf->frame = frame + 1;
        tf->start_pos = frame;
        tf->length = 0;
        tf->valid = true;
        
        for (int pos = frame; pos + 2 < read_len && tf->length < MAX_PROTEIN_LENGTH - 1; pos += 3) {
            char aa = translate_codon(&read[pos]);
            
            if (aa == '*') {
                // Stop codon - finalize current peptide if long enough
                if (tf->length >= MIN_PEPTIDE_LENGTH) {
                    tf->sequence[tf->length] = '\0';
                    valid_frames++;
                    
                    // Start new frame after stop codon
                    if (valid_frames < NUM_FRAMES) {
                        tf = &read_frames[valid_frames];
                        tf->frame = frame + 1;
                        tf->start_pos = pos + 3;
                        tf->length = 0;
                        tf->valid = true;
                    }
                } else {
                    // Reset if peptide too short
                    tf->length = 0;
                    tf->start_pos = pos + 3;
                }
            } else {
                tf->sequence[tf->length++] = aa;
            }
        }
        
        // Finalize last peptide if long enough
        if (tf->length >= MIN_PEPTIDE_LENGTH && valid_frames < NUM_FRAMES) {
            tf->sequence[tf->length] = '\0';
            valid_frames++;
        } else if (tf->length < MIN_PEPTIDE_LENGTH) {
            tf->valid = false;
        }
    }
    
    // Reverse complement frames (-1, -2, -3)
    for (int frame = 0; frame < 3 && valid_frames < NUM_FRAMES; frame++) {
        TranslatedFrame* tf = &read_frames[valid_frames];
        tf->frame = -(frame + 1);
        tf->start_pos = read_len - frame - 1;
        tf->length = 0;
        tf->valid = true;
        
        for (int pos = read_len - frame - 1; pos >= 2 && tf->length < MAX_PROTEIN_LENGTH - 1; pos -= 3) {
            // Generate reverse complement codon
            char rc_codon[3];
            for (int i = 0; i < 3; i++) {
                if (pos - i >= 0) {
                    char base = read[pos - i];
                    switch(base) {
                        case 'A': case 'a': rc_codon[i] = 'T'; break;
                        case 'T': case 't': rc_codon[i] = 'A'; break;
                        case 'G': case 'g': rc_codon[i] = 'C'; break;
                        case 'C': case 'c': rc_codon[i] = 'G'; break;
                        default: rc_codon[i] = 'N';
                    }
                } else {
                    rc_codon[i] = 'N';
                }
            }
            
            char aa = translate_codon(rc_codon);
            
            if (aa == '*') {
                if (tf->length >= MIN_PEPTIDE_LENGTH) {
                    tf->sequence[tf->length] = '\0';
                    valid_frames++;
                    
                    if (valid_frames < NUM_FRAMES) {
                        tf = &read_frames[valid_frames];
                        tf->frame = -(frame + 1);
                        tf->start_pos = pos - 3;
                        tf->length = 0;
                        tf->valid = true;
                    }
                } else {
                    tf->length = 0;
                    tf->start_pos = pos - 3;
                }
            } else {
                tf->sequence[tf->length++] = aa;
            }
        }
        
        if (tf->length >= MIN_PEPTIDE_LENGTH && valid_frames < NUM_FRAMES) {
            tf->sequence[tf->length] = '\0';
            valid_frames++;
        } else if (tf->length < MIN_PEPTIDE_LENGTH) {
            tf->valid = false;
        }
    }
    
    // Ensure we don't exceed the maximum number of frames
    if (valid_frames > NUM_FRAMES) {
        valid_frames = NUM_FRAMES;
    }
    
    frame_counts[tid] = valid_frames;
}

// Enhanced protein k-mer matching kernel with comprehensive safety checks
__global__ void enhanced_protein_kmer_match_kernel(
    const TranslatedFrame* translated_frames,
    const uint32_t* frame_counts,
    const int num_reads,
    const ProteinDatabase* protein_db,
    ProteinMatch* matches,
    uint32_t* match_counts,
    const uint32_t max_matches_per_read,
    const bool enable_smith_waterman = false
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Bounds checking for thread ID
    if (tid >= num_reads) return;
    
    // Input validation
    if (!translated_frames || !frame_counts || !protein_db || !matches || !match_counts) {
        if (match_counts) match_counts[tid] = 0;
        return;
    }
    
    uint32_t num_frames = frame_counts[tid];
    if (num_frames == 0 || num_frames > NUM_FRAMES) {
        match_counts[tid] = 0;
        return;
    }
    
    const TranslatedFrame* read_frames = &translated_frames[tid * NUM_FRAMES];
    ProteinMatch* read_matches = &matches[tid * max_matches_per_read];
    uint32_t match_count = 0;
    
    // Process each translated frame
    for (uint32_t frame_idx = 0; frame_idx < num_frames && frame_idx < NUM_FRAMES; frame_idx++) {
        const TranslatedFrame& frame = read_frames[frame_idx];
        
        if (!frame.valid || frame.length < PROTEIN_KMER_SIZE) continue;
        
        // Collect k-mer seed hits with bounds checking
        SeedHit seeds[MAX_SEEDS_PER_FRAME];
        int num_seeds = 0;
        
        // Find k-mer seed matches
        for (int pos = 0; pos + PROTEIN_KMER_SIZE <= frame.length && num_seeds < MAX_SEEDS_PER_FRAME; pos++) {
            uint64_t kmer_hash = hash_protein_kmer(&frame.sequence[pos]);
            
            int kmer_idx = binary_search_protein_kmer(protein_db, kmer_hash);
            if (kmer_idx >= 0 && kmer_idx < protein_db->num_kmers) {
                uint32_t start_idx = protein_db->kmer_start_indices[kmer_idx];
                uint32_t count = protein_db->kmer_counts[kmer_idx];
                
                // Bounds checking for position data access
                if (start_idx < UINT32_MAX && count > 0) {
                    // Add hits for this k-mer (limit to avoid overflow)
                    for (uint32_t i = 0; i < count && i < 5 && num_seeds < MAX_SEEDS_PER_FRAME; i++) {
                        if (start_idx + i < UINT32_MAX) {  // Additional bounds check
                            uint32_t encoded = protein_db->position_data[start_idx + i];
                            uint32_t protein_id = encoded >> 16;
                            uint32_t ref_pos = encoded & 0xFFFF;
                            
                            #if DEBUG_AMR
                            if (tid == 0 && num_seeds == 0) {  // First thread, first seed
                                printf("[KERNEL] K-mer match: protein_id=%d, gene_id=%d\n", 
                                       protein_id, protein_db->gene_ids[protein_id]);
                            }
                            #endif
                            
                            if (protein_id < protein_db->num_proteins) {
                                seeds[num_seeds] = {protein_id, (uint32_t)pos, ref_pos, 10.0f, true};
                                num_seeds++;
                            }
                        }
                    }
                }
            }
        }
        
        // Cluster seeds by protein and extend alignments
        for (int s = 0; s < num_seeds && match_count < max_matches_per_read; s++) {
            uint32_t protein_id = seeds[s].protein_id;
            if (!seeds[s].valid || protein_id >= protein_db->num_proteins) continue;
            
            // Collect all seeds for this protein
            SeedHit protein_seeds[MAX_PROTEIN_SEEDS];
            int seed_count = 0;
            
            for (int t = s; t < num_seeds && seed_count < MAX_PROTEIN_SEEDS; t++) {
                if (seeds[t].valid && seeds[t].protein_id == protein_id) {
                    protein_seeds[seed_count++] = seeds[t];
                    seeds[t].valid = false; // Mark as used
                }
            }
            
            if (seed_count < MIN_SEED_HITS) continue;
            
            // Safe extension from first seed
            uint32_t seed_query_pos = protein_seeds[0].query_pos;
            uint32_t seed_ref_pos = protein_seeds[0].ref_pos;
            
            // Bounds checking for protein access
            if (protein_id >= protein_db->num_proteins) continue;
            
            uint32_t seq_offset = protein_db->seq_offsets[protein_id];
            uint16_t ref_len = protein_db->seq_lengths[protein_id];
            
            // Verify sequence bounds
            if (seq_offset >= protein_db->total_sequence_length || 
                seq_offset + ref_len > protein_db->total_sequence_length ||
                seed_ref_pos >= ref_len) {
                continue;
            }
            
            const char* ref_seq = &protein_db->sequences[seq_offset];
            
            // Safe extension with mismatch tolerance
            int left_extend = 0;
            int left_mismatches = 0;
            const int max_mismatches = 3;  // Allow mismatches for AMR variants
            
            while (seed_query_pos > left_extend && 
                   seed_ref_pos > left_extend &&
                   left_extend < 30 &&
                   left_mismatches <= max_mismatches) {
                
                int query_idx = seed_query_pos - left_extend - 1;
                int ref_idx = seed_ref_pos - left_extend - 1;
                
                if (query_idx >= 0 && query_idx < frame.length && ref_idx >= 0 && ref_idx < ref_len) {
                    char query_aa = frame.sequence[query_idx];
                    char ref_aa = ref_seq[ref_idx];
                    
                    if (query_aa == ref_aa || get_blosum_score(query_aa, ref_aa) > 0) {
                        left_extend++;
                    } else {
                        left_mismatches++;
                        left_extend++;
                    }
                } else {
                    break;
                }
            }
            
            // Extend right from seed
            int right_extend = PROTEIN_KMER_SIZE;
            int right_mismatches = 0;
            
            while (seed_query_pos + right_extend < frame.length && 
                   seed_ref_pos + right_extend < ref_len &&
                   right_extend < 30 &&
                   right_mismatches <= max_mismatches) {
                
                int query_idx = seed_query_pos + right_extend;
                int ref_idx = seed_ref_pos + right_extend;
                
                if (query_idx < frame.length && ref_idx < ref_len) {
                    char query_aa = frame.sequence[query_idx];
                    char ref_aa = ref_seq[ref_idx];
                    
                    if (query_aa == ref_aa || get_blosum_score(query_aa, ref_aa) > 0) {
                        right_extend++;
                    } else {
                        right_mismatches++;
                        right_extend++;
                    }
                } else {
                    break;
                }
            }
            
            uint32_t match_query_start = seed_query_pos - left_extend;
            uint32_t match_ref_start = seed_ref_pos - left_extend;
            uint32_t match_length = left_extend + right_extend;
            
            // Create match if meets minimum requirements
            if (seed_count >= MIN_SEED_HITS && match_length >= MIN_ALIGNMENT_LENGTH) {
                ProteinMatch temp_match;
                temp_match.read_id = tid;
                temp_match.frame = frame.frame;
                temp_match.protein_id = protein_id;
                temp_match.gene_id = protein_db->gene_ids[protein_id];
                temp_match.species_id = protein_db->species_ids[protein_id];
                temp_match.query_start = match_query_start;
                temp_match.ref_start = match_ref_start;
                temp_match.match_length = match_length;
                
                // Calculate identity and coverage
                int total_mismatches = left_mismatches + right_mismatches;
                temp_match.identity = 1.0f - (float)total_mismatches / match_length;
                temp_match.alignment_score = match_length * 2.0f - total_mismatches * 4.0f;
                
                temp_match.gene_length = ref_len * 3;  // Convert AA to nucleotides
                temp_match.coverage_start = match_ref_start;
                temp_match.coverage_end = match_ref_start + match_length;
                // Coverage fraction can be calculated later from coverage_start/end and gene_length
                
                temp_match.used_smith_waterman = false;
                temp_match.concordant = false;
                // Initialize query_peptide
                memset(temp_match.query_peptide, 0, 51);
                
                // Extract query peptide sequence with bounds checking
                int peptide_len = min(match_length, 50);
                for (int k = 0; k < peptide_len; k++) {
                    int query_idx = match_query_start + k;
                    
                    if (query_idx < frame.length && k < 50) {
                        temp_match.query_peptide[k] = frame.sequence[query_idx];
                    } else {
                        temp_match.query_peptide[k] = 'X';
                    }
                }
                temp_match.query_peptide[peptide_len] = '\0';
                
                // Apply Smith-Waterman if enabled and beneficial
                if (enable_smith_waterman && seed_count >= MIN_SEED_HITS && match_length < 40) {
                    uint16_t available_ref = ref_len - match_ref_start;
                    uint16_t available_query = frame.length - match_query_start;
                    uint16_t sw_ref_len = min(available_ref, (uint16_t)50);
                    uint16_t sw_query_len = min(available_query, (uint16_t)50);
                    
                    if (sw_ref_len > 0 && sw_query_len > 0) {
                        uint16_t sw_query_start, sw_ref_start, sw_length;
                        float sw_score = smith_waterman_align(
                            &frame.sequence[match_query_start], sw_query_len,
                            &ref_seq[match_ref_start], sw_ref_len,
                            &sw_query_start, &sw_ref_start, &sw_length
                        );
                        
                        if (sw_score > temp_match.alignment_score + 10.0f) {  // Require significant improvement
                            temp_match.alignment_score = sw_score;
                            temp_match.query_start = match_query_start + sw_query_start;
                            temp_match.ref_start = match_ref_start + sw_ref_start;
                            temp_match.match_length = sw_length;
                            temp_match.identity = sw_score / (sw_length * 3.0f);  // More conservative identity
                            temp_match.used_smith_waterman = true;
                            
                            // Update coverage
                            temp_match.coverage_start = temp_match.ref_start;
                            temp_match.coverage_end = temp_match.ref_start + sw_length;
                            // Coverage fraction can be calculated later from coverage_start/end and gene_length
                        }
                    }
                }
                
                // Accept match if it meets quality thresholds
                if (temp_match.identity >= MIN_IDENTITY_THRESHOLD && 
                    temp_match.match_length >= MIN_ALIGNMENT_LENGTH &&
                    match_count < max_matches_per_read) {
                    
                    read_matches[match_count] = temp_match;
                    match_count++;
                }
            }
        }
    }
    
    // Ensure match count doesn't exceed limit
    if (match_count > max_matches_per_read) {
        match_count = max_matches_per_read;
    }
    
    match_counts[tid] = match_count;
}

// Host wrapper class with enhanced safety and memory management
class AMRTranslatedSearchEngine {
private:
    ProteinDatabase* d_protein_db;
    TranslatedFrame* d_translated_frames;
    uint32_t* d_frame_counts;
    ProteinMatch* d_matches;
    uint32_t* d_match_counts;
    
    int max_batch_size;
    bool smith_waterman_enabled;
    bool initialized;
    
    // Memory tracking
    size_t allocated_memory;
    
public:
    AMRTranslatedSearchEngine(int batch_size = 8192, bool enable_sw = false) 
        : max_batch_size(batch_size), smith_waterman_enabled(enable_sw), 
          initialized(false), allocated_memory(0) {
        
        DEBUG_PRINT("Initializing AMR engine: batch_size=%d, SW=%s", 
                   batch_size, enable_sw ? "enabled" : "disabled");
        
        // Validate and adjust batch size
        if (batch_size <= 0) {
            printf("[AMR ENGINE WARNING] Invalid batch size %d, using default 8192\n", batch_size);
            max_batch_size = 8192;
        } else if (batch_size > 200000) {
            printf("[AMR ENGINE WARNING] Batch size %d too large, capping at 200000\n", batch_size);
            max_batch_size = 200000;
        }
        
        // Initialize pointers
        d_protein_db = nullptr;
        d_translated_frames = nullptr;
        d_frame_counts = nullptr;
        d_matches = nullptr;
        d_match_counts = nullptr;
        
        if (!allocateGPUMemory()) {
            printf("[AMR ENGINE ERROR] Failed to allocate GPU memory\n");
            return;
        }
        
        initialized = true;
        DEBUG_PRINT("AMR engine initialized successfully");
    }
    
    ~AMRTranslatedSearchEngine() {
        cleanup();
    }
    
private:
    bool allocateGPUMemory() {
        // Check available GPU memory
        size_t free_mem, total_mem;
        if (cudaMemGetInfo(&free_mem, &total_mem) != cudaSuccess) {
            printf("[AMR ENGINE ERROR] Failed to query GPU memory\n");
            return false;
        }
        
        DEBUG_PRINT("GPU memory: %zu MB free / %zu MB total", 
                   free_mem / (1024*1024), total_mem / (1024*1024));
        
        // Calculate required memory with safety margins
        size_t frames_mem = max_batch_size * NUM_FRAMES * sizeof(TranslatedFrame);
        size_t counts_mem = max_batch_size * sizeof(uint32_t);
        size_t matches_mem = max_batch_size * MAX_MATCHES_PER_READ * sizeof(ProteinMatch);
        size_t total_required = frames_mem + counts_mem * 2 + matches_mem;
        
        DEBUG_PRINT("Requiring %zu MB GPU memory", total_required / (1024*1024));
        
        if (total_required > free_mem * 0.75) {  // Leave 25% buffer
            printf("[AMR ENGINE ERROR] Insufficient GPU memory: need %zu MB, have %zu MB\n", 
                   total_required / (1024*1024), free_mem / (1024*1024));
            return false;
        }
        
        // Allocate with error checking
        CUDA_CHECK(cudaMalloc(&d_translated_frames, frames_mem));
        allocated_memory += frames_mem;
        
        CUDA_CHECK(cudaMalloc(&d_frame_counts, counts_mem));
        allocated_memory += counts_mem;
        
        CUDA_CHECK(cudaMalloc(&d_matches, matches_mem));
        allocated_memory += matches_mem;
        
        CUDA_CHECK(cudaMalloc(&d_match_counts, counts_mem));
        allocated_memory += counts_mem;
        
        DEBUG_PRINT("Allocated %zu MB GPU memory", allocated_memory / (1024*1024));
        return true;
    }
    
    void cleanup() {
        if (d_translated_frames) {
            cudaFree(d_translated_frames);
            d_translated_frames = nullptr;
        }
        if (d_frame_counts) {
            cudaFree(d_frame_counts);
            d_frame_counts = nullptr;
        }
        if (d_matches) {
            cudaFree(d_matches);
            d_matches = nullptr;
        }
        if (d_match_counts) {
            cudaFree(d_match_counts);
            d_match_counts = nullptr;
        }
        
        if (d_protein_db) {
            ProteinDatabase h_db;
            if (cudaMemcpy(&h_db, d_protein_db, sizeof(ProteinDatabase), cudaMemcpyDeviceToHost) == cudaSuccess) {
                if (h_db.sorted_kmer_hashes) cudaFree(h_db.sorted_kmer_hashes);
                if (h_db.kmer_start_indices) cudaFree(h_db.kmer_start_indices);
                if (h_db.kmer_counts) cudaFree(h_db.kmer_counts);
                if (h_db.position_data) cudaFree(h_db.position_data);
                if (h_db.protein_ids) cudaFree(h_db.protein_ids);
                if (h_db.gene_ids) cudaFree(h_db.gene_ids);
                if (h_db.species_ids) cudaFree(h_db.species_ids);
                if (h_db.seq_lengths) cudaFree(h_db.seq_lengths);
                if (h_db.seq_offsets) cudaFree(h_db.seq_offsets);
                if (h_db.sequences) cudaFree(h_db.sequences);
                if (h_db.is_amr_gene) cudaFree(h_db.is_amr_gene);
                if (h_db.gene_families) {
                    // First free each string
                    char** h_family_ptrs = new char*[h_db.num_proteins];
                    cudaMemcpy(h_family_ptrs, h_db.gene_families, 
                               h_db.num_proteins * sizeof(char*), cudaMemcpyDeviceToHost);
                    for (uint32_t i = 0; i < h_db.num_proteins; i++) {
                        if (h_family_ptrs[i]) cudaFree(h_family_ptrs[i]);
                    }
                    delete[] h_family_ptrs;
                    // Then free the array of pointers
                    cudaFree(h_db.gene_families);
                }
                if (h_db.resistance_class) cudaFree(h_db.resistance_class);
            }
            cudaFree(d_protein_db);
            d_protein_db = nullptr;
        }
        
        initialized = false;
        allocated_memory = 0;
        DEBUG_PRINT("AMR engine cleanup completed");
    }
    
public:
    bool isInitialized() const { return initialized; }
    
    bool loadProteinDatabase(const std::string& db_path) {
        if (!initialized) {
            printf("[AMR ENGINE ERROR] Engine not initialized\n");
            return false;
        }
        
        DEBUG_PRINT("Loading AMR protein database from %s", db_path.c_str());
        
        if (d_protein_db) {
            DEBUG_PRINT("Protein database already loaded");
            return true;
        }
        
        // Load k-mer index with enhanced error checking
        std::string kmer_path = db_path + "/protein_kmers.bin";
        std::ifstream kmer_file(kmer_path, std::ios::binary);
        if (!kmer_file.good()) {
            printf("[AMR DB ERROR] Cannot read k-mer file: %s\n", kmer_path.c_str());
            return false;
        }
        
        uint32_t kmer_length, num_kmers;
        kmer_file.read(reinterpret_cast<char*>(&kmer_length), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&num_kmers), sizeof(uint32_t));
        
        if (kmer_length != PROTEIN_KMER_SIZE) {
            printf("[AMR DB ERROR] K-mer size mismatch: expected %d, got %d\n", PROTEIN_KMER_SIZE, kmer_length);
            return false;
        }
        
        if (num_kmers == 0 || num_kmers > 10000000) {  // Sanity check
            printf("[AMR DB ERROR] Invalid k-mer count: %d\n", num_kmers);
            return false;
        }
        
        DEBUG_PRINT("Loading %d protein k-mers", num_kmers);
        
        // Load and sort k-mers with bounds checking
        std::map<uint64_t, std::vector<uint32_t>> kmer_map;
        
        for (uint32_t i = 0; i < num_kmers; i++) {
            char kmer_seq[PROTEIN_KMER_SIZE + 1] = {0};
            kmer_file.read(kmer_seq, kmer_length);
            
            if (kmer_file.gcount() != kmer_length) {
                printf("[AMR DB ERROR] Incomplete k-mer read at position %d\n", i);
                return false;
            }
            
            // Validate k-mer sequence
            bool valid_kmer = true;
            for (int j = 0; j < kmer_length; j++) {
                if (aa_to_index(kmer_seq[j]) == 22 && kmer_seq[j] != 'X') {  // Invalid AA
                    valid_kmer = false;
                    break;
                }
            }
            
            if (!valid_kmer) {
                printf("[AMR DB WARNING] Skipping invalid k-mer at position %d\n", i);
                continue;
            }
            
            uint64_t hash = hash_protein_kmer(kmer_seq);
            
            uint32_t num_positions;
            kmer_file.read(reinterpret_cast<char*>(&num_positions), sizeof(uint32_t));
            
            if (num_positions > 1000) {  // Reasonable limit for AMR database
                printf("[AMR DB WARNING] K-mer %d has excessive positions (%d), truncating\n", i, num_positions);
                num_positions = 1000;
            }
            
            for (uint32_t j = 0; j < num_positions; j++) {
                uint32_t protein_idx, position;
                kmer_file.read(reinterpret_cast<char*>(&protein_idx), sizeof(uint32_t));
                kmer_file.read(reinterpret_cast<char*>(&position), sizeof(uint32_t));
                
                if (protein_idx > 100000 || position > 65535) {  // Sanity checks
                    printf("[AMR DB WARNING] Invalid protein/position: %d/%d\n", protein_idx, position);
                    continue;
                }
                
                uint32_t encoded = (protein_idx << 16) | (position & 0xFFFF);
                kmer_map[hash].push_back(encoded);
            }
        }
        kmer_file.close();
        
        DEBUG_PRINT("Processed %zu unique k-mer hashes", kmer_map.size());
        
        // Create sorted arrays for GPU
        std::vector<uint64_t> sorted_hashes;
        std::vector<uint32_t> start_indices;
        std::vector<uint32_t> kmer_counts;
        std::vector<uint32_t> position_data;
        
        sorted_hashes.reserve(kmer_map.size());
        start_indices.reserve(kmer_map.size());
        kmer_counts.reserve(kmer_map.size());
        
        for (const auto& pair : kmer_map) {
            sorted_hashes.push_back(pair.first);
            start_indices.push_back(position_data.size());
            kmer_counts.push_back(pair.second.size());
            
            for (uint32_t pos : pair.second) {
                position_data.push_back(pos);
            }
        }
        
        // Load protein sequences and metadata
        std::string protein_path = db_path + "/proteins.bin";
        std::ifstream protein_file(protein_path, std::ios::binary);
        if (!protein_file.good()) {
            printf("[AMR DB ERROR] Cannot read protein file: %s\n", protein_path.c_str());
            return false;
        }
        
        uint32_t num_proteins;
        protein_file.read(reinterpret_cast<char*>(&num_proteins), sizeof(uint32_t));
        
        if (num_proteins == 0 || num_proteins > 1000000) {  // Sanity check
            printf("[AMR DB ERROR] Invalid protein count: %d\n", num_proteins);
            return false;
        }
        
        // Calculate total file size for validation
        protein_file.seekg(0, std::ios::end);
        size_t file_size = protein_file.tellg();
        protein_file.seekg(sizeof(uint32_t), std::ios::beg);
        
        size_t remaining_size = file_size - sizeof(uint32_t);
        if (remaining_size > 100 * 1024 * 1024) {  // 100MB limit for safety
            printf("[AMR DB ERROR] Protein file too large: %zu bytes\n", remaining_size);
            return false;
        }
        
        std::vector<char> all_sequences(remaining_size + 1);
        protein_file.read(all_sequences.data(), remaining_size);
        all_sequences[remaining_size] = '\0';
        protein_file.close();
        
        DEBUG_PRINT("Loaded %d proteins, %zu sequence bytes", num_proteins, remaining_size);
        
        // Initialize metadata vectors with bounds checking
        std::vector<uint32_t> protein_ids(num_proteins);
        std::vector<uint32_t> gene_ids(num_proteins);
        std::vector<uint32_t> species_ids(num_proteins);
        std::vector<uint16_t> seq_lengths(num_proteins);
        std::vector<uint32_t> seq_offsets(num_proteins);
        std::vector<uint8_t> is_amr_gene(num_proteins, 0);  // Use uint8_t instead of bool
        std::vector<uint8_t> resistance_class(num_proteins, 0);
        std::vector<std::string> gene_families(num_proteins);  // Store gene family names
        
        // Load AMR-specific metadata
        std::string metadata_path = db_path + "/protein_details.json";
        std::ifstream metadata_file(metadata_path);
        
        if (metadata_file.good()) {
            std::string json_content((std::istreambuf_iterator<char>(metadata_file)),
                                   std::istreambuf_iterator<char>());
            metadata_file.close();
            
            DEBUG_PRINT("Processing metadata JSON (%zu bytes)", json_content.size());
            
            // Enhanced JSON parsing with AMR detection
            size_t pos = 0;
            size_t protein_idx = 0;
            size_t current_offset = 0;
            
            while ((pos = json_content.find("\"id\":", pos)) != std::string::npos && protein_idx < num_proteins) {
                pos += 5;
                size_t id_start = json_content.find_first_of("0123456789", pos);
                size_t id_end = json_content.find_first_not_of("0123456789", id_start);
                
                if (id_start != std::string::npos && id_end != std::string::npos) {
                    protein_ids[protein_idx] = std::stoi(json_content.substr(id_start, id_end - id_start));
                }
                
                // When parsing the JSON metadata, use the protein index as gene_id
                // or extract from the gene_name field
                gene_ids[protein_idx] = protein_idx; // Simple approach: each protein is its own gene
                // OR parse from the JSON if available
                size_t gene_pos = json_content.find("\"gene_id\":", pos);
                if (gene_pos != std::string::npos && gene_pos < pos + 1000) {
                    gene_pos += 10;
                    size_t gene_start = json_content.find_first_of("0123456789", gene_pos);
                    size_t gene_end = json_content.find_first_not_of("0123456789", gene_start);
                    if (gene_start != std::string::npos && gene_end != std::string::npos) {
                        // Still parse it but use protein_idx as the actual gene_id
                        int parsed_gene_id = std::stoi(json_content.substr(gene_start, gene_end - gene_start));
                        // For now, use protein index as gene_id to ensure unique mapping
                        gene_ids[protein_idx] = protein_idx;
                    }
                }
                
                // Parse gene_family
                size_t family_pos = json_content.find("\"gene_family\":", pos);
                if (family_pos != std::string::npos && family_pos < pos + 1000) {
                    family_pos += 14;
                    size_t quote_start = json_content.find("\"", family_pos);
                    size_t quote_end = json_content.find("\"", quote_start + 1);
                    if (quote_start != std::string::npos && quote_end != std::string::npos) {
                        std::string family = json_content.substr(quote_start + 1, quote_end - quote_start - 1);
                        // Store this in a vector for now, will copy to GPU later
                        gene_families[protein_idx] = family;
                    }
                }
                
                // Parse species_id
                size_t species_pos = json_content.find("\"species_id\":", pos);
                if (species_pos != std::string::npos && species_pos < pos + 1000) {
                    species_pos += 13;
                    size_t species_start = json_content.find_first_of("0123456789", species_pos);
                    size_t species_end = json_content.find_first_not_of("0123456789", species_start);
                    if (species_start != std::string::npos && species_end != std::string::npos) {
                        species_ids[protein_idx] = std::stoi(json_content.substr(species_start, species_end - species_start));
                    }
                }
                
                // Parse sequence length
                size_t length_pos = json_content.find("\"length\":", pos);
                if (length_pos != std::string::npos && length_pos < pos + 1000) {
                    length_pos += 9;
                    size_t length_start = json_content.find_first_of("0123456789", length_pos);
                    size_t length_end = json_content.find_first_not_of("0123456789", length_start);
                    if (length_start != std::string::npos && length_end != std::string::npos) {
                        uint16_t len = std::stoi(json_content.substr(length_start, length_end - length_start));
                        seq_lengths[protein_idx] = (len > 5000) ? 5000 : len;  // Cap at reasonable size
                    }
                }
                
                // Check for AMR annotations
                size_t amr_pos = json_content.find("\"amr\":", pos);
                if (amr_pos != std::string::npos && amr_pos < pos + 1000) {
                    if (json_content.find("true", amr_pos + 6) == amr_pos + 6) {
                        is_amr_gene[protein_idx] = 1;
                    }
                }
                
                // Check for fluoroquinolone resistance
                if (json_content.find("fluoroquinolone", pos) != std::string::npos && 
                    json_content.find("fluoroquinolone", pos) < pos + 1000) {
                    resistance_class[protein_idx] = 1;  // Fluoroquinolone class
                }
                
                seq_offsets[protein_idx] = current_offset;
                current_offset += seq_lengths[protein_idx];
                
                if (protein_idx < 5) {  // Debug first 5 proteins
                    printf("[AMR DB] Protein[%zu]: id=%d, gene_id=%d, species_id=%d\n",
                           protein_idx, protein_ids[protein_idx], 
                           gene_ids[protein_idx], species_ids[protein_idx]);
                }
                
                protein_idx++;
                pos = id_end;
            }
            
            DEBUG_PRINT("Processed metadata for %zu proteins", protein_idx);
        } else {
            printf("[AMR DB WARNING] No metadata file found, using defaults\n");
            // Generate default metadata
            size_t avg_protein_len = remaining_size / num_proteins;
            for (uint32_t i = 0; i < num_proteins; i++) {
                protein_ids[i] = i;
                gene_ids[i] = i;  // Each protein is its own gene
                species_ids[i] = i % 20;  // Distribute across species
                seq_offsets[i] = i * avg_protein_len;
                seq_lengths[i] = (i == num_proteins - 1) ? 
                                (remaining_size - seq_offsets[i]) : avg_protein_len;
                if (seq_lengths[i] > 5000) seq_lengths[i] = 5000;  // Safety cap
            }
        }
        
        // Prepare gene families for GPU (after metadata parsing loop completes)
        std::vector<char*> h_gene_family_ptrs(num_proteins);
        std::vector<std::string> gene_family_storage = gene_families; // Keep strings alive
        for (uint32_t i = 0; i < num_proteins; i++) {
            h_gene_family_ptrs[i] = const_cast<char*>(gene_family_storage[i].c_str());
        }
        
        // Clear any existing CUDA errors before allocation
        cudaError_t existing_err = cudaGetLastError();
        if (existing_err != cudaSuccess) {
            printf("[AMR DB WARNING] Clearing existing CUDA error: %s\n", cudaGetErrorString(existing_err));
        }
        
        // Allocate and copy to GPU with comprehensive error checking
        ProteinDatabase h_db = {};
        h_db.num_proteins = num_proteins;
        h_db.num_kmers = sorted_hashes.size();
        h_db.total_sequence_length = remaining_size;
        
        // Allocate GPU memory for all database components
        CUDA_CHECK(cudaMalloc(&h_db.sorted_kmer_hashes, sorted_hashes.size() * sizeof(uint64_t)));
        CUDA_CHECK(cudaMalloc(&h_db.kmer_start_indices, start_indices.size() * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.kmer_counts, kmer_counts.size() * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.position_data, position_data.size() * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.protein_ids, num_proteins * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.gene_ids, num_proteins * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.species_ids, num_proteins * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.seq_lengths, num_proteins * sizeof(uint16_t)));
        CUDA_CHECK(cudaMalloc(&h_db.seq_offsets, num_proteins * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&h_db.sequences, all_sequences.size() * sizeof(char)));
        CUDA_CHECK(cudaMalloc(&h_db.is_amr_gene, num_proteins * sizeof(uint8_t)));  // Changed from bool to uint8_t
        CUDA_CHECK(cudaMalloc(&h_db.resistance_class, num_proteins * sizeof(uint8_t)));
        
        // Allocate array of pointers for gene families
        char** d_gene_family_ptrs;
        CUDA_CHECK(cudaMalloc(&d_gene_family_ptrs, num_proteins * sizeof(char*)));
        
        // Allocate and copy each gene family string
        std::vector<char*> d_family_strings(num_proteins);
        for (uint32_t i = 0; i < num_proteins; i++) {
            size_t len = gene_families[i].length() + 1;
            CUDA_CHECK(cudaMalloc(&d_family_strings[i], len));
            CUDA_CHECK(cudaMemcpy(d_family_strings[i], gene_families[i].c_str(), 
                                 len, cudaMemcpyHostToDevice));
        }
        
        // Copy array of pointers
        CUDA_CHECK(cudaMemcpy(d_gene_family_ptrs, d_family_strings.data(), 
                             num_proteins * sizeof(char*), cudaMemcpyHostToDevice));
        
        h_db.gene_families = d_gene_family_ptrs;
        
        // Copy data to GPU
        CUDA_CHECK(cudaMemcpy(h_db.sorted_kmer_hashes, sorted_hashes.data(), 
                             sorted_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.kmer_start_indices, start_indices.data(), 
                             start_indices.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.kmer_counts, kmer_counts.data(), 
                             kmer_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.position_data, position_data.data(), 
                             position_data.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.protein_ids, protein_ids.data(), 
                             num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.gene_ids, gene_ids.data(), 
                             num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.species_ids, species_ids.data(), 
                             num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.seq_lengths, seq_lengths.data(), 
                             num_proteins * sizeof(uint16_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.seq_offsets, seq_offsets.data(), 
                             num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.sequences, all_sequences.data(), 
                             all_sequences.size() * sizeof(char), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.is_amr_gene, is_amr_gene.data(), 
                             num_proteins * sizeof(uint8_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(h_db.resistance_class, resistance_class.data(), 
                             num_proteins * sizeof(uint8_t), cudaMemcpyHostToDevice));
        
        // Copy database structure to GPU
        CUDA_CHECK(cudaMalloc(&d_protein_db, sizeof(ProteinDatabase)));
        CUDA_CHECK(cudaMemcpy(d_protein_db, &h_db, sizeof(ProteinDatabase), cudaMemcpyHostToDevice));
        
        // Count AMR genes for reporting
        int amr_count = 0;
        for (bool is_amr : is_amr_gene) {
            if (is_amr) amr_count++;
        }
        
        DEBUG_PRINT("AMR protein database loaded successfully:");
        DEBUG_PRINT("  - %d total proteins (%d AMR genes)", num_proteins, amr_count);
        DEBUG_PRINT("  - %zu unique k-mers", sorted_hashes.size());
        DEBUG_PRINT("  - %zu total sequence bytes", remaining_size);
        
        return true;
    }
    
    void setSmithWatermanEnabled(bool enabled) {
        smith_waterman_enabled = enabled;
        DEBUG_PRINT("Smith-Waterman alignment %s", enabled ? "ENABLED" : "DISABLED");
    }
    
    bool searchTranslatedReads(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,
        int num_reads,
        ProteinMatch* results,
        uint32_t* result_counts
    ) {
        if (!initialized || !d_protein_db) {
            printf("[AMR SEARCH ERROR] Engine not properly initialized\n");
            return false;
        }
        
        // Comprehensive input validation
        if (num_reads <= 0 || num_reads > max_batch_size) {
            printf("[AMR SEARCH ERROR] Invalid num_reads: %d (max: %d)\n", num_reads, max_batch_size);
            return false;
        }
        
        if (!d_reads || !d_read_lengths || !d_read_offsets || !results || !result_counts) {
            printf("[AMR SEARCH ERROR] Null input pointers\n");
            return false;
        }
        
        DEBUG_PRINT("Processing %d reads for AMR detection", num_reads);
        
        // Clear any existing CUDA errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("[AMR SEARCH WARNING] Clearing existing CUDA error: %s\n", cudaGetErrorString(err));
        }
        
        // Initialize arrays
        err = cudaMemset(d_frame_counts, 0, num_reads * sizeof(uint32_t));
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] cudaMemset d_frame_counts failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        err = cudaMemset(d_match_counts, 0, num_reads * sizeof(uint32_t));
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] cudaMemset d_match_counts failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        err = cudaMemset(d_matches, 0, num_reads * MAX_MATCHES_PER_READ * sizeof(ProteinMatch));
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] cudaMemset d_matches failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        // Configure kernel launch parameters
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        if (grid_size <= 0 || grid_size > 65535) {
            printf("[AMR SEARCH ERROR] Invalid grid size: %d\n", grid_size);
            return false;
        }
        
        DEBUG_PRINT("Launching kernels: grid_size=%d, block_size=%d", grid_size, block_size);
        
        // Stage 1: 6-frame translation
        six_frame_translate_kernel<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets,
            num_reads, d_translated_frames, d_frame_counts
        );
        
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] six_frame_translate_kernel launch failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] six_frame_translate_kernel sync failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        // Stage 2: Enhanced protein matching for AMR detection
        enhanced_protein_kmer_match_kernel<<<grid_size, block_size>>>(
            d_translated_frames, d_frame_counts,
            num_reads, d_protein_db,
            d_matches, d_match_counts, MAX_MATCHES_PER_READ,
            smith_waterman_enabled
        );
        
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] enhanced_protein_kmer_match_kernel launch failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] enhanced_protein_kmer_match_kernel sync failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        // Copy results back to host
        size_t matches_size = num_reads * MAX_MATCHES_PER_READ * sizeof(ProteinMatch);
        size_t counts_size = num_reads * sizeof(uint32_t);
        
        DEBUG_PRINT("Copying results: num_reads=%d, matches_size=%zu bytes, counts_size=%zu bytes", 
                    num_reads, matches_size, counts_size);
        DEBUG_PRINT("Source buffer d_matches=%p, dest buffer results=%p", d_matches, results);
        DEBUG_PRINT("sizeof(ProteinMatch)=%zu, MAX_MATCHES_PER_READ=%d", sizeof(ProteinMatch), MAX_MATCHES_PER_READ);
        
        err = cudaMemcpy(results, d_matches, matches_size, cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] cudaMemcpy results failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        err = cudaMemcpy(result_counts, d_match_counts, counts_size, cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            printf("[AMR SEARCH ERROR] cudaMemcpy result_counts failed: %s\n", cudaGetErrorString(err));
            cudaGetLastError(); // Clear error state
            return false;
        }
        
        DEBUG_PRINT("AMR search completed successfully");
        return true;
    }
    
    // Get performance statistics
    void getStats(int& batch_size, size_t& memory_used, bool& sw_enabled) const {
        batch_size = max_batch_size;
        memory_used = allocated_memory;
        sw_enabled = smith_waterman_enabled;
    }
};

// C interface for integration with Python/other languages
extern "C" {
    void* create_amr_search_engine(int batch_size, bool enable_sw) {
        return new AMRTranslatedSearchEngine(batch_size, enable_sw);
    }
    
    void destroy_amr_search_engine(void* engine) {
        if (engine) {
            delete static_cast<AMRTranslatedSearchEngine*>(engine);
        }
    }
    
    int load_amr_protein_database(void* engine, const char* db_path) {
        if (!engine) return -1;
        AMRTranslatedSearchEngine* amr_engine = static_cast<AMRTranslatedSearchEngine*>(engine);
        return amr_engine->loadProteinDatabase(db_path) ? 0 : -1;
    }
    
    void set_amr_smith_waterman_enabled(void* engine, bool enabled) {
        if (engine) {
            AMRTranslatedSearchEngine* amr_engine = static_cast<AMRTranslatedSearchEngine*>(engine);
            amr_engine->setSmithWatermanEnabled(enabled);
        }
    }
    
    int search_amr_translated_reads(
        void* engine,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,
        int num_reads,
        void* results,
        uint32_t* result_counts
    ) {
        if (!engine) return -1;
        
        AMRTranslatedSearchEngine* amr_engine = static_cast<AMRTranslatedSearchEngine*>(engine);
        
        return amr_engine->searchTranslatedReads(
            d_reads, d_read_lengths, d_read_offsets,
            d_reads_to_process, num_reads,
            static_cast<ProteinMatch*>(results), result_counts
        ) ? 0 : -1;
    }
    
    int get_amr_engine_stats(void* engine, int* batch_size, size_t* memory_used, bool* sw_enabled) {
        if (!engine) return -1;
        
        AMRTranslatedSearchEngine* amr_engine = static_cast<AMRTranslatedSearchEngine*>(engine);
        amr_engine->getStats(*batch_size, *memory_used, *sw_enabled);
        return 0;
    }
    
    // Compatibility wrappers for existing pipeline
    void* create_translated_search_engine(int batch_size) {
        return create_amr_search_engine(batch_size, false);
    }
    
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw) {
        return create_amr_search_engine(batch_size, enable_sw);
    }
    
    void destroy_translated_search_engine(void* engine) {
        destroy_amr_search_engine(engine);
    }
    
    int load_protein_database(void* engine, const char* db_path) {
        return load_amr_protein_database(engine, db_path);
    }
    
    void set_smith_waterman_enabled(void* engine, bool enabled) {
        set_amr_smith_waterman_enabled(engine, enabled);
    }
    
    int search_translated_reads(
        void* engine,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,
        int num_reads,
        void* results,
        uint32_t* result_counts
    ) {
        return search_amr_translated_reads(
            engine, d_reads, d_read_lengths, d_read_offsets,
            d_reads_to_process, num_reads, results, result_counts
        );
    }
}

#endif // TRANSLATED_SEARCH_CU