// translated_search.cu
// GPU-accelerated 6-frame translation and protein search for resistance detection
// Enhanced with 5-mer seeding, clustering, extension, and optional Smith-Waterman

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

namespace cg = cooperative_groups;

// Debug macros
#define DEBUG_TRANS 0
#define DEBUG_PRINT(fmt, ...) if(DEBUG_TRANS) { printf("[TRANS DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

// Enhanced constants for 5-mer approach
#define CODON_SIZE 3
#define NUM_FRAMES 6
#define MAX_PROTEIN_LENGTH 200
#define PROTEIN_KMER_SIZE 5         // Increased from 3 to 5 amino acids
#define MIN_PEPTIDE_LENGTH 20       // Minimum peptide length to consider
#define MIN_SEED_HITS 1            // Require only one k-mer hit for extension
#define EXTENSION_THRESHOLD 15     // Minimum amino acids for valid match
#define MIN_IDENTITY_THRESHOLD 0.6f // Minimum 60% identity for valid protein match (lowered for more hits)
#define SW_SCORE_THRESHOLD 40.0f   // Lowered threshold for Smith-Waterman alignment to detect mutations
#define AA_ALPHABET_SIZE 24

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

// BLOSUM62 matrix (simplified, key values)
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

// Amino acid to index mapping
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

// helper function to check if alignment covers QRDR regions
// Updated to use correct gene_ids from the database
__device__ bool covers_qrdr_region(uint32_t gene_id, uint16_t ref_start, uint16_t match_length) {
    uint16_t ref_end = ref_start + match_length;
    
    // Based on protein_details.json:
    // gene_id 1 = gyrA
    // gene_id 2 = parC
    // gene_id 0 = grlA (S. aureus equivalent of gyrA)
    
    if (gene_id == 1) {  // gyrA
        // gyrA QRDR positions (0-based): S83 (pos 82) and D87 (pos 86)
        return (ref_start <= 82 && ref_end > 82) || (ref_start <= 86 && ref_end > 86);
    } else if (gene_id == 2) {  // parC
        // parC QRDR positions (0-based): S80 (pos 79) and E84 (pos 83)
        return (ref_start <= 79 && ref_end > 79) || (ref_start <= 83 && ref_end > 83);
    } else if (gene_id == 0) {  // grlA (S. aureus)
        // grlA QRDR positions similar to gyrA: S84 (pos 83) and E88 (pos 87)
        return (ref_start <= 83 && ref_end > 83) || (ref_start <= 87 && ref_end > 87);
    }
    
    return false;
}


// Get BLOSUM score
__device__ inline float get_blosum_score(char aa1, char aa2) {
    int idx1 = aa_to_index(aa1);
    int idx2 = aa_to_index(aa2);
    return BLOSUM62_SCORES[idx1 * 24 + idx2];
}

// Base encoding for translation
__device__ inline int base_to_index(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;  // Invalid base
    }
}

// Translate codon to amino acid
__device__ inline char translate_codon(const char* codon) {
    int idx = 0;
    for (int i = 0; i < 3; i++) {
        int base = base_to_index(codon[i]);
        if (base < 0) return 'X';  // Unknown
        idx = (idx << 2) | base;
    }
    return GENETIC_CODE[idx];
}

// Optimized hash function for 5-mer protein k-mers
__device__ inline uint64_t hash_protein_5mer(const char* kmer) {
    uint64_t hash = 0;
    const uint64_t prime = 31;
    
    for (int i = 0; i < PROTEIN_KMER_SIZE; i++) {
        // Map amino acid to value preserving chemical properties
        uint8_t aa_val = aa_to_index(kmer[i]);
        hash = hash * prime + aa_val;
    }
    
    return hash;
}

// Structure for translated read
struct TranslatedFrame {
    char sequence[MAX_PROTEIN_LENGTH];
    uint16_t length;
    int8_t frame;  // -3, -2, -1, +1, +2, +3
    uint16_t start_pos;  // Start position in original read
};

// Enhanced structure for protein match with mutation details
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
    uint8_t num_mutations;
    uint8_t mutation_positions[10];  // Up to 10 mutations
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
    bool used_smith_waterman;  // Flag indicating if SW was used
    char query_peptide[51];  // Store aligned peptide sequence (up to 50 AA + null terminator)
    bool is_qrdr_alignment;  // Flag for QRDR region alignment
};

// Enhanced protein database structure with sorted k-mer index
struct ProteinDatabase {
    uint32_t num_proteins;
    uint32_t num_kmers;
    
    // Sorted k-mer index for binary search
    uint64_t* sorted_kmer_hashes;
    uint32_t* kmer_start_indices;  // Start index in position array
    uint32_t* kmer_counts;         // Number of positions per k-mer
    uint32_t* position_data;       // Encoded protein ID + position
    
    // Protein metadata
    uint32_t* protein_ids;
    uint32_t* gene_ids;
    uint32_t* species_ids;
    uint16_t* seq_lengths;
    
    // Reference sequences for Smith-Waterman
    char* sequences;
    uint32_t* seq_offsets;
};

// Helper structure for seed clustering
struct SeedHit {
    uint32_t protein_id;
    uint32_t query_pos;
    uint32_t ref_pos;
    float score;
};

// Binary search for protein k-mer
__device__ inline int binary_search_protein_kmer(
    const ProteinDatabase* db,
    uint64_t target_hash
) {
    int left = 0;
    int right = db->num_kmers - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
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

// Smith-Waterman local alignment (simplified for GPU)
__device__ float smith_waterman_align(
    const char* query, uint16_t query_len,
    const char* ref, uint16_t ref_len,
    uint16_t* best_query_start,
    uint16_t* best_ref_start,
    uint16_t* best_length,
    uint8_t* mutations,
    char* ref_mutations,
    char* query_mutations,
    float* blosum_mutations,
    uint8_t* num_mutations
) {
    // Simplified banded Smith-Waterman for GPU (more lenient scoring)
    const int GAP_OPEN = -2;
    const int GAP_EXTEND = -1;
    const int BAND_WIDTH = 10;
    
    // Use local arrays for alignments - increased to handle full 150bp reads (50 AA)
    // For banded alignment, we only need a strip of the matrix
    const int MAX_ALIGN_LEN = 60; // Support up to 60 AA (180 bp)
    float H[MAX_ALIGN_LEN][MAX_ALIGN_LEN] = {0}; // Score matrix
    
    if (query_len > MAX_ALIGN_LEN-1 || ref_len > MAX_ALIGN_LEN-1) {
        // Fall back to simple scoring for very large sequences
        float score = 0.0f;
        int matches = 0;
        int aligned_len = min(query_len, ref_len);
        
        for (int i = 0; i < aligned_len; i++) {
            float blosum = get_blosum_score(query[i], ref[i]);
            score += blosum;
            if (blosum > 0) matches++;
        }
        
        *best_query_start = 0;
        *best_ref_start = 0;
        *best_length = aligned_len;
        *num_mutations = 0;
        
        return score;
    }
    
    float max_score = 0.0f;
    int max_i = 0, max_j = 0;
    
    // Fill scoring matrix with banding optimization
    // Since we expect high similarity, only compute cells near the diagonal
    for (int i = 1; i <= query_len; i++) {
        // Banded alignment: only fill cells within BAND_WIDTH of diagonal
        int j_start = max(1, i - BAND_WIDTH);
        int j_end = min((int)ref_len, i + BAND_WIDTH);
        
        for (int j = j_start; j <= j_end; j++) {
            float match = H[i-1][j-1] + get_blosum_score(query[i-1], ref[j-1]);
            float delete_gap = (j > 1) ? H[i-1][j] + GAP_OPEN : GAP_OPEN;
            float insert_gap = (i > 1) ? H[i][j-1] + GAP_OPEN : GAP_OPEN;
            
            H[i][j] = fmaxf(0.0f, fmaxf(match, fmaxf(delete_gap, insert_gap)));
            
            if (H[i][j] > max_score) {
                max_score = H[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback to find alignment
    int i = max_i, j = max_j;
    int align_len = 0;
    uint8_t mut_count = 0;
    
    while (i > 0 && j > 0 && H[i][j] > 0 && align_len < MAX_ALIGN_LEN-1) {
        if (i > 0 && j > 0 && H[i][j] == H[i-1][j-1] + get_blosum_score(query[i-1], ref[j-1])) {
            // Match or mismatch
            if (query[i-1] != ref[j-1] && mut_count < 10) {
                mutations[mut_count] = j-1; // Reference position
                ref_mutations[mut_count] = ref[j-1];
                query_mutations[mut_count] = query[i-1];
                blosum_mutations[mut_count] = get_blosum_score(query[i-1], ref[j-1]);
                mut_count++;
            }
            i--; j--;
        } else if (i > 0 && H[i][j] == H[i-1][j] + GAP_OPEN) {
            // Deletion in reference
            i--;
        } else if (j > 0 && H[i][j] == H[i][j-1] + GAP_OPEN) {
            // Insertion in reference  
            j--;
        } else {
            break;
        }
        align_len++;
    }
    
    *best_query_start = i;
    *best_ref_start = j;
    *best_length = align_len;
    *num_mutations = mut_count;
    
    return max_score;
}

// 6-frame translation kernel (unchanged)
__global__ void six_frame_translate_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const int num_reads,
    TranslatedFrame* translated_frames,
    uint32_t* frame_counts
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const char* read = reads + read_offsets[tid];
    const int read_len = read_lengths[tid];
    
    TranslatedFrame* read_frames = &translated_frames[tid * NUM_FRAMES];
    int valid_frames = 0;
    
    // Forward frames (+1, +2, +3)
    for (int frame = 0; frame < 3; frame++) {
        TranslatedFrame& tf = read_frames[valid_frames];
        tf.frame = frame + 1;
        tf.start_pos = frame;
        tf.length = 0;
        
        for (int pos = frame; pos + 2 < read_len; pos += 3) {
            if (tf.length >= MAX_PROTEIN_LENGTH - 1) break;
            
            char aa = translate_codon(&read[pos]);
            if (aa == '*') {
                if (tf.length >= MIN_PEPTIDE_LENGTH) {
                    tf.sequence[tf.length] = '\0';
                    valid_frames++;
                    
                    if (valid_frames < NUM_FRAMES) {
                        tf = read_frames[valid_frames];
                        tf.frame = frame + 1;
                        tf.start_pos = pos + 3;
                        tf.length = 0;
                    }
                } else {
                    tf.length = 0;
                    tf.start_pos = pos + 3;
                }
            } else {
                tf.sequence[tf.length++] = aa;
            }
        }
        
        if (tf.length >= MIN_PEPTIDE_LENGTH) {
            tf.sequence[tf.length] = '\0';
            valid_frames++;
        }
    }
    
    // Reverse complement frames (-1, -2, -3)
    for (int frame = 0; frame < 3 && valid_frames < NUM_FRAMES; frame++) {
        TranslatedFrame& tf = read_frames[valid_frames];
        tf.frame = -(frame + 1);
        tf.start_pos = read_len - frame - 1;
        tf.length = 0;
        
        for (int pos = read_len - frame - 1; pos >= 2; pos -= 3) {
            if (tf.length >= MAX_PROTEIN_LENGTH - 1) break;
            
            char rc_codon[3];
            for (int i = 0; i < 3; i++) {
                char base = read[pos - i];
                switch(base) {
                    case 'A': case 'a': rc_codon[i] = 'T'; break;
                    case 'T': case 't': rc_codon[i] = 'A'; break;
                    case 'G': case 'g': rc_codon[i] = 'C'; break;
                    case 'C': case 'c': rc_codon[i] = 'G'; break;
                    default: rc_codon[i] = 'N';
                }
            }
            
            char aa = translate_codon(rc_codon);
            if (aa == '*') {
                if (tf.length >= MIN_PEPTIDE_LENGTH) {
                    tf.sequence[tf.length] = '\0';
                    valid_frames++;
                    
                    if (valid_frames < NUM_FRAMES) {
                        tf = read_frames[valid_frames];
                        tf.frame = -(frame + 1);
                        tf.start_pos = pos - 3;
                        tf.length = 0;
                    }
                } else {
                    tf.length = 0;
                    tf.start_pos = pos - 3;
                }
            } else {
                tf.sequence[tf.length++] = aa;
            }
        }
        
        if (tf.length >= MIN_PEPTIDE_LENGTH) {
            tf.sequence[tf.length] = '\0';
            valid_frames++;
        }
    }
    
    frame_counts[tid] = valid_frames;
    
    if (tid == 0 && DEBUG_TRANS) {
        DEBUG_PRINT("Read 0: %d valid frames from %d bp", valid_frames, read_len);
    }
}

// Enhanced protein k-mer matching kernel with 5-mer seeding and extension
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
    if (tid >= num_reads) return;
    
    uint32_t num_frames = frame_counts[tid];
    if (num_frames == 0) {
        match_counts[tid] = 0;
        return;
    }
    
    const TranslatedFrame* read_frames = &translated_frames[tid * NUM_FRAMES];
    ProteinMatch* read_matches = &matches[tid * max_matches_per_read];
    uint32_t match_count = 0;
    
    // For each translated frame
    for (uint32_t frame_idx = 0; frame_idx < num_frames; frame_idx++) {
        const TranslatedFrame& frame = read_frames[frame_idx];
        
        if (frame.length < PROTEIN_KMER_SIZE) continue;
        
        // Collect 5-mer seed hits
        SeedHit seeds[100];
        int num_seeds = 0;
        
        // Find 5-mer seed matches
        for (int pos = 0; pos + PROTEIN_KMER_SIZE <= frame.length && num_seeds < 100; pos++) {
            uint64_t kmer_hash = hash_protein_5mer(&frame.sequence[pos]);
            
            int kmer_idx = binary_search_protein_kmer(protein_db, kmer_hash);
            if (kmer_idx >= 0) {
                uint32_t start_idx = protein_db->kmer_start_indices[kmer_idx];
                uint32_t count = protein_db->kmer_counts[kmer_idx];
                
                // Add hits for this k-mer (limit to avoid overflow)
                for (uint32_t i = 0; i < count && i < 5 && num_seeds < 100; i++) {
                    uint32_t encoded = protein_db->position_data[start_idx + i];
                    uint32_t protein_id = encoded >> 16;
                    uint32_t ref_pos = encoded & 0xFFFF;
                    
                    seeds[num_seeds] = {protein_id, (uint32_t)pos, ref_pos, 10.0f};
                    num_seeds++;
                }
            }
        }
        
        // Cluster seeds by protein and extend
        for (int s = 0; s < num_seeds && match_count < max_matches_per_read; s++) {
            uint32_t protein_id = seeds[s].protein_id;
            if (protein_id == UINT32_MAX) continue; // Already processed
            
            // Collect all seeds for this protein
            SeedHit protein_seeds[20];
            int seed_count = 0;
            
            for (int t = s; t < num_seeds && seed_count < 20; t++) {
                if (seeds[t].protein_id == protein_id) {
                    protein_seeds[seed_count++] = seeds[t];
                    seeds[t].protein_id = UINT32_MAX; // Mark as used
                }
            }
            
            if (seed_count < MIN_SEED_HITS) continue;
            
            // Safe extension from first seed - extend while maintaining alignment
            uint32_t seed_query_pos = protein_seeds[0].query_pos;
            uint32_t seed_ref_pos = protein_seeds[0].ref_pos;
            
            // Get reference sequence for this protein
            const char* ref_seq = &protein_db->sequences[protein_db->seq_offsets[protein_id]];
            uint16_t ref_len = protein_db->seq_lengths[protein_id];
            
            // Extend left from seed (allow some mismatches)
            int left_extend = 0;
            int left_mismatches = 0;
            const int max_mismatches = 2;  // Allow up to 2 mismatches during extension
            
            // Arrays to store mutation information
            uint8_t mutation_positions[10];
            char ref_aas[10];
            char query_aas[10];
            int num_mutations = 0;
            
            while (seed_query_pos - left_extend > 0 && 
                   seed_ref_pos - left_extend > 0 &&
                   left_extend < 20 &&
                   left_mismatches <= max_mismatches) {
                char query_aa = frame.sequence[seed_query_pos - left_extend - 1];
                char ref_aa = ref_seq[seed_ref_pos - left_extend - 1];
                
                if (query_aa == ref_aa) {
                    left_extend++;
                } else {
                    // Record mutation details
                    if (num_mutations < 10) {
                        mutation_positions[num_mutations] = left_extend - 1;
                        ref_aas[num_mutations] = ref_aa;
                        query_aas[num_mutations] = query_aa;
                        num_mutations++;
                    }
                    left_mismatches++;
                    left_extend++;  // Continue extending despite mismatch
                }
            }
            
            // Extend right from seed (allow some mismatches)
            int right_extend = PROTEIN_KMER_SIZE;  // Start after the k-mer
            int right_mismatches = 0;
            
            while (seed_query_pos + right_extend < frame.length && 
                   seed_ref_pos + right_extend < ref_len &&
                   right_extend < 50 &&
                   right_mismatches <= max_mismatches) {
                char query_aa = frame.sequence[seed_query_pos + right_extend];
                char ref_aa = ref_seq[seed_ref_pos + right_extend];
                
                if (query_aa == ref_aa) {
                    right_extend++;
                } else {
                    // Record mutation details
                    if (num_mutations < 10) {
                        mutation_positions[num_mutations] = seed_query_pos - (seed_query_pos - left_extend) + right_extend;
                        ref_aas[num_mutations] = ref_aa;
                        query_aas[num_mutations] = query_aa;
                        num_mutations++;
                    }
                    right_mismatches++;
                    right_extend++;  // Continue extending despite mismatch
                }
            }
            
            uint32_t min_query = seed_query_pos - left_extend;
            uint32_t min_ref = seed_ref_pos - left_extend;
            uint32_t query_span = left_extend + right_extend;
            uint32_t ref_span = query_span;  // Same span due to exact extension
            
            // Always proceed if we have enough seed hits
            if (seed_count >= MIN_SEED_HITS) {
                
                // Create temporary match to check identity
                ProteinMatch temp_match;
                temp_match.read_id = tid;
                temp_match.frame = frame.frame;
                temp_match.protein_id = protein_id;
                temp_match.gene_id = protein_db->gene_ids[protein_id];
                temp_match.species_id = protein_db->species_ids[protein_id];
                temp_match.query_start = min_query;
                temp_match.ref_start = min_ref;
                temp_match.match_length = query_span;  // Use the extended length
                
                // Calculate identity based on mismatches found during extension
                int total_mismatches = left_mismatches + right_mismatches;
                temp_match.identity = 1.0f - (float)total_mismatches / query_span;
                temp_match.alignment_score = query_span * 2.0f - total_mismatches * 4.0f;  // Penalize mismatches
                temp_match.num_mutations = num_mutations;
                temp_match.used_smith_waterman = false;
                
                // Copy mutation details
                for (int i = 0; i < num_mutations && i < 10; i++) {
                    temp_match.mutation_positions[i] = mutation_positions[i];
                    temp_match.ref_aas[i] = ref_aas[i];
                    temp_match.query_aas[i] = query_aas[i];
                }
                
                // Initialize new fields
                temp_match.is_qrdr_alignment = false;  // Will be set on host side
                // Extract peptide sequence from translated frame (extended region)
                int peptide_len = min(query_span, 50);
                for (int k = 0; k < peptide_len && k < 50; k++) {
                    if (min_query + k < frame.length) {
                        temp_match.query_peptide[k] = frame.sequence[min_query + k];
                    } else {
                        temp_match.query_peptide[k] = 'X';
                    }
                }
                temp_match.query_peptide[peptide_len] = '\0';
                
                // Debug: Print scoring info for first few matches
                if (tid == 0 && match_count < 3 && DEBUG_TRANS) {
                    DEBUG_PRINT("Match %d: seed_count=%d, score=%.1f, threshold=%.1f, SW enabled=%s", 
                               match_count, seed_count, temp_match.alignment_score, SW_SCORE_THRESHOLD, enable_smith_waterman ? "YES" : "NO");
                }
                
                // Apply Smith-Waterman to extend alignments and find more mutations
                // Run SW for all matches to get better extension (not just high-scoring ones)
                if (enable_smith_waterman && seed_count >= MIN_SEED_HITS) {
                    if (tid == 0 && match_count < 3 && DEBUG_TRANS) {
                        DEBUG_PRINT("Triggering Smith-Waterman: protein_id=%d, score=%.1f", protein_id, temp_match.alignment_score);
                    }
                    
                    // Get reference sequence
                    const char* ref_seq = &protein_db->sequences[protein_db->seq_offsets[protein_id]];
                    uint16_t ref_len = protein_db->seq_lengths[protein_id];
                    
                    if (temp_match.ref_start < ref_len) {
                        // Use the already extended region for Smith-Waterman
                        // Allow SW to see the full remaining sequence (up to 60 AA)
                        uint16_t available_ref = ref_len - temp_match.ref_start;
                        uint16_t available_query = frame.length - temp_match.query_start;
                        uint16_t sw_ref_len = min(available_ref, (uint16_t)60);
                        uint16_t sw_query_len = min(available_query, (uint16_t)60);
                        
                        uint16_t sw_query_start, sw_ref_start, sw_length;
                        uint8_t sw_num_mutations;
                        
                        // Debug: print sequences being aligned (first few chars)
                        if (tid == 0 && DEBUG_TRANS) {
                            DEBUG_PRINT("SW input: query_len=%d, ref_len=%d, protein_id=%d", 
                                    sw_query_len, sw_ref_len, protein_id);
                            DEBUG_PRINT("Query offset=%d: %.15s", temp_match.query_start, &frame.sequence[temp_match.query_start]);
                            DEBUG_PRINT("Ref offset=%d: %.15s", temp_match.ref_start, &ref_seq[temp_match.ref_start]);
                            DEBUG_PRINT("Ref full seq (first 20): %.20s", ref_seq);
                        }
                        
                        float sw_score = smith_waterman_align(
                            &frame.sequence[temp_match.query_start], sw_query_len,
                            &ref_seq[temp_match.ref_start], sw_ref_len,
                            &sw_query_start, &sw_ref_start, &sw_length,
                            temp_match.mutation_positions,
                            temp_match.ref_aas,
                            temp_match.query_aas,
                            temp_match.blosum_scores,
                            &sw_num_mutations
                        );
                        
                        // ADD DETAILED DEBUG OUTPUT HERE FOR MUTATION DETECTION
                        if (sw_score > temp_match.alignment_score) {
                            temp_match.alignment_score = sw_score;
                            temp_match.query_start += sw_query_start;
                            temp_match.ref_start += sw_ref_start;
                            temp_match.match_length = sw_length;
                            temp_match.num_mutations = sw_num_mutations;
                            temp_match.identity = (float)(sw_length - sw_num_mutations) / sw_length;
                            temp_match.used_smith_waterman = true;
                            
                            // DEBUG: Add detailed mutation analysis
                            if (tid < 5 && sw_num_mutations > 0) {  // First 5 reads with mutations
                                printf("\n[MUTATION DEBUG] Read %d, Frame %d, Gene %d:\n", 
                                    tid, frame.frame, temp_match.gene_id);
                                printf("  Alignment: Query pos %d-%d → Ref pos %d-%d\n",
                                    temp_match.query_start, temp_match.query_start + sw_length,
                                    temp_match.ref_start, temp_match.ref_start + sw_length);
                                printf("  Identity: %.1f%% (%d mutations in %d AA)\n",
                                    temp_match.identity * 100, sw_num_mutations, sw_length);
                                
                                // Check for QRDR positions
                                // gyrA: positions 83, 87
                                // parC: positions 80, 84
                                bool has_qrdr_mutation = false;
                                
                                for (int i = 0; i < sw_num_mutations; i++) {
                                    int mutation_pos_in_protein = temp_match.ref_start + temp_match.mutation_positions[i];
                                    
                                    printf("  Mutation %d: Position %d, %c → %c (BLOSUM: %.1f)\n",
                                        i + 1, mutation_pos_in_protein,
                                        temp_match.ref_aas[i], temp_match.query_aas[i],
                                        temp_match.blosum_scores[i]);
                                    
                                    // Check if this is a QRDR position
                                    if (temp_match.gene_id == 1) {  // gyrA (corrected gene_id)
                                        if (mutation_pos_in_protein == 82) {  // 0-based, so 83-1
                                            printf("    *** QRDR MUTATION: gyrA S83%c ***\n", temp_match.query_aas[i]);
                                            has_qrdr_mutation = true;
                                        } else if (mutation_pos_in_protein == 86) {  // 87-1
                                            printf("    *** QRDR MUTATION: gyrA D87%c ***\n", temp_match.query_aas[i]);
                                            has_qrdr_mutation = true;
                                        }
                                    } else if (temp_match.gene_id == 2) {  // parC (corrected gene_id)
                                        if (mutation_pos_in_protein == 79) {  // 80-1
                                            printf("    *** QRDR MUTATION: parC S80%c ***\n", temp_match.query_aas[i]);
                                            has_qrdr_mutation = true;
                                        } else if (mutation_pos_in_protein == 83) {  // 84-1
                                            printf("    *** QRDR MUTATION: parC E84%c ***\n", temp_match.query_aas[i]);
                                            has_qrdr_mutation = true;
                                        }
                                    } else if (temp_match.gene_id == 0) {  // grlA (S. aureus)
                                        if (mutation_pos_in_protein == 83) {  // S84 (0-based)
                                            printf("    *** QRDR MUTATION: grlA S84%c ***\n", temp_match.query_aas[i]);
                                            has_qrdr_mutation = true;
                                        } else if (mutation_pos_in_protein == 87) {  // E88
                                            printf("    *** QRDR MUTATION: grlA E88%c ***\n", temp_match.query_aas[i]);
                                            has_qrdr_mutation = true;
                                        }
                                    }
                                }
                                
                                if (!has_qrdr_mutation && sw_num_mutations > 0) {
                                    printf("  Note: Mutations found but NOT in QRDR positions\n");
                                }
                                
                                // Show the actual alignment region
                                printf("  Query peptide: ");
                                for (int j = 0; j < min(sw_length, 30); j++) {
                                    printf("%c", frame.sequence[temp_match.query_start + j]);
                                }
                                printf("\n");
                                
                                printf("  Ref peptide:   ");
                                for (int j = 0; j < min(sw_length, 30); j++) {
                                    printf("%c", ref_seq[temp_match.ref_start + j]);
                                }
                                printf("\n\n");
                            }
                        }
                    }
                }
                                
                // Only accept matches with sufficient identity
                if (temp_match.identity >= MIN_IDENTITY_THRESHOLD) {
                    read_matches[match_count] = temp_match;
                    match_count++;
                }
            }
        }
    }
    
    // Keep only the best scoring match (unless there are ties)
    if (match_count > 1) {
        float max_score = read_matches[0].alignment_score;
        int best_idx = 0;
        
        // Find the maximum score
        for (uint32_t i = 1; i < match_count; i++) {
            if (read_matches[i].alignment_score > max_score) {
                max_score = read_matches[i].alignment_score;
                best_idx = i;
            }
        }
        
        // Count how many matches have the same max score
        int ties = 0;
        for (uint32_t i = 0; i < match_count; i++) {
            if (read_matches[i].alignment_score == max_score) {
                ties++;
            }
        }
        
        // If there's a single best match, keep only that one
        if (ties == 1) {
            if (best_idx != 0) {
                read_matches[0] = read_matches[best_idx];
            }
            match_count = 1;
        } else {
            // Keep all matches with the max score
            uint32_t kept = 0;
            for (uint32_t i = 0; i < match_count; i++) {
                if (read_matches[i].alignment_score == max_score) {
                    if (kept != i) {
                        read_matches[kept] = read_matches[i];
                    }
                    kept++;
                }
            }
            match_count = kept;
        }
    }
    
    match_counts[tid] = match_count;
    
    if (tid == 0 && match_count > 0 && DEBUG_TRANS) {
        DEBUG_PRINT("Read 0: %d protein matches found (SW enabled: %s)", 
                   match_count, enable_smith_waterman ? "YES" : "NO");
    }
    
    // Add summary statistics at the end of the kernel
    if (tid == 0) {
        // Count how many reads had QRDR coverage
        __shared__ int qrdr_coverage_count;
        __shared__ int qrdr_mutation_count;
        
        if (threadIdx.x == 0) {
            qrdr_coverage_count = 0;
            qrdr_mutation_count = 0;
        }
        __syncthreads();
        
        // Each thread checks its matches
        for (uint32_t i = 0; i < match_count; i++) {
            if (read_matches[i].used_smith_waterman) {
                if (covers_qrdr_region(read_matches[i].gene_id, 
                                      read_matches[i].ref_start, 
                                      read_matches[i].match_length)) {
                    atomicAdd(&qrdr_coverage_count, 1);
                    
                    // Check if any mutations are in QRDR positions
                    for (int j = 0; j < read_matches[i].num_mutations; j++) {
                        int pos = read_matches[i].ref_start + read_matches[i].mutation_positions[j];
                        // Check for mutations at key QRDR positions (works for any gene)
                        // gyrA: S83 (pos 82) and D87 (pos 86)
                        // parC: S80 (pos 79) and E84 (pos 83)
                        if (pos == 82 || pos == 86 || pos == 79 || pos == 83) {
                            atomicAdd(&qrdr_mutation_count, 1);
                            break;
                        }
                    }
                }
            }
        }
        
        __syncthreads();
        if (threadIdx.x == 0) {
            printf("[SUMMARY] QRDR coverage: %d alignments, QRDR mutations: %d\n", 
                   qrdr_coverage_count, qrdr_mutation_count);
        }
    }
}

// Host wrapper class
class TranslatedSearchEngine {
private:
    ProteinDatabase* d_protein_db;
    TranslatedFrame* d_translated_frames;
    uint32_t* d_frame_counts;
    ProteinMatch* d_matches;
    uint32_t* d_match_counts;
    
    int max_batch_size;
    bool smith_waterman_enabled;
    
public:
    TranslatedSearchEngine(int batch_size = 10000, bool enable_sw = false) 
        : max_batch_size(batch_size), smith_waterman_enabled(enable_sw) {
        
        DEBUG_PRINT("Initializing TranslatedSearchEngine (batch=%d, SW=%s)", 
                   batch_size, enable_sw ? "enabled" : "disabled");
        
        cudaError_t err;
        err = cudaMalloc(&d_translated_frames, batch_size * NUM_FRAMES * sizeof(TranslatedFrame));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate d_translated_frames: %s", cudaGetErrorString(err));
        }
        err = cudaMalloc(&d_frame_counts, batch_size * sizeof(uint32_t));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate d_frame_counts: %s", cudaGetErrorString(err));
        }
        err = cudaMalloc(&d_matches, batch_size * 32 * sizeof(ProteinMatch));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate d_matches: %s", cudaGetErrorString(err));
        }
        err = cudaMalloc(&d_match_counts, batch_size * sizeof(uint32_t));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate d_match_counts: %s", cudaGetErrorString(err));
        }
        
        d_protein_db = nullptr;
    }
    
    ~TranslatedSearchEngine() {
        if (d_translated_frames) cudaFree(d_translated_frames);
        if (d_frame_counts) cudaFree(d_frame_counts);
        if (d_matches) cudaFree(d_matches);
        if (d_match_counts) cudaFree(d_match_counts);
        
        if (d_protein_db) {
            ProteinDatabase h_db;
            cudaMemcpy(&h_db, d_protein_db, sizeof(ProteinDatabase), cudaMemcpyDeviceToHost);
            
            if (h_db.sorted_kmer_hashes) cudaFree(h_db.sorted_kmer_hashes);
            if (h_db.kmer_start_indices) cudaFree(h_db.kmer_start_indices);
            if (h_db.kmer_counts) cudaFree(h_db.kmer_counts);
            if (h_db.position_data) cudaFree(h_db.position_data);
            if (h_db.protein_ids) cudaFree(h_db.protein_ids);
            if (h_db.gene_ids) cudaFree(h_db.gene_ids);
            if (h_db.species_ids) cudaFree(h_db.species_ids);
            if (h_db.seq_lengths) cudaFree(h_db.seq_lengths);
            if (h_db.sequences) cudaFree(h_db.sequences);
            if (h_db.seq_offsets) cudaFree(h_db.seq_offsets);
            
            cudaFree(d_protein_db);
        }
    }
    
    bool loadProteinDatabase(const std::string& db_path) {
        DEBUG_PRINT("Loading enhanced protein database from %s", db_path.c_str());
        
        if (d_protein_db) {
            DEBUG_PRINT("Protein database already loaded");
            return true;
        }
        
        // Load k-mer index
        std::string kmer_path = db_path + "/protein_kmers.bin";
        std::ifstream kmer_file(kmer_path, std::ios::binary);
        if (!kmer_file.good()) {
            DEBUG_PRINT("ERROR: Cannot read k-mer file: %s", kmer_path.c_str());
            return false;
        }
        
        uint32_t kmer_length, num_kmers;
        kmer_file.read(reinterpret_cast<char*>(&kmer_length), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&num_kmers), sizeof(uint32_t));
        
        if (kmer_length != PROTEIN_KMER_SIZE) {
            DEBUG_PRINT("ERROR: K-mer size mismatch: expected %d, got %d", PROTEIN_KMER_SIZE, kmer_length);
            return false;
        }
        
        DEBUG_PRINT("Loading %d protein 5-mers", num_kmers);
        
        // Load and sort k-mers
        std::map<uint64_t, std::vector<uint32_t>> kmer_map;
        
        for (uint32_t i = 0; i < num_kmers; i++) {
            char kmer_seq[6] = {0};
            kmer_file.read(kmer_seq, kmer_length);
            
            uint64_t hash = 0;
            const uint64_t prime = 31;
            for (int j = 0; j < kmer_length; j++) {
                hash = hash * prime + aa_to_index(kmer_seq[j]);
            }
            
            uint32_t num_positions;
            kmer_file.read(reinterpret_cast<char*>(&num_positions), sizeof(uint32_t));
            
            for (uint32_t j = 0; j < num_positions; j++) {
                uint32_t protein_idx, position;
                kmer_file.read(reinterpret_cast<char*>(&protein_idx), sizeof(uint32_t));
                kmer_file.read(reinterpret_cast<char*>(&position), sizeof(uint32_t));
                
                uint32_t encoded = (protein_idx << 16) | (position & 0xFFFF);
                kmer_map[hash].push_back(encoded);
            }
        }
        kmer_file.close();
        
        // Create sorted arrays
        std::vector<uint64_t> sorted_hashes;
        std::vector<uint32_t> start_indices;
        std::vector<uint32_t> kmer_counts;
        std::vector<uint32_t> position_data;
        
        for (const auto& pair : kmer_map) {
            sorted_hashes.push_back(pair.first);
            start_indices.push_back(position_data.size());
            kmer_counts.push_back(pair.second.size());
            
            for (uint32_t pos : pair.second) {
                position_data.push_back(pos);
            }
        }
        
        // Load protein sequences for Smith-Waterman
        std::string protein_path = db_path + "/proteins.bin";
        std::ifstream protein_file(protein_path, std::ios::binary);
        if (!protein_file.good()) {
            DEBUG_PRINT("ERROR: Cannot read protein file: %s", protein_path.c_str());
            return false;
        }
        
        uint32_t num_proteins;
        protein_file.read(reinterpret_cast<char*>(&num_proteins), sizeof(uint32_t));
        
        DEBUG_PRINT("Reading %d proteins from database", num_proteins);
        
        // Read all remaining data as one big sequence block
        protein_file.seekg(0, std::ios::end);
        size_t file_size = protein_file.tellg();
        protein_file.seekg(sizeof(uint32_t), std::ios::beg); // Skip num_proteins
        
        size_t remaining_size = file_size - sizeof(uint32_t);
        std::vector<char> all_sequences(remaining_size + 1);
        protein_file.read(all_sequences.data(), remaining_size);
        all_sequences[remaining_size] = '\0';
        
        DEBUG_PRINT("Read %zu bytes of sequence data", remaining_size);
        printf("[DATABASE DEBUG] Loading protein database from: %s\n", db_path.c_str());
        printf("[DATABASE DEBUG] Number of proteins: %d\n", num_proteins);
        printf("[DATABASE DEBUG] Total sequence bytes: %zu\n", remaining_size);
        
        // Load protein metadata from protein_details.json
        std::vector<uint32_t> protein_ids(num_proteins);
        std::vector<uint32_t> gene_ids(num_proteins);
        std::vector<uint32_t> species_ids(num_proteins);
        std::vector<uint16_t> seq_lengths(num_proteins);
        std::vector<uint32_t> seq_offsets(num_proteins);
        
        DEBUG_PRINT("Loading protein metadata for %d proteins", num_proteins);
        
        // Try to load metadata from protein_details.json
        std::string metadata_path = db_path + "/protein_details.json";
        std::ifstream metadata_file(metadata_path);
        
        if (metadata_file.good()) {
            // Parse JSON metadata
            std::string json_content((std::istreambuf_iterator<char>(metadata_file)),
                                   std::istreambuf_iterator<char>());
            metadata_file.close();
            
            // Simple JSON parsing for the array of protein objects
            size_t pos = 0;
            size_t protein_idx = 0;
            size_t current_offset = 0;
            
            while ((pos = json_content.find("\"id\":", pos)) != std::string::npos && protein_idx < num_proteins) {
                // Parse protein ID
                pos += 5;
                size_t id_start = json_content.find_first_of("0123456789", pos);
                size_t id_end = json_content.find_first_not_of("0123456789", id_start);
                protein_ids[protein_idx] = std::stoi(json_content.substr(id_start, id_end - id_start));
                
                // Parse gene_id
                size_t gene_pos = json_content.find("\"gene_id\":", pos);
                if (gene_pos != std::string::npos && gene_pos < pos + 500) {
                    gene_pos += 10;
                    size_t gene_start = json_content.find_first_of("0123456789", gene_pos);
                    size_t gene_end = json_content.find_first_not_of("0123456789", gene_start);
                    gene_ids[protein_idx] = std::stoi(json_content.substr(gene_start, gene_end - gene_start));
                }
                
                // Parse species_id
                size_t species_pos = json_content.find("\"species_id\":", pos);
                if (species_pos != std::string::npos && species_pos < pos + 500) {
                    species_pos += 13;
                    size_t species_start = json_content.find_first_of("0123456789", species_pos);
                    size_t species_end = json_content.find_first_not_of("0123456789", species_start);
                    species_ids[protein_idx] = std::stoi(json_content.substr(species_start, species_end - species_start));
                }
                
                // Parse length
                size_t length_pos = json_content.find("\"length\":", pos);
                if (length_pos != std::string::npos && length_pos < pos + 500) {
                    length_pos += 9;
                    size_t length_start = json_content.find_first_of("0123456789", length_pos);
                    size_t length_end = json_content.find_first_not_of("0123456789", length_start);
                    seq_lengths[protein_idx] = std::stoi(json_content.substr(length_start, length_end - length_start));
                }
                
                // Set offset
                seq_offsets[protein_idx] = current_offset;
                current_offset += seq_lengths[protein_idx];
                
                protein_idx++;
                pos = id_end;
            }
            
            DEBUG_PRINT("Loaded metadata for %d proteins from JSON", protein_idx);
            
            // Print mappings for debugging
            for (uint32_t i = 0; i < std::min((uint32_t)5, num_proteins); i++) {
                printf("[GENE MAPPING] Protein %d: gene_id=%d, species_id=%d, length=%d, offset=%d\n",
                       i, gene_ids[i], species_ids[i], seq_lengths[i], seq_offsets[i]);
            }
        } else {
            // Fallback: estimate from sequence data
            printf("[WARNING] Could not load protein_details.json, using estimated mappings\n");
            size_t avg_protein_len = remaining_size / num_proteins;
            for (uint32_t i = 0; i < num_proteins; i++) {
                protein_ids[i] = i;
                gene_ids[i] = 1;  // Default to gyrA (gene_id 1)
                species_ids[i] = i % 6;  // Cycle through species
                seq_offsets[i] = i * avg_protein_len;
                seq_lengths[i] = (i == num_proteins - 1) ? 
                                (remaining_size - seq_offsets[i]) : avg_protein_len;
            }
        }
        
        // Debug: Print first few proteins
        if (DEBUG_TRANS) {
            for (int i = 0; i < 3 && i < num_proteins; i++) {
                const char* seq_start = all_sequences.data() + seq_offsets[i];
                DEBUG_PRINT("Protein %d: offset=%d, len=%d, seq=%.15s...", 
                           i, seq_offsets[i], seq_lengths[i], seq_start);
            }
        }
        
        // Use the all_sequences vector as our sequences
        std::vector<char> sequences = std::move(all_sequences);
        protein_file.close();
        
        DEBUG_PRINT("Loaded %d proteins, total sequence length: %zu", num_proteins, sequences.size());
        
        // Allocate and copy to GPU
        ProteinDatabase h_db;
        h_db.num_proteins = num_proteins;
        h_db.num_kmers = sorted_hashes.size();
        
        cudaError_t err;
        
        // K-mer data
        err = cudaMalloc(&h_db.sorted_kmer_hashes, sorted_hashes.size() * sizeof(uint64_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.kmer_start_indices, start_indices.size() * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.kmer_counts, kmer_counts.size() * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.position_data, position_data.size() * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        
        // Protein metadata
        err = cudaMalloc(&h_db.protein_ids, num_proteins * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.gene_ids, num_proteins * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.species_ids, num_proteins * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.seq_lengths, num_proteins * sizeof(uint16_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.seq_offsets, num_proteins * sizeof(uint32_t));
        if (err != cudaSuccess) return false;
        err = cudaMalloc(&h_db.sequences, sequences.size() * sizeof(char));
        if (err != cudaSuccess) return false;
        
        // Copy data
        cudaMemcpy(h_db.sorted_kmer_hashes, sorted_hashes.data(), sorted_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.kmer_start_indices, start_indices.data(), start_indices.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.kmer_counts, kmer_counts.data(), kmer_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.position_data, position_data.data(), position_data.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.protein_ids, protein_ids.data(), num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.gene_ids, gene_ids.data(), num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.species_ids, species_ids.data(), num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.seq_lengths, seq_lengths.data(), num_proteins * sizeof(uint16_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.seq_offsets, seq_offsets.data(), num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(h_db.sequences, sequences.data(), sequences.size() * sizeof(char), cudaMemcpyHostToDevice);
        
        // Copy database structure
        err = cudaMalloc(&d_protein_db, sizeof(ProteinDatabase));
        if (err != cudaSuccess) return false;
        err = cudaMemcpy(d_protein_db, &h_db, sizeof(ProteinDatabase), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) return false;
        
        DEBUG_PRINT("Enhanced protein database loaded: %d proteins, %d unique 5-mers, SW=%s", 
                   num_proteins, (int)sorted_hashes.size(), smith_waterman_enabled ? "enabled" : "disabled");
        return true;
    }
    
    void setSmithWatermanEnabled(bool enabled) {
        smith_waterman_enabled = enabled;
        DEBUG_PRINT("Smith-Waterman alignment %s", enabled ? "ENABLED" : "DISABLED");
    }
    
    void searchTranslatedReads(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,
        int num_reads,
        ProteinMatch* results,
        uint32_t* result_counts
    ) {
        if (num_reads > max_batch_size) {
            DEBUG_PRINT("ERROR: num_reads (%d) exceeds max_batch_size (%d)", num_reads, max_batch_size);
            return;
        }
        
        if (!d_protein_db) {
            DEBUG_PRINT("ERROR: Protein database not loaded");
            return;
        }
        
        cudaError_t err;
        err = cudaMemset(d_frame_counts, 0, num_reads * sizeof(uint32_t));
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Failed to reset frame counts: %s", cudaGetErrorString(err));
            return;
        }
        err = cudaMemset(d_match_counts, 0, num_reads * sizeof(uint32_t));
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Failed to reset match counts: %s", cudaGetErrorString(err));
            return;
        }
        
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        // Stage 1: 6-frame translation
        six_frame_translate_kernel<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets,
            num_reads, d_translated_frames, d_frame_counts
        );
        
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Translation kernel failed: %s", cudaGetErrorString(err));
            return;
        }
        
        // Stage 2: Enhanced 5-mer protein matching with optional Smith-Waterman
        enhanced_protein_kmer_match_kernel<<<grid_size, block_size>>>(
            d_translated_frames, d_frame_counts,
            num_reads, d_protein_db,
            d_matches, d_match_counts, 32,
            smith_waterman_enabled
        );
        
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Enhanced protein matching kernel failed: %s", cudaGetErrorString(err));
            return;
        }
        
        // Copy results
        err = cudaMemcpy(results, d_matches, num_reads * 32 * sizeof(ProteinMatch), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Failed to copy protein matches: %s", cudaGetErrorString(err));
            return;
        }
        err = cudaMemcpy(result_counts, d_match_counts, num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Failed to copy match counts: %s", cudaGetErrorString(err));
            return;
        }
    }
};

// C interface for integration
extern "C" {
    void* create_translated_search_engine(int batch_size) {
        TranslatedSearchEngine* engine = new TranslatedSearchEngine(batch_size, false); // SW disabled by default
        return engine;
    }
    
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw) {
        TranslatedSearchEngine* engine = new TranslatedSearchEngine(batch_size, enable_sw);
        return engine;
    }
    
    void destroy_translated_search_engine(void* engine) {
        if (engine) {
            delete static_cast<TranslatedSearchEngine*>(engine);
        }
    }
    
    int load_protein_database(void* engine, const char* db_path) {
        if (!engine) return -1;
        TranslatedSearchEngine* te = static_cast<TranslatedSearchEngine*>(engine);
        return te->loadProteinDatabase(db_path) ? 0 : -1;
    }
    
    void set_smith_waterman_enabled(void* engine, bool enabled) {
        if (engine) {
            TranslatedSearchEngine* te = static_cast<TranslatedSearchEngine*>(engine);
            te->setSmithWatermanEnabled(enabled);
        }
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
        if (!engine) return -1;
        TranslatedSearchEngine* te = static_cast<TranslatedSearchEngine*>(engine);
        te->searchTranslatedReads(
            d_reads, d_read_lengths, d_read_offsets,
            d_reads_to_process, num_reads,
            static_cast<ProteinMatch*>(results), result_counts
        );
        return 0;
    }
}

#endif // TRANSLATED_SEARCH_CU