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
#define DEBUG_TRANS 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG_TRANS) { printf("[TRANS DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

// Enhanced constants for 5-mer approach
#define CODON_SIZE 3
#define NUM_FRAMES 6
#define MAX_PROTEIN_LENGTH 200
#define PROTEIN_KMER_SIZE 5         // Increased from 3 to 5 amino acids
#define MIN_PEPTIDE_LENGTH 20       // Minimum peptide length to consider
#define MIN_SEED_HITS 2            // Require multiple k-mer hits for extension
#define EXTENSION_THRESHOLD 15     // Minimum amino acids for valid match
#define MIN_IDENTITY_THRESHOLD 0.9f // Minimum 90% identity for valid protein match
#define SW_SCORE_THRESHOLD 60.0f   // Threshold for Smith-Waterman alignment
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
    // Simplified banded Smith-Waterman for GPU
    const int GAP_OPEN = -10;
    const int GAP_EXTEND = -1;
    const int BAND_WIDTH = 10;
    
    // Use local arrays for small alignments
    float H[50][50] = {0}; // Score matrix (limited size for GPU)
    
    if (query_len > 49 || ref_len > 49) {
        // Fall back to simple scoring for large sequences
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
    
    // Fill scoring matrix
    for (int i = 1; i <= query_len; i++) {
        for (int j = 1; j <= ref_len; j++) {
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
    
    // Traceback to find alignment
    int i = max_i, j = max_j;
    int align_len = 0;
    uint8_t mut_count = 0;
    
    while (i > 0 && j > 0 && H[i][j] > 0 && align_len < 49) {
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
            
            // Calculate match span
            uint32_t min_query = protein_seeds[0].query_pos;
            uint32_t max_query = protein_seeds[0].query_pos + PROTEIN_KMER_SIZE;
            uint32_t min_ref = protein_seeds[0].ref_pos;
            uint32_t max_ref = protein_seeds[0].ref_pos + PROTEIN_KMER_SIZE;
            
            for (int i = 1; i < seed_count; i++) {
                min_query = min(min_query, protein_seeds[i].query_pos);
                max_query = max(max_query, protein_seeds[i].query_pos + PROTEIN_KMER_SIZE);
                min_ref = min(min_ref, protein_seeds[i].ref_pos);
                max_ref = max(max_ref, protein_seeds[i].ref_pos + PROTEIN_KMER_SIZE);
            }
            
            uint32_t query_span = max_query - min_query;
            uint32_t ref_span = max_ref - min_ref;
            
            // Quality filter: reasonable span and length
            if (query_span >= EXTENSION_THRESHOLD && 
                query_span <= frame.length && 
                abs((int)query_span - (int)ref_span) <= 5) {
                
                // Create temporary match to check identity
                ProteinMatch temp_match;
                temp_match.read_id = tid;
                temp_match.frame = frame.frame;
                temp_match.protein_id = protein_id;
                temp_match.gene_id = protein_db->gene_ids[protein_id];
                temp_match.species_id = protein_db->species_ids[protein_id];
                temp_match.query_start = min_query;
                temp_match.ref_start = min_ref;
                temp_match.match_length = query_span;
                temp_match.identity = (float)(seed_count * PROTEIN_KMER_SIZE) / query_span;
                temp_match.alignment_score = seed_count * 10.0f;
                temp_match.num_mutations = 0;
                temp_match.used_smith_waterman = false;
                
                // Optional Smith-Waterman for high-scoring matches
                if (enable_smith_waterman && temp_match.alignment_score >= SW_SCORE_THRESHOLD) {
                    // Get reference sequence
                    const char* ref_seq = &protein_db->sequences[protein_db->seq_offsets[protein_id]];
                    uint16_t ref_len = protein_db->seq_lengths[protein_id];
                    
                    if (temp_match.ref_start < ref_len) {
                        uint16_t available_ref = ref_len - temp_match.ref_start;
                        uint16_t available_query = frame.length - temp_match.query_start;
                        uint16_t sw_ref_len = min(available_ref, (uint16_t)50);
                        uint16_t sw_query_len = min(available_query, (uint16_t)50);
                        
                        uint16_t sw_query_start, sw_ref_start, sw_length;
                        uint8_t sw_num_mutations;
                        
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
                        
                        if (sw_score > temp_match.alignment_score) {
                            temp_match.alignment_score = sw_score;
                            temp_match.query_start += sw_query_start;
                            temp_match.ref_start += sw_ref_start;
                            temp_match.match_length = sw_length;
                            temp_match.num_mutations = sw_num_mutations;
                            temp_match.identity = (float)(sw_length - sw_num_mutations) / sw_length;
                            temp_match.used_smith_waterman = true;
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
    
    match_counts[tid] = match_count;
    
    if (tid == 0 && match_count > 0 && DEBUG_TRANS) {
        DEBUG_PRINT("Read 0: %d protein matches found (SW enabled: %s)", 
                   match_count, enable_smith_waterman ? "YES" : "NO");
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
        
        std::vector<uint32_t> protein_ids(num_proteins);
        std::vector<uint32_t> gene_ids(num_proteins);
        std::vector<uint32_t> species_ids(num_proteins);
        std::vector<uint16_t> seq_lengths(num_proteins);
        std::vector<uint32_t> seq_offsets(num_proteins);
        std::vector<char> sequences;
        
        for (uint32_t i = 0; i < num_proteins; i++) {
            protein_file.read(reinterpret_cast<char*>(&protein_ids[i]), sizeof(uint32_t));
            protein_file.read(reinterpret_cast<char*>(&gene_ids[i]), sizeof(uint32_t));
            protein_file.read(reinterpret_cast<char*>(&species_ids[i]), sizeof(uint32_t));
            
            uint16_t seq_len;
            protein_file.read(reinterpret_cast<char*>(&seq_len), sizeof(uint16_t));
            seq_lengths[i] = seq_len;
            seq_offsets[i] = sequences.size();
            
            std::vector<char> seq(seq_len + 1);
            protein_file.read(seq.data(), seq_len);
            seq[seq_len] = '\0';
            
            sequences.insert(sequences.end(), seq.begin(), seq.end());
        }
        protein_file.close();
        
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