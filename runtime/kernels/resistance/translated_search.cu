// translated_search.cu
// GPU-accelerated 6-frame translation and protein search for resistance detection

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

namespace cg = cooperative_groups;

// Debug macros
#define DEBUG_TRANS 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG_TRANS) { printf("[TRANS DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

// Constants
#define CODON_SIZE 3
#define NUM_FRAMES 6
#define MAX_PROTEIN_LENGTH 200  // Max length for translated peptides
#define PROTEIN_KMER_SIZE 3     // K-mer size for proteins
#define MIN_PEPTIDE_LENGTH 20   // Minimum peptide length to consider
#define AA_ALPHABET_SIZE 24     // 20 standard + B, Z, X, *

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

// BLOSUM62 matrix for GPU (flattened)
__constant__ float BLOSUM62_GPU[24*24];  // To be loaded from file

// Amino acid to index mapping
__device__ inline int aa_to_index(char aa) {
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
    return BLOSUM62_GPU[idx1 * 24 + idx2];
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

// Structure for translated read
struct TranslatedFrame {
    char sequence[MAX_PROTEIN_LENGTH];
    uint16_t length;
    int8_t frame;  // -3, -2, -1, +1, +2, +3
    uint16_t start_pos;  // Start position in original read
};

// Structure for protein match
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
};

// Protein database structure
struct ProteinDatabase {
    uint32_t num_proteins;
    uint32_t* protein_ids;
    uint32_t* gene_ids;
    uint32_t* species_ids;
    char* sequences;         // Concatenated sequences
    uint32_t* seq_offsets;   // Start position of each sequence
    uint16_t* seq_lengths;   // Length of each sequence
    
    // K-mer index
    uint32_t num_kmers;
    uint64_t* kmer_hashes;   // Hash of each k-mer
    uint32_t* kmer_positions; // Protein ID and position encoded
};

// 6-frame translation kernel
__global__ void six_frame_translate_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const int num_reads,
    TranslatedFrame* translated_frames,  // 6 frames per read
    uint32_t* frame_counts              // Number of valid frames per read
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const char* read = reads + read_offsets[tid];
    const int read_len = read_lengths[tid];
    
    // Output base for this read's frames
    TranslatedFrame* read_frames = &translated_frames[tid * NUM_FRAMES];
    int valid_frames = 0;
    
    // Forward frames (+1, +2, +3)
    for (int frame = 0; frame < 3; frame++) {
        TranslatedFrame& tf = read_frames[valid_frames];
        tf.frame = frame + 1;
        tf.start_pos = frame;
        tf.length = 0;
        
        // Translate forward frame
        for (int pos = frame; pos + 2 < read_len; pos += 3) {
            if (tf.length >= MAX_PROTEIN_LENGTH - 1) break;
            
            char aa = translate_codon(&read[pos]);
            if (aa == '*') {
                // Stop codon - save peptide if long enough
                if (tf.length >= MIN_PEPTIDE_LENGTH) {
                    tf.sequence[tf.length] = '\0';
                    valid_frames++;
                    
                    // Start new peptide
                    if (valid_frames < NUM_FRAMES) {
                        tf = read_frames[valid_frames];
                        tf.frame = frame + 1;
                        tf.start_pos = pos + 3;
                        tf.length = 0;
                    }
                }
                else {
                    // Reset for next ORF
                    tf.length = 0;
                    tf.start_pos = pos + 3;
                }
            }
            else {
                tf.sequence[tf.length++] = aa;
            }
        }
        
        // Save final peptide if long enough
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
        
        // Translate reverse complement
        for (int pos = read_len - frame - 1; pos >= 2; pos -= 3) {
            if (tf.length >= MAX_PROTEIN_LENGTH - 1) break;
            
            // Get reverse complement codon
            char rc_codon[3];
            for (int i = 0; i < 3; i++) {
                char base = read[pos - i];
                // Complement
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
                // Stop codon
                if (tf.length >= MIN_PEPTIDE_LENGTH) {
                    tf.sequence[tf.length] = '\0';
                    valid_frames++;
                    
                    // Start new peptide
                    if (valid_frames < NUM_FRAMES) {
                        tf = read_frames[valid_frames];
                        tf.frame = -(frame + 1);
                        tf.start_pos = pos - 3;
                        tf.length = 0;
                    }
                }
                else {
                    tf.length = 0;
                    tf.start_pos = pos - 3;
                }
            }
            else {
                tf.sequence[tf.length++] = aa;
            }
        }
        
        // Save final peptide
        if (tf.length >= MIN_PEPTIDE_LENGTH) {
            tf.sequence[tf.length] = '\0';
            valid_frames++;
        }
    }
    
    frame_counts[tid] = valid_frames;
    
    // Debug first read
    if (tid == 0 && DEBUG_TRANS) {
        DEBUG_PRINT("Read 0: %d valid frames from %d bp", valid_frames, read_len);
        for (int i = 0; i < valid_frames; i++) {
            DEBUG_PRINT("  Frame %+d: %d aa, start=%d", 
                       read_frames[i].frame, read_frames[i].length, read_frames[i].start_pos);
        }
    }
}

// Protein k-mer matching kernel
__global__ void protein_kmer_match_kernel(
    const TranslatedFrame* translated_frames,
    const uint32_t* frame_counts,
    const int num_reads,
    const ProteinDatabase* protein_db,
    ProteinMatch* matches,
    uint32_t* match_counts,
    const uint32_t max_matches_per_read
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
        
        // Collect all k-mer hits for this frame, then cluster them by protein
        struct KmerHit {
            uint32_t protein_id;
            uint32_t query_pos;
            uint32_t ref_pos;
        };
        KmerHit frame_hits[200];  // Max hits per frame
        int num_hits = 0;
        
        // Extract k-mers and find all matches
        for (int pos = 0; pos + PROTEIN_KMER_SIZE <= frame.length && num_hits < 200; pos++) {
            // Hash the k-mer using same algorithm as loader
            uint64_t kmer_hash = 0;
            for (int i = 0; i < PROTEIN_KMER_SIZE; i++) {
                kmer_hash = kmer_hash * 31 + (uint8_t)frame.sequence[pos + i];
            }
            
            // Find first match for this k-mer (to avoid too many hits)
            for (uint32_t i = 0; i < protein_db->num_kmers && num_hits < 200; i++) {
                if (protein_db->kmer_hashes[i] == kmer_hash) {
                    uint32_t encoded = protein_db->kmer_positions[i];
                    uint32_t protein_id = encoded >> 16;
                    uint32_t ref_pos = encoded & 0xFFFF;
                    
                    frame_hits[num_hits].protein_id = protein_id;
                    frame_hits[num_hits].query_pos = pos;
                    frame_hits[num_hits].ref_pos = ref_pos;
                    num_hits++;
                    break;  // Only take first match per k-mer
                }
            }
        }
        
        // Cluster hits by protein and create consolidated matches
        for (int h = 0; h < num_hits && match_count < max_matches_per_read; h++) {
            uint32_t protein_id = frame_hits[h].protein_id;
            
            // Check if we already have a match for this protein in this frame
            bool found_existing = false;
            for (uint32_t m = 0; m < match_count; m++) {
                if (read_matches[m].protein_id == protein_id && read_matches[m].frame == frame.frame) {
                    // Update existing match - extend it
                    read_matches[m].match_length += PROTEIN_KMER_SIZE;
                    read_matches[m].alignment_score += 6.0f;  // Add score for additional k-mer
                    found_existing = true;
                    break;
                }
            }
            
            // Create new match if this is the first hit for this protein
            if (!found_existing) {
                ProteinMatch& match = read_matches[match_count];
                match.read_id = tid;
                match.frame = frame.frame;
                match.protein_id = protein_id;
                match.gene_id = 0;  // Not loaded yet
                match.species_id = 0;  // Not loaded yet
                match.query_start = frame_hits[h].query_pos;
                match.ref_start = frame_hits[h].ref_pos;
                match.match_length = PROTEIN_KMER_SIZE;
                match.identity = 1.0f;  // Exact k-mer match
                match.alignment_score = 6.0f;  // Base score for k-mer match
                match.num_mutations = 0;
                
                match_count++;
            }
        }
    }
    
    match_counts[tid] = match_count;
    
    if (tid == 0 && match_count > 0 && DEBUG_TRANS) {
        DEBUG_PRINT("Read 0: %d protein matches found", match_count);
    }
}

// Calculate resistance confidence scores
__global__ void calculate_resistance_confidence_kernel(
    const ProteinMatch* matches,
    const uint32_t* match_counts,
    const int num_reads,
    const float* known_resistance_positions,  // Pre-loaded resistance position data
    float* confidence_scores,
    uint32_t* resistance_flags
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    uint32_t num_matches = match_counts[tid];
    if (num_matches == 0) {
        confidence_scores[tid] = 0.0f;
        resistance_flags[tid] = 0;
        return;
    }
    
    const ProteinMatch* read_matches = &matches[tid * 32];  // max_matches_per_read
    
    float max_confidence = 0.0f;
    uint32_t resistance_found = 0;
    
    for (uint32_t i = 0; i < num_matches; i++) {
        const ProteinMatch& match = read_matches[i];
        
        // Check each mutation
        for (int m = 0; m < match.num_mutations; m++) {
            uint16_t pos = match.mutation_positions[m];
            
            // Check if this is a known resistance position
            // (Simplified - would need proper lookup)
            bool is_resistance_pos = false;
            
            // Special check for key positions
            if (match.gene_id == 0) {  // gyrA
                if (pos == 83 || pos == 87) is_resistance_pos = true;
            } else if (match.gene_id == 1) {  // parC
                if (pos == 80 || pos == 84) is_resistance_pos = true;
            }
            
            if (is_resistance_pos) {
                // Calculate confidence based on:
                // 1. BLOSUM score (normalized)
                // 2. Identity of surrounding region
                // 3. Match quality
                
                float blosum_factor = (match.blosum_scores[m] + 4.0f) / 8.0f;
                float identity_factor = match.identity;
                float length_factor = min(1.0f, match.match_length / 30.0f);
                
                float confidence = blosum_factor * identity_factor * length_factor;
                
                if (confidence > max_confidence) {
                    max_confidence = confidence;
                }
                
                resistance_found = 1;
            }
        }
    }
    
    confidence_scores[tid] = max_confidence;
    resistance_flags[tid] = resistance_found;
}

// Host wrapper class
class TranslatedSearchEngine {
private:
    ProteinDatabase* d_protein_db;
    float* d_blosum_matrix;
    TranslatedFrame* d_translated_frames;
    uint32_t* d_frame_counts;
    ProteinMatch* d_matches;
    uint32_t* d_match_counts;
    
    int max_batch_size;
    
public:
    TranslatedSearchEngine(int batch_size = 10000) : max_batch_size(batch_size) {
        // Allocate GPU memory with error checking
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
        d_blosum_matrix = nullptr;
    }
    
    ~TranslatedSearchEngine() {
        if (d_translated_frames) cudaFree(d_translated_frames);
        if (d_frame_counts) cudaFree(d_frame_counts);
        if (d_matches) cudaFree(d_matches);
        if (d_match_counts) cudaFree(d_match_counts);
        
        // Free protein database and its k-mer data
        if (d_protein_db) {
            // First copy the database structure to get pointers to k-mer data
            ProteinDatabase h_db;
            cudaMemcpy(&h_db, d_protein_db, sizeof(ProteinDatabase), cudaMemcpyDeviceToHost);
            
            // Free k-mer data if it exists
            if (h_db.kmer_hashes) cudaFree(h_db.kmer_hashes);
            if (h_db.kmer_positions) cudaFree(h_db.kmer_positions);
            
            // Free the database structure itself
            cudaFree(d_protein_db);
        }
        
        if (d_blosum_matrix) cudaFree(d_blosum_matrix);
    }
    
    bool loadProteinDatabase(const std::string& db_path) {
        // Load protein database from binary files
        DEBUG_PRINT("Loading protein database from %s", db_path.c_str());
        
        if (d_protein_db) {
            DEBUG_PRINT("Protein database already loaded");
            return true;
        }
        
        // Load metadata for validation
        std::string metadata_path = db_path + "/metadata.json";
        std::ifstream metadata_file(metadata_path);
        if (!metadata_file.good()) {
            DEBUG_PRINT("ERROR: Cannot read metadata file: %s", metadata_path.c_str());
            return false;
        }
        
        // Load k-mer index
        std::string kmer_path = db_path + "/protein_kmers.bin";
        std::ifstream kmer_file(kmer_path, std::ios::binary);
        if (!kmer_file.good()) {
            DEBUG_PRINT("ERROR: Cannot read k-mer file: %s", kmer_path.c_str());
            return false;
        }
        
        // Read k-mer header
        uint32_t kmer_length, num_kmers;
        kmer_file.read(reinterpret_cast<char*>(&kmer_length), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&num_kmers), sizeof(uint32_t));
        
        if (kmer_length != PROTEIN_KMER_SIZE) {
            DEBUG_PRINT("ERROR: K-mer size mismatch: expected %d, got %d", PROTEIN_KMER_SIZE, kmer_length);
            return false;
        }
        
        DEBUG_PRINT("Loading %d protein k-mers", num_kmers);
        
        // Allocate host memory for k-mer data
        std::vector<uint64_t> h_kmer_hashes(num_kmers);
        std::vector<uint32_t> h_kmer_positions;
        
        // Read k-mers and their positions
        for (uint32_t i = 0; i < num_kmers; i++) {
            // Read k-mer sequence (3 bytes)
            char kmer_seq[4] = {0};
            kmer_file.read(kmer_seq, kmer_length);
            
            // Convert to hash (simple hash for now)
            uint64_t hash = 0;
            for (int j = 0; j < kmer_length; j++) {
                hash = hash * 31 + (uint8_t)kmer_seq[j];
            }
            h_kmer_hashes[i] = hash;
            
            // Read number of positions
            uint32_t num_positions;
            kmer_file.read(reinterpret_cast<char*>(&num_positions), sizeof(uint32_t));
            
            // Read positions
            for (uint32_t j = 0; j < num_positions; j++) {
                uint32_t protein_idx, position;
                kmer_file.read(reinterpret_cast<char*>(&protein_idx), sizeof(uint32_t));
                kmer_file.read(reinterpret_cast<char*>(&position), sizeof(uint32_t));
                
                // Encode protein index and position into single uint32
                uint32_t encoded = (protein_idx << 16) | (position & 0xFFFF);
                h_kmer_positions.push_back(encoded);
            }
        }
        kmer_file.close();
        
        // Load protein sequences
        std::string protein_path = db_path + "/proteins.bin";
        std::ifstream protein_file(protein_path, std::ios::binary);
        if (!protein_file.good()) {
            DEBUG_PRINT("ERROR: Cannot read protein file: %s", protein_path.c_str());
            return false;
        }
        
        // Read protein header
        uint32_t num_proteins;
        protein_file.read(reinterpret_cast<char*>(&num_proteins), sizeof(uint32_t));
        
        DEBUG_PRINT("Loading %d proteins", num_proteins);
        
        // For now, create a minimal database with k-mer index only
        ProteinDatabase h_db;
        h_db.num_proteins = num_proteins;
        h_db.num_kmers = num_kmers;
        h_db.protein_ids = nullptr;
        h_db.gene_ids = nullptr;
        h_db.species_ids = nullptr;
        h_db.sequences = nullptr;
        h_db.seq_offsets = nullptr;
        h_db.seq_lengths = nullptr;
        
        // Allocate and copy k-mer data to GPU
        cudaError_t err;
        uint64_t* d_kmer_hashes;
        uint32_t* d_kmer_positions;
        
        err = cudaMalloc(&d_kmer_hashes, num_kmers * sizeof(uint64_t));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate k-mer hashes: %s", cudaGetErrorString(err));
            return false;
        }
        
        err = cudaMalloc(&d_kmer_positions, h_kmer_positions.size() * sizeof(uint32_t));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate k-mer positions: %s", cudaGetErrorString(err));
            cudaFree(d_kmer_hashes);
            return false;
        }
        
        err = cudaMemcpy(d_kmer_hashes, h_kmer_hashes.data(), num_kmers * sizeof(uint64_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to copy k-mer hashes: %s", cudaGetErrorString(err));
            cudaFree(d_kmer_hashes);
            cudaFree(d_kmer_positions);
            return false;
        }
        
        err = cudaMemcpy(d_kmer_positions, h_kmer_positions.data(), h_kmer_positions.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to copy k-mer positions: %s", cudaGetErrorString(err));
            cudaFree(d_kmer_hashes);
            cudaFree(d_kmer_positions);
            return false;
        }
        
        h_db.kmer_hashes = d_kmer_hashes;
        h_db.kmer_positions = d_kmer_positions;
        
        // Allocate and copy database structure to GPU
        err = cudaMalloc(&d_protein_db, sizeof(ProteinDatabase));
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to allocate protein database: %s", cudaGetErrorString(err));
            cudaFree(d_kmer_hashes);
            cudaFree(d_kmer_positions);
            return false;
        }
        
        err = cudaMemcpy(d_protein_db, &h_db, sizeof(ProteinDatabase), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            DEBUG_PRINT("Failed to copy protein database: %s", cudaGetErrorString(err));
            cudaFree(d_protein_db);
            cudaFree(d_kmer_hashes);
            cudaFree(d_kmer_positions);
            d_protein_db = nullptr;
            return false;
        }
        
        DEBUG_PRINT("Protein database loaded successfully: %d proteins, %d k-mers", num_proteins, num_kmers);
        return true;
    }
    
    void searchTranslatedReads(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,  // Which reads passed earlier filters
        int num_reads,
        ProteinMatch* results,
        uint32_t* result_counts
    ) {
        // Bounds checking
        if (num_reads > max_batch_size) {
            DEBUG_PRINT("ERROR: num_reads (%d) exceeds max_batch_size (%d)", num_reads, max_batch_size);
            return;
        }
        
        if (!d_protein_db) {
            DEBUG_PRINT("ERROR: Protein database not loaded");
            return;
        }
        
        // Reset counts with error checking
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
        
        // Stage 1: 6-frame translation
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        six_frame_translate_kernel<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets,
            num_reads, d_translated_frames, d_frame_counts
        );
        
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Translation kernel failed: %s", cudaGetErrorString(err));
            return;
        }
        
        // Stage 2: Protein k-mer matching
        protein_kmer_match_kernel<<<grid_size, block_size>>>(
            d_translated_frames, d_frame_counts,
            num_reads, d_protein_db,
            d_matches, d_match_counts, 32
        );
        
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            DEBUG_PRINT("ERROR: Protein matching kernel failed: %s", cudaGetErrorString(err));
            return;
        }
        
        // Copy results with error checking
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
        TranslatedSearchEngine* engine = new TranslatedSearchEngine(batch_size);
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