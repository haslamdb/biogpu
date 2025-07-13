// amr_detection_kernels.cu
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <cub/cub.cuh>
#include "amr_detection_pipeline.h"
#include "ncbi_amr_database_loader.h"

namespace cg = cooperative_groups;

// Constants
__constant__ char GENETIC_CODE[64];  // Will be initialized with codon table
__constant__ int BLOOM_HASH_FUNCS = 3;
constexpr int MAX_MATCHES_PER_READ = 32;  // Must match the value in amr_detection_pipeline.h

// Hash functions for bloom filter and minimizers
__device__ __forceinline__ uint64_t murmur_hash64(uint64_t key, uint64_t seed) {
    key ^= seed;
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

// Device-compatible strlen
__device__ __forceinline__ int device_strlen(const char* str) {
    int len = 0;
    while (str[len] != '\0') len++;
    return len;
}

// Convert DNA character to 2-bit encoding
__device__ __forceinline__ uint8_t dna_to_bits(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 'u': case 'U': return 3;
        default: return 255;  // N or invalid
    }
}

// Reverse complement of 2-bit encoded DNA
__device__ __forceinline__ uint64_t reverse_complement_kmer(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        uint8_t base = (kmer >> (2 * i)) & 3;
        uint8_t comp = 3 - base;  // A<->T, C<->G
        rc = (rc << 2) | comp;
    }
    return rc;
}

// Build bloom filter from AMR nucleotide sequences
__global__ void build_bloom_filter_kernel(
    const char* amr_sequences,
    const uint32_t* sequence_offsets,
    const uint32_t* sequence_lengths,
    uint64_t* bloom_filter,
    const uint32_t num_sequences,
    const uint32_t bloom_size,
    const int k
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x * blockDim.x;
    
    for (int seq_idx = tid; seq_idx < num_sequences; seq_idx += stride) {
        const uint32_t offset = sequence_offsets[seq_idx];
        const uint32_t length = sequence_lengths[seq_idx];
        const char* sequence = &amr_sequences[offset];
        
        // Process all k-mers in this sequence
        if (length >= k) {
            uint64_t kmer = 0;
            uint64_t mask = (1ULL << (2 * k)) - 1;
            int valid_bases = 0;
            
            for (uint32_t i = 0; i < length; i++) {
                uint8_t base = dna_to_bits(sequence[i]);
                
                if (base < 4) {
                    kmer = ((kmer << 2) | base) & mask;
                    valid_bases++;
                    
                    if (valid_bases >= k) {
                        // Add k-mer and its reverse complement to bloom filter
                        uint64_t rc_kmer = reverse_complement_kmer(kmer, k);
                        uint64_t canonical = min(kmer, rc_kmer);
                        
                        // Add to bloom filter with multiple hash functions
                        for (int h = 0; h < BLOOM_HASH_FUNCS; h++) {
                            uint64_t hash = murmur_hash64(canonical, h);
                            uint64_t bit_pos = hash % (bloom_size * 64);
                            uint64_t word_idx = bit_pos / 64;
                            uint64_t bit_idx = bit_pos % 64;
                            atomicOr((unsigned long long*)&bloom_filter[word_idx], 1ULL << bit_idx);
                        }
                    }
                } else {
                    valid_bases = 0;
                    kmer = 0;
                }
            }
        }
    }
}

// Generate minimizers from reads
__global__ void generate_minimizers_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    Minimizer* minimizers,
    uint32_t* minimizer_counts,
    uint32_t* minimizer_offsets,
    const int num_reads,
    const int k,  // minimizer k-mer size
    const int w   // window size
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const int offset = read_offsets[tid];
    const int length = read_lengths[tid];
    const char* read = &reads[offset];
    
    const int max_minimizers_per_read = 1000;  // Adjust as needed
    const int output_offset = minimizer_offsets[tid];
    
    uint32_t num_minimizers = 0;
    
    if (length >= k) {
        // Sliding window to find minimizers
        for (int window_start = 0; window_start <= length - k - w + 1; window_start++) {
            uint64_t min_hash = UINT64_MAX;
            int min_pos = -1;
            bool min_is_reverse = false;
            
            // Find minimum k-mer in window
            for (int pos = window_start; pos < window_start + w && pos <= length - k; pos++) {
                uint64_t kmer = 0;
                uint64_t mask = (1ULL << (2 * k)) - 1;
                bool valid = true;
                
                // Build k-mer
                for (int i = 0; i < k; i++) {
                    uint8_t base = dna_to_bits(read[pos + i]);
                    if (base >= 4) {
                        valid = false;
                        break;
                    }
                    kmer = ((kmer << 2) | base) & mask;
                }
                
                if (valid) {
                    // Get canonical k-mer
                    uint64_t rc_kmer = reverse_complement_kmer(kmer, k);
                    bool is_reverse = rc_kmer < kmer;
                    uint64_t canonical = is_reverse ? rc_kmer : kmer;
                    
                    // Hash the k-mer
                    uint64_t hash = murmur_hash64(canonical, 0);
                    
                    if (hash < min_hash) {
                        min_hash = hash;
                        min_pos = pos;
                        min_is_reverse = is_reverse;
                    }
                }
            }
            
            // Store minimizer if found
            if (min_pos >= 0 && num_minimizers < max_minimizers_per_read) {
                minimizers[output_offset + num_minimizers] = {
                    min_hash,
                    static_cast<uint32_t>(min_pos),
                    min_is_reverse
                };
                num_minimizers++;
            }
        }
    }
    
    minimizer_counts[tid] = num_minimizers;
}

// Screen minimizers against bloom filter
__global__ void screen_minimizers_kernel(
    const Minimizer* minimizers,
    const uint32_t* minimizer_counts,
    const uint32_t* minimizer_offsets,
    const uint64_t* bloom_filter,
    bool* read_passes_filter,
    const int num_reads,
    const uint64_t bloom_size
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const uint32_t count = minimizer_counts[tid];
    const uint32_t offset = minimizer_offsets[tid];
    
    bool passes = false;
    
    // Check if any minimizer is in bloom filter
    for (uint32_t i = 0; i < count && !passes; i++) {
        uint64_t hash = minimizers[offset + i].hash;
        
        bool in_bloom = true;
        for (int h = 0; h < BLOOM_HASH_FUNCS; h++) {
            uint64_t hash_val = murmur_hash64(hash, h);
            uint64_t bit_pos = hash_val % (bloom_size * 64);
            uint64_t word_idx = bit_pos / 64;
            uint64_t bit_idx = bit_pos % 64;
            
            if (!(bloom_filter[word_idx] & (1ULL << bit_idx))) {
                in_bloom = false;
                break;
            }
        }
        
        if (in_bloom) {
            passes = true;
        }
    }
    
    read_passes_filter[tid] = passes;
}

// Banded Smith-Waterman alignment kernel
__device__ int banded_smith_waterman(
    const char* query,
    const int query_len,
    const char* reference,
    const int ref_len,
    const int band_width,
    int& ref_start,
    int& ref_end,
    int& query_start,
    int& query_end
) {
    // Simplified banded SW implementation
    // In practice, you'd want a more optimized version
    
    const int MATCH = 2;
    const int MISMATCH = -1;
    const int GAP = -1;
    
    extern __shared__ int shared_mem[];
    int* scores = shared_mem;  // Dynamic shared memory
    
    int max_score = 0;
    int max_i = 0, max_j = 0;
    
    // Initialize first row and column
    for (int i = 0; i <= band_width && i <= query_len; i++) {
        scores[i * (2 * band_width + 1) + band_width] = 0;
    }
    
    // Fill scoring matrix with banding
    for (int i = 1; i <= query_len; i++) {
        int j_start = max(1, i - band_width);
        int j_end = min(ref_len, i + band_width);
        
        for (int j = j_start; j <= j_end; j++) {
            int diag_offset = j - i + band_width;
            
            if (diag_offset < 0 || diag_offset > 2 * band_width) continue;
            
            int match = (query[i-1] == reference[j-1]) ? MATCH : MISMATCH;
            
            int diag_score = scores[(i-1) * (2 * band_width + 1) + diag_offset] + match;
            int up_score = (diag_offset > 0) ? 
                scores[(i-1) * (2 * band_width + 1) + diag_offset - 1] + GAP : 0;
            int left_score = (diag_offset < 2 * band_width) ? 
                scores[(i-1) * (2 * band_width + 1) + diag_offset + 1] + GAP : 0;
            
            int score = max(0, max(diag_score, max(up_score, left_score)));
            scores[i * (2 * band_width + 1) + diag_offset] = score;
            
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback to find alignment boundaries
    ref_end = max_j;
    query_end = max_i;
    
    // Simplified traceback - in practice would be more complex
    ref_start = max(0, max_j - query_len);
    query_start = max(0, max_i - ref_len);
    
    return max_score;
}

// Main translated alignment kernel
__global__ void translated_alignment_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    const bool* read_passes_filter,
    const char* amr_proteins,
    const uint32_t* protein_offsets,
    const uint32_t* protein_lengths,
    const AMRGeneEntry* gene_entries,
    AMRHit* hits,
    uint32_t* hit_counts,
    const int num_reads,
    const int num_proteins,
    const AMRDetectionConfig config
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    // Skip reads that didn't pass bloom filter
    if (!read_passes_filter[tid]) {
        hit_counts[tid] = 0;
        return;
    }
    
    const int offset = read_offsets[tid];
    const int length = read_lengths[tid];
    const char* read = &reads[offset];
    
    // Allocate space for 6-frame translation in shared memory
    extern __shared__ char shared_translations[];
    const int max_protein_len = (length / 3) + 1;
    char* translations = &shared_translations[threadIdx.x * max_protein_len * 6];
    
    // Perform 6-frame translation
    for (int frame = 0; frame < 6; frame++) {
        char* frame_translation = &translations[frame * max_protein_len];
        int aa_pos = 0;
        
        bool is_reverse = frame >= 3;
        int start_pos = frame % 3;
        
        if (!is_reverse) {
            // Forward frames
            for (int pos = start_pos; pos + 2 < length; pos += 3) {
                uint8_t codon_idx = 0;
                bool valid = true;
                
                for (int i = 0; i < 3; i++) {
                    uint8_t base = dna_to_bits(read[pos + i]);
                    if (base >= 4) {
                        valid = false;
                        break;
                    }
                    codon_idx = (codon_idx << 2) | base;
                }
                
                frame_translation[aa_pos++] = valid ? GENETIC_CODE[codon_idx] : 'X';
            }
        } else {
            // Reverse frames - implement reverse complement translation
            // ... (similar logic but with RC)
        }
        
        frame_translation[aa_pos] = '\0';
    }
    
    // Search against AMR proteins
    uint32_t num_hits = 0;
    const int max_hits_per_read = 10;
    
    for (int prot_idx = 0; prot_idx < num_proteins && num_hits < max_hits_per_read; prot_idx++) {
        const uint32_t prot_offset = protein_offsets[prot_idx];
        const uint32_t prot_length = protein_lengths[prot_idx];
        const char* protein = &amr_proteins[prot_offset];
        const AMRGeneEntry& gene = gene_entries[prot_idx];
        
        // Try each translation frame
        for (int frame = 0; frame < 6; frame++) {
            char* frame_translation = &translations[frame * max_protein_len];
            int trans_len = device_strlen(frame_translation);
            
            if (trans_len < config.min_alignment_length) continue;
            
            // Perform banded alignment
            int ref_start, ref_end, query_start, query_end;
            int score = banded_smith_waterman(
                frame_translation, trans_len,
                protein, prot_length,
                config.band_width,
                ref_start, ref_end, query_start, query_end
            );
            
            int alignment_length = ref_end - ref_start;
            if (alignment_length < config.min_alignment_length) continue;
            
            // Calculate identity
            int matches = 0;
            for (int i = 0; i < alignment_length; i++) {
                if (frame_translation[query_start + i] == protein[ref_start + i]) {
                    matches++;
                }
            }
            
            float identity = (float)matches / alignment_length;
            float coverage = (float)alignment_length / prot_length;
            
            if (identity >= config.min_identity && coverage >= config.min_coverage) {
                // Found a hit!
                AMRHit& hit = hits[tid * max_hits_per_read + num_hits];
                hit.read_id = tid;
                hit.gene_id = prot_idx;
                hit.ref_start = ref_start;
                hit.ref_end = ref_end;
                hit.read_start = query_start;
                hit.read_end = query_end;
                hit.identity = identity;
                hit.coverage = coverage;
                hit.frame = (frame < 3) ? frame + 1 : -(frame - 2);
                hit.is_complete_gene = (coverage >= 0.95f);
                
                // Copy gene metadata
                memcpy(hit.gene_name, gene.gene_name, 64);
                memcpy(hit.drug_class, gene.class_, 32);
                
                num_hits++;
                break;  // Only take best frame per protein
            }
        }
    }
    
    hit_counts[tid] = num_hits;
}

// Extend alignments using minimizer information
__global__ void extend_alignments_kernel(
    const char* reads,
    const int* read_offsets,
    const int* read_lengths,
    const Minimizer* minimizers,
    const uint32_t* minimizer_counts,
    const uint32_t* minimizer_offsets,
    AMRHit* hits,
    const uint32_t* hit_counts,
    const char* amr_proteins,
    const uint32_t* protein_offsets,
    const uint32_t* protein_lengths,
    const int num_reads
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const uint32_t num_hits = hit_counts[tid];
    if (num_hits == 0) return;
    
    const uint32_t num_minimizers = minimizer_counts[tid];
    const uint32_t minimizer_offset = minimizer_offsets[tid];
    
    // For each hit, try to extend using nearby minimizers
    for (uint32_t h = 0; h < num_hits; h++) {
        AMRHit& hit = hits[tid * 10 + h];
        
        // Find minimizers near the alignment boundaries
        for (uint32_t m = 0; m < num_minimizers; m++) {
            const Minimizer& minimizer = minimizers[minimizer_offset + m];
            
            // Check if minimizer is near alignment boundary
            int frame_adjusted_pos = minimizer.pos;
            if (hit.frame < 0) {
                frame_adjusted_pos = read_lengths[tid] - minimizer.pos - 15;  // k-mer size
            }
            
            int translated_pos = frame_adjusted_pos / 3;
            
            // If minimizer is within extension range
            if (abs(translated_pos - hit.read_end) < 20 ||
                abs(translated_pos - hit.read_start) < 20) {
                
                // Attempt to extend alignment
                // This is a simplified version - you'd want more sophisticated extension
                if (translated_pos > hit.read_end && hit.ref_end < protein_lengths[hit.gene_id]) {
                    // Try to extend to the right
                    int extension_length = min(10, (int)protein_lengths[hit.gene_id] - hit.ref_end);
                    hit.ref_end += extension_length;
                    hit.read_end += extension_length;
                    hit.coverage = (float)(hit.ref_end - hit.ref_start) / protein_lengths[hit.gene_id];
                }
            }
        }
    }
}

// Update coverage statistics
__global__ void update_coverage_stats_kernel(
    const AMRHit* hits,
    const uint32_t* hit_counts,
    AMRCoverageStats* coverage_stats,
    const int num_reads,
    const int num_genes
) {
    // Use a grid-stride loop to process all genes
    for (int gene_idx = blockIdx.x * blockDim.x + threadIdx.x;
         gene_idx < num_genes;
         gene_idx += gridDim.x * blockDim.x)
    {
        AMRCoverageStats& stats = coverage_stats[gene_idx];

        // Iterate through all reads to find hits for this gene
        for (int read_idx = 0; read_idx < num_reads; read_idx++) {
            const uint32_t num_hits = hit_counts[read_idx];
            for (uint32_t h = 0; h < num_hits; h++) {
                const AMRHit& hit = hits[read_idx * MAX_MATCHES_PER_READ + h]; // Assumes MAX_MATCHES_PER_READ is known or passed

                if (hit.gene_id == gene_idx) {
                    // Atomically update the aggregate statistics
                    atomicAdd(&stats.total_reads, 1);
                    atomicAdd(&stats.total_bases_mapped, hit.ref_end - hit.ref_start);

                    // Atomically update the per-position coverage counts
                    for (uint16_t pos = hit.ref_start; pos < hit.ref_end && pos < stats.gene_length; pos++) {
                        atomicAdd(&stats.position_counts[pos], 1);
                    }
                }
            }
        }
    }
}

// Calculate final coverage statistics
__global__ void finalize_coverage_stats_kernel(
    AMRCoverageStats* coverage_stats,
    const AMRGeneEntry* gene_entries,
    const int num_genes
) {
    const int gene_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (gene_idx >= num_genes) return;

    AMRCoverageStats& stats = coverage_stats[gene_idx];
    if (stats.total_reads == 0) return; // Skip genes with no reads

    const uint16_t gene_length = gene_entries[gene_idx].protein_length;
    stats.gene_length = gene_length;

    // Calculate total covered positions
    uint32_t total_coverage = 0;
    uint16_t covered_positions = 0;
    for (uint16_t pos = 0; pos < gene_length; pos++) {
        uint32_t depth = stats.position_counts[pos];
        if (depth > 0) {
            covered_positions++;
            total_coverage += depth;
        }
    }
    stats.covered_positions = covered_positions;

    // Calculate final metrics
    if (covered_positions > 0) {
        stats.percent_coverage = (float)covered_positions / gene_length * 100.0f;
        stats.mean_coverage = (float)total_coverage / gene_length; // Avg coverage over all positions
        stats.mean_depth = (float)total_coverage / covered_positions; // Avg depth over covered positions

        // Calculate uniformity (coefficient of variation)
        float sum_squared_diff = 0.0f;
        for (uint16_t pos = 0; pos < gene_length; pos++) {
            float diff = (float)stats.position_counts[pos] - stats.mean_coverage;
            sum_squared_diff += diff * diff;
        }
        float variance = sum_squared_diff / gene_length;
        float std_dev = sqrtf(variance);
        stats.coverage_uniformity = 1.0f - (std_dev / (stats.mean_coverage + 1e-6f)); // Add epsilon for stability
    } else {
        stats.percent_coverage = 0.0f;
        stats.mean_coverage = 0.0f;
        stats.mean_depth = 0.0f;
        stats.coverage_uniformity = 0.0f;
    }
}
