// enhanced_mutation_detection.cu
// Improved mutation detection that compares against wild-type references
// regardless of what's in the protein database

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdint.h>
#include <cstdio>

#ifdef __NVCC__
#ifndef __CUDA_ARCH__
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>
#include <iomanip>
#endif
#endif

// Wild-type reference sequences for QRDR regions
__constant__ char WILDTYPE_GYRA_QRDR[54] = "TGCCTATGACATCGCAGTCGAGCGAACCTTATTGAAGATGTTTCCAAGCTGAA"; // Example gyrA QRDR
__constant__ char WILDTYPE_PARC_QRDR[55] = "TTCAGCACAGAGCGCTTGCTGGAAGATGTTTCTAAGCTGAAAGCGTTGGATCGC"; // Example parC QRDR

// Known resistance positions and their wild-type amino acids
struct ResistancePosition {
    uint32_t gene_id;
    uint16_t position;      // 1-based amino acid position
    char wildtype_aa;
    char resistant_aas[8];  // Common resistant amino acids
    int num_resistant;
    float resistance_score;
};

__constant__ ResistancePosition RESISTANCE_POSITIONS[16] = {
    // gyrA
    {0, 83, 'S', {'L', 'F', 'W', 'A'}, 4, 0.9f},
    {0, 87, 'D', {'N', 'G', 'Y', 'H'}, 4, 0.7f},
    
    // parC
    {1, 80, 'S', {'I', 'R', 'F'}, 3, 0.6f},
    {1, 84, 'E', {'V', 'K', 'G', 'A'}, 4, 0.4f},
    
    // gyrB
    {2, 426, 'G', {'D', 'N'}, 2, 0.5f},
    {2, 447, 'A', {'V', 'T'}, 2, 0.3f},
    
    // parE
    {3, 416, 'G', {'D', 'S'}, 2, 0.4f},
    {3, 420, 'A', {'V', 'T'}, 2, 0.3f},
    
    // Terminator
    {0, 0, 0, {0}, 0, 0.0f}
};

// Enhanced protein match structure with detailed mutation info
struct EnhancedProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    
    // Enhanced mutation details
    uint8_t num_resistance_mutations;
    uint8_t num_total_mutations;
    uint16_t resistance_positions[10];
    char wildtype_aas[10];
    char observed_aas[10];
    float resistance_scores[10];
    bool is_known_resistance[10];
    
    // Alignment quality
    float blosum_total_score;
    bool used_smith_waterman;
    bool covers_qrdr;
    
    // Translated sequences (for debugging)
    char query_sequence[100];
    char ref_sequence[100];
};

// Genetic code for translation
static __constant__ char GENETIC_CODE[64] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S',
    'W', 'C', '*', 'C', 'L', 'F', 'L', 'F'
};

// Base to index conversion
__device__ inline int base_to_index(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}

// Translate codon to amino acid
__device__ inline char translate_codon_device(const char* codon) {
    int idx = 0;
    for (int i = 0; i < 3; i++) {
        int base = base_to_index(codon[i]);
        if (base < 0) return 'X';
        idx = (idx << 2) | base;
    }
    return GENETIC_CODE[idx];
}

// Get BLOSUM62 score (simplified version)
__device__ inline float get_blosum_score_simple(char aa1, char aa2) {
    if (aa1 == aa2) return 5.0f;  // Match
    
    // Simplified BLOSUM-like scoring
    if ((aa1 == 'L' && aa2 == 'I') || (aa1 == 'I' && aa2 == 'L')) return 2.0f;
    if ((aa1 == 'K' && aa2 == 'R') || (aa1 == 'R' && aa2 == 'K')) return 3.0f;
    if ((aa1 == 'F' && aa2 == 'Y') || (aa1 == 'Y' && aa2 == 'F')) return 3.0f;
    if ((aa1 == 'S' && aa2 == 'T') || (aa1 == 'T' && aa2 == 'S')) return 1.0f;
    
    return -2.0f;  // Mismatch
}

// Check if amino acid is a known resistance mutation
__device__ inline bool is_resistance_mutation(uint32_t gene_id, uint16_t position, char observed_aa, float* resistance_score) {
    for (int i = 0; RESISTANCE_POSITIONS[i].position != 0; i++) {
        const ResistancePosition& rp = RESISTANCE_POSITIONS[i];
        
        if (rp.gene_id == gene_id && rp.position == position) {
            for (int j = 0; j < rp.num_resistant; j++) {
                if (rp.resistant_aas[j] == observed_aa) {
                    *resistance_score = rp.resistance_score;
                    return true;
                }
            }
            break;
        }
    }
    *resistance_score = 0.0f;
    return false;
}

// Enhanced mutation detection kernel
__global__ void enhanced_mutation_detection_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const int num_reads,
    const char* protein_sequences,     // Database protein sequences
    const uint32_t* protein_offsets,   // Offsets in protein_sequences
    const uint32_t* gene_ids,          // Gene ID for each protein
    const uint32_t* protein_starts,    // Start position of each protein in gene
    const uint32_t num_proteins,
    EnhancedProteinMatch* matches,
    uint32_t* match_counts,
    const uint32_t max_matches_per_read
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const char* read = reads + read_offsets[tid];
    const int read_len = read_lengths[tid];
    
    EnhancedProteinMatch* read_matches = &matches[tid * max_matches_per_read];
    uint32_t match_count = 0;
    
    // Try all 6 reading frames
    for (int frame = -3; frame <= 3; frame++) {
        if (frame == 0) continue;  // No frame 0
        
        // Translate read in current frame
        char translated_sequence[200];
        int translated_length = 0;
        int start_pos, step;
        
        if (frame > 0) {
            start_pos = frame - 1;
            step = 3;
        } else {
            start_pos = read_len + frame;
            step = -3;
        }
        
        // Translate sequence
        for (int pos = start_pos; 
             (frame > 0 ? pos + 2 < read_len : pos >= 2) && translated_length < 199; 
             pos += step) {
            
            char codon[4];
            if (frame > 0) {
                codon[0] = read[pos];
                codon[1] = read[pos + 1];
                codon[2] = read[pos + 2];
            } else {
                // Reverse complement
                char bases[3] = {read[pos], read[pos - 1], read[pos - 2]};
                for (int i = 0; i < 3; i++) {
                    switch(bases[i]) {
                        case 'A': case 'a': codon[i] = 'T'; break;
                        case 'T': case 't': codon[i] = 'A'; break;
                        case 'G': case 'g': codon[i] = 'C'; break;
                        case 'C': case 'c': codon[i] = 'G'; break;
                        default: codon[i] = 'N'; break;
                    }
                }
            }
            codon[3] = '\0';
            
            char aa = translate_codon_device(codon);
            if (aa == '*') {
                if (translated_length >= 20) {  // Minimum peptide length
                    break;  // End of ORF
                }
                translated_length = 0;  // Reset for next ORF
                continue;
            }
            
            translated_sequence[translated_length++] = aa;
        }
        
        if (translated_length < 20) continue;  // Too short
        translated_sequence[translated_length] = '\0';
        
        // Compare against protein database with wild-type mutation detection
        for (uint32_t protein_idx = 0; protein_idx < num_proteins && match_count < max_matches_per_read; protein_idx++) {
            const char* protein_seq = protein_sequences + protein_offsets[protein_idx];
            uint32_t gene_id = gene_ids[protein_idx];
            uint32_t protein_start = protein_starts[protein_idx];
            
            // Simple local alignment (sliding window)
            int best_score = 0;
            int best_query_pos = 0;
            int best_ref_pos = 0;
            int best_length = 0;
            
            for (int ref_pos = 0; ref_pos < 100; ref_pos++) {  // Limit search space
                for (int query_pos = 0; query_pos < translated_length - 10; query_pos++) {
                    int score = 0;
                    int length = 0;
                    
                    // Score alignment
                    while (query_pos + length < translated_length && 
                           ref_pos + length < 100 &&  // Assume max protein length 100
                           protein_seq[ref_pos + length] != '\0') {
                        
                        char query_aa = translated_sequence[query_pos + length];
                        char ref_aa = protein_seq[ref_pos + length];
                        
                        score += (int)get_blosum_score_simple(query_aa, ref_aa);
                        length++;
                        
                        if (length >= 15 && score > best_score) {
                            best_score = score;
                            best_query_pos = query_pos;
                            best_ref_pos = ref_pos;
                            best_length = length;
                        }
                    }
                }
            }
            
            // Check if alignment is good enough
            if (best_score >= 30 && best_length >= 15) {
                EnhancedProteinMatch& match = read_matches[match_count];
                
                match.read_id = tid;
                match.frame = frame;
                match.protein_id = protein_idx;
                match.gene_id = gene_id;
                match.species_id = 0;  // Simplified
                match.query_start = best_query_pos;
                match.ref_start = best_ref_pos;
                match.match_length = best_length;
                match.alignment_score = best_score;
                match.blosum_total_score = best_score;
                match.used_smith_waterman = false;
                
                // Enhanced mutation detection against wild-type
                match.num_resistance_mutations = 0;
                match.num_total_mutations = 0;
                match.covers_qrdr = false;
                
                // Copy sequences for debugging
                int copy_len = min(best_length, 99);
                for (int i = 0; i < copy_len; i++) {
                    match.query_sequence[i] = translated_sequence[best_query_pos + i];
                    match.ref_sequence[i] = protein_seq[best_ref_pos + i];
                }
                match.query_sequence[copy_len] = '\0';
                match.ref_sequence[copy_len] = '\0';
                
                // Check each position for mutations
                for (int pos = 0; pos < best_length && match.num_total_mutations < 10; pos++) {
                    char query_aa = translated_sequence[best_query_pos + pos];
                    char ref_aa = protein_seq[best_ref_pos + pos];
                    
                    int global_position = protein_start + best_ref_pos + pos + 1; // 1-based
                    
                    // Check if this position is a known resistance site
                    float resistance_score;
                    bool is_resistance_site = false;
                    char wildtype_aa = ref_aa;  // Default to reference
                    
                    // Override with known wild-type if this is a resistance position
                    for (int i = 0; RESISTANCE_POSITIONS[i].position != 0; i++) {
                        const ResistancePosition& rp = RESISTANCE_POSITIONS[i];
                        if (rp.gene_id == gene_id && rp.position == global_position) {
                            wildtype_aa = rp.wildtype_aa;
                            match.covers_qrdr = true;
                            
                            // Check if observed AA is a resistance mutation
                            if (is_resistance_mutation(gene_id, global_position, query_aa, &resistance_score)) {
                                match.resistance_positions[match.num_resistance_mutations] = global_position;
                                match.wildtype_aas[match.num_resistance_mutations] = wildtype_aa;
                                match.observed_aas[match.num_resistance_mutations] = query_aa;
                                match.resistance_scores[match.num_resistance_mutations] = resistance_score;
                                match.is_known_resistance[match.num_resistance_mutations] = true;
                                match.num_resistance_mutations++;
                                is_resistance_site = true;
                            }
                            break;
                        }
                    }
                    
                    // Record any mutation (resistance or not)
                    if (query_aa != wildtype_aa) {
                        if (match.num_total_mutations < 10) {
                            if (!is_resistance_site) {
                                match.resistance_positions[match.num_total_mutations] = global_position;
                                match.wildtype_aas[match.num_total_mutations] = wildtype_aa;
                                match.observed_aas[match.num_total_mutations] = query_aa;
                                match.resistance_scores[match.num_total_mutations] = 0.0f;
                                match.is_known_resistance[match.num_total_mutations] = false;
                            }
                            match.num_total_mutations++;
                        }
                    }
                }
                
                // Calculate identity based on wild-type comparison
                int matches = 0;
                for (int pos = 0; pos < best_length; pos++) {
                    char query_aa = translated_sequence[best_query_pos + pos];
                    int global_position = protein_start + best_ref_pos + pos + 1;
                    
                    // Get wild-type amino acid for this position
                    char wildtype_aa = protein_seq[best_ref_pos + pos];  // Default
                    for (int i = 0; RESISTANCE_POSITIONS[i].position != 0; i++) {
                        const ResistancePosition& rp = RESISTANCE_POSITIONS[i];
                        if (rp.gene_id == gene_id && rp.position == global_position) {
                            wildtype_aa = rp.wildtype_aa;
                            break;
                        }
                    }
                    
                    if (query_aa == wildtype_aa) matches++;
                }
                
                match.identity = (float)matches / best_length;
                
                match_count++;
            }
        }
    }
    
    match_counts[tid] = match_count;
}

// Host wrapper function
extern "C" {
    int enhanced_mutation_detection(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        const char* d_protein_sequences,
        const uint32_t* d_protein_offsets,
        const uint32_t* d_gene_ids,
        const uint32_t* d_protein_starts,
        uint32_t num_proteins,
        EnhancedProteinMatch* d_matches,
        uint32_t* d_match_counts,
        uint32_t max_matches_per_read
    ) {
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        enhanced_mutation_detection_kernel<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets, num_reads,
            d_protein_sequences, d_protein_offsets, d_gene_ids, d_protein_starts, num_proteins,
            d_matches, d_match_counts, max_matches_per_read
        );
        
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            return -1;
        }
        
        cudaDeviceSynchronize();
        return 0;
    }
}

#ifdef __NVCC__
#ifndef __CUDA_ARCH__
// Enhanced mutation analysis for host-side processing
class EnhancedMutationAnalyzer {
private:
    struct KnownMutation {
        uint32_t gene_id;
        int position;
        char wildtype_aa;
        std::vector<char> resistant_aas;
        float resistance_level;
        std::string description;
    };
    
    std::vector<KnownMutation> known_mutations = {
        {0, 83, 'S', {'L', 'F', 'W', 'A'}, 0.9f, "gyrA S83L/F/W/A - High-level FQ resistance"},
        {0, 87, 'D', {'N', 'G', 'Y', 'H'}, 0.7f, "gyrA D87N/G/Y/H - Moderate FQ resistance"},
        {1, 80, 'S', {'I', 'R', 'F'}, 0.6f, "parC S80I/R/F - Moderate FQ resistance"},
        {1, 84, 'E', {'V', 'K', 'G', 'A'}, 0.4f, "parC E84V/K/G/A - Low-moderate FQ resistance"},
        {2, 426, 'G', {'D', 'N'}, 0.5f, "gyrB G426D/N - FQ resistance"},
        {2, 447, 'A', {'V', 'T'}, 0.3f, "gyrB A447V/T - Low FQ resistance"},
        {3, 416, 'G', {'D', 'S'}, 0.4f, "parE G416D/S - FQ resistance"},
        {3, 420, 'A', {'V', 'T'}, 0.3f, "parE A420V/T - Low FQ resistance"}
    };
    
public:
    struct MutationReport {
        std::string gene_name;
        int position;
        char wildtype_aa;
        char observed_aa;
        float resistance_score;
        std::string description;
        bool is_qrdr;
        bool is_known_resistance;
    };
    
    std::vector<MutationReport> analyzeMutations(const EnhancedProteinMatch& match) {
        std::vector<MutationReport> reports;
        
        std::map<uint32_t, std::string> gene_names = {
            {0, "gyrA"}, {1, "parC"}, {2, "gyrB"}, {3, "parE"}
        };
        
        // Analyze resistance mutations
        for (int i = 0; i < match.num_resistance_mutations; i++) {
            MutationReport report;
            report.gene_name = gene_names.count(match.gene_id) ? 
                              gene_names[match.gene_id] : "Unknown";
            report.position = match.resistance_positions[i];
            report.wildtype_aa = match.wildtype_aas[i];
            report.observed_aa = match.observed_aas[i];
            report.resistance_score = match.resistance_scores[i];
            report.is_known_resistance = match.is_known_resistance[i];
            report.is_qrdr = true;  // All resistance mutations are in QRDR by definition
            
            // Find description
            for (const auto& km : known_mutations) {
                if (km.gene_id == match.gene_id && km.position == report.position) {
                    if (std::find(km.resistant_aas.begin(), km.resistant_aas.end(), 
                                 report.observed_aa) != km.resistant_aas.end()) {
                        report.description = km.description;
                        break;
                    }
                }
            }
            
            if (report.description.empty()) {
                report.description = "Novel mutation at resistance position";
            }
            
            reports.push_back(report);
        }
        
        return reports;
    }
    
    std::string generateMutationSummary(const EnhancedProteinMatch& match) {
        if (match.num_resistance_mutations == 0) {
            if (match.covers_qrdr) {
                return "QRDR covered - No resistance mutations detected";
            } else {
                return "No QRDR coverage in this alignment";
            }
        }
        
        std::stringstream summary;
        summary << "RESISTANCE DETECTED: " << (int)match.num_resistance_mutations << " mutations";
        
        auto reports = analyzeMutations(match);
        for (const auto& report : reports) {
            summary << "\n  " << report.gene_name << " " << report.wildtype_aa 
                   << report.position << report.observed_aa 
                   << " (score: " << std::fixed << std::setprecision(2) << report.resistance_score << ")";
        }
        
        return summary.str();
    }
    
    float calculateOverallResistanceScore(const EnhancedProteinMatch& match) {
        float max_score = 0.0f;
        for (int i = 0; i < match.num_resistance_mutations; i++) {
            if (match.resistance_scores[i] > max_score) {
                max_score = match.resistance_scores[i];
            }
        }
        return max_score;
    }
};
#endif // __CUDA_ARCH__
#endif // __NVCC__

// Debugging function to print alignment details
extern "C" {
    void print_enhanced_match_debug(const EnhancedProteinMatch* match) {
        printf("=== Enhanced Match Debug ===\n");
        printf("Read %u, Frame %d, Gene %u\n", match->read_id, match->frame, match->gene_id);
        printf("Position: %u-%u, Length: %u\n", match->ref_start, 
               match->ref_start + match->match_length - 1, match->match_length);
        printf("Score: %.1f, Identity: %.1f%%\n", match->alignment_score, match->identity * 100);
        printf("Covers QRDR: %s\n", match->covers_qrdr ? "YES" : "NO");
        printf("Resistance mutations: %u/%u total mutations\n", 
               match->num_resistance_mutations, match->num_total_mutations);
        
        if (match->num_resistance_mutations > 0) {
            printf("RESISTANCE MUTATIONS:\n");
            for (int i = 0; i < match->num_resistance_mutations; i++) {
                printf("  Position %u: %c -> %c (score: %.2f, known: %s)\n",
                       match->resistance_positions[i],
                       match->wildtype_aas[i],
                       match->observed_aas[i],
                       match->resistance_scores[i],
                       match->is_known_resistance[i] ? "YES" : "NO");
            }
        }
        
        printf("Query:  %s\n", match->query_sequence);
        printf("Ref:    %s\n", match->ref_sequence);
        printf("========================\n\n");
    }
}