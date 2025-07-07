// simplified_translated_search.cu
// Simplified GPU-accelerated 6-frame translation and protein search
// Focus: Basic alignment and coverage statistics without complex k-mer indexing

#ifndef SIMPLIFIED_TRANSLATED_SEARCH_CU
#define SIMPLIFIED_TRANSLATED_SEARCH_CU

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstring>

// Simple constants
#define MAX_READ_LENGTH 300
#define MAX_PROTEIN_LENGTH 500
#define MIN_ALIGNMENT_LENGTH 20
#define MIN_IDENTITY 0.80f
#define MAX_PROTEINS 50000
#define MAX_MATCHES_PER_READ 10

// Genetic code table
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

// Simple structures
struct SimpleProtein {
    char sequence[MAX_PROTEIN_LENGTH];
    char name[64];
    char drug_class[32];
    uint16_t length;
    uint32_t gene_id;
};

struct SimpleAlignment {
    uint32_t read_id;
    uint32_t protein_id;
    uint32_t gene_id;
    int8_t frame;
    uint16_t read_start;
    uint16_t read_end;
    uint16_t protein_start;
    uint16_t protein_end;
    float identity;
    bool valid;
};

struct SimpleCoverage {
    uint32_t gene_id;
    uint32_t total_reads;
    uint32_t covered_positions;
    uint16_t gene_length;
    float percent_coverage;
    float mean_depth;
    uint32_t* position_counts;  // Will be allocated separately
};

// Simple protein database
struct SimpleProteinDB {
    SimpleProtein* proteins;
    uint32_t num_proteins;
    bool loaded;
};

// Device functions
__device__ inline int base_to_index(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}

__device__ inline char translate_codon(const char* codon) {
    int idx = 0;
    for (int i = 0; i < 3; i++) {
        int base = base_to_index(codon[i]);
        if (base < 0) return 'X';
        idx = (idx << 2) | base;
    }
    return GENETIC_CODE[idx];
}

__device__ inline char reverse_complement_base(char base) {
    switch(base) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'G': case 'g': return 'C';
        case 'C': case 'c': return 'G';
        default: return 'N';
    }
}

// Simple alignment scoring
__device__ float score_alignment(const char* query, const char* ref, int length) {
    int matches = 0;
    for (int i = 0; i < length; i++) {
        if (query[i] == ref[i]) {
            matches++;
        }
    }
    return (float)matches / length;
}

// Find best local alignment using simple sliding window
__device__ void find_best_alignment(
    const char* query, int query_len,
    const char* ref, int ref_len,
    int* best_query_start, int* best_ref_start, 
    int* best_length, float* best_identity
) {
    *best_identity = 0.0f;
    *best_length = 0;
    
    // Try all possible alignments
    for (int q_start = 0; q_start < query_len; q_start++) {
        for (int r_start = 0; r_start < ref_len; r_start++) {
            int max_len = min(query_len - q_start, ref_len - r_start);
            if (max_len < MIN_ALIGNMENT_LENGTH) continue;
            
            // Try different alignment lengths
            for (int len = MIN_ALIGNMENT_LENGTH; len <= max_len; len++) {
                float identity = score_alignment(
                    &query[q_start], &ref[r_start], len
                );
                
                if (identity > *best_identity) {
                    *best_identity = identity;
                    *best_length = len;
                    *best_query_start = q_start;
                    *best_ref_start = r_start;
                }
            }
        }
    }
}

// Main kernel: translate reads and align to proteins
__global__ void simple_translate_and_align_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    int num_reads,
    const SimpleProteinDB* protein_db,
    SimpleAlignment* alignments,
    uint32_t* alignment_counts
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    const char* read = reads + read_offsets[tid];
    int read_len = read_lengths[tid];
    
    if (read_len < MIN_ALIGNMENT_LENGTH * 3) {
        alignment_counts[tid] = 0;
        return;
    }
    
    SimpleAlignment* read_alignments = &alignments[tid * MAX_MATCHES_PER_READ];
    uint32_t align_count = 0;
    
    // Translate in 6 frames
    for (int frame = -3; frame <= 3; frame++) {
        if (frame == 0) continue;
        
        char translated[MAX_READ_LENGTH];
        int trans_len = 0;
        
        if (frame > 0) {
            // Forward frames
            for (int pos = frame - 1; pos + 2 < read_len; pos += 3) {
                translated[trans_len++] = translate_codon(&read[pos]);
                if (trans_len >= MAX_READ_LENGTH - 1) break;
            }
        } else {
            // Reverse frames
            char rc_read[MAX_READ_LENGTH];
            for (int i = 0; i < read_len; i++) {
                rc_read[i] = reverse_complement_base(read[read_len - 1 - i]);
            }
            
            for (int pos = (-frame) - 1; pos + 2 < read_len; pos += 3) {
                translated[trans_len++] = translate_codon(&rc_read[pos]);
                if (trans_len >= MAX_READ_LENGTH - 1) break;
            }
        }
        
        if (trans_len < MIN_ALIGNMENT_LENGTH) continue;
        translated[trans_len] = '\0';
        
        // Align against all proteins
        for (uint32_t p = 0; p < protein_db->num_proteins && align_count < MAX_MATCHES_PER_READ; p++) {
            const SimpleProtein& protein = protein_db->proteins[p];
            
            int best_query_start, best_ref_start, best_length;
            float best_identity;
            
            find_best_alignment(
                translated, trans_len,
                protein.sequence, protein.length,
                &best_query_start, &best_ref_start,
                &best_length, &best_identity
            );
            
            if (best_identity >= MIN_IDENTITY && best_length >= MIN_ALIGNMENT_LENGTH) {
                SimpleAlignment& align = read_alignments[align_count];
                align.read_id = tid;
                align.protein_id = p;
                align.gene_id = protein.gene_id;
                align.frame = frame;
                align.read_start = best_query_start;
                align.read_end = best_query_start + best_length;
                align.protein_start = best_ref_start;
                align.protein_end = best_ref_start + best_length;
                align.identity = best_identity;
                align.valid = true;
                
                align_count++;
            }
        }
    }
    
    alignment_counts[tid] = align_count;
}

// Coverage calculation kernel
__global__ void calculate_coverage_kernel(
    const SimpleAlignment* alignments,
    const uint32_t* alignment_counts,
    int num_reads,
    SimpleCoverage* coverage_stats,
    const SimpleProteinDB* protein_db
) {
    int gene_id = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Find the maximum gene_id first
    uint32_t max_gene_id = 0;
    for (uint32_t p = 0; p < protein_db->num_proteins; p++) {
        if (protein_db->proteins[p].gene_id > max_gene_id) {
            max_gene_id = protein_db->proteins[p].gene_id;
        }
    }
    
    if (gene_id > max_gene_id) return;
    
    SimpleCoverage& stats = coverage_stats[gene_id];
    stats.gene_id = gene_id;
    stats.total_reads = 0;
    stats.covered_positions = 0;
    stats.gene_length = 0;
    
    // Find gene length from any protein with this gene_id
    for (uint32_t p = 0; p < protein_db->num_proteins; p++) {
        if (protein_db->proteins[p].gene_id == gene_id) {
            stats.gene_length = protein_db->proteins[p].length;
            break;
        }
    }
    
    if (stats.gene_length == 0) return;
    
    // Count alignments for this gene
    for (int r = 0; r < num_reads; r++) {
        uint32_t count = alignment_counts[r];
        const SimpleAlignment* read_aligns = &alignments[r * MAX_MATCHES_PER_READ];
        
        for (uint32_t a = 0; a < count; a++) {
            if (read_aligns[a].gene_id == gene_id && read_aligns[a].valid) {
                atomicAdd(&stats.total_reads, 1);
                
                // Mark covered positions
                for (uint16_t pos = read_aligns[a].protein_start; 
                     pos < read_aligns[a].protein_end && pos < stats.gene_length; 
                     pos++) {
                    if (stats.position_counts) {
                        atomicAdd(&stats.position_counts[pos], 1);
                    }
                }
            }
        }
    }
    
    // Calculate coverage percentage
    if (stats.position_counts) {
        uint32_t covered = 0;
        uint32_t total_depth = 0;
        
        for (uint16_t pos = 0; pos < stats.gene_length; pos++) {
            if (stats.position_counts[pos] > 0) {
                covered++;
                total_depth += stats.position_counts[pos];
            }
        }
        
        stats.covered_positions = covered;
        stats.percent_coverage = (float)covered / stats.gene_length * 100.0f;
        stats.mean_depth = covered > 0 ? (float)total_depth / covered : 0.0f;
    }
}

// Host class for simple translated search
class SimpleTranslatedSearch {
private:
    SimpleProteinDB* d_protein_db;
    SimpleProteinDB h_protein_db;
    
    SimpleAlignment* d_alignments;
    uint32_t* d_alignment_counts;
    SimpleCoverage* d_coverage_stats;
    
    int max_batch_size;
    int max_genes;
    bool initialized;
    
public:
    SimpleTranslatedSearch(int batch_size = 10000) 
        : max_batch_size(batch_size), max_genes(1000), initialized(false) {
        
        // Allocate GPU memory
        cudaMalloc(&d_alignments, batch_size * MAX_MATCHES_PER_READ * sizeof(SimpleAlignment));
        cudaMalloc(&d_alignment_counts, batch_size * sizeof(uint32_t));
        cudaMalloc(&d_coverage_stats, max_genes * sizeof(SimpleCoverage));
        cudaMalloc(&d_protein_db, sizeof(SimpleProteinDB));
        
        h_protein_db.proteins = nullptr;
        h_protein_db.num_proteins = 0;
        h_protein_db.loaded = false;
        
        initialized = true;
    }
    
    ~SimpleTranslatedSearch() {
        if (d_alignments) cudaFree(d_alignments);
        if (d_alignment_counts) cudaFree(d_alignment_counts);
        if (d_coverage_stats) cudaFree(d_coverage_stats);
        
        if (h_protein_db.proteins) {
            // Free proteins on GPU
            cudaFree(h_protein_db.proteins);
        }
        
        if (d_protein_db) cudaFree(d_protein_db);
    }
    
    bool loadProteins(const std::string& fasta_path) {
        printf("Loading proteins from %s\n", fasta_path.c_str());
        
        std::ifstream file(fasta_path);
        if (!file.is_open()) {
            printf("Error: Cannot open %s\n", fasta_path.c_str());
            return false;
        }
        
        std::vector<SimpleProtein> proteins;
        std::string line, header, sequence;
        uint32_t gene_counter = 0;
        
        while (std::getline(file, line)) {
            if (line[0] == '>') {
                // Process previous protein
                if (!header.empty() && !sequence.empty()) {
                    SimpleProtein protein = {};
                    
                    // Copy sequence
                    int len = std::min((int)sequence.length(), MAX_PROTEIN_LENGTH - 1);
                    strncpy(protein.sequence, sequence.c_str(), len);
                    protein.sequence[len] = '\0';
                    protein.length = len;
                    
                    // Extract name from header
                    size_t space_pos = header.find(' ');
                    std::string name = (space_pos != std::string::npos) ? 
                                     header.substr(1, space_pos - 1) : header.substr(1);
                    strncpy(protein.name, name.c_str(), 63);
                    protein.name[63] = '\0';
                    
                    // Simple drug class assignment
                    if (header.find("bla") != std::string::npos) {
                        strcpy(protein.drug_class, "BETA_LACTAM");
                    } else if (header.find("van") != std::string::npos) {
                        strcpy(protein.drug_class, "GLYCOPEPTIDE");
                    } else if (header.find("tet") != std::string::npos) {
                        strcpy(protein.drug_class, "TETRACYCLINE");
                    } else if (header.find("qnr") != std::string::npos || 
                              header.find("gyr") != std::string::npos ||
                              header.find("par") != std::string::npos) {
                        strcpy(protein.drug_class, "FLUOROQUINOLONE");
                    } else {
                        strcpy(protein.drug_class, "UNKNOWN");
                    }
                    
                    protein.gene_id = gene_counter++;
                    proteins.push_back(protein);
                }
                
                header = line;
                sequence.clear();
            } else {
                // Remove whitespace
                for (char c : line) {
                    if (!isspace(c)) {
                        sequence += toupper(c);
                    }
                }
            }
        }
        
        // Process last protein
        if (!header.empty() && !sequence.empty()) {
            SimpleProtein protein = {};
            
            int len = std::min((int)sequence.length(), MAX_PROTEIN_LENGTH - 1);
            strncpy(protein.sequence, sequence.c_str(), len);
            protein.sequence[len] = '\0';
            protein.length = len;
            
            size_t space_pos = header.find(' ');
            std::string name = (space_pos != std::string::npos) ? 
                             header.substr(1, space_pos - 1) : header.substr(1);
            strncpy(protein.name, name.c_str(), 63);
            protein.name[63] = '\0';
            
            strcpy(protein.drug_class, "UNKNOWN");
            protein.gene_id = gene_counter++;
            proteins.push_back(protein);
        }
        
        file.close();
        
        if (proteins.empty()) {
            printf("No proteins loaded\n");
            return false;
        }
        
        // Copy to GPU
        SimpleProtein* d_proteins;
        cudaMalloc(&d_proteins, proteins.size() * sizeof(SimpleProtein));
        cudaMemcpy(d_proteins, proteins.data(), 
                   proteins.size() * sizeof(SimpleProtein), 
                   cudaMemcpyHostToDevice);
        
        h_protein_db.proteins = d_proteins;
        h_protein_db.num_proteins = proteins.size();
        h_protein_db.loaded = true;
        
        cudaMemcpy(d_protein_db, &h_protein_db, sizeof(SimpleProteinDB), 
                   cudaMemcpyHostToDevice);
        
        printf("Loaded %zu proteins\n", proteins.size());
        return true;
    }
    
    bool searchReads(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        std::vector<SimpleAlignment>& alignments,
        std::vector<SimpleCoverage>& coverage
    ) {
        if (!initialized || !h_protein_db.loaded) {
            printf("Error: Not initialized or proteins not loaded\n");
            return false;
        }
        
        if (num_reads > max_batch_size) {
            printf("Error: num_reads (%d) exceeds batch size (%d)\n", 
                   num_reads, max_batch_size);
            return false;
        }
        
        // Clear previous results
        cudaMemset(d_alignment_counts, 0, num_reads * sizeof(uint32_t));
        cudaMemset(d_coverage_stats, 0, max_genes * sizeof(SimpleCoverage));
        
        // Launch alignment kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        simple_translate_and_align_kernel<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets, num_reads,
            d_protein_db, d_alignments, d_alignment_counts
        );
        
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("Alignment kernel error: %s\n", cudaGetErrorString(err));
            return false;
        }
        
        cudaDeviceSynchronize();
        
        // Launch coverage kernel
        int coverage_grid_size = (max_genes + block_size - 1) / block_size;
        
        calculate_coverage_kernel<<<coverage_grid_size, block_size>>>(
            d_alignments, d_alignment_counts, num_reads,
            d_coverage_stats, d_protein_db
        );
        
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("Coverage kernel error: %s\n", cudaGetErrorString(err));
            return false;
        }
        
        cudaDeviceSynchronize();
        
        // Copy results back
        std::vector<uint32_t> counts(num_reads);
        cudaMemcpy(counts.data(), d_alignment_counts, 
                   num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::vector<SimpleAlignment> all_alignments(num_reads * MAX_MATCHES_PER_READ);
        cudaMemcpy(all_alignments.data(), d_alignments,
                   all_alignments.size() * sizeof(SimpleAlignment), 
                   cudaMemcpyDeviceToHost);
        
        std::vector<SimpleCoverage> all_coverage(max_genes);
        cudaMemcpy(all_coverage.data(), d_coverage_stats,
                   max_genes * sizeof(SimpleCoverage), cudaMemcpyDeviceToHost);
        
        // Extract valid alignments
        alignments.clear();
        for (int i = 0; i < num_reads; i++) {
            for (uint32_t j = 0; j < counts[i]; j++) {
                const SimpleAlignment& align = all_alignments[i * MAX_MATCHES_PER_READ + j];
                if (align.valid) {
                    alignments.push_back(align);
                }
            }
        }
        
        // Extract coverage stats with data
        coverage.clear();
        for (int i = 0; i < max_genes; i++) {
            if (all_coverage[i].total_reads > 0) {
                coverage.push_back(all_coverage[i]);
            }
        }
        
        printf("Found %zu alignments, %zu genes with coverage\n", 
               alignments.size(), coverage.size());
        
        return true;
    }
    
    void writeResults(const std::string& output_prefix,
                     const std::vector<SimpleAlignment>& alignments,
                     const std::vector<SimpleCoverage>& coverage) {
        
        // Write alignments
        std::string align_file = output_prefix + "_alignments.tsv";
        std::ofstream align_out(align_file);
        align_out << "read_id\tprotein_id\tgene_id\tframe\tidentity\tread_start\tread_end\tprotein_start\tprotein_end\n";
        
        for (const auto& align : alignments) {
            align_out << align.read_id << "\t"
                     << align.protein_id << "\t"
                     << align.gene_id << "\t"
                     << (int)align.frame << "\t"
                     << align.identity << "\t"
                     << align.read_start << "\t"
                     << align.read_end << "\t"
                     << align.protein_start << "\t"
                     << align.protein_end << "\n";
        }
        align_out.close();
        
        // Write coverage
        std::string cov_file = output_prefix + "_coverage.tsv";
        std::ofstream cov_out(cov_file);
        cov_out << "gene_id\ttotal_reads\tcovered_positions\tgene_length\tpercent_coverage\tmean_depth\n";
        
        for (const auto& cov : coverage) {
            cov_out << cov.gene_id << "\t"
                   << cov.total_reads << "\t"
                   << cov.covered_positions << "\t"
                   << cov.gene_length << "\t"
                   << cov.percent_coverage << "\t"
                   << cov.mean_depth << "\n";
        }
        cov_out.close();
        
        printf("Results written to %s and %s\n", align_file.c_str(), cov_file.c_str());
    }
};

// C interface
extern "C" {
    void* create_simple_translated_search(int batch_size) {
        return new SimpleTranslatedSearch(batch_size);
    }
    
    void destroy_simple_translated_search(void* engine) {
        if (engine) {
            delete static_cast<SimpleTranslatedSearch*>(engine);
        }
    }
    
    int load_simple_proteins(void* engine, const char* fasta_path) {
        if (!engine) return -1;
        SimpleTranslatedSearch* sts = static_cast<SimpleTranslatedSearch*>(engine);
        return sts->loadProteins(fasta_path) ? 0 : -1;
    }
    
    int search_simple_reads(
        void* engine,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        void* alignments_out,
        void* coverage_out
    ) {
        if (!engine) return -1;
        
        SimpleTranslatedSearch* sts = static_cast<SimpleTranslatedSearch*>(engine);
        std::vector<SimpleAlignment> alignments;
        std::vector<SimpleCoverage> coverage;
        
        bool success = sts->searchReads(d_reads, d_read_lengths, d_read_offsets, 
                                       num_reads, alignments, coverage);
        
        // Note: You'll need to handle copying results back to the caller
        // This is a simplified interface
        
        return success ? 0 : -1;
    }
}

#endif // SIMPLIFIED_TRANSLATED_SEARCH_CU
