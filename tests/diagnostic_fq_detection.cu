// diagnostic_fq_detection.cu
// Diagnostic version to troubleshoot why FQ resistance mutations aren't being detected

#include <cuda_runtime.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>

// Genetic code for translation
__device__ const char CODON_TABLE[64] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
};

// Base to index mapping
__device__ int base_to_index(char base) {
    switch(base) {
        case 'T': case 't': return 0;
        case 'C': case 'c': return 1;
        case 'A': case 'a': return 2;
        case 'G': case 'g': return 3;
        default: return -1;
    }
}

// Translate codon to amino acid
__device__ char translate_codon(const char* codon) {
    int idx1 = base_to_index(codon[0]);
    int idx2 = base_to_index(codon[1]);
    int idx3 = base_to_index(codon[2]);
    
    if (idx1 < 0 || idx2 < 0 || idx3 < 0) return 'X'; // Unknown
    
    int codon_index = idx1 * 16 + idx2 * 4 + idx3;
    return CODON_TABLE[codon_index];
}

// Structure for QRDR positions we care about
struct QRDRPosition {
    const char* gene;
    int codon_position;      // 1-based codon position
    int gene_start;          // 0-based start position in reference
    char wild_type_aa;
    const char* resistant_aas;
};

// Key positions for E. coli fluoroquinolone resistance
__constant__ QRDRPosition QRDR_POSITIONS[] = {
    {"gyrA", 83, 247, 'S', "LFW"},   // Ser83 -> Leu/Phe/Trp
    {"gyrA", 87, 259, 'D', "NGY"},   // Asp87 -> Asn/Gly/Tyr
    {"parC", 80, 238, 'S', "IR"},    // Ser80 -> Ile/Arg
    {"parC", 84, 250, 'E', "VKG"},   // Glu84 -> Val/Lys/Gly
};

// Diagnostic information structure
struct DiagnosticInfo {
    int read_id;
    int alignment_pos;
    int gene_match;
    char sequence[300];
    char translation[100];
    int frame_used;
    bool spans_qrdr;
    int mutations_checked;
    int mutations_found;
    char debug_msg[256];
};

// Check if read spans a QRDR position
__device__ bool check_spans_position(const char* read_seq, int read_len, 
                                   int align_start, int qrdr_pos) {
    int read_end = align_start + read_len;
    // Need at least 3 bases to cover a codon
    return (align_start <= qrdr_pos - 2) && (read_end >= qrdr_pos + 3);
}

// Extract and translate the region around a QRDR position
__device__ void extract_and_translate(const char* read_seq, int read_len,
                                    int align_start, int qrdr_pos,
                                    char* codon_out, char* aa_out,
                                    DiagnosticInfo* diag) {
    // Calculate position in read
    int codon_start_in_ref = qrdr_pos;
    int offset_in_read = codon_start_in_ref - align_start;
    
    // Debug info
    sprintf(diag->debug_msg + strlen(diag->debug_msg), 
            "Align:%d QRDRPos:%d Offset:%d ", 
            align_start, qrdr_pos, offset_in_read);
    
    // Check bounds
    if (offset_in_read < 0 || offset_in_read + 3 > read_len) {
        *aa_out = '?';
        return;
    }
    
    // Extract codon
    codon_out[0] = read_seq[offset_in_read];
    codon_out[1] = read_seq[offset_in_read + 1];
    codon_out[2] = read_seq[offset_in_read + 2];
    codon_out[3] = '\0';
    
    // Translate
    *aa_out = translate_codon(codon_out);
    
    // Add to debug
    sprintf(diag->debug_msg + strlen(diag->debug_msg), 
            "Codon:%s AA:%c ", codon_out, *aa_out);
}

// Main diagnostic kernel
__global__ void diagnose_fq_resistance(
    const char* reads,          // Concatenated reads
    const int* read_lengths,    // Length of each read
    const int* read_offsets,    // Start position in concatenated array
    const int* align_positions, // Where each read aligns to reference
    const int* gene_ids,        // Which gene each read maps to (0=gyrA, 1=parC, etc)
    DiagnosticInfo* diagnostics,
    int num_reads,
    int num_qrdr_positions
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_reads) return;
    
    // Get read info
    int read_len = read_lengths[tid];
    int read_offset = read_offsets[tid];
    const char* read_seq = reads + read_offset;
    int align_pos = align_positions[tid];
    int gene_id = gene_ids[tid];
    
    // Initialize diagnostic info
    DiagnosticInfo* diag = &diagnostics[tid];
    diag->read_id = tid;
    diag->alignment_pos = align_pos;
    diag->gene_match = gene_id;
    diag->mutations_checked = 0;
    diag->mutations_found = 0;
    diag->spans_qrdr = false;
    
    // Copy sequence for debugging (truncate if needed)
    int copy_len = min(read_len, 299);
    memcpy(diag->sequence, read_seq, copy_len);
    diag->sequence[copy_len] = '\0';
    
    // Clear debug message
    diag->debug_msg[0] = '\0';
    
    // Try all three reading frames
    sprintf(diag->debug_msg, "Gene%d ", gene_id);
    
    // Check each QRDR position
    for (int i = 0; i < num_qrdr_positions; i++) {
        QRDRPosition qrdr = QRDR_POSITIONS[i];
        
        // Skip if wrong gene
        if ((gene_id == 0 && strcmp(qrdr.gene, "gyrA") != 0) ||
            (gene_id == 1 && strcmp(qrdr.gene, "parC") != 0)) {
            continue;
        }
        
        // Calculate actual position in reference
        int qrdr_dna_pos = qrdr.gene_start + (qrdr.codon_position - 1) * 3;
        
        // Check if read spans this position
        if (check_spans_position(read_seq, read_len, align_pos, qrdr_dna_pos)) {
            diag->spans_qrdr = true;
            diag->mutations_checked++;
            
            char codon[4];
            char aa;
            extract_and_translate(read_seq, read_len, align_pos, 
                                qrdr_dna_pos, codon, &aa, diag);
            
            // Check if it's wild type
            if (aa == qrdr.wild_type_aa) {
                sprintf(diag->debug_msg + strlen(diag->debug_msg), 
                        "WT@%d ", qrdr.codon_position);
            } 
            // Check if it's a known resistance mutation
            else if (strchr(qrdr.resistant_aas, aa) != NULL) {
                diag->mutations_found++;
                sprintf(diag->debug_msg + strlen(diag->debug_msg), 
                        "MUT@%d:%c>%c! ", qrdr.codon_position, 
                        qrdr.wild_type_aa, aa);
            }
            // Unknown mutation
            else if (aa != '?' && aa != 'X') {
                sprintf(diag->debug_msg + strlen(diag->debug_msg), 
                        "Novel@%d:%c>%c ", qrdr.codon_position, 
                        qrdr.wild_type_aa, aa);
            }
        }
    }
    
    // Translate first 30 codons for debugging
    for (int frame = 0; frame < 3; frame++) {
        if (frame > 0) continue; // Just do frame 0 for now
        
        diag->frame_used = frame;
        int aa_count = 0;
        for (int i = frame; i + 2 < read_len && aa_count < 30; i += 3) {
            char codon[3] = {read_seq[i], read_seq[i+1], read_seq[i+2]};
            diag->translation[aa_count++] = translate_codon(codon);
        }
        diag->translation[aa_count] = '\0';
    }
}

// Host-side diagnostic function
void run_diagnostic_test() {
    // Create synthetic test data
    std::vector<std::string> test_reads = {
        // Wild type gyrA around codon 83 (S)
        "ATGGCTGAATTCCGCCATGCGCGTACTTTCGCCGAGCTGCGCGATCAGGTCTTCCGCTCATGGCGGATCCAGCAGGTTTTCCAGCTTCGATGCGA",
        
        // Mutant gyrA with S83L (TCT->TTG)
        "ATGGCTGAATTCCGCCATGCGCGTACTTTCGCCGAGCTGCGCGATCAGGTCTTGCGCTCATGGCGGATCCAGCAGGTTTTCCAGCTTCGATGCGA",
        
        // Wild type parC around codon 80
        "ATGAGCGAAATCGCGCCGTGCGTGACCGCCTGGAAGCGGTAGATCAGGAACGCTCGATGGCGGTACAGCAGGTTATCCACCTTCGACGCGGCGT",
        
        // Mutant parC with S80I (AGC->ATC)
        "ATGAGCGAAATCGCGCCGTGCGTGACCGCCTGGAAGCGGTAGATCAGGAACGCATCATGGCGGTACAGCAGGTTATCCACCTTCGACGCGGCGT"
    };
    
    // Alignment positions (where reads start in reference)
    std::vector<int> align_positions = {200, 200, 190, 190};
    std::vector<int> gene_ids = {0, 0, 1, 1}; // 0=gyrA, 1=parC
    
    // Prepare data for GPU
    std::string concatenated;
    std::vector<int> lengths, offsets;
    
    for (const auto& read : test_reads) {
        offsets.push_back(concatenated.length());
        lengths.push_back(read.length());
        concatenated += read;
    }
    
    // Allocate GPU memory
    char* d_reads;
    int *d_lengths, *d_offsets, *d_aligns, *d_genes;
    DiagnosticInfo* d_diagnostics;
    
    cudaMalloc(&d_reads, concatenated.size());
    cudaMalloc(&d_lengths, lengths.size() * sizeof(int));
    cudaMalloc(&d_offsets, offsets.size() * sizeof(int));
    cudaMalloc(&d_aligns, align_positions.size() * sizeof(int));
    cudaMalloc(&d_genes, gene_ids.size() * sizeof(int));
    cudaMalloc(&d_diagnostics, test_reads.size() * sizeof(DiagnosticInfo));
    
    // Copy to GPU
    cudaMemcpy(d_reads, concatenated.c_str(), concatenated.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lengths, lengths.data(), lengths.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, offsets.data(), offsets.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_aligns, align_positions.data(), align_positions.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_genes, gene_ids.data(), gene_ids.size() * sizeof(int), cudaMemcpyHostToDevice);
    
    // Run diagnostic kernel
    int blockSize = 256;
    int gridSize = (test_reads.size() + blockSize - 1) / blockSize;
    diagnose_fq_resistance<<<gridSize, blockSize>>>(
        d_reads, d_lengths, d_offsets, d_aligns, d_genes,
        d_diagnostics, test_reads.size(), 4
    );
    
    // Get results
    std::vector<DiagnosticInfo> diagnostics(test_reads.size());
    cudaMemcpy(diagnostics.data(), d_diagnostics, 
               test_reads.size() * sizeof(DiagnosticInfo), cudaMemcpyDeviceToHost);
    
    // Print diagnostic info
    printf("\n=== FLUOROQUINOLONE RESISTANCE DIAGNOSTIC REPORT ===\n\n");
    
    for (int i = 0; i < test_reads.size(); i++) {
        const auto& diag = diagnostics[i];
        printf("Read %d:\n", i);
        printf("  Gene: %s (id=%d)\n", diag.gene_match == 0 ? "gyrA" : "parC", diag.gene_match);
        printf("  Alignment position: %d\n", diag.alignment_pos);
        printf("  Spans QRDR: %s\n", diag.spans_qrdr ? "YES" : "NO");
        printf("  Mutations checked: %d\n", diag.mutations_checked);
        printf("  Mutations found: %d\n", diag.mutations_found);
        printf("  Debug: %s\n", diag.debug_msg);
        printf("  First 30 AA: %.30s\n", diag.translation);
        printf("  Sequence: %.60s...\n\n", diag.sequence);
    }
    
    // Cleanup
    cudaFree(d_reads);
    cudaFree(d_lengths);
    cudaFree(d_offsets);
    cudaFree(d_aligns);
    cudaFree(d_genes);
    cudaFree(d_diagnostics);
}

// Function to create synthetic resistant reads
void generate_test_fastq(const char* filename) {
    FILE* fp = fopen(filename, "w");
    
    // E. coli gyrA reference around codon 83 (includes upstream for context)
    // This is the QRDR region - codons 67-106
    const char* gyrA_ref = 
        "ATGAGCGAAGAACAAAACAGATCCGGGAAGCGCGCCGCAAAGCCCTGTACGACCTGGTC"
        "GAGCAACTGCCGGAAACCGAATTCGGCAACGGTAAGCGCATCGAAGTGATGACTGGACA"
        "AGAAATCTCGCACGGTGATCCTGCGCGACGGCAGCTATGAGCGTATCAACGATATGCGT"
        "ATGCGCTCAGCGTATGCGATTTCCGTTGGGCAATGAACAGTGGAACAACCTGTTCCGCT"
        "ATCTGCCGGTGTACGACCCGCAGGATCTGGCGATCCGTCACCGTCGCCACGGCGATATT"; // S83 is at TCG
    
    // Write wild type
    fprintf(fp, "@read1_gyrA_wildtype\n%s\n+\n", gyrA_ref);
    for (int i = 0; i < strlen(gyrA_ref); i++) fprintf(fp, "I");
    fprintf(fp, "\n");
    
    // Create S83L mutant (TCG -> TTG)
    std::string mutant = gyrA_ref;
    size_t s83_pos = mutant.find("ATCTCGCAC"); // Find context before S83
    if (s83_pos != std::string::npos) {
        s83_pos += 9; // Move to codon 83
        mutant[s83_pos] = 'T';   // C->T
        mutant[s83_pos+1] = 'T'; // C->T  
        mutant[s83_pos+2] = 'G'; // G stays G
    }
    
    fprintf(fp, "@read2_gyrA_S83L_mutant\n%s\n+\n", mutant.c_str());
    for (int i = 0; i < mutant.length(); i++) fprintf(fp, "I");
    fprintf(fp, "\n");
    
    fclose(fp);
    printf("Created test FASTQ: %s\n", filename);
}

int main() {
    printf("Running Fluoroquinolone Resistance Diagnostic Test\n");
    printf("==================================================\n\n");
    
    // Generate test data
    generate_test_fastq("test_fq_resistance.fastq");
    
    // Run diagnostic
    run_diagnostic_test();
    
    return 0;
}