// integrated_amr_pipeline.cpp
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <map>
#include "direct_translated_amr_search.cu"
#include "ncbi_amr_database_loader.h"

class IntegratedAMRPipeline {
private:
    // Existing components
    NCBIAMRDatabaseLoader amr_db;
    
    // GPU memory
    char* d_reads;
    int* d_read_offsets;
    int* d_read_lengths;
    char* d_translated_reads;
    int* d_frame_offsets;
    int* d_frame_lengths;
    AMRMatch* d_matches;
    int* d_match_counts;
    
    // Configuration
    int max_reads_per_batch = 100000;
    int max_read_length = 300;  // For paired-end reads
    int threads_per_block = 256;
    
public:
    // Main entry point combining both approaches
    void processMetagenomicReads(
        const std::vector<std::string>& reads,
        bool use_nucleotide_search = true,
        bool use_protein_search = true
    ) {
        std::cout << "Processing " << reads.size() << " reads for AMR detection" << std::endl;
        
        if (use_protein_search) {
            // Direct translated search for all AMR proteins
            performTranslatedSearch(reads);
        }
        
        if (use_nucleotide_search) {
            // Your existing mutation detection for quinolone resistance
            performTargetedMutationDetection(reads);
        }
        
        // Combine and report results
        generateCombinedReport();
    }
    
    void performTranslatedSearch(const std::vector<std::string>& reads) {
        // Prepare read data for GPU
        prepareReadsForGPU(reads);
        
        // Allocate GPU memory for translations
        size_t translation_size = reads.size() * max_read_length * 6; // 6 frames
        cudaMalloc(&d_translated_reads, translation_size);
        cudaMalloc(&d_frame_offsets, reads.size() * 6 * sizeof(int));
        cudaMalloc(&d_frame_lengths, reads.size() * 6 * sizeof(int));
        
        // Launch translation kernel
        int blocks = (reads.size() + threads_per_block - 1) / threads_per_block;
        translateReadsKernel<<<blocks, threads_per_block>>>(
            d_reads, d_read_offsets, d_read_lengths,
            d_translated_reads, d_frame_offsets, d_frame_lengths,
            reads.size(), max_read_length
        );
        
        // Allocate memory for matches
        int max_matches_per_read = 10;
        cudaMalloc(&d_matches, reads.size() * max_matches_per_read * sizeof(AMRMatch));
        cudaMalloc(&d_match_counts, reads.size() * sizeof(int));
        cudaMemset(d_match_counts, 0, reads.size() * sizeof(int));
        
        // Launch search kernel
        int search_threads = reads.size() * 6; // One thread per frame
        int search_blocks = (search_threads + threads_per_block - 1) / threads_per_block;
        
        searchAMRProteinsKernel<<<search_blocks, threads_per_block, 4096>>>(
            d_translated_reads, d_frame_offsets, d_frame_lengths,
            amr_db.getGPUProteinSequences(),
            amr_db.getProteinOffsets(),
            amr_db.getProteinLengths(),
            amr_db.getIdentityThresholds(),
            d_matches, d_match_counts,
            reads.size(), amr_db.getNumGenes()
        );
        
        // Copy results back
        std::vector<AMRMatch> matches(reads.size() * max_matches_per_read);
        std::vector<int> match_counts(reads.size());
        
        cudaMemcpy(matches.data(), d_matches, 
                   matches.size() * sizeof(AMRMatch), cudaMemcpyDeviceToHost);
        cudaMemcpy(match_counts.data(), d_match_counts,
                   match_counts.size() * sizeof(int), cudaMemcpyDeviceToHost);
        
        // Process results
        processTranslatedSearchResults(matches, match_counts);
    }
    
    void performTargetedMutationDetection(const std::vector<std::string>& reads) {
        // Your existing code for gyrA, gyrB, parC, parE mutation detection
        // This can run in parallel with translated search
        
        // Target specific positions known for fluoroquinolone resistance
        std::vector<std::string> target_genes = {"gyrA", "gyrB", "parC", "parE"};
        std::map<std::string, std::vector<int>> target_positions = {
            {"gyrA", {83, 87}},
            {"gyrB", {426, 447}},
            {"parC", {80, 84}},
            {"parE", {445, 458}}
        };
        
        // Use your existing CUDA kernel for targeted mutation detection
        // This is more efficient for known resistance positions
    }
    
    void processTranslatedSearchResults(
        const std::vector<AMRMatch>& matches,
        const std::vector<int>& match_counts
    ) {
        // Group matches by AMR gene family
        std::map<std::string, int> gene_family_counts;
        std::map<std::string, float> gene_family_coverage;
        
        for (size_t i = 0; i < match_counts.size(); i++) {
            for (int j = 0; j < match_counts[i]; j++) {
                const AMRMatch& match = matches[i * 10 + j];
                
                // Get gene information from database
                std::string gene_name = amr_db.getGeneName(match.protein_id);
                std::string drug_class = amr_db.getDrugClass(match.protein_id);
                
                // Track high-confidence matches
                if (match.identity >= 0.95 && match.is_complete_protein) {
                    gene_family_counts[gene_name]++;
                    
                    // Calculate coverage depth
                    if (gene_family_coverage.find(gene_name) == gene_family_coverage.end()) {
                        gene_family_coverage[gene_name] = 0;
                    }
                    gene_family_coverage[gene_name] += 1.0f / match_counts.size();
                }
            }
        }
        
        // Report findings
        std::cout << "\n=== AMR Genes Detected (Translated Search) ===" << std::endl;
        for (const auto& [gene, count] : gene_family_counts) {
            std::cout << gene << ": " << count << " reads (coverage: " 
                      << gene_family_coverage[gene] * 100 << "%)" << std::endl;
        }
    }
    
    void generateCombinedReport() {
        // Combine results from both search methods
        std::cout << "\n=== Combined AMR Detection Report ===" << std::endl;
        
        // Priority findings for fluoroquinolones
        std::cout << "\nFluoroquinolone Resistance:" << std::endl;
        // Report mutations in gyrA, gyrB, parC, parE
        
        // Other resistance mechanisms
        std::cout << "\nOther AMR Genes Detected:" << std::endl;
        // Report beta-lactamases, aminoglycoside resistance, etc.
        
        // Clinical recommendations
        std::cout << "\nClinical Interpretation:" << std::endl;
        // Suggest alternative antibiotics based on resistance profile
    }
    
private:
    void prepareReadsForGPU(const std::vector<std::string>& reads) {
        // Convert reads to GPU-friendly format
        std::vector<char> concatenated_reads;
        std::vector<int> offsets;
        std::vector<int> lengths;
        
        for (const auto& read : reads) {
            offsets.push_back(concatenated_reads.size());
            lengths.push_back(read.length());
            concatenated_reads.insert(concatenated_reads.end(), 
                                    read.begin(), read.end());
        }
        
        // Allocate and copy to GPU
        cudaMalloc(&d_reads, concatenated_reads.size());
        cudaMalloc(&d_read_offsets, offsets.size() * sizeof(int));
        cudaMalloc(&d_read_lengths, lengths.size() * sizeof(int));
        
        cudaMemcpy(d_reads, concatenated_reads.data(), 
                   concatenated_reads.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, offsets.data(), 
                   offsets.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, lengths.data(), 
                   lengths.size() * sizeof(int), cudaMemcpyHostToDevice);
    }
};

// Example usage
int main() {
    IntegratedAMRPipeline pipeline;
    
    // Load your AMR database (proteins only)
    pipeline.loadAMRDatabase("amr_proteins.fasta");
    
    // Process reads
    std::vector<std::string> reads = loadReadsFromFastq("patient_sample.fastq");
    
    // Run both nucleotide and protein-based detection
    pipeline.processMetagenomicReads(reads, true, true);
    
    return 0;
}
