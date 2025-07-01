#include "ncbi_amr_database_loader.h"
#include <iostream>
#include <iomanip>
#include <cuda_runtime.h>

int main() {
    NCBIAMRDatabaseLoader loader;
    
    std::cout << "=== Testing NCBI AMR Database Loader ===" << std::endl;
    
    // Load the FASTA files
    bool success = loader.loadFromFastaFiles(
        "amr_nucl.fasta",
        "amr_prot.fasta"
    );
    
    if (!success) {
        std::cerr << "Failed to load AMR database!" << std::endl;
        return 1;
    }
    
    // Get the gene entries from CPU memory to examine them
    // We need to copy back from GPU to examine the data
    uint32_t num_genes = loader.getNumGenes();
    AMRGeneEntry* gpu_genes = loader.getGPUGeneEntries();
    
    // Allocate CPU memory for gene entries
    AMRGeneEntry* cpu_genes = new AMRGeneEntry[num_genes];
    
    // Copy gene entries from GPU to CPU
    cudaMemcpy(cpu_genes, gpu_genes, num_genes * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    // Get the sequence data pointers
    char* gpu_dna = loader.getGPUDNASequences();
    char* gpu_proteins = loader.getGPUProteinSequences();
    
    std::cout << "\n=== First 20 Gene Entries ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Accession" 
              << std::setw(20) << "Gene Name" 
              << std::setw(15) << "Gene Family"
              << std::setw(20) << "Resistance Class"
              << std::setw(10) << "DNA Len"
              << std::setw(10) << "Prot Len"
              << std::setw(15) << "Has Protein"
              << std::setw(10) << "Identity"
              << std::setw(10) << "Coverage"
              << std::endl;
    std::cout << std::string(150, '-') << std::endl;
    
    // Display first 20 genes
    for (int i = 0; i < std::min(20, (int)num_genes); i++) {
        AMRGeneEntry& gene = cpu_genes[i];
        
        std::cout << std::left << std::setw(20) << gene.accession
                  << std::setw(20) << gene.gene_name
                  << std::setw(15) << gene.gene_family
                  << std::setw(20) << gene.class_
                  << std::setw(10) << gene.dna_length
                  << std::setw(10) << gene.protein_length
                  << std::setw(15) << (gene.has_protein_match ? "Yes" : "No")
                  << std::setw(10) << std::fixed << std::setprecision(2) << gene.identity_threshold
                  << std::setw(10) << std::fixed << std::setprecision(2) << gene.coverage_threshold
                  << std::endl;
    }
    
    // Show some with protein matches
    std::cout << "\n=== Genes with Protein Matches (First 10) ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Accession" 
              << std::setw(20) << "Gene Name" 
              << std::setw(15) << "Gene Family"
              << std::setw(20) << "Resistance Class"
              << std::setw(10) << "DNA Len"
              << std::setw(10) << "Prot Len"
              << "Description (truncated)"
              << std::endl;
    std::cout << std::string(150, '-') << std::endl;
    
    int count = 0;
    for (int i = 0; i < num_genes && count < 10; i++) {
        AMRGeneEntry& gene = cpu_genes[i];
        if (gene.has_protein_match) {
            std::cout << std::left << std::setw(20) << gene.accession
                      << std::setw(20) << gene.gene_name
                      << std::setw(15) << gene.gene_family
                      << std::setw(20) << gene.class_
                      << std::setw(10) << gene.dna_length
                      << std::setw(10) << gene.protein_length
                      << std::string(gene.description).substr(0, 50) << "..."
                      << std::endl;
            count++;
        }
    }
    
    // Clean up
    delete[] cpu_genes;
    
    std::cout << "\nDetailed test completed successfully!" << std::endl;
    
    return 0;
}