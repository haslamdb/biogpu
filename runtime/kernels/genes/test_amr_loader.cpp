#include "ncbi_amr_database_loader.h"
#include <iostream>
#include <iomanip>

int main() {
    NCBIAMRDatabaseLoader loader;
    
    std::cout << "=== Testing NCBI AMR Database Loader ===" << std::endl;
    std::cout << "Working directory: " << __FILE__ << std::endl;
    
    // Load the FASTA files
    bool success = loader.loadFromFastaFiles(
        "AMR_CDS.fa",
        "AMRProt.fa"
    );
    
    if (!success) {
        std::cerr << "Failed to load AMR database!" << std::endl;
        return 1;
    }
    
    std::cout << "\n=== Loading Summary ===" << std::endl;
    std::cout << "Successfully loaded " << loader.getNumGenes() << " AMR genes" << std::endl;
    
    // Print detailed database statistics
    loader.printDatabaseStats();
    
    // Get GPU pointers for verification
    AMRGeneEntry* gpu_genes = loader.getGPUGeneEntries();
    char* gpu_dna = loader.getGPUDNASequences();
    char* gpu_proteins = loader.getGPUProteinSequences();
    
    std::cout << "\n=== GPU Memory Allocation ===" << std::endl;
    std::cout << "GPU gene entries: " << (gpu_genes ? "Allocated" : "NULL") << std::endl;
    std::cout << "GPU DNA sequences: " << (gpu_dna ? "Allocated" : "NULL") << std::endl;
    std::cout << "GPU protein sequences: " << (gpu_proteins ? "Allocated" : "NULL") << std::endl;
    
    std::cout << "\nDatabase loader test completed successfully!" << std::endl;
    
    return 0;
}