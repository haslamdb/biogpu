// build_arg_database.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "arg_database.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <ncbi_arg.fasta> <output_dir> <protein_db_path>" << std::endl;
        return 1;
    }
    
    std::string fasta_path = argv[1];
    std::string output_dir = argv[2];
    std::string protein_db_path = argv[3];
    
    // Create ARG database
    ARGDatabase arg_db;
    
    // Load from NCBI ARG database
    if (!arg_db.loadFromNCBI(fasta_path, "")) {
        std::cerr << "Failed to load NCBI ARG database" << std::endl;
        return 1;
    }
    
    // Build protein sequences and k-mer index
    std::cout << "Building protein k-mer index for ARGs..." << std::endl;
    
    // This integrates with your existing protein database builder
    // The ARGs become additional proteins in your database
    // but with special metadata marking them as ARGs
    
    std::string arg_proteins_path = output_dir + "/arg_proteins.fasta";
    std::string arg_metadata_path = output_dir + "/arg_metadata.json";
    
    arg_db.exportProteinSequences(arg_proteins_path);
    arg_db.exportMetadata(arg_metadata_path);
    
    // Run your existing protein indexer
    std::string cmd = "./build_protein_index " + arg_proteins_path + " " + 
                     output_dir + " --append-to " + protein_db_path;
    system(cmd.c_str());
    
    // Create mapping file
    arg_db.createProteinToARGMapping(output_dir + "/protein_to_arg.bin");
    
    std::cout << "ARG database built successfully" << std::endl;
    return 0;
}
