// Example usage in main pipeline
int main(int argc, char** argv) {
    // Your existing arguments plus ARG database
    std::string nucleotide_index = argv[1];
    std::string protein_db = argv[2];
    std::string arg_db_path = argv[3];  // New: path to ARG database
    std::string r1_path = argv[4];
    std::string r2_path = argv[5];
    std::string output_prefix = argv[6];
    
    // Create enhanced pipeline with ARG detection
    ARGResistancePipeline pipeline(true, true, 5, 0);
    
    // Load databases
    pipeline.loadDatabases(nucleotide_index, protein_db, "");
    pipeline.loadARGDatabase(arg_db_path + "/sequences.fasta", 
                           arg_db_path + "/metadata.json");
    
    // Process reads
    pipeline.processReads(r1_path, r2_path, output_prefix);
    
    // Generate reports
    pipeline.generateCoverageReport(output_prefix + "_arg_coverage.txt");
    
    return 0;
}
