// Test program to verify fixed_database_loader functionality
#include "runtime/kernels/resistance/fixed_database_loader.cpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <database_path>" << std::endl;
        std::cerr << "\nExamples:" << std::endl;
        std::cerr << "  JSON approach:      " << argv[0] << " data/integrated_fq_resistance_database_v2/protein_database.json" << std::endl;
        std::cerr << "  Directory approach: " << argv[0] << " data/integrated_fq_resistance_database_v2" << std::endl;
        return 1;
    }
    
    SpeciesGeneDatabaseLoader db_loader;
    
    std::string db_path = argv[1];
    std::cout << "\n=== Testing Fixed Database Loader ===" << std::endl;
    std::cout << "Database path: " << db_path << std::endl;
    
    bool success = false;
    
    // Try JSON first
    if (db_path.find(".json") != std::string::npos) {
        std::cout << "\nAttempting to load as JSON format..." << std::endl;
        success = db_loader.loadFromJSON(db_path);
    } else {
        std::cout << "\nAttempting to load as directory format (TSV/FASTA)..." << std::endl;
        std::string seq_file = db_path + "/wildtype_proteins.fasta";
        std::string map_file = db_path + "/species_gene_to_id.tsv";
        success = db_loader.loadFromTSV(seq_file, map_file);
    }
    
    if (!success) {
        std::cerr << "ERROR: Failed to load database!" << std::endl;
        return 1;
    }
    
    // Print summary
    db_loader.printSummary();
    
    // Test some lookups
    std::cout << "\n=== Testing Database Lookups ===" << std::endl;
    
    // Test species-gene ID lookup
    int ecoli_gyrA_id = db_loader.getSpeciesGeneId("Escherichia_coli", "gyrA");
    if (ecoli_gyrA_id >= 0) {
        std::cout << "E. coli gyrA ID: " << ecoli_gyrA_id << std::endl;
        
        // Get protein info
        const SpeciesGeneProtein* protein = db_loader.getProtein(ecoli_gyrA_id);
        if (protein) {
            std::cout << "  Sequence length: " << protein->length << " aa" << std::endl;
            std::cout << "  First 50 aa: " << protein->sequence.substr(0, 50) << "..." << std::endl;
        }
    }
    
    // Test reverse lookup
    auto species_gene = db_loader.getSpeciesGene(0);
    std::cout << "\nID 0 maps to: " << species_gene.first << " " << species_gene.second << std::endl;
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    return 0;
}