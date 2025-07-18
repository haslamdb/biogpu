// build_amr_protein_database.cpp
// Builds protein database for clinical diagnostic pipeline
// Compatible with TranslatedSearchEngine from translated_search_amr.cu

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cstdint>

struct ProteinEntry {
    std::string accession;
    std::string gene_name;
    std::string gene_family;
    std::string drug_class;
    std::string sequence;
    uint32_t gene_id;
    uint16_t length;
};

// Parse AMRProt.fa header format
// Format: >0|WP_033990638.1|1|1|blaPDC-11|blaPDC|hydrolase|2|CEPHALOSPORIN|BETA-LACTAM|class_C_beta-lactamase_PDC-11
void parseAMRProtHeader(const std::string& header, ProteinEntry& entry) {
    std::vector<std::string> parts;
    std::stringstream ss(header);
    std::string part;
    
    while (std::getline(ss, part, '|')) {
        parts.push_back(part);
    }
    
    // Based on actual format:
    // >0|WP_033990638.1|1|1|blaPDC-11|blaPDC|hydrolase|2|CEPHALOSPORIN|BETA-LACTAM|class_C_beta-lactamase_PDC-11
    if (parts.size() >= 10) {
        entry.accession = parts[1];        // WP_033990638.1
        entry.gene_id = std::stoi(parts[0]); // 0
        entry.gene_name = parts[4];        // blaPDC-11
        entry.gene_family = parts[5];      // blaPDC
        entry.drug_class = parts[9];       // BETA-LACTAM
        
        // Alternative drug class in position 8
        if (entry.drug_class.empty() && parts.size() > 8) {
            entry.drug_class = parts[8];   // CEPHALOSPORIN
        }
    } else if (parts.size() >= 5) {
        entry.accession = parts[1];
        // Use the index (first field) as gene_id for now
        entry.gene_id = std::stoi(parts[0]);
        entry.gene_name = parts[4]; // e.g., "blaPDC-11"
        
        // Extract gene family from parts[5] if available
        if (parts.size() > 5) {
            entry.gene_family = parts[5]; // e.g., "blaPDC"
        } else {
            // Fallback to extracting from gene name
            size_t dash_pos = entry.gene_name.find('-');
            if (dash_pos != std::string::npos) {
                entry.gene_family = entry.gene_name.substr(0, dash_pos);
            } else {
                entry.gene_family = entry.gene_name;
            }
        }
        
        // Extract drug class from later fields
        if (parts.size() > 9) {
            entry.drug_class = parts[9]; // e.g., "BETA-LACTAM"
        } else if (parts.size() > 8) {
            entry.drug_class = parts[8]; // e.g., "CEPHALOSPORIN"
        } else {
            // Fallback to inferring from gene name
            std::string gene_lower = entry.gene_name;
            std::transform(gene_lower.begin(), gene_lower.end(), gene_lower.begin(), ::tolower);
            
            if (gene_lower.find("bla") == 0) {
                entry.drug_class = "BETA-LACTAM";
            } else if (gene_lower.find("tet") == 0) {
                entry.drug_class = "TETRACYCLINE";
            } else {
                entry.drug_class = "UNKNOWN";
            }
        }
    }
}

#ifdef TEST_GENE_FAMILY_EXTRACTION
void testGeneFamilyExtraction() {
    std::cout << "\n=== Testing Gene Family Extraction ===" << std::endl;
    
    std::vector<std::pair<std::string, std::string>> test_cases = {
        {"blaKPC-2", "blaKPC"},
        {"qnrA1", "qnrA1"},
        {"vanA", "vanA"},
        {"aac(6')-Ib", "aac(6')"},
        {"tet(M)", "tet(M)"},
        {"mecA", "mecA"},
        {"blaOXA-48", "blaOXA"},
        {"blaKPC-3", "blaKPC"},
        {"blaNDM-1", "blaNDM"},
        {"armA", "armA"}
    };
    
    std::cout << "Test cases:" << std::endl;
    for (const auto& [gene_name, expected_family] : test_cases) {
        // Extract gene family using same logic as parseAMRProtHeader
        std::string name = gene_name;
        size_t dash_pos = name.find('-');
        std::string actual_family;
        if (dash_pos != std::string::npos) {
            actual_family = name.substr(0, dash_pos);
        } else {
            actual_family = name;
        }
        
        bool passed = (actual_family == expected_family);
        std::cout << "  " << gene_name << " -> " << actual_family 
                  << " (expected: " << expected_family << ") " 
                  << (passed ? "PASS" : "FAIL") << std::endl;
    }
    
    std::cout << "\nTest complete!" << std::endl;
}
#endif

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <AMRProt.fa> <output_dir>" << std::endl;
        std::cerr << "Builds protein reference database for clinical AMR gene detection" << std::endl;
        #ifdef TEST_GENE_FAMILY_EXTRACTION
        std::cerr << "       " << argv[0] << " --test  (run gene family extraction tests)" << std::endl;
        #endif
        return 1;
    }
    
    #ifdef TEST_GENE_FAMILY_EXTRACTION
    if (std::string(argv[1]) == "--test") {
        testGeneFamilyExtraction();
        return 0;
    }
    #endif
    
    std::string input_file = argv[1];
    std::string output_dir = argv[2];
    
    // Create output directory if needed
    std::string mkdir_cmd = "mkdir -p " + output_dir;
    system(mkdir_cmd.c_str());
    
    std::cout << "Building AMR protein database for clinical diagnostics..." << std::endl;
    std::cout << "Input: " << input_file << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    
    // Parse FASTA file
    std::vector<ProteinEntry> proteins;
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << input_file << std::endl;
        return 1;
    }
    
    std::string line;
    ProteinEntry current_entry;
    bool in_sequence = false;
    
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (in_sequence && !current_entry.sequence.empty()) {
                current_entry.length = current_entry.sequence.length();
                proteins.push_back(current_entry);
            }
            
            current_entry = ProteinEntry();
            parseAMRProtHeader(line.substr(1), current_entry);
            in_sequence = true;
        } else if (in_sequence) {
            // Remove whitespace and convert to uppercase
            for (char c : line) {
                if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
                    current_entry.sequence += toupper(c);
                }
            }
        }
    }
    
    if (in_sequence && !current_entry.sequence.empty()) {
        current_entry.length = current_entry.sequence.length();
        proteins.push_back(current_entry);
    }
    file.close();
    
    std::cout << "Loaded " << proteins.size() << " AMR protein sequences" << std::endl;
    
    
    // Statistics
    std::map<std::string, int> drug_class_counts;
    for (const auto& p : proteins) {
        drug_class_counts[p.drug_class]++;
    }
    
    std::cout << "\nDrug class distribution:" << std::endl;
    for (const auto& pair : drug_class_counts) {
        std::cout << "  " << pair.first << ": " << pair.second << " genes" << std::endl;
    }
    
    // Build k-mer index
    const int KMER_SIZE = 8;
    std::map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> kmer_to_positions;
    
    for (uint32_t i = 0; i < proteins.size(); i++) {
        const std::string& seq = proteins[i].sequence;
        
        for (size_t pos = 0; pos + KMER_SIZE <= seq.length(); pos++) {
            std::string kmer = seq.substr(pos, KMER_SIZE);
            kmer_to_positions[kmer].push_back({i, pos});
        }
    }
    
    std::cout << "\nBuilt k-mer index with " << kmer_to_positions.size() << " unique " 
              << KMER_SIZE << "-mers" << std::endl;
    
    // Write protein sequences (all concatenated)
    std::string protein_file = output_dir + "/proteins.bin";
    std::ofstream prot_out(protein_file, std::ios::binary);
    
    uint32_t num_proteins = proteins.size();
    prot_out.write(reinterpret_cast<const char*>(&num_proteins), sizeof(uint32_t));
    
    // Write all sequences concatenated
    for (const auto& protein : proteins) {
        prot_out.write(protein.sequence.c_str(), protein.sequence.length());
    }
    prot_out.close();
    
    // Write metadata JSON
    std::string metadata_file = output_dir + "/protein_details.json";
    std::ofstream meta_out(metadata_file);
    
    meta_out << "[" << std::endl;
    size_t offset = 0;
    for (size_t i = 0; i < proteins.size(); i++) {
        const auto& p = proteins[i];
        meta_out << "  {" << std::endl;
        meta_out << "    \"id\": " << i << "," << std::endl;
        meta_out << "    \"accession\": \"" << p.accession << "\"," << std::endl;
        meta_out << "    \"gene_name\": \"" << p.gene_name << "\"," << std::endl;
        meta_out << "    \"gene_family\": \"" << p.gene_family << "\"," << std::endl;
        meta_out << "    \"drug_class\": \"" << p.drug_class << "\"," << std::endl;
        meta_out << "    \"gene_id\": " << p.gene_id << "," << std::endl;
        meta_out << "    \"length\": " << p.length << "," << std::endl;
        meta_out << "    \"offset\": " << offset << std::endl;
        meta_out << "  }" << (i < proteins.size() - 1 ? "," : "") << std::endl;
        offset += p.length;
    }
    meta_out << "]" << std::endl;
    meta_out.close();
    
    // Write k-mer index (format expected by TranslatedSearchEngine)
    std::string kmer_file = output_dir + "/protein_kmers.bin";
    std::ofstream kmer_out(kmer_file, std::ios::binary);
    
    uint32_t kmer_size = KMER_SIZE;
    uint32_t num_kmers = kmer_to_positions.size();
    kmer_out.write(reinterpret_cast<const char*>(&kmer_size), sizeof(uint32_t));
    kmer_out.write(reinterpret_cast<const char*>(&num_kmers), sizeof(uint32_t));
    
    // Write k-mers in sorted order (map is already sorted)
    for (const auto& [kmer, positions] : kmer_to_positions) {
        // Write k-mer sequence
        kmer_out.write(kmer.c_str(), KMER_SIZE);
        
        // Write number of positions
        uint32_t num_positions = positions.size();
        kmer_out.write(reinterpret_cast<const char*>(&num_positions), sizeof(uint32_t));
        
        // Write positions
        for (const auto& [protein_idx, pos] : positions) {
            kmer_out.write(reinterpret_cast<const char*>(&protein_idx), sizeof(uint32_t));
            kmer_out.write(reinterpret_cast<const char*>(&pos), sizeof(uint32_t));
        }
    }
    kmer_out.close();
    
    std::cout << "\nDatabase files created:" << std::endl;
    std::cout << "  " << protein_file << " (" << num_proteins << " proteins)" << std::endl;
    std::cout << "  " << metadata_file << " (metadata for clinical reporting)" << std::endl;
    std::cout << "  " << kmer_file << " (" << num_kmers << " k-mers indexed)" << std::endl;
    
    // Create a summary file for verification
    std::string summary_file = output_dir + "/database_summary.txt";
    std::ofstream summary(summary_file);
    summary << "AMR Protein Database Summary" << std::endl;
    summary << "===========================" << std::endl;
    summary << "Total proteins: " << proteins.size() << std::endl;
    summary << "Total k-mers: " << kmer_to_positions.size() << std::endl;
    summary << "K-mer size: " << KMER_SIZE << std::endl;
    summary << "\nDrug classes:" << std::endl;
    for (const auto& [drug_class, count] : drug_class_counts) {
        summary << "  " << drug_class << ": " << count << std::endl;
    }
    summary << "\nExample proteins:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(10), proteins.size()); i++) {
        summary << "  " << proteins[i].gene_name << " (" << proteins[i].drug_class 
                << "): " << proteins[i].length << " aa" << std::endl;
    }
    summary.close();
    
    std::cout << "\nDatabase build complete!" << std::endl;
    return 0;
}