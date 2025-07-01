#ifndef NCBI_AMR_DATABASE_LOADER_H
#define NCBI_AMR_DATABASE_LOADER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <regex>
#include <cctype>
#include <cstring>
#include <cuda_runtime.h>
#include <json/json.h>

// Simplified AMR gene entry for FASTA-only data
struct AMRGeneEntry {
    // Basic identifiers (extracted from FASTA headers)
    char accession[64];          // Full accession from header
    char gene_name[64];          // Gene name if parseable from header
    char description[256];       // Full description from FASTA header
    
    // Sequence data for translated alignment
    char* dna_sequence;         // Nucleotide sequence
    uint32_t dna_length;
    char* protein_sequence;     // Amino acid sequence
    uint32_t protein_length;
    
    // Detection parameters (set to reasonable defaults)
    float identity_threshold;    // 0.90 default
    float coverage_threshold;    // 0.80 default
    
    // Optional fields (can be filled later or left empty)
    char gene_family[32];       // Extracted from gene_name if possible
    char class_[32];            // Can be inferred later
    char organisms[128];        // From header if present
    char mge_type[32];          // Mobile genetic element type
    
    // GPU memory layout optimization
    uint32_t gpu_offset_dna;    // Offset in GPU sequence array
    uint32_t gpu_offset_protein; // Offset in GPU protein array
    
    // Matching info (for linking DNA and protein sequences)
    bool has_protein_match;     // Whether we found matching protein
    uint32_t protein_index;     // Index of matching protein sequence
};

// Mobile genetic elements and associated genes
struct MGEElement {
    char element_type[32];      // transposon, integron, plasmid
    char element_name[64];      // Tn4401, In37, IncF
    char* sequence;
    uint32_t length;
    char* protein_sequence;     // If coding
    uint32_t protein_length;
    char associated_amr_genes[512]; // Commonly linked AMR genes
    
    // GPU layout
    uint32_t gpu_offset_dna;
    uint32_t gpu_offset_protein;
};

class NCBIAMRDatabaseLoader {
private:
    std::vector<AMRGeneEntry> amr_genes;
    
    // For protein sequence matching
    std::unordered_map<std::string, std::string> protein_map;
    std::unordered_map<std::string, std::string> protein_header_map;
    
    // Sequence storage for GPU transfer
    std::vector<char> concatenated_dna_sequences;
    std::vector<char> concatenated_protein_sequences;
    
    // GPU memory pointers
    char* d_dna_sequences;
    char* d_protein_sequences;
    AMRGeneEntry* d_gene_entries;
    MGEElement* d_mge_elements;
    char* d_flanking_sequences;
    
    // Additional storage
    std::vector<char> concatenated_flanking_sequences;
    std::vector<MGEElement> mge_elements;
    std::unordered_map<std::string, std::vector<std::string>> class_to_genes;

public:
    NCBIAMRDatabaseLoader() : 
        d_dna_sequences(nullptr),
        d_protein_sequences(nullptr), 
        d_gene_entries(nullptr),
        d_mge_elements(nullptr),
        d_flanking_sequences(nullptr) {}
    
    ~NCBIAMRDatabaseLoader() {
        cleanup_gpu_memory();
    }
    
    // Main loading function for your FASTA files
    bool loadFromFastaFiles(const std::string& dna_fasta_path, const std::string& protein_fasta_path) {
        std::cout << "Loading AMR database from FASTA files..." << std::endl;
        std::cout << "DNA sequences: " << dna_fasta_path << std::endl;
        std::cout << "Protein sequences: " << protein_fasta_path << std::endl;
        
        // Load DNA sequences
        if (!loadDNASequences(dna_fasta_path)) {
            std::cerr << "Failed to load DNA sequences from " << dna_fasta_path << std::endl;
            return false;
        }
        
        // Load protein sequences and match with DNA
        if (!loadProteinSequences(protein_fasta_path)) {
            std::cerr << "Failed to load protein sequences from " << protein_fasta_path << std::endl;
            return false;
        }
        
        // Match DNA and protein sequences
        matchDNAandProteinSequences();
        
        // Set reasonable detection thresholds
        setDefaultThresholds();
        
        // Organize for GPU efficiency
        organizeForGPU();
        
        // Transfer to GPU
        if (!transferToGPU()) {
            std::cerr << "Failed to transfer database to GPU" << std::endl;
            return false;
        }
        
        std::cout << "Successfully loaded " << amr_genes.size() << " AMR genes" << std::endl;
        printLoadingStats();
        
        return true;
    }
    
private:
    bool loadDNASequences(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }
        
        std::string line, current_header, current_sequence;
        while (std::getline(file, line)) {
            if (line[0] == '>') {
                // Process previous sequence
                if (!current_header.empty() && !current_sequence.empty()) {
                    processDNASequence(current_header, current_sequence);
                }
                current_header = line.substr(1); // Remove '>'
                current_sequence.clear();
            } else {
                // Remove any whitespace and convert to uppercase
                for (char c : line) {
                    if (!isspace(c)) {
                        current_sequence += toupper(c);
                    }
                }
            }
        }
        
        // Process last sequence
        if (!current_header.empty() && !current_sequence.empty()) {
            processDNASequence(current_header, current_sequence);
        }
        
        file.close();
        std::cout << "Loaded " << amr_genes.size() << " DNA sequences" << std::endl;
        return true;
    }
    
    bool loadProteinSequences(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }
        
        std::string line, current_header, current_sequence;
        std::unordered_map<std::string, std::string> protein_sequences;
        std::unordered_map<std::string, std::string> protein_headers;
        
        while (std::getline(file, line)) {
            if (line[0] == '>') {
                // Process previous sequence
                if (!current_header.empty() && !current_sequence.empty()) {
                    std::string accession = extractAccessionFromHeader(current_header);
                    protein_sequences[accession] = current_sequence;
                    protein_headers[accession] = current_header;
                }
                current_header = line.substr(1);
                current_sequence.clear();
            } else {
                // Remove whitespace and convert to uppercase
                for (char c : line) {
                    if (!isspace(c)) {
                        current_sequence += toupper(c);
                    }
                }
            }
        }
        
        // Process last sequence
        if (!current_header.empty() && !current_sequence.empty()) {
            std::string accession = extractAccessionFromHeader(current_header);
            protein_sequences[accession] = current_sequence;
            protein_headers[accession] = current_header;
        }
        
        file.close();
        
        // Store protein sequences for matching
        this->protein_map = protein_sequences;
        this->protein_header_map = protein_headers;
        
        std::cout << "Loaded " << protein_sequences.size() << " protein sequences" << std::endl;
        return true;
    }
    
    void processDNASequence(const std::string& header, const std::string& sequence) {
        AMRGeneEntry gene = {};
        
        // Extract accession and basic info from header
        std::string accession = extractAccessionFromHeader(header);
        strncpy(gene.accession, accession.c_str(), 63);
        gene.accession[63] = '\0';
        
        // Store full header as description
        strncpy(gene.description, header.c_str(), 255);
        gene.description[255] = '\0';
        
        // Try to extract gene name from header
        std::string gene_name = extractGeneNameFromHeader(header);
        strncpy(gene.gene_name, gene_name.c_str(), 63);
        gene.gene_name[63] = '\0';
        
        // Store DNA sequence
        gene.dna_length = sequence.length();
        gene.dna_sequence = new char[gene.dna_length + 1];
        strcpy(gene.dna_sequence, sequence.c_str());
        
        // Initialize other fields
        gene.identity_threshold = 0.90f;  // Default threshold
        gene.coverage_threshold = 0.80f;  // Default threshold
        gene.has_protein_match = false;
        gene.protein_sequence = nullptr;
        gene.protein_length = 0;
        
        // Try to extract gene family from gene name
        extractGeneFamilyFromName(gene_name, gene.gene_family);
        
        amr_genes.push_back(gene);
    }
    
    std::string extractAccessionFromHeader(const std::string& header) {
        // Handle various FASTA header formats
        // Common formats:
        // >accession description
        // >gi|number|ref|accession| description
        // >accession.version description
        
        std::istringstream iss(header);
        std::string first_token;
        iss >> first_token;
        
        // Handle RefSeq format: gi|number|ref|accession|
        if (first_token.find('|') != std::string::npos) {
            std::vector<std::string> parts = splitString(first_token, '|');
            if (parts.size() >= 4) {
                return parts[3]; // Extract accession part
            }
        }
        
        // Handle simple format: accession description
        return first_token;
    }
    
    std::string extractGeneNameFromHeader(const std::string& header) {
        // Try to extract gene name from various header formats
        // Look for patterns like blaKPC, vanA, mecA, etc.
        
        std::string lower_header = header;
        std::transform(lower_header.begin(), lower_header.end(), lower_header.begin(), ::tolower);
        
        // Common gene name patterns
        std::vector<std::string> patterns = {
            "bla[a-z0-9_-]+",     // Beta-lactamases
            "van[a-z]",           // Vancomycin resistance
            "mec[a-z]",           // Methicillin resistance
            "tet[a-z0-9]+",       // Tetracycline resistance
            "sul[0-9]+",          // Sulfonamide resistance
            "aac\\([0-9']+\\)-[a-z0-9]+", // Aminoglycoside resistance
            "ant\\([0-9\"]+\\)-[a-z0-9]+",
            "aph\\([0-9']+\\)-[a-z0-9]+"
        };
        
        for (const auto& pattern : patterns) {
            std::regex regex_pattern(pattern);
            std::smatch match;
            if (std::regex_search(lower_header, match, regex_pattern)) {
                return match.str();
            }
        }
        
        // If no specific pattern found, try to extract first meaningful word
        std::istringstream iss(header);
        std::string token;
        while (iss >> token) {
            // Skip common prefixes
            if (token != ">" && token.find("gi|") != 0 && token.find("ref|") != 0) {
                // Take first part before space or punctuation
                size_t end = token.find_first_of(" \t.,;:|");
                if (end != std::string::npos) {
                    token = token.substr(0, end);
                }
                if (token.length() > 2) {  // Must be at least 3 characters
                    return token;
                }
            }
        }
        
        return "unknown";
    }
    
    void extractGeneFamilyFromName(const std::string& gene_name, char* gene_family) {
        // Extract gene family from gene name
        // e.g., blaKPC-2 -> blaKPC, vanA -> van, mecA -> mec
        
        std::string name = gene_name;
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        
        // Remove numbers and special characters from the end
        std::regex family_regex("([a-z]+).*");
        std::smatch match;
        
        if (std::regex_match(name, match, family_regex)) {
            std::string family = match[1].str();
            strncpy(gene_family, family.c_str(), 31);
            gene_family[31] = '\0';
        } else {
            strcpy(gene_family, "unknown");
        }
    }
    
    void matchDNAandProteinSequences() {
        // Try to match DNA sequences with protein sequences
        int matches = 0;
        
        for (auto& gene : amr_genes) {
            std::string dna_accession(gene.accession);
            
            // Try exact match first
            auto it = protein_map.find(dna_accession);
            if (it != protein_map.end()) {
                attachProteinSequence(gene, it->second);
                matches++;
                continue;
            }
            
            // Try fuzzy matching (remove version numbers, etc.)
            std::string base_accession = removeVersionFromAccession(dna_accession);
            for (const auto& prot_pair : protein_map) {
                std::string prot_base = removeVersionFromAccession(prot_pair.first);
                if (base_accession == prot_base) {
                    attachProteinSequence(gene, prot_pair.second);
                    matches++;
                    break;
                }
            }
            
            // Try matching by gene name
            if (!gene.has_protein_match) {
                std::string gene_name(gene.gene_name);
                for (const auto& prot_pair : protein_header_map) {
                    if (prot_pair.second.find(gene_name) != std::string::npos) {
                        auto seq_it = protein_map.find(prot_pair.first);
                        if (seq_it != protein_map.end()) {
                            attachProteinSequence(gene, seq_it->second);
                            matches++;
                            break;
                        }
                    }
                }
            }
        }
        
        std::cout << "Matched " << matches << " DNA sequences with protein sequences" << std::endl;
    }
    
    void attachProteinSequence(AMRGeneEntry& gene, const std::string& protein_seq) {
        gene.protein_length = protein_seq.length();
        gene.protein_sequence = new char[gene.protein_length + 1];
        strcpy(gene.protein_sequence, protein_seq.c_str());
        gene.has_protein_match = true;
    }
    
    std::string removeVersionFromAccession(const std::string& accession) {
        // Remove version number (e.g., NP_123456.1 -> NP_123456)
        size_t dot_pos = accession.find('.');
        if (dot_pos != std::string::npos) {
            return accession.substr(0, dot_pos);
        }
        return accession;
    }
    
    void setDefaultThresholds() {
        // Set reasonable detection thresholds based on gene family
        for (auto& gene : amr_genes) {
            std::string family(gene.gene_family);
            
            if (family.find("bla") == 0) {
                // Beta-lactamases: high stringency
                gene.identity_threshold = 0.95f;
                gene.coverage_threshold = 0.90f;
                strcpy(gene.class_, "BETA_LACTAM");
            } else if (family.find("van") == 0) {
                // Vancomycin: very high stringency
                gene.identity_threshold = 0.98f;
                gene.coverage_threshold = 0.95f;
                strcpy(gene.class_, "GLYCOPEPTIDE");
            } else if (family.find("tet") == 0) {
                // Tetracycline: moderate stringency
                gene.identity_threshold = 0.92f;
                gene.coverage_threshold = 0.85f;
                strcpy(gene.class_, "TETRACYCLINE");
            } else if (family.find("aac") == 0 || family.find("ant") == 0 || family.find("aph") == 0) {
                // Aminoglycosides: moderate stringency
                gene.identity_threshold = 0.92f;
                gene.coverage_threshold = 0.85f;
                strcpy(gene.class_, "AMINOGLYCOSIDE");
            } else {
                // Default thresholds
                gene.identity_threshold = 0.90f;
                gene.coverage_threshold = 0.80f;
                strcpy(gene.class_, "UNKNOWN");
            }
        }
    }
    
    void printLoadingStats() {
        std::cout << "\n=== Database Loading Statistics ===" << std::endl;
        std::cout << "Total genes loaded: " << amr_genes.size() << std::endl;
        
        int with_protein = 0;
        std::unordered_map<std::string, int> family_counts;
        std::unordered_map<std::string, int> class_counts;
        
        for (const auto& gene : amr_genes) {
            if (gene.has_protein_match) with_protein++;
            family_counts[gene.gene_family]++;
            class_counts[gene.class_]++;
        }
        
        std::cout << "Genes with protein matches: " << with_protein << std::endl;
        std::cout << "Genes without protein matches: " << (amr_genes.size() - with_protein) << std::endl;
        
        std::cout << "\nTop gene families:" << std::endl;
        for (const auto& pair : family_counts) {
            if (pair.second > 10) {  // Only show families with >10 genes
                std::cout << "  " << pair.first << ": " << pair.second << std::endl;
            }
        }
        
        std::cout << "\nResistance classes:" << std::endl;
        for (const auto& pair : class_counts) {
            std::cout << "  " << pair.first << ": " << pair.second << std::endl;
        }
    }
    
    std::vector<std::string> splitTSVLine(const std::string& line) {
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }
        
        return fields;
    }
    
    std::vector<std::string> splitString(const std::string& str, char delimiter) {
        std::vector<std::string> tokens;
        std::stringstream ss(str);
        std::string token;
        
        while (std::getline(ss, token, delimiter)) {
            tokens.push_back(token);
        }
        
        return tokens;
    }
    
    void setDetectionThresholds(AMRGeneEntry& gene) {
        // Set thresholds based on gene family and resistance class
        std::string gene_family(gene.gene_family);
        std::string resistance_class(gene.class_);
        
        if (resistance_class == "BETA-LACTAM") {
            gene.identity_threshold = 0.95f;  // High stringency for beta-lactamases
            gene.coverage_threshold = 0.90f;
        } else if (resistance_class == "GLYCOPEPTIDE") {
            gene.identity_threshold = 0.98f;  // Very high for vancomycin resistance
            gene.coverage_threshold = 0.95f;
        } else if (resistance_class == "AMINOGLYCOSIDE") {
            gene.identity_threshold = 0.92f;  // Moderate stringency
            gene.coverage_threshold = 0.85f;
        } else {
            gene.identity_threshold = 0.90f;  // Default
            gene.coverage_threshold = 0.80f;
        }
    }
    
    void extractMGEContext(AMRGeneEntry& gene, const std::string& context_field) {
        // Parse mobile genetic element information
        if (context_field.find("plasmid") != std::string::npos) {
            strncpy(gene.mge_type, "plasmid", 31);
        } else if (context_field.find("transposon") != std::string::npos) {
            strncpy(gene.mge_type, "transposon", 31);
        } else if (context_field.find("integron") != std::string::npos) {
            strncpy(gene.mge_type, "integron", 31);
        } else {
            strncpy(gene.mge_type, "chromosomal", 31);
        }
    }
    
    void organizeForGPU() {
        // Sort genes by resistance class for coalesced memory access
        std::sort(amr_genes.begin(), amr_genes.end(),
            [](const AMRGeneEntry& a, const AMRGeneEntry& b) {
                return std::string(a.class_) < std::string(b.class_);
            });
        
        // Concatenate sequences for GPU memory efficiency
        concatenateSequencesForGPU();
        
        // Build class boundaries for GPU kernels
        buildClassBoundaries();
    }
    
    void concatenateSequencesForGPU() {
        concatenated_dna_sequences.clear();
        concatenated_protein_sequences.clear();
        
        for (auto& gene : amr_genes) {
            // DNA sequences
            gene.gpu_offset_dna = concatenated_dna_sequences.size();
            for (uint32_t i = 0; i < gene.dna_length; i++) {
                concatenated_dna_sequences.push_back(gene.dna_sequence[i]);
            }
            concatenated_dna_sequences.push_back('\0'); // Null terminator
            
            // Protein sequences (if available)
            if (gene.has_protein_match && gene.protein_sequence) {
                gene.gpu_offset_protein = concatenated_protein_sequences.size();
                for (uint32_t i = 0; i < gene.protein_length; i++) {
                    concatenated_protein_sequences.push_back(gene.protein_sequence[i]);
                }
                concatenated_protein_sequences.push_back('\0');
            } else {
                gene.gpu_offset_protein = 0; // No protein sequence
            }
        }
        
        std::cout << "Concatenated " << concatenated_dna_sequences.size() << " DNA bytes and " 
                  << concatenated_protein_sequences.size() << " protein bytes for GPU" << std::endl;
    }
    
    void buildClassBoundaries() {
        // Create index of start/end positions for each resistance class
        class_to_genes.clear();
        std::string current_class = "";
        
        for (size_t i = 0; i < amr_genes.size(); i++) {
            std::string gene_class(amr_genes[i].class_);
            if (gene_class != current_class) {
                if (!current_class.empty()) {
                    class_to_genes[current_class].push_back(std::to_string(i));
                }
                current_class = gene_class;
                class_to_genes[current_class].push_back(std::to_string(i));
            }
        }
        // Add final boundary
        if (!current_class.empty()) {
            class_to_genes[current_class].push_back(std::to_string(amr_genes.size()));
        }
    }
    
    bool transferToGPU() {
        // Allocate GPU memory for sequences
        size_t dna_size = concatenated_dna_sequences.size();
        size_t protein_size = concatenated_protein_sequences.size();
        size_t flanking_size = concatenated_flanking_sequences.size();
        
        if (dna_size > 0) {
            cudaMalloc(&d_dna_sequences, dna_size);
            cudaMemcpy(d_dna_sequences, concatenated_dna_sequences.data(), 
                      dna_size, cudaMemcpyHostToDevice);
        }
        
        if (protein_size > 0) {
            cudaMalloc(&d_protein_sequences, protein_size);
            cudaMemcpy(d_protein_sequences, concatenated_protein_sequences.data(),
                      protein_size, cudaMemcpyHostToDevice);
        }
        
        if (flanking_size > 0) {
            cudaMalloc(&d_flanking_sequences, flanking_size);
            cudaMemcpy(d_flanking_sequences, concatenated_flanking_sequences.data(),
                      flanking_size, cudaMemcpyHostToDevice);
        }
        
        // Transfer gene metadata
        size_t genes_size = amr_genes.size() * sizeof(AMRGeneEntry);
        cudaMalloc(&d_gene_entries, genes_size);
        cudaMemcpy(d_gene_entries, amr_genes.data(), genes_size, cudaMemcpyHostToDevice);
        
        // Transfer MGE elements if any
        if (!mge_elements.empty()) {
            size_t mge_size = mge_elements.size() * sizeof(MGEElement);
            cudaMalloc(&d_mge_elements, mge_size);
            cudaMemcpy(d_mge_elements, mge_elements.data(), mge_size, cudaMemcpyHostToDevice);
        }
        
        return true;
    }
    
    void cleanup_gpu_memory() {
        if (d_dna_sequences) cudaFree(d_dna_sequences);
        if (d_protein_sequences) cudaFree(d_protein_sequences);
        if (d_gene_entries) cudaFree(d_gene_entries);
        if (d_mge_elements) cudaFree(d_mge_elements);
        if (d_flanking_sequences) cudaFree(d_flanking_sequences);
    }

public:
    // Accessors for GPU data
    AMRGeneEntry* getGPUGeneEntries() { return d_gene_entries; }
    char* getGPUDNASequences() { return d_dna_sequences; }
    char* getGPUProteinSequences() { return d_protein_sequences; }
    
    uint32_t getNumGenes() const { return amr_genes.size(); }
    
    void printDatabaseStats() {
        std::cout << "\n=== NCBI AMR Database Statistics ===" << std::endl;
        std::cout << "Total AMR genes loaded: " << amr_genes.size() << std::endl;
        std::cout << "Total DNA sequence data: " << concatenated_dna_sequences.size() << " bytes" << std::endl;
        std::cout << "Total protein sequence data: " << concatenated_protein_sequences.size() << " bytes" << std::endl;
        
        // Print resistance class distribution
        std::unordered_map<std::string, int> class_counts;
        for (const auto& gene : amr_genes) {
            class_counts[gene.class_]++;
        }
        
        std::cout << "\nResistance class distribution:" << std::endl;
        for (const auto& pair : class_counts) {
            std::cout << "  " << pair.first << ": " << pair.second << " genes" << std::endl;
        }
    }
};

#endif // NCBI_AMR_DATABASE_LOADER_H
