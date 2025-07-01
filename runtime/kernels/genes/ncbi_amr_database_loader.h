#ifndef NCBI_AMR_DATABASE_LOADER_H
#define NCBI_AMR_DATABASE_LOADER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <json/json.h>

// AMR gene entry supporting translated alignment
struct AMRGeneEntry {
    // Identifiers
    char accession[32];          // NG_047606.1
    char gene_symbol[32];        // blaKPC-2, vanA, mecA
    char gene_family[32];        // blaKPC, van, mec
    char allele[16];            // -2, -3, etc.
    
    // Classification
    char class_[32];            // BETA-LACTAM, GLYCOPEPTIDE
    char subclass[32];          // CARBAPENEM, VANCOMYCIN
    char mechanism[64];         // HYDROLYSIS, TARGET_ALTERATION
    
    // Sequence data for translated alignment
    char* dna_sequence;         // Nucleotide sequence
    uint32_t dna_length;
    char* protein_sequence;     // Amino acid sequence
    uint32_t protein_length;
    
    // Detection parameters
    float identity_threshold;    // 0.95 for strict, 0.85 for permissive
    float coverage_threshold;    // 0.90 for complete gene
    
    // Organism associations
    char organisms[256];        // Pipe-delimited species list
    float host_specificity;     // 0.0=broad, 1.0=narrow
    
    // Mobile genetic element context
    char mge_type[32];          // plasmid, transposon, integron, chromosomal
    char mge_name[64];          // Tn4401, IncF, etc.
    
    // Flanking regions (if available)
    char* upstream_flanking;    // 500bp upstream
    char* downstream_flanking;  // 500bp downstream
    uint32_t flanking_length;
    bool has_flanking;
    
    // GPU memory layout optimization
    uint32_t gpu_offset_dna;    // Offset in GPU sequence array
    uint32_t gpu_offset_protein; // Offset in GPU protein array
    uint32_t gpu_offset_flanking; // Offset in GPU flanking array
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
    std::vector<MGEElement> mge_elements;
    
    // Sequence storage for GPU transfer
    std::vector<char> concatenated_dna_sequences;
    std::vector<char> concatenated_protein_sequences;
    std::vector<char> concatenated_flanking_sequences;
    
    // Metadata mappings
    std::unordered_map<std::string, std::string> accession_to_class;
    std::unordered_map<std::string, std::vector<std::string>> class_to_genes;
    
    // GPU memory pointers
    char* d_dna_sequences;
    char* d_protein_sequences;
    char* d_flanking_sequences;
    AMRGeneEntry* d_gene_entries;
    MGEElement* d_mge_elements;
    uint32_t* d_class_boundaries;

public:
    NCBIAMRDatabaseLoader() : 
        d_dna_sequences(nullptr),
        d_protein_sequences(nullptr), 
        d_flanking_sequences(nullptr),
        d_gene_entries(nullptr),
        d_mge_elements(nullptr),
        d_class_boundaries(nullptr) {}
    
    ~NCBIAMRDatabaseLoader() {
        cleanup_gpu_memory();
    }
    
    // Main loading function
    bool loadNCBIDatabase(const std::string& base_path, bool include_flanking = false) {
        std::cout << "Loading NCBI AMR Database from: " << base_path << std::endl;
        
        // Load core NCBI files
        if (!loadAMRSequences(base_path + "/AMR_CDS.fa")) {
            std::cerr << "Failed to load AMR_CDS.fa" << std::endl;
            return false;
        }
        
        if (!loadAMRProteins(base_path + "/AMRProt.fa")) {
            std::cerr << "Failed to load AMRProt.fa" << std::endl;
            return false;
        }
        
        if (!loadAMRMetadata(base_path + "/AMR.tsv")) {
            std::cerr << "Failed to load AMR.tsv" << std::endl;
            return false;
        }
        
        // Optional: load flanking regions
        if (include_flanking) {
            loadFlankingRegions(base_path + "/flanking_regions.fa");
        }
        
        // Load mobile genetic elements
        loadMGEDatabase(base_path + "/mge_database.fa");
        
        // Organize for GPU efficiency
        organizeForGPU();
        
        // Transfer to GPU
        if (!transferToGPU()) {
            std::cerr << "Failed to transfer database to GPU" << std::endl;
            return false;
        }
        
        std::cout << "Successfully loaded " << amr_genes.size() << " AMR genes and " 
                  << mge_elements.size() << " MGE elements" << std::endl;
        
        return true;
    }
    
private:
    bool loadAMRSequences(const std::string& filename) {
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
                    processAMRSequence(current_header, current_sequence);
                }
                current_header = line.substr(1); // Remove '>'
                current_sequence.clear();
            } else {
                current_sequence += line;
            }
        }
        
        // Process last sequence
        if (!current_header.empty() && !current_sequence.empty()) {
            processAMRSequence(current_header, current_sequence);
        }
        
        file.close();
        return true;
    }
    
    bool loadAMRProteins(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }
        
        std::string line, current_header, current_sequence;
        std::unordered_map<std::string, std::string> protein_sequences;
        
        while (std::getline(file, line)) {
            if (line[0] == '>') {
                // Process previous sequence
                if (!current_header.empty() && !current_sequence.empty()) {
                    // Extract accession from header
                    std::string accession = extractAccession(current_header);
                    protein_sequences[accession] = current_sequence;
                }
                current_header = line.substr(1);
                current_sequence.clear();
            } else {
                current_sequence += line;
            }
        }
        
        // Process last sequence
        if (!current_header.empty() && !current_sequence.empty()) {
            std::string accession = extractAccession(current_header);
            protein_sequences[accession] = current_sequence;
        }
        
        // Match proteins to DNA sequences
        for (auto& gene : amr_genes) {
            auto it = protein_sequences.find(gene.accession);
            if (it != protein_sequences.end()) {
                gene.protein_length = it->second.length();
                gene.protein_sequence = new char[gene.protein_length + 1];
                strcpy(gene.protein_sequence, it->second.c_str());
            }
        }
        
        file.close();
        return true;
    }
    
    bool loadAMRMetadata(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }
        
        std::string line;
        bool header_skipped = false;
        
        while (std::getline(file, line)) {
            if (!header_skipped) {
                header_skipped = true;
                continue; // Skip header
            }
            
            std::vector<std::string> fields = splitTSVLine(line);
            if (fields.size() < 10) continue; // Ensure minimum fields
            
            // Find corresponding gene entry
            std::string accession = fields[0]; // Assuming first field is accession
            auto gene_it = std::find_if(amr_genes.begin(), amr_genes.end(),
                [&accession](const AMRGeneEntry& gene) {
                    return std::string(gene.accession) == accession;
                });
            
            if (gene_it != amr_genes.end()) {
                // Populate metadata fields
                if (fields.size() > 1) strncpy(gene_it->gene_symbol, fields[1].c_str(), 31);
                if (fields.size() > 2) strncpy(gene_it->class_, fields[2].c_str(), 31);
                if (fields.size() > 3) strncpy(gene_it->subclass, fields[3].c_str(), 31);
                if (fields.size() > 4) strncpy(gene_it->mechanism, fields[4].c_str(), 63);
                if (fields.size() > 5) strncpy(gene_it->gene_family, fields[5].c_str(), 31);
                if (fields.size() > 6) strncpy(gene_it->organisms, fields[6].c_str(), 255);
                
                // Set detection thresholds based on gene family
                setDetectionThresholds(*gene_it);
                
                // Extract MGE context if available
                if (fields.size() > 7) extractMGEContext(*gene_it, fields[7]);
            }
        }
        
        file.close();
        return true;
    }
    
    bool loadFlankingRegions(const std::string& filename) {
        // Optional: Load flanking regions for species attribution
        // This would require additional data preparation but is feasible
        std::cout << "Flanking regions loading not yet implemented" << std::endl;
        return true;
    }
    
    bool loadMGEDatabase(const std::string& filename) {
        // Load transposons, integrons, plasmid sequences
        // These help identify horizontal gene transfer context
        std::cout << "MGE database loading not yet implemented" << std::endl;
        return true;
    }
    
    void processAMRSequence(const std::string& header, const std::string& sequence) {
        AMRGeneEntry gene = {};
        
        // Extract accession from header
        std::string accession = extractAccession(header);
        strncpy(gene.accession, accession.c_str(), 31);
        
        // Store DNA sequence
        gene.dna_length = sequence.length();
        gene.dna_sequence = new char[gene.dna_length + 1];
        strcpy(gene.dna_sequence, sequence.c_str());
        
        // Initialize other fields
        gene.identity_threshold = 0.90f;  // Default threshold
        gene.coverage_threshold = 0.80f;  // Default threshold
        gene.has_flanking = false;
        gene.upstream_flanking = nullptr;
        gene.downstream_flanking = nullptr;
        gene.flanking_length = 0;
        
        amr_genes.push_back(gene);
    }
    
    std::string extractAccession(const std::string& header) {
        // Extract accession number from FASTA header
        // Format may vary: >accession description or >gi|number|ref|accession|
        size_t start = 0;
        size_t end = header.find(' ');
        if (end == std::string::npos) end = header.length();
        
        std::string accession = header.substr(start, end - start);
        
        // Handle RefSeq format: gi|number|ref|accession|
        if (accession.find('|') != std::string::npos) {
            std::vector<std::string> parts = splitString(accession, '|');
            if (parts.size() >= 4) {
                accession = parts[3]; // Extract accession part
            }
        }
        
        return accession;
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
        concatenated_flanking_sequences.clear();
        
        for (auto& gene : amr_genes) {
            // DNA sequences
            gene.gpu_offset_dna = concatenated_dna_sequences.size();
            for (uint32_t i = 0; i < gene.dna_length; i++) {
                concatenated_dna_sequences.push_back(gene.dna_sequence[i]);
            }
            concatenated_dna_sequences.push_back('\0'); // Null terminator
            
            // Protein sequences
            if (gene.protein_sequence) {
                gene.gpu_offset_protein = concatenated_protein_sequences.size();
                for (uint32_t i = 0; i < gene.protein_length; i++) {
                    concatenated_protein_sequences.push_back(gene.protein_sequence[i]);
                }
                concatenated_protein_sequences.push_back('\0');
            }
            
            // Flanking sequences (if available)
            if (gene.has_flanking && gene.upstream_flanking) {
                gene.gpu_offset_flanking = concatenated_flanking_sequences.size();
                for (uint32_t i = 0; i < gene.flanking_length; i++) {
                    concatenated_flanking_sequences.push_back(gene.upstream_flanking[i]);
                }
                concatenated_flanking_sequences.push_back('|'); // Separator
                for (uint32_t i = 0; i < gene.flanking_length; i++) {
                    concatenated_flanking_sequences.push_back(gene.downstream_flanking[i]);
                }
                concatenated_flanking_sequences.push_back('\0');
            }
        }
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
        if (d_flanking_sequences) cudaFree(d_flanking_sequences);
        if (d_gene_entries) cudaFree(d_gene_entries);
        if (d_mge_elements) cudaFree(d_mge_elements);
        if (d_class_boundaries) cudaFree(d_class_boundaries);
    }

public:
    // Accessors for GPU data
    AMRGeneEntry* getGPUGeneEntries() { return d_gene_entries; }
    char* getGPUDNASequences() { return d_dna_sequences; }
    char* getGPUProteinSequences() { return d_protein_sequences; }
    char* getGPUFlankingSequences() { return d_flanking_sequences; }
    MGEElement* getGPUMGEElements() { return d_mge_elements; }
    
    uint32_t getNumGenes() const { return amr_genes.size(); }
    uint32_t getNumMGEElements() const { return mge_elements.size(); }
    
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
