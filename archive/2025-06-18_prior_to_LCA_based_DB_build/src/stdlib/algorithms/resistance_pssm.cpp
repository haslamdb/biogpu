// resistance_pssm.cpp
// Position-Specific Scoring Matrix for resistance detection
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iomanip>

// Nucleotide encoding
const int ALPHABET_SIZE = 4;
const std::map<char, int> BASE_TO_IDX = {
    {'A', 0}, {'a', 0},
    {'C', 1}, {'c', 1},
    {'G', 2}, {'g', 2},
    {'T', 3}, {'t', 3}
};

const std::vector<char> IDX_TO_BASE = {'A', 'C', 'G', 'T'};

// Pseudocount for PSSM construction
const float PSEUDOCOUNT = 0.1;

// Genetic code for mutation annotation
std::map<std::string, char> create_codon_table() {
    return {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
}

// Utility function to normalize a sequence to uppercase
std::string normalize_sequence(const std::string& seq) {
    std::string normalized;
    for (char c : seq) {
        normalized += std::toupper(c);
    }
    return normalized;
}

// PSSM for a specific gene region
class ResistancePSSM {
private:
    std::vector<std::vector<float>> scores;  // Position x Base
    std::vector<std::vector<float>> frequencies;  // For building PSSM
    int window_size;
    int center_position;  // Position of interest (e.g., codon 83)
    std::string gene_name;
    float resistance_threshold;
    
    // Background frequencies (E. coli default)
    std::vector<float> background_freq = {0.25, 0.25, 0.25, 0.25};  // Simplified
    
public:
    // Default constructor for map insertion
    ResistancePSSM() : window_size(0), center_position(0), resistance_threshold(0.0) {}
    
    ResistancePSSM(const std::string& gene, int center, int window) 
        : gene_name(gene), center_position(center), window_size(window) {
        scores.resize(window_size, std::vector<float>(ALPHABET_SIZE, 0.0));
        frequencies.resize(window_size, std::vector<float>(ALPHABET_SIZE, PSEUDOCOUNT));
        resistance_threshold = 0.0;  // Will be calculated after training
    }
    
    // Add a training sequence (window around mutation site)
    void add_sequence(const std::string& seq, bool is_resistant) {
        if (seq.length() != window_size) {
            std::cerr << "Warning: Sequence length mismatch. Expected " 
                      << window_size << ", got " << seq.length() << "\n";
            return;
        }
        
        std::string norm_seq = normalize_sequence(seq);
        
        // Add counts to frequency matrix
        for (int i = 0; i < window_size; i++) {
            if (BASE_TO_IDX.count(norm_seq[i])) {
                int base_idx = BASE_TO_IDX.at(norm_seq[i]);
                // Weight resistant sequences more heavily
                float weight = is_resistant ? 2.0 : 1.0;
                frequencies[i][base_idx] += weight;
            }
        }
    }
    
    // Build PSSM from frequency counts
    void build_pssm() {
        // Calculate position-specific scores
        for (int pos = 0; pos < window_size; pos++) {
            float total = 0;
            for (int base = 0; base < ALPHABET_SIZE; base++) {
                total += frequencies[pos][base];
            }
            
            // Convert to log-odds scores
            for (int base = 0; base < ALPHABET_SIZE; base++) {
                float freq = frequencies[pos][base] / total;
                scores[pos][base] = std::log2(freq / background_freq[base]);
            }
        }
        
        // Calculate resistance threshold (simplified)
        resistance_threshold = calculate_threshold();
    }
    
    // Score a sequence
    float score_sequence(const std::string& seq) const {
        if (seq.length() != window_size) {
            return -std::numeric_limits<float>::infinity();
        }
        
        std::string norm_seq = normalize_sequence(seq);
        float total_score = 0.0;
        
        for (int i = 0; i < window_size; i++) {
            if (BASE_TO_IDX.count(norm_seq[i])) {
                int base_idx = BASE_TO_IDX.at(norm_seq[i]);
                total_score += scores[i][base_idx];
            }
        }
        
        return total_score;
    }
    
    // Check if sequence indicates resistance
    bool is_resistant(const std::string& seq) const {
        return score_sequence(seq) > resistance_threshold;
    }
    
    // Get conservation at each position
    std::vector<float> get_conservation() const {
        std::vector<float> conservation(window_size);
        
        for (int pos = 0; pos < window_size; pos++) {
            // Calculate information content at this position
            float info_content = 0.0;
            for (int base = 0; base < ALPHABET_SIZE; base++) {
                float freq = std::exp2(scores[pos][base]) * background_freq[base];
                if (freq > 0) {
                    info_content += freq * std::log2(freq);
                }
            }
            conservation[pos] = 2.0 + info_content;  // Max is 2 bits
        }
        
        return conservation;
    }
    
    // Print PSSM as a matrix
    void print_pssm() const {
        std::cout << "PSSM for " << gene_name << " (center: " << center_position << ")\n";
        std::cout << "Pos\tA\tC\tG\tT\tCons\n";
        
        auto conservation = get_conservation();
        
        for (int pos = 0; pos < window_size; pos++) {
            std::cout << pos << "\t";
            for (int base = 0; base < ALPHABET_SIZE; base++) {
                std::cout << std::fixed << std::setprecision(2) 
                          << scores[pos][base] << "\t";
            }
            std::cout << conservation[pos] << "\n";
        }
        
        std::cout << "Resistance threshold: " << resistance_threshold << "\n";
    }
    
private:
    float calculate_threshold() {
        // Simplified threshold calculation
        // In practice, would use ROC analysis on training data
        return window_size * 0.1;  // Placeholder
    }
};

// Collection of PSSMs for resistance genes
class ResistancePSSMCollection {
private:
    std::map<std::string, ResistancePSSM> pssms;
    std::map<std::string, char> codon_table;
    
    // Known resistance mutations for training
    struct ResistanceMutation {
        std::string gene;
        int codon_position;
        std::string wild_type_codon;
        std::string mutant_codon;
        char wild_type_aa;
        char mutant_aa;
    };
    
    std::vector<ResistanceMutation> known_mutations;
    
public:
    ResistancePSSMCollection() {
        codon_table = create_codon_table();
        initialize_fluoroquinolone_pssms();
    }
    
    void initialize_fluoroquinolone_pssms() {
        // Define key positions and their surrounding regions
        // gyrA codon 83 (typically bases 247-249 in E. coli)
        add_pssm("gyrA_83", "gyrA", 247, 21);  // 10 bases on each side
        
        // gyrA codon 87 (typically bases 259-261)
        add_pssm("gyrA_87", "gyrA", 259, 21);
        
        // parC codon 80 (typically bases 238-240)
        add_pssm("parC_80", "parC", 238, 21);
        
        // parC codon 84 (typically bases 250-252)
        add_pssm("parC_84", "parC", 250, 21);
        
        // Add known mutations for training
        add_known_mutation("gyrA", 83, "TCT", "TTG", 'S', 'L');  // S83L
        add_known_mutation("gyrA", 83, "TCT", "TTC", 'S', 'F');  // S83F
        add_known_mutation("gyrA", 87, "GAC", "AAC", 'D', 'N');  // D87N
        add_known_mutation("gyrA", 87, "GAC", "GGC", 'D', 'G');  // D87G
        add_known_mutation("parC", 80, "AGC", "ATC", 'S', 'I');  // S80I
        add_known_mutation("parC", 80, "AGC", "TTG", 'S', 'L');  // S80L
        add_known_mutation("parC", 84, "GAA", "GTA", 'E', 'V');  // E84V
        add_known_mutation("parC", 84, "GAA", "AAA", 'E', 'K');  // E84K
    }
    
    void add_pssm(const std::string& name, const std::string& gene, 
                  int center_base_position, int window) {
        pssms.emplace(name, ResistancePSSM(gene, center_base_position, window));
    }
    
    void add_known_mutation(const std::string& gene, int codon_pos,
                           const std::string& wt_codon, const std::string& mut_codon,
                           char wt_aa, char mut_aa) {
        known_mutations.push_back({gene, codon_pos, wt_codon, mut_codon, wt_aa, mut_aa});
    }
    
    // Scan a sequence for resistance mutations
    struct ResistanceHit {
        std::string pssm_name;
        float score;
        bool is_resistant;
        std::string sequence_context;
        std::string predicted_mutation;
    };
    
    std::vector<ResistanceHit> scan_sequence(const std::string& sequence, 
                                           const std::string& gene_hint = "") {
        std::vector<ResistanceHit> hits;
        
        for (auto& [name, pssm] : pssms) {
            // Skip if gene hint doesn't match
            if (!gene_hint.empty() && name.find(gene_hint) == std::string::npos) {
                continue;
            }
            
            // Scan with sliding window
            int window = 21;  // PSSM window size
            
            for (size_t i = 0; i <= sequence.length() - window; i++) {
                std::string subseq = sequence.substr(i, window);
                float score = pssm.score_sequence(subseq);
                
                if (pssm.is_resistant(subseq)) {
                    ResistanceHit hit;
                    hit.pssm_name = name;
                    hit.score = score;
                    hit.is_resistant = true;
                    hit.sequence_context = subseq;
                    
                    // Try to identify the mutation
                    hit.predicted_mutation = identify_mutation(subseq, name);
                    
                    hits.push_back(hit);
                }
            }
        }
        
        return hits;
    }
    
    // Train PSSMs with example sequences
    void train_with_examples() {
        // Add wild-type sequences
        if (pssms.find("gyrA_83") != pssms.end()) {
            pssms["gyrA_83"].add_sequence("ATCGATCGTCTGATCGATCGA", false);  // Wild-type S (TCT)
            pssms["gyrA_83"].add_sequence("ATCGATCGTTGGATCGATCGA", true);   // Mutant L (TTG)
            pssms["gyrA_83"].add_sequence("ATCGATCGTTCGATCGATCGA", true);   // Mutant F (TTC)
        }
        
        if (pssms.find("gyrA_87") != pssms.end()) {
            pssms["gyrA_87"].add_sequence("ATCGATCGGACGATCGATCGA", false);  // Wild-type D (GAC)
            pssms["gyrA_87"].add_sequence("ATCGATCGAACGATCGATCGA", true);   // Mutant N (AAC)
        }
        
        // Build all PSSMs
        for (auto& [name, pssm] : pssms) {
            pssm.build_pssm();
        }
    }
    
    // Visualize conservation across resistance sites
    void visualize_conservation() {
        std::cout << "\nConservation Analysis of Resistance Sites:\n";
        std::cout << "=========================================\n\n";
        
        for (const auto& [name, pssm] : pssms) {
            std::cout << name << ":\n";
            auto conservation = pssm.get_conservation();
            
            // Create ASCII visualization
            for (int i = 0; i < conservation.size(); i++) {
                std::cout << std::setw(3) << i << ": ";
                int bars = (int)(conservation[i] * 20);  // Scale to 0-40 chars
                for (int j = 0; j < bars; j++) {
                    std::cout << "█";
                }
                std::cout << " " << std::fixed << std::setprecision(2) 
                          << conservation[i] << "\n";
            }
            
            // Mark the center position
            std::cout << "     ";
            for (int i = 0; i < conservation.size(); i++) {
                if (i == conservation.size() / 2) {
                    std::cout << "^";
                } else if (i == conservation.size() / 2 - 1 || 
                           i == conservation.size() / 2 + 1) {
                    std::cout << "-";
                } else {
                    std::cout << " ";
                }
            }
            std::cout << " <- Critical codon\n\n";
        }
    }
    
private:
    std::string identify_mutation(const std::string& sequence, const std::string& pssm_name) {
        // Extract the center codon and identify mutation
        int center = sequence.length() / 2;
        std::string codon = sequence.substr(center - 1, 3);
        
        // Map to amino acid
        if (codon_table.count(codon)) {
            char aa = codon_table[codon];
            
            // Check against known mutations
            for (const auto& mut : known_mutations) {
                if (pssm_name.find(std::to_string(mut.codon_position)) != std::string::npos) {
                    if (codon == mut.mutant_codon) {
                        return std::string(1, mut.wild_type_aa) + 
                               std::to_string(mut.codon_position) + 
                               std::string(1, mut.mutant_aa);
                    }
                }
            }
            
            return "Unknown mutation at " + pssm_name;
        }
        
        return "Invalid codon";
    }
};

// Forward declaration
class ResistancePSSMCollection;

// GPU-optimized PSSM scanner (preparation for CUDA implementation)
class GPUPSSMScanner {
private:
    // Flattened data structure for GPU transfer
    struct GPUPSSMData {
        float* scores;          // All PSSM scores concatenated
        int* window_sizes;      // Window size for each PSSM
        float* thresholds;      // Resistance threshold for each PSSM
        int num_pssms;
        int max_window_size;
    };
    
    GPUPSSMData gpu_data;
    
public:
    void prepare_for_gpu(const ResistancePSSMCollection& collection) {
        // This would prepare data for GPU transfer
        // Flatten all PSSMs into contiguous memory
        std::cout << "Preparing PSSM data for GPU acceleration...\n";
    }
    
    // Simulated GPU kernel for batch scanning
    std::vector<std::vector<ResistancePSSMCollection::ResistanceHit>> 
    batch_scan(const std::vector<std::string>& sequences) {
        std::vector<std::vector<ResistancePSSMCollection::ResistanceHit>> results;
        
        // In CUDA, this would be:
        // __global__ void pssm_scan_kernel(char* sequences, GPUPSSMData* pssms, ...)
        
        #pragma omp parallel for
        for (size_t i = 0; i < sequences.size(); i++) {
            // Simulated GPU processing
            std::vector<ResistancePSSMCollection::ResistanceHit> seq_hits;
            results.push_back(seq_hits);
        }
        
        return results;
    }
};

// Example usage and testing
#ifdef TEST
int main() {
    // Create and train PSSM collection
    ResistancePSSMCollection pssm_collection;
    pssm_collection.train_with_examples();
    
    // Visualize conservation
    pssm_collection.visualize_conservation();
    
    // Test sequences
    std::vector<std::string> test_sequences = {
        // gyrA S83L mutation context
        "ATCGATCGTTGGATCGATCGA",  // Contains TTG (Leu) instead of TCT (Ser)
        
        // Wild-type gyrA
        "ATCGATCGTCTGATCGATCGA",  // Contains TCT (Ser)
        
        // gyrA D87N mutation context  
        "ATCGATCGAACGATCGATCGA",  // Contains AAC (Asn) instead of GAC (Asp)
    };
    
    std::cout << "\nScanning test sequences for resistance mutations:\n";
    std::cout << "================================================\n\n";
    
    for (size_t i = 0; i < test_sequences.size(); i++) {
        std::cout << "Sequence " << i << ": " << test_sequences[i] << "\n";
        
        auto hits = pssm_collection.scan_sequence(test_sequences[i]);
        
        if (hits.empty()) {
            std::cout << "  No resistance mutations detected\n";
        } else {
            for (const auto& hit : hits) {
                std::cout << "  RESISTANCE DETECTED!\n";
                std::cout << "    PSSM: " << hit.pssm_name << "\n";
                std::cout << "    Score: " << hit.score << "\n";
                std::cout << "    Mutation: " << hit.predicted_mutation << "\n";
            }
        }
        std::cout << "\n";
    }
    
    // Demonstrate GPU preparation
    GPUPSSMScanner gpu_scanner;
    gpu_scanner.prepare_for_gpu(pssm_collection);
    
    std::cout << "\nBatch processing " << test_sequences.size() 
              << " sequences (GPU simulation)...\n";
    auto batch_results = gpu_scanner.batch_scan(test_sequences);
    std::cout << "Batch processing complete.\n";
    
    return 0;
}
#endif // TEST