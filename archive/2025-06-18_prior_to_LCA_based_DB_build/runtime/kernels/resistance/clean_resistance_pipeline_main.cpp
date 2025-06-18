// clean_resistance_pipeline_main.cpp
// Clean integrated pipeline for GPU-accelerated fluoroquinolone resistance detection
// Uses existing CUDA kernels for nucleotide k-mer matching and translated protein search
// Now with configurable bloom filter and Smith-Waterman alignment
// INCLUDES: Batch processing, hardcoded FQ mutations, AND fixed allele frequency tracking

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include <cuda_runtime.h>
#include <zlib.h>
#include <unordered_map>
#include <climits>

#include "fq_mutation_detector.cuh"
#include "hdf5_alignment_writer.h"
#include "global_fq_resistance_mapper.h"
#include "fq_resistance_positions.h"
#include "fq_mutations_hardcoded.h"
#include "sample_csv_parser.h"

// External global FQ resistance database
extern FQResistanceDatabase* g_fq_resistance_db;

// FQ mutation reporter declarations
extern "C" {
    void init_fq_resistance_database(const char* csv_path);
    void cleanup_fq_resistance_database();
    void* create_fq_mutation_reporter(const char* output_path);
    void destroy_fq_mutation_reporter(void* reporter);
    void set_fq_reporter_gene_mapping(void* reporter, uint32_t id, const char* name);
    void set_fq_reporter_species_mapping(void* reporter, uint32_t id, const char* name);
    void process_protein_match_for_fq_report(void* reporter, const void* match);
    void generate_fq_mutation_report(void* reporter);
}

// Clinical report generator declarations
extern "C" {
    void* create_clinical_fq_report_generator(const char* output_path);
    void destroy_clinical_report_generator(void* generator);
    void process_match_for_clinical_report(void* generator, uint32_t read_id,
                                         const char* species_name, const char* gene_name,
                                         uint32_t gene_id, uint32_t species_id,
                                         int num_mutations, const int* positions,
                                         const char* ref_aas, const char* mut_aas,
                                         float alignment_score, float identity,
                                         int match_length, int ref_start, int ref_end);
    void update_clinical_report_stats(void* generator, int total_reads, int reads_with_matches);
    void update_clinical_report_performance(void* generator, double processing_seconds, double reads_per_sec);
    void generate_clinical_report(void* generator);
    void add_allele_frequency_to_clinical_report(void* generator, 
                                                  const char* species, const char* gene,
                                                  uint16_t position, uint32_t total_depth,
                                                  char wildtype_aa, uint32_t wildtype_count, 
                                                  float wildtype_frequency,
                                                  char dominant_mutant_aa, uint32_t dominant_mutant_count, 
                                                  float dominant_mutant_frequency, 
                                                  float total_resistant_frequency,
                                                  bool has_resistance_mutation, 
                                                  const char* mutation_summary);
}

// External CUDA kernel declarations
extern "C" {
    // Bloom filter functions
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int bloom_filter_screen_reads_with_rc(void* filter, const char* d_reads, const int* d_read_lengths,
                                          const int* d_read_offsets, int num_reads, bool* d_read_passes,
                                          int* d_kmers_found, int min_kmers_threshold, bool check_rc,
                                          int* d_debug_stats);
    
    // Translated search functions
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
    void set_smith_waterman_enabled(void* engine, bool enabled);
    int search_translated_reads(void* engine, const char* d_reads, const int* d_read_lengths,
                                const int* d_read_offsets, const bool* d_reads_to_process,
                                int num_reads, void* results, uint32_t* result_counts);
    
    // K-mer filtering functions
    void launch_kmer_filter(const char* d_reads, const int* d_read_lengths, const int* d_read_offsets,
                           FQMutationDetectorCUDA& detector, int num_reads, CandidateMatch* d_candidates,
                           uint32_t* d_candidate_counts, bool check_reverse_complement);
    
    void launch_position_weighted_alignment(const char* d_reads, const int* d_read_lengths,
                                          const int* d_read_offsets, const CandidateMatch* d_candidates,
                                          const uint32_t* d_candidate_counts, FQMutationDetectorCUDA& detector,
                                          int num_reads, AlignmentResult* d_results, uint32_t* d_result_count);
}

// Structure definitions (matching your existing code)
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    uint8_t num_mutations;
    uint8_t mutation_positions[10];
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
    bool used_smith_waterman;
    char query_peptide[51];
    bool is_qrdr_alignment;
};

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

// Allele frequency tracking structure
struct AlleleFrequencyData {
    uint32_t species_id;
    uint32_t gene_id;
    std::string species_name;
    std::string gene_name;
    uint16_t position;           // 1-based position
    char wildtype_aa;
    std::map<char, uint32_t> aa_counts;  // amino acid -> count
    uint32_t total_depth;
    float wildtype_frequency;
    float dominant_mutant_frequency;
    char dominant_mutant_aa;
    std::vector<char> known_resistant_aas;
    float total_resistant_frequency;
    bool has_resistance_mutation;
    std::string mutation_summary;  // e.g., "S83L(45%),S83F(12%)"
};

// CUDA error checking
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "[CUDA ERROR] %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

class CleanResistancePipeline {
private:
    // FIXED Allele Frequency Analyzer class with proper wildtype tracking
    class AlleleFrequencyAnalyzer {
    private:
        // Pileup data structure: [species_gene_key][position][amino_acid] = count
        std::unordered_map<std::string, std::unordered_map<uint16_t, std::unordered_map<char, uint32_t>>> pileup_data;
        std::unordered_map<std::string, std::unordered_map<uint16_t, uint32_t>> depth_data;
        
        // FIXED: Store reference amino acids observed from actual alignments
        // [species_gene_key][position] = reference_amino_acid
        std::unordered_map<std::string, std::unordered_map<uint16_t, char>> reference_aas;
        
        const std::map<uint32_t, std::string>& gene_id_to_name;
        const std::map<uint32_t, std::string>& species_id_to_name;
        
    public:
        AlleleFrequencyAnalyzer(const std::map<uint32_t, std::string>& genes,
                               const std::map<uint32_t, std::string>& species)
            : gene_id_to_name(genes), species_id_to_name(species) {}
        
        void addProteinAlignment(const ProteinMatch& match) {
            // Get species and gene names
            std::string species_name = species_id_to_name.count(match.species_id) ? 
                                      species_id_to_name.at(match.species_id) : 
                                      "Species" + std::to_string(match.species_id);
            std::string gene_name = gene_id_to_name.count(match.gene_id) ? 
                                   gene_id_to_name.at(match.gene_id) : 
                                   "Gene" + std::to_string(match.gene_id);
            
            std::string species_gene_key = species_name + "_" + gene_name;
            
            // FIXED: First, store reference amino acids from mutations
            // This gives us the wildtype at mutated positions
            for (uint8_t i = 0; i < match.num_mutations && i < 10; i++) {
                uint16_t global_pos = match.ref_start + match.mutation_positions[i] + 1; // Convert to 1-based
                char ref_aa = match.ref_aas[i];
                
                // Store the reference amino acid if we haven't seen it before
                if (reference_aas[species_gene_key].find(global_pos) == reference_aas[species_gene_key].end()) {
                    reference_aas[species_gene_key][global_pos] = ref_aa;
                }
            }
            
            // FIXED: Process ALL positions in the alignment using query_peptide
            // This ensures we count both wildtype and mutant amino acids correctly
            for (uint16_t i = 0; i < match.match_length && i < 50; i++) { // query_peptide max length is 50
                char query_aa = match.query_peptide[i];
                if (query_aa == '\0') break; // End of peptide
                
                uint16_t global_pos = match.ref_start + i + 1; // Convert to 1-based
                
                // Add the observed amino acid to pileup
                if (query_aa != 'X' && query_aa != '*') { // Skip unknown/stop codons
                    pileup_data[species_gene_key][global_pos][query_aa]++;
                    depth_data[species_gene_key][global_pos]++;
                }
            }
        }
        
        std::vector<AlleleFrequencyData> generateAlleleFrequencies(uint32_t min_depth = 5) {
            std::vector<AlleleFrequencyData> results;
            
            for (const auto& species_gene_pair : pileup_data) {
                const std::string& species_gene_key = species_gene_pair.first;
                const auto& positions = species_gene_pair.second;
                
                // Parse species and gene names
                size_t underscore_pos = species_gene_key.find_last_of('_');
                if (underscore_pos == std::string::npos) continue;
                
                std::string species_name = species_gene_key.substr(0, underscore_pos);
                std::string gene_name = species_gene_key.substr(underscore_pos + 1);
                
                // Get species and gene IDs
                uint32_t species_id = getSpeciesId(species_name);
                uint32_t gene_id = getGeneId(gene_name);
                
                for (const auto& pos_pair : positions) {
                    uint16_t position = pos_pair.first;
                    const auto& aa_counts = pos_pair.second;
                    
                    uint32_t total_depth = depth_data[species_gene_key][position];
                    if (total_depth < min_depth) continue;
                    
                    // NEW: Only report known FQ resistance positions
                    if (!isKnownFQResistancePosition(species_name, gene_name, position)) {
                        continue;
                    }
                    
                    AlleleFrequencyData freq_data;
                    freq_data.species_id = species_id;
                    freq_data.gene_id = gene_id;
                    freq_data.species_name = species_name;
                    freq_data.gene_name = gene_name;
                    freq_data.position = position;
                    // Convert unordered_map to map
                    for (const auto& pair : aa_counts) {
                        freq_data.aa_counts[pair.first] = pair.second;
                    }
                    freq_data.total_depth = total_depth;
                    
                    // Get wildtype from hardcoded FQ mutations or observed references
                    freq_data.wildtype_aa = getWildtypeFromFQMutations(species_name, gene_name, position);
                    if (freq_data.wildtype_aa == 'X' && reference_aas[species_gene_key].find(position) != reference_aas[species_gene_key].end()) {
                        freq_data.wildtype_aa = reference_aas[species_gene_key][position];
                    }
                    freq_data.known_resistant_aas = getKnownResistantAAs(species_name, gene_name, position);
                    
                    // Calculate frequencies
                    uint32_t wt_count = aa_counts.count(freq_data.wildtype_aa) ? 
                                       aa_counts.at(freq_data.wildtype_aa) : 0;
                    freq_data.wildtype_frequency = (float)wt_count / total_depth;
                    
                    // Find dominant mutant
                    uint32_t max_mutant_count = 0;
                    char dominant_mutant = 'X';
                    uint32_t total_resistant_count = 0;
                    std::vector<std::string> mutation_details;
                    
                    for (const auto& aa_pair : aa_counts) {
                        char aa = aa_pair.first;
                        uint32_t count = aa_pair.second;
                        
                        if (aa != freq_data.wildtype_aa) {
                            if (count > max_mutant_count) {
                                max_mutant_count = count;
                                dominant_mutant = aa;
                            }
                            
                            // Check if this is a known resistance mutation
                            bool is_resistant = std::find(freq_data.known_resistant_aas.begin(),
                                                         freq_data.known_resistant_aas.end(),
                                                         aa) != freq_data.known_resistant_aas.end();
                            if (is_resistant) {
                                total_resistant_count += count;
                                freq_data.has_resistance_mutation = true;
                            }
                            
                            // Add to mutation summary
                            float aa_freq = (float)count / total_depth;
                            if (aa_freq >= 0.05) {  // Only include variants ≥5%
                                std::string change = std::string(1, freq_data.wildtype_aa) + 
                                                    std::to_string(position) + aa;
                                mutation_details.push_back(change + "(" + 
                                                          std::to_string((int)(aa_freq * 100)) + "%)");
                            }
                        }
                    }
                    
                    freq_data.dominant_mutant_aa = dominant_mutant;
                    freq_data.dominant_mutant_frequency = (float)max_mutant_count / total_depth;
                    freq_data.total_resistant_frequency = (float)total_resistant_count / total_depth;
                    
                    // Create mutation summary string
                    if (!mutation_details.empty()) {
                        freq_data.mutation_summary = "";
                        for (size_t i = 0; i < mutation_details.size(); i++) {
                            if (i > 0) freq_data.mutation_summary += ",";
                            freq_data.mutation_summary += mutation_details[i];
                        }
                    } else {
                        freq_data.mutation_summary = "None";
                    }
                    
                    results.push_back(freq_data);
                }
            }
            
            // Sort by resistance frequency (descending)
            std::sort(results.begin(), results.end(),
                      [](const AlleleFrequencyData& a, const AlleleFrequencyData& b) {
                          return a.total_resistant_frequency > b.total_resistant_frequency;
                      });
            
            return results;
        }
        
    private:
        char getWildtypeAA(const std::string& species, const std::string& gene, uint16_t position) {
            // Use global FQ resistance database if available
            if (g_fq_resistance_db) {
                char wt = g_fq_resistance_db->getWildtypeAA(gene, position);
                if (wt != 'X') return wt;
            }
            
            // Fallback to known wildtype amino acids for common FQ resistance positions
            static const std::map<std::string, std::map<uint16_t, char>> wildtype_map = {
                {"gyrA", {{83, 'S'}, {87, 'D'}}},
                {"parC", {{80, 'S'}, {84, 'E'}}},
                {"gyrB", {{426, 'D'}, {447, 'K'}}},
                {"parE", {{416, 'S'}, {420, 'L'}}},
                {"grlA", {{80, 'S'}, {84, 'E'}}}  // S. aureus
            };
            
            auto gene_it = wildtype_map.find(gene);
            if (gene_it != wildtype_map.end()) {
                auto pos_it = gene_it->second.find(position);
                if (pos_it != gene_it->second.end()) {
                    return pos_it->second;
                }
            }
            
            return 'X';  // Unknown
        }
        
        std::vector<char> getKnownResistantAAs(const std::string& species, const std::string& gene, uint16_t position) {
            // Use global FQ resistance database if available
            if (g_fq_resistance_db) {
                std::vector<char> resistant_aas;
                // Query the database for known resistant amino acids at this position
                // This would require adding a method to FQResistanceDatabase
                // For now, fall back to hardcoded values
            }
            
            // Known resistant amino acids by species/gene/position
            static const std::map<std::string, std::map<std::string, std::map<uint16_t, std::vector<char>>>> resistance_map = {
                {"Escherichia_coli", {
                    {"gyrA", {{83, {'L', 'F', 'W'}}, {87, {'N', 'G', 'Y', 'H'}}}},
                    {"parC", {{80, {'I', 'R', 'F'}}, {84, {'V', 'K', 'G', 'A'}}}}
                }},
                {"Enterococcus_faecium", {
                    {"gyrA", {{83, {'R', 'I', 'Y', 'N'}}, {87, {'G', 'K'}}}},
                    {"parC", {{80, {'I', 'R'}}, {84, {'K', 'A'}}}}
                }},
                {"Staphylococcus_aureus", {
                    {"gyrA", {{84, {'L'}}, {88, {'Y'}}}},
                    {"grlA", {{80, {'F', 'Y'}}, {84, {'L', 'V'}}}}
                }}
            };
            
            auto species_it = resistance_map.find(species);
            if (species_it != resistance_map.end()) {
                auto gene_it = species_it->second.find(gene);
                if (gene_it != species_it->second.end()) {
                    auto pos_it = gene_it->second.find(position);
                    if (pos_it != gene_it->second.end()) {
                        return pos_it->second;
                    }
                }
            }
            
            return std::vector<char>();
        }
        
        uint32_t getSpeciesId(const std::string& species_name) {
            for (const auto& pair : species_id_to_name) {
                if (pair.second == species_name) return pair.first;
            }
            return UINT32_MAX;
        }
        
        uint32_t getGeneId(const std::string& gene_name) {
            for (const auto& pair : gene_id_to_name) {
                if (pair.second == gene_name) return pair.first;
            }
            return UINT32_MAX;
        }
        
        // Check if this is a known FQ resistance position based on hardcoded mutations
        bool isKnownFQResistancePosition(const std::string& species, const std::string& gene, uint16_t position) {
            // Check if this species/gene/position combination exists in our hardcoded mutations
            for (const auto& mutation : HARDCODED_FQ_MUTATIONS) {
                if (mutation.species == species && mutation.gene == gene && mutation.position == position) {
                    return true;
                }
            }
            return false;
        }
        
        // Get wildtype AA from hardcoded FQ mutations
        char getWildtypeFromFQMutations(const std::string& species, const std::string& gene, uint16_t position) {
            // Find the wildtype for this species/gene/position
            for (const auto& mutation : HARDCODED_FQ_MUTATIONS) {
                if (mutation.species == species && mutation.gene == gene && mutation.position == position) {
                    return mutation.wildtype_aa;
                }
            }
            return 'X'; // Unknown
        }
    };

    // GPU memory for reads
    char* d_reads_r1;
    char* d_reads_r2;
    int* d_lengths_r1;
    int* d_lengths_r2;
    int* d_offsets_r1;
    int* d_offsets_r2;
    
    // Bloom filter
    void* bloom_filter;
    bool* d_bloom_passes_r1;
    bool* d_bloom_passes_r2;
    int* d_bloom_kmers_r1;
    int* d_bloom_kmers_r2;
    
    // Nucleotide search results
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    AlignmentResult* d_results;
    uint32_t* d_result_count;
    
    // Protein search
    void* translated_search_engine;
    ProteinMatch* d_protein_matches;
    uint32_t* d_protein_match_counts;
    
    // Components
    FQMutationDetectorCUDA detector;
    HDF5AlignmentWriter* hdf5_writer;
    void* fq_mutation_reporter;
    void* clinical_report_generator;
    AlleleFrequencyAnalyzer* allele_analyzer;
    std::string protein_db_path;
    
    // Mappings
    std::map<uint32_t, std::string> gene_id_to_name;
    std::map<uint32_t, std::string> species_id_to_name;
    
    // Configuration
    const int batch_size = 10000;
    const int max_read_length = 300;
    const int bloom_min_kmers = 3;
    
    // Configuration flags
    bool use_bloom_filter = true;        // Enable bloom filtering by default
    bool use_smith_waterman = true;      // Enable Smith-Waterman by default
    uint32_t min_allele_depth = 5;       // Minimum depth for allele frequency analysis
    uint32_t min_report_depth = 0;       // Minimum depth for reporting (0 = no filter)
    const int kmer_length = 15;
    const int max_candidates_per_read = 64;
    const int max_protein_matches_per_read = 32;
    
    // Statistics
    struct Stats {
        int total_reads = 0;
        int bloom_passed = 0;
        int nucleotide_matches = 0;
        int protein_matches = 0;
        int resistance_mutations = 0;
        int fq_resistance_mutations = 0;
    } stats;

public:
    CleanResistancePipeline(bool enable_bloom = true, bool enable_sw = true,
                            uint32_t min_allele_d = 5, uint32_t min_report_d = 0) 
        : use_bloom_filter(enable_bloom), use_smith_waterman(enable_sw),
          min_allele_depth(min_allele_d), min_report_depth(min_report_d) {
        
        std::cout << "Initializing Clean Resistance Detection Pipeline\n";
        std::cout << "  Bloom filter: " << (use_bloom_filter ? "ENABLED" : "DISABLED") << "\n";
        std::cout << "  Smith-Waterman: " << (use_smith_waterman ? "ENABLED" : "DISABLED") << "\n";
        std::cout << "  Min allele depth: " << min_allele_depth << "\n";
        std::cout << "  Min report depth: " << min_report_depth << "\n";
        
        // Only create bloom filter if enabled
        if (use_bloom_filter) {
            bloom_filter = create_bloom_filter(kmer_length);
            if (!bloom_filter) {
                throw std::runtime_error("Failed to create bloom filter");
            }
        } else {
            bloom_filter = nullptr;
        }
        
        // Create translated search engine with configurable Smith-Waterman
        translated_search_engine = create_translated_search_engine_with_sw(batch_size, use_smith_waterman);
        if (!translated_search_engine) {
            throw std::runtime_error("Failed to create translated search engine");
        }
        
        // Allocate GPU memory
        size_t read_buffer_size = batch_size * max_read_length;
        CUDA_CHECK(cudaMalloc(&d_reads_r1, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_reads_r2, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_lengths_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_lengths_r2, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r2, batch_size * sizeof(int)));
        
        // Bloom filter results
        CUDA_CHECK(cudaMalloc(&d_bloom_passes_r1, batch_size * sizeof(bool)));
        CUDA_CHECK(cudaMalloc(&d_bloom_passes_r2, batch_size * sizeof(bool)));
        CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r2, batch_size * sizeof(int)));
        
        // Nucleotide search results
        CUDA_CHECK(cudaMalloc(&d_candidates, batch_size * max_candidates_per_read * sizeof(CandidateMatch)));
        CUDA_CHECK(cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_results, batch_size * max_candidates_per_read * sizeof(AlignmentResult)));
        CUDA_CHECK(cudaMalloc(&d_result_count, sizeof(uint32_t)));
        
        // Protein search results
        CUDA_CHECK(cudaMalloc(&d_protein_matches, batch_size * max_protein_matches_per_read * sizeof(ProteinMatch)));
        CUDA_CHECK(cudaMalloc(&d_protein_match_counts, batch_size * sizeof(uint32_t)));
        
        hdf5_writer = nullptr;
        fq_mutation_reporter = nullptr;
        clinical_report_generator = nullptr;
        allele_analyzer = nullptr;
    }
    
    ~CleanResistancePipeline() {
        // DISABLED: FQ mutation reporter
        // if (fq_mutation_reporter) {
        //     generate_fq_mutation_report(fq_mutation_reporter);
        //     destroy_fq_mutation_reporter(fq_mutation_reporter);
        // }
        
        // Clean up FQ resistance database
        cleanup_fq_resistance_database();
        
        // Clean up
        if (bloom_filter) destroy_bloom_filter(bloom_filter);
        if (translated_search_engine) destroy_translated_search_engine(translated_search_engine);
        if (hdf5_writer) delete hdf5_writer;
        if (allele_analyzer) delete allele_analyzer;
        
        // Free GPU memory
        cudaFree(d_reads_r1);
        cudaFree(d_reads_r2);
        cudaFree(d_lengths_r1);
        cudaFree(d_lengths_r2);
        cudaFree(d_offsets_r1);
        cudaFree(d_offsets_r2);
        cudaFree(d_bloom_passes_r1);
        cudaFree(d_bloom_passes_r2);
        cudaFree(d_bloom_kmers_r1);
        cudaFree(d_bloom_kmers_r2);
        cudaFree(d_candidates);
        cudaFree(d_candidate_counts);
        cudaFree(d_results);
        cudaFree(d_result_count);
        cudaFree(d_protein_matches);
        cudaFree(d_protein_match_counts);
    }
    
    void setBloomFilterEnabled(bool enabled) {
        if (enabled && !bloom_filter) {
            bloom_filter = create_bloom_filter(kmer_length);
            if (!bloom_filter) {
                std::cerr << "Warning: Failed to create bloom filter\n";
                use_bloom_filter = false;
                return;
            }
        } else if (!enabled && bloom_filter) {
            destroy_bloom_filter(bloom_filter);
            bloom_filter = nullptr;
        }
        use_bloom_filter = enabled;
        std::cout << "Bloom filter " << (enabled ? "ENABLED" : "DISABLED") << "\n";
    }
    
    void setSmithWatermanEnabled(bool enabled) {
        use_smith_waterman = enabled;
        if (translated_search_engine) {
            set_smith_waterman_enabled(translated_search_engine, enabled);
        }
        std::cout << "Smith-Waterman " << (enabled ? "ENABLED" : "DISABLED") << "\n";
    }
    
    void loadDatabases(const std::string& nucleotide_index_path, 
                      const std::string& protein_db_path,
                      const std::string& fq_csv_path) {
        std::cout << "Loading databases...\n";
        
        // Initialize FQ resistance database (if path provided)
        if (!fq_csv_path.empty()) {
            init_fq_resistance_database(fq_csv_path.c_str());
        }
        
        // Initialize global FQ mapper (will use hardcoded data if no CSV provided)
        GlobalFQResistanceMapper& mapper = GlobalFQResistanceMapper::getInstance();
        const char* csv_path_ptr = fq_csv_path.empty() ? nullptr : fq_csv_path.c_str();
        if (init_global_fq_mapper(csv_path_ptr, protein_db_path.c_str()) != 0) {
            std::cerr << "WARNING: Failed to initialize global FQ resistance mapper\n";
        } else {
            std::cout << "Global FQ resistance mapper initialized successfully\n";
            mapper.printSummary();
        }
        
        // Store protein DB path for later use
        this->protein_db_path = protein_db_path;
        
        // Load nucleotide k-mer index
        detector.loadIndex(nucleotide_index_path.c_str());
        
        // Build bloom filter from k-mer index (only if enabled)
        if (use_bloom_filter && bloom_filter && detector.d_kmer_sorted && detector.num_kmers > 0) {
            if (build_bloom_filter_from_index(bloom_filter, detector.d_kmer_sorted, 
                                            detector.num_kmers) == 0) {
                std::cout << "Bloom filter built with " << detector.num_kmers << " k-mers\n";
            }
        }
        
        // Load protein database
        if (load_protein_database(translated_search_engine, protein_db_path.c_str()) == 0) {
            std::cout << "Protein database loaded from " << protein_db_path << "\n";
        } else {
            std::cerr << "Warning: Failed to load protein database\n";
        }
    }
    
    void loadMappingsForReporter() {
        // Always load mappings since we need them for the main pipeline too
        if (protein_db_path.empty()) return;
        
        // Load metadata from protein database
        std::string metadata_path = protein_db_path + "/metadata.json";
        std::ifstream metadata_file(metadata_path);
        
        if (!metadata_file.good()) {
            std::cerr << "Warning: Could not load metadata from " << metadata_path << std::endl;
            return;
        }
        
        gene_id_to_name.clear();
        species_id_to_name.clear();
        
        std::string line;
        bool in_species_map = false;
        bool in_gene_map = false;
        
        while (std::getline(metadata_file, line)) {
            // Check for section markers
            if (line.find("\"species_map\"") != std::string::npos) {
                in_species_map = true;
                in_gene_map = false;
                continue;
            }
            if (line.find("\"gene_map\"") != std::string::npos) {
                in_gene_map = true;
                in_species_map = false;
                continue;
            }
            
            // Parse "id": "name" entries
            if ((in_species_map || in_gene_map) && line.find(":") != std::string::npos) {
                size_t first_quote = line.find("\"");
                size_t second_quote = line.find("\"", first_quote + 1);
                size_t third_quote = line.find("\"", second_quote + 1);
                size_t fourth_quote = line.find("\"", third_quote + 1);
                
                if (first_quote != std::string::npos && fourth_quote != std::string::npos) {
                    std::string id_str = line.substr(first_quote + 1, second_quote - first_quote - 1);
                    std::string name = line.substr(third_quote + 1, fourth_quote - third_quote - 1);
                    
                    try {
                        uint32_t id = std::stoi(id_str);
                        if (in_species_map) {
                            species_id_to_name[id] = name;
                            // Skip setting to fq_mutation_reporter since it's disabled
                            // set_fq_reporter_species_mapping(fq_mutation_reporter, id, name.c_str());
                        } else if (in_gene_map) {
                            gene_id_to_name[id] = name;
                            // Skip setting to fq_mutation_reporter since it's disabled
                            // set_fq_reporter_gene_mapping(fq_mutation_reporter, id, name.c_str());
                        }
                    } catch (...) {
                        // Skip invalid entries
                    }
                }
            }
            
            // Check for section end
            if ((in_species_map || in_gene_map) && line.find("}") != std::string::npos) {
                in_species_map = false;
                in_gene_map = false;
            }
        }
        
        metadata_file.close();
        std::cout << "Loaded " << gene_id_to_name.size() << " gene mappings and " 
                  << species_id_to_name.size() << " species mappings\n";
        
        // Initialize allele frequency analyzer with loaded mappings
        if (allele_analyzer) {
            delete allele_analyzer;
        }
        allele_analyzer = new AlleleFrequencyAnalyzer(gene_id_to_name, species_id_to_name);
    }
    
    void processReads(const std::string& r1_path, const std::string& r2_path,
                     const std::string& output_path) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Initialize HDF5 output
        std::string hdf5_path = output_path + ".h5";
        hdf5_writer = new HDF5AlignmentWriter(hdf5_path);
        hdf5_writer->initialize("", r1_path, r2_path);
        
        // DISABLED: FQ mutation reporter to prevent interference with downstream analysis
        // fq_mutation_reporter = create_fq_mutation_reporter(output_path.c_str());
        fq_mutation_reporter = nullptr;
        
        // Initialize clinical report generator
        clinical_report_generator = create_clinical_fq_report_generator(output_path.c_str());
        
        // Load and set mappings from protein database metadata
        loadMappingsForReporter();
        
        // Open FASTQ files
        gzFile gz_r1 = gzopen(r1_path.c_str(), "r");
        gzFile gz_r2 = gzopen(r2_path.c_str(), "r");
        
        if (!gz_r1 || !gz_r2) {
            throw std::runtime_error("Failed to open input files");
        }
        
        // Open output JSON
        std::ofstream json_output(output_path + ".json");
        json_output << "{\n";
        json_output << "  \"pipeline\": \"clean_resistance_detection\",\n";
        json_output << "  \"version\": \"1.0\",\n";
        json_output << "  \"configuration\": {\n";
        json_output << "    \"bloom_filter\": " << (use_bloom_filter ? "true" : "false") << ",\n";
        json_output << "    \"smith_waterman\": " << (use_smith_waterman ? "true" : "false") << "\n";
        json_output << "  },\n";
        json_output << "  \"results\": [\n";
        
        bool first_result = true;
        
        // Open CSV for all protein matches
        std::ofstream csv_output(output_path + "_protein_matches.csv");
        csv_output << "read_id,read_number,species_name,gene_name,protein_id,gene_id,species_id,"
                   << "frame,query_start,query_end,ref_start,ref_end,match_length,"
                   << "alignment_score,identity,num_mutations,used_smith_waterman,"
                   << "query_peptide,mutations,is_qrdr_alignment\n";
        
        // Process in batches
        while (true) {
            std::vector<FastqRecord> batch_r1, batch_r2;
            
            if (!readBatch(gz_r1, gz_r2, batch_r1, batch_r2, batch_size)) {
                break;
            }
            
            if (batch_r1.empty()) break;
            
            // Process batch
            processBatch(batch_r1, batch_r2, json_output, csv_output, first_result);
            
            // Progress update
            if (stats.total_reads % 100000 == 0) {
                std::cout << "Processed " << stats.total_reads << " reads...\n";
            }
        }
        
        // Close files
        gzclose(gz_r1);
        gzclose(gz_r2);
        
        // Generate allele frequency analysis before finalizing JSON
        std::string allele_frequency_json = "";
        std::vector<AlleleFrequencyData> allele_frequencies;  // Declare outside if block
        if (allele_analyzer) {
            allele_frequencies = allele_analyzer->generateAlleleFrequencies(min_allele_depth);
            
            // Write allele frequencies to CSV
            std::string allele_csv_path = output_path + "_allele_frequencies.csv";
            
            // Extract sample name from output path
            std::string sample_name = output_path;
            size_t last_slash = sample_name.find_last_of("/");
            if (last_slash != std::string::npos) {
                sample_name = sample_name.substr(last_slash + 1);
            }
            
            writeAlleleFrequenciesToCSV(allele_frequencies, allele_csv_path, sample_name, min_report_depth);
            
            // Prepare allele frequency summary for JSON
            int total_resistance_positions = 0;
            int positions_with_resistance = 0;
            float max_resistance_freq = 0.0f;
            
            for (const auto& freq : allele_frequencies) {
                if (freq.total_depth >= 5) {  // Only count positions with sufficient depth
                    total_resistance_positions++;
                    if (freq.has_resistance_mutation) {
                        positions_with_resistance++;
                    }
                    if (freq.total_resistant_frequency > max_resistance_freq) {
                        max_resistance_freq = freq.total_resistant_frequency;
                    }
                }
            }
            
            // Build allele frequency JSON section
            std::stringstream af_json;
            af_json << "    \"allele_frequency_analysis\": {\n";
            af_json << "      \"total_positions_analyzed\": " << total_resistance_positions << ",\n";
            af_json << "      \"positions_with_resistance\": " << positions_with_resistance << ",\n";
            af_json << "      \"max_resistance_frequency\": " << std::fixed << std::setprecision(4) << max_resistance_freq << ",\n";
            af_json << "      \"allele_frequency_csv\": \"" << sample_name << "_allele_frequencies.csv\"\n";
            af_json << "    },\n";
            allele_frequency_json = af_json.str();
        }
        
        // Finalize output
        json_output << "\n  ],\n";
        json_output << "  \"summary\": {\n";
        json_output << "    \"total_reads\": " << stats.total_reads << ",\n";
        json_output << "    \"bloom_passed\": " << stats.bloom_passed << ",\n";
        json_output << "    \"nucleotide_matches\": " << stats.nucleotide_matches << ",\n";
        json_output << "    \"protein_matches\": " << stats.protein_matches << ",\n";
        json_output << "    \"resistance_mutations\": " << stats.resistance_mutations << ",\n";
        json_output << "    \"fq_resistance_mutations\": " << stats.fq_resistance_mutations << ",\n";
        json_output << allele_frequency_json;
        json_output << "    \"processing_completed\": true\n";
        json_output << "  }\n";
        json_output << "}\n";
        json_output.close();
        
        csv_output.close();
        
        hdf5_writer->finalize(output_path + "_summary.json");
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Calculate reads per second
        double reads_per_second = stats.total_reads / static_cast<double>(duration.count());
        
        // Update clinical report with final statistics and performance
        if (clinical_report_generator) {
            update_clinical_report_stats(clinical_report_generator, stats.total_reads, stats.protein_matches);
            update_clinical_report_performance(clinical_report_generator, static_cast<double>(duration.count()), reads_per_second);
            
            // Add allele frequency data to clinical report
            if (allele_analyzer) {
                for (const auto& freq : allele_frequencies) {
                    // Calculate counts
                    uint32_t wt_count = freq.aa_counts.count(freq.wildtype_aa) ? 
                                       freq.aa_counts.at(freq.wildtype_aa) : 0;
                    uint32_t dominant_count = freq.aa_counts.count(freq.dominant_mutant_aa) ? 
                                             freq.aa_counts.at(freq.dominant_mutant_aa) : 0;
                    
                    add_allele_frequency_to_clinical_report(
                        clinical_report_generator,
                        freq.species_name.c_str(),
                        freq.gene_name.c_str(),
                        freq.position,
                        freq.total_depth,
                        freq.wildtype_aa,
                        wt_count,
                        freq.wildtype_frequency,
                        freq.dominant_mutant_aa,
                        dominant_count,
                        freq.dominant_mutant_frequency,
                        freq.total_resistant_frequency,
                        freq.has_resistance_mutation,
                        freq.mutation_summary.c_str()
                    );
                }
            }
            
            // Generate the clinical report now that we have all the data
            generate_clinical_report(clinical_report_generator);
            destroy_clinical_report_generator(clinical_report_generator);
            clinical_report_generator = nullptr;  // Set to nullptr after destruction
        }
        
        // Print summary
        std::cout << "\n=== PROCESSING COMPLETE ===\n";
        std::cout << "Total reads: " << stats.total_reads << "\n";
        std::cout << "Bloom filter passed: " << stats.bloom_passed << " ("
                  << std::fixed << std::setprecision(1) 
                  << (100.0 * stats.bloom_passed / stats.total_reads) << "%)\n";
        std::cout << "Nucleotide matches: " << stats.nucleotide_matches << "\n";
        std::cout << "Protein matches: " << stats.protein_matches << "\n";
        std::cout << "Total mutations detected: " << stats.resistance_mutations << "\n";
        std::cout << "FQ resistance mutations: " << stats.fq_resistance_mutations << "\n";
        std::cout << "Processing time: " << duration.count() << " seconds\n";
        std::cout << "Performance: " << std::fixed << std::setprecision(0) 
                  << reads_per_second << " reads/second\n";
        std::cout << "Output: " << output_path << ".json and " << hdf5_path << "\n";
    }

private:
    void writeAlleleFrequenciesToCSV(const std::vector<AlleleFrequencyData>& allele_data,
                                    const std::string& csv_path,
                                    const std::string& sample_name,
                                    uint32_t min_report_depth = 0) {
        std::ofstream csv_file(csv_path);
        if (!csv_file.good()) {
            std::cerr << "ERROR: Failed to create allele frequency CSV file: " << csv_path << std::endl;
            return;
        }
        
        // Write CSV header
        csv_file << "sample_name,species,gene,position,total_depth,wildtype_aa,wildtype_count,"
                 << "wildtype_frequency,dominant_mutant_aa,dominant_mutant_count,dominant_mutant_frequency,"
                 << "total_resistant_frequency,has_resistance_mutation,mutation_summary,all_amino_acids\n";
        
        for (const auto& freq : allele_data) {
            // Skip positions below minimum reporting depth
            if (freq.total_depth < min_report_depth) {
                continue;
            }
            
            // Calculate wildtype count
            uint32_t wt_count = freq.aa_counts.count(freq.wildtype_aa) ? 
                               freq.aa_counts.at(freq.wildtype_aa) : 0;
            
            // Calculate dominant mutant count
            uint32_t dominant_count = freq.aa_counts.count(freq.dominant_mutant_aa) ? 
                                     freq.aa_counts.at(freq.dominant_mutant_aa) : 0;
            
            // Create amino acid summary
            std::string aa_summary = "";
            for (const auto& aa_pair : freq.aa_counts) {
                if (!aa_summary.empty()) aa_summary += ";";
                aa_summary += std::string(1, aa_pair.first) + ":" + std::to_string(aa_pair.second);
            }
            
            csv_file << sample_name << ","
                     << freq.species_name << ","
                     << freq.gene_name << ","
                     << freq.position << ","
                     << freq.total_depth << ","
                     << freq.wildtype_aa << ","
                     << wt_count << ","
                     << std::fixed << std::setprecision(4) << freq.wildtype_frequency << ","
                     << freq.dominant_mutant_aa << ","
                     << dominant_count << ","
                     << std::fixed << std::setprecision(4) << freq.dominant_mutant_frequency << ","
                     << std::fixed << std::setprecision(4) << freq.total_resistant_frequency << ","
                     << (freq.has_resistance_mutation ? "TRUE" : "FALSE") << ","
                     << freq.mutation_summary << ","
                     << aa_summary << "\n";
        }
        
        csv_file.close();
        
        // Count how many positions were actually written
        size_t positions_written = 0;
        for (const auto& freq : allele_data) {
            if (freq.total_depth >= min_report_depth) {
                positions_written++;
            }
        }
        
        std::cout << "Wrote allele frequencies for " << positions_written << " positions";
        if (min_report_depth > 0) {
            std::cout << " (filtered from " << allele_data.size() << " total, min depth=" << min_report_depth << ")";
        }
        std::cout << " to: " << csv_path << std::endl;
    }

    bool readBatch(gzFile gz_r1, gzFile gz_r2, std::vector<FastqRecord>& batch_r1,
                  std::vector<FastqRecord>& batch_r2, int max_size) {
        char buffer[1024];
        
        for (int i = 0; i < max_size; i++) {
            FastqRecord rec1, rec2;
            
            // Read R1
            if (gzgets(gz_r1, buffer, 1024) == NULL) return i > 0;
            rec1.header = std::string(buffer);
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
            rec1.sequence = std::string(buffer);
            rec1.sequence.pop_back(); // Remove newline
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false; // +
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
            rec1.quality = std::string(buffer);
            rec1.quality.pop_back();
            
            // Read R2
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.header = std::string(buffer);
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.sequence = std::string(buffer);
            rec2.sequence.pop_back();
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false; // +
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.quality = std::string(buffer);
            rec2.quality.pop_back();
            
            batch_r1.push_back(rec1);
            batch_r2.push_back(rec2);
        }
        
        return true;
    }
    
    ReadBatch prepareBatch(const std::vector<FastqRecord>& records) {
        ReadBatch batch;
        batch.num_reads = records.size();
        batch.total_bases = 0;
        
        for (const auto& rec : records) {
            batch.total_bases += rec.sequence.length();
        }
        
        batch.sequences = new char[batch.total_bases];
        batch.lengths = new int[batch.num_reads];
        batch.offsets = new int[batch.num_reads];
        
        int offset = 0;
        for (int i = 0; i < batch.num_reads; i++) {
            const std::string& seq = records[i].sequence;
            memcpy(batch.sequences + offset, seq.c_str(), seq.length());
            batch.lengths[i] = seq.length();
            batch.offsets[i] = offset;
            offset += seq.length();
        }
        
        return batch;
    }
    
    // Add validation function for protein matches
    bool validateProteinMatch(const ProteinMatch& match, 
                             const std::string& species_name,
                             const std::string& gene_name) {
        // Basic sanity checks
        
        // 1. Check alignment score and identity
        if (match.alignment_score < 30.0f || match.identity < 0.70f) {
            return false;  // Too low quality
        }
        
        // 2. Check alignment length
        if (match.match_length < 20) {
            return false;  // Too short
        }
        
        // 3. For FQ resistance genes, just check if alignment is in reasonable QRDR range
        // Removed hard-coded position requirements to allow flexibility across species
        
        // 4. Basic peptide quality check
        std::string peptide(match.query_peptide);
        
        // Check for obvious problems like all X's or all stops
        int x_count = 0;
        int stop_count = 0;
        for (char c : peptide) {
            if (c == 'X') x_count++;
            if (c == '*') stop_count++;
        }
        
        // Reject if more than 50% X's or stops
        if (x_count > peptide.length() / 2 || stop_count > peptide.length() / 2) {
            return false;
        }
        
        // 5. Species-specific validation removed - was causing issues with non-E.coli samples
        
        return true;
    }

    void processBatch(const std::vector<FastqRecord>& batch_r1,
                     const std::vector<FastqRecord>& batch_r2,
                     std::ofstream& json_output,
                     std::ofstream& csv_output,
                     bool& first_result) {
        int num_reads = batch_r1.size();
        stats.total_reads += num_reads;
        
        // Prepare batches
        ReadBatch r1_batch = prepareBatch(batch_r1);
        ReadBatch r2_batch = prepareBatch(batch_r2);
        
        // Transfer to GPU
        CUDA_CHECK(cudaMemcpy(d_reads_r1, r1_batch.sequences, r1_batch.total_bases, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_lengths_r1, r1_batch.lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_offsets_r1, r1_batch.offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        
        CUDA_CHECK(cudaMemcpy(d_reads_r2, r2_batch.sequences, r2_batch.total_bases, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_lengths_r2, r2_batch.lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_offsets_r2, r2_batch.offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        
        // Reset counters
        CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
        CUDA_CHECK(cudaMemset(d_candidate_counts, 0, num_reads * sizeof(uint32_t)));
        CUDA_CHECK(cudaMemset(d_protein_match_counts, 0, num_reads * sizeof(uint32_t)));
        
        // Stage 1: Bloom filter screening
        if (use_bloom_filter && bloom_filter) {
            int bloom_result_r1 = bloom_filter_screen_reads_with_rc(
                bloom_filter, d_reads_r1, d_lengths_r1, d_offsets_r1,
                num_reads, d_bloom_passes_r1, d_bloom_kmers_r1,
                bloom_min_kmers, false, nullptr
            );
            
            int bloom_result_r2 = bloom_filter_screen_reads_with_rc(
                bloom_filter, d_reads_r2, d_lengths_r2, d_offsets_r2,
                num_reads, d_bloom_passes_r2, d_bloom_kmers_r2,
                bloom_min_kmers, true, nullptr
            );
        
            // Get bloom results
            std::vector<uint8_t> h_bloom_r1(num_reads), h_bloom_r2(num_reads);
            CUDA_CHECK(cudaMemcpy(h_bloom_r1.data(), d_bloom_passes_r1, num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_bloom_r2.data(), d_bloom_passes_r2, num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
            
            // Count bloom passes
            for (int i = 0; i < num_reads; i++) {
                if (h_bloom_r1[i] || h_bloom_r2[i]) {
                    stats.bloom_passed++;
                }
            }
        } else {
            // If bloom filtering is disabled, mark all reads as passed
            CUDA_CHECK(cudaMemset(d_bloom_passes_r1, 1, num_reads * sizeof(bool)));
            CUDA_CHECK(cudaMemset(d_bloom_passes_r2, 1, num_reads * sizeof(bool)));
            stats.bloom_passed += num_reads;
        }
        
        // Stage 2: Nucleotide k-mer matching
        launch_kmer_filter(d_reads_r1, d_lengths_r1, d_offsets_r1,
                          detector, num_reads, d_candidates, d_candidate_counts, false);
        
        launch_position_weighted_alignment(d_reads_r1, d_lengths_r1, d_offsets_r1,
                                         d_candidates, d_candidate_counts,
                                         detector, num_reads, d_results, d_result_count);
        
        // Get nucleotide results
        uint32_t num_results;
        CUDA_CHECK(cudaMemcpy(&num_results, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        
        if (num_results > 0) {
            std::vector<AlignmentResult> h_results(num_results);
            CUDA_CHECK(cudaMemcpy(h_results.data(), d_results, 
                                 num_results * sizeof(AlignmentResult), cudaMemcpyDeviceToHost));
            
            stats.nucleotide_matches += num_results;
            
            // Add to HDF5
            hdf5_writer->addAlignmentBatch(h_results.data(), num_results, 
                                         stats.total_reads - num_reads);
        }
        
        // Stage 3: Protein search
        search_translated_reads(translated_search_engine,
                               d_reads_r1, d_lengths_r1, d_offsets_r1,
                               d_bloom_passes_r1, num_reads,
                               d_protein_matches, d_protein_match_counts);
        
        // Get protein results
        std::vector<uint32_t> h_protein_counts(num_reads);
        CUDA_CHECK(cudaMemcpy(h_protein_counts.data(), d_protein_match_counts,
                             num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
        
        int total_protein_matches = 0;
        for (int i = 0; i < num_reads; i++) {
            total_protein_matches += h_protein_counts[i];
        }
        
        if (total_protein_matches > 0) {
            std::vector<ProteinMatch> h_protein_matches(num_reads * max_protein_matches_per_read);
            CUDA_CHECK(cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                                 num_reads * max_protein_matches_per_read * sizeof(ProteinMatch),
                                 cudaMemcpyDeviceToHost));
            
            // Process protein matches with validation and FQ resistance checking
            GlobalFQResistanceMapper& fq_mapper = GlobalFQResistanceMapper::getInstance();
            
            int validated_matches = 0;
            int rejected_matches = 0;
            
            for (int i = 0; i < num_reads; i++) {
                for (uint32_t j = 0; j < h_protein_counts[i]; j++) {
                    ProteinMatch& match = h_protein_matches[i * max_protein_matches_per_read + j];
                    
                    // Get gene and species names
                    std::string gene_name = "unknown";
                    std::string species_name = "unknown";
                    
                    auto gene_it = gene_id_to_name.find(match.gene_id);
                    if (gene_it != gene_id_to_name.end()) {
                        gene_name = gene_it->second;
                    }
                    
                    auto species_it = species_id_to_name.find(match.species_id);
                    if (species_it != species_id_to_name.end()) {
                        species_name = species_it->second;
                    }
                    
                    // VALIDATE THE MATCH FIRST
                    if (!validateProteinMatch(match, species_name, gene_name)) {
                        rejected_matches++;
                        continue;  // Skip this match entirely
                    }
                    
                    validated_matches++;
                    stats.protein_matches++;
                    
                    // Add to allele frequency analysis
                    if (allele_analyzer) {
                        allele_analyzer->addProteinAlignment(match);
                    }
                    
                    // Check if this alignment covers QRDR region
                    // Use the FQ resistance database to check each mutation position
                    bool is_qrdr_alignment = false;
                    if (g_fq_resistance_db) {
                        // Check if any position in the alignment range is a QRDR position
                        for (int pos = match.ref_start; pos < match.ref_start + match.match_length; pos++) {
                            if (g_fq_resistance_db->isQRDRPosition(gene_name, pos)) {
                                is_qrdr_alignment = true;
                                break;
                            }
                        }
                        
                        // Also check specifically if any mutations are at QRDR positions
                        for (int k = 0; k < match.num_mutations; k++) {
                            int mutation_pos = match.ref_start + match.mutation_positions[k];
                            if (g_fq_resistance_db->isQRDRPosition(gene_name, mutation_pos)) {
                                is_qrdr_alignment = true;
                                match.is_qrdr_alignment = true;  // Update the match structure
                                break;
                            }
                        }
                    }
                    
                    // Override the GPU-set value with our host-side QRDR detection
                    match.is_qrdr_alignment = is_qrdr_alignment;
                    
                    // Write to CSV - ALL protein matches
                    csv_output << (stats.total_reads - num_reads + i) << ","
                              << i << ","
                              << species_name << ","
                              << gene_name << ","
                              << match.protein_id << ","
                              << match.gene_id << ","
                              << match.species_id << ","
                              << (int)match.frame << ","
                              << match.query_start << ","
                              << (match.query_start + match.match_length - 1) << ","
                              << match.ref_start << ","
                              << (match.ref_start + match.match_length - 1) << ","
                              << match.match_length << ","
                              << match.alignment_score << ","
                              << match.identity << ","
                              << (int)match.num_mutations << ","
                              << (match.used_smith_waterman ? "true" : "false") << ","
                              << "\"" << match.query_peptide << "\",\"";
                    
                    // Add mutations to CSV
                    for (int k = 0; k < match.num_mutations; k++) {
                        if (k > 0) csv_output << ";";
                        uint16_t global_position = match.ref_start + match.mutation_positions[k] + 1;
                        csv_output << match.ref_aas[k] << global_position << match.query_aas[k];
                    }
                    csv_output << "\","
                              << (match.is_qrdr_alignment ? "true" : "false") << "\n";
                    
                    // DISABLED: FQ mutation reporter was filtering out non-QRDR mutations
                    // This was interfering with downstream mutation calling
                    // All mutation reporting is now handled by the CSV output
                    // if (fq_mutation_reporter) {
                    //     process_protein_match_for_fq_report(fq_mutation_reporter, &match);
                    // }
                    
                    // Process match for clinical report
                    if (clinical_report_generator && match.num_mutations > 0) {
                        // Prepare arrays for C interface
                        int positions[10];
                        char ref_aas[10];
                        char mut_aas[10];
                        
                        for (int k = 0; k < match.num_mutations; k++) {
                            positions[k] = match.mutation_positions[k];
                            ref_aas[k] = match.ref_aas[k];
                            mut_aas[k] = match.query_aas[k];
                        }
                        
                        process_match_for_clinical_report(
                            clinical_report_generator,
                            stats.total_reads - num_reads + i,  // global read ID
                            species_name.c_str(),
                            gene_name.c_str(),
                            match.gene_id,
                            match.species_id,
                            match.num_mutations,
                            positions,
                            ref_aas,
                            mut_aas,
                            match.alignment_score,
                            match.identity,
                            match.match_length,
                            match.ref_start,
                            match.ref_start + match.match_length - 1
                        );
                    }
                    
                    // Debug output for matches with mutations
                    if (match.num_mutations > 0 && stats.protein_matches <= 100) {
                        std::cout << "DEBUG: Protein match " << stats.protein_matches << ": "
                                  << species_name << " " << gene_name
                                  << ", mutations=" << (int)match.num_mutations
                                  << ", score=" << match.alignment_score
                                  << ", identity=" << match.identity << "\n";
                    }
                    
                    // Check each mutation for FQ resistance
                    int fq_mutations_in_match = 0;
                    bool has_fq_resistance = false;
                    
                    for (int k = 0; k < match.num_mutations; k++) {
                        // Calculate global position (1-based)
                        uint16_t global_position = match.ref_start + match.mutation_positions[k] + 1;
                        char ref_aa = match.ref_aas[k];
                        char query_aa = match.query_aas[k];
                        
                        // Debug output for each mutation
                        if (stats.protein_matches <= 100) {
                            std::cout << "  Mutation " << (k+1) << ": " << ref_aa << global_position 
                                      << query_aa << " in " << gene_name;
                        }
                        
                        // Check if this is a known FQ resistance mutation
                        bool is_fq = fq_mapper.isResistanceMutation(
                            species_name, gene_name, global_position, ref_aa, query_aa
                        );
                        
                        if (is_fq) {
                            fq_mutations_in_match++;
                            has_fq_resistance = true;
                            stats.fq_resistance_mutations++;
                            
                            if (stats.protein_matches <= 100) {
                                std::cout << " -> FQ RESISTANCE!";
                            }
                        } else {
                            // Check if this position is known to have resistance mutations
                            char expected_wt = fq_mapper.getWildtypeAA(
                                species_name, gene_name, global_position
                            );
                            
                            if (expected_wt != 'X' && expected_wt != ref_aa) {
                                if (stats.protein_matches <= 100) {
                                    std::cout << " -> Note: Expected WT=" << expected_wt;
                                }
                            }
                        }
                        
                        if (stats.protein_matches <= 100) {
                            std::cout << "\n";
                        }
                    }
                    
                    // Always count mutations (not just FQ resistance)
                    stats.resistance_mutations += match.num_mutations;
                    
                    // Output FQ resistance to JSON
                    if (has_fq_resistance) {
                        if (!first_result) json_output << ",\n";
                        json_output << "    {\n";
                        json_output << "      \"type\": \"FQ_RESISTANCE\",\n";
                        json_output << "      \"read_id\": " << (stats.total_reads - num_reads + i) << ",\n";
                        json_output << "      \"gene_name\": \"" << gene_name << "\",\n";
                        json_output << "      \"species_name\": \"" << species_name << "\",\n";
                        json_output << "      \"gene_id\": " << match.gene_id << ",\n";
                        json_output << "      \"species_id\": " << match.species_id << ",\n";
                        json_output << "      \"frame\": " << (int)match.frame << ",\n";
                        json_output << "      \"alignment_score\": " << match.alignment_score << ",\n";
                        json_output << "      \"identity\": " << match.identity << ",\n";
                        json_output << "      \"num_mutations\": " << (int)match.num_mutations << ",\n";
                        json_output << "      \"num_fq_mutations\": " << fq_mutations_in_match << ",\n";
                        json_output << "      \"mutations\": [\n";
                        
                        bool first_mut = true;
                        for (int k = 0; k < match.num_mutations; k++) {
                            uint16_t global_position = match.ref_start + match.mutation_positions[k] + 1;
                            char ref_aa = match.ref_aas[k];
                            char query_aa = match.query_aas[k];
                            
                            bool is_fq = fq_mapper.isResistanceMutation(
                                species_name, gene_name, global_position, ref_aa, query_aa
                            );
                            
                            if (!first_mut) json_output << ",\n";
                            json_output << "        {\n";
                            json_output << "          \"position\": " << global_position << ",\n";
                            json_output << "          \"ref_aa\": \"" << ref_aa << "\",\n";
                            json_output << "          \"query_aa\": \"" << query_aa << "\",\n";
                            json_output << "          \"is_fq_resistance\": " << (is_fq ? "true" : "false") << "\n";
                            json_output << "        }";
                            first_mut = false;
                        }
                        
                        json_output << "\n      ],\n";
                        json_output << "      \"peptide\": \"" << match.query_peptide << "\"\n";
                        json_output << "    }";
                        first_result = false;
                    }
                }
            }
            
            // Print validation summary
            std::cout << "Protein match validation: " << validated_matches << " accepted, " 
                      << rejected_matches << " rejected" << std::endl;
            
            // Add to HDF5
            hdf5_writer->addTranslatedResults(h_protein_matches.data(), h_protein_counts.data(),
                                            num_reads, stats.total_reads - num_reads);
        }
        
        // Clean up batch memory
        delete[] r1_batch.sequences;
        delete[] r1_batch.lengths;
        delete[] r1_batch.offsets;
        delete[] r2_batch.sequences;
        delete[] r2_batch.lengths;
        delete[] r2_batch.offsets;
    }
};

// Function to process a single sample
int processSingleSample(const BioGPU::SampleInfo& sample, 
                       const std::string& output_prefix,
                       const std::string& nucleotide_index,
                       const std::string& protein_db,
                       const std::string& fq_csv_path,
                       bool use_bloom,
                       bool use_sw,
                       uint32_t min_allele_depth,
                       uint32_t min_report_depth) {
    try {
        std::cout << "\n=== Processing Sample: " << sample.sample_name << " ===" << std::endl;
        std::cout << "Directory: " << sample.file_path << std::endl;
        std::cout << "R1 file: " << sample.read1_filename << " => " << sample.read1_path << std::endl;
        if (!sample.read2_filename.empty()) {
            std::cout << "R2 file: " << sample.read2_filename << " => " << sample.read2_path << std::endl;
        }
        std::cout << "Output: " << output_prefix << std::endl;
        
        // Create and run pipeline
        CleanResistancePipeline pipeline(use_bloom, use_sw, min_allele_depth, min_report_depth);
        pipeline.loadDatabases(nucleotide_index, protein_db, fq_csv_path);
        pipeline.processReads(sample.read1_path, sample.read2_path, output_prefix);
        
        std::cout << "Sample " << sample.sample_name << " completed successfully!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR processing sample " << sample.sample_name << ": " << e.what() << std::endl;
        return 1;
    }
}

// Print usage information
void printUsage(const char* program_name) {
    std::cerr << "Usage:\n";
    std::cerr << "\nSingle sample mode:\n";
    std::cerr << "  " << program_name << " <nucleotide_index> <protein_db> <reads_r1.fastq.gz> <reads_r2.fastq.gz> [options]\n";
    std::cerr << "\nBatch mode:\n";
    std::cerr << "  " << program_name << " <nucleotide_index> <protein_db> --csv <samples.csv> [options]\n";
    std::cerr << "\nRequired:\n";
    std::cerr << "  nucleotide_index: Directory containing k-mer index\n";
    std::cerr << "  protein_db: Directory containing protein database\n";
    std::cerr << "\nFor single sample:\n";
    std::cerr << "  reads_r1.fastq.gz: Forward reads\n";
    std::cerr << "  reads_r2.fastq.gz: Reverse reads\n";
    std::cerr << "\nFor batch mode:\n";
    std::cerr << "  --csv <samples.csv>: CSV file with sample information\n";
    std::cerr << "    CSV columns: SampleName, FilePath, R1 file, R2 file (optional)\n";
    std::cerr << "\nOptions:\n";
    std::cerr << "  --output-prefix <prefix>: Output file prefix (single sample only)\n";
    std::cerr << "  --output-dir <dir>: Base output directory (default: results)\n";
    std::cerr << "  --fq-csv <path>: Path to FQ resistance mutations CSV (optional; uses hardcoded data if not provided)\n";
    std::cerr << "  --no-bloom: Disable bloom filter pre-screening\n";
    std::cerr << "  --no-sw: Disable Smith-Waterman alignment\n";
    std::cerr << "  --min-allele-depth <N>: Minimum depth for allele frequency analysis (default: 5)\n";
    std::cerr << "  --min-report-depth <N>: Minimum depth for reporting polymorphisms (default: 0)\n";
    std::cerr << "  --stop-on-error: Stop batch processing on first error (batch mode only)\n";
    std::cerr << "  --dry-run: Validate CSV and show what would be processed without running (batch mode only)\n";
    std::cerr << "\nExample CSV format:\n";
    std::cerr << "  SampleName,FilePath,R1 file,R2 file\n";
    std::cerr << "  Sample1,~/data/,sample1_R1.fastq.gz,sample1_R2.fastq.gz\n";
    std::cerr << "  Sample2,~/data/,sample2_R1.fastq.gz,sample2_R2.fastq.gz\n";
}

// Main function
int main(int argc, char** argv) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }
    
    // Required arguments
    std::string nucleotide_index = argv[1];
    std::string protein_db = argv[2];
    
    // Default values
    std::string output_dir = "results";
    std::string fq_csv_path = "";  // Empty by default - will use hardcoded data
    std::string csv_input_path = "";
    std::string single_output_prefix = "";
    bool use_bloom = true;
    bool use_sw = true;
    uint32_t min_allele_depth = 5;
    uint32_t min_report_depth = 0;
    bool stop_on_error = false;
    bool dry_run = false;
    bool batch_mode = false;
    
    // For single sample mode
    std::string r1_path = "";
    std::string r2_path = "";
    
    // Parse command line arguments
    int i = 3;
    while (i < argc) {
        std::string arg(argv[i]);
        
        if (arg == "--csv" && i + 1 < argc) {
            csv_input_path = argv[++i];
            batch_mode = true;
        } else if (arg == "--output-prefix" && i + 1 < argc) {
            single_output_prefix = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "--fq-csv" && i + 1 < argc) {
            fq_csv_path = argv[++i];
        } else if (arg == "--no-bloom") {
            use_bloom = false;
        } else if (arg == "--no-sw") {
            use_sw = false;
        } else if (arg == "--min-allele-depth" && i + 1 < argc) {
            min_allele_depth = std::stoi(argv[++i]);
        } else if (arg == "--min-report-depth" && i + 1 < argc) {
            min_report_depth = std::stoi(argv[++i]);
        } else if (arg == "--stop-on-error") {
            stop_on_error = true;
        } else if (arg == "--dry-run") {
            dry_run = true;
        } else if (arg[0] != '-' && r1_path.empty()) {
            // First non-option argument is R1 (single sample mode)
            r1_path = arg;
        } else if (arg[0] != '-' && r2_path.empty()) {
            // Second non-option argument is R2 (single sample mode)
            r2_path = arg;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
        i++;
    }
    
    // Validate arguments
    if (batch_mode && !r1_path.empty()) {
        std::cerr << "Error: Cannot specify both --csv and individual read files" << std::endl;
        printUsage(argv[0]);
        return 1;
    }
    
    if (!batch_mode && (r1_path.empty() || r2_path.empty())) {
        std::cerr << "Error: Must specify either --csv for batch mode or both R1 and R2 files" << std::endl;
        printUsage(argv[0]);
        return 1;
    }
    
    // Check CUDA device
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "Memory: " << (prop.totalGlobalMem / (1024*1024*1024)) << " GB" << std::endl;
    
    // Print configuration
    std::cout << "\n=== Clean Fluoroquinolone Resistance Detection Pipeline ===" << std::endl;
    std::cout << "Mode: " << (batch_mode ? "BATCH" : "SINGLE SAMPLE") << std::endl;
    std::cout << "Nucleotide index: " << nucleotide_index << std::endl;
    std::cout << "Protein database: " << protein_db << std::endl;
    if (!fq_csv_path.empty()) {
        std::cout << "FQ resistance CSV: " << fq_csv_path << std::endl;
    } else {
        std::cout << "FQ resistance data: Using hardcoded mutations" << std::endl;
    }
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Bloom filter: " << (use_bloom ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "  Smith-Waterman: " << (use_sw ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "  Min allele depth: " << min_allele_depth << std::endl;
    std::cout << "  Min report depth: " << min_report_depth << std::endl;
    
    if (batch_mode) {
        // Batch processing mode
        std::cout << "\nBatch input CSV: " << csv_input_path << std::endl;
        if (stop_on_error) {
            std::cout << "Stop on error: ENABLED" << std::endl;
        }
        if (dry_run) {
            std::cout << "Mode: DRY RUN (validation only)" << std::endl;
        }
        
        // Create batch processor
        BioGPU::BatchProcessor batch_processor(output_dir, true);
        
        // Load samples from CSV
        if (!batch_processor.loadSamples(csv_input_path)) {
            std::cerr << "Failed to load samples from CSV" << std::endl;
            return 1;
        }
        
        // Print summary
        batch_processor.getParser().printSummary();
        
        // Validate samples
        auto errors = batch_processor.getParser().getValidationErrors();
        if (!errors.empty()) {
            std::cerr << "\nValidation errors found:" << std::endl;
            for (const auto& error : errors) {
                std::cerr << "  " << error << std::endl;
            }
            if (stop_on_error) {
                return 1;
            }
        }
        
        if (dry_run) {
            std::cout << "\nDry run mode - showing what would be processed:" << std::endl;
            batch_processor.getParser().printDetailed();
            return 0;
        }
        
        // Process batch
        auto process_func = [&](const BioGPU::SampleInfo& sample, const std::string& output_prefix) {
            return processSingleSample(sample, output_prefix, nucleotide_index, protein_db,
                                     fq_csv_path, use_bloom, use_sw, min_allele_depth, min_report_depth);
        };
        
        int failed_count = batch_processor.processBatch(process_func, stop_on_error);
        
        return (failed_count > 0) ? 1 : 0;
        
    } else {
        // Single sample mode
        std::cout << "\nR1 reads: " << r1_path << std::endl;
        std::cout << "R2 reads: " << r2_path << std::endl;
        
        // Create sample info
        BioGPU::SampleInfo sample;
        sample.read1_path = r1_path;
        sample.read2_path = r2_path;
        
        // Extract sample name if not provided
        if (single_output_prefix.empty()) {
            size_t split_pos = r1_path.find("_R1.fastq.gz");
            if (split_pos != std::string::npos) {
                size_t dir_pos = r1_path.find_last_of("/\\");
                size_t start_pos = (dir_pos != std::string::npos) ? dir_pos + 1 : 0;
                sample.sample_name = r1_path.substr(start_pos, split_pos - start_pos);
            } else {
                size_t dir_pos = r1_path.find_last_of("/\\");
                size_t start_pos = (dir_pos != std::string::npos) ? dir_pos + 1 : 0;
                size_t dot_pos = r1_path.find_last_of(".");
                size_t end_pos = (dot_pos != std::string::npos) ? dot_pos : r1_path.length();
                sample.sample_name = r1_path.substr(start_pos, end_pos - start_pos);
            }
            
            single_output_prefix = output_dir + "/" + sample.sample_name + "/" + sample.sample_name;
            
            // Create output directory
            std::string mkdir_cmd = "mkdir -p " + output_dir + "/" + sample.sample_name;
            system(mkdir_cmd.c_str());
        } else {
            // Extract sample name from output prefix
            size_t last_slash = single_output_prefix.find_last_of("/");
            if (last_slash != std::string::npos) {
                sample.sample_name = single_output_prefix.substr(last_slash + 1);
            } else {
                sample.sample_name = single_output_prefix;
            }
        }
        
        std::cout << "Sample name: " << sample.sample_name << std::endl;
        std::cout << "Output prefix: " << single_output_prefix << std::endl;
        
        // Process single sample
        return processSingleSample(sample, single_output_prefix, nucleotide_index, protein_db,
                                 fq_csv_path, use_bloom, use_sw, min_allele_depth, min_report_depth);
    }
}