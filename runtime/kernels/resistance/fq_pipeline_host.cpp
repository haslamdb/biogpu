// enhanced_fq_pipeline_host.cpp
// Enhanced main host code with integrated diagnostic reporting and improved mutation detection
// Added QRDR alignments CSV export functionality
// Added optional bloom filter and k-mer matching flags

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iomanip>
#include <cstring>
#include <zlib.h>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cuda_runtime.h>
#include "fq_mutation_detector.cuh"
#include "hdf5_alignment_writer.h"

// Include Bloom filter declarations
extern "C" {
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int bloom_filter_screen_reads(void* filter, const char* d_reads, const int* d_read_lengths,
                                  const int* d_read_offsets, int num_reads, bool* d_read_passes,
                                  int* d_kmers_found, int min_kmers_threshold);
    int save_bloom_filter(void* filter, const char* filename);
    int load_bloom_filter(void* filter, const char* filename);
    
    // RC-aware functions
    int bloom_filter_screen_reads_with_rc(
        void* filter,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        bool* d_read_passes,
        int* d_kmers_found,
        int min_kmers_threshold,
        bool check_rc,
        int* d_debug_stats
    );
    
    bool bloom_filter_has_rc(void* filter);
}

// Enhanced translated search declarations
extern "C" {
    void* create_translated_search_engine(int batch_size);
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
    void set_smith_waterman_enabled(void* engine, bool enabled);
    int search_translated_reads(
        void* engine,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,
        int num_reads,
        void* results,
        uint32_t* result_counts
    );
}

// Include diagnostic reporter
extern "C" {
    void* create_diagnostic_reporter(const char* output_path);
    void destroy_diagnostic_reporter(void* reporter);
    void update_pipeline_statistics(void* reporter, int total_reads, int after_bloom, 
                                  int after_kmer, int protein_alignments, int mutations_found);
    void add_protein_match_to_report(void* reporter, const ProteinMatch* match);
    void generate_diagnostic_report(void* reporter);
}

// Enhanced ProteinMatch structure (ensure consistency)
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
    char query_peptide[51];  // Store aligned peptide sequence (up to 50 AA + null terminator)
    bool is_qrdr_alignment;  // Flag for QRDR region alignment
};

// Debug macros
#define DEBUG 0
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { fprintf(stderr, "[DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); fflush(stderr); }

// CUDA error checking
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "[CUDA ERROR] %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// FASTQ record structure
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Read batch for GPU processing
struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

// CSV writer for QRDR alignments
void writeQRDRAlignmentsToCSV(const std::vector<ProteinMatch>& qrdr_alignments,
                              const std::map<uint32_t, std::string>& gene_id_to_name,
                              const std::map<uint32_t, std::string>& species_id_to_name,
                              const std::string& csv_path,
                              int total_reads_processed) {
    std::ofstream csv_file(csv_path);
    if (!csv_file.good()) {
        std::cerr << "ERROR: Failed to create CSV file: " << csv_path << std::endl;
        return;
    }
    
    // Write CSV header
    csv_file << "read_id,gene_name,species_name,frame,alignment_score,identity_percent,"
             << "ref_start,ref_end,query_start,query_end,match_length,"
             << "num_mutations,mutation_positions,mutation_changes,"
             << "peptide_sequence,covers_qrdr_positions\n";
    
    // Helper to identify which QRDR positions are covered
    auto getQRDRPositions = [](uint32_t gene_id, uint16_t ref_start, uint16_t match_length) -> std::string {
        std::vector<int> covered_positions;
        uint16_t ref_end = ref_start + match_length;
        
        if (gene_id == 0) {  // gyrA
            if (ref_start <= 82 && ref_end > 82) covered_positions.push_back(83);  // 0-based to 1-based
            if (ref_start <= 86 && ref_end > 86) covered_positions.push_back(87);
        } else if (gene_id == 1) {  // parC
            if (ref_start <= 79 && ref_end > 79) covered_positions.push_back(80);
            if (ref_start <= 83 && ref_end > 83) covered_positions.push_back(84);
        } else if (gene_id == 2) {  // gyrB
            if (ref_start <= 425 && ref_end > 425) covered_positions.push_back(426);
            if (ref_start <= 446 && ref_end > 446) covered_positions.push_back(447);
        } else if (gene_id == 3) {  // parE
            if (ref_start <= 415 && ref_end > 415) covered_positions.push_back(416);
            if (ref_start <= 419 && ref_end > 419) covered_positions.push_back(420);
        }
        
        std::string result;
        for (size_t i = 0; i < covered_positions.size(); i++) {
            if (i > 0) result += ";";
            result += std::to_string(covered_positions[i]);
        }
        return result.empty() ? "none" : result;
    };
    
    // Write each QRDR alignment
    for (const auto& match : qrdr_alignments) {
        // Get gene and species names
        std::string gene_name = gene_id_to_name.count(match.gene_id) ? 
                               gene_id_to_name.at(match.gene_id) : "Gene" + std::to_string(match.gene_id);
        std::string species_name = species_id_to_name.count(match.species_id) ? 
                                  species_id_to_name.at(match.species_id) : "Species" + std::to_string(match.species_id);
        
        // Format mutation positions and changes
        std::string mutation_positions;
        std::string mutation_changes;
        for (int i = 0; i < match.num_mutations; i++) {
            if (i > 0) {
                mutation_positions += ";";
                mutation_changes += ";";
            }
            int global_pos = match.ref_start + match.mutation_positions[i] + 1;  // Convert to 1-based
            mutation_positions += std::to_string(global_pos);
            mutation_changes += std::string(1, match.ref_aas[i]) + std::to_string(global_pos) + 
                               std::string(1, match.query_aas[i]);
        }
        
        if (mutation_positions.empty()) {
            mutation_positions = "none";
            mutation_changes = "none";
        }
        
        // Get covered QRDR positions
        std::string qrdr_positions = getQRDRPositions(match.gene_id, match.ref_start, match.match_length);
        
        // Write row
        csv_file << match.read_id << ","
                 << gene_name << ","
                 << species_name << ","
                 << (int)match.frame << ","
                 << std::fixed << std::setprecision(1) << match.alignment_score << ","
                 << std::fixed << std::setprecision(1) << (match.identity * 100) << ","
                 << (match.ref_start + 1) << ","  // Convert to 1-based
                 << (match.ref_start + match.match_length) << ","  // End position (1-based)
                 << (match.query_start + 1) << ","  // Convert to 1-based
                 << (match.query_start + match.match_length) << ","
                 << match.match_length << ","
                 << (int)match.num_mutations << ","
                 << mutation_positions << ","
                 << mutation_changes << ","
                 << match.query_peptide << ","
                 << qrdr_positions << "\n";
    }
    
    csv_file.close();
    
    std::cout << "\nWrote " << qrdr_alignments.size() << " QRDR alignments to: " << csv_path << std::endl;
}

// Enhanced mutation detection logic with species-aware resistance patterns
class MutationAnalyzer {
private:
    // Species and gene mappings loaded from database metadata
    std::map<uint32_t, std::string> species_id_to_name;
    std::map<uint32_t, std::string> gene_id_to_name;
    
    // Species-aware resistance mutation definition
    struct SpeciesAwareMutation {
        uint32_t species_id;
        uint32_t gene_id;
        int position;
        char wildtype_aa;
        std::vector<char> resistant_aas;
        float resistance_level;
        std::string description;
    };
    
    // Species-specific resistance mutations based on published data
    std::vector<SpeciesAwareMutation> species_mutations = {
        // Enterococcus faecium (species_id=0) mutations
        {0, 0, 83, 'S', {'R', 'I', 'Y', 'N'}, 0.9f, "E.faecium gyrA S83R/I/Y/N"},
        {0, 0, 87, 'E', {'G', 'K'}, 0.8f, "E.faecium gyrA E87G/K"},
        {0, 1, 80, 'S', {'I', 'R'}, 0.7f, "E.faecium parC S80I/R"},
        {0, 1, 84, 'E', {'K', 'A'}, 0.6f, "E.faecium parC E84K/A"},
        
        // Escherichia coli (species_id=1) mutations
        {1, 0, 83, 'S', {'L', 'F', 'W'}, 0.9f, "E.coli gyrA S83L/F/W"},
        {1, 0, 87, 'D', {'N', 'G', 'Y', 'H'}, 0.7f, "E.coli gyrA D87N/G/Y/H"},
        {1, 1, 80, 'S', {'I', 'R', 'F'}, 0.6f, "E.coli parC S80I/R/F"},
        {1, 1, 84, 'E', {'V', 'K', 'G', 'A'}, 0.4f, "E.coli parC E84V/K/G/A"},
        {1, 2, 529, 'I', {'L'}, 0.5f, "E.coli parE I529L"},
    };
    
public:
    // Load species and gene mappings from database metadata
    void loadDatabaseMappings(const std::string& db_path) {
        std::string metadata_path = db_path + "/metadata.json";
        std::ifstream metadata_file(metadata_path);
        
        if (!metadata_file.good()) {
            std::cerr << "WARNING: Could not read metadata for species-aware detection: " << metadata_path << std::endl;
            return;
        }
        
        std::string line;
        bool in_species_map = false;
        bool in_gene_map = false;
        
        while (std::getline(metadata_file, line)) {
            // Parse species_map section
            if (line.find("\"species_map\"") != std::string::npos) {
                in_species_map = true;
                continue;
            }
            if (line.find("\"gene_map\"") != std::string::npos) {
                in_gene_map = true;
                in_species_map = false;
                continue;
            }
            
            // End of current section
            if ((in_species_map || in_gene_map) && line.find("}") != std::string::npos) {
                in_species_map = false;
                in_gene_map = false;
                continue;
            }
            
            // Parse mappings: "0": "Enterococcus_faecium",
            if ((in_species_map || in_gene_map) && line.find("\":") != std::string::npos) {
                size_t id_start = line.find("\"") + 1;
                size_t id_end = line.find("\"", id_start);
                size_t name_start = line.find("\"", id_end + 1) + 1;
                size_t name_end = line.find("\"", name_start);
                
                if (id_end != std::string::npos && name_end != std::string::npos) {
                    uint32_t id = std::stoi(line.substr(id_start, id_end - id_start));
                    std::string name = line.substr(name_start, name_end - name_start);
                    
                    if (in_species_map) {
                        species_id_to_name[id] = name;
                        std::cout << "[SPECIES-AWARE] Loaded species mapping: " << id << " -> " << name << std::endl;
                    } else if (in_gene_map) {
                        gene_id_to_name[id] = name;
                        std::cout << "[SPECIES-AWARE] Loaded gene mapping: " << id << " -> " << name << std::endl;
                    }
                }
            }
        }
        
        metadata_file.close();
        std::cout << "[SPECIES-AWARE] Loaded " << species_id_to_name.size() << " species and " 
                  << gene_id_to_name.size() << " gene mappings" << std::endl;
    }
    // Species-aware resistance mutation detection
    int detectResistanceMutations(ProteinMatch& match) {
        int resistance_mutations = 0;
        
        if (match.num_mutations == 0) {
            return 0; // No mutations detected by Smith-Waterman
        }
        
        // Get species and gene names for reporting
        std::string species_name = species_id_to_name.count(match.species_id) ? 
                                  species_id_to_name[match.species_id] : "Unknown";
        std::string gene_name = gene_id_to_name.count(match.gene_id) ? 
                               gene_id_to_name[match.gene_id] : "Unknown";
        
        // Debug: Print all mutations found by Smith-Waterman for first few matches
        static int debug_count = 0;
        // Species-aware debug output disabled for production
        /*
        if (debug_count < 5) {
            printf("[SPECIES-AWARE DEBUG] Read %d, %s %s (species_id=%d, gene_id=%d): %d mutations detected by SW\n", 
                   match.read_id, species_name.c_str(), gene_name.c_str(), 
                   match.species_id, match.gene_id, match.num_mutations);
            for (int i = 0; i < match.num_mutations; i++) {
                int global_pos = match.ref_start + match.mutation_positions[i];
                printf("  Mutation %d: position %d, %c->%c\n", 
                       i, global_pos, match.ref_aas[i], match.query_aas[i]);
            }
            
            // Show species-specific resistance patterns we're looking for
            if (debug_count == 0) {
                printf("[SPECIES-AWARE PATTERNS] Looking for these resistance mutations:\n");
                for (const auto& mut : species_mutations) {
                    std::string sp_name = species_id_to_name.count(mut.species_id) ? 
                                         species_id_to_name[mut.species_id] : "Unknown";
                    std::string gn_name = gene_id_to_name.count(mut.gene_id) ? 
                                         gene_id_to_name[mut.gene_id] : "Unknown";
                    printf("  %s %s: position %d, %c->", sp_name.c_str(), gn_name.c_str(), 
                           mut.position, mut.wildtype_aa);
                    for (char aa : mut.resistant_aas) printf("%c", aa);
                    printf(" (%s)\n", mut.description.c_str());
                }
            }
            debug_count++;
        }
        */
        
        // Species-aware mutation classification
        for (int i = 0; i < match.num_mutations; i++) {
            int global_pos = match.ref_start + match.mutation_positions[i];
            char observed_aa = match.query_aas[i];
            char reference_aa = match.ref_aas[i];
            
            // Look for species+gene specific resistance mutations
            for (const auto& resistance_mut : species_mutations) {
                // Must match species, gene, and position
                if (resistance_mut.species_id != match.species_id || 
                    resistance_mut.gene_id != match.gene_id || 
                    resistance_mut.position != global_pos) {
                    continue;
                }
                
                // Mutation detection logging disabled for production
                /*
                printf("[SPECIES-AWARE MATCH] Found mutation at %s %s resistance position %d: %c->%c\n",
                       species_name.c_str(), gene_name.c_str(), global_pos, reference_aa, observed_aa);
                */
                
                // Check if wildtype amino acid matches expected
                if (reference_aa != resistance_mut.wildtype_aa) {
                    // Warning logging disabled for production
                    /*
                    printf("[WARNING] Expected wildtype %c at position %d, but found %c in reference\n",
                           resistance_mut.wildtype_aa, global_pos, reference_aa);
                    */
                }
                
                // Check if observed AA is a known resistance variant for this species
                if (std::find(resistance_mut.resistant_aas.begin(), 
                             resistance_mut.resistant_aas.end(), 
                             observed_aa) != resistance_mut.resistant_aas.end()) {
                    resistance_mutations++;
                    // Resistance detection logging disabled for production
                    /*
                    printf("[🚨 RESISTANCE DETECTED] %s: %c%d%c (resistance level: %.1f)\n",
                           resistance_mut.description.c_str(), 
                           reference_aa, global_pos, observed_aa, resistance_mut.resistance_level);
                    */
                } else {
                    // Variant logging disabled for production
                    /*
                    printf("[VARIANT] %s %s position %d: %c->%c (not known resistance pattern: ",
                           species_name.c_str(), gene_name.c_str(), global_pos, reference_aa, observed_aa);
                    for (char aa : resistance_mut.resistant_aas) printf("%c", aa);
                    printf(")\n");
                    */
                }
                // Continue checking other resistance patterns - don't break!
                // This allows detection of cross-species resistance patterns
            }
        }
        
        return resistance_mutations;
    }
    
    // Check if alignment covers key resistance positions (species-aware)
    bool coversResistanceRegion(const ProteinMatch& match) {
        for (const auto& resistance_mut : species_mutations) {
            // Check all species patterns, not just exact species/gene match
            // This allows cross-species detection and ranking by alignment score
            if (resistance_mut.gene_id != match.gene_id) {
                continue; // Only check same gene, but allow different species
            }
            
            int rel_pos = resistance_mut.position - match.ref_start;
            if (rel_pos >= 0 && rel_pos < match.match_length) {
                return true; // Alignment covers a resistance position
            }
        }
        return false;
    }
};

class FastqReader {
private:
    gzFile file;
    bool is_open;
    
public:
    FastqReader(const std::string& filename) : is_open(false) {
        DEBUG_PRINT("Opening FASTQ file: %s", filename.c_str());
        file = gzopen(filename.c_str(), "r");
        is_open = (file != NULL);
        if (!is_open) {
            DEBUG_PRINT("ERROR: Failed to open FASTQ file: %s", filename.c_str());
        } else {
            DEBUG_PRINT("Successfully opened FASTQ file");
        }
    }
    
    ~FastqReader() {
        if (is_open) {
            gzclose(file);
        }
    }
    
    bool readRecord(FastqRecord& record) {
        if (!is_open) return false;
        
        const int buffer_size = 1024;
        char buffer[buffer_size];
        
        // Read header
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.header = std::string(buffer);
        if (record.header.empty() || record.header[0] != '@') return false;
        
        // Read sequence
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.sequence = std::string(buffer);
        record.sequence.pop_back(); // Remove newline
        
        // Read plus line
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        
        // Read quality
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.quality = std::string(buffer);
        record.quality.pop_back(); // Remove newline
        
        return true;
    }
    
    bool isOpen() const { return is_open; }
};

class EnhancedFQResistancePipeline {
private:
    // GPU memory for reads
    char* d_reads_r1;
    char* d_reads_r2;
    int* d_lengths_r1;
    int* d_lengths_r2;
    int* d_offsets_r1;
    int* d_offsets_r2;
    
    // Bloom filter memory
    void* bloom_filter;
    bool* d_bloom_passes_r1;
    bool* d_bloom_passes_r2;
    int* d_bloom_kmers_r1;
    int* d_bloom_kmers_r2;
    int* d_bloom_debug_stats;
    
    // Results for nucleotide search
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    AlignmentResult* d_results;
    uint32_t* d_result_count;
    
    // Enhanced translated search engine and results
    void* translated_search_engine;
    ProteinMatch* d_protein_matches;
    uint32_t* d_protein_match_counts;
    bool enable_translated_search;
    bool enable_smith_waterman;
    bool enable_diagnostic_reporting;
    bool enable_bloom_filter;  // NEW: Optional bloom filter
    bool enable_kmer_match;    // NEW: Optional k-mer matching
    std::string protein_db_path;
    
    // Diagnostic reporting
    void* diagnostic_reporter;
    MutationAnalyzer mutation_analyzer;
    
    // HDF5 writer
    HDF5AlignmentWriter* hdf5_writer;
    
    // Batch parameters
    const int batch_size = 10000;
    const int max_read_length = 300;
    const int bloom_min_kmers = 3;
    const int kmer_length = 15;
    const int max_protein_matches_per_read = 32;
    
    FQMutationDetectorCUDA detector;
    std::string current_index_path;
    
    // Enhanced statistics tracking
    struct EnhancedStats {
        int total_reads_processed = 0;
        int reads_passed_bloom = 0;
        int reads_with_candidates = 0;
        int reads_with_kmer_hits = 0;  // Reads retained after k-mer enrichment
        int total_protein_matches = 0;
        int reads_with_sw_alignments = 0;  // Reads with Smith-Waterman alignments
        int resistance_mutations_found = 0;
        int qrdr_alignments = 0;
        int high_confidence_matches = 0;
        int total_candidates_found = 0;
        int total_mutations_found = 0;
        int total_protein_resistance_found = 0;
        int total_smith_waterman_alignments = 0;
    } enhanced_stats;
    
    // Store top Smith-Waterman results for analysis
    std::vector<ProteinMatch> top_sw_matches;
    std::vector<ProteinMatch> qrdr_alignments;  // Store all QRDR alignments for ranking
    std::set<int> reads_with_sw_results;  // Track unique reads with SW results
    
    // Dynamic gene ID to name mapping loaded from protein database metadata
    std::map<uint32_t, std::string> gene_id_to_name;
    
    // Dynamic species ID to name mapping loaded from protein database metadata
    std::map<uint32_t, std::string> species_id_to_name;
    
    // QRDR position definitions for key genes
    std::map<uint32_t, std::pair<int, int>> qrdr_regions = {
        {0, {75, 95}},   // gyrA QRDR (approximately positions 75-95)
        {1, {75, 90}},   // parC QRDR (approximately positions 75-90)
        {2, {420, 450}}, // gyrB QRDR
        {3, {410, 425}}  // parE QRDR
    };
    
    // Helper function to check if alignment overlaps with QRDR
    bool isQRDRAlignment(const ProteinMatch& match) {
        auto it = qrdr_regions.find(match.gene_id);
        if (it == qrdr_regions.end()) return false;
        
        int qrdr_start = it->second.first;
        int qrdr_end = it->second.second;
        int match_start = match.ref_start;
        int match_end = match.ref_start + match.match_length;
        
        // Check for overlap
        return !(match_end <= qrdr_start || match_start >= qrdr_end);
    }
    
    // Helper function to extract peptide sequence from translated frame
    void extractPeptideSequence(ProteinMatch& match, const char* translated_sequence, int frame_length) {
        // Initialize peptide sequence
        memset(match.query_peptide, 0, sizeof(match.query_peptide));
        
        // Ensure we don't exceed bounds
        int start_pos = match.query_start;
        int length = std::min((int)match.match_length, 50);
        
        if (start_pos >= 0 && start_pos + length <= frame_length) {
            strncpy(match.query_peptide, translated_sequence + start_pos, length);
            match.query_peptide[length] = '\0';
        } else {
            strcpy(match.query_peptide, "SEQUENCE_UNAVAILABLE");
        }
    }
    
public:
    EnhancedFQResistancePipeline(bool use_translated_search = false, 
                                bool use_smith_waterman = false, 
                                bool use_diagnostic_reporting = false,
                                bool use_bloom_filter = true,      // NEW: Default on
                                bool use_kmer_match = true)         // NEW: Default on
        : enable_translated_search(use_translated_search), 
          enable_smith_waterman(use_smith_waterman), 
          enable_diagnostic_reporting(use_diagnostic_reporting),
          enable_bloom_filter(use_bloom_filter),
          enable_kmer_match(use_kmer_match) {
        
        DEBUG_PRINT("Initializing Enhanced FQ Resistance Pipeline");
        DEBUG_PRINT("Batch size: %d, Max read length: %d", batch_size, max_read_length);
        DEBUG_PRINT("Translated search: %s", enable_translated_search ? "ENABLED" : "DISABLED");
        DEBUG_PRINT("Smith-Waterman: %s", enable_smith_waterman ? "ENABLED" : "DISABLED");
        DEBUG_PRINT("Diagnostic reporting: %s", enable_diagnostic_reporting ? "ENABLED" : "DISABLED");
        DEBUG_PRINT("Bloom filter: %s", enable_bloom_filter ? "ENABLED" : "DISABLED");
        DEBUG_PRINT("K-mer matching: %s", enable_kmer_match ? "ENABLED" : "DISABLED");
        
        // Initialize diagnostic reporting
        diagnostic_reporter = nullptr;
        
        // Initialize HDF5 writer to null
        hdf5_writer = nullptr;
        
        // Create Bloom filter only if enabled
        bloom_filter = nullptr;
        if (enable_bloom_filter) {
            bloom_filter = create_bloom_filter(kmer_length);
            if (!bloom_filter) {
                std::cerr << "ERROR: Failed to create Bloom filter" << std::endl;
                exit(1);
            }
            DEBUG_PRINT("Bloom filter created successfully");
        }
        
        // Initialize enhanced translated search engine
        translated_search_engine = nullptr;
        d_protein_matches = nullptr;
        d_protein_match_counts = nullptr;
        
        if (enable_translated_search) {
            translated_search_engine = create_translated_search_engine_with_sw(batch_size, enable_smith_waterman);
            if (!translated_search_engine) {
                std::cerr << "ERROR: Failed to create enhanced translated search engine" << std::endl;
                exit(1);
            }
            DEBUG_PRINT("Enhanced translated search engine created (SW: %s)", 
                       enable_smith_waterman ? "enabled" : "disabled");
            
            // Allocate memory for protein search results
            CUDA_CHECK(cudaMalloc(&d_protein_matches, 
                                  batch_size * max_protein_matches_per_read * sizeof(ProteinMatch)));
            CUDA_CHECK(cudaMalloc(&d_protein_match_counts, batch_size * sizeof(uint32_t)));
        }
        
        // Allocate GPU memory for batch processing
        size_t read_buffer_size = batch_size * max_read_length;
        DEBUG_PRINT("Allocating GPU memory: read buffer size = %zu bytes", read_buffer_size);
        
        CUDA_CHECK(cudaMalloc(&d_reads_r1, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_reads_r2, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_lengths_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_lengths_r2, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r2, batch_size * sizeof(int)));
        
        // Allocate Bloom filter results only if enabled
        if (enable_bloom_filter) {
            CUDA_CHECK(cudaMalloc(&d_bloom_passes_r1, batch_size * sizeof(bool)));
            CUDA_CHECK(cudaMalloc(&d_bloom_passes_r2, batch_size * sizeof(bool)));
            CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r1, batch_size * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r2, batch_size * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&d_bloom_debug_stats, 4 * sizeof(int)));
        }
        
        // Allocate memory for nucleotide search results only if k-mer matching is enabled
        if (enable_kmer_match) {
            size_t candidates_size = batch_size * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch);
            size_t results_size = batch_size * MAX_CANDIDATES_PER_READ * sizeof(AlignmentResult);
            DEBUG_PRINT("Allocating GPU memory: candidates = %zu bytes, results = %zu bytes", 
                        candidates_size, results_size);
            
            CUDA_CHECK(cudaMalloc(&d_candidates, candidates_size));
            CUDA_CHECK(cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t)));
            CUDA_CHECK(cudaMalloc(&d_results, results_size));
            CUDA_CHECK(cudaMalloc(&d_result_count, sizeof(uint32_t)));
        }
        
        DEBUG_PRINT("GPU memory allocation completed successfully");
    }
    
    ~EnhancedFQResistancePipeline() {
        // Destroy diagnostic reporter
        if (diagnostic_reporter) {
            destroy_diagnostic_reporter(diagnostic_reporter);
        }
        
        // Destroy HDF5 writer
        if (hdf5_writer) {
            delete hdf5_writer;
        }
        
        // Destroy enhanced translated search engine
        if (translated_search_engine) {
            destroy_translated_search_engine(translated_search_engine);
        }
        
        // Destroy Bloom filter
        if (bloom_filter) {
            destroy_bloom_filter(bloom_filter);
        }
        
        // Free GPU memory
        cudaFree(d_reads_r1);
        cudaFree(d_reads_r2);
        cudaFree(d_lengths_r1);
        cudaFree(d_lengths_r2);
        cudaFree(d_offsets_r1);
        cudaFree(d_offsets_r2);
        
        if (enable_bloom_filter) {
            cudaFree(d_bloom_passes_r1);
            cudaFree(d_bloom_passes_r2);
            cudaFree(d_bloom_kmers_r1);
            cudaFree(d_bloom_kmers_r2);
            cudaFree(d_bloom_debug_stats);
        }
        
        if (enable_kmer_match) {
            cudaFree(d_candidates);
            cudaFree(d_candidate_counts);
            cudaFree(d_results);
            cudaFree(d_result_count);
        }
        
        if (d_protein_matches) cudaFree(d_protein_matches);
        if (d_protein_match_counts) cudaFree(d_protein_match_counts);
    }
    
    void initializeDiagnosticReporting(const std::string& output_path) {
        std::string diagnostic_path = output_path.substr(0, output_path.find_last_of('.')) + "_diagnostic.txt";
        diagnostic_reporter = create_diagnostic_reporter(diagnostic_path.c_str());
        
        if (diagnostic_reporter) {
            std::cout << "Diagnostic reporting enabled: " << diagnostic_path << std::endl;
        }
    }
    
    void setProteinDatabase(const std::string& db_path) {
        protein_db_path = db_path;
        if (translated_search_engine && !protein_db_path.empty()) {
            DEBUG_PRINT("Loading enhanced protein database from: %s", protein_db_path.c_str());
            
            // Load gene mappings from metadata
            loadGeneIdToNameMapping(protein_db_path);
            
            // Load species-aware resistance mutation mappings
            mutation_analyzer.loadDatabaseMappings(protein_db_path);
            
            int result = load_protein_database(translated_search_engine, protein_db_path.c_str());
            if (result != 0) {
                std::cerr << "WARNING: Failed to load enhanced protein database" << std::endl;
                enable_translated_search = false;
            } else {
                std::cout << "Enhanced protein database loaded successfully (5-mer k-mers)" << std::endl;
                if (enable_smith_waterman) {
                    std::cout << "Smith-Waterman alignment enabled for high-scoring matches" << std::endl;
                }
            }
        }
    }
    
    void loadGeneIdToNameMapping(const std::string& db_path) {
        std::string metadata_path = db_path + "/metadata.json";
        std::ifstream metadata_file(metadata_path);
        
        if (!metadata_file.good()) {
            std::cerr << "WARNING: Could not read metadata.json from " << metadata_path << std::endl;
            std::cerr << "Using default gene ID mappings" << std::endl;
            
            // Fallback mappings
            gene_id_to_name[0] = "gene0";
            gene_id_to_name[1] = "gene1";
            gene_id_to_name[2] = "gene2";
            gene_id_to_name[3] = "gene3";
            species_id_to_name[0] = "Species0";
            species_id_to_name[1] = "Species1";
            return;
        }
        
        std::string line;
        bool in_species_map = false;
        bool in_gene_map = false;
        
        while (std::getline(metadata_file, line)) {
            // Look for "species_map" section
            if (line.find("\"species_map\"") != std::string::npos) {
                in_species_map = true;
                in_gene_map = false;
                continue;
            }
            
            // Look for "gene_map" section
            if (line.find("\"gene_map\"") != std::string::npos) {
                in_gene_map = true;
                in_species_map = false;
                continue;
            }
            
            // End of current section
            if ((in_species_map || in_gene_map) && line.find("}") != std::string::npos) {
                in_species_map = false;
                in_gene_map = false;
                continue;
            }
            
            // Parse mappings: "0": "name",
            if ((in_species_map || in_gene_map) && line.find("\":") != std::string::npos) {
                size_t id_start = line.find("\"") + 1;
                size_t id_end = line.find("\"", id_start);
                size_t name_start = line.find("\"", id_end + 1) + 1;
                size_t name_end = line.find("\"", name_start);
                
                if (id_end != std::string::npos && name_end != std::string::npos) {
                    uint32_t id = std::stoi(line.substr(id_start, id_end - id_start));
                    std::string name = line.substr(name_start, name_end - name_start);
                    
                    if (in_species_map) {
                        species_id_to_name[id] = name;
                        std::cout << "Loaded species mapping: " << id << " -> " << name << std::endl;
                    } else if (in_gene_map) {
                        gene_id_to_name[id] = name;
                        std::cout << "Loaded gene mapping: " << id << " -> " << name << std::endl;
                    }
                }
            }
        }
        
        metadata_file.close();
        
        if (gene_id_to_name.empty()) {
            std::cerr << "WARNING: No gene mappings found in metadata.json" << std::endl;
            // Fallback
            gene_id_to_name[0] = "gene0";
            gene_id_to_name[1] = "gene1";
        }
        
        if (species_id_to_name.empty()) {
            std::cerr << "WARNING: No species mappings found in metadata.json" << std::endl;
            // Fallback
            species_id_to_name[0] = "Species0";
            species_id_to_name[1] = "Species1";
        }
    }
    
    void setSmithWatermanEnabled(bool enabled) {
        enable_smith_waterman = enabled;
        if (translated_search_engine) {
            set_smith_waterman_enabled(translated_search_engine, enabled);
            DEBUG_PRINT("Smith-Waterman alignment %s", enabled ? "ENABLED" : "DISABLED");
        }
    }
    
    void loadIndex(const std::string& index_path) {
        DEBUG_PRINT("Loading index from: %s", index_path.c_str());
        
        // Store index path for HDF5 initialization
        current_index_path = index_path;
        
        // Load k-mer index only if k-mer matching is enabled
        if (enable_kmer_match) {
            detector.loadIndex(index_path.c_str());
            DEBUG_PRINT("K-mer index loaded");
        }
        
        // Build or load Bloom filter only if enabled
        if (enable_bloom_filter && enable_kmer_match) {
            std::string bloom_path = index_path + "/bloom_filter.bin";
            
            // Try to load existing Bloom filter
            if (load_bloom_filter(bloom_filter, bloom_path.c_str()) == 0) {
                std::cout << "Loaded pre-built Bloom filter from: " << bloom_path << std::endl;
                
                // Check if it has RC k-mers
                bool has_rc = bloom_filter_has_rc(bloom_filter);
                std::cout << "Bloom filter contains RC k-mers: " << (has_rc ? "YES" : "NO") << std::endl;
                
                if (!has_rc) {
                    std::cout << "WARNING: Bloom filter does not contain RC k-mers. Rebuilding..." << std::endl;
                    // Force rebuild with RC
                    if (detector.d_kmer_sorted && detector.num_kmers > 0) {
                        if (build_bloom_filter_from_index(bloom_filter, detector.d_kmer_sorted, detector.num_kmers) == 0) {
                            std::cout << "Bloom filter rebuilt with RC k-mers" << std::endl;
                            save_bloom_filter(bloom_filter, bloom_path.c_str());
                        }
                    }
                }
            } else {
                // Build Bloom filter from k-mer index
                std::cout << "Building Bloom filter from k-mer index (with RC k-mers)..." << std::endl;
                
                if (detector.d_kmer_sorted && detector.num_kmers > 0) {
                    if (build_bloom_filter_from_index(bloom_filter, detector.d_kmer_sorted, detector.num_kmers) == 0) {
                        std::cout << "Bloom filter built successfully with " << detector.num_kmers << " k-mers (including RC)" << std::endl;
                        
                        // Save for future use
                        save_bloom_filter(bloom_filter, bloom_path.c_str());
                        std::cout << "Saved Bloom filter to: " << bloom_path << std::endl;
                    } else {
                        std::cerr << "WARNING: Failed to build Bloom filter from index" << std::endl;
                    }
                } else {
                    std::cerr << "WARNING: No k-mers available to build Bloom filter" << std::endl;
                }
            }
        }
        
        DEBUG_PRINT("Index loading completed");
    }
    
    ReadBatch prepareBatch(const std::vector<FastqRecord>& records) {
        ReadBatch batch;
        batch.num_reads = records.size();
        
        // Calculate total size needed
        batch.total_bases = 0;
        for (const auto& record : records) {
            batch.total_bases += record.sequence.length();
        }
        
        // Allocate host memory
        batch.sequences = new char[batch.total_bases];
        batch.lengths = new int[batch.num_reads];
        batch.offsets = new int[batch.num_reads];
        
        // Copy sequences
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
    
    void processPairedReads(const std::string& r1_path, const std::string& r2_path, 
                           const std::string& output_path) {
        DEBUG_PRINT("Starting enhanced paired read processing with diagnostic reporting");
        
        // Extract sample name from r1_path by splitting on "_R1.fastq.gz"
        std::string sample_name;
        size_t split_pos = r1_path.find("_R1.fastq.gz");
        if (split_pos != std::string::npos) {
            // Find the last directory separator to get just the filename
            size_t dir_pos = r1_path.find_last_of("/\\");
            size_t start_pos = (dir_pos != std::string::npos) ? dir_pos + 1 : 0;
            sample_name = r1_path.substr(start_pos, split_pos - start_pos);
        } else {
            // Fallback: use base filename without extension
            size_t dir_pos = r1_path.find_last_of("/\\");
            size_t start_pos = (dir_pos != std::string::npos) ? dir_pos + 1 : 0;
            size_t dot_pos = r1_path.find_last_of(".");
            size_t end_pos = (dot_pos != std::string::npos) ? dot_pos : r1_path.length();
            sample_name = r1_path.substr(start_pos, end_pos - start_pos);
        }
        
        // Create report directory
        std::string report_dir = sample_name + ".fq_report";
        std::string mkdir_cmd = "mkdir -p " + report_dir;
        system(mkdir_cmd.c_str());
        
        // Update output paths to use the report directory
        std::string new_output_path = report_dir + "/" + sample_name + ".output.json";
        
        // Initialize diagnostic reporting
        if (enable_diagnostic_reporting) {
            initializeDiagnosticReporting(new_output_path);
        }
        
        // Create HDF5 output file in report directory
        std::string hdf5_path = report_dir + "/" + sample_name + ".h5";
        hdf5_writer = new HDF5AlignmentWriter(hdf5_path);
        hdf5_writer->initialize(current_index_path, r1_path, r2_path);
        DEBUG_PRINT("HDF5 output initialized: %s", hdf5_path.c_str());
        
        FastqReader reader1(r1_path);
        FastqReader reader2(r2_path);
        
        if (!reader1.isOpen() || !reader2.isOpen()) {
            std::cerr << "ERROR: Failed to open input files" << std::endl;
            DEBUG_PRINT("Reader1 open: %s, Reader2 open: %s", 
                       reader1.isOpen() ? "YES" : "NO",
                       reader2.isOpen() ? "YES" : "NO");
            return;
        }
        
        DEBUG_PRINT("Both FASTQ files opened successfully");
        
        std::ofstream output(new_output_path);
        output << "{\n";
        output << "  \"sample\": \"" << sample_name << "\",\n";
        output << "  \"hdf5_output\": \"" << hdf5_path << "\",\n";
        output << "  \"pipeline_version\": \"0.4.0-enhanced-diagnostics\",\n";
        output << "  \"diagnostic_reporting\": " << (enable_diagnostic_reporting ? "true" : "false") << ",\n";
        output << "  \"bloom_filter_enabled\": " << (enable_bloom_filter ? "true" : "false") << ",\n";
        output << "  \"kmer_match_enabled\": " << (enable_kmer_match ? "true" : "false") << ",\n";
        output << "  \"translated_search_enabled\": " << (enable_translated_search ? "true" : "false") << ",\n";
        output << "  \"smith_waterman_enabled\": " << (enable_smith_waterman ? "true" : "false") << ",\n";
        output << "  \"protein_kmer_size\": 5,\n";
        if (enable_translated_search && !protein_db_path.empty()) {
            output << "  \"protein_database\": \"" << protein_db_path << "\",\n";
        }
        output << "  \"mutations\": [\n";
        
        bool first_mutation = true;
        auto pipeline_start = std::chrono::high_resolution_clock::now();
        
        while (true) {
            // Read batch of paired reads
            std::vector<FastqRecord> batch_r1, batch_r2;
            
            DEBUG_PRINT("Reading batch %d (up to %d reads)", enhanced_stats.total_reads_processed/batch_size + 1, batch_size);
            
            for (int i = 0; i < batch_size; i++) {
                FastqRecord rec1, rec2;
                bool got1 = reader1.readRecord(rec1);
                bool got2 = reader2.readRecord(rec2);
                
                if (!got1 || !got2) {
                    if (i == 0) {
                        DEBUG_PRINT("No more reads available");
                    } else {
                        DEBUG_PRINT("Batch partially filled: %d reads", i);
                    }
                    break;
                }
                
                batch_r1.push_back(rec1);
                batch_r2.push_back(rec2);
            }
            
            if (batch_r1.empty()) {
                DEBUG_PRINT("Empty batch, ending processing");
                break;
            }
            
            DEBUG_PRINT("Batch contains %zu read pairs", batch_r1.size());
            
            enhanced_stats.total_reads_processed += batch_r1.size();
            
            // Prepare batches
            ReadBatch batch1 = prepareBatch(batch_r1);
            ReadBatch batch2 = prepareBatch(batch_r2);
            
            // Transfer to GPU
            DEBUG_PRINT("Transferring batch to GPU: %d reads, %d total bases (R1)", 
                       batch1.num_reads, batch1.total_bases);
            
            CUDA_CHECK(cudaMemcpy(d_reads_r1, batch1.sequences, batch1.total_bases, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lengths_r1, batch1.lengths, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_offsets_r1, batch1.offsets, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            
            CUDA_CHECK(cudaMemcpy(d_reads_r2, batch2.sequences, batch2.total_bases, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lengths_r2, batch2.lengths, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_offsets_r2, batch2.offsets, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            
            // Reset result counters
            if (enable_kmer_match) {
                CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
                CUDA_CHECK(cudaMemset(d_candidate_counts, 0, batch1.num_reads * sizeof(uint32_t)));
            }
            if (enable_translated_search) {
                CUDA_CHECK(cudaMemset(d_protein_match_counts, 0, batch1.num_reads * sizeof(uint32_t)));
            }
            
            DEBUG_PRINT("GPU transfer completed");
            
            // Initialize read pair pass status - default to all pass if no filtering
            std::vector<bool> h_pair_passes(batch1.num_reads, true);
            int reads_passed_bloom_combined = batch1.num_reads;
            
            // =========================
            // Stage 0: Bloom Filter Pre-screening (OPTIONAL)
            // =========================
            if (enable_bloom_filter) {
                auto bloom_start = std::chrono::high_resolution_clock::now();
                DEBUG_PRINT("Stage 0: Bloom filter pre-screening for %d reads", batch1.num_reads);
                
                // Reset debug statistics
                CUDA_CHECK(cudaMemset(d_bloom_debug_stats, 0, 4 * sizeof(int)));
                
                // Screen R1 reads (forward orientation expected)
                int bloom_result = bloom_filter_screen_reads_with_rc(
                    bloom_filter,
                    d_reads_r1, d_lengths_r1, d_offsets_r1,
                    batch1.num_reads,
                    d_bloom_passes_r1, d_bloom_kmers_r1,
                    bloom_min_kmers,
                    false,  // R1: don't check RC
                    d_bloom_debug_stats
                );
                
                if (bloom_result != 0) {
                    DEBUG_PRINT("WARNING: Bloom filter screening failed with code %d", bloom_result);
                }
                
                // Get R1 Bloom filter results
                std::vector<char> h_bloom_passes_raw(batch1.num_reads);
                std::vector<int> h_bloom_kmers(batch1.num_reads);
                CUDA_CHECK(cudaMemcpy(h_bloom_passes_raw.data(), d_bloom_passes_r1, 
                                     batch1.num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
                std::vector<bool> h_bloom_passes(batch1.num_reads);
                for (int i = 0; i < batch1.num_reads; i++) {
                    h_bloom_passes[i] = h_bloom_passes_raw[i];
                }
                CUDA_CHECK(cudaMemcpy(h_bloom_kmers.data(), d_bloom_kmers_r1, 
                                     batch1.num_reads * sizeof(int), cudaMemcpyDeviceToHost));
                
                // Get R1 statistics
                int r1_stats[4];
                CUDA_CHECK(cudaMemcpy(r1_stats, d_bloom_debug_stats, 4 * sizeof(int), cudaMemcpyDeviceToHost));
                
                // Count R1 reads that passed
                int reads_passed_bloom_r1 = 0;
                for (int i = 0; i < batch1.num_reads; i++) {
                    if (h_bloom_passes[i]) reads_passed_bloom_r1++;
                }
                
                DEBUG_PRINT("R1 Bloom filter results: %d/%d reads passed (%.1f%%)", 
                           reads_passed_bloom_r1, batch1.num_reads, 
                           100.0 * reads_passed_bloom_r1 / batch1.num_reads);
                DEBUG_PRINT("  R1 stats - Forward hits: %d, RC hits: %d", r1_stats[1], r1_stats[2]);
                
                // Screen R2 reads with RC support
                DEBUG_PRINT("Stage 0b: Bloom filter pre-screening for R2 reads (with RC)");
                
                // Reset debug statistics
                CUDA_CHECK(cudaMemset(d_bloom_debug_stats, 0, 4 * sizeof(int)));
                
                bloom_result = bloom_filter_screen_reads_with_rc(
                    bloom_filter,
                    d_reads_r2, d_lengths_r2, d_offsets_r2,
                    batch2.num_reads,
                    d_bloom_passes_r2, d_bloom_kmers_r2,
                    bloom_min_kmers,
                    true,   // R2: CHECK RC!
                    d_bloom_debug_stats
                );
                
                // Get R2 Bloom results
                std::vector<char> h_bloom_passes_r2_raw(batch2.num_reads);
                CUDA_CHECK(cudaMemcpy(h_bloom_passes_r2_raw.data(), d_bloom_passes_r2, 
                                     batch2.num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
                std::vector<bool> h_bloom_passes_r2(batch2.num_reads);
                for (int i = 0; i < batch2.num_reads; i++) {
                    h_bloom_passes_r2[i] = h_bloom_passes_r2_raw[i];
                }
                
                // Get R2 statistics
                int r2_stats[4];
                CUDA_CHECK(cudaMemcpy(r2_stats, d_bloom_debug_stats, 4 * sizeof(int), cudaMemcpyDeviceToHost));
                
                int reads_passed_bloom_r2 = 0;
                for (int i = 0; i < batch2.num_reads; i++) {
                    if (h_bloom_passes_r2[i]) reads_passed_bloom_r2++;
                }
                
                DEBUG_PRINT("R2 Bloom filter results: %d/%d reads passed (%.1f%%)", 
                           reads_passed_bloom_r2, batch2.num_reads, 
                           100.0 * reads_passed_bloom_r2 / batch2.num_reads);
                DEBUG_PRINT("  R2 stats - Forward hits: %d, RC hits: %d", r2_stats[1], r2_stats[2]);
                
                auto bloom_end = std::chrono::high_resolution_clock::now();
                auto bloom_time = std::chrono::duration_cast<std::chrono::microseconds>(bloom_end - bloom_start).count();
                
                // Combine R1 and R2 Bloom results - a read pair passes if EITHER read passes
                reads_passed_bloom_combined = 0;
                for (int i = 0; i < batch1.num_reads; i++) {
                    h_pair_passes[i] = h_bloom_passes[i] || h_bloom_passes_r2[i];
                    if (h_pair_passes[i]) reads_passed_bloom_combined++;
                }
                
                DEBUG_PRINT("Combined Bloom results: %d/%d read pairs passed (%.1f ms)", 
                           reads_passed_bloom_combined, batch1.num_reads, bloom_time / 1000.0);
                
                enhanced_stats.reads_passed_bloom += reads_passed_bloom_combined;
            } else {
                // No bloom filter - all reads pass
                DEBUG_PRINT("Bloom filter disabled - all %d reads proceed to next stage", batch1.num_reads);
            }
            
            // Only proceed if some read pairs passed Bloom filter (or bloom is disabled)
            if (reads_passed_bloom_combined > 0) {
                // We'll process both R1 and R2, then combine results
                std::map<int, AlignmentResult> best_results_per_read_pair;
                
                // =========================
                // Process R1/R2 reads with nucleotide search (OPTIONAL)
                // =========================
                if (enable_kmer_match) {
                    DEBUG_PRINT("Processing R1 reads with k-mer matching...");
                    
                    // Stage 1: K-mer Filtering for R1
                    auto kmer_start = std::chrono::high_resolution_clock::now();
                    DEBUG_PRINT("Stage 1: K-mer filtering for R1 reads");
                    
                    launch_kmer_filter(
                        d_reads_r1, d_lengths_r1, d_offsets_r1,
                        detector, batch1.num_reads,
                        d_candidates, d_candidate_counts,
                        false  // R1: don't check RC in k-mer screening
                    );
                    CUDA_CHECK(cudaDeviceSynchronize());
                    
                    // Check candidate counts for R1
                    std::vector<uint32_t> h_candidate_counts_r1(batch1.num_reads);
                    CUDA_CHECK(cudaMemcpy(h_candidate_counts_r1.data(), d_candidate_counts, 
                                         batch1.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                    
                    int total_candidates_r1 = 0;
                    int reads_with_kmer_hits_r1 = 0;
                    for (int i = 0; i < batch1.num_reads; i++) {
                        if (h_pair_passes[i] && h_candidate_counts_r1[i] > 0) {
                            total_candidates_r1 += h_candidate_counts_r1[i];
                            reads_with_kmer_hits_r1++;
                        }
                    }
                    
                    DEBUG_PRINT("R1 k-mer filtering: %d total candidates", total_candidates_r1);
                    
                    // Stage 2: Alignment for R1
                    if (total_candidates_r1 > 0) {
                        DEBUG_PRINT("Stage 2: Position-weighted alignment for R1");
                        
                        CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
                        
                        launch_position_weighted_alignment(
                            d_reads_r1, d_lengths_r1, d_offsets_r1,
                            d_candidates, d_candidate_counts,
                            detector, batch1.num_reads,
                            d_results, d_result_count
                        );
                        CUDA_CHECK(cudaDeviceSynchronize());
                        
                        // Get R1 results
                        uint32_t num_results_r1;
                        CUDA_CHECK(cudaMemcpy(&num_results_r1, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        
                        if (num_results_r1 > 0) {
                            std::vector<AlignmentResult> results_r1(num_results_r1);
                            cudaMemcpy(results_r1.data(), d_results, num_results_r1 * sizeof(AlignmentResult), 
                                      cudaMemcpyDeviceToHost);
                            
                            // Add to HDF5
                            hdf5_writer->addAlignmentBatch(results_r1.data(), num_results_r1, 
                                                         enhanced_stats.total_reads_processed - batch1.num_reads);
                            
                            // Store R1 results for JSON output
                            for (const auto& result : results_r1) {
                                if (result.num_mutations_detected > 0 && h_pair_passes[result.read_id]) {
                                    int read_pair_id = result.read_id;
                                    
                                    if (best_results_per_read_pair.find(read_pair_id) == best_results_per_read_pair.end() ||
                                        result.alignment_score > best_results_per_read_pair[read_pair_id].alignment_score) {
                                        best_results_per_read_pair[read_pair_id] = result;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Process R2 reads with nucleotide search (with RC support)
                    DEBUG_PRINT("Processing R2 reads (with RC support)...");
                    
                    // Reset candidate counts for R2
                    CUDA_CHECK(cudaMemset(d_candidate_counts, 0, batch2.num_reads * sizeof(uint32_t)));
                    
                    // Stage 1: K-mer Filtering for R2 with RC
                    DEBUG_PRINT("Stage 1: K-mer filtering for R2 reads (with RC)");
                    
                    launch_kmer_filter(
                        d_reads_r2, d_lengths_r2, d_offsets_r2,
                        detector, batch2.num_reads,
                        d_candidates, d_candidate_counts,
                        true   // R2: CHECK RC in k-mer screening!
                    );
                    CUDA_CHECK(cudaDeviceSynchronize());
                    
                    // Check candidate counts for R2
                    std::vector<uint32_t> h_candidate_counts_r2(batch2.num_reads);
                    CUDA_CHECK(cudaMemcpy(h_candidate_counts_r2.data(), d_candidate_counts, 
                                         batch2.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                    
                    int total_candidates_r2 = 0;
                    int reads_with_kmer_hits_r2 = 0;
                    for (int i = 0; i < batch2.num_reads; i++) {
                        if (h_pair_passes[i] && h_candidate_counts_r2[i] > 0) {
                            total_candidates_r2 += h_candidate_counts_r2[i];
                            reads_with_kmer_hits_r2++;
                        }
                    }
                    
                    DEBUG_PRINT("R2 k-mer filtering: %d total candidates", total_candidates_r2);
                    
                    // Stage 2: Alignment for R2
                    if (total_candidates_r2 > 0) {
                        DEBUG_PRINT("Stage 2: Position-weighted alignment for R2");
                        
                        CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
                        
                        launch_position_weighted_alignment(
                            d_reads_r2, d_lengths_r2, d_offsets_r2,
                            d_candidates, d_candidate_counts,
                            detector, batch2.num_reads,
                            d_results, d_result_count
                        );
                        CUDA_CHECK(cudaDeviceSynchronize());
                        
                        // Get R2 results
                        uint32_t num_results_r2;
                        CUDA_CHECK(cudaMemcpy(&num_results_r2, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        
                        if (num_results_r2 > 0) {
                            std::vector<AlignmentResult> results_r2(num_results_r2);
                            cudaMemcpy(results_r2.data(), d_results, num_results_r2 * sizeof(AlignmentResult), 
                                      cudaMemcpyDeviceToHost);
                            
                            // Add to HDF5
                            hdf5_writer->addAlignmentBatch(results_r2.data(), num_results_r2, 
                                                         enhanced_stats.total_reads_processed - batch2.num_reads);
                            
                            // Combine with R1 results for JSON output, preferring higher scores
                            for (const auto& result : results_r2) {
                                if (result.num_mutations_detected > 0 && h_pair_passes[result.read_id]) {
                                    int read_pair_id = result.read_id;
                                    
                                    // If R2 has better score than R1, use R2
                                    if (best_results_per_read_pair.find(read_pair_id) == best_results_per_read_pair.end() ||
                                        result.alignment_score > best_results_per_read_pair[read_pair_id].alignment_score) {
                                        best_results_per_read_pair[read_pair_id] = result;
                                    }
                                }
                            }
                        }
                    }
                    
                    auto kmer_end = std::chrono::high_resolution_clock::now();
                    auto kmer_time = std::chrono::duration_cast<std::chrono::microseconds>(kmer_end - kmer_start).count();
                    
                    DEBUG_PRINT("Combined nucleotide results: %zu read pairs with mutations (%.1f ms)", 
                               best_results_per_read_pair.size(), kmer_time / 1000.0);
                    
                    enhanced_stats.total_candidates_found += total_candidates_r1 + total_candidates_r2;
                    enhanced_stats.reads_with_candidates += reads_passed_bloom_combined;
                    
                    // Count unique read pairs with k-mer hits (either R1 or R2 has hits)
                    int read_pairs_with_kmer_hits = 0;
                    for (int i = 0; i < batch1.num_reads; i++) {
                        if ((h_pair_passes[i] && h_candidate_counts_r1[i] > 0) || 
                            (h_pair_passes[i] && h_candidate_counts_r2[i] > 0)) {
                            read_pairs_with_kmer_hits++;
                        }
                    }
                    enhanced_stats.reads_with_kmer_hits += read_pairs_with_kmer_hits;
                } else {
                    DEBUG_PRINT("K-mer matching disabled - skipping nucleotide search");
                }
                
                // =========================
                // Stage 3: Enhanced 5-mer Translated Search (if enabled)
                // =========================
                if (enable_translated_search && translated_search_engine) {
                    auto trans_start = std::chrono::high_resolution_clock::now();
                    DEBUG_PRINT("Stage 3: Enhanced 6-frame translated search (5-mer k-mers, SW: %s)", 
                               enable_smith_waterman ? "enabled" : "disabled");
                    
                    // Create bool array for reads to process based on h_pair_passes
                    bool* d_reads_to_process;
                    CUDA_CHECK(cudaMalloc(&d_reads_to_process, batch1.num_reads * sizeof(bool)));
                    
                    if (enable_bloom_filter) {
                        // Copy the bloom filter results
                        CUDA_CHECK(cudaMemcpy(d_reads_to_process, d_bloom_passes_r1, 
                                             batch1.num_reads * sizeof(bool), cudaMemcpyDeviceToDevice));
                    } else {
                        // All reads pass if bloom filter is disabled
                        std::vector<uint8_t> all_pass(batch1.num_reads, 1);
                        CUDA_CHECK(cudaMemcpy(d_reads_to_process, all_pass.data(), 
                                             batch1.num_reads * sizeof(bool), cudaMemcpyHostToDevice));
                    }
                    
                    // Search both R1 and R2 for protein matches
                    int trans_result = search_translated_reads(
                        translated_search_engine,
                        d_reads_r1, d_lengths_r1, d_offsets_r1,
                        d_reads_to_process,  // Only search reads that passed bloom
                        batch1.num_reads,
                        d_protein_matches,
                        d_protein_match_counts
                    );
                    
                    if (trans_result == 0) {
                        // Get protein match counts
                        std::vector<uint32_t> h_protein_match_counts(batch1.num_reads);
                        CUDA_CHECK(cudaMemcpy(h_protein_match_counts.data(), d_protein_match_counts,
                                             batch1.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        
                        // Process protein matches with enhanced mutation detection
                        if (std::accumulate(h_protein_match_counts.begin(), h_protein_match_counts.end(), 0) > 0) {
                            std::vector<ProteinMatch> h_protein_matches(batch1.num_reads * max_protein_matches_per_read);
                            CUDA_CHECK(cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                                                 batch1.num_reads * max_protein_matches_per_read * sizeof(ProteinMatch),
                                                 cudaMemcpyDeviceToHost));
                            
                            // Analyze each protein match
                            for (int i = 0; i < batch1.num_reads; i++) {
                                if (h_pair_passes[i] && h_protein_match_counts[i] > 0) {
                                    for (uint32_t j = 0; j < h_protein_match_counts[i]; j++) {
                                        ProteinMatch& pm = h_protein_matches[i * max_protein_matches_per_read + j];
                                        
                                        if (pm.used_smith_waterman) {
                                            enhanced_stats.total_protein_matches++;
                                        }
                                        
                                        // Enhanced mutation analysis
                                        int resistance_muts = mutation_analyzer.detectResistanceMutations(pm);
                                        enhanced_stats.resistance_mutations_found += resistance_muts;
                                        
                                        // Check if covers resistance region (Smith-Waterman only)
                                        if (pm.used_smith_waterman && mutation_analyzer.coversResistanceRegion(pm)) {
                                            enhanced_stats.qrdr_alignments++;
                                        }
                                        
                                        // Track high confidence matches (Smith-Waterman only)
                                        if (pm.used_smith_waterman && pm.identity >= 0.95f && pm.alignment_score >= 50.0f) {
                                            enhanced_stats.high_confidence_matches++;
                                        }
                                        
                                        if (pm.used_smith_waterman) {
                                            enhanced_stats.total_smith_waterman_alignments++;
                                            
                                            // Check if this is a QRDR alignment
                                            pm.is_qrdr_alignment = isQRDRAlignment(pm);
                                            
                                            // Store for top display
                                            top_sw_matches.push_back(pm);
                                            
                                            // Store QRDR alignments separately for ranking
                                            if (pm.is_qrdr_alignment) {
                                                qrdr_alignments.push_back(pm);
                                            }
                                            
                                            // Track unique reads with SW results (R1)
                                            reads_with_sw_results.insert(pm.read_id + enhanced_stats.total_reads_processed - batch1.num_reads);
                                        }
                                        
                                        // Add to diagnostic report (all matches for analysis)
                                        if (enable_diagnostic_reporting && diagnostic_reporter) {
                                            add_protein_match_to_report(diagnostic_reporter, &pm);
                                        }
                                        
                                        // Check if any mutations are at resistance positions
                                        bool has_resistance = false;
                                        for (int m = 0; m < pm.num_mutations; m++) {
                                            if (pm.gene_id == 0 && (pm.mutation_positions[m] == 83 || pm.mutation_positions[m] == 87)) {
                                                has_resistance = true;
                                                break;
                                            }
                                            if (pm.gene_id == 1 && (pm.mutation_positions[m] == 80 || pm.mutation_positions[m] == 84)) {
                                                has_resistance = true;
                                                break;
                                            }
                                        }
                                        if (has_resistance) {
                                            enhanced_stats.total_protein_resistance_found++;
                                        }
                                    }
                                }
                            }
                            
                            // Add to HDF5
                            hdf5_writer->addTranslatedResults(h_protein_matches.data(), h_protein_match_counts.data(),
                                                            batch1.num_reads, enhanced_stats.total_reads_processed - batch1.num_reads);
                        }
                    }
                    
                    // Also search R2 (optional, could skip for speed)
                    CUDA_CHECK(cudaMemset(d_protein_match_counts, 0, batch2.num_reads * sizeof(uint32_t)));
                    
                    if (enable_bloom_filter) {
                        CUDA_CHECK(cudaMemcpy(d_reads_to_process, d_bloom_passes_r2, 
                                             batch2.num_reads * sizeof(bool), cudaMemcpyDeviceToDevice));
                    }
                    
                    trans_result = search_translated_reads(
                        translated_search_engine,
                        d_reads_r2, d_lengths_r2, d_offsets_r2,
                        d_reads_to_process,
                        batch2.num_reads,
                        d_protein_matches,
                        d_protein_match_counts
                    );
                    
                    if (trans_result == 0) {
                        std::vector<uint32_t> h_protein_match_counts(batch2.num_reads);
                        CUDA_CHECK(cudaMemcpy(h_protein_match_counts.data(), d_protein_match_counts,
                                             batch2.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        
                        if (std::accumulate(h_protein_match_counts.begin(), h_protein_match_counts.end(), 0) > 0) {
                            std::vector<ProteinMatch> h_protein_matches(batch2.num_reads * max_protein_matches_per_read);
                            CUDA_CHECK(cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                                                 batch2.num_reads * max_protein_matches_per_read * sizeof(ProteinMatch),
                                                 cudaMemcpyDeviceToHost));
                            
                            // Process R2 protein matches
                            for (int i = 0; i < batch2.num_reads; i++) {
                                if (h_pair_passes[i] && h_protein_match_counts[i] > 0) {
                                    for (uint32_t j = 0; j < h_protein_match_counts[i]; j++) {
                                        ProteinMatch& pm = h_protein_matches[i * max_protein_matches_per_read + j];
                                        
                                        if (pm.used_smith_waterman) {
                                            enhanced_stats.total_protein_matches++;
                                        }
                                        
                                        // Enhanced mutation analysis
                                        int resistance_muts = mutation_analyzer.detectResistanceMutations(pm);
                                        enhanced_stats.resistance_mutations_found += resistance_muts;
                                        
                                        // Check if covers resistance region (Smith-Waterman only)
                                        if (pm.used_smith_waterman && mutation_analyzer.coversResistanceRegion(pm)) {
                                            enhanced_stats.qrdr_alignments++;
                                        }
                                        
                                        // Track high confidence matches (Smith-Waterman only)
                                        if (pm.used_smith_waterman && pm.identity >= 0.95f && pm.alignment_score >= 50.0f) {
                                            enhanced_stats.high_confidence_matches++;
                                        }
                                        
                                        if (pm.used_smith_waterman) {
                                            enhanced_stats.total_smith_waterman_alignments++;
                                            
                                            // Check if this is a QRDR alignment
                                            pm.is_qrdr_alignment = isQRDRAlignment(pm);
                                            
                                            // Store for top display
                                            top_sw_matches.push_back(pm);
                                            
                                            // Store QRDR alignments separately for ranking
                                            if (pm.is_qrdr_alignment) {
                                                qrdr_alignments.push_back(pm);
                                            }
                                            
                                            // Track unique reads with SW results (R2)
                                            reads_with_sw_results.insert(pm.read_id + enhanced_stats.total_reads_processed - batch2.num_reads);
                                        }
                                        
                                        // Add to diagnostic report
                                        if (enable_diagnostic_reporting && diagnostic_reporter) {
                                            add_protein_match_to_report(diagnostic_reporter, &pm);
                                        }
                                        
                                        // Check for resistance
                                        bool has_resistance = false;
                                        for (int m = 0; m < pm.num_mutations; m++) {
                                            // Convert relative position to global position (1-based)
                                            int global_pos = pm.ref_start + pm.mutation_positions[m] + 1;
                                            
                                            // Based on the wildtype_protein_db mapping:
                                            // gene_id 0 = parE (positions 416, 420)
                                            // gene_id 1 = gyrA (positions 83, 87)
                                            if (pm.gene_id == 1 && (global_pos == 83 || global_pos == 87)) {
                                                has_resistance = true;
                                                break;
                                            }
                                            if (pm.gene_id == 0 && (global_pos == 416 || global_pos == 420)) {
                                                has_resistance = true;
                                                break;
                                            }
                                        }
                                        if (has_resistance) {
                                            enhanced_stats.total_protein_resistance_found++;
                                        }
                                    }
                                }
                            }
                            
                            hdf5_writer->addTranslatedResults(h_protein_matches.data(), h_protein_match_counts.data(),
                                                            batch2.num_reads, enhanced_stats.total_reads_processed - batch2.num_reads);
                        }
                    }
                    
                    CUDA_CHECK(cudaFree(d_reads_to_process));
                    
                    auto trans_end = std::chrono::high_resolution_clock::now();
                    auto trans_time = std::chrono::duration_cast<std::chrono::microseconds>(trans_end - trans_start).count();
                    DEBUG_PRINT("Enhanced translated search completed (%.1f ms)", trans_time / 1000.0);
                }
                
                // Output combined results to JSON
                if (enable_kmer_match) {
                    for (const auto& pair : best_results_per_read_pair) {
                        const auto& result = pair.second;
                        if (!first_mutation) output << ",\n";
                        output << "    {\n";
                        output << "      \"read_pair\": " << (result.read_id + enhanced_stats.total_reads_processed - batch1.num_reads) << ",\n";
                        output << "      \"gene_id\": " << result.gene_id << ",\n";
                        output << "      \"species_id\": " << result.species_id << ",\n";
                        output << "      \"alignment_score\": " << result.alignment_score << ",\n";
                        output << "      \"identity\": " << result.identity << ",\n";
                        output << "      \"mutations_detected\": " << (int)result.num_mutations_detected << ",\n";
                        output << "      \"source\": \"" << (result.start_pos < 1000 ? "R1" : "R2") << "\"\n";
                        output << "    }";
                        first_mutation = false;
                        enhanced_stats.total_mutations_found++;
                    }
                }
            }
            
            // Update diagnostic statistics periodically
            if (enable_diagnostic_reporting && diagnostic_reporter && enhanced_stats.total_reads_processed % 50000 == 0) {
                update_pipeline_statistics(
                    diagnostic_reporter,
                    enhanced_stats.total_reads_processed,
                    enhanced_stats.reads_passed_bloom,
                    enhanced_stats.reads_with_candidates,
                    enhanced_stats.total_protein_matches,
                    enhanced_stats.resistance_mutations_found
                );
            }
            
            // Cleanup batch memory
            delete[] batch1.sequences;
            delete[] batch1.lengths;
            delete[] batch1.offsets;
            delete[] batch2.sequences;
            delete[] batch2.lengths;
            delete[] batch2.offsets;
            
            // Progress update disabled for production
            /*
            if (enhanced_stats.total_reads_processed % 100000 == 0) {
                std::cout << "Processed " << enhanced_stats.total_reads_processed << " read pairs..." << std::endl;
                if (enable_bloom_filter) {
                    std::cout << "  Bloom filter pass rate: " << (100.0 * enhanced_stats.reads_passed_bloom / enhanced_stats.total_reads_processed) << "%" << std::endl;
                }
                if (enable_translated_search) {
                    std::cout << "  Protein matches found: " << enhanced_stats.total_protein_matches << std::endl;
                    if (enable_smith_waterman) {
                        std::cout << "  Smith-Waterman alignments: " << enhanced_stats.total_smith_waterman_alignments << std::endl;
                    }
                }
            }
            */
        }
        
        auto pipeline_end = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(pipeline_end - pipeline_start).count();
        
        // Update final Smith-Waterman read count
        enhanced_stats.reads_with_sw_alignments = reads_with_sw_results.size();
        
        // Sort top Smith-Waterman matches by alignment score
        if (!top_sw_matches.empty()) {
            std::sort(top_sw_matches.begin(), top_sw_matches.end(),
                     [](const ProteinMatch& a, const ProteinMatch& b) {
                         return a.alignment_score > b.alignment_score;
                     });
        }
        
        // Create QRDR alignments CSV file
        if (!qrdr_alignments.empty()) {
            std::string csv_path = report_dir + "/" + sample_name + "_qrdr_alignments.csv";
            writeQRDRAlignmentsToCSV(qrdr_alignments, gene_id_to_name, species_id_to_name, 
                                    csv_path, enhanced_stats.total_reads_processed);
        }
        
        // Finalize HDF5 output
        hdf5_writer->finalize(output_path);
        
        // Finalize output
        output << "\n  ],\n";
        output << "  \"enhanced_summary\": {\n";
        output << "    \"total_reads\": " << enhanced_stats.total_reads_processed << ",\n";
        if (enable_bloom_filter) {
            output << "    \"reads_passed_bloom\": " << enhanced_stats.reads_passed_bloom << ",\n";
        }
        if (enable_kmer_match) {
            output << "    \"reads_with_candidates\": " << enhanced_stats.reads_with_candidates << ",\n";
            output << "    \"reads_with_kmer_hits\": " << enhanced_stats.reads_with_kmer_hits << ",\n";
        }
        if (enable_translated_search) {
            output << "    \"protein_matches\": " << enhanced_stats.total_protein_matches << ",\n";
            output << "    \"reads_with_sw_alignments\": " << enhanced_stats.reads_with_sw_alignments << ",\n";
            output << "    \"resistance_mutations_found\": " << enhanced_stats.resistance_mutations_found << ",\n";
            output << "    \"qrdr_alignments\": " << enhanced_stats.qrdr_alignments << ",\n";
            output << "    \"high_confidence_matches\": " << enhanced_stats.high_confidence_matches << ",\n";
        }
        if (enable_bloom_filter) {
            output << "    \"bloom_retention_rate\": " << (double)enhanced_stats.reads_passed_bloom / enhanced_stats.total_reads_processed << ",\n";
        }
        if (enable_translated_search && enhanced_stats.reads_passed_bloom > 0) {
            output << "    \"protein_hit_rate\": " << (double)enhanced_stats.total_protein_matches / enhanced_stats.reads_passed_bloom << ",\n";
        }
        if (enable_kmer_match) {
            output << "    \"total_candidates\": " << enhanced_stats.total_candidates_found << ",\n";
            output << "    \"mutations_found\": " << enhanced_stats.total_mutations_found << ",\n";
        }
        if (enable_translated_search) {
            output << "    \"protein_resistance_found\": " << enhanced_stats.total_protein_resistance_found << ",\n";
            output << "    \"protein_kmer_size\": 5,\n";
            if (enable_smith_waterman) {
                output << "    \"smith_waterman_alignments\": " << enhanced_stats.total_smith_waterman_alignments << ",\n";
            }
        }
        output << "    \"processing_time_seconds\": " << total_time << "\n";
        output << "  }\n";
        output << "}\n";
        output.close();
        
        // Generate final diagnostic report
        if (enable_diagnostic_reporting && diagnostic_reporter) {
            update_pipeline_statistics(
                diagnostic_reporter,
                enhanced_stats.total_reads_processed,
                enhanced_stats.reads_passed_bloom,
                enhanced_stats.reads_with_candidates,
                enhanced_stats.total_protein_matches,
                enhanced_stats.resistance_mutations_found
            );
            
            generate_diagnostic_report(diagnostic_reporter);
            std::cout << "\nDiagnostic report generated. Check *_diagnostic.txt file for detailed analysis." << std::endl;
        }
        
        // Enhanced final summary disabled for production
        /*
        std::cout << "\n=== ENHANCED PROCESSING COMPLETE ===" << std::endl;
        std::cout << "Total reads: " << enhanced_stats.total_reads_processed << std::endl;
        
        if (enable_bloom_filter || enable_kmer_match) {
            std::cout << "Filter performance:" << std::endl;
            if (enable_bloom_filter) {
                std::cout << "  Bloom filter retention: " << enhanced_stats.reads_passed_bloom 
                          << " (" << (100.0 * enhanced_stats.reads_passed_bloom / enhanced_stats.total_reads_processed) << "%)" << std::endl;
            }
            if (enable_kmer_match) {
                std::cout << "  Reads with k-mer hits: " << enhanced_stats.reads_with_kmer_hits << std::endl;
            }
        }
        
        if (enable_translated_search) {
            std::cout << "Protein analysis:" << std::endl;
            std::cout << "  Total protein matches: " << enhanced_stats.total_protein_matches << std::endl;
            std::cout << "  Reads with Smith-Waterman alignments: " << enhanced_stats.reads_with_sw_alignments << std::endl;
            std::cout << "  QRDR region alignments: " << enhanced_stats.qrdr_alignments << std::endl;
            std::cout << "  High confidence matches: " << enhanced_stats.high_confidence_matches << std::endl;
            std::cout << "  Resistance mutations detected: " << enhanced_stats.resistance_mutations_found << std::endl;
        }
        
        // Display top 20 QRDR alignments ranked by alignment score
        if (!qrdr_alignments.empty() && enable_smith_waterman) {
            // Sort QRDR alignments by alignment score (descending)
            std::sort(qrdr_alignments.begin(), qrdr_alignments.end(),
                     [](const ProteinMatch& a, const ProteinMatch& b) {
                         return a.alignment_score > b.alignment_score;
                     });
            
            std::cout << "\nTop 20 QRDR alignments (sorted by alignment score):" << std::endl;
            std::cout << "These represent alignments within Quinolone Resistance Determining Regions:" << std::endl;
            
            int num_qrdr_to_show = std::min(20, (int)qrdr_alignments.size());
            for (int i = 0; i < num_qrdr_to_show; i++) {
                const auto& match = qrdr_alignments[i];
                
                // Map gene IDs to names using dynamic mapping from metadata
                std::string gene_name;
                if (gene_id_to_name.count(match.gene_id)) {
                    gene_name = gene_id_to_name[match.gene_id];
                } else {
                    gene_name = "Gene" + std::to_string(match.gene_id);
                }
                
                // Map species IDs to names using dynamic mapping from metadata
                std::string species_name;
                if (species_id_to_name.count(match.species_id)) {
                    species_name = species_id_to_name[match.species_id];
                } else {
                    species_name = "Species" + std::to_string(match.species_id);
                }
                
                std::cout << "  " << (i + 1) << ". Read " << match.read_id 
                         << ", " << species_name << " " << gene_name 
                         << ", Frame " << (int)match.frame
                         << ", Score: " << std::fixed << std::setprecision(1) << match.alignment_score
                         << ", Identity: " << std::fixed << std::setprecision(1) << (match.identity * 100) << "%"
                         << ", Position: " << match.ref_start << "-" << (match.ref_start + match.match_length - 1)
                         << ", Length: " << match.match_length
                         << ", Peptide: " << match.query_peptide << std::endl;
            }
            
            std::cout << "\nTotal QRDR alignments found: " << qrdr_alignments.size() << std::endl;
            std::cout << "Full QRDR alignment details exported to: " << sample_name << "_qrdr_alignments.csv" << std::endl;
            std::cout << "NOTE: These peptide sequences should be verified using BLAST to confirm alignment accuracy." << std::endl;
        } else if (enable_smith_waterman) {
            std::cout << "\nNo QRDR alignments detected." << std::endl;
            std::cout << "This may indicate:" << std::endl;
            std::cout << "  - Reads don't contain QRDR sequences" << std::endl;
            std::cout << "  - Alignment score thresholds are too stringent" << std::endl;
            std::cout << "  - QRDR region definitions need adjustment" << std::endl;
        }
        
        if (enhanced_stats.resistance_mutations_found > 0) {
            std::cout << "\n🚨 RESISTANCE DETECTED: " << enhanced_stats.resistance_mutations_found << " mutations found" << std::endl;
            std::cout << "   Recommendation: Consider alternative antimicrobials" << std::endl;
        } else if (enhanced_stats.qrdr_alignments > 0) {
            std::cout << "\n✅ QRDR regions covered but no resistance mutations detected" << std::endl;
            std::cout << "   Note: Check diagnostic report for alignment quality assessment" << std::endl;
        } else {
            std::cout << "\n❓ No QRDR coverage detected" << std::endl;
            if (enable_diagnostic_reporting) {
                std::cout << "   Check diagnostic report for troubleshooting advice" << std::endl;
            }
        }
        
        std::cout << "\nTotal time: " << total_time << " seconds" << std::endl;
        std::cout << "Throughput: " << (enhanced_stats.total_reads_processed / (double)total_time) << " reads/second" << std::endl;
        std::cout << "\nHDF5 output: " << hdf5_path << std::endl;
        if (enable_diagnostic_reporting) {
            std::cout << "Check diagnostic report for detailed analysis and troubleshooting" << std::endl;
        }
        */ // End of summary output commenting
    }
};

// Updated main function with new flags
int main(int argc, char** argv) {
    DEBUG_PRINT("=== Enhanced FQ Pipeline GPU Starting (v0.4.0) ===");
    DEBUG_PRINT("Features: Optional Bloom filter + Optional K-mer enrichment + Enhanced protein search + Diagnostic reporting");
    
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <index_path> <reads_R1.fastq.gz> <reads_R2.fastq.gz> [options]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --enable-bloom-filter          Enable bloom filter pre-screening (default: on)\n";
        std::cerr << "  --disable-bloom-filter         Disable bloom filter pre-screening\n";
        std::cerr << "  --enable-kmer-match            Enable k-mer matching (default: on)\n";
        std::cerr << "  --disable-kmer-match           Disable k-mer matching\n";
        std::cerr << "  --enable-translated-search     Enable 6-frame translated search with 5-mer k-mers\n";
        std::cerr << "  --enable-smith-waterman        Enable Smith-Waterman for high-scoring protein matches\n";
        std::cerr << "  --enable-diagnostic-reporting  Enable detailed diagnostic reporting\n";
        std::cerr << "  --protein-db <path>            Path to protein resistance database\n";
        std::cerr << "\nNote: Output files will be automatically created based on sample name extracted from input files\n";
        return 1;
    }
    
    std::string index_path = argv[1];
    std::string r1_path = argv[2];
    std::string r2_path = argv[3];
    std::string output_path = "";  // No longer used, kept for compatibility
    
    // Parse enhanced command line options
    bool enable_bloom_filter = true;       // Default: enabled
    bool enable_kmer_match = true;         // Default: enabled
    bool enable_translated_search = false;
    bool enable_smith_waterman = false;
    bool enable_diagnostic_reporting = false;
    std::string protein_db_path;
    
    for (int i = 4; i < argc; i++) {
        if (std::string(argv[i]) == "--enable-bloom-filter") {
            enable_bloom_filter = true;
        } else if (std::string(argv[i]) == "--disable-bloom-filter") {
            enable_bloom_filter = false;
        } else if (std::string(argv[i]) == "--enable-kmer-match") {
            enable_kmer_match = true;
        } else if (std::string(argv[i]) == "--disable-kmer-match") {
            enable_kmer_match = false;
        } else if (std::string(argv[i]) == "--enable-translated-search") {
            enable_translated_search = true;
        } else if (std::string(argv[i]) == "--enable-smith-waterman") {
            enable_smith_waterman = true;
        } else if (std::string(argv[i]) == "--enable-diagnostic-reporting") {
            enable_diagnostic_reporting = true;
        } else if (std::string(argv[i]) == "--protein-db" && i + 1 < argc) {
            protein_db_path = argv[i + 1];
            i++; // Skip next argument
        }
    }
    
    // Display configuration
    std::cout << "=== Enhanced Fluoroquinolone Resistance Detection Pipeline (v0.4.0) ===" << std::endl;
    std::cout << "Features:";
    if (enable_bloom_filter) std::cout << " Bloom filter";
    if (enable_kmer_match) std::cout << " + K-mer enrichment";
    if (enable_translated_search) {
        std::cout << " + Enhanced 5-mer protein search";
        if (enable_smith_waterman) {
            std::cout << " + Smith-Waterman";
        }
    }
    if (enable_diagnostic_reporting) std::cout << " + Diagnostic reporting";
    std::cout << std::endl;
    
    std::cout << "Index: " << index_path << std::endl;
    std::cout << "R1 reads: " << r1_path << std::endl;
    std::cout << "R2 reads: " << r2_path << std::endl;
    std::cout << "Bloom filter: " << (enable_bloom_filter ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "K-mer matching: " << (enable_kmer_match ? "ENABLED" : "DISABLED") << std::endl;
    
    if (enable_translated_search) {
        std::cout << "Enhanced translated search: ENABLED (5-mer k-mers)" << std::endl;
        if (enable_smith_waterman) {
            std::cout << "Smith-Waterman alignment: ENABLED (for high-scoring matches)" << std::endl;
        }
        if (!protein_db_path.empty()) {
            std::cout << "Protein database: " << protein_db_path << std::endl;
        }
    }
    
    // Check input files exist
    DEBUG_PRINT("Checking input files...");
    std::ifstream idx_check(index_path);
    if (!idx_check.good()) {
        std::cerr << "ERROR: Index file not found: " << index_path << std::endl;
        return 1;
    }
    idx_check.close();
    
    std::ifstream r1_check(r1_path);
    if (!r1_check.good()) {
        std::cerr << "ERROR: R1 file not found: " << r1_path << std::endl;
        return 1;
    }
    r1_check.close();
    
    std::ifstream r2_check(r2_path);
    if (!r2_check.good()) {
        std::cerr << "ERROR: R2 file not found: " << r2_path << std::endl;
        return 1;
    }
    r2_check.close();
    
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
    std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "Memory: " << prop.totalGlobalMem / (1024*1024*1024) << " GB" << std::endl;
    std::cout << std::endl;
    
    // Run enhanced pipeline
    try {
        DEBUG_PRINT("Creating enhanced pipeline instance");
        EnhancedFQResistancePipeline pipeline(enable_translated_search, enable_smith_waterman, 
                                            enable_diagnostic_reporting, enable_bloom_filter, 
                                            enable_kmer_match);
        
        DEBUG_PRINT("Loading index");
        pipeline.loadIndex(index_path);
        
        // Load protein database if provided
        if (enable_translated_search && !protein_db_path.empty()) {
            pipeline.setProteinDatabase(protein_db_path);
        }
        
        DEBUG_PRINT("Starting enhanced paired read processing");
        pipeline.processPairedReads(r1_path, r2_path, output_path);
        
        std::cout << "\n=== PIPELINE COMPLETED SUCCESSFULLY ===" << std::endl;
        std::cout << "Output files are generated in the sample-specific directory" << std::endl;
        std::cout << "based on the input filename (extracted from R1 filename)" << std::endl;
        
        DEBUG_PRINT("=== Enhanced FQ Pipeline GPU Completed Successfully ===");
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Exception caught: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "ERROR: Unknown exception caught" << std::endl;
        return 1;
    }
    
    return 0;
}