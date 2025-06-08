// biogpu_core.h - Core data structures for BioGPU DSL
#ifndef BIOGPU_CORE_H
#define BIOGPU_CORE_H

#include <cstdint>
#include <vector>
#include <map>
#include <string>

namespace biogpu {

// =============================================================================
// SEQUENCE READ DATA STRUCTURES
// =============================================================================

// 2-bit encoding for DNA sequences (A=00, C=01, G=10, T=11)
struct CompressedReadBatch {
    // Core sequence data
    uint32_t* packed_sequences;      // 16 bases per uint32_t
    uint8_t* quality_scores;          // Binned Phred scores (0-40 -> 0-255)
    uint32_t* read_metadata;          // Packed: read_id(20b) + length(12b)
    
    // Paired-end support
    uint32_t* pair_offsets;           // Offset to mate pair
    uint8_t* pair_orientation;        // FR, RF, FF, RR
    
    // Pre-computed features for fast filtering/searching
    uint64_t* minimizer_sketches;     // Pre-computed minimizers for each read
    uint16_t* gc_content;             // GC percentage * 100 (0-10000)
    uint32_t* kmer_hashes;            // Key k-mers for indexing
    
    // Batch metadata
    size_t num_reads;
    size_t max_read_length;
    size_t batch_id;
    
    // GPU memory management
    bool on_device;
    int device_id;
};

// Wavelet tree for efficient pattern matching across reads
struct WaveletReadTree {
    struct Level {
        uint64_t* bitvector;
        size_t size;
        uint32_t* rank_samples;    // Sampled rank values for O(1) query
    };
    
    std::vector<Level> levels;
    std::map<char, int> alphabet_mapping;
    size_t total_length;
    
    // GPU kernels for operations
    __device__ int rank(char symbol, size_t position);
    __device__ int select(char symbol, int occurrence);
    __device__ void extract_pattern_matches(const char* pattern, int* matches);
};

// Memory-efficient storage for very large datasets
struct MemoryMappedReads {
    void* mmap_ptr;
    size_t file_size;
    
    struct ReadIndex {
        uint64_t* read_offsets;      // File offset for each read
        uint32_t* read_lengths;
        uint64_t num_reads;
        
        // Optional indices for fast access
        std::map<uint64_t, std::vector<uint32_t>> minimizer_to_reads;
        std::map<std::string, std::vector<uint32_t>> sample_to_reads;
    } index;
    
    // Load specific batch into GPU memory
    CompressedReadBatch load_batch(size_t start_read, size_t num_reads);
};

// =============================================================================
// RESISTANCE DATABASE STRUCTURES
// =============================================================================

// Taxonomic information for organisms
struct TaxonomyInfo {
    uint32_t species_id;
    uint32_t genus_id;
    uint32_t family_id;
    std::string species_name;
    std::string strain;
};

// Node representing an organism in the resistance graph
struct OrganismNode {
    TaxonomyInfo taxonomy;
    std::vector<uint32_t> reference_genome_ids;
    uint64_t typical_resistance_profile;    // Bit vector of common resistances
    float clinical_prevalence;              // 0-1, based on clinical data
    
    // Organism-specific metadata
    float typical_generation_time;          // Hours
    bool is_pathogen;
    std::vector<std::string> common_infections;
};

// Node representing a resistance gene
struct GeneNode {
    uint32_t gene_id;
    std::string gene_name;                  // e.g., "gyrA"
    std::string functional_category;        // e.g., "DNA_GYRASE_SUBUNIT_A"
    
    // Sequence information
    std::vector<std::string> sequence_variants;
    std::vector<std::pair<int, int>> conserved_regions;  // Start, end positions
    std::vector<std::pair<int, int>> qrdr_regions;      // QRDR boundaries
    
    // For alignment
    float* position_weight_matrix;          // PWM for this gene family
    int alignment_length;
};

// Specific mutation information
struct MutationNode {
    uint32_t mutation_id;
    uint32_t gene_id;
    
    // Position and change
    int codon_position;                     // Which codon (1-based)
    int nucleotide_position;                // Position in gene (1-based)
    std::string wild_type_codon;            // e.g., "TCT" (Ser)
    std::string mutant_codon;               // e.g., "TTG" (Leu)
    char wild_type_aa;                      // e.g., 'S'
    char mutant_aa;                         // e.g., 'L'
    
    // Resistance information
    struct ResistancePhenotype {
        std::vector<std::string> drugs_affected;  // e.g., ["levofloxacin", "ciprofloxacin"]
        float mic_fold_change;                    // e.g., 8.0 for 8-fold increase
        float confidence_score;                   // 0-1, based on evidence
        std::vector<uint32_t> pubmed_ids;        // Supporting literature
    } phenotype;
    
    // Organism associations
    std::vector<uint32_t> found_in_species;      // Which species have this mutation
    std::map<uint32_t, float> species_frequency; // How common in each species
};

// Main resistance graph structure
struct ResistanceGraph {
    // Nodes
    std::map<uint32_t, OrganismNode> organisms;
    std::map<uint32_t, GeneNode> genes;
    std::map<uint32_t, MutationNode> mutations;
    
    // Edges for efficient traversal
    std::map<uint32_t, std::vector<uint32_t>> organism_to_genes;
    std::map<uint32_t, std::vector<uint32_t>> gene_to_mutations;
    std::map<uint32_t, std::vector<uint32_t>> mutation_to_organisms;
    
    // Fast lookup structures
    std::map<std::string, uint32_t> gene_name_to_id;
    std::map<std::pair<uint32_t, std::string>, uint32_t> gene_mutation_to_id;
};

// GPU-optimized mutation scanning index
struct MutationIndex {
    // Bloom filter cascade for progressive filtering
    struct BloomFilter {
        uint64_t* bit_array;
        size_t size;
        int num_hash_functions;
        int kmer_size;
    };
    std::vector<BloomFilter> bloom_cascade;  // Different k-mer sizes
    
    // Perfect hash for known mutations
    struct MutationSignature {
        uint32_t gene_id;
        uint64_t sequence_hash;
        int position;
    };
    
    // GPU-friendly hash table
    MutationSignature* signatures;
    MutationNode** mutation_pointers;
    size_t table_size;
    
    // Position-specific scoring matrices for each gene
    struct PSSM {
        float* scores;              // 4 x length matrix (ACGT)
        int length;
        float threshold;
    };
    std::map<uint32_t, PSSM> gene_pssms;
    
    // Codon-aware indexing
    struct CodonIndex {
        uint32_t gene_id;
        int* codon_starts;          // Start position of each codon
        int num_codons;
        int* important_codons;      // Which codons to check for resistance
        int num_important;
    };
    std::vector<CodonIndex> codon_indices;
};

// =============================================================================
// ANALYSIS RESULTS STRUCTURES
// =============================================================================

// Results from mutation detection
struct MutationHit {
    uint32_t read_id;
    uint32_t mutation_id;
    float quality_score;            // Confidence in the call
    int position_in_read;           // Where mutation was found
    bool is_perfect_match;
    
    // For novel variants
    std::string observed_sequence;
    float similarity_score;
};

// Organism abundance profile
struct CommunityProfile {
    std::map<uint32_t, float> species_abundance;     // Species ID -> relative abundance
    std::map<uint32_t, float> species_coverage;      // Species ID -> genome coverage
    std::map<uint32_t, uint64_t> species_read_count; // Raw read counts
    
    float total_reads;
    float shannon_diversity;
    float evenness;
};

// Resistance profile for a sample
struct ResistanceProfile {
    // Per-organism resistance
    struct OrganismResistance {
        uint32_t species_id;
        std::vector<MutationHit> mutations_found;
        float resistance_score;              // 0-1, overall resistance level
        std::map<std::string, float> drug_resistance;  // Drug -> predicted MIC change
    };
    std::vector<OrganismResistance> organism_profiles;
    
    // Summary statistics
    float overall_fq_resistance;             // Fraction of FQ-resistant organisms
    std::map<uint32_t, float> mutation_frequencies;  // How common each mutation
    
    // Clinical interpretation
    std::vector<std::string> resistance_warnings;
    std::vector<std::string> treatment_recommendations;
};

// =============================================================================
// GPU KERNEL DECLARATIONS
// =============================================================================

// Core sequence operations
__global__ void compress_reads_kernel(const char* raw_sequences, CompressedReadBatch* output);
__global__ void build_minimizer_index(CompressedReadBatch* reads, uint64_t* minimizers);
__global__ void compute_kmer_hashes(CompressedReadBatch* reads, int k, uint32_t* hashes);

// Mutation detection
__global__ void scan_mutations_bloom(CompressedReadBatch* reads, MutationIndex* index, bool* potential_hits);
__global__ void verify_mutations_exact(CompressedReadBatch* reads, MutationIndex* index, MutationHit* results);
__global__ void detect_novel_variants(CompressedReadBatch* reads, PSSM* gene_pssms, MutationHit* novel);

// Community profiling
__global__ void map_reads_to_species(CompressedReadBatch* reads, ReferenceIndex* refs, uint32_t* assignments);
__global__ void calculate_abundances(uint32_t* assignments, CommunityProfile* profile);

// Resistance attribution
__global__ void link_mutations_to_organisms(MutationHit* mutations, uint32_t* read_species, ResistanceProfile* profile);

} // namespace biogpu

#endif // BIOGPU_CORE_H