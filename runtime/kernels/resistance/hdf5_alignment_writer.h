// hdf5_alignment_writer.h
// HDF5 writer for BioGPU pipeline with translated search support

#ifndef HDF5_ALIGNMENT_WRITER_H
#define HDF5_ALIGNMENT_WRITER_H

#include <H5Cpp.h>
#include <vector>
#include <string>
#include "fq_mutation_detector.cuh"

// Forward declaration for ProteinMatch structure from translated_search.cu
struct ProteinMatch;

class HDF5AlignmentWriter {
private:
    H5::H5File* file;
    std::string filename;
    
    // Buffers for nucleotide alignment batch writing
    std::vector<uint32_t> read_ids;
    std::vector<uint32_t> gene_ids;
    std::vector<uint32_t> species_ids;
    std::vector<uint32_t> seq_ids;
    std::vector<float> alignment_scores;
    std::vector<float> identities;
    std::vector<uint16_t> start_positions;
    std::vector<uint16_t> match_counts;
    std::vector<uint8_t> num_mutations;
    std::vector<uint32_t> mutation_positions;  // Flattened array of mutation positions
    std::vector<uint32_t> mutation_offsets;    // Offset into mutation_positions for each alignment
    
    // NEW: Buffers for translated search results
    std::vector<uint32_t> trans_read_ids;
    std::vector<int8_t> trans_frames;
    std::vector<uint32_t> trans_protein_ids;
    std::vector<uint32_t> trans_gene_ids;
    std::vector<uint32_t> trans_species_ids;
    std::vector<uint16_t> trans_query_starts;
    std::vector<uint16_t> trans_ref_starts;
    std::vector<uint16_t> trans_match_lengths;
    std::vector<uint32_t> trans_codon_starts;
    std::vector<float> trans_alignment_scores;
    std::vector<float> trans_identities;
    std::vector<uint8_t> trans_num_mutations;
    std::vector<float> trans_confidence_scores;
    std::vector<bool> trans_is_resistance;
    
    // Mutation details for translated search
    std::vector<uint16_t> trans_mutation_positions;
    std::vector<char> trans_ref_aas;
    std::vector<char> trans_query_aas;
    std::vector<float> trans_blosum_scores;
    std::vector<uint32_t> trans_mutation_offsets;
    
    // Metadata
    uint64_t total_alignments_written;
    uint64_t total_reads_processed;
    uint64_t total_translated_alignments_written;
    
public:
    HDF5AlignmentWriter(const std::string& output_path);
    ~HDF5AlignmentWriter();
    
    // Initialize HDF5 file structure
    void initialize(const std::string& index_path, 
                   const std::string& r1_path, 
                   const std::string& r2_path);
    
    // Add nucleotide alignment results from GPU
    void addAlignmentBatch(const AlignmentResult* h_results, 
                          size_t num_results,
                          size_t batch_offset);
    
    // Add k-mer screening results (for ML features)
    void addKmerScreeningResults(const CandidateMatch* h_candidates,
                                const uint32_t* h_candidate_counts,
                                size_t num_reads,
                                size_t batch_offset);
    
    // NEW: Add translated search results
    void addTranslatedResults(const ProteinMatch* h_matches,
                             const uint32_t* h_match_counts,
                             size_t num_reads,
                             size_t batch_offset);
    
    // Write all buffered data to disk
    void flush();
    
    // NEW: Flush translated results
    void flushTranslatedResults();
    
    // Finalize and close file
    void finalize(const std::string& json_summary_path);
    
    // Export to SAM format (optional)
    void exportToSAM(const std::string& sam_path);
    
private:
    void createDatasets();
    void createTranslatedDatasets();  // NEW
    void writeMetadata();
    void writeSummaryStats();
    
    // NEW: Helper functions for translated search
    float calculateResistanceConfidence(const ProteinMatch& match);
    bool checkIfResistanceMutation(const ProteinMatch& match);
};

#endif // HDF5_ALIGNMENT_WRITER_H