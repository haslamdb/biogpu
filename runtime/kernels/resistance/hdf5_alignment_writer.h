#ifndef HDF5_ALIGNMENT_WRITER_H
#define HDF5_ALIGNMENT_WRITER_H

#include <H5Cpp.h>
#include <vector>
#include <string>
#include "fq_mutation_detector.cuh"

class HDF5AlignmentWriter {
private:
    H5::H5File* file;
    std::string filename;
    
    // Buffers for batch writing
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
    
    // Metadata
    uint64_t total_alignments_written;
    uint64_t total_reads_processed;
    
public:
    HDF5AlignmentWriter(const std::string& output_path);
    ~HDF5AlignmentWriter();
    
    // Initialize HDF5 file structure
    void initialize(const std::string& index_path, 
                   const std::string& r1_path, 
                   const std::string& r2_path);
    
    // Add alignment results from GPU
    void addAlignmentBatch(const AlignmentResult* h_results, 
                          size_t num_results,
                          size_t batch_offset);
    
    // Add k-mer screening results (for ML features)
    void addKmerScreeningResults(const CandidateMatch* h_candidates,
                                const uint32_t* h_candidate_counts,
                                size_t num_reads,
                                size_t batch_offset);
    
    // Write all buffered data to disk
    void flush();
    
    // Finalize and close file
    void finalize(const std::string& json_summary_path);
    
    // Export to SAM format (optional)
    void exportToSAM(const std::string& sam_path);
    
private:
    void createDatasets();
    void writeMetadata();
    void writeSummaryStats();
};

#endif // HDF5_ALIGNMENT_WRITER_H