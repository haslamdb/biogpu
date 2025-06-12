// minimizer_io.h - Fast I/O for minimizer data
#ifndef MINIMIZER_IO_H
#define MINIMIZER_IO_H

#include <string>
#include <vector>
#include <fstream>
#include "minimizer_common.h"

namespace biogpu {

// Binary format header
struct MinimizerFileHeader {
    uint32_t magic = 0x4D494E49;  // "MINI" in hex
    uint32_t version = 1;
    uint32_t k_mer_size;
    uint32_t window_size;
    uint64_t num_reads;
    uint64_t total_minimizers;
};

// Read record in file - stores metadata for each read
struct ReadRecord {
    uint32_t read_id;
    uint32_t num_minimizers;
    uint32_t read_length;
    uint32_t reserved;  // For alignment
};

class MinimizerWriter {
private:
    std::ofstream file;
    MinimizerFileHeader header;
    uint64_t reads_written = 0;
    uint64_t minimizers_written = 0;
    
public:
    MinimizerWriter(const std::string& filename, int k, int m);
    ~MinimizerWriter();
    
    // Write minimizers for a batch of reads
    void write_batch(const std::vector<std::string>& read_ids,
                     const std::vector<std::vector<Minimizer>>& minimizers,
                     const std::vector<uint32_t>& read_lengths);
    
    // Finalize and close file
    void close();
};

class MinimizerReader {
private:
    std::ifstream file;
    MinimizerFileHeader header;
    
public:
    MinimizerReader(const std::string& filename);
    ~MinimizerReader();
    
    // Get file metadata
    const MinimizerFileHeader& get_header() const { return header; }
    
    // Read all minimizers for a specific read
    bool read_next(ReadRecord& record, std::vector<Minimizer>& minimizers);
    
    // Read entire file into memory (for smaller files)
    void read_all(std::vector<ReadRecord>& records,
                  std::vector<std::vector<Minimizer>>& all_minimizers);
};

// Utility function to save minimizers from a FASTQ file
void save_minimizers_from_fastq(const std::string& fastq_file,
                               const std::string& output_file,
                               int k = 31,
                               int m = 15,
                               size_t batch_size = 10000,
                               int gpu_threads = 2);

} // namespace biogpu

#endif // MINIMIZER_IO_H