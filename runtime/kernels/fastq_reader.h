// fastq_reader.h

#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <string>
#include <vector>
#include <zlib.h>

// Represents a single FASTQ record
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Represents a batch of reads ready for GPU transfer
struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

// Reads a batch of paired-end records from two gzipped FASTQ files.
// Returns true if a full or partial batch was read.
bool read_batch_paired(gzFile r1_file, gzFile r2_file, 
                       std::vector<FastqRecord>& batch_r1, 
                       std::vector<FastqRecord>& batch_r2, 
                       int max_batch_size);

// Reads a batch of single-end records from a gzipped FASTQ file.
// Returns true if a full or partial batch was read.
bool read_batch_single(gzFile file, 
                       std::vector<FastqRecord>& batch, 
                       int max_batch_size);

// Prepares a vector of FastqRecord structs into a contiguous memory block for the GPU.
// The caller is responsible for deleting the allocated memory in the returned ReadBatch.
ReadBatch prepare_batch_for_gpu(const std::vector<FastqRecord>& records);

#endif // FASTQ_READER_H
