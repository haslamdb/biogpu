#include "fastq_reader.h"
#include <cstring> // For memcpy

bool read_batch_paired(gzFile r1_file, gzFile r2_file, 
                       std::vector<FastqRecord>& batch_r1, 
                       std::vector<FastqRecord>& batch_r2, 
                       int max_batch_size) {
    char buffer[1024];
    
    for (int i = 0; i < max_batch_size; ++i) {
        FastqRecord rec1, rec2;
        
        // Read R1 record
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return i > 0; // End of file
        rec1.header = std::string(buffer);
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return false; // Truncated file
        rec1.sequence = std::string(buffer);
        rec1.sequence.pop_back(); // Remove newline
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return false; // '+' line
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return false; // Quality
        rec1.quality = std::string(buffer);
        rec1.quality.pop_back();

        // Read R2 record
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false; // Paired file mismatch
        rec2.header = std::string(buffer);
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        rec2.sequence = std::string(buffer);
        rec2.sequence.pop_back();
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        rec2.quality = std::string(buffer);
        rec2.quality.pop_back();
        
        batch_r1.push_back(rec1);
        batch_r2.push_back(rec2);
    }
    
    return true;
}

bool read_batch_single(gzFile file, 
                       std::vector<FastqRecord>& batch, 
                       int max_batch_size) {
    char buffer[1024];
    
    for (int i = 0; i < max_batch_size; ++i) {
        FastqRecord rec;
        
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return i > 0; // End of file
        rec.header = std::string(buffer);
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return false; // Truncated
        rec.sequence = std::string(buffer);
        rec.sequence.pop_back(); // Remove newline
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return false; // '+' line
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return false; // Quality
        rec.quality = std::string(buffer);
        rec.quality.pop_back();
        
        batch.push_back(rec);
    }
    
    return true;
}

ReadBatch prepare_batch_for_gpu(const std::vector<FastqRecord>& records) {
    ReadBatch batch;
    batch.num_reads = records.size();
    batch.total_bases = 0;
    
    for (const auto& rec : records) {
        batch.total_bases += rec.sequence.length();
    }
    
    batch.sequences = new char[batch.total_bases];
    batch.lengths = new int[batch.num_reads];
    batch.offsets = new int[batch.num_reads];
    
    int current_offset = 0;
    for (int i = 0; i < batch.num_reads; ++i) {
        const std::string& seq = records[i].sequence;
        memcpy(batch.sequences + current_offset, seq.c_str(), seq.length());
        batch.lengths[i] = seq.length();
        batch.offsets[i] = current_offset;
        current_offset += seq.length();
    }
    
    return batch;
}