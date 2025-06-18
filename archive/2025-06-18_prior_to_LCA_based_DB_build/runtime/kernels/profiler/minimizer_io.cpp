// minimizer_io.cpp - Implementation of minimizer I/O functionality
#include "minimizer_io.h"
#include "fastq_processing.h"
#include <iostream>
#include <cstring>

namespace biogpu {

MinimizerWriter::MinimizerWriter(const std::string& filename, int k, int m) 
    : file(filename, std::ios::binary) {
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }
    
    header.k_mer_size = k;
    header.window_size = m;
    
    // Write header placeholder - will update when closing
    file.write(reinterpret_cast<const char*>(&header), sizeof(header));
}

MinimizerWriter::~MinimizerWriter() {
    if (file.is_open()) {
        close();
    }
}

void MinimizerWriter::write_batch(const std::vector<std::string>& read_ids,
                                 const std::vector<std::vector<Minimizer>>& minimizers,
                                 const std::vector<uint32_t>& read_lengths) {
    if (read_ids.size() != minimizers.size() || read_ids.size() != read_lengths.size()) {
        throw std::runtime_error("Mismatched batch sizes");
    }
    
    for (size_t i = 0; i < read_ids.size(); i++) {
        ReadRecord record;
        record.read_id = reads_written + i;
        record.num_minimizers = minimizers[i].size();
        record.read_length = read_lengths[i];
        record.reserved = 0;
        
        // Write read record
        file.write(reinterpret_cast<const char*>(&record), sizeof(record));
        
        // Write minimizers for this read
        if (record.num_minimizers > 0) {
            file.write(reinterpret_cast<const char*>(minimizers[i].data()), 
                      record.num_minimizers * sizeof(Minimizer));
        }
        
        minimizers_written += record.num_minimizers;
    }
    
    reads_written += read_ids.size();
}

void MinimizerWriter::close() {
    if (!file.is_open()) return;
    
    // Update header with final counts
    header.num_reads = reads_written;
    header.total_minimizers = minimizers_written;
    
    // Seek to beginning and rewrite header
    file.seekp(0);
    file.write(reinterpret_cast<const char*>(&header), sizeof(header));
    
    file.close();
    
    std::cout << "Wrote " << reads_written << " reads with " 
              << minimizers_written << " total minimizers\n";
}

MinimizerReader::MinimizerReader(const std::string& filename) 
    : file(filename, std::ios::binary) {
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open input file: " + filename);
    }
    
    // Read header
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    
    if (header.magic != 0x4D494E49) {
        throw std::runtime_error("Invalid minimizer file format");
    }
    
    if (header.version != 1) {
        throw std::runtime_error("Unsupported file version");
    }
}

MinimizerReader::~MinimizerReader() {
    if (file.is_open()) {
        file.close();
    }
}

bool MinimizerReader::read_next(ReadRecord& record, std::vector<Minimizer>& minimizers) {
    if (!file.read(reinterpret_cast<char*>(&record), sizeof(record))) {
        return false;  // End of file
    }
    
    minimizers.resize(record.num_minimizers);
    if (record.num_minimizers > 0) {
        file.read(reinterpret_cast<char*>(minimizers.data()), 
                  record.num_minimizers * sizeof(Minimizer));
    }
    
    return true;
}

void MinimizerReader::read_all(std::vector<ReadRecord>& records,
                              std::vector<std::vector<Minimizer>>& all_minimizers) {
    records.clear();
    all_minimizers.clear();
    
    ReadRecord record;
    std::vector<Minimizer> minimizers;
    
    while (read_next(record, minimizers)) {
        records.push_back(record);
        all_minimizers.push_back(std::move(minimizers));
    }
}

void save_minimizers_from_fastq(const std::string& fastq_file,
                               const std::string& output_file,
                               int k, int m,
                               size_t batch_size,
                               int gpu_threads) {
    // Create pipeline and writer
    GPUMinimizerPipeline pipeline(k, m, batch_size, gpu_threads);
    MinimizerWriter writer(output_file, k, m);
    
    // Process file and save minimizers
    std::cout << "Extracting and saving minimizers from " << fastq_file << "...\n";
    
    size_t total_reads = 0;
    pipeline.process_file(fastq_file, 
        [&writer, &total_reads](const ReadBatch& batch, 
                               const std::vector<std::vector<Minimizer>>& minimizers) {
            // Extract read lengths
            std::vector<uint32_t> read_lengths;
            read_lengths.reserve(batch.size());
            for (const auto& seq : batch.sequences) {
                read_lengths.push_back(seq.length());
            }
            
            // Write batch to file
            writer.write_batch(batch.headers, minimizers, read_lengths);
            
            total_reads += batch.size();
            if (total_reads % 100000 == 0) {
                std::cout << "Processed and saved " << total_reads << " reads...\n";
            }
        });
    
    writer.close();
    
    auto stats = pipeline.get_statistics();
    std::cout << "\nMinimizer extraction complete:\n";
    std::cout << "  Total reads: " << stats.total_reads << "\n";
    std::cout << "  Processing time: " << stats.processing_time_ms / 1000.0 << " seconds\n";
    std::cout << "  Output file: " << output_file << "\n";
}

} // namespace biogpu