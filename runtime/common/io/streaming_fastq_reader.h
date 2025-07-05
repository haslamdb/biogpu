// runtime/common/io/streaming_fastq_reader.h
#ifndef STREAMING_FASTQ_READER_H
#define STREAMING_FASTQ_READER_H

#include <fstream>
#include <memory>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <zlib.h>
#include "sequence_batch.h"

namespace BioGPU {

class StreamingFastqReader {
private:
    // File handles
    std::ifstream file_stream;
    gzFile gz_file = nullptr;
    bool is_gzipped = false;
    bool is_open = false;
    
    // Paired-end support
    std::ifstream file_stream_r2;
    gzFile gz_file_r2 = nullptr;
    bool has_r2 = false;
    
    // Configuration
    size_t batch_size;
    size_t max_queue_size;
    
    // Async reading
    std::thread reader_thread;
    std::queue<std::shared_ptr<SequenceBatch>> batch_queue;
    mutable std::mutex queue_mutex;
    std::condition_variable queue_cv;
    bool stop_reading = false;
    bool reading_complete = false;
    
    // Statistics
    size_t total_reads = 0;
    size_t total_bases = 0;
    
public:
    StreamingFastqReader(size_t batch_size = 100000, size_t max_queue_size = 10)
        : batch_size(batch_size), max_queue_size(max_queue_size) {}
    
    ~StreamingFastqReader() {
        close();
    }
    
    // Open single or paired-end files
    bool open(const std::string& filename);
    bool openPaired(const std::string& r1_file, const std::string& r2_file);
    void close();
    
    // Get next batch (blocking)
    std::shared_ptr<SequenceBatch> getNextBatch();
    
    // Check if more batches available
    bool hasNext() const {
        std::lock_guard<std::mutex> lock(queue_mutex);
        return !batch_queue.empty() || !reading_complete;
    }
    
    // Get statistics
    size_t getTotalReads() const { return total_reads; }
    size_t getTotalBases() const { return total_bases; }
    
private:
    void readerThreadFunc();
    bool readBatchFromFile(SequenceBatch& batch);
    bool readLine(std::string& line);
    bool readLineR2(std::string& line);
};

} // namespace BioGPU

#endif // STREAMING_FASTQ_READER_H