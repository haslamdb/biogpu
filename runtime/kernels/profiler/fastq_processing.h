// fastq_processing.h
#ifndef FASTQ_PROCESSING_H
#define FASTQ_PROCESSING_H

#include <string>
#include <memory>
#include <functional>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <zlib.h>
#include "minimizer_common.h"
#include "minimizer_extractor.h"

namespace biogpu {

// Thread-safe queue for batch processing
template<typename T>
class ThreadSafeQueue {
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
public:
    ThreadSafeQueue();
    ~ThreadSafeQueue();
    
    void push(T item);
    bool pop(T& item);
    void set_finished();
};

// FASTQ Reader class
class FastqReader {
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
public:
    FastqReader(const std::string& filename);
    ~FastqReader();
    
    bool read_batch(ReadBatch& batch, size_t batch_size);
    bool is_open() const;
    void close();
};

// GPU Pipeline for minimizer extraction
class GPUMinimizerPipeline {
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
public:
    GPUMinimizerPipeline(int k = 31, int m = 15, 
                         size_t batch_size = 10000, 
                         int gpu_threads = 1);
    ~GPUMinimizerPipeline();
    
    // Process a FASTQ file with a callback for each batch
    void process_file(const std::string& fastq_file, 
                      std::function<void(const ReadBatch&, 
                                       const std::vector<std::vector<Minimizer>>&)> callback);
    
    // Process a single batch
    std::vector<std::vector<Minimizer>> process_batch(const ReadBatch& batch);
    
    // Get pipeline statistics
    struct Statistics {
        size_t total_reads;
        size_t total_batches;
        double processing_time_ms;
    };
    
    Statistics get_statistics() const;
};

// Statistics collector for minimizer analysis
class MinimizerStats {
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
public:
    MinimizerStats();
    ~MinimizerStats();
    
    void process_batch(const ReadBatch& batch, 
                       const std::vector<std::vector<Minimizer>>& minimizers);
    
    void print_summary() const;
    
    // Get specific statistics
    size_t get_total_reads() const;
    size_t get_total_minimizers() const;
    size_t get_unique_minimizers() const;
    double get_average_minimizers_per_read() const;
    std::vector<std::pair<uint64_t, uint32_t>> get_top_minimizers(size_t n = 10) const;
};

} // namespace biogpu

#endif // FASTQ_PROCESSING_H