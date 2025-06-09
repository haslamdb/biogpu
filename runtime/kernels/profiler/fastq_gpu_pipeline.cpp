// fastq_gpu_pipeline.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <zlib.h>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <unordered_map>
#include <algorithm>
#include "minimizer_common.h"

// Forward declarations
class MinimizerExtractor {
public:
    MinimizerExtractor(int k_size = 31, int window_size = 15);
    ~MinimizerExtractor();
    std::vector<std::vector<Minimizer>> extract_minimizers(const std::vector<std::string>& sequences);
private:
    void allocate_device_memory(size_t num_reads, size_t total_sequence_length);
    int k;
    int m;
    size_t allocated_reads;
    void* d_sequences;
    void* d_sequence_offsets;
    void* d_sequence_lengths;
    void* d_minimizers;
    void* d_minimizer_counts;
};

// Thread-safe queue for batches
template<typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue;
    mutable std::mutex mutex;
    std::condition_variable cond;
    bool finished = false;
    
public:
    void push(T item) {
        {
            std::lock_guard<std::mutex> lock(mutex);
            queue.push(std::move(item));
        }
        cond.notify_one();
    }
    
    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mutex);
        cond.wait(lock, [this] { return !queue.empty() || finished; });
        
        if (queue.empty()) return false;
        
        item = std::move(queue.front());
        queue.pop();
        return true;
    }
    
    void set_finished() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            finished = true;
        }
        cond.notify_all();
    }
};

// FASTQ Reader class
class FastqReader {
private:
    std::string filename;
    gzFile gz_file;
    bool is_gzipped;
    std::ifstream text_file;
    
public:
    FastqReader(const std::string& fname) : filename(fname) {
        // Check if gzipped
        is_gzipped = (filename.size() > 3 && 
                      filename.substr(filename.size() - 3) == ".gz");
        
        if (is_gzipped) {
            gz_file = gzopen(filename.c_str(), "r");
            if (!gz_file) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
        } else {
            text_file.open(filename);
            if (!text_file.is_open()) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
        }
    }
    
    ~FastqReader() {
        if (is_gzipped && gz_file) {
            gzclose(gz_file);
        }
    }
    
    bool read_batch(ReadBatch& batch, size_t batch_size) {
        batch.clear();
        batch.reserve(batch_size);
        
        std::string line;
        size_t count = 0;
        
        while (count < batch_size) {
            // Read header
            if (!getline(line)) break;
            if (line.empty() || line[0] != '@') continue;
            
            std::string header = line.substr(1);  // Remove @
            
            // Read sequence
            if (!getline(line)) break;
            std::string sequence = line;
            
            // Read plus line
            if (!getline(line)) break;
            
            // Read quality
            if (!getline(line)) break;
            std::string quality = line;
            
            batch.headers.push_back(std::move(header));
            batch.sequences.push_back(std::move(sequence));
            batch.qualities.push_back(std::move(quality));
            
            count++;
        }
        
        return batch.size() > 0;
    }
    
private:
    bool getline(std::string& line) {
        if (is_gzipped) {
            char buffer[4096];
            if (gzgets(gz_file, buffer, sizeof(buffer)) == NULL) {
                return false;
            }
            line = buffer;
            // Remove newline
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            return true;
        } else {
            return std::getline(text_file, line).good();
        }
    }
};

// GPU Pipeline class
class GPUMinimizerPipeline {
private:
    MinimizerExtractor extractor;
    ThreadSafeQueue<std::shared_ptr<ReadBatch>> input_queue;
    ThreadSafeQueue<std::pair<std::shared_ptr<ReadBatch>, 
                              std::vector<std::vector<Minimizer>>>> output_queue;
    
    std::vector<std::thread> gpu_workers;
    size_t batch_size;
    int num_gpu_threads;
    
public:
    GPUMinimizerPipeline(int k = 31, int m = 15, 
                         size_t batch_sz = 10000, 
                         int gpu_threads = 1) 
        : extractor(k, m), batch_size(batch_sz), num_gpu_threads(gpu_threads) {
    }
    
    // Process FASTQ file
    void process_file(const std::string& fastq_file, 
                      std::function<void(const ReadBatch&, 
                                       const std::vector<std::vector<Minimizer>>&)> callback) {
        // Start GPU workers
        for (int i = 0; i < num_gpu_threads; i++) {
            gpu_workers.emplace_back([this, i] {
                cudaSetDevice(i % 2);  // Alternate between GPUs
                gpu_worker_thread();
            });
        }
        
        // Start result processor
        std::thread result_processor([this, callback] {
            process_results(callback);
        });
        
        // Read file and queue batches
        FastqReader reader(fastq_file);
        size_t total_reads = 0;
        
        while (true) {
            auto batch = std::make_shared<ReadBatch>();
            if (!reader.read_batch(*batch, batch_size)) {
                break;
            }
            
            total_reads += batch->size();
            input_queue.push(batch);
            
            if (total_reads % 100000 == 0) {
                std::cout << "Queued " << total_reads << " reads...\n";
            }
        }
        
        // Signal completion
        input_queue.set_finished();
        
        // Wait for GPU workers
        for (auto& worker : gpu_workers) {
            worker.join();
        }
        
        output_queue.set_finished();
        result_processor.join();
        
        std::cout << "Processed " << total_reads << " total reads\n";
    }
    
private:
    void gpu_worker_thread() {
        std::shared_ptr<ReadBatch> batch;
        
        while (input_queue.pop(batch)) {
            // Extract minimizers on GPU
            auto minimizers = extractor.extract_minimizers(batch->sequences);
            
            // Queue results
            output_queue.push({batch, std::move(minimizers)});
        }
    }
    
    void process_results(std::function<void(const ReadBatch&, 
                                          const std::vector<std::vector<Minimizer>>&)> callback) {
        std::pair<std::shared_ptr<ReadBatch>, 
                  std::vector<std::vector<Minimizer>>> result;
        
        while (output_queue.pop(result)) {
            callback(*result.first, result.second);
        }
    }
};

// Example: Statistics collector
class MinimizerStats {
private:
    std::mutex mutex;
    size_t total_reads = 0;
    size_t total_minimizers = 0;
    std::unordered_map<uint64_t, uint32_t> minimizer_counts;
    
public:
    void process_batch(const ReadBatch& batch, 
                       const std::vector<std::vector<Minimizer>>& minimizers) {
        std::lock_guard<std::mutex> lock(mutex);
        
        total_reads += batch.size();
        
        for (const auto& read_minimizers : minimizers) {
            total_minimizers += read_minimizers.size();
            
            for (const auto& m : read_minimizers) {
                minimizer_counts[m.hash]++;
            }
        }
    }
    
    void print_summary() {
        std::cout << "\n=== Minimizer Statistics ===\n";
        std::cout << "Total reads: " << total_reads << "\n";
        std::cout << "Total minimizers: " << total_minimizers << "\n";
        std::cout << "Unique minimizers: " << minimizer_counts.size() << "\n";
        std::cout << "Average minimizers per read: " 
                  << (double)total_minimizers / total_reads << "\n";
        
        // Find most common minimizers
        std::vector<std::pair<uint64_t, uint32_t>> top_minimizers;
        for (const auto& [hash, count] : minimizer_counts) {
            top_minimizers.push_back({hash, count});
        }
        
        std::sort(top_minimizers.begin(), top_minimizers.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        std::cout << "\nTop 10 most frequent minimizers:\n";
        for (size_t i = 0; i < std::min(size_t(10), top_minimizers.size()); i++) {
            std::cout << "  Hash " << top_minimizers[i].first 
                      << ": " << top_minimizers[i].second << " occurrences\n";
        }
    }
};

// Main function
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fastq_file>\n";
        return 1;
    }
    
    std::string fastq_file = argv[1];
    
    // Create pipeline with k=31, m=15, batch_size=10000
    GPUMinimizerPipeline pipeline(31, 15, 10000, 2);  // Use 2 GPU threads
    
    // Create statistics collector
    MinimizerStats stats;
    
    // Process file
    std::cout << "Processing " << fastq_file << "...\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    
    pipeline.process_file(fastq_file, 
        [&stats](const ReadBatch& batch, 
                 const std::vector<std::vector<Minimizer>>& minimizers) {
            stats.process_batch(batch, minimizers);
        });
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    stats.print_summary();
    
    std::cout << "\nTotal processing time: " << duration.count() / 1000.0 << " seconds\n";
    
    return 0;
}

// Compilation command:
// nvcc -O3 -std=c++17 fastq_gpu_pipeline.cpp minimizer_extraction.cu -o minimizer_pipeline -lz -lpthread