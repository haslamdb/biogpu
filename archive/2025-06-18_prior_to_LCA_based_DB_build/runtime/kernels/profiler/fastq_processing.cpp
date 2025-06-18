// fastq_processing.cpp
#include "fastq_processing.h"
#include <iostream>
#include <fstream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <algorithm>
#include <atomic>
#include <cuda_runtime.h>

namespace biogpu {

// ThreadSafeQueue implementation
template<typename T>
class ThreadSafeQueue<T>::Impl {
public:
    std::queue<T> queue;
    mutable std::mutex mutex;
    std::condition_variable cond;
    bool finished = false;
};

template<typename T>
ThreadSafeQueue<T>::ThreadSafeQueue() : pImpl(std::make_unique<Impl>()) {}

template<typename T>
ThreadSafeQueue<T>::~ThreadSafeQueue() = default;

template<typename T>
void ThreadSafeQueue<T>::push(T item) {
    {
        std::lock_guard<std::mutex> lock(pImpl->mutex);
        pImpl->queue.push(std::move(item));
    }
    pImpl->cond.notify_one();
}

template<typename T>
bool ThreadSafeQueue<T>::pop(T& item) {
    std::unique_lock<std::mutex> lock(pImpl->mutex);
    pImpl->cond.wait(lock, [this] { return !pImpl->queue.empty() || pImpl->finished; });
    
    if (pImpl->queue.empty()) return false;
    
    item = std::move(pImpl->queue.front());
    pImpl->queue.pop();
    return true;
}

template<typename T>
void ThreadSafeQueue<T>::set_finished() {
    {
        std::lock_guard<std::mutex> lock(pImpl->mutex);
        pImpl->finished = true;
    }
    pImpl->cond.notify_all();
}

// Explicit instantiations for types we use
template class ThreadSafeQueue<std::shared_ptr<ReadBatch>>;
template class ThreadSafeQueue<std::pair<std::shared_ptr<ReadBatch>, 
                                       std::vector<std::vector<Minimizer>>>>;

// FastqReader implementation
class FastqReader::Impl {
public:
    std::string filename;
    gzFile gz_file = nullptr;
    bool is_gzipped = false;
    std::ifstream text_file;
    
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

FastqReader::FastqReader(const std::string& filename) : pImpl(std::make_unique<Impl>()) {
    pImpl->filename = filename;
    
    // Check if gzipped
    pImpl->is_gzipped = (filename.size() > 3 && 
                        filename.substr(filename.size() - 3) == ".gz");
    
    if (pImpl->is_gzipped) {
        pImpl->gz_file = gzopen(filename.c_str(), "r");
        if (!pImpl->gz_file) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    } else {
        pImpl->text_file.open(filename);
        if (!pImpl->text_file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    }
}

FastqReader::~FastqReader() {
    close();
}

void FastqReader::close() {
    if (pImpl->is_gzipped && pImpl->gz_file) {
        gzclose(pImpl->gz_file);
        pImpl->gz_file = nullptr;
    }
    if (pImpl->text_file.is_open()) {
        pImpl->text_file.close();
    }
}

bool FastqReader::is_open() const {
    if (pImpl->is_gzipped) {
        return pImpl->gz_file != nullptr;
    } else {
        return pImpl->text_file.is_open();
    }
}

bool FastqReader::read_batch(ReadBatch& batch, size_t batch_size) {
    batch.clear();
    batch.reserve(batch_size);
    
    std::string line;
    size_t count = 0;
    
    while (count < batch_size) {
        // Read header
        if (!pImpl->getline(line)) break;
        if (line.empty() || line[0] != '@') continue;
        
        std::string header = line.substr(1);  // Remove @
        
        // Read sequence
        if (!pImpl->getline(line)) break;
        std::string sequence = line;
        
        // Read plus line
        if (!pImpl->getline(line)) break;
        
        // Read quality
        if (!pImpl->getline(line)) break;
        std::string quality = line;
        
        batch.headers.push_back(std::move(header));
        batch.sequences.push_back(std::move(sequence));
        batch.qualities.push_back(std::move(quality));
        
        count++;
    }
    
    return batch.size() > 0;
}

// GPUMinimizerPipeline implementation
class GPUMinimizerPipeline::Impl {
public:
    // Parameters for creating extractors
    int k_mer_size;
    int minimizer_window;
    
    ThreadSafeQueue<std::shared_ptr<ReadBatch>> input_queue;
    ThreadSafeQueue<std::pair<std::shared_ptr<ReadBatch>, 
                              std::vector<std::vector<Minimizer>>>> output_queue;
    
    std::vector<std::thread> gpu_workers;
    size_t batch_size;
    int num_gpu_threads;
    
    // Statistics
    std::atomic<size_t> total_reads{0};
    std::atomic<size_t> total_batches{0};
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
    
    Impl(int k, int m, size_t batch_sz, int gpu_threads)
        : k_mer_size(k), minimizer_window(m), batch_size(batch_sz), num_gpu_threads(gpu_threads) {}
    
    void gpu_worker_thread(int thread_id) {
        // Set GPU device for this thread
        int num_devices;
        cudaGetDeviceCount(&num_devices);
        if (num_devices > 1) {
            cudaSetDevice(thread_id % num_devices);
        }
        
        // Create thread-local extractor to avoid race conditions
        MinimizerExtractor extractor(k_mer_size, minimizer_window);
        
        std::shared_ptr<ReadBatch> batch;
        
        while (input_queue.pop(batch)) {
            try {
                // Extract minimizers on GPU
                auto minimizers = extractor.extract_minimizers(batch->sequences);
                
                // Queue results
                output_queue.push({batch, std::move(minimizers)});
                
                total_batches++;
            } catch (const std::exception& e) {
                std::cerr << "GPU worker " << thread_id << " error: " << e.what() << std::endl;
                // Continue processing other batches
            }
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

GPUMinimizerPipeline::GPUMinimizerPipeline(int k, int m, size_t batch_size, int gpu_threads)
    : pImpl(std::make_unique<Impl>(k, m, batch_size, gpu_threads)) {}

GPUMinimizerPipeline::~GPUMinimizerPipeline() = default;

void GPUMinimizerPipeline::process_file(const std::string& fastq_file,
                                        std::function<void(const ReadBatch&, 
                                                         const std::vector<std::vector<Minimizer>>&)> callback) {
    pImpl->start_time = std::chrono::steady_clock::now();
    pImpl->total_reads = 0;
    pImpl->total_batches = 0;
    
    // Start GPU workers
    for (int i = 0; i < pImpl->num_gpu_threads; i++) {
        pImpl->gpu_workers.emplace_back(&Impl::gpu_worker_thread, pImpl.get(), i);
    }
    
    // Start result processor
    std::thread result_processor(&Impl::process_results, pImpl.get(), callback);
    
    // Read file and queue batches
    FastqReader reader(fastq_file);
    
    while (true) {
        auto batch = std::make_shared<ReadBatch>();
        if (!reader.read_batch(*batch, pImpl->batch_size)) {
            break;
        }
        
        pImpl->total_reads += batch->size();
        pImpl->input_queue.push(batch);
        
        if (pImpl->total_reads % 100000 == 0) {
            std::cout << "Queued " << pImpl->total_reads << " reads...\n";
        }
    }
    
    // Signal completion
    pImpl->input_queue.set_finished();
    
    // Wait for GPU workers
    for (auto& worker : pImpl->gpu_workers) {
        worker.join();
    }
    
    pImpl->output_queue.set_finished();
    result_processor.join();
    
    pImpl->end_time = std::chrono::steady_clock::now();
    
    std::cout << "Processed " << pImpl->total_reads << " total reads in " 
              << pImpl->total_batches << " batches\n";
}

std::vector<std::vector<Minimizer>> GPUMinimizerPipeline::process_batch(const ReadBatch& batch) {
    // Create a temporary extractor for single batch processing
    MinimizerExtractor extractor(pImpl->k_mer_size, pImpl->minimizer_window);
    return extractor.extract_minimizers(batch.sequences);
}

GPUMinimizerPipeline::Statistics GPUMinimizerPipeline::get_statistics() const {
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        pImpl->end_time - pImpl->start_time);
    
    return {
        pImpl->total_reads.load(),
        pImpl->total_batches.load(),
        static_cast<double>(duration.count())
    };
}

// MinimizerStats implementation
class MinimizerStats::Impl {
public:
    std::mutex mutex;
    size_t total_reads = 0;
    size_t total_minimizers = 0;
    std::unordered_map<uint64_t, uint32_t> minimizer_counts;
};

MinimizerStats::MinimizerStats() : pImpl(std::make_unique<Impl>()) {}

MinimizerStats::~MinimizerStats() = default;

void MinimizerStats::process_batch(const ReadBatch& batch, 
                                  const std::vector<std::vector<Minimizer>>& minimizers) {
    std::lock_guard<std::mutex> lock(pImpl->mutex);
    
    pImpl->total_reads += batch.size();
    
    for (const auto& read_minimizers : minimizers) {
        pImpl->total_minimizers += read_minimizers.size();
        
        for (const auto& m : read_minimizers) {
            pImpl->minimizer_counts[m.hash]++;
        }
    }
}

void MinimizerStats::print_summary() const {
    std::lock_guard<std::mutex> lock(pImpl->mutex);
    
    std::cout << "\n=== Minimizer Statistics ===\n";
    std::cout << "Total reads: " << pImpl->total_reads << "\n";
    std::cout << "Total minimizers: " << pImpl->total_minimizers << "\n";
    std::cout << "Unique minimizers: " << pImpl->minimizer_counts.size() << "\n";
    
    if (pImpl->total_reads > 0) {
        std::cout << "Average minimizers per read: " 
                  << static_cast<double>(pImpl->total_minimizers) / pImpl->total_reads << "\n";
    }
    
    // Skip top minimizers display - not needed for production use
    // This was causing hangs due to mutex/threading issues
}

size_t MinimizerStats::get_total_reads() const {
    std::lock_guard<std::mutex> lock(pImpl->mutex);
    return pImpl->total_reads;
}

size_t MinimizerStats::get_total_minimizers() const {
    std::lock_guard<std::mutex> lock(pImpl->mutex);
    return pImpl->total_minimizers;
}

size_t MinimizerStats::get_unique_minimizers() const {
    std::lock_guard<std::mutex> lock(pImpl->mutex);
    return pImpl->minimizer_counts.size();
}

double MinimizerStats::get_average_minimizers_per_read() const {
    std::lock_guard<std::mutex> lock(pImpl->mutex);
    if (pImpl->total_reads == 0) return 0.0;
    return static_cast<double>(pImpl->total_minimizers) / pImpl->total_reads;
}

std::vector<std::pair<uint64_t, uint32_t>> MinimizerStats::get_top_minimizers(size_t n) const {
    // Copy data while holding the lock
    std::vector<std::pair<uint64_t, uint32_t>> all_minimizers;
    {
        std::lock_guard<std::mutex> lock(pImpl->mutex);
        all_minimizers.reserve(pImpl->minimizer_counts.size());
        for (const auto& [hash, count] : pImpl->minimizer_counts) {
            all_minimizers.push_back({hash, count});
        }
    }
    // Release lock before sorting
    
    std::sort(all_minimizers.begin(), all_minimizers.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    if (all_minimizers.size() > n) {
        all_minimizers.resize(n);
    }
    
    return all_minimizers;
}

} // namespace biogpu