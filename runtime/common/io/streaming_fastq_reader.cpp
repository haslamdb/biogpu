// runtime/common/io/streaming_fastq_reader.cpp
#include "streaming_fastq_reader.h"
#include <iostream>
#include <algorithm>
#include <cstring>

namespace BioGPU {

bool StreamingFastqReader::open(const std::string& filename) {
    close();  // Close any existing files
    
    // Check if file is gzipped
    is_gzipped = (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz");
    
    if (is_gzipped) {
        gz_file = gzopen(filename.c_str(), "rb");
        if (!gz_file) {
            std::cerr << "Error: Cannot open gzipped file " << filename << std::endl;
            return false;
        }
        // Set buffer size for better performance
        gzbuffer(gz_file, 256 * 1024);  // 256KB buffer
    } else {
        file_stream.open(filename, std::ios::binary);
        if (!file_stream.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }
    }
    
    is_open = true;
    has_r2 = false;
    stop_reading = false;
    reading_complete = false;
    
    // Start the reader thread
    reader_thread = std::thread(&StreamingFastqReader::readerThreadFunc, this);
    
    return true;
}

bool StreamingFastqReader::openPaired(const std::string& r1_file, const std::string& r2_file) {
    close();  // Close any existing files
    
    // Open R1
    is_gzipped = (r1_file.size() > 3 && r1_file.substr(r1_file.size() - 3) == ".gz");
    
    if (is_gzipped) {
        gz_file = gzopen(r1_file.c_str(), "rb");
        gz_file_r2 = gzopen(r2_file.c_str(), "rb");
        
        if (!gz_file || !gz_file_r2) {
            std::cerr << "Error: Cannot open gzipped files " << r1_file << " or " << r2_file << std::endl;
            if (gz_file) gzclose(gz_file);
            if (gz_file_r2) gzclose(gz_file_r2);
            gz_file = nullptr;
            gz_file_r2 = nullptr;
            return false;
        }
        gzbuffer(gz_file, 256 * 1024);
        gzbuffer(gz_file_r2, 256 * 1024);
    } else {
        file_stream.open(r1_file, std::ios::binary);
        file_stream_r2.open(r2_file, std::ios::binary);
        
        if (!file_stream.is_open() || !file_stream_r2.is_open()) {
            std::cerr << "Error: Cannot open files " << r1_file << " or " << r2_file << std::endl;
            file_stream.close();
            file_stream_r2.close();
            return false;
        }
    }
    
    is_open = true;
    has_r2 = true;
    stop_reading = false;
    reading_complete = false;
    
    // Start the reader thread
    reader_thread = std::thread(&StreamingFastqReader::readerThreadFunc, this);
    
    return true;
}

void StreamingFastqReader::close() {
    // Signal reader thread to stop
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        stop_reading = true;
    }
    queue_cv.notify_all();
    
    // Wait for reader thread to finish
    if (reader_thread.joinable()) {
        reader_thread.join();
    }
    
    // Close files
    if (is_gzipped) {
        if (gz_file) {
            gzclose(gz_file);
            gz_file = nullptr;
        }
        if (gz_file_r2) {
            gzclose(gz_file_r2);
            gz_file_r2 = nullptr;
        }
    } else {
        if (file_stream.is_open()) {
            file_stream.close();
        }
        if (file_stream_r2.is_open()) {
            file_stream_r2.close();
        }
    }
    
    // Clear queue
    std::queue<std::shared_ptr<SequenceBatch>> empty;
    std::swap(batch_queue, empty);
    
    is_open = false;
    has_r2 = false;
}

std::shared_ptr<SequenceBatch> StreamingFastqReader::getNextBatch() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    
    // Wait until there's a batch or reading is complete
    queue_cv.wait(lock, [this] {
        return !batch_queue.empty() || reading_complete || stop_reading;
    });
    
    if (!batch_queue.empty()) {
        auto batch = batch_queue.front();
        batch_queue.pop();
        queue_cv.notify_all();  // Notify reader thread there's space
        return batch;
    }
    
    return nullptr;
}

void StreamingFastqReader::readerThreadFunc() {
    while (!stop_reading) {
        auto batch = std::make_shared<SequenceBatch>();
        batch->reserve(batch_size);
        
        bool success = readBatchFromFile(*batch);
        
        if (!success || batch->empty()) {
            // End of file or error
            std::lock_guard<std::mutex> lock(queue_mutex);
            reading_complete = true;
            queue_cv.notify_all();
            break;
        }
        
        // Update statistics
        total_reads += batch->size();
        total_bases += batch->getTotalBases();
        
        // Add batch to queue
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            
            // Wait if queue is full
            queue_cv.wait(lock, [this] {
                return batch_queue.size() < max_queue_size || stop_reading;
            });
            
            if (stop_reading) break;
            
            batch_queue.push(batch);
        }
        queue_cv.notify_all();
    }
}

bool StreamingFastqReader::readBatchFromFile(SequenceBatch& batch) {
    std::string line1, line2, line3, line4;  // FASTQ format: header, seq, +, qual
    std::string line1_r2, line2_r2, line3_r2, line4_r2;  // For R2
    
    batch.clear();
    
    while (batch.size() < batch_size) {
        // Read R1
        if (!readLine(line1) || !readLine(line2) || !readLine(line3) || !readLine(line4)) {
            break;  // End of file
        }
        
        // Validate FASTQ format
        if (line1.empty() || line1[0] != '@' || line3.empty() || line3[0] != '+') {
            std::cerr << "Warning: Invalid FASTQ format, skipping read" << std::endl;
            continue;
        }
        
        // If paired-end, read R2
        if (has_r2) {
            if (!readLineR2(line1_r2) || !readLineR2(line2_r2) || 
                !readLineR2(line3_r2) || !readLineR2(line4_r2)) {
                std::cerr << "Error: R2 file shorter than R1" << std::endl;
                return false;
            }
            
            if (line1_r2.empty() || line1_r2[0] != '@' || line3_r2.empty() || line3_r2[0] != '+') {
                std::cerr << "Warning: Invalid FASTQ format in R2, skipping read" << std::endl;
                continue;
            }
            
            // Verify read names match (ignoring /1 and /2 suffixes)
            std::string name1 = line1.substr(1);
            std::string name2 = line1_r2.substr(1);
            size_t space1 = name1.find(' ');
            size_t space2 = name2.find(' ');
            if (space1 != std::string::npos) name1 = name1.substr(0, space1);
            if (space2 != std::string::npos) name2 = name2.substr(0, space2);
            
            // Remove /1 or /2 suffix if present
            if (name1.size() > 2 && name1[name1.size()-2] == '/') {
                name1 = name1.substr(0, name1.size()-2);
            }
            if (name2.size() > 2 && name2[name2.size()-2] == '/') {
                name2 = name2.substr(0, name2.size()-2);
            }
            
            if (name1 != name2) {
                std::cerr << "Warning: Read names don't match: " << name1 << " vs " << name2 << std::endl;
            }
        }
        
        // Add to batch
        if (has_r2) {
            // For paired-end, concatenate R1 and R2 with a separator
            // Or store them separately based on your needs
            batch.addRead(line1, line2, line4);
            batch.addRead(line1_r2, line2_r2, line4_r2);
            batch.is_paired_end = true;
        } else {
            batch.addRead(line1, line2, line4);
            batch.is_paired_end = false;
        }
    }
    
    return !batch.empty();
}

bool StreamingFastqReader::readLine(std::string& line) {
    line.clear();
    
    if (is_gzipped) {
        if (!gz_file || gzeof(gz_file)) {
            return false;
        }
        
        char buffer[8192];
        if (gzgets(gz_file, buffer, sizeof(buffer)) == nullptr) {
            return false;
        }
        
        line = buffer;
        // Remove newline
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    } else {
        if (!std::getline(file_stream, line)) {
            return false;
        }
        // Remove carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    }
    
    return true;
}

bool StreamingFastqReader::readLineR2(std::string& line) {
    line.clear();
    
    if (is_gzipped) {
        if (!gz_file_r2 || gzeof(gz_file_r2)) {
            return false;
        }
        
        char buffer[8192];
        if (gzgets(gz_file_r2, buffer, sizeof(buffer)) == nullptr) {
            return false;
        }
        
        line = buffer;
        // Remove newline
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    } else {
        if (!std::getline(file_stream_r2, line)) {
            return false;
        }
        // Remove carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    }
    
    return true;
}

} // namespace BioGPU