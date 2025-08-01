#ifndef UNIFIED_BATCH_PROCESSOR_H
#define UNIFIED_BATCH_PROCESSOR_H

#include "detector_interface.h"
#include "../../common/io/streaming_fastq_reader.h"
#include "../../common/io/sequence_batch.h"
#include <vector>
#include <memory>
#include <functional>

namespace BioGPU {

// Convert SequenceBatch to ReadBatch for detector interface
class BatchConverter {
public:
    static ReadBatch convertToReadBatch(const SequenceBatch& seq_batch) {
        ReadBatch read_batch;
        
        if (seq_batch.is_paired_end) {
            // Paired-end: StreamingFastqReader stores R1 and R2 alternately
            size_t num_pairs = seq_batch.size() / 2;
            read_batch.num_reads = num_pairs;
            
            read_batch.read_ids.reserve(num_pairs);
            read_batch.sequences_r1.reserve(num_pairs);
            read_batch.sequences_r2.reserve(num_pairs);
            read_batch.qualities_r1.reserve(num_pairs);
            read_batch.qualities_r2.reserve(num_pairs);
            
            for (size_t i = 0; i < num_pairs; ++i) {
                size_t r1_idx = i * 2;
                size_t r2_idx = i * 2 + 1;
                
                // Extract read ID from R1 header (before first space)
                std::string header = seq_batch.headers[r1_idx];
                size_t space_pos = header.find(' ');
                std::string read_id = (space_pos != std::string::npos) ? 
                                      header.substr(0, space_pos) : header;
                
                // Remove @ prefix if present
                if (!read_id.empty() && read_id[0] == '@') {
                    read_id = read_id.substr(1);
                }
                
                read_batch.read_ids.push_back(read_id);
                read_batch.sequences_r1.push_back(seq_batch.sequences[r1_idx]);
                read_batch.sequences_r2.push_back(seq_batch.sequences[r2_idx]);
                read_batch.qualities_r1.push_back(seq_batch.qualities[r1_idx]);
                read_batch.qualities_r2.push_back(seq_batch.qualities[r2_idx]);
            }
        } else {
            // Single-end reads
            size_t num_reads = seq_batch.size();
            read_batch.num_reads = num_reads;
            
            read_batch.read_ids.reserve(num_reads);
            read_batch.sequences_r1.reserve(num_reads);
            read_batch.qualities_r1.reserve(num_reads);
            
            for (size_t i = 0; i < num_reads; ++i) {
                std::string header = seq_batch.headers[i];
                size_t space_pos = header.find(' ');
                std::string read_id = (space_pos != std::string::npos) ? 
                                      header.substr(0, space_pos) : header;
                
                // Remove @ prefix if present
                if (!read_id.empty() && read_id[0] == '@') {
                    read_id = read_id.substr(1);
                }
                
                read_batch.read_ids.push_back(read_id);
                read_batch.sequences_r1.push_back(seq_batch.sequences[i]);
                read_batch.qualities_r1.push_back(seq_batch.qualities[i]);
            }
        }
        
        return read_batch;
    }
    
    static ReadBatch convertToReadBatchSingle(const SequenceBatch& seq_batch) {
        ReadBatch read_batch;
        
        size_t num_reads = seq_batch.size();
        read_batch.num_reads = num_reads;
        
        read_batch.read_ids.reserve(num_reads);
        read_batch.sequences_r1.reserve(num_reads);
        read_batch.qualities_r1.reserve(num_reads);
        
        for (size_t i = 0; i < num_reads; ++i) {
            std::string header = seq_batch.headers[i];
            size_t space_pos = header.find(' ');
            std::string read_id = (space_pos != std::string::npos) ? 
                                  header.substr(0, space_pos) : header;
            
            read_batch.read_ids.push_back(read_id);
            read_batch.sequences_r1.push_back(seq_batch.sequences[i]);
            read_batch.qualities_r1.push_back(seq_batch.qualities[i]);
        }
        
        return read_batch;
    }
};

// Unified batch processor that feeds multiple detectors
class UnifiedBatchProcessor {
public:
    UnifiedBatchProcessor(size_t batch_size = 50000, size_t queue_depth = 10) 
        : batch_size_(batch_size), queue_depth_(queue_depth) {}
    
    // Add a detector to the pipeline
    void addDetector(std::shared_ptr<DetectorInterface> detector) {
        detectors_.push_back(detector);
    }
    
    // Process paired-end FASTQ files
    bool processPairedEnd(const std::string& r1_path, const std::string& r2_path,
                          ProgressCallback progress_callback = nullptr) {
        
        // Create a single streaming reader for paired-end files
        StreamingFastqReader reader(batch_size_, queue_depth_);
        
        if (!reader.openPaired(r1_path, r2_path)) {
            std::cerr << "Failed to open paired-end files: " << r1_path << " and " << r2_path << std::endl;
            return false;
        }
        
        // Initialize all detectors
        for (auto& detector : detectors_) {
            if (!detector->initialize(detector->config)) {
                std::cerr << "Failed to initialize detector: " << detector->getName() << std::endl;
                return false;
            }
        }
        
        // Process batches asynchronously
        int batch_num = 0;
        size_t total_reads = 0;
        
        while (reader.hasNext()) {
            auto batch = reader.getNextBatch();
            
            if (!batch) break;
            
            // Convert to ReadBatch format
            ReadBatch read_batch = BatchConverter::convertToReadBatch(*batch);
            total_reads += read_batch.num_reads;
            
            // Feed batch to all detectors
            for (auto& detector : detectors_) {
                detector->processBatch(read_batch);
            }
            
            batch_num++;
            
            // Report progress
            if (progress_callback) {
                int progress = std::min(90, batch_num * 10);  // Estimate progress
                progress_callback("batch_processing", progress, 
                    "Processed " + std::to_string(total_reads) + " reads in " + 
                    std::to_string(batch_num) + " batches");
            }
        }
        
        // Finalize all detectors
        for (auto& detector : detectors_) {
            detector->finalize();
        }
        
        // Get statistics
        total_reads_processed_ = reader.getTotalReads();
        total_bases_processed_ = reader.getTotalBases();
        
        return true;
    }
    
    // Process single-end FASTQ file
    bool processSingleEnd(const std::string& fastq_path,
                          ProgressCallback progress_callback = nullptr) {
        
        StreamingFastqReader reader(batch_size_, queue_depth_);
        
        if (!reader.open(fastq_path)) {
            std::cerr << "Failed to open FASTQ file: " << fastq_path << std::endl;
            return false;
        }
        
        // Initialize all detectors
        for (auto& detector : detectors_) {
            if (!detector->initialize(detector->config)) {
                std::cerr << "Failed to initialize detector: " << detector->getName() << std::endl;
                return false;
            }
        }
        
        // Process batches
        int batch_num = 0;
        size_t total_reads = 0;
        
        while (reader.hasNext()) {
            auto batch = reader.getNextBatch();
            if (!batch) break;
            
            ReadBatch read_batch = BatchConverter::convertToReadBatchSingle(*batch);
            total_reads += read_batch.num_reads;
            
            // Feed batch to all detectors
            for (auto& detector : detectors_) {
                detector->processBatch(read_batch);
            }
            
            batch_num++;
            
            if (progress_callback) {
                int progress = std::min(90, batch_num * 10);
                progress_callback("batch_processing", progress, 
                    "Processed " + std::to_string(total_reads) + " reads");
            }
        }
        
        // Finalize all detectors
        for (auto& detector : detectors_) {
            detector->finalize();
        }
        
        total_reads_processed_ = reader.getTotalReads();
        total_bases_processed_ = reader.getTotalBases();
        
        return true;
    }
    
    // Get processing statistics
    size_t getTotalReadsProcessed() const { return total_reads_processed_; }
    size_t getTotalBasesProcessed() const { return total_bases_processed_; }
    
    // Get results from all detectors
    void getAllResults(Json::Value& combined_results) {
        combined_results = Json::Value(Json::objectValue);
        
        for (auto& detector : detectors_) {
            Json::Value detector_results;
            detector->getResults(detector_results);
            combined_results[detector->getName()] = detector_results;
        }
        
        // Add processing statistics
        combined_results["statistics"]["total_reads"] = (Json::Value::UInt64)total_reads_processed_;
        combined_results["statistics"]["total_bases"] = (Json::Value::UInt64)total_bases_processed_;
    }
    
    // Write output files for all detectors
    void writeAllOutputs(const std::string& output_dir) {
        for (auto& detector : detectors_) {
            detector->writeOutputFiles(output_dir);
        }
    }
    
private:
    std::vector<std::shared_ptr<DetectorInterface>> detectors_;
    size_t batch_size_;
    size_t queue_depth_;
    size_t total_reads_processed_ = 0;
    size_t total_bases_processed_ = 0;
};

} // namespace BioGPU

#endif // UNIFIED_BATCH_PROCESSOR_H