#ifndef DETECTOR_INTERFACE_H
#define DETECTOR_INTERFACE_H

#include <string>
#include <vector>
#include <memory>
#include <json/json.h>

// Forward declarations
struct ReadBatch {
    std::vector<std::string> read_ids;
    std::vector<std::string> sequences_r1;
    std::vector<std::string> sequences_r2;
    std::vector<std::string> qualities_r1;
    std::vector<std::string> qualities_r2;
    size_t num_reads;
    
    ReadBatch() : num_reads(0) {}
    
    void clear() {
        read_ids.clear();
        sequences_r1.clear();
        sequences_r2.clear();
        qualities_r1.clear();
        qualities_r2.clear();
        num_reads = 0;
    }
};

// Base configuration for all detectors
struct DetectorConfig {
    int gpu_device;
    int batch_size;
    bool use_bloom_filter;
    double min_identity;
    double min_coverage;
    int threads;
    std::string output_dir;
    std::string sample_id;
    
    DetectorConfig() : 
        gpu_device(0), 
        batch_size(50000), 
        use_bloom_filter(false),
        min_identity(0.90),
        min_coverage(0.80),
        threads(8) {}
};

// Progress callback function type
typedef std::function<void(const std::string& stage, int percentage, const std::string& message)> ProgressCallback;

// Abstract base class for detection modules
class DetectorInterface {
public:
    virtual ~DetectorInterface() = default;
    
    // Initialize the detector with configuration
    virtual bool initialize(const DetectorConfig& config) = 0;
    
    // Process a batch of reads
    virtual void processBatch(const ReadBatch& batch) = 0;
    
    // Finalize processing and generate results
    virtual void finalize() = 0;
    
    // Get results in JSON format
    virtual void getResults(Json::Value& results) = 0;
    
    // Generate output files
    virtual void writeOutputFiles(const std::string& output_dir) = 0;
    
    // Set progress callback
    virtual void setProgressCallback(ProgressCallback callback) { 
        progress_callback = callback; 
    }
    
    // Get detector name for logging
    virtual std::string getName() const = 0;
    
protected:
    ProgressCallback progress_callback;
    DetectorConfig config;
    
    // Helper to report progress
    void reportProgress(const std::string& stage, int percentage, const std::string& message) {
        if (progress_callback) {
            progress_callback(stage, percentage, message);
        }
    }
};

// Factory function type for creating detectors
typedef std::unique_ptr<DetectorInterface> (*DetectorFactory)();

#endif // DETECTOR_INTERFACE_H