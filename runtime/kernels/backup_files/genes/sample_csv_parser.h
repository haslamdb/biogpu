// sample_csv_parser.h
// CSV parser for batch sample processing in BioGPU pipelines
// Parses CSV files containing sample names and file paths

#ifndef SAMPLE_CSV_PARSER_H
#define SAMPLE_CSV_PARSER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <functional>

namespace BioGPU {

// Structure to hold information about a single sample
struct SampleInfo {
    std::string sample_name;
    std::string file_path;      // Base directory path
    std::string read1_filename; // R1 filename (without path)
    std::string read2_filename; // R2 filename (without path)
    std::string read1_path;     // Full path (file_path + read1_filename)
    std::string read2_path;     // Full path (file_path + read2_filename)
    std::map<std::string, std::string> metadata;  // Additional optional columns
    
    bool isPairedEnd() const {
        return !read2_path.empty() && !read2_filename.empty();
    }
    
    bool isValid() const {
        return !sample_name.empty() && !read1_path.empty() && !file_path.empty();
    }
};

// Main CSV parser class
class SampleCSVParser {
private:
    std::vector<SampleInfo> samples;
    std::vector<std::string> header_columns;
    bool has_header;
    char delimiter;
    
    // Recognized column name variations
    const std::vector<std::string> SAMPLE_NAME_VARIANTS = {
        "Sample Name", "SampleName", "Sample ID", "SampleID", 
        "sample_name", "samplename", "sample_id", "sampleid",
        "Sample", "sample", "ID", "id", "Name", "name"
    };
    
    const std::vector<std::string> FILE_PATH_VARIANTS = {
        "FilePath", "File Path", "filepath", "file_path",
        "Path", "path", "Directory", "directory", "Dir", "dir",
        "BaseDir", "basedir", "Base Dir", "base_dir"
    };
    
    const std::vector<std::string> READ1_VARIANTS = {
        "R1 file", "R1file", "R1_file", "R1", "r1",
        "Read1", "read1", "Forward", "forward",
        "File1", "file1", "FASTQ1", "fastq1"
    };
    
    const std::vector<std::string> READ2_VARIANTS = {
        "R2 file", "R2file", "R2_file", "R2", "r2",
        "Read2", "read2", "Reverse", "reverse",
        "File2", "file2", "FASTQ2", "fastq2"
    };
    
    // Helper functions
    std::string trim(const std::string& str);
    std::vector<std::string> split(const std::string& line, char delim);
    bool matchesColumnName(const std::string& column, const std::vector<std::string>& variants);
    std::string expandPath(const std::string& path) const;
    std::string stripPath(const std::string& filepath) const;
    std::string combinePath(const std::string& dir, const std::string& filename) const;
    bool fileExists(const std::string& path) const;
    void detectDelimiter(const std::string& first_line);
    
public:
    SampleCSVParser();
    
    // Parse CSV file
    bool parseFile(const std::string& csv_path, bool validate_paths = true);
    
    // Get parsed samples
    const std::vector<SampleInfo>& getSamples() const { return samples; }
    size_t getSampleCount() const { return samples.size(); }
    
    // Get specific sample
    const SampleInfo* getSample(size_t index) const {
        if (index < samples.size()) {
            return &samples[index];
        }
        return nullptr;
    }
    
    // Find sample by name
    const SampleInfo* findSampleByName(const std::string& name) const;
    
    // Validation
    bool validateSample(const SampleInfo& sample) const;
    std::vector<std::string> getValidationErrors() const;
    
    // Utility functions
    void printSummary() const;
    void printDetailed() const;
    
    // Configuration
    void setDelimiter(char delim) { delimiter = delim; }
    void setHasHeader(bool has_hdr) { has_header = has_hdr; }
};

// Utility function for batch processing
class BatchProcessor {
private:
    SampleCSVParser parser;
    std::string output_base_dir;
    bool create_sample_dirs;
    
public:
    BatchProcessor(const std::string& output_dir = "results", bool create_dirs = true);
    
    // Load samples from CSV
    bool loadSamples(const std::string& csv_path);
    
    // Get output path for a sample
    std::string getOutputPath(const SampleInfo& sample, const std::string& suffix = "") const;
    
    // Process function callback type
    typedef std::function<int(const SampleInfo&, const std::string&)> ProcessFunc;
    
    // Run batch processing
    int processBatch(ProcessFunc process_func, bool stop_on_error = false);
    
    // Get parser for direct access
    SampleCSVParser& getParser() { return parser; }
};

} // namespace BioGPU

#endif // SAMPLE_CSV_PARSER_H