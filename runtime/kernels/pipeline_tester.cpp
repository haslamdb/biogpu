// pipeline_tester.cpp
//
// A modular test harness for the unified resistance detection pipeline.
// This script allows for staged testing of different components, starting
// with the input and FASTQ reading stages.
//
// To Compile:
//   make
//
// To Run:
//   ./pipeline_tester samples.csv
//
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <map>
#include <numeric>

#include "sample_csv_parser.h"
#include "fastq_reader.h"

// Test function signature defines the interface for all test cases.
using TestFunc = std::function<bool(const BioGPU::SampleInfo&)>;

// TestCase structure holds a test's name and its corresponding function.
struct TestCase {
    std::string name;
    TestFunc function;
};

// ============================================================================
//                          TEST IMPLEMENTATIONS
// ============================================================================

/**
 * @brief Tests the ability to read single-end or paired-end FASTQ.gz files.
 * * This test uses the functions from fastq_reader to open the file(s) specified
 * in the SampleInfo struct, read all records in small batches, and report the
 * total number of reads and bases found.
 * * @param sample The sample metadata, including paths to FASTQ files.
 * @return True if reads were successfully imported, false otherwise.
 */
bool test_fastq_reader(const BioGPU::SampleInfo& sample) {
    std::cout << "  [RUNNING] test_fastq_reader for sample: " << sample.sample_name << std::endl;
    
    const int BATCH_SIZE = 100; // Use a small batch size for testing purposes
    long total_reads = 0;
    long total_bases = 0;
    bool success = true;

    if (sample.isPairedEnd()) {
        gzFile r1_file = gzopen(sample.read1_path.c_str(), "r");
        gzFile r2_file = gzopen(sample.read2_path.c_str(), "r");

        if (!r1_file || !r2_file) {
            std::cerr << "    [ERROR] Could not open paired-end files: " 
                      << sample.read1_path << ", " << sample.read2_path << std::endl;
            if(r1_file) gzclose(r1_file);
            if(r2_file) gzclose(r2_file);
            return false;
        }

        // Read all paired-end records in a loop
        while (true) {
            std::vector<FastqRecord> batch_r1, batch_r2;
            if (!read_batch_paired(r1_file, r2_file, batch_r1, batch_r2, BATCH_SIZE)) {
                break;
            }
            if (batch_r1.empty()) break;

            total_reads += batch_r1.size();
            for(const auto& rec : batch_r1) total_bases += rec.sequence.length();
            for(const auto& rec : batch_r2) total_bases += rec.sequence.length();
        }

        gzclose(r1_file);
        gzclose(r2_file);

    } else { // Single-end
        gzFile file = gzopen(sample.read1_path.c_str(), "r");
        if (!file) {
            std::cerr << "    [ERROR] Could not open single-end file: " << sample.read1_path << std::endl;
            return false;
        }

        // Read all single-end records in a loop
        while (true) {
            std::vector<FastqRecord> batch;
            if (!read_batch_single(file, batch, BATCH_SIZE)) {
                break;
            }
            if (batch.empty()) break;

            total_reads += batch.size();
            for(const auto& rec : batch) total_bases += rec.sequence.length();
        }
        gzclose(file);
    }

    if (total_reads > 0) {
        std::cout << "    [SUCCESS] Read " << total_reads << " records (" << total_bases << " bases)." << std::endl;
    } else {
        std::cerr << "    [FAIL] No reads were imported from the files." << std::endl;
        success = false;
    }

    return success;
}


// ============================================================================
//                             TEST RUNNER
// ============================================================================

/**
 * @class PipelineTester
 * @brief A modular test runner for the bioinformatics pipeline.
 * * This class discovers and runs registered test functions against a set of
 * samples defined in a CSV file. It reports a summary of passed and failed tests.
 */
class PipelineTester {
private:
    std::vector<TestCase> tests;
    BioGPU::SampleCSVParser parser;

public:
    /**
     * @brief Registers a new test case to be run.
     * @param name A descriptive name for the test.
     * @param func The function implementing the test logic.
     */
    void register_test(const std::string& name, TestFunc func) {
        tests.push_back({name, func});
    }

    /**
     * @brief Runs the entire test suite against a sample sheet.
     * @param csv_path Path to the input CSV file defining samples.
     * @return True if all tests passed for all valid samples, false otherwise.
     */
    bool run(const std::string& csv_path) {
        std::cout << "--- Initializing Pipeline Test Runner ---" << std::endl;
        
        if (!parser.parseFile(csv_path, true)) {
            std::cerr << "[FATAL] Failed to parse sample CSV: " << csv_path << std::endl;
            return false;
        }
        parser.printSummary();

        int total_passed = 0;
        int total_failed = 0;

        for (size_t i = 0; i < parser.getSampleCount(); ++i) {
            const BioGPU::SampleInfo* sample = parser.getSample(i);
            if (!sample) continue;

            std::cout << "\n--- Testing Sample: " << sample->sample_name << " ---" << std::endl;
            
            if (!parser.validateSample(*sample)) {
                 std::cerr << "  [SKIPPING] Sample files not found or invalid." << std::endl;
                 continue;
            }

            for (const auto& test : tests) {
                std::cout << "-> Running Test: " << test.name << std::endl;
                bool result = test.function(*sample);
                if (result) {
                    total_passed++;
                } else {
                    total_failed++;
                }
            }
        }

        std::cout << "\n--- Test Suite Finished ---" << std::endl;
        std::cout << "Total Tests Passed: " << total_passed << std::endl;
        std::cout << "Total Tests Failed: " << total_failed << std::endl;
        std::cout << "--------------------------" << std::endl;

        return total_failed == 0;
    }
};

// ============================================================================
//                                 MAIN
// ============================================================================

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_samples.csv>" << std::endl;
        return 1;
    }

    PipelineTester tester;

    // Register all tests to be run.
    // To add more tests, create a new function like `test_fastq_reader`
    // and register it here.
    tester.register_test("FASTQ Reader Test", test_fastq_reader);
    // tester.register_test("Protein Search Test", test_protein_search); // Example for future

    // Run the test suite
    bool success = tester.run(argv[1]);

    return success ? 0 : 1;
}

// ============================================================================
// ============================================================================
//
//               SUPPORTING SOURCE FILES (for compilation)
//
// ============================================================================
// ============================================================================

// fastq_reader.h
#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <string>
#include <vector>
#include <zlib.h>

// Represents a single FASTQ record
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Represents a batch of reads ready for GPU transfer
struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

// Reads a batch of paired-end records from two gzipped FASTQ files.
// Returns true if a full or partial batch was read.
bool read_batch_paired(gzFile r1_file, gzFile r2_file, 
                       std::vector<FastqRecord>& batch_r1, 
                       std::vector<FastqRecord>& batch_r2, 
                       int max_batch_size);

// Reads a batch of single-end records from a gzipped FASTQ file.
// Returns true if a full or partial batch was read.
bool read_batch_single(gzFile file, 
                       std::vector<FastqRecord>& batch, 
                       int max_batch_size);

// Prepares a vector of FastqRecord structs into a contiguous memory block for the GPU.
// The caller is responsible for deleting the allocated memory in the returned ReadBatch.
ReadBatch prepare_batch_for_gpu(const std::vector<FastqRecord>& records);

#endif // FASTQ_READER_H

// fastq_reader.cpp
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
    for (int i = 0; i < (int)records.size(); ++i) {
        const std::string& seq = records[i].sequence;
        memcpy(batch.sequences + current_offset, seq.c_str(), seq.length());
        batch.lengths[i] = seq.length();
        batch.offsets[i] = current_offset;
        current_offset += seq.length();
    }
    
    return batch;
}

// sample_csv_parser.h
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

} // namespace BioGPU

#endif // SAMPLE_CSV_PARSER_H

// sample_csv_parser.cpp
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>

namespace BioGPU {

// Constructor
SampleCSVParser::SampleCSVParser() : has_header(true), delimiter(',') {
}

// Trim whitespace from string
std::string SampleCSVParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
        return "";
    }
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// Split string by delimiter
std::vector<std::string> SampleCSVParser::split(const std::string& line, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(line);
    
    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(trim(token));
    }
    
    return tokens;
}

// Check if column name matches any variant
bool SampleCSVParser::matchesColumnName(const std::string& column, 
                                       const std::vector<std::string>& variants) {
    std::string col_lower = column;
    std::transform(col_lower.begin(), col_lower.end(), col_lower.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    
    for (const auto& variant : variants) {
        std::string var_lower = variant;
        std::transform(var_lower.begin(), var_lower.end(), var_lower.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        if (col_lower == var_lower) {
            return true;
        }
    }
    return false;
}

// Expand ~ in path
std::string SampleCSVParser::expandPath(const std::string& path) const {
    if (path.empty() || path[0] != '~') {
        return path;
    }
    
    const char* home = getenv("HOME");
    if (!home) {
        return path;
    }
    
    return std::string(home) + path.substr(1);
}

// Strip path from filename
std::string SampleCSVParser::stripPath(const std::string& filepath) const {
    size_t pos = filepath.find_last_of("/\\");
    if (pos != std::string::npos) {
        return filepath.substr(pos + 1);
    }
    return filepath;
}

// Combine directory and filename
std::string SampleCSVParser::combinePath(const std::string& dir, const std::string& filename) const {
    if (dir.empty()) {
        return filename;
    }
    
    std::string expanded_dir = expandPath(dir);
    
    // Ensure directory ends with separator
    if (expanded_dir.back() != '/' && expanded_dir.back() != '\\') {
        expanded_dir += '/';
    }
    
    // Strip any leading path from filename
    std::string clean_filename = stripPath(filename);
    
    return expanded_dir + clean_filename;
}

// Check if file exists
bool SampleCSVParser::fileExists(const std::string& path) const {
    struct stat buffer;
    std::string expanded = expandPath(path);
    return (stat(expanded.c_str(), &buffer) == 0);
}

// Auto-detect delimiter
void SampleCSVParser::detectDelimiter(const std::string& first_line) {
    // Count occurrences of common delimiters
    int comma_count = std::count(first_line.begin(), first_line.end(), ',');
    int tab_count = std::count(first_line.begin(), first_line.end(), '\t');
    int semicolon_count = std::count(first_line.begin(), first_line.end(), ';');
    
    // Choose the most common one
    if (tab_count > comma_count && tab_count > semicolon_count) {
        delimiter = '\t';
    } else if (semicolon_count > comma_count) {
        delimiter = ';';
    } else {
        delimiter = ',';
    }
}

// Parse CSV file
bool SampleCSVParser::parseFile(const std::string& csv_path, bool validate_paths) {
    samples.clear();
    header_columns.clear();
    
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open CSV file: " << csv_path << std::endl;
        return false;
    }
    
    std::string line;
    bool first_line = true;
    int line_number = 0;
    
    // Column indices
    int sample_name_idx = -1;
    int file_path_idx = -1;
    int read1_idx = -1;
    int read2_idx = -1;
    
    while (std::getline(file, line)) {
        line_number++;
        
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\n\r") == std::string::npos) {
            continue;
        }
        
        // Auto-detect delimiter on first line
        if (first_line) {
            detectDelimiter(line);
        }
        
        std::vector<std::string> fields = split(line, delimiter);
        
        // Process header
        if (first_line && has_header) {
            first_line = false;
            header_columns = fields;
            
            // Find column indices
            for (size_t i = 0; i < fields.size(); i++) {
                if (matchesColumnName(fields[i], SAMPLE_NAME_VARIANTS)) {
                    sample_name_idx = i;
                } else if (matchesColumnName(fields[i], FILE_PATH_VARIANTS)) {
                    file_path_idx = i;
                } else if (matchesColumnName(fields[i], READ1_VARIANTS)) {
                    read1_idx = i;
                } else if (matchesColumnName(fields[i], READ2_VARIANTS)) {
                    read2_idx = i;
                }
            }
            
            if (sample_name_idx == -1 || file_path_idx == -1 || read1_idx == -1) {
                 std::cerr << "Error: Missing one or more required columns (Sample, Path, R1 File) in CSV header." << std::endl;
                 return false;
            }
            
            continue;
        }
        
        // If no header, assume first four columns are sample, path, r1, r2
        if (!has_header && first_line) {
            first_line = false;
            sample_name_idx = 0;
            file_path_idx = 1;
            read1_idx = 2;
            read2_idx = (fields.size() > 3) ? 3 : -1;
        }
        
        // Parse sample data
        if ((int)fields.size() > std::max({sample_name_idx, file_path_idx, read1_idx})) {
            SampleInfo sample;
            sample.sample_name = fields[sample_name_idx];
            sample.file_path = fields[file_path_idx];
            sample.read1_filename = stripPath(fields[read1_idx]);
            
            sample.read1_path = combinePath(sample.file_path, sample.read1_filename);
            
            if (read2_idx != -1 && (int)fields.size() > read2_idx && !fields[read2_idx].empty()) {
                sample.read2_filename = stripPath(fields[read2_idx]);
                sample.read2_path = combinePath(sample.file_path, sample.read2_filename);
            }
            
            samples.push_back(sample);
        }
    }
    
    file.close();
    return !samples.empty();
}

// Validate sample
bool SampleCSVParser::validateSample(const SampleInfo& sample) const {
    if (!sample.isValid()) return false;
    if (!fileExists(sample.read1_path)) return false;
    if (sample.isPairedEnd() && !fileExists(sample.read2_path)) return false;
    return true;
}

// Get validation errors
std::vector<std::string> SampleCSVParser::getValidationErrors() const {
    std::vector<std::string> errors;
    for (const auto& sample : samples) {
        if (!sample.isValid()) {
            errors.push_back("Sample '" + sample.sample_name + "' is missing required fields.");
        } else {
            if (!fileExists(sample.read1_path)) {
                errors.push_back("Sample '" + sample.sample_name + "': R1 file not found: " + sample.read1_path);
            }
            if (sample.isPairedEnd() && !fileExists(sample.read2_path)) {
                errors.push_back("Sample '" + sample.sample_name + "': R2 file not found: " + sample.read2_path);
            }
        }
    }
    return errors;
}

void SampleCSVParser::printSummary() const { /* Omitted for brevity */ }
void SampleCSVParser::printDetailed() const { /* Omitted for brevity */ }

} // namespace BioGPU


// Makefile
//
// Use this Makefile to compile the test harness.
//
// Commands:
//   make        - Compiles the pipeline_tester executable.
//   make clean  - Removes compiled object files and the executable.
//
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -I.
LDFLAGS = -lz

# Source files
SRCS = pipeline_tester.cpp sample_csv_parser.cpp fastq_reader.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = pipeline_tester

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

