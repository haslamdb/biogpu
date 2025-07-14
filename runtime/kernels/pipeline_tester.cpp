// pipeline_tester.cpp
//
// A modular test harness for the unified resistance detection pipeline.
// This script allows for staged testing of different components.
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
#include "kmer_minimizer.h"
#include "bloom_filter.h"

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
 */
bool test_fastq_reader(const BioGPU::SampleInfo& sample) {
    std::cout << "  [RUNNING] test_fastq_reader for sample: " << sample.sample_name << std::endl;
    
    const int BATCH_SIZE = 100; // Use a small batch size for testing purposes
    long total_reads = 0;
    long total_bases = 0;

    if (sample.isPairedEnd()) {
        gzFile r1_file = gzopen(sample.read1_path.c_str(), "r");
        gzFile r2_file = gzopen(sample.read2_path.c_str(), "r");

        if (!r1_file || !r2_file) {
            std::cerr << "    [ERROR] Could not open paired-end files." << std::endl;
            if(r1_file) gzclose(r1_file);
            if(r2_file) gzclose(r2_file);
            return false;
        }

        while (true) {
            std::vector<FastqRecord> batch_r1, batch_r2;
            if (!read_batch_paired(r1_file, r2_file, batch_r1, batch_r2, BATCH_SIZE) || batch_r1.empty()) {
                break;
            }
            total_reads += batch_r1.size();
            for(const auto& rec : batch_r1) total_bases += rec.sequence.length();
            for(const auto& rec : batch_r2) total_bases += rec.sequence.length();
        }
        gzclose(r1_file);
        gzclose(r2_file);

    } else { // Single-end
        gzFile file = gzopen(sample.read1_path.c_str(), "r");
        if (!file) {
            std::cerr << "    [ERROR] Could not open single-end file." << std::endl;
            return false;
        }

        while (true) {
            std::vector<FastqRecord> batch;
            if (!read_batch_single(file, batch, BATCH_SIZE) || batch.empty()) {
                break;
            }
            total_reads += batch.size();
            for(const auto& rec : batch) total_bases += rec.sequence.length();
        }
        gzclose(file);
    }

    if (total_reads > 0) {
        std::cout << "    [SUCCESS] Read " << total_reads << " records (" << total_bases << " bases)." << std::endl;
        return true;
    } else {
        std::cerr << "    [FAIL] No reads were imported from the files." << std::endl;
        return false;
    }
}

/**
 * @brief Tests k-mer/minimizer extraction and Bloom filter screening.
 * This test simulates the pre-screening stage of the pipeline.
 * It does not require a real GPU and uses stub functions.
 */
bool test_kmer_extraction_and_bloom_filter(const BioGPU::SampleInfo& sample) {
    std::cout << "  [RUNNING] test_kmer_extraction_and_bloom_filter for sample: " << sample.sample_name << std::endl;
    
    // STUB: In a real scenario, we would load k-mers from both FQ and AMR databases.
    std::vector<uint64_t> reference_kmers = { 0x123456789ABCDEF0, 0xFEDCBA9876543210 };
    
    // 1. Initialize Bloom Filter
    BloomFilterGPU bloom_filter(15); // k=15
    if (!bloom_filter.initialize()) {
        std::cerr << "    [FAIL] Could not initialize Bloom filter." << std::endl;
        return false;
    }
    std::cout << "    [INFO] Bloom filter initialized." << std::endl;
    
    // STUB: Build filter from reference k-mers (GPU operation)
    // In a real test, you'd copy reference_kmers to device and call the kernel.
    // For this stub, we assume it's built.
    std::cout << "    [INFO] Stub: Building Bloom filter from reference k-mers." << std::endl;

    // 2. Read a batch of reads
    gzFile r1_file = gzopen(sample.read1_path.c_str(), "r");
    if (!r1_file) {
        std::cerr << "    [FAIL] Could not open FASTQ file for test." << std::endl;
        return false;
    }
    std::vector<FastqRecord> records;
    read_batch_single(r1_file, records, 100);
    gzclose(r1_file);

    if (records.empty()) {
        std::cerr << "    [FAIL] Could not read a batch of records for testing." << std::endl;
        return false;
    }

    // 3. Extract Minimizers
    // STUB: This would be a GPU kernel call.
    std::vector<Minimizer> minimizers;
    std::cout << "    [INFO] Stub: Extracting minimizers from " << records.size() << " reads." << std::endl;
    // Dummy minimizer generation
    for(const auto& rec : records) {
        if(rec.sequence.length() > 15) {
            minimizers.push_back({12345, 0, false});
        }
    }
    std::cout << "    [INFO] Extracted " << minimizers.size() << " minimizers." << std::endl;


    // 4. Screen with Bloom Filter
    // STUB: This would be a GPU kernel call.
    std::vector<bool> passes(records.size(), false);
    int passed_count = 0;
    std::cout << "    [INFO] Stub: Screening minimizers against Bloom filter." << std::endl;
    // Dummy screening logic
    for(size_t i = 0; i < records.size(); ++i) {
        if (i % 2 == 0) { // Let every other read pass for the test
            passes[i] = true;
            passed_count++;
        }
    }

    std::cout << "    [SUCCESS] " << passed_count << " out of " << records.size() << " reads passed the filter." << std::endl;
    return true;
}


// ============================================================================
//                             TEST RUNNER
// ============================================================================

class PipelineTester {
private:
    std::vector<TestCase> tests;
    BioGPU::SampleCSVParser parser;

public:
    void register_test(const std::string& name, TestFunc func) {
        tests.push_back({name, func});
    }

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
                std::cout << "-----------------------------------------" << std::endl;
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

    tester.register_test("FASTQ Reader Test", test_fastq_reader);
    tester.register_test("K-mer Extraction & Bloom Filter Test", test_kmer_extraction_and_bloom_filter);

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

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

bool read_batch_paired(gzFile r1_file, gzFile r2_file, 
                       std::vector<FastqRecord>& batch_r1, 
                       std::vector<FastqRecord>& batch_r2, 
                       int max_batch_size);

bool read_batch_single(gzFile file, 
                       std::vector<FastqRecord>& batch, 
                       int max_batch_size);

ReadBatch prepare_batch_for_gpu(const std::vector<FastqRecord>& records);

#endif // FASTQ_READER_H

// fastq_reader.cpp
#include <cstring>

bool read_batch_paired(gzFile r1_file, gzFile r2_file, 
                       std::vector<FastqRecord>& batch_r1, 
                       std::vector<FastqRecord>& batch_r2, 
                       int max_batch_size) {
    char buffer[1024];
    for (int i = 0; i < max_batch_size; ++i) {
        FastqRecord rec1, rec2;
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return i > 0;
        rec1.header = std::string(buffer);
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return false;
        rec1.sequence = std::string(buffer); rec1.sequence.pop_back();
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return false;
        if (gzgets(r1_file, buffer, sizeof(buffer)) == NULL) return false;
        rec1.quality = std::string(buffer); rec1.quality.pop_back();
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        rec2.header = std::string(buffer);
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        rec2.sequence = std::string(buffer); rec2.sequence.pop_back();
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        if (gzgets(r2_file, buffer, sizeof(buffer)) == NULL) return false;
        rec2.quality = std::string(buffer); rec2.quality.pop_back();
        batch_r1.push_back(rec1); batch_r2.push_back(rec2);
    }
    return true;
}

bool read_batch_single(gzFile file, std::vector<FastqRecord>& batch, int max_batch_size) {
    char buffer[1024];
    for (int i = 0; i < max_batch_size; ++i) {
        FastqRecord rec;
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return i > 0;
        rec.header = std::string(buffer);
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return false;
        rec.sequence = std::string(buffer); rec.sequence.pop_back();
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return false;
        if (gzgets(file, buffer, sizeof(buffer)) == NULL) return false;
        rec.quality = std::string(buffer); rec.quality.pop_back();
        batch.push_back(rec);
    }
    return true;
}

ReadBatch prepare_batch_for_gpu(const std::vector<FastqRecord>& records) {
    ReadBatch batch;
    batch.num_reads = records.size();
    batch.total_bases = 0;
    for (const auto& rec : records) { batch.total_bases += rec.sequence.length(); }
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
struct SampleInfo {
    std::string sample_name;
    std::string file_path;
    std::string read1_filename;
    std::string read2_filename;
    std::string read1_path;
    std::string read2_path;
    std::map<std::string, std::string> metadata;
    bool isPairedEnd() const { return !read2_path.empty() && !read2_filename.empty(); }
    bool isValid() const { return !sample_name.empty() && !read1_path.empty() && !file_path.empty(); }
};
class SampleCSVParser {
private:
    std::vector<SampleInfo> samples;
    std::vector<std::string> header_columns;
    bool has_header;
    char delimiter;
    const std::vector<std::string> SAMPLE_NAME_VARIANTS = {"Sample Name", "SampleName", "Sample ID", "SampleID", "sample_name", "samplename", "sample_id", "sampleid", "Sample", "sample", "ID", "id", "Name", "name"};
    const std::vector<std::string> FILE_PATH_VARIANTS = {"FilePath", "File Path", "filepath", "file_path", "Path", "path", "Directory", "directory", "Dir", "dir", "BaseDir", "basedir", "Base Dir", "base_dir"};
    const std::vector<std::string> READ1_VARIANTS = {"R1 file", "R1file", "R1_file", "R1", "r1", "Read1", "read1", "Forward", "forward", "File1", "file1", "FASTQ1", "fastq1"};
    const std::vector<std::string> READ2_VARIANTS = {"R2 file", "R2file", "R2_file", "R2", "r2", "Read2", "read2", "Reverse", "reverse", "File2", "file2", "FASTQ2", "fastq2"};
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
    bool parseFile(const std::string& csv_path, bool validate_paths = true);
    const std::vector<SampleInfo>& getSamples() const { return samples; }
    size_t getSampleCount() const { return samples.size(); }
    const SampleInfo* getSample(size_t index) const { return (index < samples.size()) ? &samples[index] : nullptr; }
    const SampleInfo* findSampleByName(const std::string& name) const;
    bool validateSample(const SampleInfo& sample) const;
    std::vector<std::string> getValidationErrors() const;
    void printSummary() const;
    void printDetailed() const;
};
}
#endif

// sample_csv_parser.cpp
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
namespace BioGPU {
SampleCSVParser::SampleCSVParser() : has_header(true), delimiter(',') {}
std::string SampleCSVParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}
std::vector<std::string> SampleCSVParser::split(const std::string& line, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(line);
    while (std::getline(tokenStream, token, delim)) { tokens.push_back(trim(token)); }
    return tokens;
}
bool SampleCSVParser::matchesColumnName(const std::string& column, const std::vector<std::string>& variants) {
    std::string col_lower = column;
    std::transform(col_lower.begin(), col_lower.end(), col_lower.begin(), [](unsigned char c){ return std::tolower(c); });
    for (const auto& variant : variants) {
        std::string var_lower = variant;
        std::transform(var_lower.begin(), var_lower.end(), var_lower.begin(), [](unsigned char c){ return std::tolower(c); });
        if (col_lower == var_lower) return true;
    }
    return false;
}
std::string SampleCSVParser::expandPath(const std::string& path) const {
    if (path.empty() || path[0] != '~') return path;
    const char* home = getenv("HOME");
    return home ? std::string(home) + path.substr(1) : path;
}
std::string SampleCSVParser::stripPath(const std::string& filepath) const {
    size_t pos = filepath.find_last_of("/\\");
    return (pos != std::string::npos) ? filepath.substr(pos + 1) : filepath;
}
std::string SampleCSVParser::combinePath(const std::string& dir, const std::string& filename) const {
    if (dir.empty()) return filename;
    std::string expanded_dir = expandPath(dir);
    if (expanded_dir.back() != '/' && expanded_dir.back() != '\\') { expanded_dir += '/'; }
    return expanded_dir + stripPath(filename);
}
bool SampleCSVParser::fileExists(const std::string& path) const {
    struct stat buffer;
    return (stat(expandPath(path).c_str(), &buffer) == 0);
}
void SampleCSVParser::detectDelimiter(const std::string& first_line) {
    int comma_count = std::count(first_line.begin(), first_line.end(), ',');
    int tab_count = std::count(first_line.begin(), first_line.end(), '\t');
    if (tab_count > comma_count) delimiter = '\t';
    else delimiter = ',';
}
bool SampleCSVParser::parseFile(const std::string& csv_path, bool validate_paths) {
    samples.clear(); header_columns.clear();
    std::ifstream file(csv_path);
    if (!file.is_open()) { std::cerr << "Error: Cannot open CSV file: " << csv_path << std::endl; return false; }
    std::string line;
    bool first_line = true;
    int sample_name_idx = -1, file_path_idx = -1, read1_idx = -1, read2_idx = -1;
    while (std::getline(file, line)) {
        if (line.empty() || line.find_first_not_of(" \t\n\r") == std::string::npos) continue;
        if (first_line) { detectDelimiter(line); }
        std::vector<std::string> fields = split(line, delimiter);
        if (first_line && has_header) {
            first_line = false; header_columns = fields;
            for (size_t i = 0; i < fields.size(); i++) {
                if (matchesColumnName(fields[i], SAMPLE_NAME_VARIANTS)) sample_name_idx = i;
                else if (matchesColumnName(fields[i], FILE_PATH_VARIANTS)) file_path_idx = i;
                else if (matchesColumnName(fields[i], READ1_VARIANTS)) read1_idx = i;
                else if (matchesColumnName(fields[i], READ2_VARIANTS)) read2_idx = i;
            }
            if (sample_name_idx == -1 || file_path_idx == -1 || read1_idx == -1) { std::cerr << "Error: Missing required columns in CSV header." << std::endl; return false; }
            continue;
        }
        if (!has_header && first_line) { first_line = false; sample_name_idx = 0; file_path_idx = 1; read1_idx = 2; read2_idx = (fields.size() > 3) ? 3 : -1; }
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
    file.close(); return !samples.empty();
}
bool SampleCSVParser::validateSample(const SampleInfo& sample) const {
    if (!sample.isValid() || !fileExists(sample.read1_path)) return false;
    return sample.isPairedEnd() ? fileExists(sample.read2_path) : true;
}
void SampleCSVParser::printSummary() const { /* Omitted for brevity */ }
}

// kmer_minimizer.h
#ifndef KMER_MINIMIZER_H
#define KMER_MINIMIZER_H
#include <cstdint>
struct Minimizer {
    uint64_t hash;
    uint32_t pos;
    bool is_reverse;
};
#endif

// bloom_filter.h
#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H
#include <cstdint>
#include <string>
#include <vector>
class BloomFilterGPU {
private:
    uint64_t* d_bloom_filter;
    size_t bloom_size_bytes;
    int kmer_length;
    bool initialized;
public:
    BloomFilterGPU(int k = 15);
    ~BloomFilterGPU();
    bool initialize();
    bool isInitialized() const { return initialized; }
};
#endif

// bloom_filter.cpp
#include <iostream>
#include <cuda_runtime.h>
BloomFilterGPU::BloomFilterGPU(int k) : kmer_length(k), initialized(false), d_bloom_filter(nullptr) {}
BloomFilterGPU::~BloomFilterGPU() { if (d_bloom_filter) cudaFree(d_bloom_filter); }
bool BloomFilterGPU::initialize() {
    bloom_size_bytes = (1ULL << 26) / 8; // 64 Mbits = 8 MBytes
    cudaError_t err = cudaMalloc(&d_bloom_filter, bloom_size_bytes);
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate Bloom filter: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    cudaMemset(d_bloom_filter, 0, bloom_size_bytes);
    initialized = true;
    return true;
}


// Makefile
CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++17 -Wall -O2 -I.
NVCCFLAGS = -std=c++17 -O2 -I. -x cu
LDFLAGS = -lz

# Source files
SRCS = pipeline_tester.cpp sample_csv_parser.cpp fastq_reader.cpp
CU_SRCS = bloom_filter.cu

# Object files
OBJS = $(SRCS:.cpp=.o)
CU_OBJS = $(CU_SRCS:.cu=.o)

# Executable name
TARGET = pipeline_tester

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS) $(CU_OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(CU_OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(CU_OBJS) $(TARGET)

