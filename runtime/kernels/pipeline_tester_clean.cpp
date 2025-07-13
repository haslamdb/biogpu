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
 * This test uses the functions from fastq_reader to open the file(s) specified
 * in the SampleInfo struct, read all records in small batches, and report the
 * total number of reads and bases found.
 * @param sample The sample metadata, including paths to FASTQ files.
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
 * This class discovers and runs registered test functions against a set of
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