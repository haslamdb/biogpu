/*
 * test_database_builder.cu
 *
 * An integration test for the GPUKrakenDatabaseBuilder class.
 * This test simulates a real build process with synthetic data to ensure
 * the end-to-end pipeline is stable and functional.
 *
 * How to Compile and Run:
 * 1. Add the following to your CMakeLists.txt:
 * add_executable(test_db_builder test_database_builder.cu)
 * target_link_libraries(test_db_builder biogpu) # If you have a library
 * # Or link sources directly if not using a library
 *
 * 2. Compile:
 * cmake .
 * make test_db_builder
 *
 * 3. Run:
 * ./test_db_builder
 */

#include "gpu_kraken_database_builder.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>

// Helper function to create a dummy FASTA file for testing
bool create_dummy_genome_file(const std::string& filepath, const std::string& sequence) {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Failed to create dummy file: " << filepath << std::endl;
        return false;
    }
    // Header format includes taxid for easy parsing by the builder
    outfile << ">taxid_12345|dummy_genome\n";
    outfile << sequence << "\n";
    outfile.close();
    return true;
}

int main() {
    std::cout << "--- Running GPUKrakenDatabaseBuilder Integration Test ---" << std::endl;

    // 1. Setup temporary directories for the test
    const std::string temp_genome_dir = "temp_test_genomes";
    const std::string temp_db_dir = "temp_test_db_output";

    try {
        std::filesystem::create_directory(temp_genome_dir);
        std::filesystem::create_directory(temp_db_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Failed to create temporary directories: " << e.what() << std::endl;
        return 1;
    }

    // 2. Create synthetic genome files
    std::cout << "Creating synthetic genome data..." << std::endl;
    std::string genome_seq_1 = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
    std::string genome_seq_2 = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
    
    if (!create_dummy_genome_file(temp_genome_dir + "/genome1_taxid_1001.fasta", genome_seq_1)) return 1;
    if (!create_dummy_genome_file(temp_genome_dir + "/genome2_taxid_1002.fasta", genome_seq_2)) return 1;

    // 3. Instantiate and run the database builder
    bool success = false;
    try {
        std::cout << "Instantiating GPUKrakenDatabaseBuilder..." << std::endl;
        ClassificationParams params;
        params.k = 31; // Use standard params for testing
        params.ell = 15;

        GPUKrakenDatabaseBuilder builder(temp_db_dir, params);

        std::cout << "Initializing CUDA context..." << std::endl;
        if (!builder.initialize_cuda_context()) {
            std::cerr << "TEST FAILED: CUDA context initialization failed." << std::endl;
            throw std::runtime_error("CUDA init failed");
        }

        std::cout << "Starting database build..." << std::endl;
        if (builder.build_database_from_genomes(temp_genome_dir)) {
            std::cout << "Database build process completed." << std::endl;
            success = true;
        } else {
            std::cerr << "TEST FAILED: build_database_from_genomes returned false." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "TEST FAILED: An exception was thrown during the build process: " << e.what() << std::endl;
        success = false;
    }

    // 4. Verify results (simple check for output files)
    if (success) {
        bool hash_exists = std::filesystem::exists(temp_db_dir + "/hash_table.k2d");
        bool tax_exists = std::filesystem::exists(temp_db_dir + "/taxonomy.tsv");

        if (hash_exists && tax_exists) {
            std::cout << "VERIFICATION PASSED: Output files were created." << std::endl;
        } else {
            std::cerr << "VERIFICATION FAILED: Output files are missing." << std::endl;
            success = false;
        }
    }

    // 5. Cleanup
    std::cout << "Cleaning up temporary files and directories..." << std::endl;
    try {
        std::filesystem::remove_all(temp_genome_dir);
        std::filesystem::remove_all(temp_db_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Cleanup failed: " << e.what() << std::endl;
        // Don't fail the test for cleanup issues
    }

    if (success) {
        std::cout << "\n--- TEST PASSED ---" << std::endl;
        return 0;
    } else {
        std::cout << "\n--- TEST FAILED ---" << std::endl;
        return 1;
    }
}
