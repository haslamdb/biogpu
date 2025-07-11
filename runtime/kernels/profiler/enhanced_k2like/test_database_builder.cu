/*
 * test_database_builder.cu
 *
 * An integration test for the GPUKrakenDatabaseBuilder class.
 * This test simulates a real build process with a specified directory of
 * reference genomes to ensure the end-to-end pipeline is stable and functional.
 *
 * How to Compile and Run:
 * 1. Ensure you have a directory with reference genomes (e.g., data/type_strain_reference_genomes).
 * 2. Add this file to your CMakeLists.txt:
 * add_executable(test_db_builder test_database_builder.cu)
 * target_link_libraries(test_db_builder biogpu) # Or link sources directly
 *
 * 3. Compile:
 * cmake .
 * make test_db_builder
 *
 * 4. Run:
 * ./test_db_builder
 */

#include "core/gpu_database_builder_core.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>


int main() {
    std::cout << "--- Running GPUKrakenDatabaseBuilder Integration Test ---" << std::endl;

    // 1. Define paths for the test
    const std::string concatenated_fna_file = "test_clean.fna";
    const std::string temp_db_dir = "temp_real_db_output";

    // **Pre-flight Check**: Ensure the concatenated FNA file exists
    if (!std::filesystem::exists(concatenated_fna_file)) {
        std::cerr << "TEST FAILED: The concatenated FNA file does not exist." << std::endl;
        std::cerr << "Please ensure '" << concatenated_fna_file << "' exists." << std::endl;
        return 1;
    }
    
    std::cout << "Using concatenated FNA file: " << concatenated_fna_file << std::endl;

    // Setup temporary output directory
    try {
        if (std::filesystem::exists(temp_db_dir)) {
            std::filesystem::remove_all(temp_db_dir);
        }
        std::filesystem::create_directory(temp_db_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Failed to create temporary output directory: " << e.what() << std::endl;
        return 1;
    }

    // 2. Instantiate and run the database builder
    bool success = false;
    try {
        std::cout << "Instantiating GPUKrakenDatabaseBuilder..." << std::endl;
        DatabaseBuildConfig config = GPUKrakenDatabaseBuilder::create_default_config();
        config.k_value = 35;  // Now testing with k=35
        config.ell_value = 31;
        config.memory_config.auto_scale_enabled = true;
        config.memory_config.max_memory_fraction = 80;  // Use 80% of GPU memory

        GPUKrakenDatabaseBuilder builder(temp_db_dir, config);

        // The new builder initializes CUDA context internally, so the explicit call is no longer needed.
        std::cout << "Starting database build from concatenated FNA file..." << std::endl;
        if (builder.build_database_from_streaming_fna(concatenated_fna_file)) {
            std::cout << "Database build process completed." << std::endl;
            success = true;
        } else {
            std::cerr << "TEST FAILED: build_database_from_genomes returned false." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "TEST FAILED: An exception was thrown during the build process: " << e.what() << std::endl;
        success = false;
    }

    // 3. Verify results (simple check for output files)
    if (success) {
        bool hash_exists = std::filesystem::exists(temp_db_dir + "/hash_table.k2d");
        bool tax_exists = std::filesystem::exists(temp_db_dir + "/taxonomy.tsv");

        if (hash_exists && tax_exists) {
            std::cout << "VERIFICATION PASSED: Output files were created successfully." << std::endl;
        } else {
            std::cerr << "VERIFICATION FAILED: One or more output database files are missing." << std::endl;
            if (!hash_exists) std::cerr << "Missing: " << temp_db_dir << "/hash_table.k2d" << std::endl;
            if (!tax_exists) std::cerr << "Missing: " << temp_db_dir << "/taxonomy.tsv" << std::endl;
            success = false;
        }
    }

    // 4. Cleanup the temporary output directory (only on success)
    if (success) {
        std::cout << "Cleaning up temporary output directory..." << std::endl;
        try {
            std::filesystem::remove_all(temp_db_dir);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Cleanup of output directory failed: " << e.what() << std::endl;
        }
    } else {
        std::cout << "Keeping output directory for debugging: " << temp_db_dir << std::endl;
    }

    if (success) {
        std::cout << "\n--- TEST PASSED ---" << std::endl;
        return 0;
    } else {
        std::cout << "\n--- TEST FAILED ---" << std::endl;
        return 1;
    }
}
