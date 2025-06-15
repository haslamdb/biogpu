# BioGPU Root CMakeLists.txt
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(BioGPU LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set CUDA standard 
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Find required packages
find_package(Threads REQUIRED)

# Check CUDA compute capability
include(CheckLanguage)
check_language(CUDA)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 60 70 75 80 86)  # Support multiple GPU generations
endif()

# Set compilation flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo -use_fast_math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

# Debug vs Release
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -G")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
else()
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/runtime/kernels)
include_directories(${CUDA_INCLUDE_DIRS})

# Add subdirectories
add_subdirectory(runtime/kernels)

# Optional: Add tests if you want them
option(BUILD_TESTS "Build test programs" OFF)
if(BUILD_TESTS)
    enable_testing()
    # add_subdirectory(tests)  # Commented out until tests directory exists
endif()

# Optional: Add examples
option(BUILD_EXAMPLES "Build example programs" OFF)
if(BUILD_EXAMPLES)
    # add_subdirectory(examples)  # Commented out until examples directory exists
endif()

# Print configuration summary
message(STATUS "BioGPU Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "  CUDA flags: ${CMAKE_CUDA_FLAGS}")
message(STATUS "  CXX flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "  Build tests: ${BUILD_TESTS}")
message(STATUS "  Build examples: ${BUILD_EXAMPLES}")

# Create config header (commented out until config.h.in exists)
# configure_file(
#     ${CMAKE_SOURCE_DIR}/include/biogpu/config.h.in
#     ${CMAKE_BINARY_DIR}/include/biogpu/config.h
# )

# Install configuration
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
        
# install(FILES ${CMAKE_BINARY_DIR}/include/biogpu/config.h
#         DESTINATION include/biogpu)

# ===========================
# Database Building Targets
# ===========================

# Build clean resistance database (position-only, no drug info)
add_custom_target(build_clean_resistance_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_clean_resistance_db.py
            --fasta-dir ${CMAKE_SOURCE_DIR}/data/fq_genes
            --csv ${CMAKE_SOURCE_DIR}/data/quinolone_resistance_mutation_table.csv
            --output-dir ${CMAKE_SOURCE_DIR}/data/clean_resistance_db
            --position-only
    COMMENT "Building clean resistance database with position information only"
)

# Build integrated resistance database (combines all data sources)
add_custom_target(build_integrated_resistance_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_integrated_resistance_db.py
            --fasta-dir ${CMAKE_SOURCE_DIR}/data/fq_genes
            --csv ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            --output-dir ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db
            --add-manual
    COMMENT "Building integrated resistance database with nucleotide, protein, and mutation data"
)

# Build individual database components
add_custom_target(build_nucleotide_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py 
            ${CMAKE_SOURCE_DIR}/data/fq_genes
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/nucleotide
            --kmer-length 15
    COMMENT "Building nucleotide k-mer index"
    DEPENDS build_integrated_resistance_db
)

add_custom_target(build_protein_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_protein_resistance_db.py
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/protein
            --kmer-length 5
            --enable-minimizers
    COMMENT "Building protein database with 5-mer index and minimizers"
    DEPENDS build_integrated_resistance_db
)

add_custom_target(build_resistance_catalog
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/compile_resistance_catalog.py
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/resistance_db.json
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/resistance_catalog.bin
            --include-pmqr
            --include-efflux
    COMMENT "Compiling comprehensive resistance mutation catalog"
    DEPENDS build_integrated_resistance_db
)

# ===========================
# Validation and Testing Targets
# ===========================

# Add custom target to validate FQ database
add_custom_target(validate_fq_database
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/resistance/test_fq_mapper
            ${CMAKE_SOURCE_DIR}/data/quinolone_resistance_mutation_table.csv
            ${CMAKE_SOURCE_DIR}/data/integrated_clean_db/protein
    DEPENDS runtime/kernels/resistance/test_fq_mapper
    COMMENT "Validating FQ resistance database integration"
)

# Test clean resistance pipeline
add_custom_target(test_clean_pipeline
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/resistance/clean_resistance_pipeline
            ${CMAKE_SOURCE_DIR}/data/clean_resistance_db
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R2.fastq.gz
            ${CMAKE_BINARY_DIR}/clean_pipeline_output
    DEPENDS runtime/kernels/resistance/clean_resistance_pipeline build_clean_resistance_db
    COMMENT "Testing clean resistance pipeline with position-only database"
)

# Validate integrated database
add_custom_target(validate_integrated_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/validate_integrated_db.py
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db
            --check-consistency
            --test-queries
            --generate-report
    COMMENT "Validating integrated resistance database consistency"
    DEPENDS build_integrated_resistance_db
)

# Test integrated pipeline with synthetic data
add_custom_target(test_integrated_pipeline
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/resistance/integrated_resistance_pipeline
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/nucleotide
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/protein
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R2.fastq.gz
            ${CMAKE_BINARY_DIR}/integrated_test_output
    DEPENDS runtime/kernels/resistance/integrated_resistance_pipeline build_integrated_resistance_db
    COMMENT "Testing integrated pipeline with synthetic resistance data"
)

# Benchmark integrated vs individual pipelines
add_custom_target(benchmark_integrated
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/benchmark_integrated.py
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db
            ${CMAKE_SOURCE_DIR}/data/test_fastq
            --compare-methods
            --plot-results
            --output ${CMAKE_BINARY_DIR}/benchmark_report.html
    COMMENT "Benchmarking integrated pipeline performance"
    DEPENDS runtime/kernels/resistance/integrated_resistance_pipeline
)

# ===========================
# Kraken2-style Taxonomic Profiling Targets
# ===========================

# Build Kraken2-style database from genomes
add_custom_target(build_kraken_database
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_pipeline
            build
            --genome-dir ${CMAKE_SOURCE_DIR}/data/reference_genomes
            --output ${CMAKE_SOURCE_DIR}/data/kraken_database
            --k 35 --minimizer-len 31 --spaces 7
    DEPENDS runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_pipeline
    COMMENT "Building Kraken2-style taxonomic database from reference genomes"
)

# Test Kraken2-style classification
add_custom_target(test_kraken_classification
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_pipeline
            classify
            --database ${CMAKE_SOURCE_DIR}/data/kraken_database
            --reads ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            --output ${CMAKE_BINARY_DIR}/kraken_results.txt
            --report ${CMAKE_BINARY_DIR}/kraken_summary.txt
            --confidence 0.1 --batch-size 10000
    DEPENDS runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_pipeline build_kraken_database
    COMMENT "Testing Kraken2-style taxonomic classification"
)

# Complete Kraken2-style pipeline (build + classify)
add_custom_target(test_complete_kraken_pipeline
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_pipeline
            pipeline
            --genome-dir ${CMAKE_SOURCE_DIR}/data/reference_genomes
            --reads ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            --output ${CMAKE_BINARY_DIR}/complete_kraken_test
            --confidence 0.1
    DEPENDS runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_pipeline
    COMMENT "Testing complete Kraken2-style pipeline (build database + classify reads)"
)

# Test paired-end Kraken2-style classification
add_custom_target(test_kraken_paired_classification
    COMMAND ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_paired_classifier
            --database ${CMAKE_SOURCE_DIR}/data/kraken_database
            --reads1 ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            --reads2 ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R2.fastq.gz
            --output ${CMAKE_BINARY_DIR}/kraken_paired_results.txt
            --confidence 0.1 --paired-weight 2.0
    DEPENDS runtime/kernels/profiler/k2like_from_scratch/gpu_kraken_paired_classifier build_kraken_database
    COMMENT "Testing paired-end Kraken2-style classification with concordance scoring"
)

# ===========================
# Analysis and Reporting Targets
# ===========================

# Generate resistance profile report
add_custom_target(analyze_resistance_profile
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/resistance_profile_analyzer.py
            ${CMAKE_BINARY_DIR}/integrated_test_output.h5
            ${CMAKE_BINARY_DIR}/integrated_test_output.json
            --generate-html
            --include-alignments
            --output ${CMAKE_BINARY_DIR}/resistance_report.html
    COMMENT "Generating comprehensive resistance profile report"
)

# Create clinical interpretation
add_custom_target(clinical_interpretation
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/clinical_interpreter.py
            ${CMAKE_BINARY_DIR}/integrated_test_output.json
            --drug-database ${CMAKE_SOURCE_DIR}/data/drug_resistance_db.json
            --guidelines CLSI
            --output ${CMAKE_BINARY_DIR}/clinical_report.pdf
    COMMENT "Generating clinical interpretation report"
)

# Generate taxonomic profiling report
add_custom_target(analyze_taxonomic_profile
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/taxonomic_profile_analyzer.py
            ${CMAKE_BINARY_DIR}/kraken_results.txt
            ${CMAKE_BINARY_DIR}/kraken_summary.txt
            --generate-html
            --include-lineages
            --output ${CMAKE_BINARY_DIR}/taxonomic_report.html
    COMMENT "Generating comprehensive taxonomic profiling report"
)

# ===========================
# Installation Rules
# ===========================

# Install Python scripts
install(PROGRAMS
    ${CMAKE_SOURCE_DIR}/src/python/build_integrated_resistance_db.py
    ${CMAKE_SOURCE_DIR}/src/python/build_clean_resistance_db.py
    ${CMAKE_SOURCE_DIR}/src/python/validate_integrated_db.py
    ${CMAKE_SOURCE_DIR}/src/python/resistance_profile_analyzer.py
    ${CMAKE_SOURCE_DIR}/src/python/clinical_interpreter.py
    ${CMAKE_SOURCE_DIR}/src/python/benchmark_integrated.py
    ${CMAKE_SOURCE_DIR}/src/python/read_hdf5_alignments.py
    ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py
    ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
    ${CMAKE_SOURCE_DIR}/src/python/build_protein_resistance_db.py
    ${CMAKE_SOURCE_DIR}/src/python/diagnostic_analyzer.py
    ${CMAKE_SOURCE_DIR}/src/python/taxonomic_profile_analyzer.py
    DESTINATION bin
)

# Install data processing scripts
install(DIRECTORY
    ${CMAKE_SOURCE_DIR}/scripts/
    DESTINATION share/biogpu/scripts
    FILES_MATCHING PATTERN "*.py" PATTERN "*.sh"
)

# ===========================
# Configuration Summary
# ===========================

message(STATUS "")
message(STATUS "BioGPU Modular Pipeline Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86")
message(STATUS "")
message(STATUS "Pipeline Modules:")
message(STATUS "  ✓ Resistance Detection (runtime/kernels/resistance)")
message(STATUS "  ✓ Profiling Pipeline (runtime/kernels/profiler)")
message(STATUS "  ✓ Kraken2-style Classification (runtime/kernels/profiler/k2like_from_scratch)")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  === Resistance Detection ===")
message(STATUS "  build_clean_resistance_db     - Build simplified resistance database")
message(STATUS "  build_integrated_resistance_db - Build unified resistance database")
message(STATUS "  test_clean_pipeline           - Test simplified pipeline")
message(STATUS "  test_integrated_pipeline      - Test with synthetic data")
message(STATUS "  validate_fq_database          - Validate FQ resistance database")
message(STATUS "  benchmark_integrated          - Performance comparison")
message(STATUS "  === Taxonomic Profiling ===")
message(STATUS "  build_kraken_database         - Build Kraken2-style taxonomic database")
message(STATUS "  test_kraken_classification    - Test single-end classification")
message(STATUS "  test_kraken_paired_classification - Test paired-end classification")
message(STATUS "  test_complete_kraken_pipeline - Complete build + classify pipeline")
message(STATUS "  === Analysis and Reporting ===")
message(STATUS "  analyze_resistance_profile    - Generate HTML resistance reports")
message(STATUS "  analyze_taxonomic_profile     - Generate HTML taxonomic reports")
message(STATUS "  clinical_interpretation       - Create clinical reports")
message(STATUS "")
message(STATUS "Pipeline executables will be built in:")
message(STATUS "  ${CMAKE_BINARY_DIR}/runtime/kernels/resistance/")
message(STATUS "  ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/")
message(STATUS "  ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/k2like_from_scratch/")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  1. Build everything: make")
message(STATUS "  2. Build only resistance: make -C runtime/kernels/resistance")
message(STATUS "  3. Build only profiler: make -C runtime/kernels/profiler")
message(STATUS "  4. Build Kraken database: make build_kraken_database")
message(STATUS "  5. Test classification: make test_kraken_classification")
message(STATUS "  6. Complete pipeline: make test_complete_kraken_pipeline")
message(STATUS "")