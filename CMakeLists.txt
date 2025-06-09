cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.6.0 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Project include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# CUDA include directories
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA library directories
link_directories(/usr/local/cuda/lib64)
link_directories(/usr/local/cuda/targets/x86_64-linux/lib)

# ===========================
# Add Subdirectories
# ===========================

# Add resistance detection pipeline
add_subdirectory(runtime/kernels/resistance)

# Add profiling pipeline (when ready)
add_subdirectory(runtime/kernels/profiler)

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
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  build_clean_resistance_db     - Build simplified resistance database")
message(STATUS "  build_integrated_resistance_db - Build unified resistance database")
message(STATUS "  test_clean_pipeline           - Test simplified pipeline")
message(STATUS "  test_integrated_pipeline      - Test with synthetic data")
message(STATUS "  validate_fq_database          - Validate FQ resistance database")
message(STATUS "  benchmark_integrated          - Performance comparison")
message(STATUS "  analyze_resistance_profile    - Generate HTML reports")
message(STATUS "  clinical_interpretation       - Create clinical reports")
message(STATUS "")
message(STATUS "Pipeline executables will be built in:")
message(STATUS "  ${CMAKE_BINARY_DIR}/runtime/kernels/resistance/")
message(STATUS "  ${CMAKE_BINARY_DIR}/runtime/kernels/profiler/")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  1. Build everything: make")
message(STATUS "  2. Build only resistance: make -C runtime/kernels/resistance")
message(STATUS "  3. Build only profiler: make -C runtime/kernels/profiler")
message(STATUS "")