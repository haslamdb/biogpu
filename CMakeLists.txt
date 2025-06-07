cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.5.0 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS C CXX)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_VERSION}")
    include_directories(${HDF5_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "HDF5 not found. Please install HDF5 development libraries.")
endif()

# Find jsoncpp
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)
if(JSONCPP_FOUND)
    message(STATUS "Found jsoncpp: ${JSONCPP_VERSION}")
    include_directories(${JSONCPP_INCLUDE_DIRS})
    link_directories(${JSONCPP_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "jsoncpp not found. Please install jsoncpp development libraries (libjsoncpp-dev).")
endif()

# Project include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/runtime/kernels/resistance)

# CUDA include directories
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA library directories
link_directories(/usr/local/cuda/lib64)
link_directories(/usr/local/cuda/targets/x86_64-linux/lib)

# ===========================
# Core Libraries
# ===========================

# Bloom filter library
add_library(bloom_filter STATIC
    runtime/kernels/resistance/bloom_filter.cu
)

set_target_properties(bloom_filter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Enhanced HDF5 alignment writer library
add_library(hdf5_writer STATIC
    runtime/kernels/resistance/hdf5_alignment_writer.cpp
)

target_link_libraries(hdf5_writer
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Enhanced translated search library
add_library(translated_search STATIC
    runtime/kernels/resistance/translated_search_revised.cu
)

set_target_properties(translated_search PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# K-mer screening library
add_library(kmer_screening STATIC
    runtime/kernels/resistance/kmer_screening.cu
)

set_target_properties(kmer_screening PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Diagnostic reporting library
add_library(diagnostic_reporter STATIC
    runtime/kernels/resistance/diagnostic_report.cpp
)

# Enhanced mutation detection library
add_library(enhanced_mutation_detection STATIC
    runtime/kernels/resistance/enhanced_mutation_detection.cu
)

set_target_properties(enhanced_mutation_detection PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# ===========================
# NEW: Clean Resistance Loader Library
# ===========================

# Add the fix_pipeline_resistance_loader as a library
add_library(clean_resistance_loader STATIC
    runtime/kernels/resistance/fix_pipeline_resistance_loader.cpp
)

# No special dependencies needed for this loader

# ===========================
# UPDATED: Integrated Resistance Pipeline Libraries
# ===========================

# GPU-accelerated resistance detection library
add_library(resistance_detection_gpu STATIC
    runtime/kernels/resistance/resistance_detection_gpu.cu
)

set_target_properties(resistance_detection_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Keep the original fixed_database_loader for compatibility
add_library(fixed_database_loader STATIC
    runtime/kernels/resistance/fixed_database_loader.cpp
)

target_link_libraries(fixed_database_loader
    ${JSONCPP_LIBRARIES}
)

# ===========================
# Executables
# ===========================

# Enhanced FQ pipeline (original)
add_executable(fq_pipeline_gpu
    runtime/kernels/resistance/fq_pipeline_host.cpp
    runtime/kernels/resistance/fq_mutation_detector.cu
)

set_target_properties(fq_pipeline_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(fq_pipeline_gpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(fq_pipeline_gpu
    bloom_filter
    hdf5_writer
    translated_search
    kmer_screening
    diagnostic_reporter
    enhanced_mutation_detection
    clean_resistance_loader    # Add the new loader
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# ===========================
# UPDATED: Integrated Resistance Pipeline Executable
# ===========================

add_executable(integrated_resistance_pipeline
    runtime/kernels/resistance/integrate_resistance_pipeline.cpp
    runtime/kernels/resistance/fq_mutation_detector.cu
)

set_target_properties(integrated_resistance_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(integrated_resistance_pipeline PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(integrated_resistance_pipeline
    bloom_filter
    hdf5_writer
    translated_search
    kmer_screening
    resistance_detection_gpu
    fixed_database_loader        # Keep original for species-gene tracking
    clean_resistance_loader      # Add the new simplified loader
    diagnostic_reporter
    enhanced_mutation_detection
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${JSONCPP_LIBRARIES}
)

target_compile_definitions(integrated_resistance_pipeline PRIVATE ${HDF5_DEFINITIONS})

# ===========================
# NEW: Clean Resistance Pipeline Executable
# ===========================

# Create a new executable that uses the clean resistance loader
add_executable(clean_resistance_pipeline
    runtime/kernels/resistance/clean_resistance_pipeline_main.cpp
    runtime/kernels/resistance/fq_mutation_detector.cu
)

set_target_properties(clean_resistance_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(clean_resistance_pipeline PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(clean_resistance_pipeline
    bloom_filter
    hdf5_writer
    translated_search
    kmer_screening
    resistance_detection_gpu
    clean_resistance_loader      # Use only the clean loader
    diagnostic_reporter
    enhanced_mutation_detection
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# ===========================
# Test Executables
# ===========================


# Test program for resistance detection GPU
add_executable(test_resistance_gpu
    tests/test_resistance_detection_gpu.cu
)

set_target_properties(test_resistance_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(test_resistance_gpu
    resistance_detection_gpu
    cudart
    cuda
)

# Other test executables (unchanged)
add_executable(test_bloom_filter
    runtime/kernels/resistance/test_bloom_filter.cu
)

set_target_properties(test_bloom_filter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(test_bloom_filter
    bloom_filter
    cudart
    cuda
)

add_executable(test_hdf5_reader
    tests/test_hdf5_reader.cpp
)

target_link_libraries(test_hdf5_reader
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

add_executable(test_translated_search
    tests/test_translated_search.cu
)

set_target_properties(test_translated_search PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(test_translated_search
    translated_search
    cudart
    cuda
)

# ===========================
# NEW: Clean Database Building Targets
# ===========================

# Build clean resistance database (position-only, no drug info)
add_custom_target(build_clean_resistance_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_clean_resistance_db.py
            --fasta-dir ${CMAKE_SOURCE_DIR}/data/fq_genes
            --csv ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            --output-dir ${CMAKE_SOURCE_DIR}/data/clean_resistance_db
            --position-only
    COMMENT "Building clean resistance database with position information only"
)

# ===========================
# UPDATED: Integrated Database Building Targets
# ===========================

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

# Test clean resistance pipeline
add_custom_target(test_clean_pipeline
    COMMAND ./clean_resistance_pipeline
            ${CMAKE_SOURCE_DIR}/data/clean_resistance_db
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R2.fastq.gz
            ${CMAKE_BINARY_DIR}/clean_pipeline_output
    DEPENDS clean_resistance_pipeline build_clean_resistance_db
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
    COMMAND ./integrated_resistance_pipeline
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/nucleotide
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db/protein
            ${CMAKE_SOURCE_DIR}/data/integrated_resistance_db
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R1.fastq.gz
            ${CMAKE_SOURCE_DIR}/data/test_fastq/synthetic_R2.fastq.gz
            ${CMAKE_BINARY_DIR}/integrated_test_output
    DEPENDS integrated_resistance_pipeline build_integrated_resistance_db
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
    DEPENDS integrated_resistance_pipeline
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

# Install executables
install(TARGETS 
    fq_pipeline_gpu
    integrated_resistance_pipeline
    clean_resistance_pipeline      # NEW
    test_resistance_gpu
    test_bloom_filter 
    test_hdf5_reader 
    test_translated_search
    RUNTIME DESTINATION bin
)

# Install libraries
install(TARGETS 
    bloom_filter 
    hdf5_writer 
    translated_search
    kmer_screening
    diagnostic_reporter
    enhanced_mutation_detection
    resistance_detection_gpu
    fixed_database_loader
    clean_resistance_loader        # NEW
    ARCHIVE DESTINATION lib
)

# Install Python scripts
install(PROGRAMS
    ${CMAKE_SOURCE_DIR}/src/python/build_integrated_resistance_db.py
    ${CMAKE_SOURCE_DIR}/src/python/build_clean_resistance_db.py      # NEW
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
message(STATUS "BioGPU Integrated Resistance Pipeline Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86")
message(STATUS "  HDF5 Version: ${HDF5_VERSION}")
message(STATUS "")
message(STATUS "NEW: Clean Resistance Pipeline Components:")
message(STATUS "  ✓ Simplified resistance loader (position-only)")
message(STATUS "  ✓ Clean resistance pipeline executable")
message(STATUS "  ✓ Position-focused mutation detection")
message(STATUS "")
message(STATUS "Integrated Pipeline Components:")
message(STATUS "  ✓ Multi-stage resistance detection")
message(STATUS "  ✓ Unified resistance database")
message(STATUS "  ✓ Population-level variant calling")
message(STATUS "  ✓ Clinical interpretation module")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  clean_resistance_pipeline     - Simplified pipeline using position-only DB")
message(STATUS "  integrated_resistance_pipeline - Complete 3-stage detection pipeline")
message(STATUS "  build_clean_resistance_db     - Build simplified resistance database")
message(STATUS "  build_integrated_resistance_db - Build unified resistance database")
message(STATUS "  test_clean_pipeline           - Test simplified pipeline")
message(STATUS "  test_integrated_pipeline      - Test with synthetic data")
message(STATUS "  benchmark_integrated          - Performance comparison")
message(STATUS "  analyze_resistance_profile    - Generate HTML reports")
message(STATUS "  clinical_interpretation       - Create clinical reports")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  1. Build clean database: make build_clean_resistance_db")
message(STATUS "  2. Run clean pipeline:")
message(STATUS "     ./clean_resistance_pipeline \\")
message(STATUS "       data/clean_resistance_db \\")
message(STATUS "       reads_R1.fastq.gz reads_R2.fastq.gz output")
message(STATUS "")
message(STATUS "  3. Build integrated database: make build_integrated_resistance_db")
message(STATUS "  4. Run integrated pipeline:")
message(STATUS "     ./integrated_resistance_pipeline \\")
message(STATUS "       data/integrated_resistance_db/nucleotide \\")
message(STATUS "       data/integrated_resistance_db/protein \\")
message(STATUS "       data/integrated_resistance_db \\")
message(STATUS "       reads_R1.fastq.gz reads_R2.fastq.gz output")
message(STATUS "")