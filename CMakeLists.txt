cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.4.0 LANGUAGES CXX CUDA)

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
# Libraries
# ===========================

# Bloom filter library
add_library(bloom_filter STATIC
    runtime/kernels/resistance/bloom_filter.cu
)

set_target_properties(bloom_filter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# HDF5 alignment writer library (updated with translated search support)
add_library(hdf5_writer STATIC
    runtime/kernels/resistance/hdf5_alignment_writer.cpp
)

target_link_libraries(hdf5_writer
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# NEW: Translated search library
add_library(translated_search STATIC
    runtime/kernels/resistance/translated_search.cu
)

set_target_properties(translated_search PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# ===========================
# FQ Resistance Pipeline (GPU)
# ===========================

# FQ pipeline CUDA sources - includes ALL necessary kernels
set(FQ_CUDA_SOURCES
    runtime/kernels/resistance/kmer_screening.cu
    runtime/kernels/resistance/fq_mutation_detector.cu  # Contains identity calculation fix
)

# FQ pipeline host sources
set(FQ_HOST_SOURCES
    runtime/kernels/resistance/fq_pipeline_host.cpp
)

# Create FQ pipeline executable
add_executable(fq_pipeline_gpu
    ${FQ_HOST_SOURCES}
    ${FQ_CUDA_SOURCES}
)

# Set CUDA properties for FQ pipeline
set_target_properties(fq_pipeline_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"  # Include Titan Xp (61)
)

# CUDA-specific compile options
target_compile_options(fq_pipeline_gpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
    >
)

# Link libraries for FQ pipeline (added translated_search)
target_link_libraries(fq_pipeline_gpu
    bloom_filter         # Bloom filter library
    hdf5_writer         # HDF5 writer library
    translated_search   # NEW: Translated search library
    z                   # for zlib (FASTQ.gz reading)
    cudart
    cuda
    ${HDF5_LIBRARIES}      # HDF5 libraries
    ${HDF5_CXX_LIBRARIES}  # HDF5 C++ libraries
)

# Add HDF5 definitions
target_compile_definitions(fq_pipeline_gpu PRIVATE ${HDF5_DEFINITIONS})

# ===========================
# Test and Debug Executables
# ===========================

# Debug tool for gene/species ID analysis
add_executable(debug_ids
    tests/debug_gene_species_ids.cpp
)
target_link_libraries(debug_ids cudart cuda)

# Test program for Bloom filter
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

# Test program for HDF5 alignment reader
add_executable(test_hdf5_reader
    tests/test_hdf5_reader.cpp
)

target_link_libraries(test_hdf5_reader
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# NEW: Test program for translated search
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
# Development Targets
# ===========================

# Build k-mer index with correct paths
add_custom_target(build_kmer_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py 
            ${CMAKE_SOURCE_DIR}/output/GeneFiles
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            --kmer-length 15
    COMMENT "Building enhanced k-mer index"
)

# NEW: Build protein resistance database
add_custom_target(build_protein_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_protein_resistance_db.py
            ${CMAKE_SOURCE_DIR}/data/fq_genes
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/protein_resistance_db
    COMMENT "Building protein resistance database for translated search"
)

# Test Bloom filter functionality
add_custom_target(test_bloom
    COMMAND ./test_bloom_filter
    DEPENDS test_bloom_filter
    COMMENT "Testing Bloom filter implementation"
)

# Debug gene/species ID assignments
add_custom_target(debug_index
    COMMAND ./debug_ids ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
    DEPENDS debug_ids
    COMMENT "Debug gene/species ID assignments in k-mer index"
)

# Validate existing index
add_custom_target(validate_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            --create-synthetic
            --test-data ${CMAKE_SOURCE_DIR}/data/test_fastq
    COMMENT "Validating k-mer index with synthetic data"
)

# Test HDF5 output reading
add_custom_target(test_hdf5
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/read_hdf5_alignments.py
            ${CMAKE_BINARY_DIR}/test_output.h5
            --summarize --plot
    COMMENT "Testing HDF5 output reading and visualization"
)

# NEW: Test translated search functionality
add_custom_target(test_translation
    COMMAND ./test_translated_search
    DEPENDS test_translated_search
    COMMENT "Testing 6-frame translation and protein search"
)

# ===========================
# Installation Rules
# ===========================

# Install executables (added test_translated_search)
install(TARGETS fq_pipeline_gpu debug_ids test_bloom_filter test_hdf5_reader test_translated_search
    RUNTIME DESTINATION bin
)

# Install libraries (added translated_search)
install(TARGETS bloom_filter hdf5_writer translated_search
    ARCHIVE DESTINATION lib
)

# Install Python scripts (added protein DB builder)
install(PROGRAMS
    ${CMAKE_SOURCE_DIR}/src/python/read_hdf5_alignments.py
    ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py
    ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
    ${CMAKE_SOURCE_DIR}/src/python/build_protein_resistance_db.py
    DESTINATION bin
)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ===========================
# Configuration Summary
# ===========================

message(STATUS "")
message(STATUS "BioGPU Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION} (added translated search support)")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86 (includes Titan Xp)")
message(STATUS "  HDF5 Version: ${HDF5_VERSION}")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  fq_pipeline_gpu    - GPU-accelerated FQ resistance detection with translated search")
message(STATUS "  debug_ids          - Debug gene/species ID assignments")
message(STATUS "  test_bloom_filter  - Test Bloom filter implementation")
message(STATUS "  test_hdf5_reader   - Test HDF5 alignment reading")
message(STATUS "  test_translated_search - Test 6-frame translation and protein search")
message(STATUS "")
message(STATUS "Development targets:")
message(STATUS "  build_kmer_index - Build k-mer index from gene files")
message(STATUS "  build_protein_db - Build protein database for translated search")
message(STATUS "  validate_index   - Validate index with synthetic data")
message(STATUS "  test_bloom       - Run Bloom filter tests")
message(STATUS "  debug_index      - Debug gene/species ID assignments")
message(STATUS "  test_hdf5        - Test HDF5 output reading and visualization")
message(STATUS "  test_translation - Test 6-frame translation functionality")
message(STATUS "")