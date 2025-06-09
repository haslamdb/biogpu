cmake_minimum_required(VERSION 3.20)
project(BioGPU_Profiling VERSION 0.1.0 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard  
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Project include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# CUDA include directories
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA library directories
link_directories(/usr/local/cuda/lib64)
link_directories(/usr/local/cuda/targets/x86_64-linux/lib)

# Find OpenMP for parallel processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP: ${OpenMP_CXX_VERSION}")
else()
    message(WARNING "OpenMP not found - some features may be slower")
endif()

# Find zlib for compressed FASTQ support
find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    message(STATUS "Found zlib: ${ZLIB_VERSION_STRING}")
else()
    message(FATAL_ERROR "zlib not found. Please install zlib development libraries.")
endif()

# ===========================
# Core Profiling Libraries
# ===========================

# Minimizer extraction library
add_library(minimizer_extraction STATIC
    minimizer_extraction.cu
    minimizer_common.h
)

set_target_properties(minimizer_extraction PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(minimizer_extraction PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

# FASTQ processing library
add_library(fastq_processing STATIC
    fastq_gpu_pipeline.cpp
)

target_link_libraries(fastq_processing
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fastq_processing OpenMP::OpenMP_CXX)
endif()

# Kraken2 database interface library
add_library(kraken_interface STATIC
    kraken_database_loader.cpp
    kraken_database_loader.h
)

target_link_libraries(kraken_interface
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(kraken_interface OpenMP::OpenMP_CXX)
endif()

# Hybrid CPU-GPU classification library
add_library(hybrid_classification STATIC
    hybrid_classification.cu
    hybrid_classification.h
)

set_target_properties(hybrid_classification PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(hybrid_classification PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(hybrid_classification
    kraken_interface
    minimizer_extraction
    cudart
)

# ===========================
# Test Executables
# ===========================

# Test minimizer extraction
add_executable(test_minimizer
    test_minimizer.cpp
)

target_link_libraries(test_minimizer
    minimizer_extraction
    cudart
    cuda
)

# Test FASTQ pipeline
add_executable(test_fastq_pipeline
    test_fastq_main.cpp
)

target_link_libraries(test_fastq_pipeline
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(test_fastq_pipeline OpenMP::OpenMP_CXX)
endif()

# Test Kraken database loading
add_executable(test_kraken_loader
    test_kraken_loader.cpp
)

target_link_libraries(test_kraken_loader
    kraken_interface
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(test_kraken_loader OpenMP::OpenMP_CXX)
endif()

# ===========================
# Main Pipeline Executables
# ===========================

# Hybrid profiler pipeline (CPU database + GPU minimizers)
add_executable(hybrid_profiler_pipeline
    hybrid_profiler_main.cpp
)

target_link_libraries(hybrid_profiler_pipeline
    hybrid_classification
    kraken_interface
    minimizer_extraction
    fastq_processing
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(hybrid_profiler_pipeline OpenMP::OpenMP_CXX)
endif()

# GPU minimizer profiler (standalone)
add_executable(gpu_minimizer_profiler
    gpu_profiler_main.cpp
)

target_link_libraries(gpu_minimizer_profiler
    minimizer_extraction
    fastq_processing
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(gpu_minimizer_profiler OpenMP::OpenMP_CXX)
endif()

# ===========================
# Custom Targets for Testing
# ===========================

# Test with synthetic data
add_custom_target(test_synthetic
    COMMAND ./test_minimizer
    DEPENDS test_minimizer
    COMMENT "Testing minimizer extraction with synthetic sequences"
)

# Test with actual FASTQ (user must provide FASTQ_FILE)
add_custom_target(test_real_fastq
    COMMAND ./test_fastq_pipeline $ENV{FASTQ_FILE}
    DEPENDS test_fastq_pipeline
    COMMENT "Testing FASTQ pipeline - set FASTQ_FILE environment variable"
)

# Test Kraken database loading
add_custom_target(test_kraken_db
    COMMAND ./test_kraken_loader $ENV{KRAKEN_DB_PATH}
    DEPENDS test_kraken_loader
    COMMENT "Testing Kraken database loading - set KRAKEN_DB_PATH environment variable"
)

# Benchmark minimizer performance
add_custom_target(benchmark_minimizers
    COMMAND ./test_minimizer --benchmark
    DEPENDS test_minimizer
    COMMENT "Benchmarking GPU minimizer extraction performance"
)

# ===========================
# Installation Rules
# ===========================

# Install executables
install(TARGETS 
    test_minimizer
    test_fastq_pipeline
    test_kraken_loader
    hybrid_profiler_pipeline
    gpu_minimizer_profiler
    RUNTIME DESTINATION bin
)

# Install libraries
install(TARGETS 
    minimizer_extraction
    fastq_processing
    kraken_interface
    hybrid_classification
    ARCHIVE DESTINATION lib
)

# Install headers
install(FILES 
    minimizer_common.h
    kraken_database_loader.h
    hybrid_classification.h
    DESTINATION include/biogpu/profiler
)

# ===========================
# Configuration Summary
# ===========================

message(STATUS "")
message(STATUS "BioGPU Microbiome Profiling Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86")
message(STATUS "  OpenMP Support: ${OpenMP_CXX_FOUND}")
message(STATUS "  zlib Version: ${ZLIB_VERSION_STRING}")
message(STATUS "")
message(STATUS "Profiling Pipeline Components:")
message(STATUS "  ✓ GPU minimizer extraction")
message(STATUS "  ✓ FASTQ processing with compression support")
message(STATUS "  ✓ Kraken2 database interface")
message(STATUS "  ✓ Hybrid CPU-GPU classification")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  test_minimizer           - Test GPU minimizer extraction")
message(STATUS "  test_fastq_pipeline      - Test FASTQ processing")
message(STATUS "  test_kraken_loader       - Test Kraken database loading")
message(STATUS "  hybrid_profiler_pipeline   - Main profiling pipeline")
message(STATUS "  gpu_minimizer_profiler   - Standalone GPU profiler")
message(STATUS "")
message(STATUS "Testing targets:")
message(STATUS "  make test_synthetic      - Test with synthetic data")
message(STATUS "  make test_real_fastq     - Test with real FASTQ (set FASTQ_FILE)")
message(STATUS "  make test_kraken_db      - Test Kraken DB (set KRAKEN_DB_PATH)")
message(STATUS "  make benchmark_minimizers - Performance benchmarking")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  1. Test minimizer extraction:")
message(STATUS "     ./test_minimizer")
message(STATUS "")
message(STATUS "  2. Test with your FASTQ file:")
message(STATUS "     FASTQ_FILE=/path/to/reads.fastq.gz make test_real_fastq")
message(STATUS "")
message(STATUS "  3. Run hybrid profiling pipeline:")
message(STATUS "     ./hybrid_profiler_pipeline /path/to/kraken_db reads.fastq.gz output.txt")
message(STATUS "")
message(STATUS "  4. Test Kraken database loading:")
message(STATUS "     KRAKEN_DB_PATH=/path/to/kraken_db make test_kraken_db")
message(STATUS "")