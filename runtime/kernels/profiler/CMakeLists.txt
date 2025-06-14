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
# Core Libraries
# ===========================

# Minimizer extraction library
add_library(minimizer_extraction STATIC
    minimizer_extraction.cu
    minimizer_extractor.h
    minimizer_common.h
)

set_target_properties(minimizer_extraction PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
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
    fastq_processing.cpp
    fastq_processing.h
)

target_link_libraries(fastq_processing
    minimizer_extraction
    ${ZLIB_LIBRARIES}
    cuda
    cudart
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fastq_processing OpenMP::OpenMP_CXX)
endif()

# GPU k-mer database library
add_library(gpu_kmer_database STATIC
    gpu_kmer_database.cu
    gpu_kmer_database.h
)

set_target_properties(gpu_kmer_database PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    POSITION_INDEPENDENT_CODE ON
)

# Hierarchical GPU database library - ARCHIVED
# See archived_hierarchical_db/ directory for implementation
# Removed from build as we're moving to CPU RAM streaming approach

# ===========================
# Main Executables
# ===========================

# GPU profiler pipeline
add_executable(gpu_profiler_pipeline k2_like/gpu_profiler_pipeline.cu)
target_link_libraries(gpu_profiler_pipeline 
    minimizer_extraction 
    fastq_processing
    gpu_kmer_database
    ${CUDA_LIBRARIES} 
    ${ZLIB_LIBRARIES}
)

set_target_properties(gpu_profiler_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Hierarchical profiler pipeline - ARCHIVED
# See archived_hierarchical_db/ directory

# Hybrid profiler pipeline
add_executable(hybrid_profiler_pipeline
    k2_like/hybrid_profiler_pipeline.cu
)

set_target_properties(hybrid_profiler_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(hybrid_profiler_pipeline
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(hybrid_profiler_pipeline OpenMP::OpenMP_CXX)
endif()

# Hybrid GPU/CPU profiler
add_executable(hybrid_gpu_cpu_profiler
    hybrid_gpu_cpu_profiler.cu
)

set_target_properties(hybrid_gpu_cpu_profiler PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(hybrid_gpu_cpu_profiler
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(hybrid_gpu_cpu_profiler OpenMP::OpenMP_CXX)
endif()

# Optimization for hybrid profiler
target_compile_options(hybrid_gpu_cpu_profiler PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        -O3
        -DNDEBUG
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -O3
        -march=native
        -DNDEBUG
    >
)

# ===========================
# Database Building Tools
# ===========================

# Build test database tool
add_executable(build_test_database
    k2_like/build_test_database.cpp
)

target_link_libraries(build_test_database
    gpu_kmer_database
    cudart
    cuda
)

# Build database from k-mer list
add_executable(build_db_from_kmers
    k2_like/build_db_from_kmers.cpp
)

target_link_libraries(build_db_from_kmers
    gpu_kmer_database
    cudart
    cuda
)

# Hierarchical database builder - ARCHIVED
# See archived_hierarchical_db/ directory

# Create hybrid GPU/CPU database
add_executable(create_hybrid_gpu_cpu_db
    gpu_cpu_hybrid/create_hybrid_gpu_cpu_db.cpp
)

target_link_libraries(create_hybrid_gpu_cpu_db
    gpu_kmer_database
    cudart
    cuda
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
    fastq_pipeline_main.cpp
)

set_target_properties(test_fastq_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(test_fastq_pipeline
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

# Debug minimizer program
add_executable(debug_minimizer
    debug_minimizer.cpp
)

set_target_properties(debug_minimizer PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(debug_minimizer
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

# Test stats program
add_executable(test_stats
    test_stats.cpp
)

target_link_libraries(test_stats
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

# Debug pipeline program
add_executable(fastq_pipeline_debug
    fastq_pipeline_debug.cpp
)

target_link_libraries(fastq_pipeline_debug
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)


# Debug version of hybrid GPU/CPU profiler
add_executable(hybrid_gpu_cpu_profiler_debug
    hybrid_gpu_cpu_profiler.cu
)

set_target_properties(hybrid_gpu_cpu_profiler_debug PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(hybrid_gpu_cpu_profiler_debug
    fastq_processing
    minimizer_extraction
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(hybrid_gpu_cpu_profiler_debug OpenMP::OpenMP_CXX)
endif()

# Debug compile options
target_compile_options(hybrid_gpu_cpu_profiler_debug PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        -g
        -G
        -O0
        -DDEBUG
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -g
        -O0
        -DDEBUG
    >
)

# ===========================
# Minimized Directory Executables
# ===========================

# Database analyzer
add_executable(analyze_database
    minimized/analyze_database.cu
)

set_target_properties(analyze_database PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(analyze_database
    cudart
    cuda
)

# Build minimizer database tool
add_executable(build_minimizer_db
    minimized/build_minimizer_db.cpp
)

target_link_libraries(build_minimizer_db
    ${ZLIB_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(build_minimizer_db OpenMP::OpenMP_CXX)
endif()

# GPU community profiler
add_executable(gpu_community_profiler
    minimized/gpu_community_profiler.cu
)

set_target_properties(gpu_community_profiler PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(gpu_community_profiler
    fastq_processing
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

# Optimize compile options for GPU community profiler
target_compile_options(gpu_community_profiler PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        -O3
        -DNDEBUG
        --maxrregcount=64
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -O3
        -march=native
        -DNDEBUG
    >
)

# Streaming GPU profiler
add_executable(streaming_gpu_profiler
    minimized/streaming_gpu_profiler.cu
)

set_target_properties(streaming_gpu_profiler PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(streaming_gpu_profiler
    fastq_processing
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

# Optimize compile options for streaming profiler
target_compile_options(streaming_gpu_profiler PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        -O3
        -DNDEBUG
        --maxrregcount=64
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -O3
        -march=native
        -DNDEBUG
    >
)

# Debug version of GPU community profiler
add_executable(gpu_community_profiler_debug
    minimized/gpu_community_profiler.cu
)

set_target_properties(gpu_community_profiler_debug PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(gpu_community_profiler_debug
    fastq_processing
    cudart
    cuda
    ${ZLIB_LIBRARIES}
)

# Debug compile options
target_compile_options(gpu_community_profiler_debug PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        -g
        -G
        -O0
        -DDEBUG
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -g
        -O0
        -DDEBUG
    >
)

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
    COMMAND ./test_minimizer $ENV{FASTQ_FILE}
    DEPENDS test_minimizer
    COMMENT "Testing with real FASTQ - set FASTQ_FILE environment variable"
)

# Test FASTQ pipeline
add_custom_target(test_fastq_pipeline_run
    COMMAND ./test_fastq_pipeline $ENV{FASTQ_FILE}
    DEPENDS test_fastq_pipeline
    COMMENT "Testing FASTQ pipeline - set FASTQ_FILE environment variable"
)

# Test database building and profiling
add_custom_target(test_database
    COMMAND ./build_test_database test_microbiome.db
    COMMAND ./gpu_profiler_pipeline test_microbiome.db test_microbiome.db.test.fastq
    DEPENDS build_test_database gpu_profiler_pipeline
    COMMENT "Building test database and running profiler"
)



# Test minimizer database building and GPU profiling
add_custom_target(test_minimizer_pipeline
    COMMAND ./build_minimizer_db test_genomes test_microbes_minimizer.db
    COMMAND ./gpu_community_profiler test_microbes_minimizer.db test_sample.fastq test_results
    DEPENDS build_minimizer_db gpu_community_profiler
    COMMENT "Building minimizer database and running GPU profiler"
)

# Performance benchmark
add_custom_target(benchmark_gpu_profiler
    COMMAND ./gpu_community_profiler test_microbes_minimizer.db large_sample.fastq benchmark_results
    DEPENDS gpu_community_profiler
    COMMENT "Benchmarking GPU community profiler performance"
)

# Compare CPU vs GPU profiler performance (if CPU profiler exists)
add_custom_target(compare_profilers
    COMMAND echo "Note: CPU profiler not yet implemented in minimized/"
    COMMAND echo "Running GPU profiler..."
    COMMAND time ./gpu_community_profiler test_microbes_minimizer.db test_sample.fastq gpu_results
    DEPENDS gpu_community_profiler
    COMMENT "Running GPU profiler (CPU profiler not yet available)"
)

# ===========================
# Installation Rules
# ===========================

# Install executables
install(TARGETS 
    test_minimizer
    test_fastq_pipeline
    hybrid_profiler_pipeline
    gpu_profiler_pipeline
    build_test_database
    build_db_from_kmers
    create_hybrid_gpu_cpu_db
    hybrid_gpu_cpu_profiler
    hybrid_gpu_cpu_profiler_debug
    # Minimized executables
    analyze_database
    build_minimizer_db
    gpu_community_profiler
    streaming_gpu_profiler
    gpu_community_profiler_debug
    RUNTIME DESTINATION bin
)

# Install libraries
install(TARGETS 
    minimizer_extraction
    fastq_processing
    gpu_kmer_database
    ARCHIVE DESTINATION lib
)

# Install headers
install(FILES 
    minimizer_common.h
    minimizer_extractor.h
    fastq_processing.h
    gpu_kmer_database.h
    kraken_database_loader.h
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
message(STATUS "Currently Available Components:")
message(STATUS "  ✓ GPU minimizer extraction")
message(STATUS "  ✓ GPU k-mer database")
message(STATUS "  ✓ GPU profiler pipeline")
message(STATUS "  ✓ Database building tools")
message(STATUS "  ✓ Test programs for all components")
message(STATUS "  ✓ FASTQ GPU pipeline with zlib support")
message(STATUS "  ✓ Hybrid profiler pipeline (needs Kraken2 database)")
message(STATUS "  ✓ Hybrid GPU/CPU profiler with optimizations")
message(STATUS "  ✓ Hybrid database creation tools")
message(STATUS "  ✓ NEW: Minimizer-based GPU community profiler")
message(STATUS "  ✓ NEW: Streaming GPU profiler for large databases")
message(STATUS "  ✓ NEW: Database analyzer tool")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  === K2-like Implementation ===")
message(STATUS "  gpu_profiler_pipeline       - Main GPU-accelerated profiler")
message(STATUS "  build_test_database         - Build test microbial database")
message(STATUS "  build_db_from_kmers         - Build database from k-mer list")
message(STATUS "  hybrid_profiler_pipeline    - Hybrid CPU-GPU profiler (needs Kraken2 DB)")
message(STATUS "  === Hybrid GPU/CPU Implementation ===")
message(STATUS "  hybrid_gpu_cpu_profiler     - Optimized hybrid GPU/CPU profiler")
message(STATUS "  create_hybrid_gpu_cpu_db    - Create hybrid GPU/CPU database")
message(STATUS "  === Minimizer-based Implementation ===")
message(STATUS "  gpu_community_profiler      - GPU-accelerated community profiler")
message(STATUS "  streaming_gpu_profiler      - Streaming profiler for large databases")
message(STATUS "  build_minimizer_db          - Build minimizer database")
message(STATUS "  analyze_database            - Analyze database statistics")
message(STATUS "  === Test and Debug Tools ===")
message(STATUS "  test_minimizer              - Test GPU minimizer extraction")
message(STATUS "  test_fastq_pipeline         - Full FASTQ processing pipeline")
message(STATUS "  gpu_community_profiler_debug - Debug version of community profiler")
message(STATUS "  hybrid_gpu_cpu_profiler_debug - Debug version of hybrid profiler")
message(STATUS "")
message(STATUS "Testing targets:")
message(STATUS "  make test_synthetic         - Test with synthetic data")
message(STATUS "  make test_real_fastq        - Test minimizer with real FASTQ")
message(STATUS "  make test_fastq_pipeline_run - Test full pipeline with FASTQ")
message(STATUS "  make test_database          - Build test DB and run profiler")
message(STATUS "  make test_minimizer_pipeline - Test minimizer DB building and profiling")
message(STATUS "  make benchmark_gpu_profiler - Benchmark GPU community profiler")
message(STATUS "  make compare_profilers      - Compare profiler implementations")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  1. Build:")
message(STATUS "     mkdir build && cd build")
message(STATUS "     cmake ..")
message(STATUS "     make")
message(STATUS "")
message(STATUS "  2. Quick test with synthetic database:")
message(STATUS "     make test_database")
message(STATUS "")
message(STATUS "  3. Build database from k-mer list:")
message(STATUS "     ./build_db_from_kmers kmers.txt microbiome.db")
message(STATUS "")
message(STATUS "  4. Run profiler on FASTQ:")
message(STATUS "     ./gpu_profiler_pipeline microbiome.db reads.fastq")
message(STATUS "")
message(STATUS "  5. Test with your FASTQ file:")
message(STATUS "     ./test_minimizer /path/to/reads.fastq")
message(STATUS "     # or")
message(STATUS "     FASTQ_FILE=/path/to/reads.fastq make test_real_fastq")
message(STATUS "")
message(STATUS "  6. Test hybrid profiler (requires Kraken2 database):")
message(STATUS "     ./hybrid_profiler_pipeline /path/to/kraken_db /path/to/reads.fastq")
message(STATUS "")
message(STATUS "  7. Use optimized hybrid GPU/CPU profiler:")
message(STATUS "     ./create_hybrid_gpu_cpu_db kmers.txt hybrid_gpu_cpu.db")
message(STATUS "     ./hybrid_gpu_cpu_profiler hybrid_gpu_cpu.db reads.fastq")
message(STATUS "")
message(STATUS "  8. Use minimizer-based GPU community profiler:")
message(STATUS "     ./build_minimizer_db /path/to/genomes microbes_minimizer.db")
message(STATUS "     ./gpu_community_profiler microbes_minimizer.db sample.fastq results")
message(STATUS "     ./analyze_database microbes_minimizer.db")
message(STATUS "")
message(STATUS "  9. Stream large databases with GPU acceleration:")
message(STATUS "     ./streaming_gpu_profiler large_microbes.db sample.fastq results")
message(STATUS "")