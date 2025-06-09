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

# Standalone FASTQ pipeline test (if we build it separately)
# Uncomment when fastq_gpu_pipeline.cpp has its own main()
# add_executable(test_fastq_pipeline
#     fastq_gpu_pipeline.cpp
# )
# 
# target_link_libraries(test_fastq_pipeline
#     minimizer_extraction
#     cudart
#     cuda
#     ${ZLIB_LIBRARIES}
# )

# Hybrid profiler pipeline - currently just a test
# Uncomment when we have proper separation of concerns
# add_executable(hybrid_profiler_pipeline
#     hybrid_profiler_profile.cpp
# )
# 
# target_link_libraries(hybrid_profiler_pipeline
#     minimizer_extraction
#     cudart
#     cuda
#     ${ZLIB_LIBRARIES}
# )
# 
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(hybrid_profiler_pipeline OpenMP::OpenMP_CXX)
# endif()

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

# ===========================
# Installation Rules
# ===========================

# Install executables
install(TARGETS 
    test_minimizer
    RUNTIME DESTINATION bin
)

# Install libraries
install(TARGETS 
    minimizer_extraction
    ARCHIVE DESTINATION lib
)

# Install headers
install(FILES 
    minimizer_common.h
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
message(STATUS "  ✓ Test program for minimizers")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  test_minimizer - Test GPU minimizer extraction")
message(STATUS "")
message(STATUS "Testing targets:")
message(STATUS "  make test_synthetic - Test with synthetic data")
message(STATUS "  make test_real_fastq - Test with real FASTQ (set FASTQ_FILE)")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  1. Build:")
message(STATUS "     mkdir build && cd build")
message(STATUS "     cmake ..")
message(STATUS "     make")
message(STATUS "")
message(STATUS "  2. Test with synthetic data:")
message(STATUS "     ./test_minimizer")
message(STATUS "")
message(STATUS "  3. Test with your FASTQ file:")
message(STATUS "     ./test_minimizer /path/to/reads.fastq")
message(STATUS "     # or")
message(STATUS "     FASTQ_FILE=/path/to/reads.fastq make test_real_fastq")
message(STATUS "")