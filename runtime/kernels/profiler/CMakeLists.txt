cmake_minimum_required(VERSION 3.20)
project(BioGPU_Enhanced_K2like VERSION 0.1.0 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard  
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Project include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/../../../include)

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

# Find Thrust (usually comes with CUDA)
find_package(Thrust CONFIG)
if(Thrust_FOUND)
    message(STATUS "Found Thrust: ${Thrust_VERSION}")
    if(Thrust_VERSION VERSION_GREATER_EQUAL "1.15.0")
        message(STATUS "Modern Thrust detected - async operations available")
        add_compile_definitions(THRUST_ASYNC_AVAILABLE)
    else()
        message(STATUS "Older Thrust detected - using synchronous operations")
    endif()
else()
    message(WARNING "Thrust not found - GPU profiler may have reduced functionality")
endif()

# Find CUB for advanced GPU algorithms
find_path(CUB_INCLUDE_DIR cub/cub.cuh
    PATHS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    NO_DEFAULT_PATH
)
if(CUB_INCLUDE_DIR)
    message(STATUS "Found CUB: ${CUB_INCLUDE_DIR}")
    include_directories(${CUB_INCLUDE_DIR})
    add_compile_definitions(CUB_AVAILABLE)
else()
    message(WARNING "CUB not found - some GPU algorithms may be slower")
endif()

# ===========================
# Core Libraries (Required)
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

target_compile_options(fastq_processing PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:
        -O3
        -march=native
    >
)

target_link_libraries(fastq_processing
    ${ZLIB_LIBRARIES}
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

target_compile_options(gpu_kmer_database PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -O3
    >
)

# ========================================
# Enhanced K2-like Implementation
# ========================================

# GPU minimizer extraction for Kraken
add_library(gpu_kraken_minimizer_extraction STATIC
    enhanced_k2like/gpu_minimizer_extraction.cu
)

set_target_properties(gpu_kraken_minimizer_extraction PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

target_compile_options(gpu_kraken_minimizer_extraction PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        -O3
        --maxrregcount=64
        -DTHRUST_IGNORE_CUB_VERSION_CHECK
        --extended-lambda
    >
)

# GPU Kraken classifier
add_library(gpu_kraken_classifier STATIC
    enhanced_k2like/gpu_kraken_classifier.cu
)

set_target_properties(gpu_kraken_classifier PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

target_compile_options(gpu_kraken_classifier PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -O3
        --maxrregcount=64
        -DTHRUST_IGNORE_CUB_VERSION_CHECK
        --extended-lambda
    >
)

target_link_libraries(gpu_kraken_classifier
    gpu_kraken_minimizer_extraction
    fastq_processing
    ${ZLIB_LIBRARIES}
    cuda
    cudart
)

if(Thrust_FOUND)
    target_link_libraries(gpu_kraken_classifier Thrust::Thrust)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(gpu_kraken_classifier OpenMP::OpenMP_CXX)
endif()

# GPU Kraken database builder
add_library(gpu_kraken_database_builder STATIC
    enhanced_k2like/gpu_kraken_database_builder.cu
)

set_target_properties(gpu_kraken_database_builder PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

target_compile_options(gpu_kraken_database_builder PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -O3
        --maxrregcount=64
        -DTHRUST_IGNORE_CUB_VERSION_CHECK
        --extended-lambda
    >
)

target_link_libraries(gpu_kraken_database_builder
    gpu_kraken_classifier
    gpu_kraken_minimizer_extraction
    fastq_processing
    ${ZLIB_LIBRARIES}
    cuda
    cudart
    stdc++fs
)

if(Thrust_FOUND)
    target_link_libraries(gpu_kraken_database_builder Thrust::Thrust)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(gpu_kraken_database_builder OpenMP::OpenMP_CXX)
endif()

# ========================================
# Enhanced K2-like Classifier with Phylogeny
# ========================================

# Enhanced classification with phylogeny library
add_library(gpu_enhanced_classification_phylo STATIC
    enhanced_k2like/phase1_enhanced_classification.cu
)

set_target_properties(gpu_enhanced_classification_phylo PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

target_compile_options(gpu_enhanced_classification_phylo PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -O3
        --maxrregcount=64
        -DTHRUST_IGNORE_CUB_VERSION_CHECK
        --extended-lambda
    >
    $<$<COMPILE_LANGUAGE:CXX>:
        -O3
        -march=native
        -DNDEBUG
    >
)

# Enhanced Kraken pipeline executable
add_executable(enhanced_kraken_pipeline
    enhanced_k2like/kraken_pipeline_main.cu
    enhanced_k2like/sample_csv_parser.cpp
)

set_target_properties(enhanced_kraken_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/enhanced_k2like
)

target_link_libraries(enhanced_kraken_pipeline
    gpu_enhanced_classification_phylo
    gpu_kraken_database_builder
    gpu_kraken_classifier
    gpu_kraken_minimizer_extraction
    fastq_processing
    ${ZLIB_LIBRARIES}
    cuda
    cudart
    stdc++fs
)

if(Thrust_FOUND)
    target_link_libraries(enhanced_kraken_pipeline Thrust::Thrust)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(enhanced_kraken_pipeline OpenMP::OpenMP_CXX)
endif()

# Compact taxonomy builder tool (if it exists)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tools/build_compact_taxonomy.cpp)
    add_executable(build_compact_taxonomy
        tools/build_compact_taxonomy.cpp
    )

    set_target_properties(build_compact_taxonomy PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tools
    )

    target_include_directories(build_compact_taxonomy PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/tools
    )

    target_link_libraries(build_compact_taxonomy
        cuda
        cudart
        stdc++fs
    )
endif()

# Status messages
message(STATUS "")
message(STATUS "Enhanced K2-like Classifier Configuration:")
message(STATUS "  Main executable: enhanced_kraken_pipeline")
message(STATUS "  Libraries:")
message(STATUS "    - gpu_enhanced_classification_phylo")
message(STATUS "    - gpu_kraken_database_builder")
message(STATUS "    - gpu_kraken_classifier")
message(STATUS "    - gpu_kraken_minimizer_extraction")
message(STATUS "  Output directory: ${CMAKE_CURRENT_BINARY_DIR}/enhanced_k2like")
message(STATUS "")
message(STATUS "Build with: cmake --build . --target enhanced_kraken_pipeline")
message(STATUS "")