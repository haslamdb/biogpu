cmake_minimum_required(VERSION 3.18)
project(EnhancedK2Like LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find packages
find_package(CUDA REQUIRED)
find_package(Threads REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${CMAKE_CURRENT_SOURCE_DIR}/processing
    ${CMAKE_CURRENT_SOURCE_DIR}/gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/taxonomy
    ${CMAKE_CURRENT_SOURCE_DIR}/output
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
)

# Core module
set(CORE_SOURCES
    core/gpu_database_builder_core.cu
)

# Memory management module
set(MEMORY_SOURCES
    memory/gpu_memory_manager.cu
)

# File processing module
set(PROCESSING_SOURCES
    processing/genome_file_processor.cu
)

# GPU kernels module
set(GPU_SOURCES
    gpu/gpu_database_kernels.cu
)

# Taxonomy processing module
set(TAXONOMY_SOURCES
    taxonomy/taxonomy_processor.cu
)

# Output serialization module
set(OUTPUT_SOURCES
    output/database_serializer_temp.cu
)

# Combine all sources
set(ALL_SOURCES
    ${CORE_SOURCES}
    ${MEMORY_SOURCES}
    ${PROCESSING_SOURCES}
    ${GPU_SOURCES}
    ${TAXONOMY_SOURCES}
    ${OUTPUT_SOURCES}
)

# Create the main library
add_library(enhanced_k2like ${ALL_SOURCES})

# Set CUDA properties
set_property(TARGET enhanced_k2like PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET enhanced_k2like PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# CUDA compilation flags
target_compile_options(enhanced_k2like PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --gpu-architecture=sm_60
        --extended-lambda
        --expt-relaxed-constexpr
    >
)

# Link libraries
target_link_libraries(enhanced_k2like
    ${CUDA_LIBRARIES}
    Threads::Threads
    stdc++fs  # For filesystem support
)

# Main executable
add_executable(kraken_gpu_builder kraken_pipeline_main.cu)
target_link_libraries(kraken_gpu_builder enhanced_k2like)

# Example programs
add_executable(example_basic examples/basic_usage.cu)
target_link_libraries(example_basic enhanced_k2like)

add_executable(example_phylo examples/phylogenetic_analysis.cu)
target_link_libraries(example_phylo enhanced_k2like)

add_executable(example_high_memory examples/high_memory_config.cu)
target_link_libraries(example_high_memory enhanced_k2like)

add_executable(example_streaming examples/streaming_processing.cu)
target_link_libraries(example_streaming enhanced_k2like)

# Optional: Add other executables if needed
# add_executable(test_simple simple_test.cu)
# target_link_libraries(test_simple enhanced_k2like)