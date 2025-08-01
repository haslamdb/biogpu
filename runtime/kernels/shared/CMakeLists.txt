# Shared components library for BioGPU pipelines
cmake_minimum_required(VERSION 3.18)
project(biogpu_shared CUDA CXX)

# Set CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 75 80 86)

# Find CUDA
find_package(CUDA REQUIRED)
find_package(Threads REQUIRED)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

# Create shared library
add_library(biogpu_shared STATIC
    sample_csv_parser.cpp
    bloom_filter.cu
)

# Set properties for the library
set_target_properties(biogpu_shared PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
)

# Link libraries
target_link_libraries(biogpu_shared
    ${CUDA_LIBRARIES}
    Threads::Threads
)

# Export include directory for other projects
target_include_directories(biogpu_shared PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)