
cmake_minimum_required(VERSION 3.20)
project(PipelineComparison LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find packages
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(ZLIB REQUIRED)

# Include directories from main project
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(../include)
include_directories(../runtime/kernels/resistance)

# CUDA include directories
include_directories(/usr/local/cuda/include)
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA settings
set(CMAKE_CUDA_ARCHITECTURES "61;70;75;80;86")
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

# Add the comparison test executable
add_executable(pipeline_comparison_test
    ../pipeline_comparison_test.cpp
    ../runtime/kernels/resistance/fq_mutation_detector.cu
    ../runtime/kernels/resistance/bloom_filter.cu
    ../runtime/kernels/resistance/kmer_screening.cu
    ../runtime/kernels/resistance/translated_search.cu
    ../runtime/kernels/resistance/enhanced_mutation_detection.cu
    ../runtime/kernels/resistance/diagnostic_report.cpp
    ../runtime/kernels/resistance/hdf5_alignment_writer.cpp
)

set_target_properties(pipeline_comparison_test PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Link libraries
target_link_libraries(pipeline_comparison_test
    ${ZLIB_LIBRARIES}
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)
