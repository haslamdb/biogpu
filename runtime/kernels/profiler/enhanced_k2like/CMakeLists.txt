cmake_minimum_required(VERSION 3.18)
project(BioGPU LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find CUDA
find_package(CUDA REQUIRED)
if(NOT CUDA_FOUND)
    message(FATAL_ERROR "CUDA not found. Please install CUDA toolkit.")
endif()

# Check CUDA compute capability
if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "11.0")
    message(WARNING "CUDA version ${CMAKE_CUDA_COMPILER_VERSION} may not support all features")
endif()

# Organize source files by module
file(GLOB_RECURSE CORE_SOURCES "core/*.cu" "core/*.cpp")
file(GLOB_RECURSE MEMORY_SOURCES "memory/*.cu" "memory/*.cpp")
file(GLOB_RECURSE GPU_SOURCES "gpu/*.cu" "gpu/*.cpp")
file(GLOB_RECURSE TAXONOMY_SOURCES "taxonomy/*.cu" "taxonomy/*.cpp")
file(GLOB_RECURSE PROCESSING_SOURCES "processing/*.cu" "processing/*.cpp")
file(GLOB_RECURSE OUTPUT_SOURCES "output/*.cu" "output/*.cpp")
file(GLOB_RECURSE TOOLS_SOURCES "tools/*.cu" "tools/*.cpp")
file(GLOB ROOT_SOURCES "*.cu" "*.cpp")

# Exclude main files and test files
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER MEMORY_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER GPU_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER TAXONOMY_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER PROCESSING_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER OUTPUT_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER TOOLS_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*main.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*test.*")
# Exclude problematic files temporarily
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_kraken_database_builder.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_kraken_classifier.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*phase1_enhanced_classification.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_minimizer_extraction.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_database_kernels.*")

# Combine all sources
set(ALL_SOURCES ${CORE_SOURCES} ${MEMORY_SOURCES} ${GPU_SOURCES} ${TAXONOMY_SOURCES} ${PROCESSING_SOURCES} ${OUTPUT_SOURCES} ${TOOLS_SOURCES} ${ROOT_SOURCES})

# Create main library
add_library(biogpu ${ALL_SOURCES})

# Set CUDA properties
set_property(TARGET biogpu PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET biogpu PROPERTY CUDA_ARCHITECTURES 70 75 80)

# CUDA compile options
target_compile_options(biogpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    $<$<COMPILE_LANGUAGE:CUDA>:--relocatable-device-code=true>
    $<$<COMPILE_LANGUAGE:CUDA>:--maxrregcount=32>
)

# Include directories
target_include_directories(biogpu PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/tools
)

# Link CUDA libraries
target_link_libraries(biogpu ${CUDA_LIBRARIES})

# Simple test executable
add_executable(test_simple test_simple.cu)
target_link_libraries(test_simple biogpu)

# Core module test
add_executable(test_core_only test_core_only.cu)
target_link_libraries(test_core_only biogpu)

# Taxonomy module test
add_executable(test_taxonomy_only test_taxonomy_only.cu)
target_link_libraries(test_taxonomy_only biogpu)

# Database serializer test
add_executable(test_database_serializer test_database_serializer.cu)
target_link_libraries(test_database_serializer biogpu)
set_property(TARGET test_database_serializer PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_database_serializer PROPERTY CUDA_ARCHITECTURES 70 75 80)

# Add to test suite
enable_testing()
add_test(NAME DatabaseSerializerTest COMMAND test_database_serializer)

# Print configuration info
message(STATUS "CUDA Version: ${CMAKE_CUDA_COMPILER_VERSION}")
message(STATUS "CUDA Architectures: 70 75 80")
message(STATUS "Source files found: ${ALL_SOURCES}")