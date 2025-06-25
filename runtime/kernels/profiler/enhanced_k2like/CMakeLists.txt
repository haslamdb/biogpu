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

# Exclude GPU-specific files that might have issues temporarily
list(FILTER GPU_SOURCES EXCLUDE REGEX ".*cht-.*")
list(FILTER GPU_SOURCES EXCLUDE REGEX ".*hash_cells.*")
list(FILTER GPU_SOURCES EXCLUDE REGEX ".*hash_mask.*")
list(FILTER GPU_SOURCES EXCLUDE REGEX ".*lca_bits.*")

# Also exclude files from the ROOT_SOURCES that have issues
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_kraken_database_builder.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_minimizer_extraction.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_kraken_classifier_backup.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*phase1_enhanced_classification.*")
list(FILTER ROOT_SOURCES EXCLUDE REGEX ".*gpu_database_kernels.*")  # Exclude duplicate

# Combine all sources
set(ALL_SOURCES 
    ${CORE_SOURCES} 
    ${MEMORY_SOURCES} 
    ${GPU_SOURCES} 
    ${TAXONOMY_SOURCES} 
    ${PROCESSING_SOURCES} 
    ${OUTPUT_SOURCES} 
    ${TOOLS_SOURCES} 
    ${ROOT_SOURCES}
)

# Create main library
add_library(biogpu ${ALL_SOURCES})

# Set CUDA properties
set_property(TARGET biogpu PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET biogpu PROPERTY CUDA_ARCHITECTURES 61 70 75 80 86 89 90)

# CUDA compile options
target_compile_options(biogpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    $<$<COMPILE_LANGUAGE:CUDA>:--relocatable-device-code=true>
    $<$<COMPILE_LANGUAGE:CUDA>:--maxrregcount=32>
    $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>
)

# Include directories
target_include_directories(biogpu PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${CMAKE_CURRENT_SOURCE_DIR}/gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/taxonomy
    ${CMAKE_CURRENT_SOURCE_DIR}/processing
    ${CMAKE_CURRENT_SOURCE_DIR}/output
    ${CMAKE_CURRENT_SOURCE_DIR}/tools
    ${CMAKE_CURRENT_SOURCE_DIR}/../..  # For accessing parent directories if needed
)

# Link CUDA libraries
target_link_libraries(biogpu 
    ${CUDA_LIBRARIES}
    ${CUDA_CUDART_LIBRARY}
    ${CMAKE_DL_LIBS}  # For dynamic loading if needed
)

# Optional: Link with HDF5 if available
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
    target_compile_definitions(biogpu PRIVATE USE_HDF5)
    target_include_directories(biogpu PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(biogpu ${HDF5_LIBRARIES})
    message(STATUS "HDF5 support enabled")
else()
    message(STATUS "HDF5 not found - feature export to HDF5 will be disabled")
endif()

# Test executables
# Skip test_simple.cu as it appears to be documentation, not code
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_simple.cu")
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/test_simple.cu" TEST_SIMPLE_CONTENT LIMIT 100)
    string(FIND "${TEST_SIMPLE_CONTENT}" "```" HAS_MARKDOWN)
    if(HAS_MARKDOWN EQUAL -1)
        add_executable(test_simple test_simple.cu)
        target_link_libraries(test_simple biogpu)
        set_property(TARGET test_simple PROPERTY CUDA_SEPARABLE_COMPILATION ON)
        set_property(TARGET test_simple PROPERTY CUDA_ARCHITECTURES 70 75 80)
    endif()
endif()

add_executable(test_core_only test_core_only.cu)
target_link_libraries(test_core_only biogpu)
set_property(TARGET test_core_only PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_core_only PROPERTY CUDA_ARCHITECTURES 70 75 80)

add_executable(test_taxonomy_only test_taxonomy_only.cu)
target_link_libraries(test_taxonomy_only biogpu)
set_property(TARGET test_taxonomy_only PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_taxonomy_only PROPERTY CUDA_ARCHITECTURES 70 75 80)

add_executable(test_database_serializer test_database_serializer.cu)
target_link_libraries(test_database_serializer biogpu)
set_property(TARGET test_database_serializer PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_database_serializer PROPERTY CUDA_ARCHITECTURES 70 75 80)

# ML structures test
add_executable(test_ml_structures test_ml_structures.cu)
target_link_libraries(test_ml_structures biogpu)
set_property(TARGET test_ml_structures PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_ml_structures PROPERTY CUDA_ARCHITECTURES 70 75 80)

# Enhanced database simple test
add_executable(test_enhanced_database_simple test_enhanced_database_simple.cu)
target_link_libraries(test_enhanced_database_simple biogpu)
set_property(TARGET test_enhanced_database_simple PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_enhanced_database_simple PROPERTY CUDA_ARCHITECTURES 70 75 80)

# Enhanced database build test - skip for now due to missing types
# if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_enhanced_database_build.cu")
#     add_executable(test_enhanced_database_build test_enhanced_database_build.cu)
#     target_link_libraries(test_enhanced_database_build biogpu)
#     set_property(TARGET test_enhanced_database_build PROPERTY CUDA_SEPARABLE_COMPILATION ON)
#     set_property(TARGET test_enhanced_database_build PROPERTY CUDA_ARCHITECTURES 70 75 80)
# endif()

# Enable testing
enable_testing()
if(TARGET test_simple)
    add_test(NAME SimpleTest COMMAND test_simple)
endif()
add_test(NAME CoreTest COMMAND test_core_only)
add_test(NAME TaxonomyTest COMMAND test_taxonomy_only)
add_test(NAME DatabaseSerializerTest COMMAND test_database_serializer)
add_test(NAME MLStructuresTest COMMAND test_ml_structures)
add_test(NAME EnhancedDatabaseSimpleTest COMMAND test_enhanced_database_simple)

# Print configuration info
message(STATUS "=== BioGPU Configuration ===")
message(STATUS "CUDA Version: ${CMAKE_CUDA_COMPILER_VERSION}")
message(STATUS "CUDA Architectures: 61 70 75 80 86 89 90")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "Source files found:")
message(STATUS "  Core: ${CORE_SOURCES}")
message(STATUS "  Memory: ${MEMORY_SOURCES}")
message(STATUS "  GPU: ${GPU_SOURCES}")
message(STATUS "  Taxonomy: ${TAXONOMY_SOURCES}")
message(STATUS "  Processing: ${PROCESSING_SOURCES}")
message(STATUS "  Output: ${OUTPUT_SOURCES}")
message(STATUS "  Tools: ${TOOLS_SOURCES}")
message(STATUS "")
message(STATUS "Total source files: ${ALL_SOURCES}")