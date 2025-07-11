cmake_minimum_required(VERSION 3.18)
project(BioGPU LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find CUDA
find_package(CUDA REQUIRED)

# Set CUDA compilation flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")

# --- Source File Organization ---
# We now explicitly list the core modules to be included.

set(BIOGPU_SOURCES
    # Core Orchestrator
    core/gpu_database_builder_core.cu

    # GPU Kernels
    gpu/gpu_database_kernels.cu

    # Memory Management
    memory/gpu_memory_manager.cu

    # Processing Modules
    processing/genome_file_processor.cu
    processing/minimizer_feature_extractor.cu
    processing/contamination_detector.cu
    processing/feature_exporter.cu

    # Taxonomy Module
    taxonomy/taxonomy_processor_stub.cu  # Using stub for now

    # Output Module
    output/database_serializer.cu
)

# --- Library Definition ---
add_library(biogpu ${BIOGPU_SOURCES})

# --- Target Properties ---
set_property(TARGET biogpu PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET biogpu PROPERTY CUDA_ARCHITECTURES 61 70 75 80 86) # Including Pascal (61) for Titan Xp
set_property(TARGET biogpu PROPERTY POSITION_INDEPENDENT_CODE ON)

# --- Include Directories ---
target_include_directories(biogpu PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# --- Link Libraries ---
target_link_libraries(biogpu
    ${CUDA_LIBRARIES}
    ${CUDA_CUDART_LIBRARY}
)

# --- Main Executable ---
# Temporarily disabled - needs classifier refactoring
# add_executable(kraken_pipeline kraken_pipeline_main.cu)
# target_link_libraries(kraken_pipeline biogpu)
# set_property(TARGET kraken_pipeline PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# --- Test Executable ---
add_executable(test_db_builder test_database_builder.cu)
target_link_libraries(test_db_builder biogpu)
set_property(TARGET test_db_builder PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET test_db_builder PROPERTY CUDA_ARCHITECTURES 61 70 75 80 86)

# --- Testing ---
enable_testing()
add_test(NAME DatabaseBuildTest COMMAND test_db_builder)

# --- Print Configuration ---
message(STATUS "Configured BioGPU with a modular architecture.")
message(STATUS "Main executable: kraken_pipeline")
message(STATUS "Test executable: test_db_builder")