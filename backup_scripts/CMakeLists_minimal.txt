# Enhanced CMakeLists.txt for BioGPU project
cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.2.0 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Project include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/runtime/kernels/resistance)

# CUDA include directories
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA library directories
link_directories(/usr/local/cuda/lib64)
link_directories(/usr/local/cuda/targets/x86_64-linux/lib)

# Enhanced FQ pipeline with fixed k-mer screening
set(ENHANCED_CUDA_SOURCES
    runtime/kernels/resistance/fixed_kmer_screening.cu
    # runtime/kernels/resistance/fq_mutation_detector.cu  # Temporarily disabled - uses HDF5 format
)

# Host code for pipeline
set(ENHANCED_HOST_SOURCES
    runtime/kernels/resistance/fq_pipeline_host.cpp
)

# Create enhanced FQ pipeline executable
add_executable(fq_pipeline_gpu_enhanced
    ${ENHANCED_HOST_SOURCES}
    ${ENHANCED_CUDA_SOURCES}
)

# Set CUDA properties for enhanced pipeline
set_target_properties(fq_pipeline_gpu_enhanced PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"  # Include Titan Xp (61)
)

# CUDA-specific compile options for enhanced pipeline
target_compile_options(fq_pipeline_gpu_enhanced PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
    >
)

# Link libraries for enhanced pipeline
target_link_libraries(fq_pipeline_gpu_enhanced
    z          # for zlib (FASTQ.gz reading)
    cudart
    cuda
)

# Test executable for k-mer screening validation
# Commented out until test file is created
# add_executable(test_kmer_screening
#     tests/test_kmer_screening.cpp
#     ${ENHANCED_CUDA_SOURCES}
# )
# 
# set_target_properties(test_kmer_screening PROPERTIES
#     CUDA_SEPARABLE_COMPILATION ON
#     CUDA_ARCHITECTURES "61;70;75;80;86"
# )
# 
# target_link_libraries(test_kmer_screening
#     cudart
#     cuda
# )

# Convenience targets for development
add_custom_target(test_enhanced_pipeline
    COMMAND ./fq_pipeline_gpu_enhanced
    DEPENDS fq_pipeline_gpu_enhanced
    COMMENT "Testing enhanced FQ resistance pipeline"
)

add_custom_target(validate_kmer_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py 
            ${CMAKE_SOURCE_DIR}/data/sequences 
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_db/fq_mutations.csv
            ${CMAKE_SOURCE_DIR}/build/enhanced_index
    COMMENT "Building and validating enhanced k-mer index"
)

add_custom_target(test_synthetic_data
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
            ${CMAKE_SOURCE_DIR}/build/enhanced_index
            --create-synthetic
    DEPENDS validate_kmer_index
    COMMENT "Testing with synthetic data"
)

add_custom_target(run_complete_test
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/test_pipeline.py
            --sequences-dir ${CMAKE_SOURCE_DIR}/data/sequences
            --mutations-csv ${CMAKE_SOURCE_DIR}/data/fq_resistance_db/fq_mutations.csv
            --test-data ${CMAKE_SOURCE_DIR}/data/test_fastq
            --output-dir ${CMAKE_SOURCE_DIR}/build/pipeline_test_results
    DEPENDS fq_pipeline_gpu_enhanced
    COMMENT "Running complete pipeline test"
)

# Installation rules
install(TARGETS fq_pipeline_gpu_enhanced
    RUNTIME DESTINATION bin
)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Print configuration summary
message(STATUS "")
message(STATUS "Enhanced BioGPU Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86 (includes Titan Xp)")
message(STATUS "  Enhanced Features: YES")
message(STATUS "")

message(STATUS "Available targets:")
message(STATUS "  fq_pipeline_gpu_enhanced - Main enhanced pipeline")
message(STATUS "  test_kmer_screening - K-mer screening validation")
message(STATUS "  validate_kmer_index - Build and validate index")
message(STATUS "  test_synthetic_data - Test with synthetic data")
message(STATUS "  run_complete_test - Full pipeline test")
message(STATUS "")