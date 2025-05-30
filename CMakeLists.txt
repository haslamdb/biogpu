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

# ===========================
# FQ Resistance Pipeline (GPU)
# ===========================

# FQ pipeline CUDA sources - includes ALL necessary kernels
set(FQ_CUDA_SOURCES
    runtime/kernels/resistance/fixed_kmer_screening.cu
    runtime/kernels/resistance/fq_mutation_detector.cu  # Contains identity calculation fix
)

# FQ pipeline host sources
set(FQ_HOST_SOURCES
    runtime/kernels/resistance/fq_pipeline_host.cpp
)

# Create FQ pipeline executable
add_executable(fq_pipeline_gpu
    ${FQ_HOST_SOURCES}
    ${FQ_CUDA_SOURCES}
)

# Set CUDA properties for FQ pipeline
set_target_properties(fq_pipeline_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"  # Include Titan Xp (61)
)

# CUDA-specific compile options
target_compile_options(fq_pipeline_gpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
    >
)

# Link libraries for FQ pipeline
target_link_libraries(fq_pipeline_gpu
    z          # for zlib (FASTQ.gz reading)
    cudart
    cuda
)

# ===========================
# Test and Debug Executables
# ===========================

# Debug tool for gene/species ID analysis
add_executable(debug_ids
    tests/debug_gene_species_ids.cpp
)
target_link_libraries(debug_ids cudart cuda)

# ===========================
# Development Targets
# ===========================

# Build k-mer index with correct paths
add_custom_target(build_kmer_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py 
            ${CMAKE_SOURCE_DIR}/output/GeneFiles
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            --kmer-length 15
    COMMENT "Building enhanced k-mer index"
)

# Test GPU pipeline on VRE12 sample
add_custom_target(test_gpu_pipeline
    COMMAND ./fq_pipeline_gpu
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index/kmer_index.bin
            ${CMAKE_SOURCE_DIR}/data/test_fastq/VRE12_R1.fastq.gz
            ${CMAKE_SOURCE_DIR}/data/test_fastq/VRE12_R2.fastq.gz
            ${CMAKE_SOURCE_DIR}/build/gpu_results.json
    DEPENDS fq_pipeline_gpu
    COMMENT "Testing GPU pipeline on VRE12 sample with fixed index"
)

# Debug gene/species ID assignments
add_custom_target(debug_index
    COMMAND ./debug_ids ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
    DEPENDS debug_ids
    COMMENT "Debug gene/species ID assignments in k-mer index"
)

# Validate existing index
add_custom_target(validate_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            --create-synthetic
            --test-data ${CMAKE_SOURCE_DIR}/data/test_fastq
    COMMENT "Validating k-mer index with synthetic data"
)

# ===========================
# Installation Rules
# ===========================

# Install executables
install(TARGETS fq_pipeline_gpu debug_ids
    RUNTIME DESTINATION bin
)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ===========================
# Configuration Summary
# ===========================

message(STATUS "")
message(STATUS "BioGPU Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86 (includes Titan Xp)")
message(STATUS "")
message(STATUS "Build targets:")
message(STATUS "  fq_pipeline_gpu  - GPU-accelerated FQ resistance detection")
message(STATUS "  debug_ids        - Debug gene/species ID assignments")
message(STATUS "")
message(STATUS "Development targets:")
message(STATUS "  build_kmer_index - Build k-mer index from gene files")
message(STATUS "  validate_index   - Validate index with synthetic data")
message(STATUS "  test_gpu_pipeline - Test GPU pipeline on VRE12 sample")
message(STATUS "  debug_index      - Debug gene/species ID assignments")
message(STATUS "")