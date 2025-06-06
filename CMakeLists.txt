cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.4.1 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS C CXX)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_VERSION}")
    include_directories(${HDF5_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "HDF5 not found. Please install HDF5 development libraries.")
endif()

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
# Libraries
# ===========================

# Bloom filter library
add_library(bloom_filter STATIC
    runtime/kernels/resistance/bloom_filter.cu
)

set_target_properties(bloom_filter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Enhanced HDF5 alignment writer library with translated search support
add_library(hdf5_writer STATIC
    runtime/kernels/resistance/hdf5_alignment_writer.cpp
)

target_link_libraries(hdf5_writer
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Enhanced translated search library with 5-mer k-mers and Smith-Waterman
add_library(translated_search STATIC
    runtime/kernels/resistance/translated_search._revised.cu
)

set_target_properties(translated_search PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Enhanced k-mer screening library
add_library(kmer_screening STATIC
    runtime/kernels/resistance/kmer_screening.cu
)

set_target_properties(kmer_screening PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# NEW: Diagnostic reporting library
add_library(diagnostic_reporter STATIC
    runtime/kernels/resistance/diagnostic_report.cpp
)

# NEW: Enhanced mutation detection library
add_library(enhanced_mutation_detection STATIC
    runtime/kernels/resistance/enhanced_mutation_detection.cu
)

set_target_properties(enhanced_mutation_detection PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# ===========================
# Enhanced FQ Resistance Pipeline (GPU)
# ===========================

# Enhanced FQ pipeline CUDA sources with all necessary kernels
set(FQ_CUDA_SOURCES
    runtime/kernels/resistance/fq_mutation_detector.cu
)

# Enhanced FQ pipeline host sources with diagnostic support
set(FQ_HOST_SOURCES
    runtime/kernels/resistance/fq_pipeline_host.cpp
)

# Create enhanced FQ pipeline executable with diagnostic support
add_executable(fq_pipeline_gpu
    ${FQ_HOST_SOURCES}
    ${FQ_CUDA_SOURCES}
)

# Set CUDA properties for enhanced FQ pipeline
set_target_properties(fq_pipeline_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"  # Include Titan Xp (61)
)

# Enhanced CUDA-specific compile options
target_compile_options(fq_pipeline_gpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

# Link libraries for enhanced FQ pipeline with diagnostic support
target_link_libraries(fq_pipeline_gpu
    bloom_filter                # Bloom filter library
    hdf5_writer                # Enhanced HDF5 writer library
    translated_search          # Enhanced translated search library with 5-mer + SW
    kmer_screening             # Enhanced k-mer screening library
    diagnostic_reporter        # NEW: Diagnostic reporting
    enhanced_mutation_detection # NEW: Enhanced mutation detection
    z                          # for zlib (FASTQ.gz reading)
    cudart
    cuda
    ${HDF5_LIBRARIES}          # HDF5 libraries
    ${HDF5_CXX_LIBRARIES}      # HDF5 C++ libraries
)

# Add HDF5 definitions
target_compile_definitions(fq_pipeline_gpu PRIVATE ${HDF5_DEFINITIONS})


# ===========================
# Test and Debug Executables
# ===========================

# Debug tool for gene/species ID analysis
add_executable(debug_ids
    tests/debug_gene_species_ids.cpp
)
target_link_libraries(debug_ids cudart cuda)

# Enhanced test program for Bloom filter
add_executable(test_bloom_filter
    runtime/kernels/resistance/test_bloom_filter.cu
)

set_target_properties(test_bloom_filter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(test_bloom_filter
    bloom_filter
    cudart
    cuda
)

# Test program for HDF5 alignment reader
add_executable(test_hdf5_reader
    tests/test_hdf5_reader.cpp
)

target_link_libraries(test_hdf5_reader
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Enhanced test program for translated search (5-mer + Smith-Waterman)
add_executable(test_translated_search
    tests/test_translated_search.cu
)

set_target_properties(test_translated_search PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(test_translated_search
    translated_search
    cudart
    cuda
)

# ===========================
# Development Targets
# ===========================

# Build enhanced k-mer index with correct paths
add_custom_target(build_kmer_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py 
            ${CMAKE_SOURCE_DIR}/data/fq_genes
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            --kmer-length 15
    COMMENT "Building enhanced k-mer index for nucleotide search"
)

# Enhanced protein resistance database builder with 5-mer support
add_custom_target(build_protein_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_protein_resistance_db.py
            ${CMAKE_SOURCE_DIR}/data/fq_genes
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/protein_resistance_db
            --kmer-length 5
            --enable-smith-waterman
    COMMENT "Building enhanced protein resistance database for 5-mer translated search"
)

# NEW: Build wild-type protein database for proper mutation detection
add_custom_target(build_wildtype_protein_db
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/build_wildtype_protein_db.py
            ${CMAKE_SOURCE_DIR}/data/fq_genes
            ${CMAKE_SOURCE_DIR}/data/Known_Quinolone_Changes.csv
            ${CMAKE_SOURCE_DIR}/data/wildtype_protein_db
            --include-resistance-positions
            --kmer-length 5
    COMMENT "Building wild-type protein database for enhanced mutation detection"
)

# Test enhanced Bloom filter functionality
add_custom_target(test_bloom
    COMMAND ./test_bloom_filter
    DEPENDS test_bloom_filter
    COMMENT "Testing enhanced Bloom filter implementation"
)

# Debug gene/species ID assignments
add_custom_target(debug_index
    COMMAND ./debug_ids ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
    DEPENDS debug_ids
    COMMENT "Debug gene/species ID assignments in k-mer index"
)

# Validate existing index with enhanced features
add_custom_target(validate_index
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            --create-synthetic
            --test-data ${CMAKE_SOURCE_DIR}/data/test_fastq
            --enable-translated-search
    COMMENT "Validating k-mer index with synthetic data and enhanced features"
)

# Test HDF5 output reading with enhanced translated search results
add_custom_target(test_hdf5
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/read_hdf5_alignments.py
            ${CMAKE_BINARY_DIR}/test_output.h5
            --summarize --plot --include-translated
    COMMENT "Testing HDF5 output reading with enhanced translated search visualization"
)

# Test enhanced translated search functionality (5-mer + Smith-Waterman)
add_custom_target(test_translation
    COMMAND ./test_translated_search
    DEPENDS test_translated_search
    COMMENT "Testing enhanced 6-frame translation with 5-mer k-mers and Smith-Waterman"
)

# NEW: Test enhanced mutation detection
add_custom_target(test_mutation_detection
    COMMAND ./test_enhanced_mutation_detection
    DEPENDS test_enhanced_mutation_detection
    COMMENT "Testing enhanced mutation detection with wild-type comparison"
)

# NEW: Test diagnostic reporting
add_custom_target(test_diagnostics
    COMMAND ./test_diagnostic_report
    DEPENDS test_diagnostic_report
    COMMENT "Testing diagnostic report generation"
)

# Benchmark enhanced pipeline performance
add_custom_target(benchmark_pipeline
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/benchmark_pipeline.py
            ${CMAKE_SOURCE_DIR}/data/fq_resistance_index
            ${CMAKE_SOURCE_DIR}/data/test_fastq
            --enable-translated-search
            --enable-smith-waterman
            --enable-diagnostics
            --iterations 10
    COMMENT "Benchmarking enhanced pipeline performance with diagnostics"
)

# NEW: Run complete diagnostic analysis
add_custom_target(analyze_diagnostics
    COMMAND python3 ${CMAKE_SOURCE_DIR}/src/python/diagnostic_analyzer.py
            ${CMAKE_BINARY_DIR}/test_output_diagnostic.txt
            --json ${CMAKE_BINARY_DIR}/test_output.json
            --plot ${CMAKE_BINARY_DIR}/diagnostic_summary.png
    COMMENT "Analyzing diagnostic reports and generating summary plots"
)

# ===========================
# Installation Rules
# ===========================

# Install enhanced executables
install(TARGETS 
    fq_pipeline_gpu           # Enhanced version with full diagnostics
    debug_ids 
    test_bloom_filter 
    test_hdf5_reader 
    test_translated_search
    RUNTIME DESTINATION bin
)

# Install enhanced libraries
install(TARGETS 
    bloom_filter 
    hdf5_writer 
    translated_search
    kmer_screening
    diagnostic_reporter       # NEW
    enhanced_mutation_detection  # NEW
    ARCHIVE DESTINATION lib
)

# Install enhanced Python scripts
install(PROGRAMS
    ${CMAKE_SOURCE_DIR}/src/python/read_hdf5_alignments.py
    ${CMAKE_SOURCE_DIR}/src/python/enhanced_kmer_builder.py
    ${CMAKE_SOURCE_DIR}/src/python/index_validator.py
    ${CMAKE_SOURCE_DIR}/src/python/build_protein_resistance_db.py
    ${CMAKE_SOURCE_DIR}/src/python/build_wildtype_protein_db.py  # NEW
    ${CMAKE_SOURCE_DIR}/src/python/benchmark_pipeline.py
    ${CMAKE_SOURCE_DIR}/src/python/diagnostic_analyzer.py       # NEW
    DESTINATION bin
)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ===========================
# Configuration Summary
# ===========================

message(STATUS "")
message(STATUS "Enhanced BioGPU Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION} (enhanced with diagnostic reporting + mutation detection)")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86 (includes Titan Xp)")
message(STATUS "  HDF5 Version: ${HDF5_VERSION}")
message(STATUS "")
message(STATUS "Enhanced build targets:")
message(STATUS "  fq_pipeline_gpu_enhanced  - Enhanced GPU pipeline with diagnostic reporting")
message(STATUS "                             * Comprehensive mutation detection")
message(STATUS "                             * Detailed performance analysis")
message(STATUS "                             * Alignment visualization")
message(STATUS "                             * Troubleshooting recommendations")
message(STATUS "  fq_pipeline_gpu           - Original pipeline (backwards compatible)")
message(STATUS "  diagnostic_analyzer.py    - Post-processing diagnostic analysis")
message(STATUS "")
message(STATUS "Key diagnostic features:")
message(STATUS "  * Wild-type vs. mutant reference comparison")
message(STATUS "  * QRDR coverage analysis")
message(STATUS "  * Filter efficiency tracking")
message(STATUS "  * Alignment quality assessment")
message(STATUS "  * Automatic troubleshooting advice")
message(STATUS "")
message(STATUS "Usage examples:")
message(STATUS "  Enhanced pipeline with diagnostics:")
message(STATUS "    ./fq_pipeline_gpu_enhanced index/ reads_R1.fq.gz reads_R2.fq.gz output.json \\")
message(STATUS "      --enable-translated-search --protein-db protein_db/")
message(STATUS "")
message(STATUS "  Diagnostic analysis:")
message(STATUS "    python3 diagnostic_analyzer.py output_diagnostic.txt --json output.json --plot summary.png")
message(STATUS "")
message(STATUS "NEW FEATURES IN v0.4.1:")
message(STATUS "  ✓ Enhanced mutation detection with wild-type comparison")
message(STATUS "  ✓ Comprehensive diagnostic reporting")
message(STATUS "  ✓ Alignment visualization and quality metrics")
message(STATUS "  ✓ Automatic troubleshooting recommendations")
message(STATUS "  ✓ Post-processing analysis tools")
message(STATUS "")