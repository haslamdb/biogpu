# Resistance Detection Pipeline CMakeLists.txt
cmake_minimum_required(VERSION 3.18)
project(ResistanceDetectionPipeline LANGUAGES C CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required dependencies
find_package(HDF5 REQUIRED COMPONENTS C CXX)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_VERSION}")
    include_directories(${HDF5_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "HDF5 not found. Please install HDF5 development libraries.")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)
if(JSONCPP_FOUND)
    message(STATUS "Found jsoncpp: ${JSONCPP_VERSION}")
    include_directories(${JSONCPP_INCLUDE_DIRS})
    link_directories(${JSONCPP_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "jsoncpp not found. Please install jsoncpp development libraries (libjsoncpp-dev).")
endif()

# Include directories for resistance kernels
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add CUDA include directories for C++ files
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# ===========================
# Core Libraries
# ===========================

# CSV parser library (NEW)
add_library(sample_csv_parser STATIC
    sample_csv_parser.cpp
)

# Bloom filter library
add_library(bloom_filter STATIC
    bloom_filter.cu
)

set_target_properties(bloom_filter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Enhanced HDF5 alignment writer library
add_library(hdf5_writer STATIC
    hdf5_alignment_writer.cpp
)

target_link_libraries(hdf5_writer
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Enhanced translated search library
add_library(translated_search STATIC
    translated_search_revised.cu
)

set_target_properties(translated_search PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# K-mer screening library
add_library(kmer_screening STATIC
    kmer_screening.cu
)

set_target_properties(kmer_screening PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Diagnostic reporting library
add_library(diagnostic_reporter STATIC
    diagnostic_report.cpp
)

# Global FQ resistance mapper library
add_library(global_fq_mapper STATIC
    global_fq_resistance_mapper.cpp
)

target_link_libraries(global_fq_mapper
    cudart
)

# Enhanced mutation detection library
add_library(enhanced_mutation_detection STATIC
    enhanced_mutation_detection_unified.cu
)

set_target_properties(enhanced_mutation_detection PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(enhanced_mutation_detection
    global_fq_mapper
)

# Clean resistance loader library
add_library(clean_resistance_loader STATIC
    fix_pipeline_resistance_loader.cpp
)

# GPU-accelerated resistance detection library
add_library(resistance_detection_gpu STATIC
    resistance_detection_gpu.cu
)

set_target_properties(resistance_detection_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

# Fixed database loader for compatibility
add_library(fixed_database_loader STATIC
    fixed_database_loader.cpp
)

target_link_libraries(fixed_database_loader
    ${JSONCPP_LIBRARIES}
)

# ===========================
# Executables
# ===========================

# Enhanced FQ pipeline (original)
add_executable(fq_pipeline_gpu
    fq_pipeline_host.cpp
    fq_mutation_detector.cu
)

set_target_properties(fq_pipeline_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

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

target_link_libraries(fq_pipeline_gpu
    bloom_filter
    hdf5_writer
    translated_search
    kmer_screening
    diagnostic_reporter
    enhanced_mutation_detection
    clean_resistance_loader
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Integrated resistance pipeline executable
add_executable(integrated_resistance_pipeline
    integrate_resistance_pipeline.cpp
    fq_mutation_detector.cu
)

set_target_properties(integrated_resistance_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(integrated_resistance_pipeline PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(integrated_resistance_pipeline
    bloom_filter
    hdf5_writer
    translated_search
    kmer_screening
    resistance_detection_gpu
    fixed_database_loader
    clean_resistance_loader
    global_fq_mapper
    diagnostic_reporter
    enhanced_mutation_detection
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${JSONCPP_LIBRARIES}
)

target_compile_definitions(integrated_resistance_pipeline PRIVATE ${HDF5_DEFINITIONS})

# Clean resistance pipeline executable (WITH CSV PARSER SUPPORT)
add_executable(clean_resistance_pipeline
    clean_resistance_pipeline_main.cpp
    fq_mutation_detector.cu
    fq_mutation_reporter.cpp
    clinical_fq_report_generator.cpp
)

set_target_properties(clean_resistance_pipeline PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_compile_options(clean_resistance_pipeline PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(clean_resistance_pipeline
    sample_csv_parser    # Link the CSV parser library
    bloom_filter
    hdf5_writer
    translated_search
    kmer_screening
    resistance_detection_gpu
    clean_resistance_loader
    global_fq_mapper
    diagnostic_reporter
    enhanced_mutation_detection
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# ===========================
# Test Executables
# ===========================

# Test bloom filter (only test that exists in current directory)
add_executable(test_bloom_filter
    test_bloom_filter.cu
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

# Comment out test executables that reference missing files
# These can be uncommented when the test files are available

# # Test program for resistance detection GPU
# add_executable(test_resistance_gpu
#     ${CMAKE_SOURCE_DIR}/tests/test_resistance_detection_gpu.cu
# )
# 
# set_target_properties(test_resistance_gpu PROPERTIES
#     CUDA_SEPARABLE_COMPILATION ON
#     CUDA_ARCHITECTURES "61;70;75;80;86"
# )
# 
# target_link_libraries(test_resistance_gpu
#     resistance_detection_gpu
#     cudart
#     cuda
# )
# 
# # Test executable for FQ mapper
# add_executable(test_fq_mapper
#     ${CMAKE_SOURCE_DIR}/tests/test_fq_mapper.cpp
# )
# 
# target_link_libraries(test_fq_mapper
#     global_fq_mapper
#     cudart
# )
# 
# # Test HDF5 reader
# add_executable(test_hdf5_reader
#     ${CMAKE_SOURCE_DIR}/tests/test_hdf5_reader.cpp
# )
# 
# target_link_libraries(test_hdf5_reader
#     ${HDF5_LIBRARIES}
#     ${HDF5_CXX_LIBRARIES}
# )
# 
# # Test translated search
# add_executable(test_translated_search
#     ${CMAKE_SOURCE_DIR}/tests/test_translated_search.cu
# )
# 
# set_target_properties(test_translated_search PROPERTIES
#     CUDA_SEPARABLE_COMPILATION ON
#     CUDA_ARCHITECTURES "61;70;75;80;86"
# )
# 
# target_link_libraries(test_translated_search
#     translated_search
#     cudart
#     cuda
# )