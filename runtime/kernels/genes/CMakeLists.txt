# AMR Gene Detection Pipeline CMakeLists.txt
cmake_minimum_required(VERSION 3.18)
project(AMRGeneDetectionPipeline LANGUAGES C CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set CUDA architectures
set(CMAKE_CUDA_ARCHITECTURES 61 70 75 80 86)

# Find required dependencies
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

if(NOT HDF5_FOUND)
    message(FATAL_ERROR "HDF5 not found. Please install HDF5 development libraries.")
endif()

if(NOT JSONCPP_FOUND)
    message(FATAL_ERROR "jsoncpp not found. Please install jsoncpp development libraries (libjsoncpp-dev).")
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../shared)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${JSONCPP_INCLUDE_DIRS})
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Link directories
link_directories(${JSONCPP_LIBRARY_DIRS})

# Add the shared library from parent directory
add_subdirectory(../shared ${CMAKE_CURRENT_BINARY_DIR}/shared)

# ===========================
# Core Libraries
# ===========================

# AMR detection kernels library
add_library(amr_detection_kernels STATIC
    amr_detection_kernels.cu
    amr_detection_kernels_wrapper.cu
)

set_target_properties(amr_detection_kernels PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Translated search library
add_library(translated_search_amr STATIC
    translated_search_amr.cu
)

set_target_properties(translated_search_amr PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# HDF5 AMR writer library
add_library(hdf5_amr_writer STATIC
    hdf5_amr_writer.cpp
)

target_link_libraries(hdf5_amr_writer
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Clinical report generator library
add_library(clinical_amr_report_generator STATIC
    clinical_amr_report_generator.cpp
)

target_link_libraries(clinical_amr_report_generator
    ${JSONCPP_LIBRARIES}
)

# Build AMR protein database executable
add_executable(build_amr_protein_database
    build_amr_protein_database.cpp
)

target_link_libraries(build_amr_protein_database
    ${JSONCPP_LIBRARIES}
)

# ===========================
# Main Executable
# ===========================

# AMR detection main executable
add_executable(amr_detection
    amr_detection_main.cpp
    amr_detection_pipeline.cpp
)

target_compile_options(amr_detection PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        --use_fast_math
        --generate-line-info
        -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        --ptxas-options=-v
        --maxrregcount=64
        -O3
    >
)

target_link_libraries(amr_detection
    biogpu_shared  # Use shared components library
    amr_detection_kernels
    translated_search_amr
    hdf5_amr_writer
    clinical_amr_report_generator
    z
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${JSONCPP_LIBRARIES}
)