cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.1.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find LLVM
find_package(LLVM REQUIRED CONFIG)

message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")

# LLVM definitions and include directories
add_definitions(${LLVM_DEFINITIONS})
include_directories(${LLVM_INCLUDE_DIRS})
link_directories(${LLVM_LIBRARY_DIRS})

# Find CUDA
find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

# Project include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.cpp
    src/codegen/llvm_codegen.cpp
)

# Create executable
add_executable(biogpu ${SOURCES})

# LLVM components to link
llvm_map_components_to_libnames(llvm_libs 
    core 
    support 
    irreader 
    codegen 
    mc 
    mcparser 
    option 
    target
    ${LLVM_TARGETS_TO_BUILD}
)

# Link libraries
target_link_libraries(biogpu ${llvm_libs} ${CUDA_LIBRARIES})

# Enable warnings
target_compile_options(biogpu PRIVATE -Wall -Wextra -Wpedantic)
