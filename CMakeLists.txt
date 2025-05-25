cmake_minimum_required(VERSION 3.20)
project(BioGPU VERSION 0.1.0 LANGUAGES CXX CUDA)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA Standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Find LLVM
find_package(LLVM REQUIRED CONFIG)

message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")

# LLVM definitions and include directories
add_definitions(${LLVM_DEFINITIONS})
include_directories(${LLVM_INCLUDE_DIRS})
link_directories(${LLVM_LIBRARY_DIRS})

# Find ANTLR
find_package(antlr4-runtime REQUIRED)

# Project include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(/usr/include/antlr4-runtime)
include_directories(${CMAKE_SOURCE_DIR}/src/parser)  # For generated headers

# Source files
set(SOURCES
    src/compiler/compiler.cpp
    src/codegen/llvm_codegen.cpp
    src/codegen/ast_codegen.cpp
)

# Add parser source files
set(PARSER_SOURCES
    src/parser/BioGPULexer.cpp
    src/parser/BioGPUParser.cpp
    src/parser/BioGPUBaseVisitor.cpp
    src/parser/biogpu_parser.cpp
    src/parser/ast.cpp
)

# CUDA source files (if any)
set(CUDA_SOURCES
    # Add .cu files here when implemented
    # src/kernels/mutation_detection.cu
    # src/kernels/alignment.cu
)

# Combine all sources
list(APPEND SOURCES ${PARSER_SOURCES})

# Create main executable
add_executable(biogpu ${SOURCES})

# If we have CUDA sources, create a library
if(CUDA_SOURCES)
    add_library(biogpu_kernels ${CUDA_SOURCES})
    target_link_libraries(biogpu biogpu_kernels)
endif()

# Link libraries
target_link_libraries(biogpu LLVM-18 antlr4-runtime)

# Enable warnings (less strict for now)
target_compile_options(biogpu PRIVATE 
    -Wall
    $<$<CONFIG:DEBUG>:-g -O0>
    $<$<CONFIG:RELEASE>:-O3>
)

# Set CUDA architecture (adjust for your GPU)
set_target_properties(biogpu PROPERTIES 
    CUDA_ARCHITECTURES "60;70;75;80;86"
)

# Tests
# enable_testing()
# add_subdirectory(tests)

# Installation rules
install(TARGETS biogpu
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/biogpu
    DESTINATION include
)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)