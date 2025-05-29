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

# CUDA include directories
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA library directories
link_directories(/usr/local/cuda/lib64)
link_directories(/usr/local/cuda/targets/x86_64-linux/lib)

# Source files
set(SOURCES
    src/compiler/compiler.cpp
    src/codegen/llvm_codegen.cpp
    src/codegen/ast_codegen.cpp
)

# Algorithm library files
set(ALGORITHM_SOURCES
    src/stdlib/algorithms/bloom_cascade.cpp
    src/stdlib/algorithms/minimizer_classifier.cpp
    src/stdlib/algorithms/resistance_attribution.cpp
    src/stdlib/algorithms/resistance_graph.cpp
    src/stdlib/algorithms/resistance_pssm.cpp
)

# Add parser source files
set(PARSER_SOURCES
    src/parser/BioGPULexer.cpp
    src/parser/BioGPUParser.cpp
    src/parser/BioGPUBaseVisitor.cpp
    src/parser/biogpu_parser.cpp
    src/parser/ast.cpp
)

# CUDA runtime kernel files
set(CUDA_SOURCES
    runtime/kernels/resistance/fq_resistance.cu
    runtime/kernels/resistance/fq_mutation_detector.cu
    runtime/kernels/resistance/fq_pipeline_wrapper.cu
    # Add other .cu files here when implemented
    # runtime/kernels/alignment/sequence_alignment.cu
    # runtime/kernels/indexing/kmer_index.cu
)

# Runtime system files (C++ components)
set(RUNTIME_SOURCES
    # runtime/memory/gpu_memory_manager.cpp
    # runtime/io/fastq_reader.cpp
    # Add these files as you create them
)

# Combine all sources
list(APPEND SOURCES ${PARSER_SOURCES})
list(APPEND SOURCES ${ALGORITHM_SOURCES})
list(APPEND SOURCES ${RUNTIME_SOURCES})

# Create main executable
add_executable(biogpu src/main.cpp ${SOURCES})

# Create CUDA runtime library
if(CUDA_SOURCES)
    add_library(biogpu_kernels SHARED ${CUDA_SOURCES})
    
    # Set CUDA architecture specifically for Titan Xp (compute capability 6.1)
    set_target_properties(biogpu_kernels PROPERTIES 
        CUDA_ARCHITECTURES "61;70;75;80;86"  # Added 61 for Titan Xp
    )
    
    # CUDA-specific compile options
    target_compile_options(biogpu_kernels PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:
            --use_fast_math
            --generate-line-info
            -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored
        >
    )
    
    # Link CUDA libraries
    target_link_libraries(biogpu_kernels
        ${CUDA_LIBRARIES}
        cuda
        cudart
        cublas
        curand
    )
    
    # Link kernels to main executable
    target_link_libraries(biogpu biogpu_kernels)
    
    # Export symbols for Python integration
    set_target_properties(biogpu_kernels PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CUDA_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
    )
endif()

# Link libraries
target_link_libraries(biogpu LLVM-18 antlr4-runtime)

# Enable warnings (less strict for now)
target_compile_options(biogpu PRIVATE 
    -Wall
    $<$<CONFIG:DEBUG>:-g -O0>
    $<$<CONFIG:RELEASE>:-O3>
)

# Set CUDA architecture for main executable (adjust for your GPU)
set_target_properties(biogpu PROPERTIES 
    CUDA_ARCHITECTURES "61;70;75;80;86"  # Added 61 for Titan Xp
)

# Create Python extension (optional, for direct GPU access from Python)
find_package(pybind11 QUIET)
if(pybind11_FOUND AND CUDA_SOURCES AND EXISTS "${CMAKE_SOURCE_DIR}/src/python_bindings/cuda_wrapper.cpp")
    message(STATUS "Found pybind11, creating Python bindings")
    
    pybind11_add_module(biogpu_cuda 
        src/python_bindings/cuda_wrapper.cpp
    )
    
    target_link_libraries(biogpu_cuda PRIVATE 
        biogpu_kernels
        pybind11::module
    )
    
    target_compile_definitions(biogpu_cuda PRIVATE 
        VERSION_INFO=${PROJECT_VERSION}
    )
else()
    message(STATUS "pybind11 not found or cuda_wrapper.cpp missing, skipping Python bindings")
endif()

# Create test executables for development
option(BUILD_TESTS "Build test executables" ON)
if(BUILD_TESTS)
    # Test executable for CUDA kernels
    add_executable(test_cuda_kernels 
        tests/test_cuda_kernels.cpp
    )
    
    if(CUDA_SOURCES)
        target_link_libraries(test_cuda_kernels biogpu_kernels)
    endif()
    
    # Test executable for fluoroquinolone detection
    if(EXISTS "${CMAKE_SOURCE_DIR}/tests/test_fq_detection.cpp")
        add_executable(test_fq_detection
            tests/test_fq_detection.cpp
        )
        
        if(CUDA_SOURCES)
            target_link_libraries(test_fq_detection biogpu_kernels)
        endif()
    endif()
endif()

# Development convenience targets
add_custom_target(test_titanxp
    COMMAND ./test_cuda_kernels
    DEPENDS test_cuda_kernels
    COMMENT "Testing CUDA kernels on Titan Xp"
)

if(TARGET test_fq_detection)
    add_custom_target(benchmark_fq
        COMMAND ./test_fq_detection --benchmark
        DEPENDS test_fq_detection  
        COMMENT "Benchmarking fluoroquinolone resistance detection"
    )
endif()

# Add FQ mutation detection pipeline executable
add_executable(fq_pipeline_gpu
    runtime/kernels/resistance/fq_pipeline_host.cpp
    runtime/kernels/resistance/fq_mutation_detector.cu
    runtime/kernels/resistance/fq_pipeline_wrapper.cu
)

set_target_properties(fq_pipeline_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES "61;70;75;80;86"
)

target_link_libraries(fq_pipeline_gpu
    z  # for zlib (FASTQ.gz reading)
    cudart
    cuda
)

# Find HDF5
find_package(HDF5 COMPONENTS C CXX)
if(HDF5_FOUND)
    target_include_directories(fq_pipeline_gpu PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(fq_pipeline_gpu ${HDF5_LIBRARIES})
else()
    # Manual HDF5 linking for systems where find_package doesn't work well
    target_link_libraries(fq_pipeline_gpu hdf5_serial hdf5_cpp_serial)
endif()

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

# Install CUDA library
if(CUDA_SOURCES)
    install(TARGETS biogpu_kernels
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif()

# Install Python extension
if(TARGET biogpu_cuda)
    install(TARGETS biogpu_cuda
        LIBRARY DESTINATION lib/python
    )
endif()

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Print configuration summary
message(STATUS "")
message(STATUS "BioGPU Configuration Summary:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA Architectures: 61;70;75;80;86 (includes Titan Xp)")
message(STATUS "  LLVM Version: ${LLVM_PACKAGE_VERSION}")
if(CUDA_SOURCES)
    message(STATUS "  CUDA Kernels: YES")
else()
    message(STATUS "  CUDA Kernels: NO")
endif()
if(pybind11_FOUND)
    message(STATUS "  Python Bindings: YES")
else()
    message(STATUS "  Python Bindings: NO")
endif()
message(STATUS "")