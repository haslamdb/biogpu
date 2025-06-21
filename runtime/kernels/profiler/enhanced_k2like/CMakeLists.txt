cmake_minimum_required(VERSION 3.18)
project(BioGPU LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

find_package(CUDA REQUIRED)

# Simple build - combine everything into one library
file(GLOB_RECURSE SOURCES "*.cu" "*.cpp")
list(FILTER SOURCES EXCLUDE REGEX ".*main.*")  # Exclude main files

add_library(biogpu ${SOURCES})

set_property(TARGET biogpu PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET biogpu PROPERTY CUDA_ARCHITECTURES 70 75 80)

target_compile_options(biogpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    $<$<COMPILE_LANGUAGE:CUDA>:--relocatable-device-code=true>
)

# Simple test executable
add_executable(test_simple test_simple.cu)
target_link_libraries(test_simple biogpu)