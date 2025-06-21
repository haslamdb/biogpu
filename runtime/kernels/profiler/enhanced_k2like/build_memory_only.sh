#!/bin/bash

# Direct build of memory manager module only using nvcc

echo "Building GPU Memory Manager module only..."

# Create output directory
mkdir -p build_test

# Compile memory manager
nvcc -c memory/gpu_memory_manager.cu \
    -o build_test/gpu_memory_manager.o \
    -I. \
    --extended-lambda \
    --expt-relaxed-constexpr \
    -arch=sm_70 \
    -std=c++17 \
    2>&1 | tee build_test/memory_build.log

# Check if compilation succeeded
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "✓ gpu_memory_manager.cu compiled successfully!"
    
    # Try to build test executable
    echo "Building test executable..."
    nvcc test_memory_build.cu build_test/gpu_memory_manager.o \
        -o build_test/test_memory \
        -I. \
        --extended-lambda \
        --expt-relaxed-constexpr \
        -arch=sm_70 \
        -std=c++17 \
        2>&1 | tee -a build_test/memory_build.log
        
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "✓ Test executable built successfully!"
        echo "You can run: ./build_test/test_memory"
    else
        echo "✗ Test executable build failed!"
        echo "Check build_test/memory_build.log for errors"
    fi
else
    echo "✗ gpu_memory_manager.cu compilation failed!"
    echo "Check build_test/memory_build.log for errors"
    
    # Show first few errors
    echo -e "\nFirst errors:"
    grep -m 5 "error:" build_test/memory_build.log
fi