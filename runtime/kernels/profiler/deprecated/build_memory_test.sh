#!/bin/bash

# Build script for testing gpu_memory_manager module only

echo "Building GPU Memory Manager module test..."

# Create build directory
mkdir -p build_memory_test
cd build_memory_test

# Configure with the test CMakeLists
cmake -DCMAKE_BUILD_TYPE=Debug -C ../CMakeLists_memory_test.txt ..

# Build
make -j$(nproc)

# Check if build succeeded
if [ $? -eq 0 ]; then
    echo "✓ Build succeeded!"
    echo "You can run: ./test_memory_build"
else
    echo "✗ Build failed!"
    exit 1
fi