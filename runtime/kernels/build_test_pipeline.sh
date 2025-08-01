#!/bin/bash
# Build script for testing the unified pipeline framework

set -e

echo "Building BioGPU Unified Pipeline Test..."

# Create a temporary directory for the test build
TEST_DIR="unified_test_build"
if [ -d "$TEST_DIR" ]; then
    rm -rf "$TEST_DIR"
fi
mkdir -p "$TEST_DIR"

# Copy necessary files
cp CMakeLists_unified_simple.txt "$TEST_DIR/CMakeLists.txt"
cp test_unified_pipeline.cpp "$TEST_DIR/"
cp unified_amr_pipeline.cpp "$TEST_DIR/"
cp -r shared "$TEST_DIR/"
ln -s "$(pwd)/../common" "$TEST_DIR/../common"

# Build in the test directory
cd "$TEST_DIR"
mkdir build
cd build

# Configure
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Build
make -j$(nproc)

echo "Build complete!"
echo "Run: ./bio_gpu_pipeline_test"