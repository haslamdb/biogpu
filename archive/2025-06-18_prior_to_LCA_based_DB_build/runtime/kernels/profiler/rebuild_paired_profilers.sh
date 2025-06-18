#!/bin/bash
# Script to rebuild the paired-end profilers after CMakeLists.txt changes

echo "Rebuilding BioGPU Paired-End Profilers..."
echo "========================================="

# Clean build directory
if [ -d "build" ]; then
    echo "Cleaning existing build directory..."
    rm -rf build
fi

# Create new build directory
mkdir -p build
cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake ..

# Build the paired profilers
echo "Building paired-end profilers..."
make -j$(nproc) gpu_paired_profiler streaming_paired_profiler adaptive_paired_profiler

# Build debug versions if requested
if [ "$1" == "--debug" ]; then
    echo "Building debug versions..."
    make -j$(nproc) gpu_paired_profiler_debug adaptive_paired_profiler_debug
fi

# Check if builds were successful
echo ""
echo "Build Results:"
echo "=============="
if [ -f "gpu_paired_profiler" ]; then
    echo "✓ gpu_paired_profiler built successfully"
else
    echo "✗ gpu_paired_profiler build failed"
fi

if [ -f "streaming_paired_profiler" ]; then
    echo "✓ streaming_paired_profiler built successfully"
else
    echo "✗ streaming_paired_profiler build failed"
fi

if [ -f "adaptive_paired_profiler" ]; then
    echo "✓ adaptive_paired_profiler built successfully"
else
    echo "✗ adaptive_paired_profiler build failed"
fi

echo ""
echo "To test the profilers:"
echo "====================="
echo "1. Direct GPU paired-end profiler:"
echo "   ./gpu_paired_profiler ../data/microbes.db ../data/test_R1.fastq ../data/test_R2.fastq output"
echo ""
echo "2. Streaming paired-end profiler (for large databases):"
echo "   ./streaming_paired_profiler ../data/large_microbes.db ../data/test_R1.fastq ../data/test_R2.fastq output"
echo ""
echo "3. Adaptive paired profiler (auto-selects best strategy):"
echo "   ./adaptive_paired_profiler ../data/microbes.db ../data/test_R1.fastq ../data/test_R2.fastq output"
echo ""
echo "For single-end reads, just provide one FASTQ file to the adaptive profiler."