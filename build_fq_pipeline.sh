#!/bin/bash

# Build script for fluoroquinolone resistance detection pipeline
# This builds the CUDA-accelerated two-stage mutation detection system

set -e  # Exit on error

echo "Building Fluoroquinolone Resistance Detection Pipeline..."
echo "========================================================="

# Create build directory
if [ ! -d "build" ]; then
    mkdir build
fi

cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES="61;70;75;80;86" \
    -DBUILD_TESTS=ON

# Build the pipeline
echo "Building pipeline..."
make -j$(nproc) fq_pipeline_gpu

# Build the kernel library if not already built
if [ ! -f "libbiogpu_kernels.so" ]; then
    echo "Building CUDA kernels..."
    make -j$(nproc) biogpu_kernels
fi

# Check if build succeeded
if [ -f "fq_pipeline_gpu" ]; then
    echo "Build successful!"
    echo ""
    echo "Executable created: build/fq_pipeline_gpu"
    echo ""
    echo "Usage:"
    echo "  ./build/fq_pipeline_gpu <index_path> <reads_R1.fastq.gz> <reads_R2.fastq.gz> <output.json>"
    echo ""
    echo "Example:"
    echo "  ./build/fq_pipeline_gpu data/indices/fq_mutations /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz results.json"
else
    echo "Build failed!"
    exit 1
fi

# Optional: Build and run tests
echo ""
echo "Would you like to build and run tests? (y/n)"
read -r response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo "Building tests..."
    make -j$(nproc) test_cuda_kernels
    
    echo "Running CUDA kernel tests..."
    ./test_cuda_kernels
fi

cd ..
echo ""
echo "Build complete!"