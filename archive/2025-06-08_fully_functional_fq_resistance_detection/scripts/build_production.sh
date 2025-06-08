#!/bin/bash

# Production build script for reorganized BioGPU project
set -e

echo "=== Building BioGPU Production Version ==="

# Configuration
CUDA_ARCH="sm_61"  # For Titan Xp, adjust as needed
BUILD_DIR="build/production"
SRC_DIR="src"

# Create build directory
mkdir -p ${BUILD_DIR}

echo "Using CUDA architecture: ${CUDA_ARCH}"
echo "Build directory: ${BUILD_DIR}"
echo "Source directory: ${SRC_DIR}"

# Check for CUDA
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: nvcc not found. Please install CUDA toolkit."
    exit 1
fi

# Check CUDA version
echo "CUDA version:"
nvcc --version | head -n 4

echo ""
echo "=== Compiling CUDA kernels ==="

# Compile k-mer screening kernel
echo "Compiling kmer_screening.cu..."
nvcc -c ${SRC_DIR}/kernels/resistance/kmer_screening.cu \
    -o ${BUILD_DIR}/kmer_screening.o \
    -arch=${CUDA_ARCH} \
    -I${SRC_DIR}/kernels/resistance \
    -I${SRC_DIR}/common \
    -O3 \
    -Xptxas -v \
    --compiler-options -fPIC

# Compile main mutation detector
echo "Compiling fq_mutation_detector.cu..."
nvcc -c ${SRC_DIR}/kernels/resistance/fq_mutation_detector.cu \
    -o ${BUILD_DIR}/fq_mutation_detector.o \
    -arch=${CUDA_ARCH} \
    -I${SRC_DIR}/kernels/resistance \
    -I${SRC_DIR}/common \
    -O3 \
    -Xptxas -v \
    --compiler-options -fPIC

# Compile host pipeline code
echo "Compiling pipeline host code..."
nvcc -c ${SRC_DIR}/host/fq_pipeline_host.cpp \
    -o ${BUILD_DIR}/fq_pipeline_host.o \
    -arch=${CUDA_ARCH} \
    -I${SRC_DIR}/kernels/resistance \
    -I${SRC_DIR}/common \
    -O3 \
    --compiler-options -fPIC

echo ""
echo "=== Linking production executable ==="

# Link the main pipeline executable
nvcc ${BUILD_DIR}/kmer_screening.o \
     ${BUILD_DIR}/fq_mutation_detector.o \
     ${BUILD_DIR}/fq_pipeline_host.o \
     -o ${BUILD_DIR}/fq_pipeline_gpu \
     -arch=${CUDA_ARCH} \
     -lcuda \
     -lcudart \
     -lz \
     -O3

echo ""
echo "=== Build Summary ==="
echo "Production pipeline: ${BUILD_DIR}/fq_pipeline_gpu"

# Check if executable was created
if [ -f "${BUILD_DIR}/fq_pipeline_gpu" ]; then
    echo "✅ Production pipeline built successfully"
    ls -lh ${BUILD_DIR}/fq_pipeline_gpu
    echo ""
    echo "Test the build:"
    echo "  ${BUILD_DIR}/fq_pipeline_gpu --version"
    echo ""
    echo "Run with your data:"
    echo "  ${BUILD_DIR}/fq_pipeline_gpu /data/fq_resistance_index \\"
    echo "    data/test_fastq/VRE12_R1.fastq.gz \\"
    echo "    data/test_fastq/VRE12_R2.fastq.gz \\"
    echo "    VRE12_results.json"
else
    echo "❌ Production pipeline build failed"
    echo "Check the compilation errors above"
    exit 1
fi

echo ""
echo "=== Production Build Complete ==="