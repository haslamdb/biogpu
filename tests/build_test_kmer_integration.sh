#!/bin/bash

# Build script for integrated FQ resistance detection
set -e

echo "=== Building Integrated FQ Resistance Pipeline ==="

# Configuration
CUDA_ARCH="sm_61"  # For Titan Xp, adjust as needed
BUILD_DIR="build_integrated"
SRC_DIR="runtime/kernels/resistance"

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
nvcc --version | head -n 4

echo ""
echo "=== Compiling CUDA kernels ==="

# Compile the enhanced k-mer screening kernel
echo "Compiling fixed_kmer_screening.cu..."
nvcc -c ${SRC_DIR}/fixed_kmer_screening.cu \
    -o ${BUILD_DIR}/fixed_kmer_screening.o \
    -arch=${CUDA_ARCH} \
    -I${SRC_DIR} \
    -O3 \
    -Xptxas -v \
    --compiler-options -fPIC

# Compile the main mutation detector
echo "Compiling fq_mutation_detector.cu..."
nvcc -c ${SRC_DIR}/fq_mutation_detector.cu \
    -o ${BUILD_DIR}/fq_mutation_detector.o \
    -arch=${CUDA_ARCH} \
    -I${SRC_DIR} \
    -O3 \
    -Xptxas -v \
    --compiler-options -fPIC

# Compile the host pipeline code
echo "Compiling fq_pipeline_host.cpp..."
nvcc -c ${SRC_DIR}/fq_pipeline_host.cpp \
    -o ${BUILD_DIR}/fq_pipeline_host.o \
    -arch=${CUDA_ARCH} \
    -I${SRC_DIR} \
    -O3 \
    --compiler-options -fPIC

echo ""
echo "=== Linking pipeline executable ==="

# Link the main pipeline executable
nvcc ${BUILD_DIR}/fixed_kmer_screening.o \
     ${BUILD_DIR}/fq_mutation_detector.o \
     ${BUILD_DIR}/fq_pipeline_host.o \
     -o ${BUILD_DIR}/fq_pipeline_gpu \
     -arch=${CUDA_ARCH} \
     -lcuda \
     -lcudart \
     -lz \
     -O3

echo ""
echo "=== Compiling integration test ==="

# Compile the integration test
nvcc ${SRC_DIR}/test_kmer_integration.cpp \
     ${BUILD_DIR}/fixed_kmer_screening.o \
     ${BUILD_DIR}/fq_mutation_detector.o \
     -o ${BUILD_DIR}/test_integration \
     -arch=${CUDA_ARCH} \
     -I${SRC_DIR} \
     -lcuda \
     -lcudart \
     -O3

echo ""
echo "=== Build Summary ==="
echo "Main pipeline: ${BUILD_DIR}/fq_pipeline_gpu"
echo "Integration test: ${BUILD_DIR}/test_integration"

# Check if executables were created
if [ -f "${BUILD_DIR}/fq_pipeline_gpu" ]; then
    echo "✅ Main pipeline built successfully"
    ls -lh ${BUILD_DIR}/fq_pipeline_gpu
else
    echo "❌ Main pipeline build failed"
    exit 1
fi

if [ -f "${BUILD_DIR}/test_integration" ]; then
    echo "✅ Integration test built successfully"
    ls -lh ${BUILD_DIR}/test_integration
else
    echo "❌ Integration test build failed"
    exit 1
fi

echo ""
echo "=== Quick Test ==="
echo "Running integration test with default parameters..."
mkdir -p /tmp/test_index

if ${BUILD_DIR}/test_integration /tmp/test_index; then
    echo "✅ Integration test passed!"
else
    echo "❌ Integration test failed"
    exit 1
fi

echo ""
echo "=== Build Complete ==="
echo "To run the pipeline:"
echo "  ${BUILD_DIR}/fq_pipeline_gpu <index_dir> <reads_R1.fastq.gz> <reads_R2.fastq.gz> <output.json>"
echo ""
echo "To test with your own index:"
echo "  ${BUILD_DIR}/test_integration <your_index_dir>"