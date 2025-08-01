#!/bin/bash
# Build script for BioGPU Unified AMR Pipeline
# This script builds the unified pipeline that combines resistance and AMR gene detection

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}BioGPU Unified Pipeline Build Script${NC}"
echo -e "${GREEN}========================================${NC}"

# Check if we're in the right directory
if [ ! -f "unified_amr_pipeline.cpp" ]; then
    echo -e "${RED}Error: unified_amr_pipeline.cpp not found!${NC}"
    echo "Please run this script from the runtime/kernels directory"
    exit 1
fi

# Check for CUDA
if ! command -v nvcc &> /dev/null; then
    echo -e "${RED}Error: CUDA compiler (nvcc) not found!${NC}"
    echo "Please ensure CUDA toolkit is installed and in PATH"
    exit 1
fi

# Display CUDA info
echo -e "${YELLOW}CUDA Information:${NC}"
nvcc --version | head -n 4
echo ""

# Check for available GPUs
echo -e "${YELLOW}Available GPUs:${NC}"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader || echo "nvidia-smi not available"
echo ""

# Create build directory
BUILD_DIR="build_unified"
if [ -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}Removing existing build directory...${NC}"
    rm -rf "$BUILD_DIR"
fi

echo -e "${GREEN}Creating build directory...${NC}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo -e "${GREEN}Configuring with CMake...${NC}"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_ARCHITECTURES="75;80;86" \
      -DCMAKE_CXX_FLAGS="-march=native" \
      -DCMAKE_CUDA_FLAGS="--use_fast_math -O3" \
      -DCMAKE_INSTALL_PREFIX=/usr/local \
      -f ../CMakeLists_unified.txt \
      ..

# Build
echo -e "${GREEN}Building unified pipeline...${NC}"
make -j$(nproc) VERBOSE=1

# Check if build succeeded
if [ -f "bio_gpu_pipeline" ]; then
    echo -e "${GREEN}Build successful!${NC}"
    echo ""
    echo -e "${YELLOW}Executable location:${NC} $(pwd)/bio_gpu_pipeline"
    echo -e "${YELLOW}Size:${NC} $(du -h bio_gpu_pipeline | cut -f1)"
    echo ""
    
    # Test the executable
    echo -e "${GREEN}Testing executable...${NC}"
    ./bio_gpu_pipeline --help || true
    
    echo ""
    echo -e "${GREEN}Installation:${NC}"
    echo "To install system-wide to /usr/local/bin:"
    echo "  sudo make install"
    echo ""
    echo -e "${GREEN}Usage example:${NC}"
    echo "./bio_gpu_pipeline --use-multi-gpu \\"
    echo "  --r1 /data/sample_R1.fastq.gz \\"
    echo "  --r2 /data/sample_R2.fastq.gz \\"
    echo "  --output-dir /data/results \\"
    echo "  --reference-db /data/microbial_ref_db \\"
    echo "  --resistance-db /data/quinolone_resistance_db \\"
    echo "  --sample-id SAMPLE001 \\"
    echo "  --progress-json"
    
else
    echo -e "${RED}Build failed!${NC}"
    echo "Check the error messages above for details"
    exit 1
fi

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Build completed successfully!${NC}"
echo -e "${GREEN}========================================${NC}"