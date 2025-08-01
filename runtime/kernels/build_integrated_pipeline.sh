#!/bin/bash
# Build script for integrated BioGPU Unified Pipeline

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}BioGPU Integrated Pipeline Build${NC}"
echo -e "${GREEN}========================================${NC}"

# Check prerequisites
echo -e "${YELLOW}Checking prerequisites...${NC}"

# Check for required files
REQUIRED_FILES=(
    "unified_amr_pipeline_integrated.cpp"
    "resistance/CleanResistancePipeline.h"
    "resistance/CleanResistancePipelineWrapper.cpp"
    "genes/amr_detection_pipeline.h"
    "shared/sample_csv_parser.cpp"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "${RED}Error: Required file not found: $file${NC}"
        exit 1
    fi
done

# Apply memory fixes if not already applied
if [ -f "genes/amr_detection_pipeline_memory_fix.patch" ]; then
    echo -e "${YELLOW}Checking if memory fixes need to be applied...${NC}"
    if ! grep -q "database_loaded" genes/amr_detection_pipeline.cpp 2>/dev/null; then
        echo "Applying memory fixes..."
        patch -p0 < genes/amr_detection_pipeline_memory_fix.patch || true
    else
        echo "Memory fixes already applied"
    fi
fi

# Create build directory
BUILD_DIR="build_integrated"
if [ -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}Cleaning previous build...${NC}"
    rm -rf "$BUILD_DIR"
fi

echo -e "${GREEN}Creating build directory...${NC}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Copy CMakeLists.txt
cp ../CMakeLists_integrated.txt CMakeLists.txt

# Configure with CMake
echo -e "${GREEN}Configuring with CMake...${NC}"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_ARCHITECTURES="75;80;86" \
      ..

# Build with detailed output
echo -e "${GREEN}Building integrated pipeline...${NC}"
make -j$(nproc) VERBOSE=1 2>&1 | tee build.log

# Check if build succeeded
if [ -f "bio_gpu_pipeline" ]; then
    echo -e "${GREEN}Build successful!${NC}"
    echo ""
    echo -e "${YELLOW}Executable:${NC} $(pwd)/bio_gpu_pipeline"
    echo -e "${YELLOW}Size:${NC} $(du -h bio_gpu_pipeline | cut -f1)"
    
    # Test the executable
    echo ""
    echo -e "${GREEN}Testing executable...${NC}"
    ./bio_gpu_pipeline --help || true
    
    echo ""
    echo -e "${GREEN}Installation:${NC}"
    echo "To install system-wide:"
    echo "  sudo make install"
    echo ""
    echo -e "${GREEN}Example usage:${NC}"
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
    echo "Check build.log for errors"
    exit 1
fi