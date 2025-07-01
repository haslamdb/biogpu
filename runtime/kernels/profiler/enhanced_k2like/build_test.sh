#!/bin/bash
# build_test.sh - Build script for modular BioGPU tests

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
CUDA_ARCH="sm_75"  # Adjust for your GPU
BUILD_DIR="build_test"
TEST_NAME="biogpu_modular_test"

# Source directories - adjust these paths
BIOGPU_ROOT="/home/david/Documents/Code/biogpu"
INCLUDE_DIRS=(
    "-I${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like"
    "-I${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like/gpu"
    "-I${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like/features"
    "-I${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like/processing"
    "-I/usr/local/cuda/include"
)

# Source files - add new source files here
SOURCE_FILES=(
    "feature_tests.cpp"
    "${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like/gpu/gpu_database_kernels.cu"
    "${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like/processing/minimizer_feature_extractor.cu"
    "${BIOGPU_ROOT}/runtime/kernels/profiler/enhanced_k2like/processing/genome_file_processor.cpp"
    # Add more source files as needed
)

# Compiler flags
NVCC_FLAGS=(
    "-arch=${CUDA_ARCH}"
    "-std=c++17"
    "-O3"
    "-lineinfo"
    "--extended-lambda"
    "-use_fast_math"
    "--expt-relaxed-constexpr"
)

# Linker flags
LDFLAGS=(
    "-lcudart"
    "-lcuda"
    "-lstdc++fs"
)

# Debug mode
if [[ "$1" == "debug" ]]; then
    echo -e "${YELLOW}Building in DEBUG mode${NC}"
    NVCC_FLAGS+=("-g" "-G" "-DDEBUG")
else
    echo -e "${GREEN}Building in RELEASE mode${NC}"
    NVCC_FLAGS+=("-DNDEBUG")
fi

# Create build directory
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# Function to check if file exists
check_file() {
    if [[ ! -f "$1" ]]; then
        echo -e "${RED}Error: File not found: $1${NC}"
        exit 1
    fi
}

# Check all source files exist
echo "Checking source files..."
for src in "${SOURCE_FILES[@]}"; do
    check_file "$src"
done

# Build command
echo -e "${YELLOW}Building ${TEST_NAME}...${NC}"
BUILD_CMD="nvcc ${NVCC_FLAGS[@]} ${INCLUDE_DIRS[@]} -o ${TEST_NAME} ${SOURCE_FILES[@]} ${LDFLAGS[@]}"

echo "Build command:"
echo "$BUILD_CMD"
echo

# Execute build
if $BUILD_CMD 2>&1 | tee build.log; then
    echo -e "${GREEN}Build successful!${NC}"
    echo -e "Executable: ${BUILD_DIR}/${TEST_NAME}"
    
    # Print usage
    echo
    echo "Run with:"
    echo "  ./${BUILD_DIR}/${TEST_NAME}"
    echo
    echo "Options:"
    echo "  --fna FILE              Input FNA file"
    echo "  --batch-size N          Process N genomes per batch"
    echo "  --no-uniqueness         Disable uniqueness scoring"
    echo "  --no-cooccurrence       Disable co-occurrence analysis"
    echo "  --examples N            Print N example minimizers"
else
    echo -e "${RED}Build failed! Check build.log for errors${NC}"
    exit 1
fi

# Optional: Run test immediately
if [[ "$2" == "run" ]]; then
    echo
    echo -e "${YELLOW}Running test...${NC}"
    ./${TEST_NAME} "${@:3}"
fi

