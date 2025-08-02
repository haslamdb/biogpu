#!/bin/bash
# Improved build script for resistance pipeline with proper dependency paths

set -e

echo "========================================="
echo "Building Resistance Detection Pipeline"
echo "========================================="

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Check for dependencies
echo -e "${YELLOW}Checking dependencies...${NC}"

HDF5_READY=true
JSON_READY=true

if [ ! -f "/usr/include/hdf5/serial/H5Cpp.h" ]; then
    echo -e "${RED}✗ HDF5 C++ headers not found${NC}"
    echo "  Please run: ../install_dependencies.sh"
    HDF5_READY=false
else
    echo -e "${GREEN}✓ HDF5 C++ headers found${NC}"
fi

if [ ! -f "/usr/include/jsoncpp/json/json.h" ]; then
    echo -e "${RED}✗ JSON C++ headers not found${NC}"
    echo "  Please run: ../install_dependencies.sh"
    JSON_READY=false
else
    echo -e "${GREEN}✓ JSON C++ headers found${NC}"
fi

if [ "$HDF5_READY" = false ] || [ "$JSON_READY" = false ]; then
    echo ""
    echo -e "${RED}Missing dependencies. Please install them first.${NC}"
    exit 1
fi

# Compiler flags
HDF5_CFLAGS="-I/usr/include/hdf5/serial"
HDF5_LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5"
JSON_CFLAGS="-I/usr/include/jsoncpp"
JSON_LDFLAGS="-ljsoncpp"

# Create build directory
mkdir -p build
cd build

echo ""
echo -e "${GREEN}Compiling CUDA kernels...${NC}"

# Compile CUDA kernels with proper architecture
nvcc -c ../fq_mutation_detector.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../enhanced_mutation_detection_unified.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../kmer_screening.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../translated_search_revised.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../../shared/bloom_filter.cu -O3 -arch=sm_75 -I../../shared -Wno-deprecated-gpu-targets

# Fix for resistance_detection_gpu.cu - compile as C++ to avoid syntax errors
echo ""
echo -e "${GREEN}Compiling resistance detection GPU (as C++)...${NC}"
nvcc -c ../resistance_detection_gpu.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets -x cu

echo ""
echo -e "${GREEN}Compiling C++ files...${NC}"

# Compile C++ files with HDF5 and JSON includes
g++ -c ../clean_resistance_pipeline_main.cpp -std=c++17 -O3 -I.. -I../../shared $HDF5_CFLAGS
g++ -c ../fq_mutation_reporter.cpp -std=c++17 -O3
g++ -c ../diagnostic_report.cpp -std=c++17 -O3
g++ -c ../global_fq_resistance_mapper.cpp -std=c++17 -O3
g++ -c ../hdf5_alignment_writer.cpp -std=c++17 -O3 $HDF5_CFLAGS
g++ -c ../resistance_detector.cpp -std=c++17 -O3 -I../../shared $JSON_CFLAGS
g++ -c ../clinical_fq_report_generator.cpp -std=c++17 -O3
g++ -c ../fixed_database_loader.cpp -std=c++17 -O3 $JSON_CFLAGS
g++ -c ../../shared/sample_csv_parser.cpp -std=c++17 -O3 -I../../shared

echo ""
echo -e "${GREEN}Linking executable...${NC}"

# Link with all dependencies
nvcc -o clean_resistance_pipeline \
    clean_resistance_pipeline_main.o \
    fq_mutation_detector.o \
    enhanced_mutation_detection_unified.o \
    kmer_screening.o \
    resistance_detection_gpu.o \
    translated_search_revised.o \
    bloom_filter.o \
    fq_mutation_reporter.o \
    diagnostic_report.o \
    global_fq_resistance_mapper.o \
    hdf5_alignment_writer.o \
    resistance_detector.o \
    clinical_fq_report_generator.o \
    fixed_database_loader.o \
    sample_csv_parser.o \
    -lcudart -lz -lpthread \
    $HDF5_LDFLAGS $JSON_LDFLAGS \
    -O3

if [ -f "clean_resistance_pipeline" ]; then
    echo ""
    echo -e "${GREEN}=========================================${NC}"
    echo -e "${GREEN}Build successful!${NC}"
    echo -e "${GREEN}=========================================${NC}"
    echo "Executable: $(pwd)/clean_resistance_pipeline"
    
    # Copy to parent directory for easy access
    cp clean_resistance_pipeline ..
    
    echo ""
    echo "Testing executable..."
    ./clean_resistance_pipeline --help || true
    
    echo ""
    echo "The resistance pipeline is ready to integrate with the unified pipeline!"
else
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi