#!/bin/bash
# Build script for resistance pipeline executable

echo "Building resistance detection pipeline..."

# Create build directory
mkdir -p build
cd build

# Compile all necessary CUDA kernels and C++ files
echo "Compiling CUDA kernels..."
nvcc -c ../fq_mutation_detector.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../enhanced_mutation_detection_unified.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../kmer_screening.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../resistance_detection_gpu.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../translated_search_revised.cu -O3 -arch=sm_75 -Wno-deprecated-gpu-targets
nvcc -c ../../shared/bloom_filter.cu -O3 -arch=sm_75 -I../../shared -Wno-deprecated-gpu-targets

echo "Compiling C++ files..."
g++ -c ../clean_resistance_pipeline_main.cpp -std=c++17 -O3 -I.. -I../../shared
g++ -c ../fq_mutation_reporter.cpp -std=c++17 -O3
g++ -c ../diagnostic_report.cpp -std=c++17 -O3
g++ -c ../global_fq_resistance_mapper.cpp -std=c++17 -O3
g++ -c ../hdf5_alignment_writer.cpp -std=c++17 -O3
g++ -c ../resistance_detector.cpp -std=c++17 -O3
g++ -c ../clinical_fq_report_generator.cpp -std=c++17 -O3
g++ -c ../fixed_database_loader.cpp -std=c++17 -O3
g++ -c ../../shared/sample_csv_parser.cpp -std=c++17 -O3 -I../../shared

echo "Linking executable..."
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
    -lcudart -lz -lpthread -lhdf5 -O3

if [ -f "clean_resistance_pipeline" ]; then
    echo "Build successful!"
    echo "Executable: $(pwd)/clean_resistance_pipeline"
    # Copy to parent directory for easy access
    cp clean_resistance_pipeline ..
else
    echo "Build failed!"
    exit 1
fi