#!/bin/bash

# Build script for the resistance detection pipeline with allele frequency tracking
# Run this from the runtime/kernels/resistance directory

echo "Building Resistance Detection Pipeline with Allele Frequency Tracking..."
echo "=========================================================="

# Create build directory if it doesn't exist
if [ ! -d "build" ]; then
    mkdir build
fi

cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES="70;75;80;86" \
    -DCMAKE_CXX_STANDARD=17

# Check if CMake configuration was successful
if [ $? -ne 0 ]; then
    echo "CMake configuration failed!"
    exit 1
fi

# Build the clean_resistance_pipeline target
echo ""
echo "Building clean_resistance_pipeline..."
make clean_resistance_pipeline -j$(nproc)

# Check if build was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "Build successful!"
    echo "Executable: build/clean_resistance_pipeline"
    echo ""
    echo "Usage:"
    echo "  ./build/clean_resistance_pipeline <nucleotide_index> <protein_db> <reads_r1.fastq.gz> <reads_r2.fastq.gz> [output_prefix] [fq_resistance_csv] [--no-bloom] [--no-sw]"
    echo ""
    echo "Example:"
    echo "  ./build/clean_resistance_pipeline \\"
    echo "    /path/to/nucleotide/index \\"
    echo "    /path/to/protein/db \\"
    echo "    sample_R1.fastq.gz \\"
    echo "    sample_R2.fastq.gz \\"
    echo "    results/sample_output"
    echo ""
    echo "Output files will include:"
    echo "  - results/sample_output.json                          # Main results"
    echo "  - results/sample_output_protein_matches.csv           # All protein matches"
    echo "  - results/sample_output_allele_frequencies.csv        # Allele frequency data"
    echo "  - results/sample_output_clinical_fq_report.html       # Clinical HTML report with allele frequencies"
    echo "  - results/sample_output_clinical_fq_report.json       # Clinical JSON report"
    echo "  - results/sample_output_clinical_fq_report.txt        # Clinical text report"
    echo "  - results/sample_output.h5                            # HDF5 alignment data"
else
    echo ""
    echo "Build failed! Check the error messages above."
    exit 1
fi