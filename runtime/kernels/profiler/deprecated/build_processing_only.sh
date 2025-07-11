#!/bin/bash
echo "=== Building Processing Module Only ==="

# Clean previous builds
rm -f test_processing_only processing_only.o

echo "Compiling processing module..."

# Compile processing module
nvcc -c processing/genome_file_processor.cu -o processing_only.o \
     -std=c++17 --extended-lambda --expt-relaxed-constexpr \
     -I.

if [ $? -eq 0 ]; then
    echo "✓ Processing module compiled successfully"
    
    # Build test
    echo "Building processing test..."
    nvcc test_processing_only.cu processing_only.o -o test_processing_only \
         -std=c++17 --extended-lambda --expt-relaxed-constexpr \
         -I.
    
    if [ $? -eq 0 ]; then
        echo "✓ Processing test executable built successfully"
        echo "Running test..."
        ./test_processing_only
    else
        echo "✗ Processing test build failed"
        exit 1
    fi
else
    echo "✗ Processing module compilation failed"
    exit 1
fi

echo "Processing module build complete"