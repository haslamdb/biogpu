#!/bin/bash
echo "=== Building Taxonomy Module Only ==="

# Clean previous builds
rm -f test_taxonomy_only taxonomy_only.o

echo "Compiling taxonomy module..."

# Compile taxonomy module (no conditional compilation needed)
nvcc -c taxonomy/taxonomy_processor.cu -o taxonomy_only.o \
     -std=c++17 --extended-lambda --expt-relaxed-constexpr \
     -I. -I../tools

if [ $? -eq 0 ]; then
    echo "✓ Taxonomy module compiled successfully"
    
    # Build test
    echo "Building taxonomy test..."
    nvcc test_taxonomy_only.cu taxonomy_only.o -o test_taxonomy_only \
         -std=c++17 --extended-lambda --expt-relaxed-constexpr \
         -I. -I../tools
    
    if [ $? -eq 0 ]; then
        echo "✓ Taxonomy test executable built successfully"
        echo "Running test..."
        ./test_taxonomy_only
    else
        echo "✗ Taxonomy test build failed"
        exit 1
    fi
else
    echo "✗ Taxonomy module compilation failed"
    exit 1
fi

echo "Taxonomy module build complete"