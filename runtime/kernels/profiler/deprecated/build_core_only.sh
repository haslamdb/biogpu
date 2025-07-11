#!/bin/bash
echo "=== Building Core Module Only ==="

# Clean previous builds
rm -f test_core_only core_only.o

echo "Compiling core module..."

# Compile core module only
nvcc -c core/gpu_database_builder_core.cu -o core_only.o \
     -std=c++17 --extended-lambda --expt-relaxed-constexpr \
     -I. -I./memory -I./gpu

if [ $? -eq 0 ]; then
    echo "✓ Core module compiled successfully"
    
    # Try linking with test
    echo "Building core test..."
    nvcc test_core_only.cu core_only.o memory/*.cu gpu/*.cu -o test_core_only \
         -std=c++17 --extended-lambda --expt-relaxed-constexpr -I.
    
    if [ $? -eq 0 ]; then
        echo "✓ Core test executable built successfully"
        echo "Running test..."
        ./test_core_only
    else
        echo "✗ Core test build failed"
        exit 1
    fi
else
    echo "✗ Core module compilation failed"
    exit 1
fi

echo "Core module build complete"