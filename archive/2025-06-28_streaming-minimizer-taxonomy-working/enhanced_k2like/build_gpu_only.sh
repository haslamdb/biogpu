#!/bin/bash
echo "Building GPU module only..."

# Compile just the GPU module
nvcc -c gpu/gpu_database_kernels.cu -o gpu_kernels.o \
     -std=c++17 --extended-lambda --expt-relaxed-constexpr

if [ $? -eq 0 ]; then
    echo "✓ GPU module compiled successfully"
else
    echo "✗ GPU module compilation failed"
    exit 1
fi

echo "GPU module build complete"