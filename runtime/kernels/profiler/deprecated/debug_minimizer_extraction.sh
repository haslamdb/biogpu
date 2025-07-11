#!/bin/bash
# debug_minimizer_extraction.sh
# Script to help debug minimizer extraction issues

echo "=== Debugging Minimizer Extraction Pipeline ==="

# Check if test data exists
if [ ! -f "../../../../data/test_50_genomes.fna" ]; then
    echo "ERROR: Test data not found at ../../../../data/test_50_genomes.fna"
    echo "Please ensure the test data file exists"
    exit 1
fi

# Check CUDA installation
echo -e "\n1. Checking CUDA installation..."
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: nvcc not found. Please ensure CUDA is installed and in PATH"
    exit 1
fi

nvcc --version
echo ""

# Check GPU availability
echo "2. Checking GPU availability..."
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader || {
    echo "ERROR: nvidia-smi failed. Check NVIDIA driver installation"
    exit 1
}
echo ""

# Compile with debug flags
echo "3. Compiling with debug flags..."
make -f Makefile.streaming_minimizer clean
make -f Makefile.streaming_minimizer debug || {
    echo "ERROR: Compilation failed"
    exit 1
}
echo ""

# Create debug output directory
mkdir -p debug_output

# Run with cuda-memcheck for memory errors
echo "4. Running with cuda-memcheck..."
cuda-memcheck --leak-check full ./test_streaming_minimizer ../../../../data/test_50_genomes.fna 2 31 2>&1 | tee debug_output/memcheck.log

# Run with compute-sanitizer for race conditions
echo -e "\n5. Running with compute-sanitizer..."
compute-sanitizer --tool memcheck ./test_streaming_minimizer ../../../../data/test_50_genomes.fna 2 31 2>&1 | tee debug_output/sanitizer.log

# Run with nvprof for performance analysis (if available)
if command -v nvprof &> /dev/null; then
    echo -e "\n6. Running with nvprof..."
    nvprof --print-gpu-trace ./test_streaming_minimizer ../../../../data/test_50_genomes.fna 2 31 2>&1 | tee debug_output/nvprof.log
else
    echo -e "\n6. nvprof not available, skipping performance profiling"
fi

# Analyze the first few sequences from the test file
echo -e "\n7. Analyzing test data..."
echo "First 5 sequences in test file:"
awk '/^>/ {if (seq) print length(seq); seq=""; print; next} {seq=seq$0} END {if (seq) print length(seq)}' ../../../../data/test_50_genomes.fna | head -10

# Check for common issues in logs
echo -e "\n8. Checking for common issues..."

if grep -q "out of memory" debug_output/*.log 2>/dev/null; then
    echo "WARNING: Out of memory errors detected"
fi

if grep -q "invalid device" debug_output/*.log 2>/dev/null; then
    echo "WARNING: Invalid device errors detected"
fi

if grep -q "unspecified launch failure" debug_output/*.log 2>/dev/null; then
    echo "WARNING: Kernel launch failures detected"
fi

echo -e "\nDebug logs saved to debug_output/"
echo "Review the logs for detailed error information"
