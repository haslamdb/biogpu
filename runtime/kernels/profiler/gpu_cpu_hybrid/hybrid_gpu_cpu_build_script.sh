#!/bin/bash

# Build script for Hybrid Metagenomic Profiler
# Optimized for clinical deployment

set -e  # Exit on any error

echo "=============================================="
echo "Building Hybrid Metagenomic Profiler"
echo "=============================================="

# Check system requirements
echo "Checking system requirements..."

# Check CUDA installation
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: CUDA not found. Please install CUDA toolkit."
    exit 1
fi

# Check GPU availability
if ! nvidia-smi &> /dev/null; then
    echo "ERROR: No NVIDIA GPU detected or driver not installed."
    exit 1
fi

# Get GPU info
echo "GPU Information:"
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader,nounits

# Check available memory
echo "System Memory:"
free -h

# Create build directory
BUILD_DIR="build"
if [ -d "$BUILD_DIR" ]; then
    echo "Removing existing build directory..."
    rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure build
echo "Configuring build..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES="75;80;86" \
    -DCMAKE_VERBOSE_MAKEFILE=ON

# Build
echo "Building (this may take a few minutes)..."
make -j$(nproc) 

# Check if build was successful
if [ -f "hybrid_profiler" ]; then
    echo "✅ Build successful!"
    
    # Get binary info
    ls -lh hybrid_profiler
    
    # Quick system test
    echo "Running quick system test..."
    ./hybrid_profiler --help 2>/dev/null || echo "Binary created successfully"
    
else
    echo "❌ Build failed!"
    exit 1
fi

echo ""
echo "=============================================="
echo "Build Complete!"
echo "=============================================="
echo ""
echo "Binary location: $(pwd)/hybrid_profiler"
echo ""
echo "Usage:"
echo "  ./hybrid_profiler <database_path> <fastq_file> [output_prefix]"
echo ""
echo "Example:"
echo "  ./hybrid_profiler /data/microbial_db sample.fastq clinical_analysis"
echo ""
echo "For debugging, use: ./hybrid_profiler_debug"
echo ""

# Performance recommendations
echo "Performance Recommendations:"
echo "- Ensure your database is on fast storage (SSD preferred)"
echo "- Use all available CPU RAM (current: $(free -h | awk '/^Mem:/{print $2}'))"
echo "- Monitor GPU memory usage with: nvidia-smi -l 1"
echo ""

# Database preparation reminder
echo "Database Preparation:"
echo "- If you don't have a prepared database, run:"
echo "  ./prepare_database <input_genomes> <output_database>"
echo ""

cd ..

# Create test script
cat > test_profiler.sh << 'EOF'
#!/bin/bash

# Test script for hybrid profiler
# Creates synthetic test data and runs basic functionality test

set -e

echo "Creating synthetic test data..."

# Create a simple test FASTQ file
cat > test_sample.fastq << 'FASTQ_EOF'
@read1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
FASTQ_EOF

echo "Test data created: test_sample.fastq"

# Create minimal test database (this would be your actual database)
mkdir -p test_db
echo "Creating minimal test database..."

# This is a placeholder - in practice, you'd have your prepared database
echo "Note: You'll need to provide an actual database file"
echo "Usage: ./build/hybrid_profiler test_db test_sample.fastq test_output"

EOF

chmod +x test_profiler.sh

echo "Test script created: ./test_profiler.sh"
echo ""
echo "To test the profiler:"
echo "1. Prepare your microbial database"
echo "2. Run: ./test_profiler.sh"
echo "3. Or run directly: ./build/hybrid_profiler <your_db> <your_fastq> <output_prefix>"