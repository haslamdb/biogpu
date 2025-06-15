#!/bin/bash
# build_and_test.sh
# Build and test the clean Kraken2-style GPU implementation

set -e  # Exit on error

echo "🔧 Building Clean Kraken2-Style GPU Pipeline"
echo "==========================================="

# Configuration
GENOME_DIR="${1:-data/pathogen_profiler_db/genomes}"
OUTPUT_DIR="${2:-./clean_kraken_db}"
READS_FILE="${3:-sample_reads.fastq}"

# Check for CUDA
echo "Checking CUDA installation..."
if ! command -v nvcc &> /dev/null; then
    echo "❌ CUDA not found. Please install CUDA toolkit."
    exit 1
fi

CUDA_VERSION=$(nvcc --version | grep "release" | sed 's/.*release \([0-9.]*\).*/\1/')
echo "✓ Found CUDA version: $CUDA_VERSION"

# Detect GPU architecture
echo "Detecting GPU architecture..."
if command -v nvidia-smi &> /dev/null; then
    GPU_INFO=$(nvidia-smi --query-gpu=name --format=csv,noheader,nounits | head -n1)
    echo "✓ Found GPU: $GPU_INFO"
    
    # Set appropriate architecture
    if [[ $GPU_INFO == *"RTX 40"* ]] || [[ $GPU_INFO == *"RTX 30"* ]]; then
        CUDA_ARCH="sm_86"
    elif [[ $GPU_INFO == *"RTX 20"* ]] || [[ $GPU_INFO == *"Tesla V100"* ]]; then
        CUDA_ARCH="sm_70"
    elif [[ $GPU_INFO == *"GTX 10"* ]]; then
        CUDA_ARCH="sm_61"
    else
        CUDA_ARCH="sm_70"  # Default
    fi
    echo "Using CUDA architecture: $CUDA_ARCH"
else
    echo "⚠️  nvidia-smi not found, using default architecture sm_70"
    CUDA_ARCH="sm_70"
fi

# Create build directory
mkdir -p build_clean
cd build_clean

# Copy source file
cp ../kraken2_gpu.cu .

echo ""
echo "📦 Compiling Kraken2 GPU Pipeline..."
echo "======================================"

# Compile with optimizations
nvcc -std=c++17 -O3 -arch=$CUDA_ARCH \
     -lineinfo \
     --compiler-options -Wall \
     kraken2_gpu.cu -o kraken2_gpu

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful!"
else
    echo "❌ Compilation failed!"
    exit 1
fi

echo ""
echo "🧪 Testing the Pipeline..."
echo "========================="

# Check if genome directory exists
if [ ! -d "../$GENOME_DIR" ]; then
    echo "⚠️  Genome directory not found: $GENOME_DIR"
    echo "Creating test genome for demonstration..."
    
    mkdir -p test_genomes
    cat > test_genomes/test_genome.fna << 'EOF'
>test_sequence_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>test_sequence_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
EOF
    GENOME_DIR="test_genomes"
fi

# Test database building
echo "Building database from: $GENOME_DIR"
./kraken2_gpu build $GENOME_DIR $OUTPUT_DIR

if [ $? -eq 0 ]; then
    echo "✅ Database build successful!"
else
    echo "❌ Database build failed!"
    exit 1
fi

# Check if reads file exists
if [ ! -f "../$READS_FILE" ]; then
    echo "⚠️  Reads file not found: $READS_FILE"
    echo "Creating test reads for demonstration..."
    
    cat > test_reads.fastq << 'EOF'
@read_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_3
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_4
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    READS_FILE="test_reads.fastq"
fi

# Test classification
echo "Classifying reads from: $READS_FILE"
./kraken2_gpu classify $OUTPUT_DIR $READS_FILE > classification_results.txt

if [ $? -eq 0 ]; then
    echo "✅ Classification successful!"
    echo ""
    echo "📊 Results Summary:"
    echo "=================="
    
    # Show results
    echo "Classification results:"
    cat classification_results.txt | tail -n 20
    
    echo ""
    echo "Database info:"
    if [ -f "$OUTPUT_DIR/taxonomy.tsv" ]; then
        echo "Taxa in database: $(wc -l < $OUTPUT_DIR/taxonomy.tsv)"
    fi
    
    if [ -f "$OUTPUT_DIR/hash_table.k2d" ]; then
        echo "Database size: $(ls -lh $OUTPUT_DIR/hash_table.k2d | awk '{print $5}')"
    fi
    
else
    echo "❌ Classification failed!"
    exit 1
fi

echo ""
echo "🎉 All tests passed!"
echo "==================="
echo ""
echo "Usage:"
echo "  Build database: ./kraken2_gpu build <genome_dir> <output_dir>"
echo "  Classify reads: ./kraken2_gpu classify <database_dir> <reads_file>"
echo ""
echo "Performance tips:"
echo "  - Use SSD storage for better I/O performance"
echo "  - Increase GPU batch size for larger datasets"
echo "  - For production use, implement proper taxonomy loading"
echo "  - Consider using larger k-mer sizes (35-51) for better specificity"

cd ..
echo "Build completed in: build_clean/"