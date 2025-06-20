#!/bin/bash
# Simpler compile script for heap debugging without enhanced features

echo "=== COMPILING WITH HEAP CORRUPTION DETECTIVE (SIMPLE VERSION) ==="

# Clean up
rm -f *.o debug_kraken_heap_simple 2>/dev/null

# Step 1: Compile heap corruption detective
echo "Step 1: Compiling heap corruption detective..."
g++ -std=c++17 -g -O0 -fPIC -I. \
    -c heap_corruption_detective.cpp -o heap_detective.o

if [ $? -ne 0 ]; then
    echo "❌ Heap corruption detective compilation failed"
    exit 1
fi
echo "✓ heap_detective.o created"

# Step 2: Compile all CUDA files together with the heap detective
echo ""
echo "Step 2: Compiling all files together..."
nvcc -std=c++17 -g -G -O0 --extended-lambda \
     kraken_pipeline_main.cu \
     gpu_kraken_classifier.cu \
     gpu_kraken_database_builder.cu \
     gpu_minimizer_extraction.cu \
     phase1_enhanced_classification_cc.cu \
     sample_csv_parser.cpp \
     heap_corruption_detective.cpp \
     -I. -I.. -I../../../include \
     -o debug_kraken_heap_simple \
     -lcuda -lcudart -lz

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 SUCCESS! debug_kraken_heap_simple created"
    ls -la debug_kraken_heap_simple
    
    echo ""
    echo "Now you can run with heap debugging enabled!"
    echo "Example:"
    echo "  ./debug_kraken_heap_simple build --genome-dir /path/to/genomes --output /path/to/output"
else
    echo ""
    echo "❌ Compilation failed"
    exit 1
fi

# Make executable
chmod +x debug_kraken_heap_simple

echo ""
echo "=== COMPILATION COMPLETE ==="
echo "Executable: debug_kraken_heap_simple"
echo ""
echo "The heap corruption detective is now integrated!"
echo "It will automatically enable debugging when the program starts."