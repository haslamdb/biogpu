#!/bin/bash
# Compile script that includes heap corruption detective

echo "=== COMPILING WITH HEAP CORRUPTION DETECTIVE ==="

# Clean up
rm -f *.o debug_kraken_heap_debug 2>/dev/null

# Step 1: Compile heap corruption detective
echo "Step 1: Compiling heap corruption detective..."
g++ -std=c++17 -g -O0 -fPIC -I. \
    -c heap_corruption_detective.cpp -o heap_detective.o

if [ $? -ne 0 ]; then
    echo "❌ Heap corruption detective compilation failed"
    exit 1
fi
echo "✓ heap_detective.o created"

# Step 2: Compile main file
echo ""
echo "Step 2: Compiling main file..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c kraken_pipeline_main.cu -o main.o

if [ $? -ne 0 ]; then
    echo "❌ Main file compilation failed"
    exit 1
fi
echo "✓ main.o created"

# Step 3: Compile classifier
echo ""
echo "Step 3: Compiling classifier..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c gpu_kraken_classifier.cu -o classifier.o

if [ $? -ne 0 ]; then
    echo "❌ Classifier compilation failed"
    exit 1
fi
echo "✓ classifier.o created"

# Step 4: Compile database builder
echo ""
echo "Step 4: Compiling database builder..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c gpu_kraken_database_builder.cu -o db_builder.o

if [ $? -ne 0 ]; then
    echo "❌ Database builder compilation failed"
    exit 1
fi
echo "✓ db_builder.o created"

# Step 5: Compile CSV parser
echo ""
echo "Step 5: Compiling CSV parser..."
g++ -std=c++17 -g -O0 -I. -I.. -I../../../include \
    -c sample_csv_parser.cpp -o csv_parser.o

if [ $? -ne 0 ]; then
    echo "❌ CSV parser compilation failed"
    exit 1
fi
echo "✓ csv_parser.o created"

# Step 6: Compile phase1 enhanced classification
echo ""
echo "Step 6: Compiling phase1 enhanced classification..."
if [ -f "phase1_enhanced_classification_cc.cu" ]; then
    nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
         --extended-lambda \
         -c phase1_enhanced_classification_cc.cu -o phase1_enhanced.o
    
    if [ $? -eq 0 ]; then
        echo "✓ phase1_enhanced.o created"
        PHASE1_OBJ="phase1_enhanced.o"
    else
        echo "⚠️  Phase1 enhanced compilation failed, continuing without it"
        PHASE1_OBJ=""
    fi
else
    echo "⚠️  phase1_enhanced_classification_cc.cu not found"
    PHASE1_OBJ=""
fi

# Step 7: Compile minimizer extraction if exists
echo ""
echo "Step 7: Compiling minimizer extraction..."
if [ -f "gpu_minimizer_extraction.cu" ]; then
    nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
         --extended-lambda \
         -c gpu_minimizer_extraction.cu -o minimizer.o
    
    if [ $? -eq 0 ]; then
        echo "✓ minimizer.o created"
        MINIMIZER_OBJ="minimizer.o"
    else
        echo "⚠️  Minimizer compilation failed, continuing without it"
        MINIMIZER_OBJ=""
    fi
else
    echo "⚠️  gpu_minimizer_extraction.cu not found"
    MINIMIZER_OBJ=""
fi

# Step 8: Link everything including heap detective
echo ""
echo "Step 8: Linking all object files..."
echo "Object files to link:"
ls -la *.o

nvcc -std=c++17 -g -O0 \
     main.o classifier.o db_builder.o csv_parser.o heap_detective.o $PHASE1_OBJ $MINIMIZER_OBJ \
     -o debug_kraken_heap_debug \
     -lcuda -lcudart -lz -lstdc++fs

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 SUCCESS! debug_kraken_heap_debug created"
    ls -la debug_kraken_heap_debug
    
    echo ""
    echo "Now you can run with heap debugging enabled!"
    echo "Example:"
    echo "  ./debug_kraken_heap_debug build --genome-dir /path/to/genomes --output /path/to/output"
else
    echo ""
    echo "❌ Linking failed, trying without filesystem library..."
    
    # Try without std::filesystem
    nvcc -std=c++17 -g -O0 \
         main.o classifier.o db_builder.o csv_parser.o heap_detective.o $PHASE1_OBJ $MINIMIZER_OBJ \
         -o debug_kraken_heap_debug \
         -lcuda -lcudart -lz
    
    if [ $? -eq 0 ]; then
        echo "✓ Alternative linking succeeded"
    else
        echo "❌ All linking attempts failed"
        exit 1
    fi
fi

# Make executable
chmod +x debug_kraken_heap_debug

echo ""
echo "=== COMPILATION COMPLETE ==="
echo "Executable: debug_kraken_heap_debug"
echo ""
echo "The heap corruption detective is now integrated!"
echo "It will automatically enable debugging when the program starts."