#!/bin/bash
# Verify and fix the include in gpu_kraken_database_builder.cu

echo "=== CHECKING DATABASE BUILDER INCLUDE ==="

# Check what's currently included
echo "Current include in gpu_kraken_database_builder.cu:"
head -20 gpu_kraken_database_builder.cu | grep -E "#include.*gpu_kraken_classifier"

# Check if it's including the .cu file (bad) or .h file (good)
if grep -q '#include "gpu_kraken_classifier\.cu"' gpu_kraken_database_builder.cu; then
    echo "❌ FOUND PROBLEMATIC INCLUDE: gpu_kraken_classifier.cu"
    echo "Fixing it now..."
    
    # Create backup
    cp gpu_kraken_database_builder.cu gpu_kraken_database_builder.cu.backup
    
    # Replace the problematic include
    sed -i 's/#include "gpu_kraken_classifier\.cu"/#include "gpu_kraken_classifier.h"/' gpu_kraken_database_builder.cu
    
    echo "✓ Fixed: Changed to #include \"gpu_kraken_classifier.h\""
    
elif grep -q '#include "gpu_kraken_classifier\.h"' gpu_kraken_database_builder.cu; then
    echo "✓ GOOD: Already includes gpu_kraken_classifier.h"
    
else
    echo "❓ No gpu_kraken_classifier include found. Adding it..."
    
    # Add the header include after other includes
    sed -i '/^#include "gpu_minimizer_extraction\.cuh"/a #include "gpu_kraken_classifier.h"' gpu_kraken_database_builder.cu
    
    echo "✓ Added: #include \"gpu_kraken_classifier.h\""
fi

echo ""
echo "=== TESTING COMPILATION ==="

# Test compilation with separate object files
echo "Step 1: Compiling main file..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c kraken_pipeline_main.cu -o main.o 2>&1 | head -20

if [ $? -eq 0 ]; then
    echo "✓ Main file compiled successfully"
    
    echo "Step 2: Compiling classifier..."
    nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
         -c gpu_kraken_classifier.cu -o classifier.o 2>&1 | head -20
    
    if [ $? -eq 0 ]; then
        echo "✓ Classifier compiled successfully"
        
        echo "Step 3: Compiling database builder..."
        nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
             -c gpu_kraken_database_builder.cu -o db_builder.o 2>&1 | head -20
        
        if [ $? -eq 0 ]; then
            echo "✓ Database builder compiled successfully"
            
            echo "Step 4: Compiling other files..."
            
            # CSV parser
            g++ -std=c++17 -g -O0 -I. -I.. -I../../../include \
                -c sample_csv_parser.cpp -o csv_parser.o 2>/dev/null
            
            # Minimizer extraction
            nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
                 -c gpu_minimizer_extraction.cu -o minimizer.o 2>/dev/null
            
            echo "Step 5: Linking..."
            nvcc -std=c++17 -g -O0 \
                 main.o classifier.o db_builder.o csv_parser.o minimizer.o \
                 -o debug_kraken_fixed \
                 -lcuda -lcudart -lz -lstdc++fs 2>&1 | head -20
            
            if [ $? -eq 0 ]; then
                echo ""
                echo "🎉 SUCCESS! Built debug_kraken_fixed"
                echo ""
                echo "Now test the heap corruption issue:"
                echo "MALLOC_CHECK_=2 ./debug_kraken_fixed build --genome-dir YOUR_DIR --output test_db"
                echo ""
                echo "Or with step-by-step debugging:"
                echo "gdb ./debug_kraken_fixed"
                
            else
                echo "❌ Linking failed - check errors above"
            fi
        else
            echo "❌ Database builder compilation failed - check errors above"
        fi
    else
        echo "❌ Classifier compilation failed - check errors above"
    fi
else
    echo "❌ Main file compilation failed - check errors above"
fi

echo ""
echo "=== CLEANUP ==="
rm -f *.o

echo "Done!"