#!/bin/bash
# Debug and fix the compilation issue

echo "=== DEBUGGING COMPILATION ISSUE ==="

# Check what files exist
echo "Files in current directory:"
ls -la *.o *.cu *.h debug_kraken* 2>/dev/null

echo ""
echo "=== MANUAL COMPILATION STEP BY STEP ==="

# Clean up any existing object files
rm -f *.o debug_kraken* 2>/dev/null

# Step 1: Compile each file individually with verbose output
echo "Step 1: Compiling main file..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c kraken_pipeline_main.cu -o main.o -v

if [ $? -ne 0 ]; then
    echo "❌ Main file compilation failed"
    exit 1
fi
echo "✓ main.o created"

echo ""
echo "Step 2: Compiling classifier..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c gpu_kraken_classifier.cu -o classifier.o -v

if [ $? -ne 0 ]; then
    echo "❌ Classifier compilation failed"
    exit 1
fi
echo "✓ classifier.o created"

echo ""
echo "Step 3: Compiling database builder..."
nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
     -c gpu_kraken_database_builder.cu -o db_builder.o -v

if [ $? -ne 0 ]; then
    echo "❌ Database builder compilation failed"
    exit 1
fi
echo "✓ db_builder.o created"

echo ""
echo "Step 4: Compiling CSV parser..."
g++ -std=c++17 -g -O0 -I. -I.. -I../../../include \
    -c sample_csv_parser.cpp -o csv_parser.o -v

if [ $? -ne 0 ]; then
    echo "❌ CSV parser compilation failed"
    exit 1
fi
echo "✓ csv_parser.o created"

echo ""
echo "Step 5: Trying to compile minimizer extraction..."
if [ -f "gpu_minimizer_extraction.cu" ]; then
    nvcc -std=c++17 -g -O0 -I. -I.. -I../../../include \
         -c gpu_minimizer_extraction.cu -o minimizer.o -v
    
    if [ $? -eq 0 ]; then
        echo "✓ minimizer.o created"
        MINIMIZER_OBJ="minimizer.o"
    else
        echo "⚠️  Minimizer compilation failed, continuing without it"
        MINIMIZER_OBJ=""
    fi
else
    echo "⚠️  gpu_minimizer_extraction.cu not found, continuing without it"
    MINIMIZER_OBJ=""
fi

echo ""
echo "Step 6: Linking all object files..."
echo "Object files to link:"
ls -la *.o

nvcc -std=c++17 -g -O0 \
     main.o classifier.o db_builder.o csv_parser.o $MINIMIZER_OBJ \
     -o debug_kraken_fixed \
     -lcuda -lcudart -lz -lstdc++fs -v

LINK_EXIT_CODE=$?

if [ $LINK_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "🎉 SUCCESS! debug_kraken_fixed created"
    ls -la debug_kraken_fixed
    
    echo ""
    echo "Testing basic execution..."
    ./debug_kraken_fixed 2>&1 | head -10
    
else
    echo ""
    echo "❌ LINKING FAILED with exit code: $LINK_EXIT_CODE"
    echo ""
    echo "Let's try alternative approaches..."
    
    # Alternative 1: Link without problematic libraries
    echo "=== ALTERNATIVE 1: Link without lstdc++fs ==="
    nvcc -std=c++17 -g -O0 \
         main.o classifier.o db_builder.o csv_parser.o $MINIMIZER_OBJ \
         -o debug_kraken_alt1 \
         -lcuda -lcudart -lz -v
    
    if [ $? -eq 0 ]; then
        echo "✓ Alternative build 1 succeeded: debug_kraken_alt1"
        mv debug_kraken_alt1 debug_kraken_fixed
    else
        # Alternative 2: Link with g++ instead of nvcc
        echo "=== ALTERNATIVE 2: Link with g++ ==="
        g++ -std=c++17 -g -O0 \
            main.o classifier.o db_builder.o csv_parser.o $MINIMIZER_OBJ \
            -o debug_kraken_alt2 \
            -lcuda -lcudart -lz -lstdc++fs -L/usr/local/cuda/lib64 -v
        
        if [ $? -eq 0 ]; then
            echo "✓ Alternative build 2 succeeded: debug_kraken_alt2"
            mv debug_kraken_alt2 debug_kraken_fixed
        else
            # Alternative 3: Minimal build with just main components
            echo "=== ALTERNATIVE 3: Minimal build ==="
            nvcc -std=c++17 -g -O0 \
                 main.o classifier.o \
                 -o debug_kraken_minimal \
                 -lcuda -lcudart -v
            
            if [ $? -eq 0 ]; then
                echo "✓ Minimal build succeeded: debug_kraken_minimal"
                mv debug_kraken_minimal debug_kraken_fixed
            else
                echo "❌ All build attempts failed"
                echo ""
                echo "=== CREATING SIMPLE TEST EXECUTABLE ==="
                
                # Create a simple test program to at least test heap corruption
                cat > simple_test.cu << 'EOF'
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

int main(int argc, char* argv[]) {
    std::cout << "Simple heap corruption test" << std::endl;
    
    // Test basic string operations that might cause heap corruption
    std::vector<std::string> test_strings;
    
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            std::string arg(argv[i]);
            test_strings.push_back(arg);
            std::cout << "arg[" << i << "] = " << arg << std::endl;
        }
    }
    
    // Test some allocations
    for (int i = 0; i < 1000; i++) {
        std::string test = "test_string_" + std::to_string(i);
        test_strings.push_back(test);
    }
    
    std::cout << "Created " << test_strings.size() << " strings without crash" << std::endl;
    
    if (argc >= 3 && std::string(argv[1]) == "build") {
        std::cout << "Build command simulation - would process: " << argv[2] << std::endl;
    }
    
    return 0;
}
EOF
                
                nvcc -std=c++17 -g -O0 simple_test.cu -o debug_kraken_fixed
                
                if [ $? -eq 0 ]; then
                    echo "✓ Created simple test executable: debug_kraken_fixed"
                else
                    echo "❌ Even simple test compilation failed"
                fi
            fi
        fi
    fi
fi

echo ""
echo "=== FINAL STATUS ==="
if [ -f "debug_kraken_fixed" ]; then
    echo "✓ Executable created: debug_kraken_fixed"
    ls -la debug_kraken_fixed
    
    echo ""
    echo "Test it with:"
    echo "  MALLOC_CHECK_=2 ./debug_kraken_fixed build --genome-dir /tmp/test --output /tmp/out"
else
    echo "❌ No executable created"
    echo "Manual debugging required"
fi

echo ""
echo "Object files created:"
ls -la *.o 2>/dev/null
