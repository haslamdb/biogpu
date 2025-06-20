#!/bin/bash
# Fixed compilation command for debug build with AddressSanitizer

echo "=== BUILDING DEBUG VERSION WITH ADDRESSSANITIZER ==="

# Clean previous builds
rm -f debug_kraken

# Method 1: Recommended compilation with proper flags
echo "Building with AddressSanitizer..."

nvcc -std=c++17 -g -G -O0 \
     -Xcompiler "-fsanitize=address,undefined -fno-omit-frame-pointer -g3 -fstack-protector-strong" \
     kraken_pipeline_main.cu \
     gpu_kraken_classifier.cu \
     gpu_kraken_database_builder.cu \
     gpu_minimizer_extraction.cu \
     sample_csv_parser.cpp \
     -I. -I.. -I../../../include \
     -o debug_kraken \
     -lcuda -lcudart -lz -lstdc++fs \
     -Xlinker "-fsanitize=address,undefined"

# Check if build succeeded
if [ $? -eq 0 ]; then
    echo "✓ Debug build successful!"
    
    echo ""
    echo "=== RUNNING WITH ADDRESSSANITIZER ==="
    echo "Environment variables for optimal debugging:"
    echo ""
    
    # Set up environment for AddressSanitizer
    export ASAN_OPTIONS="abort_on_error=1:fast_unwind_on_malloc=0:detect_leaks=1:check_initialization_order=1:strict_init_order=1:detect_odr_violation=2"
    export UBSAN_OPTIONS="abort_on_error=1:print_stacktrace=1"
    export MALLOC_CHECK_=2
    
    echo "ASAN_OPTIONS=\"$ASAN_OPTIONS\""
    echo "UBSAN_OPTIONS=\"$UBSAN_OPTIONS\""
    echo "MALLOC_CHECK_=$MALLOC_CHECK_"
    echo ""
    
    echo "Now run your program:"
    echo "./debug_kraken build --genome-dir YOUR_GENOME_DIR --output test_db"
    echo ""
    echo "If it still crashes, try the minimal version first:"
    echo ""
    
    # Create minimal test version
    echo "=== CREATING MINIMAL TEST VERSION ==="
    
    cat > minimal_test.cu << 'EOF'
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::cout << "Minimal test starting..." << std::endl;
    
    if (argc < 2) {
        std::cout << "Need at least one argument" << std::endl;
        return 1;
    }
    
    std::cout << "argc = " << argc << std::endl;
    
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
        
        // Test string creation
        std::string arg(argv[i]);
        std::cout << "  String length: " << arg.length() << std::endl;
    }
    
    std::cout << "Minimal test completed successfully!" << std::endl;
    return 0;
}
EOF
    
    echo "Building minimal test..."
    nvcc -std=c++17 -g -O0 \
         -Xcompiler "-fsanitize=address -fno-omit-frame-pointer" \
         minimal_test.cu -o minimal_test \
         -Xlinker "-fsanitize=address"
    
    if [ $? -eq 0 ]; then
        echo "✓ Minimal test built successfully!"
        echo "Run it first to verify AddressSanitizer is working:"
        echo "./minimal_test build --genome-dir test --output test"
    fi
    
else
    echo "❌ Build failed!"
    echo ""
    echo "=== ALTERNATIVE BUILD WITHOUT ADDRESSSANITIZER ==="
    echo "If AddressSanitizer doesn't work with CUDA, try this:"
    echo ""
    
    cat > build_without_asan.sh << 'EOF'
#!/bin/bash
# Build without AddressSanitizer but with debug symbols

nvcc -std=c++17 -g -G -O0 \
     -Xcompiler "-g3 -fstack-protector-strong -fno-omit-frame-pointer" \
     kraken_pipeline_main.cu \
     gpu_kraken_classifier.cu \
     gpu_kraken_database_builder.cu \
     gpu_minimizer_extraction.cu \
     sample_csv_parser.cpp \
     -I. -I.. -I../../../include \
     -o debug_kraken_no_asan \
     -lcuda -lcudart -lz -lstdc++fs

echo "Built debug_kraken_no_asan"
echo "Run with: MALLOC_CHECK_=2 ./debug_kraken_no_asan build --genome-dir YOUR_DIR --output test_db"
echo "Or with gdb: gdb ./debug_kraken_no_asan"
EOF
    
    chmod +x build_without_asan.sh
    echo "Created build_without_asan.sh - run it if AddressSanitizer doesn't work"
fi

echo ""
echo "=== DEBUGGING TIPS ==="
echo "1. Run minimal_test first to verify AddressSanitizer works"
echo "2. If minimal_test works but main program crashes, the issue is in your code"
echo "3. If minimal_test also crashes, there's a system-level issue"
echo "4. Use 'ulimit -c unlimited' to enable core dumps"
echo "5. If you get 'AddressSanitizer failed to allocate', reduce memory usage"