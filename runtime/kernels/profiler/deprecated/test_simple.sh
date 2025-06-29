#!/bin/bash
# Test the simple executable for heap corruption

echo "=== TESTING SIMPLE EXECUTABLE FOR HEAP CORRUPTION ==="

# Set up heap debugging environment
export MALLOC_CHECK_=2
export ASAN_OPTIONS="abort_on_error=1:fast_unwind_on_malloc=0:detect_leaks=1"

echo "Environment set up for heap debugging"
echo "MALLOC_CHECK_=$MALLOC_CHECK_"

# Test 1: Basic execution
echo ""
echo "=== TEST 1: Basic execution ==="
echo "Running: ./debug_kraken_fixed"

./debug_kraken_fixed
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "✓ Basic execution successful - no heap corruption!"
elif [ $exit_code -eq 134 ]; then
    echo "❌ SIGABRT - heap corruption still present in basic execution"
elif [ $exit_code -eq 139 ]; then
    echo "❌ SIGSEGV - segmentation fault in basic execution"
else
    echo "⚠️  Exit code: $exit_code"
fi

# Test 2: With arguments (similar to the original crash)
echo ""
echo "=== TEST 2: With build arguments ==="
echo "Running: ./debug_kraken_fixed build --genome-dir /tmp/test --output /tmp/out"

./debug_kraken_fixed build --genome-dir /tmp/test --output /tmp/out
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "✓ Argument parsing successful - no heap corruption!"
elif [ $exit_code -eq 134 ]; then
    echo "❌ SIGABRT - heap corruption still present in argument parsing"
elif [ $exit_code -eq 139 ]; then
    echo "❌ SIGSEGV - segmentation fault in argument parsing"
else
    echo "⚠️  Exit code: $exit_code"
fi

# Test 3: With the exact arguments that were causing the crash
echo ""
echo "=== TEST 3: With exact problematic arguments ==="
echo "Running: ./debug_kraken_fixed build --genome-dir ../../../../biogpu/data/type_strain_genomes_processed --output type_strain_db"

./debug_kraken_fixed build --genome-dir ../../../../biogpu/data/type_strain_genomes_processed --output type_strain_db
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "✓ Exact arguments successful - no heap corruption!"
elif [ $exit_code -eq 134 ]; then
    echo "❌ SIGABRT - heap corruption still present with exact arguments"
elif [ $exit_code -eq 139 ]; then
    echo "❌ SIGSEGV - segmentation fault with exact arguments"
else
    echo "⚠️  Exit code: $exit_code"
fi

# Summary
echo ""
echo "=== SUMMARY ==="
if [ $exit_code -eq 0 ]; then
    echo "🎉 GREAT NEWS: The simple test shows no heap corruption!"
    echo "This suggests the heap corruption was related to the complex linking/multiple definitions."
    echo ""
    echo "Next step: Fix the linking issues to build the full program."
    echo "The multiple definition fix appears to have resolved the core heap corruption."
else
    echo "❌ Heap corruption still present even in simple test"
    echo "The issue is deeper than just multiple definitions."
    echo ""
    echo "Next steps:"
    echo "1. The heap corruption is in basic C++ standard library usage"
    echo "2. Could be compiler/environment issue"
    echo "3. Need to investigate with simpler test cases"
fi

echo ""
echo "=== NEXT ACTIONS ==="
if [ $exit_code -eq 0 ]; then
    echo "Since basic heap operations work, let's fix the linking issue:"
    echo ""
    echo "1. Check for missing dependencies:"
    echo "   ldd ./debug_kraken_fixed"
    echo ""
    echo "2. Try building with different library flags:"
    echo "   nvcc -std=c++17 -g main.o classifier.o db_builder.o csv_parser.o \\"
    echo "        -o full_kraken -lcuda -lcudart -lz"
    echo ""
    echo "3. Check what's causing the linking to fail:"
    echo "   nm main.o | grep -E '(undefined|U )'"
    echo "   nm classifier.o | grep -E '(undefined|U )'"
else
    echo "Since heap corruption persists even in simple cases:"
    echo ""
    echo "1. Try different compiler:"
    echo "   g++ -std=c++17 -g simple_test.cu -o test_gcc"
    echo ""
    echo "2. Check system integrity:"
    echo "   ulimit -c unlimited"
    echo "   ./debug_kraken_fixed"
    echo "   gdb ./debug_kraken_fixed core"
    echo ""
    echo "3. Try with different memory allocator:"
    echo "   LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libc_malloc_debug.so ./debug_kraken_fixed"
fi
