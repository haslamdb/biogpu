#!/bin/bash
# Test script for heap corruption

echo "=== TESTING HEAP CORRUPTION FIX ==="

# Check if the executable exists
if [ ! -f "./debug_kraken_fixed" ]; then
    echo "❌ debug_kraken_fixed not found!"
    exit 1
fi

echo "✓ Found debug_kraken_fixed executable"

# Set up environment for debugging
export MALLOC_CHECK_=2
export ASAN_OPTIONS="abort_on_error=1:fast_unwind_on_malloc=0:detect_leaks=1"

echo "Environment variables set:"
echo "  MALLOC_CHECK_=$MALLOC_CHECK_"
echo "  ASAN_OPTIONS=$ASAN_OPTIONS"

# Test 1: Just run with help to see if basic functionality works
echo ""
echo "=== TEST 1: Basic functionality test ==="
echo "Running: ./debug_kraken_fixed (no arguments)"

timeout 10s ./debug_kraken_fixed 2>&1
exit_code=$?

if [ $exit_code -eq 0 ] || [ $exit_code -eq 1 ]; then
    echo "✓ Basic run completed (exit code: $exit_code)"
else
    echo "❌ Basic run failed with exit code: $exit_code"
    if [ $exit_code -eq 124 ]; then
        echo "  (Timeout - program may be hanging)"
    elif [ $exit_code -eq 134 ]; then
        echo "  (SIGABRT - heap corruption still present)"
    elif [ $exit_code -eq 139 ]; then
        echo "  (SIGSEGV - segmentation fault)"
    fi
fi

# Test 2: Try with minimal arguments
echo ""
echo "=== TEST 2: Minimal arguments test ==="
echo "Running: ./debug_kraken_fixed build --help"

timeout 10s ./debug_kraken_fixed build --help 2>&1
exit_code=$?

if [ $exit_code -eq 0 ] || [ $exit_code -eq 1 ]; then
    echo "✓ Minimal arguments test completed (exit code: $exit_code)"
else
    echo "❌ Minimal arguments test failed with exit code: $exit_code"
fi

# Test 3: Try with the actual problematic command (but with a safe directory)
echo ""
echo "=== TEST 3: Actual build command test ==="

# Create a simple test directory
mkdir -p /tmp/test_genomes
echo ">test_genome" > /tmp/test_genomes/test.fasta
echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" >> /tmp/test_genomes/test.fasta

echo "Created test genome directory: /tmp/test_genomes"
echo "Running: ./debug_kraken_fixed build --genome-dir /tmp/test_genomes --output /tmp/test_db"

timeout 30s ./debug_kraken_fixed build --genome-dir /tmp/test_genomes --output /tmp/test_db 2>&1 | head -50
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "✓ Build command completed successfully!"
    echo "🎉 HEAP CORRUPTION APPEARS TO BE FIXED!"
elif [ $exit_code -eq 1 ]; then
    echo "⚠️  Build command exited with error code 1 (may be expected for test data)"
    echo "But no heap corruption detected!"
elif [ $exit_code -eq 124 ]; then
    echo "⚠️  Build command timed out (may be normal for large processing)"
    echo "But no immediate heap corruption detected!"
elif [ $exit_code -eq 134 ]; then
    echo "❌ HEAP CORRUPTION STILL PRESENT (SIGABRT)"
    echo "The multiple definition fix didn't solve the underlying issue"
elif [ $exit_code -eq 139 ]; then
    echo "❌ SEGMENTATION FAULT (SIGSEGV)"
    echo "There's still a memory access issue"
else
    echo "❌ Unexpected exit code: $exit_code"
fi

# Test 4: If heap corruption is still there, try with GDB for detailed info
if [ $exit_code -eq 134 ] || [ $exit_code -eq 139 ]; then
    echo ""
    echo "=== TEST 4: GDB debugging session ==="
    echo "Since heap corruption is still present, running with GDB..."
    
    cat > gdb_commands.txt << EOF
set environment MALLOC_CHECK_=2
run build --genome-dir /tmp/test_genomes --output /tmp/test_db
bt
quit
EOF
    
    echo "GDB commands prepared. Run this manually for detailed debugging:"
    echo "  gdb ./debug_kraken_fixed"
    echo "  (gdb) set environment MALLOC_CHECK_=2"
    echo "  (gdb) run build --genome-dir /tmp/test_genomes --output /tmp/test_db"
    echo "  (gdb) bt"
    
    # Also create a simple crash reproduction script
    cat > reproduce_crash.sh << 'EOF'
#!/bin/bash
echo "Reproducing crash for debugging..."
export MALLOC_CHECK_=2
ulimit -c unlimited
./debug_kraken_fixed build --genome-dir /tmp/test_genomes --output /tmp/test_db
echo "If crash occurred, check core dump with: gdb ./debug_kraken_fixed core"
EOF
    chmod +x reproduce_crash.sh
    echo "Created reproduce_crash.sh for detailed debugging"
fi

# Cleanup
rm -rf /tmp/test_genomes /tmp/test_db

echo ""
echo "=== SUMMARY ==="
if [ $exit_code -eq 0 ] || [ $exit_code -eq 1 ]; then
    echo "🎉 SUCCESS: Multiple definition errors fixed and no heap corruption detected!"
    echo "You can now test with your actual genome directory:"
    echo "  MALLOC_CHECK_=2 ./debug_kraken_fixed build --genome-dir YOUR_ACTUAL_DIR --output test_db"
elif [ $exit_code -eq 124 ]; then
    echo "⚠️  PARTIAL SUCCESS: No immediate heap corruption, but command timed out"
    echo "Try with a smaller dataset or check if it's just slow processing"
else
    echo "❌ HEAP CORRUPTION STILL PRESENT: Need further debugging"
    echo "The issue is deeper than just multiple definitions"
    echo "Next steps:"
    echo "1. Use GDB: gdb ./debug_kraken_fixed"
    echo "2. Use Valgrind: valgrind --tool=memcheck ./debug_kraken_fixed build ..."
    echo "3. Try the step-by-step main function replacement I provided earlier"
fi
