#!/bin/bash
# Build without AddressSanitizer but with debug symbols

nvcc -std=c++17 -g -G -O0 --extended-lambda \
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
