#!/bin/bash

# Test script for running FQ resistance pipeline with synthetic reads
set -e

echo "=== Testing FQ Pipeline with Synthetic Reads ==="

# Configuration
PROJECT_ROOT="$(pwd)"
DATA_DIR="../"
BUILD_DIR="../build_integrated"
RESULTS_DIR="results_synthetic"

# Input files (in parent directory)
R1_FILE="${DATA_DIR}/synthetic_reads_20250529_R1.fastq.gz"
R2_FILE="${DATA_DIR}/synthetic_reads_20250529_R2.fastq.gz"

# Output files
mkdir -p ${RESULTS_DIR}
OUTPUT_JSON="${RESULTS_DIR}/synthetic_resistance_results.json"
LOG_FILE="${RESULTS_DIR}/pipeline_debug.log"

echo "Project root: ${PROJECT_ROOT}"
echo "Data directory: ${DATA_DIR}"
echo "Results directory: ${RESULTS_DIR}"
echo ""

# Check input files exist
echo "=== Checking Input Files ==="
if [ ! -f "${R1_FILE}" ]; then
    echo "❌ ERROR: R1 file not found: ${R1_FILE}"
    echo "Current directory: $(pwd)"
    echo "Contents of parent directory:"
    ls -la ../
    exit 1
fi

if [ ! -f "${R2_FILE}" ]; then
    echo "❌ ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

echo "✅ R1 file found: ${R1_FILE}"
echo "   Size: $(du -h ${R1_FILE} | cut -f1)"
echo "✅ R2 file found: ${R2_FILE}"
echo "   Size: $(du -h ${R2_FILE} | cut -f1)"

# Check if pipeline executable exists
if [ ! -f "${BUILD_DIR}/fq_pipeline_gpu" ]; then
    echo "❌ ERROR: Pipeline executable not found: ${BUILD_DIR}/fq_pipeline_gpu"
    echo "Run ./scripts/build_integrated.sh first"
    exit 1
fi

# Quick peek at the FASTQ files
echo ""
echo "=== FASTQ File Preview ==="
echo "First few records from R1:"
zcat ${R1_FILE} | head -n 8
echo ""
echo "Read count estimate:"
zcat ${R1_FILE} | wc -l | awk '{print int($1/4) " reads (estimated)"}'

# Use the real FQ resistance index
INDEX_DIR="../data/fq_resistance_index"
echo ""
echo "=== Using Real FQ Resistance Index ==="

if [ ! -f "${INDEX_DIR}/kmer_index.bin" ]; then
    echo "❌ ERROR: FQ resistance index not found: ${INDEX_DIR}/kmer_index.bin"
    echo "Available files in ${INDEX_DIR}:"
    ls -la ${INDEX_DIR}/ 2>/dev/null || echo "Directory does not exist"
    exit 1
fi

echo "✅ Found FQ resistance index: ${INDEX_DIR}/kmer_index.bin"
echo "   Index file: $(ls -lh ${INDEX_DIR}/kmer_index.bin)"

# Check for other index files
echo "Index directory contents:"
ls -la ${INDEX_DIR}/

# Run the pipeline
echo ""
echo "=== Running FQ Resistance Pipeline ==="
echo "Command: ${BUILD_DIR}/fq_pipeline_gpu ${INDEX_DIR} ${R1_FILE} ${R2_FILE} ${OUTPUT_JSON}"
echo "Logging to: ${LOG_FILE}"
echo ""

# Run with detailed logging
${BUILD_DIR}/fq_pipeline_gpu ${INDEX_DIR} ${R1_FILE} ${R2_FILE} ${OUTPUT_JSON} 2>&1 | tee ${LOG_FILE}

# Check results
echo ""
echo "=== Pipeline Results ==="
if [ -f "${OUTPUT_JSON}" ]; then
    echo "✅ Output file created: ${OUTPUT_JSON}"
    echo "   Size: $(du -h ${OUTPUT_JSON} | cut -f1)"
    echo ""
    echo "Output content:"
    cat ${OUTPUT_JSON}
else
    echo "❌ ERROR: Output file not created"
    echo "Check the log file: ${LOG_FILE}"
    exit 1
fi

# Analyze the log for debugging info
echo ""
echo "=== Debugging Analysis ==="
echo "Checking for key pipeline stages in log..."

if grep -q "Enhanced k-mer filter started" ${LOG_FILE}; then
    echo "✅ Stage 1 (K-mer filtering) started"
else
    echo "❌ Stage 1 (K-mer filtering) not found in log"
fi

if grep -q "Enhanced k-mer filter completed" ${LOG_FILE}; then
    echo "✅ Stage 1 (K-mer filtering) completed"
else
    echo "❌ Stage 1 (K-mer filtering) did not complete"
fi

if grep -q "Simple alignment kernel started" ${LOG_FILE}; then
    echo "✅ Stage 2 (Alignment) started"
else
    echo "❌ Stage 2 (Alignment) not found in log"
fi

if grep -q "Simple alignment kernel completed" ${LOG_FILE}; then
    echo "✅ Stage 2 (Alignment) completed"
else
    echo "❌ Stage 2 (Alignment) did not complete"
fi

# Extract key metrics
echo ""
echo "=== Performance Metrics ==="
TOTAL_READS=$(grep "total_reads" ${OUTPUT_JSON} 2>/dev/null | grep -o '[0-9]*' | tail -1)
MUTATIONS_FOUND=$(grep "mutations_found" ${OUTPUT_JSON} 2>/dev/null | grep -o '[0-9]*' | tail -1)

if [ ! -z "$TOTAL_READS" ]; then
    echo "Total reads processed: ${TOTAL_READS}"
fi

if [ ! -z "$MUTATIONS_FOUND" ]; then
    echo "Mutations detected: ${MUTATIONS_FOUND}"
fi

# Look for candidate counts in the log
echo ""
echo "K-mer filtering results:"
grep "candidates" ${LOG_FILE} | head -10

echo ""
echo "Alignment results:"
grep "alignment results" ${LOG_FILE}

# Memory usage
echo ""
echo "=== GPU Memory Info ==="
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits
else
    echo "nvidia-smi not available"
fi

echo ""
echo "=== Test Complete ==="
echo "Results saved to: ${RESULTS_DIR}/"
echo "  - Pipeline output: ${OUTPUT_JSON}"
echo "  - Debug log: ${LOG_FILE}"
echo "  - Index creation log: ${RESULTS_DIR}/index_creation.log"

# Summary
if [ -f "${OUTPUT_JSON}" ] && grep -q "pipeline_ready\|mutations" ${OUTPUT_JSON}; then
    echo "✅ Pipeline test completed successfully!"
    echo ""
    echo "Next steps for debugging alignment:"
    echo "1. Check candidate counts in ${LOG_FILE}"
    echo "2. Look for CUDA kernel errors"
    echo "3. Verify k-mer matches are being found"
    echo "4. Check alignment scoring parameters"
else
    echo "❌ Pipeline test had issues - check logs"
    exit 1
fi