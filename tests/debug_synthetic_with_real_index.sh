#!/bin/bash

# Debug synthetic reads using the real FQ resistance index
set -e

echo "=== Debugging Synthetic Reads with Real FQ Index ==="

# Configuration
PROJECT_ROOT="$(pwd)"
DATA_DIR="../"
BUILD_DIR="build_integrated"
RESULTS_DIR="results_synthetic_debug"
INDEX_DIR="data/fq_resistance_index"

# Input files (in current directory)
R1_FILE="synthetic_reads_20250529_R1.fastq.gz"
R2_FILE="synthetic_reads_20250529_R2.fastq.gz"

# Output files
mkdir -p ${RESULTS_DIR}
DEBUG_LOG="${RESULTS_DIR}/alignment_debug.log"

echo "Project root: ${PROJECT_ROOT}"
echo "Data directory: ${DATA_DIR}"
echo "Results directory: ${RESULTS_DIR}"
echo "Index directory: ${INDEX_DIR}"
echo ""

# Check the real index exists
echo "=== Checking Real FQ Resistance Index ==="
if [ ! -f "${INDEX_DIR}/kmer_index.bin" ]; then
    echo "❌ ERROR: Real FQ resistance index not found: ${INDEX_DIR}/kmer_index.bin"
    echo "Available files in ${INDEX_DIR}:"
    ls -la ${INDEX_DIR}/ 2>/dev/null || echo "Directory does not exist"
    exit 1
fi

echo "✅ Found real FQ resistance index: ${INDEX_DIR}/kmer_index.bin"
INDEX_SIZE=$(du -h ${INDEX_DIR}/kmer_index.bin | cut -f1)
echo "   Index file size: ${INDEX_SIZE}"

# Show index directory contents
echo "Index directory contents:"
ls -la ${INDEX_DIR}/

# Check input files exist
echo ""
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
    echo "❌ ERROR: FQ pipeline not found: ${BUILD_DIR}/fq_pipeline_gpu"
    echo "Run build script first"
    exit 1
fi

# Quick peek at the FASTQ files
echo ""
echo "=== FASTQ File Preview ==="
echo "First few records from R1:"
zcat ${R1_FILE} | head -n 8
echo ""
echo "Read count estimate for R1:"
READ_COUNT=$(zcat ${R1_FILE} | wc -l | awk '{print int($1/4)}')
echo "${READ_COUNT} reads (estimated)"

echo "First few records from R2:"
zcat ${R2_FILE} | head -n 8

# Get some read length statistics
echo ""
echo "=== Read Length Analysis ==="
echo "Sample read lengths from R1:"
zcat ${R1_FILE} | awk 'NR%4==2 {print length($0)}' | head -10 | tr '\n' ' '
echo ""

AVERAGE_LENGTH=$(zcat ${R1_FILE} | awk 'NR%4==2 {sum+=length($0); count++} END {print int(sum/count)}' | head -1000)
echo "Average read length (first 1000 reads): ${AVERAGE_LENGTH} bp"

# Check GPU status
echo ""
echo "=== GPU Status ==="
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu --format=csv,noheader
    echo ""
else
    echo "nvidia-smi not available"
fi

# Run the FQ pipeline with detailed analysis
OUTPUT_JSON="${RESULTS_DIR}/synthetic_debug_results.json"
echo ""
echo "=== Running FQ Pipeline ==="
echo "Command: ${BUILD_DIR}/fq_pipeline_gpu ${INDEX_DIR} ${R1_FILE} ${R2_FILE} ${OUTPUT_JSON}"
echo "Logging to: ${DEBUG_LOG}"
echo ""
echo "This will analyze the reads in detail..."
echo ""

# Run with detailed logging
${BUILD_DIR}/fq_pipeline_gpu ${INDEX_DIR} ${R1_FILE} ${R2_FILE} ${OUTPUT_JSON} 2>&1 | tee ${DEBUG_LOG}

# Analyze the debug output
echo ""
echo "=== Debug Analysis Summary ==="

# Extract key metrics from the log
if [ -f "${DEBUG_LOG}" ]; then
    echo "Debug log saved to: ${DEBUG_LOG}"
    echo ""
    
    # Look for k-mer statistics
    echo "K-mer Analysis Results:"
    grep -A 5 "K-mer Content Analysis" ${DEBUG_LOG} || echo "  K-mer analysis not found in log"
    
    echo ""
    echo "K-mer Filtering Results:"
    grep -A 10 "K-mer filtering results:" ${DEBUG_LOG} || echo "  K-mer filtering results not found"
    
    echo ""
    echo "Alignment Results:"
    grep -A 10 "Alignment results:" ${DEBUG_LOG} || echo "  Alignment results not found"
    
    # Look for specific performance indicators
    echo ""
    echo "Performance Indicators:"
    KMER_TIME=$(grep "K-mer filtering completed in" ${DEBUG_LOG} | grep -o '[0-9]*' | tail -1)
    ALIGN_TIME=$(grep "Alignment completed in" ${DEBUG_LOG} | grep -o '[0-9]*' | tail -1)
    
    if [ ! -z "$KMER_TIME" ]; then
        echo "  K-mer filtering time: ${KMER_TIME} ms"
    fi
    
    if [ ! -z "$ALIGN_TIME" ]; then
        echo "  Alignment time: ${ALIGN_TIME} ms"
    fi
    
    # Check for errors
    echo ""
    echo "Error Check:"
    if grep -q "ERROR\|error\|Error" ${DEBUG_LOG}; then
        echo "  ⚠️  Errors found in log:"
        grep "ERROR\|error\|Error" ${DEBUG_LOG} | head -5
    else
        echo "  ✅ No errors found"
    fi
    
    # Success indicators
    echo ""
    echo "Success Indicators:"
    TOTAL_CANDIDATES=$(grep "Total candidates:" ${DEBUG_LOG} | grep -o '[0-9]*' | tail -1)
    TOTAL_RESULTS=$(grep "Total alignment results:" ${DEBUG_LOG} | grep -o '[0-9]*' | tail -1)
    MUTATIONS=$(grep "Results with mutations:" ${DEBUG_LOG} | grep -o '[0-9]*' | tail -1)
    
    if [ ! -z "$TOTAL_CANDIDATES" ]; then
        echo "  Total candidates found: ${TOTAL_CANDIDATES}"
    fi
    
    if [ ! -z "$TOTAL_RESULTS" ]; then
        echo "  Alignment results: ${TOTAL_RESULTS}"
    fi
    
    if [ ! -z "$MUTATIONS" ]; then
        echo "  Mutations detected: ${MUTATIONS}"
    fi
else
    echo "❌ Debug log not created"
    exit 1
fi

echo ""
echo "=== Analysis Complete ==="
echo ""
echo "Files created:"
echo "  - Debug log: ${DEBUG_LOG}"
echo ""
echo "Next steps based on results:"
echo "1. If few/no candidates found: Check if reads contain FQ resistance k-mers"
echo "2. If candidates but no alignments: Adjust alignment scoring thresholds"
echo "3. If alignments but no mutations: Check mutation detection logic"
echo "4. If everything works: Scale up to full dataset"
echo ""

# Provide specific recommendations
if [ ! -z "$TOTAL_CANDIDATES" ] && [ "$TOTAL_CANDIDATES" -gt 0 ]; then
    echo "✅ Good: K-mer filtering is finding candidates"
    
    if [ ! -z "$TOTAL_RESULTS" ] && [ "$TOTAL_RESULTS" -gt 0 ]; then
        echo "✅ Good: Alignment stage is working"
        
        if [ ! -z "$MUTATIONS" ] && [ "$MUTATIONS" -gt 0 ]; then
            echo "✅ Excellent: Mutations are being detected!"
            echo "   Ready to run on full dataset"
        else
            echo "⚠️  Issue: No mutations detected - check mutation calling logic"
        fi
    else
        echo "⚠️  Issue: No alignment results - scoring thresholds may be too strict"
    fi
else
    echo "⚠️  Issue: No candidates found - reads may not contain FQ resistance k-mers"
    echo "   Check if synthetic reads include fluoroquinolone resistance genes"
fi