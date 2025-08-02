#!/bin/bash

# Test script for unified BioGPU pipeline

echo "Testing BioGPU Unified Pipeline"
echo "================================"

# Set paths
DATA_DIR="/home/david/projects/biogpu/data"
OUTPUT_DIR="/home/david/projects/biogpu/test_output_$(date +%Y%m%d_%H%M%S)"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Test files - using smaller synthetic reads for testing
R1_FILE="$DATA_DIR/synthetic_reads_R1.fastq.gz"
R2_FILE="$DATA_DIR/synthetic_reads_R2.fastq.gz"

# Database paths
REFERENCE_DB="$DATA_DIR/pathogen_profiler_db/pathogen_subset.db"  # Using smaller subset for testing
RESISTANCE_DB="$DATA_DIR/quinolone_resistance_mutation_table.csv"
SAMPLE_ID="test_sample_001"

echo "Input files:"
echo "  R1: $R1_FILE"
echo "  R2: $R2_FILE"
echo "  Reference DB: $REFERENCE_DB"
echo "  Resistance DB: $RESISTANCE_DB"
echo "  Output: $OUTPUT_DIR"
echo ""

# Check if files exist
if [ ! -f "$R1_FILE" ]; then
    echo "Error: R1 file not found: $R1_FILE"
    exit 1
fi

if [ ! -f "$R2_FILE" ]; then
    echo "Error: R2 file not found: $R2_FILE"
    exit 1
fi

# Run the unified pipeline
echo "Starting unified pipeline..."
echo "Command: ./bio_gpu_pipeline --r1 \"$R1_FILE\" --r2 \"$R2_FILE\" --output-dir \"$OUTPUT_DIR\" --reference-db \"$REFERENCE_DB\" --resistance-db \"$RESISTANCE_DB\" --sample-id \"$SAMPLE_ID\" --use-multi-gpu --progress-json"
echo ""

./bio_gpu_pipeline \
    --r1 "$R1_FILE" \
    --r2 "$R2_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --reference-db "$REFERENCE_DB" \
    --resistance-db "$RESISTANCE_DB" \
    --sample-id "$SAMPLE_ID" \
    --use-multi-gpu \
    --progress-json

RESULT=$?

echo ""
echo "Pipeline finished with exit code: $RESULT"

if [ $RESULT -eq 0 ]; then
    echo ""
    echo "Checking output files..."
    echo "========================"
    
    # Check resistance output
    if [ -d "$OUTPUT_DIR/resistance" ]; then
        echo "Resistance output:"
        ls -la "$OUTPUT_DIR/resistance/"
    else
        echo "Warning: No resistance output directory found"
    fi
    
    echo ""
    
    # Check AMR output
    if [ -d "$OUTPUT_DIR/amr" ]; then
        echo "AMR output:"
        ls -la "$OUTPUT_DIR/amr/"
    else
        echo "Warning: No AMR output directory found"
    fi
    
    echo ""
    echo "Test completed successfully!"
else
    echo "Pipeline failed!"
fi