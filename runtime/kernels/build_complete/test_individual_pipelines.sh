#!/bin/bash

# Test individual pipelines separately

echo "Testing Individual BioGPU Pipelines"
echo "==================================="

# Set paths
DATA_DIR="/home/david/projects/biogpu/data"
OUTPUT_DIR="/home/david/projects/biogpu/test_individual_$(date +%Y%m%d_%H%M%S)"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Test with smaller files first
# Let's create a small test FASTQ file
echo "Creating small test FASTQ files..."
zcat "$DATA_DIR/synthetic_reads_R1.fastq.gz" | head -4000 > "$OUTPUT_DIR/test_R1.fastq"
zcat "$DATA_DIR/synthetic_reads_R2.fastq.gz" | head -4000 > "$OUTPUT_DIR/test_R2.fastq"

echo "Test files created:"
ls -lh "$OUTPUT_DIR"/*.fastq

echo ""
echo "1. Testing Resistance Pipeline"
echo "==============================="

# Check if resistance pipeline exists
RESISTANCE_EXEC="/home/david/projects/biogpu/runtime/kernels/resistance/build/clean_resistance_pipeline"
if [ -f "$RESISTANCE_EXEC" ]; then
    echo "Found resistance pipeline at: $RESISTANCE_EXEC"
    echo "Running resistance pipeline help..."
    "$RESISTANCE_EXEC" --help 2>&1 | head -20
else
    echo "Error: Resistance pipeline not found at $RESISTANCE_EXEC"
fi

echo ""
echo "2. Testing AMR Pipeline"
echo "======================="

# Check if AMR pipeline exists
AMR_EXEC="/home/david/projects/biogpu/runtime/kernels/amr_gpu_pipeline/gpu_amr_detection"
if [ -f "$AMR_EXEC" ]; then
    echo "Found AMR pipeline at: $AMR_EXEC"
    echo "Checking AMR pipeline..."
    # Create a simple CSV file for AMR pipeline
    echo "sample_name,fastq_path" > "$OUTPUT_DIR/test_samples.csv"
    echo "test_sample,$OUTPUT_DIR/test_R1.fastq" >> "$OUTPUT_DIR/test_samples.csv"
    
    # Test AMR pipeline with minimal data
    echo "Testing AMR detection..."
    "$AMR_EXEC" "$DATA_DIR/AMR_CDS.fa" "$OUTPUT_DIR/test_samples.csv" "$OUTPUT_DIR/amr_test" --device 0 2>&1 | head -20
else
    echo "Error: AMR pipeline not found at $AMR_EXEC"
fi

echo ""
echo "3. Testing Unified Pipeline"
echo "==========================="

UNIFIED_EXEC="./bio_gpu_pipeline"
if [ -f "$UNIFIED_EXEC" ]; then
    echo "Found unified pipeline at: $UNIFIED_EXEC"
    
    # Run with minimal test data
    echo "Running unified pipeline with test data..."
    "$UNIFIED_EXEC" \
        --r1 "$OUTPUT_DIR/test_R1.fastq" \
        --r2 "$OUTPUT_DIR/test_R2.fastq" \
        --output-dir "$OUTPUT_DIR/unified_test" \
        --reference-db "$DATA_DIR/pathogen_profiler_db/pathogen_subset.db" \
        --resistance-db "$DATA_DIR/quinolone_resistance_mutation_table.csv" \
        --sample-id "test_minimal" \
        --batch-size 1000
    
    echo ""
    echo "Unified pipeline test completed. Check output at: $OUTPUT_DIR/unified_test"
else
    echo "Error: Unified pipeline not found at $UNIFIED_EXEC"
fi

echo ""
echo "Test summary:"
echo "============="
echo "Output directory: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR/"