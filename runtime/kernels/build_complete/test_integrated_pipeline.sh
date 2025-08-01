#!/bin/bash
# Test script for integrated pipeline

echo "BioGPU Integrated Pipeline Test"
echo "==============================="

# Create test directories
mkdir -p /tmp/biogpu_test/{input,output}

# Create dummy test files (empty for validation test)
touch /tmp/biogpu_test/input/test_R1.fastq.gz
touch /tmp/biogpu_test/input/test_R2.fastq.gz

# Test help
echo -e "\n1. Testing help..."
./bio_gpu_pipeline_integrated --help

# Test validation (will fail on empty files, but tests argument parsing)
echo -e "\n2. Testing argument validation..."
./bio_gpu_pipeline_integrated \
    --r1 /tmp/biogpu_test/input/test_R1.fastq.gz \
    --r2 /tmp/biogpu_test/input/test_R2.fastq.gz \
    --output-dir /tmp/biogpu_test/output \
    --reference-db /tmp/ref_db \
    --resistance-db /tmp/res_db \
    --sample-id TEST001 \
    --progress-json 2>&1 | head -5

echo -e "\nTest complete!"
