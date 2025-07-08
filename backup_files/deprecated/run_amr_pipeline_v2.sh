#!/bin/bash

# Example script to run the refactored AMR detection pipeline

# Set environment variables for database paths
export BIOGPU_DATA="${BIOGPU_DATA:-/data/biogpu}"

# Create output directory
mkdir -p output

# Run AMR detection on test data
echo "Running AMR detection pipeline..."

# Single-end mode
if [ -f "../data/569_A_038_R1.fastq.gz" ]; then
    echo "Processing single-end FASTQ file..."
    ./test_amr_pipeline_v2 ../data/569_A_038_R1.fastq.gz
else
    echo "Test file not found: ../data/569_A_038_R1.fastq.gz"
fi

# With custom configuration
if [ -f "../../common/config/example_config.json" ]; then
    echo -e "\nProcessing with custom configuration..."
    ./test_amr_pipeline_v2 ../data/569_A_038_R1.fastq.gz ../../common/config/example_config.json
fi

echo "Done!"