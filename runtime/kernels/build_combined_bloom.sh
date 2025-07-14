#!/bin/bash
# build_combined_bloom.sh
# Builds a combined minimizer bloom filter from FQ resistance genes and AMR genes

# Default paths
FQ_RESISTANCE_FILE="../../data/resistance_db/reference_sequences.fasta"
AMR_CDS_FILE="../../data/AMR_CDS.fa"
OUTPUT_FILE="combined_minimizer_bloom.bin"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fq-resistance)
            FQ_RESISTANCE_FILE="$2"
            shift 2
            ;;
        --amr-cds)
            AMR_CDS_FILE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --fq-resistance FILE   Path to FQ resistance sequences (default: $FQ_RESISTANCE_FILE)"
            echo "  --amr-cds FILE         Path to AMR_CDS.fa file (default: $AMR_CDS_FILE)"
            echo "  --output FILE          Output bloom filter file (default: $OUTPUT_FILE)"
            echo ""
            echo "This script builds a combined minimizer bloom filter from both"
            echo "FQ resistance genes and AMR genes for efficient pre-screening."
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if files exist
if [ ! -f "$AMR_CDS_FILE" ]; then
    echo "Error: AMR_CDS.fa not found at: $AMR_CDS_FILE"
    exit 1
fi

if [ ! -f "$FQ_RESISTANCE_FILE" ]; then
    echo "Error: FQ resistance sequences not found at: $FQ_RESISTANCE_FILE"
    exit 1
fi

echo "Building combined minimizer bloom filter from:"
echo "  - FQ resistance genes: $FQ_RESISTANCE_FILE"
echo "  - AMR genes: $AMR_CDS_FILE"
echo "  - Output: $OUTPUT_FILE"
echo ""

# Build combined bloom filter from both FASTA files
./build_minimizer_bloom_filter "$OUTPUT_FILE" "$FQ_RESISTANCE_FILE" "$AMR_CDS_FILE"

if [ $? -eq 0 ]; then
    echo ""
    echo "Successfully created combined bloom filter: $OUTPUT_FILE"
    echo ""
    echo "To use this bloom filter in your pipeline:"
    echo "1. Copy to your database directory"
    echo "2. Update pipeline to load this pre-built bloom filter"
else
    echo "Error: Failed to build bloom filter"
    exit 1
fi