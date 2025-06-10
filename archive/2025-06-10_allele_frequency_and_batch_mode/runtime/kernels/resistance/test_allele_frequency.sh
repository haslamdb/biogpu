#!/bin/bash

# Test script for running the resistance pipeline with allele frequency tracking
# This uses the 1M synthetic reads dataset

echo "============================================"
echo "Testing Resistance Pipeline with Allele Frequency Tracking"
echo "============================================"

# Set paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIOGPU_ROOT="$SCRIPT_DIR/../../.."
EXECUTABLE="$SCRIPT_DIR/build/clean_resistance_pipeline"

# Input data
NUCLEOTIDE_INDEX="$BIOGPU_ROOT/data/integrated_clean_db/nucleotide"
PROTEIN_DB="$BIOGPU_ROOT/data/integrated_clean_db/protein"
READS_R1="$BIOGPU_ROOT/1M_synthetic_reads_R1.fastq.gz"
READS_R2="$BIOGPU_ROOT/1M_synthetic_reads_R2.fastq.gz"
FQ_RESISTANCE_CSV="$BIOGPU_ROOT/data/quinolone_resistance_mutation_table.csv"

# Output directory
OUTPUT_DIR="$SCRIPT_DIR/results_allele_frequency_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"
OUTPUT_PREFIX="$OUTPUT_DIR/1M_synthetic_allele_freq"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please run ./build_resistance_pipeline.sh first"
    exit 1
fi

# Check if input files exist
echo "Checking input files..."
for file in "$NUCLEOTIDE_INDEX" "$PROTEIN_DB" "$READS_R1" "$READS_R2" "$FQ_RESISTANCE_CSV"; do
    if [ ! -e "$file" ]; then
        echo "Error: Input file not found: $file"
        exit 1
    fi
done
echo "All input files found!"

# Run the pipeline
echo ""
echo "Running resistance detection pipeline..."
echo "Command: $EXECUTABLE $NUCLEOTIDE_INDEX $PROTEIN_DB $READS_R1 $READS_R2 $OUTPUT_PREFIX $FQ_RESISTANCE_CSV"
echo ""

time "$EXECUTABLE" \
    "$NUCLEOTIDE_INDEX" \
    "$PROTEIN_DB" \
    "$READS_R1" \
    "$READS_R2" \
    "$OUTPUT_PREFIX" \
    "$FQ_RESISTANCE_CSV"

# Check if run was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
    echo ""
    echo "Output files generated:"
    echo "======================================"
    
    # List all output files
    ls -lh "$OUTPUT_DIR"/*
    
    echo ""
    echo "Key files:"
    echo "----------"
    echo "1. Main results:        $(basename $OUTPUT_PREFIX).json"
    echo "2. Protein matches:     $(basename $OUTPUT_PREFIX)_protein_matches.csv"
    echo "3. Allele frequencies:  $(basename $OUTPUT_PREFIX)_allele_frequencies.csv"
    echo "4. Clinical HTML report:$(basename $OUTPUT_PREFIX)_clinical_fq_report.html"
    echo "5. HDF5 alignments:     $(basename $OUTPUT_PREFIX).h5"
    
    # Show a preview of allele frequency data if it exists
    ALLELE_FREQ_FILE="$OUTPUT_PREFIX"_allele_frequencies.csv
    if [ -f "$ALLELE_FREQ_FILE" ]; then
        echo ""
        echo "Preview of allele frequency data:"
        echo "================================="
        head -n 10 "$ALLELE_FREQ_FILE" | column -t -s,
        
        # Count lines in the file
        FREQ_COUNT=$(tail -n +2 "$ALLELE_FREQ_FILE" | wc -l)
        echo ""
        echo "Total positions analyzed: $FREQ_COUNT"
    fi
    
    # Open HTML report in browser if available
    HTML_REPORT="$OUTPUT_PREFIX"_clinical_fq_report.html
    if [ -f "$HTML_REPORT" ]; then
        echo ""
        echo "Opening HTML report in browser..."
        if command -v xdg-open > /dev/null; then
            xdg-open "$HTML_REPORT" 2>/dev/null &
        elif command -v open > /dev/null; then
            open "$HTML_REPORT"
        else
            echo "Please open $HTML_REPORT manually in your browser"
        fi
    fi
else
    echo ""
    echo "Pipeline failed! Check error messages above."
    exit 1
fi

echo ""
echo "Test complete!"