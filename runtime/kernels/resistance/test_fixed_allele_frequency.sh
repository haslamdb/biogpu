#!/bin/bash

# Test script for fixed allele frequency calculation
echo "Testing fixed allele frequency calculation..."
echo "============================================"

# Set up paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR

# Use the correct paths
NUCLEOTIDE_INDEX="/home/david/Documents/Code/biogpu/data/integrated_clean_db"
PROTEIN_DB="/home/david/Documents/Code/biogpu/data/integrated_clean_db"
TEST_DATA="/home/david/Documents/Code/biogpu/1M_synthetic_reads_R1.fastq.gz"
OUTPUT_PREFIX="results_fixed_allele_freq/test_fixed_allele_freq"
FQ_CSV="/home/david/Documents/Code/biogpu/data/quinolone_resistance_mutation_table.csv"

# Create output directory
mkdir -p results_fixed_allele_freq

# Run the pipeline with the fixed allele frequency calculation
echo "Running resistance detection pipeline..."
./build/clean_resistance_pipeline \
    $NUCLEOTIDE_INDEX \
    $PROTEIN_DB \
    $TEST_DATA \
    $TEST_DATA \
    --output-prefix $OUTPUT_PREFIX \
    --fq-csv $FQ_CSV

# Check if the allele frequency CSV was created
if [ -f "${OUTPUT_PREFIX}_allele_frequencies.csv" ]; then
    echo ""
    echo "Allele frequency CSV created. Checking results..."
    echo ""
    
    # Display the first few lines focusing on position 87
    echo "Allele frequencies at key positions (including position 87):"
    echo "------------------------------------------------------------"
    head -1 "${OUTPUT_PREFIX}_allele_frequencies.csv"
    grep -E "gyrA,87|gyrA,83|parC,80" "${OUTPUT_PREFIX}_allele_frequencies.csv" | head -10
    
    echo ""
    echo "Full results saved to: ${OUTPUT_PREFIX}_allele_frequencies.csv"
    echo "Clinical report: ${OUTPUT_PREFIX}_clinical_fq_report.html"
else
    echo "ERROR: Allele frequency CSV was not created!"
    exit 1
fi

echo ""
echo "Test complete!"