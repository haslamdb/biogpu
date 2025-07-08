#!/bin/bash
# Test script to validate gene_id and gene_family propagation

echo "=== Gene Family Validation Test ==="
echo "Testing gene_id and gene_family propagation through AMR pipeline"
echo

# Set paths
TEST_DIR="/home/david/Documents/Code/biogpu/runtime/kernels/genes"
BUILD_DIR="/home/david/Documents/Code/biogpu/build"
TEST_FASTA="$TEST_DIR/test_gene_family_validation.fasta"
TEST_DB_DIR="$TEST_DIR/test_amr_db"
TEST_OUTPUT_DIR="$TEST_DIR/test_validation_output"

# Clean up previous test runs
rm -rf "$TEST_DB_DIR" "$TEST_OUTPUT_DIR"
mkdir -p "$TEST_DB_DIR" "$TEST_OUTPUT_DIR"

# Step 1: Build the test protein database
echo "Step 1: Building test AMR protein database..."
"$BUILD_DIR/runtime/kernels/genes/build_amr_protein_database" \
    "$TEST_FASTA" \
    "$TEST_DB_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to build protein database"
    exit 1
fi

echo
echo "Database built successfully. Checking metadata..."
if [ -f "$TEST_DB_DIR/protein_details.json" ]; then
    echo "First 5 entries from protein_details.json:"
    head -n 50 "$TEST_DB_DIR/protein_details.json" | grep -E "gene_id|gene_family|gene_name" | head -20
fi

# Step 2: Create test reads (sequences from the genes with mutations)
echo
echo "Step 2: Creating test FASTQ file..."
cat > "$TEST_OUTPUT_DIR/test_reads.fastq" << 'EOF'
@read1_blaKPC-2_match
ATGTCACTGTATCGCCGTCTAGTTCTGCTGTCTTGTCTCTCATGGCCGCTGGCTGGCTTTTCTGCCACCGCGCTGACCAACCTCGTCGCGGAACCATTCGCTAAACTCGAACAGGACTTTGGCGGCTCCATCGGTGTGTACGCGATGGATACCGGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_qnrA1_match  
ATGTTCGGGCAAGAGCTGAAGTTCTACCAGAAAGTGAAAAACAGCTTTAATGGCGCAACGGTGCTGGCAGCGGTGGTGCAGCGGGTTTGCCAGAAAGTTCTGGAACATAAATATCCAAATTTTTCGCAAGAATGCATTCAGCAGCTGGAAAGTCAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3_vanA_partial
ATGAAGAAAATACTATTTGTCGGCCTGTTTCTGGCCTGTCTGCTGGCCAGCGTGAGCGCCTGGTCGCATCCGCAATTCGAAAAAGGCCTGGGTAAGATCGAAAACGCTGTGCACATTGAAGCTGGCCTGAAAGACACGACCAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read4_mecA_match
AAGAAAGACCTGCAGTTCTATCGCGATGATTCATATGAAAAAATGAAGAAACAGGTTGAAGAAATGGCTAAAACGCTGGAAGAGCAGACGCAGAAAGAGCTGGAAGGGGCCCTGAAACAGACCGTGGCGAAATATCAGAACGAACTGGAACACATT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read5_negative_control
ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Step 3: Run AMR detection pipeline
echo
echo "Step 3: Running AMR detection pipeline..."
"$BUILD_DIR/runtime/kernels/genes/amr_detection_main" \
    --reads "$TEST_OUTPUT_DIR/test_reads.fastq" \
    --protein-db "$TEST_DB_DIR" \
    --output "$TEST_OUTPUT_DIR/test_results" \
    --threads 1 \
    --batch-size 100

if [ $? -ne 0 ]; then
    echo "ERROR: AMR detection pipeline failed"
    exit 1
fi

# Step 4: Check outputs
echo
echo "Step 4: Validating outputs..."
echo

# Check TSV output
if [ -f "$TEST_OUTPUT_DIR/test_results_amr_abundance.tsv" ]; then
    echo "=== AMR Abundance TSV Output ==="
    cat "$TEST_OUTPUT_DIR/test_results_amr_abundance.tsv"
    echo
fi

# Check hits TSV
if [ -f "$TEST_OUTPUT_DIR/test_results_hits.tsv" ]; then
    echo "=== AMR Hits TSV Output ==="
    echo "First 10 lines:"
    head -10 "$TEST_OUTPUT_DIR/test_results_hits.tsv"
    echo
fi

# Check JSON output
if [ -f "$TEST_OUTPUT_DIR/test_results_clinical_amr_report.json" ]; then
    echo "=== Gene entries from JSON report ==="
    grep -A 5 -B 5 "gene_name\|gene_family\|gene_id" "$TEST_OUTPUT_DIR/test_results_clinical_amr_report.json" | head -50
    echo
fi

echo
echo "=== Validation Summary ==="
echo "Expected gene mappings:"
echo "  gene_id=1 -> gene_name=blaKPC-2, gene_family=blaKPC"
echo "  gene_id=2 -> gene_name=qnrA1, gene_family=qnrA1"
echo "  gene_id=3 -> gene_name=vanA, gene_family=vanA"
echo "  gene_id=4 -> gene_name=mecA, gene_family=mecA"
echo "  gene_id=5 -> gene_name=aac6-Ib, gene_family=aac6"
echo
echo "Check the outputs above to verify correct propagation!"