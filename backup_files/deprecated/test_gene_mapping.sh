#!/bin/bash
# test_gene_mapping.sh - Validate gene ID and family mapping

echo "=== Gene Mapping Validation Test ==="
echo "Testing gene ID and family mapping through the AMR pipeline"
echo

# Set paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="$SCRIPT_DIR"  # Use current directory for binaries

# Create test directories
TEST_DIR="$SCRIPT_DIR/test_gene_mapping"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# Create test FASTA files
cat > test_amr_cds.fa << EOF
>NC_000001|Beta-lactamase blaKPC-2
ATGTCACTGTATCGCCGTCTAGTTCTGCTGTCTTGTCTCTCATGGCCGCTGGCTGGCTTTTCTGCCACCGCGCTGACCAACCTCGTCGCGGAACCATTCGCTAAACTCGAACAGGACTTTGGCGGCTCCATCGGTGTGTACGCGATGGATACCGGCTCAGGCGCAACTGTAAGTTACCGCGCTGAGGAGCGCTTCCCACTGTGCAGCTCATTCAAGGGCTTTCTTGCTGCCGCTGTGCTGGCTCGCAGCCAGCAGCAGGCCGGCTTGCTGGACACACCCATCCGTTACGGCAAAAATGCGCTGGTTCCGTGGTCAC
>NC_000002|Fluoroquinolone resistance qnrA1
ATGGATATTATTGATAAAGTTTTTGGACAGCAAGAGATTGGCACCAAAATTATTTGCAACGAACGATGCCTGGTAGCTGCCATTAGTCAGGAGATCAACAAAAAGGTCATCCATGAAGATGCTCAGAGCAAATTCTCTGCCTGTTTTTCAGCAACAAGATGAAAACCTGGAAAATAAGGACCATGATAAGCTGTTCGATATTGAAGAACTGGCAAGACAGTACAGCGATTATCTGGAGCAGATTCTGAAACTGCATACCGGATGCAACGATCTGCGTTATCTGCAGCAGCATGCGATGTTTGGAGTTATTGGCAATGAAGGTAGAGAAGAAGATCCTTTTCATATTGATAAAAAAGAATCCCTGTCAAAGAAGGCAGATCTGGGGACAGATTTTCCAGCACTCGAAGGATTATGATACGTTTGGCGTGATTTTCAGTGATGACAAAAATGATCCGTTCCGTCGCGATCAGAAGGTTGAATAA
EOF

cat > test_amr_prot.fa << EOF
>0|WP_TEST001|1|1|blaKPC-2|Beta-lactamase KPC-2
MSLYRRLVLLSCLSWPLAGFSATALTNLVAEPFAKLEQDFGGSIGVYAMDTGSGATVSYRAEERFPLCSSFKGFLAAAVLARSQQQAGLLDTPIRYGKNALVRWSPE
>1|WP_TEST002|1|2|blaKPC-3|Beta-lactamase KPC-3  
MSLYRRLVLLSCLSWPLAGFSATALTNLVAEPFAKLEQDFGGSIGVYAMDTGSGATVSYRAEERFPLCSSFKGFLAAAVLARSQQQAGLLDTPIRYGKNALVRWSPA
>2|WP_TEST003|2|1|qnrA1|Quinolone resistance QnrA1
MFGQELKFYQKVKNSFNGATVLAAVVQRVCQKVLEHKYPNFSQECIQQLESQIQSLTNLADPTVQALRQLATGLVQVGKWFGVEPPGAPEGNPSAAELQGLLKTVFNWLASHTPLLKGLKIVGQVR
>3|WP_TEST004|3|1|vanA|Glycopeptide resistance VanA
MKKILFVGLFLACLLASVSAWSHPQFEKGLGKIENAVHIEAGLKDTTK
>4|WP_TEST005|4|1|mecA|Methicillin resistance protein
KKDLQFYRDDSYEKMKKQVEEMAKTLEEQTQKELEGALKQTVAKYQNELEHIQEKLNEQKSEKIEQEKSELKEQKNELQG
>5|WP_TEST006|5|1|aac6-Ib|Aminoglycoside acetyltransferase
MLNKVFDEFSVVQAGVSIDVDPAFVLGAARDVITLTGGLPEEEGEITLGWLQDSVIPGATVLIDGKTGKYVVSEAFSA
EOF

# Create test sample file
cat > test_samples.csv << EOF
Sample,FilePath,R1
test_sample,./,test_reads_R1.fastq
EOF

# Create test FASTQ files with reads matching our genes
cat > test_reads_R1.fastq << EOF
@read1_blaKPC-2_exact/1
ATGTCACTGTATCGCCGTCTAGTTCTGCTGTCTTGTCTCTCATGGCCGCTGGCTGGCTTTTCTGCCACCGCGCTGACCAACCTCGTCGCGGAACCATTCGCTAAACTCGAACAGGACTTTGGCGGCTCCATCGGTGTGTACGCGATGGATACCGGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_qnrA1_exact/1
ATGGATATTATTGATAAAGTTTTTGGACAGCAAGAGATTGGCACCAAAATTATTTGCAACGAACGATGCCTGGTAGCTGCCATTAGTCAGGAGATCAACAAAAAGGTCATCCATGAAGATGCTCAGAGCAAATTCTCTGCCTGTTTTTCAGCAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3_vanA_partial/1
ATGAAGAAAATACTATTTGTCGGCCTGTTTCTGGCCTGTCTGCTGGCCAGCGTGAGCGCCTGGTCGCATCCGCAATTCGAAAAAGGCCTGGGTAAGATCGAAAACGCTGTGCACATTGAAGCTGGCCTGAAAGACACGACCAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > test_reads_R2.fastq << EOF
@read1_blaKPC-2_exact/2
GTGACCACGGAACCAGCGCATTTTTGCCGTAACGGATGGGTGTGTCCAGCAAGCCGGCCTGCTGCTGGCTGCGAGCCAGCACAGCGGCAGCAAGAAAGCCCTTGAATGAGCTGCACAGTGGGAAGCGCTCCTCAGCGCGGTAACTTACAGTTGCGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_qnrA1_exact/2
GTTGCTGAAAAACAGGCAGAGAATTTGCTCTGAGCATCTTCATGGATGACCCTTTTTGTTGATCTCCTGACTAATGGCAGCTACCAGGCATCGTTCGTTGCAAATAATTTTGGTGCCAATCTCTTGCTGTCCAAAAACTTTATCAATAATATCCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3_vanA_partial/2
TTTGGTCGTGTCTTTCAGGCCAGCTTCAATGTGCACAGCGTTTTCGATCTTACCCAGGCCTTTTTCGAATTGCGGATGCGACCAGGCGCTCACGCTGGCCAGCAGACAGGCCAGAAACAGGCCGACAAATAGTATTTTTCTTCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Build test database
echo "Step 1: Building test protein database..."
"$BUILD_DIR/build_amr_protein_db" \
    test_amr_prot.fa \
    test_db/

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to build protein database"
    exit 1
fi

echo
echo "Checking database metadata..."
if [ -f "test_db/protein_details.json" ]; then
    echo "=== First few entries from protein_details.json ==="
    head -n 100 test_db/protein_details.json
fi

# Run detection with debug output
echo
echo "Step 2: Running AMR detection pipeline..."
"$BUILD_DIR/amr_detection" \
    "test_amr_cds.fa,test_amr_prot.fa" \
    "$TEST_DIR/test_samples.csv" \
    test_output \
    --protein-db test_db \
    --batch-size 100 \
    --validate > test_validation.log 2>&1

if [ $? -ne 0 ]; then
    echo "WARNING: AMR detection had issues (check log)"
fi

# Check results
echo
echo "Step 3: Checking results..."

echo
echo "=== Gene Mapping Validation from log ==="
grep -A20 "Gene Mapping Validation" test_validation.log || echo "Not found in log"

echo
echo "=== Gene entries from log ==="
grep -B2 -A2 "gene_family" test_validation.log | head -50

echo
echo "=== AMR Abundance TSV (if exists) ==="
if [ -f "test_output_amr_abundance.tsv" ]; then
    cat test_output_amr_abundance.tsv
else
    echo "File not found: test_output_amr_abundance.tsv"
fi

echo
echo "=== Gene Family Summary TSV (if exists) ==="
if [ -f "test_output_gene_family_summary.tsv" ]; then
    cat test_output_gene_family_summary.tsv
else
    echo "File not found: test_output_gene_family_summary.tsv"
fi

echo
echo "=== Clinical Report JSON (gene entries) ==="
if [ -f "test_output_clinical_amr_report.json" ]; then
    grep -A10 -B2 "gene_summaries" test_output_clinical_amr_report.json | head -50
else
    echo "File not found: test_output_clinical_amr_report.json"
fi

echo
echo "=== VALIDATION SUMMARY ==="
echo "Expected mappings:"
echo "  Gene ID 1: blaKPC-2 -> family: blaKPC"
echo "  Gene ID 1: blaKPC-3 -> family: blaKPC (same family as blaKPC-2)"
echo "  Gene ID 2: qnrA1 -> family: qnrA1"
echo "  Gene ID 3: vanA -> family: vanA"
echo "  Gene ID 4: mecA -> family: mecA"
echo "  Gene ID 5: aac6-Ib -> family: aac6"
echo
echo "Test complete. Full output in: test_validation.log"
echo "Test directory: $TEST_DIR"