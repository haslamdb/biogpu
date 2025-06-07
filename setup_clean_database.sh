#!/bin/bash
# setup_clean_database.sh
# Build complete database structure for integrated pipeline

set -e  # Exit on error

# Configuration
BASE_DIR="data"
WILDTYPE_DIR="$BASE_DIR/wildtype_protein_seqs"
MUTATIONS_CSV="$BASE_DIR/Known_Quinolone_Changes.csv"
OUTPUT_DIR="$BASE_DIR/integrated_clean_db"

echo "=== Building Clean Integrated Database ==="
echo "Source: $WILDTYPE_DIR"
echo "Output: $OUTPUT_DIR"
echo ""

# Check if wildtype directory exists
if [ ! -d "$WILDTYPE_DIR" ]; then
    echo "ERROR: Wildtype directory not found: $WILDTYPE_DIR"
    echo "Please ensure your wildtype protein sequences are in this directory"
    exit 1
fi

# Clean previous build if exists
if [ -d "$OUTPUT_DIR" ]; then
    echo "Removing previous database..."
    rm -rf "$OUTPUT_DIR"
fi

# Build the database
echo "Building database..."
python3 src/python/build_clean_dynamic_database.py \
    "$WILDTYPE_DIR" \
    "$OUTPUT_DIR" \
    --mutations-csv "$MUTATIONS_CSV" \
    --kmer-length 5

# Create empty nucleotide directory (for now)
echo ""
echo "Creating nucleotide index directory..."
mkdir -p "$OUTPUT_DIR/nucleotide"
echo "Note: Nucleotide index not implemented yet - creating empty directory"

# Create symlinks for backward compatibility
echo ""
echo "Creating symlinks..."
cd "$BASE_DIR"

# Remove old symlinks if they exist
rm -f fq_resistance_index protein_resistance_db resistance_db

# Create new symlinks
ln -s integrated_clean_db/nucleotide fq_resistance_index
ln -s integrated_clean_db/protein protein_resistance_db
ln -s integrated_clean_db/resistance resistance_db

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Database structure:"
echo "  $OUTPUT_DIR/"
echo "  ├── nucleotide/      (currently empty - k-mer index placeholder)"
echo "  ├── protein/         (protein sequences and 5-mer index)"
echo "  ├── resistance/      (mutation positions only - no drug info)"
echo "  ├── database_mappings.json"
echo "  ├── gene_mappings.tsv"
echo "  └── species_mappings.tsv"
echo ""
echo "Symlinks created:"
echo "  fq_resistance_index -> integrated_clean_db/nucleotide/"
echo "  protein_resistance_db -> integrated_clean_db/protein/"
echo "  resistance_db -> integrated_clean_db/resistance/"
echo ""
echo "To run the pipeline:"
echo "  ./integrated_resistance_pipeline \\"
echo "    $BASE_DIR/fq_resistance_index \\"
echo "    $BASE_DIR/protein_resistance_db \\"
echo "    $BASE_DIR/resistance_db \\"
echo "    test_R1.fastq.gz test_R2.fastq.gz output_prefix"