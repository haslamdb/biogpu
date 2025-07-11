#!/bin/bash

# Extract first ~25 genomes from the concatenated FNA file
# We'll get genomes from different taxids for diversity

input_file="$HOME/Data/NewK2Fasta/none/library/bacteria/library_species_taxid.fna"
output_file="test_concatenated_genomes.fna"

echo "Extracting test genomes from $input_file"

# Extract lines for first 25 different taxids
# Using line numbers from the grep output above
sed -n '1,461320p' "$input_file" > "$output_file"

# Count how many genomes we got
genome_count=$(grep -c "^>" "$output_file")
echo "Extracted $genome_count genome sequences to $output_file"

# Show file size
ls -lh "$output_file"

# Show the taxids we captured
echo "Taxids in test file:"
grep "^>" "$output_file" | sed 's/.*taxid|//' | sed 's/|.*//' | sort -u | head -20