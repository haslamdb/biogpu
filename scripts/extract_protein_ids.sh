#!/bin/bash
# Extract protein IDs from nucleotide records and download them

echo "Extracting accession numbers from DNA sequences..."
grep "^>" AMR_CDS_complete.fa | cut -d' ' -f1 | sed 's/>//' > dna_accessions.txt

echo "Found $(wc -l < dna_accessions.txt) DNA accessions"

# Get protein IDs from nucleotide records in batches
echo "Fetching protein IDs from nucleotide records..."
rm -f protein_ids.txt

batch_size=100
total=$(wc -l < dna_accessions.txt)
current=0

while IFS= read -r acc; do
    ((current++))
    echo -n "$acc" >> temp_batch.txt
    
    if [ $((current % batch_size)) -eq 0 ] || [ $current -eq $total ]; then
        echo "Processing batch ending at $current/$total..."
        # Get GenBank format to extract protein IDs
        cat temp_batch.txt | epost -db nuccore | \
        efetch -format gb | \
        grep -E "^     /protein_id=" | \
        sed 's/.*="//' | sed 's/"$//' >> protein_ids.txt
        
        rm -f temp_batch.txt
        sleep 1  # Be nice to NCBI
    else
        echo -n "," >> temp_batch.txt
    fi
done < dna_accessions.txt

echo "Found $(wc -l < protein_ids.txt) protein IDs"

# Now download the proteins
if [ -s protein_ids.txt ]; then
    echo "Downloading protein sequences..."
    cat protein_ids.txt | epost -db protein | efetch -format fasta > AMRProt_from_cds.fa
    echo "Downloaded $(grep -c '^>' AMRProt_from_cds.fa) protein sequences"
fi