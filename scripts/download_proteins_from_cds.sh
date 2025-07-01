#!/bin/bash
# Extract and download proteins from CDS records

echo "Extracting DNA accessions..."
grep "^>" AMR_CDS_complete.fa | cut -d' ' -f1 | sed 's/>//' > dna_accessions.txt
total_dna=$(wc -l < dna_accessions.txt)
echo "Found $total_dna DNA accessions"

# Process in batches to get protein IDs
echo "Extracting protein IDs from GenBank records..."
rm -f all_protein_ids.txt

batch_size=200
for ((i=0; i<total_dna; i+=batch_size)); do
    end=$((i+batch_size))
    if [ $end -gt $total_dna ]; then
        end=$total_dna
    fi
    
    echo "Processing DNA records $((i+1)) to $end..."
    
    # Get batch of accessions
    sed -n "$((i+1)),$end p" dna_accessions.txt > batch_acc.txt
    
    # Post them and fetch GenBank format to extract protein IDs
    cat batch_acc.txt | \
    epost -db nuccore -format acc | \
    efetch -format gb 2>/dev/null | \
    grep "/protein_id=" | \
    sed 's/.*="//' | sed 's/"$//' >> all_protein_ids.txt
    
    # Show progress
    current_proteins=$(wc -l < all_protein_ids.txt)
    echo "  Extracted $current_proteins protein IDs so far..."
    
    sleep 2  # Be nice to NCBI
done

# Remove duplicates
sort -u all_protein_ids.txt > unique_protein_ids.txt
total_proteins=$(wc -l < unique_protein_ids.txt)
echo "Found $total_proteins unique protein IDs"

# Download proteins in batches
echo "Downloading protein sequences..."
rm -f AMRProt_complete_from_cds.fa

batch_size=500
for ((i=0; i<total_proteins; i+=batch_size)); do
    end=$((i+batch_size))
    if [ $end -gt $total_proteins ]; then
        end=$total_proteins
    fi
    
    echo "Downloading proteins $((i+1)) to $end..."
    
    # Get batch of protein IDs
    sed -n "$((i+1)),$end p" unique_protein_ids.txt > batch_proteins.txt
    
    # Download them
    cat batch_proteins.txt | \
    epost -db protein -format acc | \
    efetch -format fasta >> AMRProt_complete_from_cds.fa
    
    # Show progress
    current_count=$(grep -c "^>" AMRProt_complete_from_cds.fa 2>/dev/null || echo 0)
    echo "  Downloaded $current_count proteins so far..."
    
    sleep 1
done

echo "Download complete!"
echo "Total protein sequences: $(grep -c '^>' AMRProt_complete_from_cds.fa)"
echo "File size: $(ls -lh AMRProt_complete_from_cds.fa | awk '{print $5}')"

# Clean up
rm -f dna_accessions.txt batch_acc.txt all_protein_ids.txt batch_proteins.txt