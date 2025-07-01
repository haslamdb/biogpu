#!/bin/bash
# Direct download of PRJNA313047 protein sequences
set -e

EMAIL="dbhaslam@gmail.com"
BATCH_SIZE=500  # Download in batches to avoid timeouts

echo "Downloading PRJNA313047 protein sequences directly..."

# Get total count
echo "Getting protein sequence count..."
COUNT=$(esearch -db protein -query "PRJNA313047" | xtract -pattern ENTREZ_DIRECT -element Count)
echo "Total protein sequences: $COUNT"

# Download in batches
echo "Downloading sequences in batches of $BATCH_SIZE..."
rm -f AMRProt_complete.fa

for start in $(seq 0 $BATCH_SIZE $COUNT); do
    echo "Downloading batch starting at $start..."
    esearch -db protein -query "PRJNA313047" | \
    efetch -format fasta -start $start -stop $((start + BATCH_SIZE - 1)) >> AMRProt_complete.fa
    
    # Show progress
    current_count=$(grep -c "^>" AMRProt_complete.fa 2>/dev/null || echo 0)
    echo "Progress: $current_count / $COUNT sequences downloaded"
    
    # Small delay to be nice to NCBI
    sleep 1
done

echo "Download complete!"
echo "Total sequences: $(grep -c '^>' AMRProt_complete.fa)"
echo "File size: $(ls -lh AMRProt_complete.fa | awk '{print $5}')"