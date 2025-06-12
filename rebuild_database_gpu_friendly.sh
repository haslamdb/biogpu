#!/bin/bash
# Rebuild the pathogen database with GPU-friendly memory settings

echo "Rebuilding pathogen database for GPU with limited memory..."

# First, let's check how many k-mers we have
KMER_COUNT=$(grep -v "^#" data/pathogen_profiler_db/database_kmers.txt | wc -l)
echo "Total k-mers in database: $KMER_COUNT"

# Calculate a reasonable table size for 12GB GPU
# Each entry is 16 bytes (8 byte hash + 4 byte taxon + 4 byte padding)
# For 10GB usable memory: 10GB / 16 bytes = ~670M entries
# With bucket_size=8, we need table_size = 670M / 8 = ~84M buckets
# Let's use 64M buckets (2^26) for safety = 8GB

# We need to modify the database building code to accept a max table size parameter
# For now, let's create a smaller k-mer file for testing

echo "Creating a test subset of k-mers..."
head -n 1000000 data/pathogen_profiler_db/database_kmers.txt > data/pathogen_profiler_db/database_kmers_subset.txt

# Build the database from the subset
echo "Building GPU database from k-mer subset..."
./runtime/kernels/profiler/build/build_db_from_kmers \
    data/pathogen_profiler_db/database_kmers_subset.txt \
    data/pathogen_profiler_db/pathogen_subset.db

echo "Done! Created subset database for testing."