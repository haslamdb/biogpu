#!/usr/bin/env python3
"""Check if synthetic reads contain sequences from resistance genes"""

import gzip
import json
import os
from collections import defaultdict

def load_gene_sequences():
    """Load gene sequences from the output directory"""
    gene_sequences = defaultdict(list)
    
    output_dir = "output/GeneFiles"
    for species in os.listdir(output_dir):
        species_dir = os.path.join(output_dir, species)
        if os.path.isdir(species_dir):
            for gene_file in os.listdir(species_dir):
                if gene_file.endswith('.json'):
                    gene_name = gene_file.replace('.json', '')
                    filepath = os.path.join(species_dir, gene_file)
                    
                    with open(filepath, 'r') as f:
                        data = json.load(f)
                        # Data is a list of entries
                        for entry in data:
                            if 'gene_features' in entry:
                                for gene_feature in entry['gene_features']:
                                    if 'nucleotide_sequence' in gene_feature:
                                        gene_sequences[gene_name].append({
                                            'species': species.replace('_', ' '),
                                            'sequence': gene_feature['nucleotide_sequence'],
                                            'length': len(gene_feature['nucleotide_sequence']),
                                            'accession': entry.get('accession', 'unknown')
                                        })
    
    return gene_sequences

def extract_kmers(sequence, k=15):
    """Extract all k-mers from a sequence"""
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(sequence[i:i+k])
    return kmers

def check_synthetic_reads(fastq_file, gene_sequences, num_reads=1000):
    """Check if synthetic reads match any gene sequences"""
    print(f"\nChecking {fastq_file}...")
    
    # Build k-mer index from gene sequences
    gene_kmers = defaultdict(set)
    total_gene_kmers = 0
    
    for gene_name, sequences in gene_sequences.items():
        for seq_info in sequences:
            seq = seq_info['sequence']
            kmers = extract_kmers(seq)
            gene_kmers[gene_name].update(kmers)
            total_gene_kmers += len(kmers)
    
    print(f"Built k-mer index: {len(gene_kmers)} genes, {total_gene_kmers} total k-mers")
    
    # Process reads
    matches = defaultdict(int)
    read_count = 0
    matching_reads = 0
    
    with gzip.open(fastq_file, 'rt') as f:
        while read_count < num_reads:
            # Read FASTQ entry
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            
            read_count += 1
            
            # Extract k-mers from read
            read_kmers = extract_kmers(sequence)
            
            # Check against gene k-mers
            read_matched = False
            for gene_name, gene_kmer_set in gene_kmers.items():
                common_kmers = read_kmers & gene_kmer_set
                if common_kmers:
                    matches[gene_name] += len(common_kmers)
                    read_matched = True
                    
                    # Show details for first few matches
                    if matching_reads < 5:
                        print(f"\nRead {read_count} matches {gene_name}:")
                        print(f"  Sequence: {sequence[:50]}...")
                        print(f"  Common k-mers: {len(common_kmers)}")
                        print(f"  Example k-mer: {list(common_kmers)[0]}")
            
            if read_matched:
                matching_reads += 1
    
    print(f"\n=== Summary ===")
    print(f"Total reads checked: {read_count}")
    print(f"Reads with matches: {matching_reads} ({matching_reads/read_count*100:.1f}%)")
    print(f"\nTop gene matches:")
    for gene, count in sorted(matches.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {gene}: {count} k-mer matches")
    
    return matches

def main():
    # Load gene sequences
    print("Loading gene sequences from output directory...")
    gene_sequences = load_gene_sequences()
    print(f"Loaded {len(gene_sequences)} genes")
    
    # Show sample genes
    print("\nSample genes loaded:")
    for gene, sequences in list(gene_sequences.items())[:5]:
        print(f"  {gene}: {len(sequences)} sequences")
        if sequences:
            print(f"    First sequence: {sequences[0]['species']}, {sequences[0]['length']}bp")
    
    # Check synthetic reads
    r1_file = "synthetic_reads_20250529_R1.fastq.gz"
    if os.path.exists(r1_file):
        matches = check_synthetic_reads(r1_file, gene_sequences)
    else:
        print(f"\nError: {r1_file} not found")
        
        # Try to find any synthetic reads
        print("\nLooking for other synthetic read files...")
        for f in os.listdir('.'):
            if 'synthetic' in f and f.endswith('.fastq.gz'):
                print(f"  Found: {f}")

if __name__ == "__main__":
    main()