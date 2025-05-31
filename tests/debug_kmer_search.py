#!/usr/bin/env python3
"""Debug k-mer search functionality"""

import json
import os
import struct
from collections import defaultdict

def load_binary_kmers(kmers_file):
    """Load k-mers from binary format"""
    kmers = {}
    with open(kmers_file, 'rb') as f:
        # Read header
        magic = f.read(4)
        version = struct.unpack('I', f.read(4))[0]
        kmer_size = struct.unpack('I', f.read(4))[0]
        num_kmers = struct.unpack('I', f.read(4))[0]
        
        print(f"Binary k-mer file info:")
        print(f"  Magic: {magic}")
        print(f"  Version: {version}")
        print(f"  K-mer size: {kmer_size}")
        print(f"  Number of k-mers: {num_kmers}")
        
        # Read k-mers
        for i in range(min(num_kmers, 1000000)):  # Limit for debug
            kmer_bytes = f.read(kmer_size)
            if not kmer_bytes:
                break
            kmer = kmer_bytes.decode('ascii', errors='ignore').rstrip('\x00')
            
            # Read associated data (gene count, etc)
            gene_count = struct.unpack('I', f.read(4))[0]
            
            if i < 5:  # Show first few
                print(f"  K-mer {i}: {kmer} (genes: {gene_count})")
            
            kmers[kmer] = {'count': gene_count}
            
            # Skip gene data for now
            f.seek(gene_count * 8, 1)  # Assuming 8 bytes per gene entry
    
    return kmers

def extract_kmers(sequence, k=15):
    """Extract all k-mers from a sequence"""
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(sequence[i:i+k])
    return kmers

def debug_kmer_search():
    # Load synthetic reads
    with open('synthetic_fq_results.json', 'r') as f:
        synthetic_data = json.load(f)
    
    # Load k-mer index from binary format
    kmers_path = 'data/gpu_resistance_db/kmers.bin'
    if not os.path.exists(kmers_path):
        print(f"Error: K-mer binary file not found at {kmers_path}")
        return
    
    print(f"Loading k-mer index from {kmers_path}...")
    kmer_index = load_binary_kmers(kmers_path)
    
    # Check index structure
    print(f"\nK-mer index type: {type(kmer_index)}")
    if isinstance(kmer_index, dict):
        print(f"Number of k-mers in index: {len(kmer_index)}")
        # Show sample k-mers
        sample_kmers = list(kmer_index.keys())[:5]
        print(f"Sample k-mers: {sample_kmers}")
        for kmer in sample_kmers:
            print(f"  {kmer}: {kmer_index[kmer]}")
    
    # Test with first few synthetic reads
    print("\n=== Testing k-mer matches for synthetic reads ===")
    for i, read_data in enumerate(synthetic_data[:5]):
        read_seq = read_data['sequence']
        print(f"\nRead {i+1} ({read_data['mutation_type']}):")
        print(f"  Species: {read_data['species']}")
        print(f"  Gene: {read_data['gene']}")
        print(f"  Length: {len(read_seq)}bp")
        print(f"  First 50bp: {read_seq[:50]}...")
        
        # Extract k-mers from read
        read_kmers = extract_kmers(read_seq)
        print(f"  Total k-mers in read: {len(read_kmers)}")
        
        # Check how many match the index
        matching_kmers = 0
        matched_genes = defaultdict(int)
        
        for kmer in read_kmers:
            if kmer in kmer_index:
                matching_kmers += 1
                # For binary format, we just have counts
                if isinstance(kmer_index[kmer], dict):
                    matched_genes['unknown'] += 1
        
        print(f"  K-mers matching index: {matching_kmers} ({matching_kmers/len(read_kmers)*100:.1f}%)")
        print(f"  Matched genes: {dict(list(matched_genes.items())[:5])}")
        
        # Show specific k-mer examples
        print(f"  First 3 k-mers from read:")
        for j, kmer in enumerate(list(read_kmers)[:3]):
            in_index = "✓" if kmer in kmer_index else "✗"
            print(f"    {kmer} {in_index}")

    # For binary format, we can't check specific genes easily
    print(f"\n=== K-mer index summary ===")
    print(f"Total unique k-mers in index: {len(kmer_index)}")

if __name__ == "__main__":
    debug_kmer_search()