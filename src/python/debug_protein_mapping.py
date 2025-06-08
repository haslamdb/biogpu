#!/usr/bin/env python3
"""
Debug protein ID mapping between k-mer index and protein sequences
"""

import struct
import json

def debug_protein_mapping(protein_db_path):
    # Load protein details
    with open(f"{protein_db_path}/protein_details.json", 'r') as f:
        protein_details = json.load(f)
    
    print(f"Number of proteins in JSON: {len(protein_details)}")
    print("\nProtein order in JSON:")
    for i, p in enumerate(protein_details):
        print(f"  Index {i}: ID={p['id']}, {p['species']} {p['gene']}, length={p['length']}")
    
    # Check k-mer index
    with open(f"{protein_db_path}/protein_kmers.bin", 'rb') as f:
        kmer_length = struct.unpack('I', f.read(4))[0]
        num_kmers = struct.unpack('I', f.read(4))[0]
        
        print(f"\nK-mer index: {num_kmers} unique {kmer_length}-mers")
        
        # Sample first few k-mers to see protein indices
        print("\nFirst 5 k-mers and their protein indices:")
        for i in range(min(5, num_kmers)):
            kmer = f.read(kmer_length).decode('ascii')
            num_positions = struct.unpack('I', f.read(4))[0]
            
            protein_indices = set()
            for j in range(num_positions):
                protein_idx, pos = struct.unpack('II', f.read(8))
                protein_indices.add(protein_idx)
            
            print(f"  K-mer '{kmer}': found in proteins {sorted(protein_indices)}")

if __name__ == "__main__":
    debug_protein_mapping("data/integrated_clean_db/protein")