#!/usr/bin/env python3
"""Debug k-mer search with proper binary format understanding"""

import json
import struct
import os

# Constants to match the CUDA code
KMER_LENGTH = 15
BASE_A = 0
BASE_C = 1
BASE_G = 2
BASE_T = 3
BASE_N = 4

def encode_base(base):
    """Encode base to 2-bit representation"""
    base = base.upper()
    if base == 'A': return BASE_A
    elif base == 'C': return BASE_C
    elif base == 'G': return BASE_G
    elif base == 'T': return BASE_T
    else: return BASE_N

def encode_kmer(seq):
    """Encode k-mer string to 64-bit integer"""
    if len(seq) != KMER_LENGTH:
        return None
    
    kmer = 0
    for base in seq:
        encoded = encode_base(base)
        if encoded == BASE_N:
            return None  # Invalid k-mer
        kmer = (kmer << 2) | encoded
    return kmer

def decode_kmer(kmer_int):
    """Decode 64-bit integer to k-mer string"""
    bases = []
    for i in range(KMER_LENGTH):
        base_code = (kmer_int >> (2 * (KMER_LENGTH - 1 - i))) & 3
        if base_code == BASE_A: bases.append('A')
        elif base_code == BASE_C: bases.append('C')
        elif base_code == BASE_G: bases.append('G')
        elif base_code == BASE_T: bases.append('T')
    return ''.join(bases)

def extract_kmers_encoded(sequence):
    """Extract all k-mers from sequence as encoded integers"""
    kmers = set()
    for i in range(len(sequence) - KMER_LENGTH + 1):
        kmer_str = sequence[i:i+KMER_LENGTH]
        encoded = encode_kmer(kmer_str)
        if encoded is not None:
            kmers.add(encoded)
    return kmers

def load_binary_index(index_path):
    """Load the binary k-mer index"""
    print(f"Loading binary k-mer index from {index_path}")
    
    kmers_dict = {}
    
    with open(index_path, 'rb') as f:
        # Read header
        num_kmers = struct.unpack('I', f.read(4))[0]
        kmer_size = struct.unpack('I', f.read(4))[0]
        
        print(f"Number of k-mers in file: {num_kmers}")
        print(f"K-mer size: {kmer_size}")
        
        # Read all k-mers
        for i in range(num_kmers):
            # Read encoded k-mer (64-bit)
            kmer_encoded = struct.unpack('Q', f.read(8))[0]
            
            # Read gene count
            gene_count = struct.unpack('I', f.read(4))[0]
            
            if i < 5:  # Show first few
                print(f"\nK-mer {i}:")
                print(f"  Encoded value: {kmer_encoded}")
                print(f"  Decoded: {decode_kmer(kmer_encoded)}")
                print(f"  Gene count: {gene_count}")
            
            gene_entries = []
            
            # Read gene data
            for g in range(gene_count):
                organism_id = struct.unpack('I', f.read(4))[0]
                gene_name_len = struct.unpack('I', f.read(4))[0]
                gene_name = f.read(gene_name_len).decode('utf-8', errors='replace')
                sequence_id = struct.unpack('I', f.read(4))[0]
                
                gene_entries.append({
                    'organism_id': organism_id,
                    'gene': gene_name,
                    'sequence_id': sequence_id
                })
                
                if i < 2 and g < 2:  # Show details for first couple
                    print(f"    Gene {g}: organism={organism_id}, gene={gene_name}, seq={sequence_id}")
            
            kmers_dict[kmer_encoded] = gene_entries
            
            if i % 10000 == 0 and i > 0:
                print(f"  Loaded {i} k-mers...")
    
    print(f"\nTotal k-mers loaded: {len(kmers_dict)}")
    return kmers_dict

def debug_kmer_matching():
    # Load synthetic reads
    with open('synthetic_fq_results.json', 'r') as f:
        data = json.load(f)
        synthetic_data = data['mutations'] if 'mutations' in data else data
    
    # Load k-mer index
    index_path = 'data/gpu_resistance_db/kmers.bin'
    if not os.path.exists(index_path):
        print(f"Error: K-mer binary file not found at {index_path}")
        return
        
    print("Loading k-mer index (this may take a moment)...")
    kmer_index = load_binary_index(index_path)
    
    print("\n=== Testing k-mer matching for synthetic reads ===")
    for i, read_data in enumerate(synthetic_data[:5]):
        read_seq = read_data['sequence']
        print(f"\nRead {i+1} ({read_data['mutation_type']}):")
        print(f"  Species: {read_data['species']}")
        print(f"  Gene: {read_data['gene']}")
        print(f"  Length: {len(read_seq)}bp")
        print(f"  First 50bp: {read_seq[:50]}...")
        
        # Extract k-mers as encoded integers
        read_kmers = extract_kmers_encoded(read_seq)
        print(f"  Total valid k-mers: {len(read_kmers)}")
        
        # Check how many match
        matching_kmers = 0
        matched_genes = {}
        
        for kmer in read_kmers:
            if kmer in kmer_index:
                matching_kmers += 1
                # Track which genes match
                for entry in kmer_index[kmer]:
                    gene_key = f"org{entry['organism_id']}:{entry['gene']}"
                    matched_genes[gene_key] = matched_genes.get(gene_key, 0) + 1
        
        print(f"  K-mers matching index: {matching_kmers} ({matching_kmers/len(read_kmers)*100:.1f}%)")
        
        # Show top matched genes
        if matched_genes:
            sorted_matches = sorted(matched_genes.items(), key=lambda x: x[1], reverse=True)
            print(f"  Top matched genes:")
            for gene, count in sorted_matches[:5]:
                print(f"    {gene}: {count} k-mer matches")
        else:
            print(f"  NO MATCHES FOUND!")
        
        # Debug: show first few k-mers
        print(f"  First 3 k-mers from read:")
        for j, kmer in enumerate(list(read_kmers)[:3]):
            kmer_str = decode_kmer(kmer)
            in_index = "✓" if kmer in kmer_index else "✗"
            print(f"    {kmer_str} -> {kmer} {in_index}")
            if kmer in kmer_index:
                print(f"      Found in: {kmer_index[kmer][:2]}...")  # Show first 2 entries

if __name__ == "__main__":
    debug_kmer_matching()