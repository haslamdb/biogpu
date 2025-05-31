#!/usr/bin/env python3
"""Debug k-mer search with real metagenomic data"""

import gzip
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
        kmer_str = sequence[i:i+KMER_LENGTH].upper()
        encoded = encode_kmer(kmer_str)
        if encoded is not None:
            kmers.add(encoded)
    return kmers

def read_fastq_sample(filename, num_reads=100):
    """Read first few reads from a FASTQ file"""
    reads = []
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    
    count = 0
    while count < num_reads:
        header = f.readline().strip()
        if not header:
            break
        seq = f.readline().strip()
        plus = f.readline().strip()
        qual = f.readline().strip()
        
        reads.append({
            'header': header,
            'sequence': seq,
            'quality': qual
        })
        count += 1
    
    f.close()
    return reads

def load_kmer_index_headers(index_path):
    """Load just k-mer headers for faster testing"""
    print(f"Loading k-mer index headers from {index_path}")
    
    kmers_set = set()
    
    with open(index_path, 'rb') as f:
        # Read header
        num_kmers = struct.unpack('I', f.read(4))[0]
        kmer_size = struct.unpack('I', f.read(4))[0]
        
        print(f"Number of k-mers in file: {num_kmers}")
        print(f"K-mer size: {kmer_size}")
        
        # Read just the k-mer values, skip gene data
        for i in range(num_kmers):
            # Read encoded k-mer (64-bit)
            kmer_encoded = struct.unpack('Q', f.read(8))[0]
            kmers_set.add(kmer_encoded)
            
            # Read gene count
            gene_count = struct.unpack('I', f.read(4))[0]
            
            # Skip gene data
            for g in range(gene_count):
                organism_id = struct.unpack('I', f.read(4))[0]
                gene_name_len = struct.unpack('I', f.read(4))[0]
                f.seek(gene_name_len, 1)  # Skip gene name
                sequence_id = struct.unpack('I', f.read(4))[0]
            
            if i % 50000 == 0 and i > 0:
                print(f"  Loaded {i} k-mers...")
    
    print(f"Total unique k-mers loaded: {len(kmers_set)}")
    return kmers_set

def debug_real_data():
    # Load k-mer index (just the k-mer values for speed)
    index_path = 'data/gpu_resistance_db/kmers.bin'
    if not os.path.exists(index_path):
        print(f"Error: K-mer index not found at {index_path}")
        return
    
    print("Loading k-mer index...")
    kmer_index = load_kmer_index_headers(index_path)
    
    # Read sample of real data
    r1_path = 'data/test_fastq/VRE12_R1.fastq.gz'
    r2_path = 'data/test_fastq/VRE12_R2.fastq.gz'
    
    print(f"\nReading sample reads from {r1_path}")
    reads_r1 = read_fastq_sample(r1_path, 50)
    
    print(f"\n=== Testing k-mer matching for real metagenomic reads ===")
    total_kmers_tested = 0
    total_matches = 0
    
    for i, read in enumerate(reads_r1[:10]):
        seq = read['sequence']
        print(f"\nRead {i+1}:")
        print(f"  Length: {len(seq)}bp")
        print(f"  First 50bp: {seq[:50]}...")
        
        # Extract k-mers
        read_kmers = extract_kmers_encoded(seq)
        print(f"  Total valid k-mers: {len(read_kmers)}")
        
        # Count matches
        matching_kmers = sum(1 for kmer in read_kmers if kmer in kmer_index)
        print(f"  K-mers matching index: {matching_kmers} ({matching_kmers/len(read_kmers)*100:.1f}%)")
        
        total_kmers_tested += len(read_kmers)
        total_matches += matching_kmers
        
        # Show first few k-mers
        print(f"  First 5 k-mers:")
        for j, kmer in enumerate(list(read_kmers)[:5]):
            kmer_str = decode_kmer(kmer)
            in_index = "✓" if kmer in kmer_index else "✗"
            print(f"    {kmer_str} {in_index}")
    
    print(f"\n=== Overall Statistics ===")
    print(f"Total k-mers tested: {total_kmers_tested}")
    print(f"Total matches: {total_matches}")
    print(f"Overall match rate: {total_matches/total_kmers_tested*100:.1f}%")

if __name__ == "__main__":
    debug_real_data()