#!/usr/bin/env python3
"""
Read the k-mer size from the protein_kmers.bin file.
The k-mer size is stored as a uint32_t in the first 4 bytes.
"""

import struct
import sys
import os

def read_kmer_size(file_path):
    """Read the k-mer size from a protein_kmers.bin file."""
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return None
    
    try:
        with open(file_path, 'rb') as f:
            # Read first 4 bytes (uint32_t)
            data = f.read(4)
            if len(data) < 4:
                print(f"Error: File too small, expected at least 4 bytes, got {len(data)}")
                return None
            
            # Unpack as little-endian unsigned 32-bit integer
            kmer_size = struct.unpack('<I', data)[0]
            return kmer_size
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def main():
    # Default path
    default_path = "data/integrated_clean_db_new/protein/protein_kmers.bin"
    
    # Use provided path or default
    file_path = sys.argv[1] if len(sys.argv) > 1 else default_path
    
    print(f"Reading k-mer size from: {file_path}")
    kmer_size = read_kmer_size(file_path)
    
    if kmer_size is not None:
        print(f"K-mer size: {kmer_size}")
        if kmer_size < 3 or kmer_size > 31:
            print(f"Warning: Unusual k-mer size detected. Common protein k-mer sizes are 3-11.")
    
    return 0 if kmer_size is not None else 1

if __name__ == "__main__":
    sys.exit(main())