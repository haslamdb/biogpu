#!/usr/bin/env python3
"""
Extract and display sample kmers from the protein database binary file
Built with kmer-length 5
"""

import struct
import json
import os
from pathlib import Path

def read_protein_kmers(db_path, kmer_length=5, num_samples=20):
    """Read and display sample kmers from the protein database"""
    
    proteins_bin = os.path.join(db_path, "protein", "proteins.bin")
    metadata_json = os.path.join(db_path, "protein", "metadata.json")
    
    print(f"Reading protein kmers from: {proteins_bin}")
    print(f"Expected k-mer length: {kmer_length}")
    
    # Check metadata
    if os.path.exists(metadata_json):
        with open(metadata_json, 'r') as f:
            metadata = json.load(f)
            print(f"\nMetadata:")
            for key, value in metadata.items():
                print(f"  {key}: {value}")
    
    # Read binary file
    if not os.path.exists(proteins_bin):
        print(f"Error: {proteins_bin} not found!")
        return
    
    file_size = os.path.getsize(proteins_bin)
    print(f"\nBinary file size: {file_size} bytes")
    
    # Let's also check the source wildtype proteins
    wildtype_dir = "/home/david/Documents/Code/biogpu/data/wildtype_protein_seqs"
    if os.path.exists(wildtype_dir):
        print(f"\nSource wildtype protein files:")
        for f in sorted(os.listdir(wildtype_dir))[:5]:
            print(f"  {f}")
    
    with open(proteins_bin, 'rb') as f:
        # Read assuming structure: header (kmer_size, num_kmers) then kmer entries
        
        # First, let's see the raw data
        sample_data = f.read(min(200, file_size))
        f.seek(0)
        
        print(f"\nFirst 100 bytes (hex):")
        for i in range(0, min(100, len(sample_data)), 16):
            hex_str = ' '.join([f'{b:02x}' for b in sample_data[i:i+16]])
            ascii_str = ''.join([chr(b) if 32 <= b <= 126 else '.' for b in sample_data[i:i+16]])
            print(f"{i:04x}: {hex_str:<48} {ascii_str}")
        
        # Try to read with kmer_length = 5
        f.seek(0)
        
        # Skip any header (try different offsets)
        offsets_to_try = [0, 4, 8, 12]
        
        for offset in offsets_to_try:
            print(f"\n\nTrying offset {offset} with k-mer length {kmer_length}:")
            f.seek(offset)
            
            try:
                # Try to read entries assuming structure: kmer(5 bytes) + gene_id(4) + species_id(4) + position(4) = 17 bytes per entry
                entry_size = kmer_length + 12  # 5 + 4 + 4 + 4 = 17 bytes
                
                print(f"Entry size: {entry_size} bytes")
                
                extracted_kmers = []
                for i in range(num_samples):
                    start_pos = f.tell()
                    entry_data = f.read(entry_size)
                    
                    if len(entry_data) < entry_size:
                        break
                    
                    # Extract kmer
                    kmer_bytes = entry_data[:kmer_length]
                    kmer_str = ''.join([chr(b) if 32 <= b <= 126 else f'\\x{b:02x}' for b in kmer_bytes])
                    
                    # Extract metadata
                    gene_id = struct.unpack('I', entry_data[kmer_length:kmer_length+4])[0]
                    species_id = struct.unpack('I', entry_data[kmer_length+4:kmer_length+8])[0]
                    position = struct.unpack('I', entry_data[kmer_length+8:kmer_length+12])[0]
                    
                    # Check if it looks valid (amino acids)
                    amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
                    clean_kmer = ''.join([c for c in kmer_str if c in amino_acids])
                    
                    if len(clean_kmer) == kmer_length:
                        extracted_kmers.append({
                            'kmer': clean_kmer,
                            'gene_id': gene_id,
                            'species_id': species_id,
                            'position': position,
                            'raw_hex': ' '.join([f'{b:02x}' for b in kmer_bytes])
                        })
                
                if extracted_kmers:
                    print(f"\nSuccessfully extracted {len(extracted_kmers)} k-mers:")
                    for i, kmer_info in enumerate(extracted_kmers[:10]):
                        print(f"{i+1:3}. K-mer: '{kmer_info['kmer']}' | Gene: {kmer_info['gene_id']:10} | Species: {kmer_info['species_id']:10} | Pos: {kmer_info['position']:10}")
                        print(f"     Hex: {kmer_info['raw_hex']}")
                    
                    # Return the extracted kmers for verification
                    return extracted_kmers
                    
            except Exception as e:
                print(f"Failed with offset {offset}: {e}")
                continue

    return None

if __name__ == "__main__":
    db_path = "/home/david/Documents/Code/biogpu/data/integrated_clean_db"
    extracted = read_protein_kmers(db_path, kmer_length=5, num_samples=50)
    
    if extracted:
        print(f"\n\nExtracted k-mers for verification:")
        print("These are the k-mers to check against the original protein sequences:")
        unique_kmers = list(set([k['kmer'] for k in extracted]))[:20]
        for i, kmer in enumerate(unique_kmers):
            print(f"  {i+1}. {kmer}")