#!/usr/bin/env python3
"""
Debug protein database gene assignments
"""
import os
import sys
import struct
import json

def read_protein_database(db_path):
    """Read and analyze protein database structure"""
    print(f"Analyzing protein database: {db_path}")
    
    # Read metadata if available
    metadata_path = os.path.join(db_path, 'metadata.json')
    if os.path.exists(metadata_path):
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        print(f"Database type: {metadata.get('database_type', 'unknown')}")
        print(f"Gene map: {metadata.get('gene_map', {})}")
        print(f"Species map: {metadata.get('species_map', {})}")
    
    # Read protein details if available
    details_path = os.path.join(db_path, 'protein_details.json')
    if os.path.exists(details_path):
        with open(details_path, 'r') as f:
            protein_details = json.load(f)
        
        print(f"\nProtein database contains {len(protein_details)} proteins:")
        print("ID | Gene     | Species                    | Length | Sequence Preview")
        print("-" * 80)
        
        for protein in protein_details[:20]:  # Show first 20
            print(f"{protein['id']:2d} | {protein['gene']:8s} | {protein['species']:25s} | {protein['length']:6d} | {protein['sequence']}")
        
        if len(protein_details) > 20:
            print(f"... and {len(protein_details) - 20} more proteins")
    
    # Read k-mer index
    kmer_path = os.path.join(db_path, 'protein_kmers.bin')
    if os.path.exists(kmer_path):
        with open(kmer_path, 'rb') as f:
            kmer_length = struct.unpack('I', f.read(4))[0]
            num_kmers = struct.unpack('I', f.read(4))[0]
            print(f"\nK-mer index: {num_kmers} unique {kmer_length}-mers")
    
    # Read binary proteins file
    proteins_path = os.path.join(db_path, 'proteins.bin')
    if os.path.exists(proteins_path):
        with open(proteins_path, 'rb') as f:
            num_proteins = struct.unpack('I', f.read(4))[0]
            print(f"Binary file: {num_proteins} proteins")
            
            # Read concatenated sequences
            remaining_data = f.read()
            print(f"Total sequence data: {len(remaining_data)} bytes")

def check_peptide_in_database(db_path, peptide, expected_gene):
    """Check where a specific peptide appears in the database"""
    details_path = os.path.join(db_path, 'protein_details.json')
    if not os.path.exists(details_path):
        print("No protein_details.json found")
        return
    
    with open(details_path, 'r') as f:
        protein_details = json.load(f)
    
    print(f"\nSearching for peptide: {peptide}")
    print(f"Expected gene: {expected_gene}")
    print("Matches found:")
    
    for protein in protein_details:
        full_seq_path = os.path.join(db_path, f"protein_{protein['id']}_full.txt")
        
        # For now, check the preview sequence
        if peptide in protein['sequence']:
            print(f"  Protein {protein['id']}: {protein['gene']} ({protein['species']}) - PREVIEW MATCH")
        
        # Also check if this matches expected gene
        if protein['gene'].lower() == expected_gene.lower():
            print(f"  Protein {protein['id']}: {protein['gene']} ({protein['species']}) - GENE MATCH")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python debug_protein_db.py <protein_db_path> [peptide] [expected_gene]")
        sys.exit(1)
    
    db_path = sys.argv[1]
    
    if not os.path.exists(db_path):
        print(f"Database path not found: {db_path}")
        sys.exit(1)
    
    read_protein_database(db_path)
    
    if len(sys.argv) >= 4:
        peptide = sys.argv[2]
        expected_gene = sys.argv[3]
        check_peptide_in_database(db_path, peptide, expected_gene)