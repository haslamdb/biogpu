#!/usr/bin/env python3
"""
Test loading and using the fluoroquinolone mutation index
"""

import h5py
import json
import numpy as np
from pathlib import Path

def test_load_index(index_dir: str = "data/indices/fq_mutations"):
    """Test loading the FQ mutation index"""
    index_path = Path(index_dir)
    
    print("=== Testing FQ Mutation Index Loading ===\n")
    
    # Load summary
    with open(index_path / "index_summary.json", 'r') as f:
        summary = json.load(f)
    
    print(f"Index Summary:")
    print(f"  Total sequences: {summary['total_sequences']}")
    print(f"  Genes indexed: {', '.join(summary['genes'])}")
    print(f"  Total k-mers: {summary['total_kmers']}")
    print(f"  K-mer size: {summary['kmer_size']}")
    print()
    
    # Load HDF5 index
    with h5py.File(index_path / "fq_mutation_index.h5", 'r') as f:
        print("HDF5 Structure:")
        
        for gene in summary['genes']:
            if gene in f:
                gene_group = f[gene]
                print(f"\n  Gene: {gene}")
                print(f"    Sequences shape: {gene_group['sequences'].shape}")
                print(f"    Mutation masks shape: {gene_group['mutation_masks'].shape}")
                print(f"    Max length: {gene_group.attrs['max_length']}")
                
                # Load metadata
                metadata = json.loads(gene_group.attrs['metadata'])
                print(f"    Number of sequences: {len(metadata['seq_ids'])}")
                print(f"    Sample species: {metadata['species'][:3]}")
                
                # Check mutation info
                total_mutations = sum(len(m) for m in metadata['mutations'])
                print(f"    Total mutations tracked: {total_mutations}")
        
        # Check k-mer index
        if 'kmer_index' in f:
            kmer_group = f['kmer_index']
            print(f"\n  K-mer Index:")
            print(f"    K-mers shape: {kmer_group['kmers'].shape}")
            print(f"    K-mer size: {kmer_group.attrs['k']}")
    
    # Load k-mer lookup
    with open(index_path / "kmer_lookup.json", 'r') as f:
        kmer_data = json.load(f)
    
    print(f"\n  K-mer Lookup:")
    print(f"    Total unique k-mers: {len(kmer_data['kmer_to_idx'])}")
    print(f"    Sample k-mers: {list(kmer_data['kmer_to_idx'].keys())[:3]}")
    
    # Memory usage estimate
    print("\n=== Memory Usage Estimate ===")
    
    # Calculate sizes
    with h5py.File(index_path / "fq_mutation_index.h5", 'r') as f:
        total_elements = 0
        for gene in summary['genes']:
            if gene in f:
                seq_array = f[gene]['sequences']
                total_elements += seq_array.shape[0] * seq_array.shape[1]
        
        # Each element is uint8 (1 byte)
        gpu_memory_mb = total_elements / (1024 * 1024)
        
    print(f"  Sequence data for GPU: ~{gpu_memory_mb:.2f} MB")
    print(f"  K-mer lookup size: ~{len(json.dumps(kmer_data)) / (1024 * 1024):.2f} MB")
    print(f"  Total estimated GPU memory: ~{gpu_memory_mb + 10:.2f} MB (with overhead)")
    
    print("\n✓ Index loading test completed successfully!")
    
    return True

def test_mutation_lookup(index_dir: str = "data/indices/fq_mutations"):
    """Test looking up specific mutations"""
    index_path = Path(index_dir)
    
    print("\n=== Testing Mutation Lookup ===\n")
    
    # Load a sample gene group and check mutations
    with h5py.File(index_path / "fq_mutation_index.h5", 'r') as f:
        # Test with gyrA (most common)
        if 'gyrA' in f:
            gene_group = f['gyrA']
            metadata = json.loads(gene_group.attrs['metadata'])
            
            # Find sequences with mutations at common positions (83, 87)
            for i, (seq_id, mutations) in enumerate(zip(metadata['seq_ids'], metadata['mutations'])):
                if i >= 3:  # Just show first 3
                    break
                    
                if mutations:
                    print(f"Sequence: {seq_id}")
                    print(f"  Species: {metadata['species'][i]}")
                    for mut in mutations:
                        print(f"  Mutation: {mut['wild_type']}{mut['codon_position']}{mut['mutant']} at position {mut['position']}")
                    print()
    
    return True

if __name__ == "__main__":
    import sys
    
    index_dir = sys.argv[1] if len(sys.argv) > 1 else "data/indices/fq_mutations"
    
    # Run tests
    test_load_index(index_dir)
    test_mutation_lookup(index_dir)