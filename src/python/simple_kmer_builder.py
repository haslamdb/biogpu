#!/usr/bin/env python3
"""
Simple k-mer index builder for FASTA files
"""

import argparse
import struct
import json
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def build_kmer_index(fasta_file, output_dir, kmer_length=15):
    """Build k-mer index from FASTA file"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load sequences
    sequences = []
    gene_map = {}
    
    logger.info(f"Loading sequences from {fasta_file}")
    for idx, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        sequences.append(str(record.seq))
        gene_map[idx] = record.id
        
    logger.info(f"Loaded {len(sequences)} sequences")
    
    # Build k-mer index
    kmer_index = defaultdict(list)
    
    for seq_idx, seq in enumerate(sequences):
        # Skip if sequence is too short
        if len(seq) < kmer_length:
            continue
            
        # Extract k-mers
        for pos in range(len(seq) - kmer_length + 1):
            kmer = seq[pos:pos + kmer_length]
            kmer_index[kmer].append((seq_idx, pos))
    
    logger.info(f"Built index with {len(kmer_index)} unique {kmer_length}-mers")
    
    # Write binary index
    index_file = output_dir / "kmer_index.bin"
    with open(index_file, 'wb') as f:
        # Header
        f.write(struct.pack('I', kmer_length))
        f.write(struct.pack('I', len(kmer_index)))
        
        # K-mers and positions
        for kmer, positions in sorted(kmer_index.items()):
            f.write(kmer.encode('ascii'))
            f.write(struct.pack('I', len(positions)))
            for seq_idx, pos in positions:
                f.write(struct.pack('II', seq_idx, pos))
    
    # Write sequences
    seq_file = output_dir / "sequences.bin"
    with open(seq_file, 'wb') as f:
        f.write(struct.pack('I', len(sequences)))
        for seq in sequences:
            f.write(struct.pack('I', len(seq)))
            f.write(seq.encode('ascii'))
    
    # Write metadata
    metadata = {
        'num_sequences': len(sequences),
        'kmer_length': kmer_length,
        'num_kmers': len(kmer_index),
        'gene_map': gene_map,
        'index_file': str(index_file),
        'sequences_file': str(seq_file)
    }
    
    with open(output_dir / 'metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"Index written to {output_dir}")
    return True

def main():
    parser = argparse.ArgumentParser(description='Build k-mer index from FASTA file')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('output_dir', help='Output directory for index')
    parser.add_argument('--kmer-length', type=int, default=15, help='K-mer length (default: 15)')
    
    args = parser.parse_args()
    
    success = build_kmer_index(args.fasta_file, args.output_dir, args.kmer_length)
    return 0 if success else 1

if __name__ == '__main__':
    exit(main())