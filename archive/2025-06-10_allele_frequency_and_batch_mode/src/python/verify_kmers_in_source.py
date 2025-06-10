#!/usr/bin/env python3
"""
Verify that extracted kmers exist in the original wildtype protein sequences
"""

import os
from pathlib import Path

def verify_kmers_in_fasta(wildtype_dir, kmers_to_check):
    """Check if kmers exist in the original FASTA files"""
    
    print(f"Checking kmers in wildtype sequences at: {wildtype_dir}")
    
    # Read all FASTA files
    protein_sequences = {}
    
    for fasta_file in Path(wildtype_dir).glob("*.fasta"):
        print(f"\nReading {fasta_file.name}...")
        with open(fasta_file, 'r') as f:
            current_header = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and current_seq:
                        seq = ''.join(current_seq)
                        protein_sequences[current_header] = seq
                        print(f"  Found sequence: {current_header[:50]}... (length: {len(seq)})")
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_header and current_seq:
                seq = ''.join(current_seq)
                protein_sequences[current_header] = seq
                print(f"  Found sequence: {current_header[:50]}... (length: {len(seq)})")
    
    print(f"\nTotal sequences loaded: {len(protein_sequences)}")
    
    # Now check each kmer
    print("\nVerifying k-mers:")
    print("-" * 80)
    
    for i, kmer in enumerate(kmers_to_check, 1):
        print(f"\n{i}. K-mer: '{kmer}'")
        found_in = []
        
        for header, sequence in protein_sequences.items():
            if kmer in sequence:
                # Find all positions
                positions = []
                start = 0
                while True:
                    pos = sequence.find(kmer, start)
                    if pos == -1:
                        break
                    positions.append(pos)
                    start = pos + 1
                
                found_in.append({
                    'header': header,
                    'positions': positions,
                    'context': []
                })
                
                # Get context around each position
                for pos in positions[:3]:  # Show first 3 occurrences
                    start = max(0, pos - 5)
                    end = min(len(sequence), pos + len(kmer) + 5)
                    context = sequence[start:end]
                    # Highlight the kmer
                    pre = sequence[start:pos]
                    post = sequence[pos + len(kmer):end]
                    found_in[-1]['context'].append(f"{pre}[{kmer}]{post}")
        
        if found_in:
            print(f"   ✓ FOUND in {len(found_in)} sequence(s):")
            for match in found_in[:3]:  # Show first 3 sequences
                print(f"     - {match['header'][:60]}...")
                print(f"       Positions: {match['positions'][:5]} {'...' if len(match['positions']) > 5 else ''}")
                for ctx in match['context'][:2]:
                    print(f"       Context: {ctx}")
        else:
            print(f"   ✗ NOT FOUND in any sequence")
    
    print("\n" + "=" * 80)
    
    # Also show some sample sequences for manual verification
    print("\nSample protein sequences (first 100 amino acids):")
    for i, (header, seq) in enumerate(list(protein_sequences.items())[:3]):
        print(f"\n{header}:")
        print(f"{seq[:100]}...")

if __name__ == "__main__":
    # K-mers extracted from the database
    kmers_to_check = [
        "LTSEM", "VARAL", "RILYG", "KKSAR", "GDYAI",
        "YRYML", "DGAAA", "TEMLR", "YDDTE", "LLVNG",
        "EADTK", "GEVLS", "YIEHQ", "EEGHY", "NDDFV",
        "DILAR", "VITLT", "FRAQR", "TSVTA", "NPDVT"
    ]
    
    wildtype_dir = "/home/david/Documents/Code/biogpu/data/wildtype_protein_seqs"
    verify_kmers_in_fasta(wildtype_dir, kmers_to_check)