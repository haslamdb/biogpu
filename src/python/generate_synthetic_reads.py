#!/usr/bin/env python3
"""
Generate synthetic paired-end reads from a FASTA reference file.
"""

import random
import gzip
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))

def generate_quality_string(length, min_qual=20, max_qual=40):
    """Generate a random quality string for FASTQ format."""
    qualities = []
    for _ in range(length):
        # Generate quality scores with some variation
        if random.random() < 0.95:  # 95% high quality
            qual = random.randint(30, max_qual)
        else:  # 5% lower quality
            qual = random.randint(min_qual, 30)
        qualities.append(chr(qual + 33))  # Convert to Phred+33
    return ''.join(qualities)

def simulate_errors(sequence, error_rate=0.00001):
    """Introduce random sequencing errors."""
    bases = ['A', 'T', 'G', 'C']
    seq_list = list(sequence)
    
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            current = seq_list[i].upper()
            if current in bases:
                # Replace with a different base
                alternatives = [b for b in bases if b != current]
                seq_list[i] = random.choice(alternatives)
    
    return ''.join(seq_list)

def generate_paired_reads(reference_file, num_reads, read_length=150, insert_size_mean=350, insert_size_std=50):
    """Generate paired-end reads from reference sequences."""
    
    # Load all sequences from the reference
    print(f"Loading reference sequences from {reference_file}...")
    sequences = []
    total_length = 0
    
    with open(reference_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seq_str = str(record.seq).upper()
            sequences.append((record.id, seq_str))
            total_length += len(seq_str)
    
    print(f"Loaded {len(sequences)} sequences with total length {total_length:,} bp")
    
    # Open output files
    r1_file = gzip.open('synthetic_reads_20250529_R1.fastq.gz', 'wt')
    r2_file = gzip.open('synthetic_reads_20250529_R2.fastq.gz', 'wt')
    
    reads_generated = 0
    
    print(f"Generating {num_reads:,} paired-end reads...")
    
    while reads_generated < num_reads:
        # Select a random sequence
        seq_id, sequence = random.choice(sequences)
        seq_length = len(sequence)
        
        # Generate insert size
        insert_size = int(random.gauss(insert_size_mean, insert_size_std))
        insert_size = max(read_length * 2, insert_size)  # Ensure insert is at least 2x read length
        
        # Skip if sequence is too short
        if seq_length < insert_size:
            continue
        
        # Select random position for the insert
        start_pos = random.randint(0, seq_length - insert_size)
        
        # Extract the insert
        insert_seq = sequence[start_pos:start_pos + insert_size]
        
        # Generate forward read (R1) from the beginning of the insert
        r1_seq = insert_seq[:read_length]
        r1_seq = simulate_errors(r1_seq)
        
        # Generate reverse read (R2) from the end of the insert (reverse complement)
        r2_seq = insert_seq[-read_length:]
        r2_seq = reverse_complement(r2_seq)
        r2_seq = simulate_errors(r2_seq)
        
        # Generate quality strings
        r1_qual = generate_quality_string(read_length)
        r2_qual = generate_quality_string(read_length)
        
        # Create read IDs
        read_id = f"SYN_{reads_generated:07d}"
        r1_header = f"@{read_id}/1"
        r2_header = f"@{read_id}/2"
        
        # Write to FASTQ files
        r1_file.write(f"{r1_header}\n{r1_seq}\n+\n{r1_qual}\n")
        r2_file.write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")
        
        reads_generated += 1
        
        if reads_generated % 100000 == 0:
            print(f"  Generated {reads_generated:,} reads...")
    
    r1_file.close()
    r2_file.close()
    
    print(f"\nCompleted! Generated {reads_generated:,} paired-end reads")
    print(f"Output files:")
    print(f"  - synthetic_reads_R1.fastq.gz")
    print(f"  - synthetic_reads_R2.fastq.gz")

if __name__ == "__main__":
    # Parameters
    reference_file = "data/test_fastq/pneumo_oxytoca_Ecoli_aureus_serratia_burk.fasta"
    num_reads = 1000000  # 1M reads for testing
    read_length = 150
    
    # Generate the reads
    generate_paired_reads(reference_file, num_reads, read_length)