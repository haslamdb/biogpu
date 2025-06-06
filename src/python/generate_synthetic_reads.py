#!/usr/bin/env python3
"""
Generate clean synthetic paired-end reads from multi-genome FASTA files.
Optimized for generating ~10M read pairs without sequence errors.
"""

import random
import gzip
import sys
import os
import argparse
from Bio import SeqIO
from pathlib import Path

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))

def generate_quality_string(length, qual=40):
    """Generate high-quality string for FASTQ format (all Q40)."""
    return chr(qual + 33) * length

def load_genome_sequences(fasta_file):
    """Load all sequences from a multi-genome FASTA file."""
    sequences = []
    total_length = 0
    
    print(f"Loading sequences from {fasta_file}...")
    
    try:
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                seq_str = str(record.seq).upper().replace('N', '')  # Remove N's
                # Filter out sequences shorter than minimum insert size
                if len(seq_str) >= 500:  # Minimum for 150bp paired reads
                    sequences.append({
                        'id': record.id,
                        'sequence': seq_str,
                        'length': len(seq_str)
                    })
                    total_length += len(seq_str)
    
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
        return [], 0
    
    return sequences, total_length

def generate_paired_reads_from_genomes(genome_files, num_reads, output_prefix, 
                                     read_length=150, insert_size_mean=350, 
                                     insert_size_std=50):
    """Generate clean paired-end reads from multiple genome files."""
    
    # Load all sequences from all genome files
    all_sequences = []
    total_genome_length = 0
    
    print("Loading genome sequences...")
    for genome_file in genome_files:
        sequences, file_length = load_genome_sequences(genome_file)
        all_sequences.extend(sequences)
        total_genome_length += file_length
        print(f"  {genome_file}: {len(sequences)} sequences, {file_length:,} bp")
    
    if not all_sequences:
        print("ERROR: No sequences loaded!")
        return False
    
    print(f"\nTotal: {len(all_sequences)} sequences from {len(genome_files)} files")
    print(f"Combined genome length: {total_genome_length:,} bp")
    
    # Create cumulative length array for weighted random selection
    cumulative_lengths = []
    running_total = 0
    for seq in all_sequences:
        running_total += seq['length']
        cumulative_lengths.append(running_total)
    
    # Open output files
    r1_filename = f"{output_prefix}_R1.fastq.gz"
    r2_filename = f"{output_prefix}_R2.fastq.gz"
    
    print(f"\nGenerating {num_reads:,} paired-end reads...")
    print(f"Output files: {r1_filename}, {r2_filename}")
    
    with gzip.open(r1_filename, 'wt') as r1_file, gzip.open(r2_filename, 'wt') as r2_file:
        reads_generated = 0
        failed_attempts = 0
        
        while reads_generated < num_reads:
            # Select sequence weighted by length
            rand_pos = random.randint(0, total_genome_length - 1)
            seq_idx = 0
            for i, cum_len in enumerate(cumulative_lengths):
                if rand_pos < cum_len:
                    seq_idx = i
                    break
            
            selected_seq = all_sequences[seq_idx]
            sequence = selected_seq['sequence']
            seq_length = selected_seq['length']
            
            # Generate insert size
            insert_size = int(random.gauss(insert_size_mean, insert_size_std))
            insert_size = max(read_length * 2 + 50, insert_size)  # Minimum reasonable insert
            insert_size = min(1000, insert_size)  # Maximum reasonable insert
            
            # Skip if sequence is too short for this insert
            if seq_length < insert_size:
                failed_attempts += 1
                if failed_attempts > num_reads * 10:  # Prevent infinite loop
                    print("ERROR: Too many failed attempts. Check sequence lengths.")
                    break
                continue
            
            # Select random position for the insert
            start_pos = random.randint(0, seq_length - insert_size)
            
            # Extract the insert
            insert_seq = sequence[start_pos:start_pos + insert_size]
            
            # Generate R1: forward read from start of insert
            r1_seq = insert_seq[:read_length]
            
            # Generate R2: reverse complement from end of insert
            r2_seq = reverse_complement(insert_seq[-read_length:])
            
            # Skip reads with N's or non-ATGC characters
            if 'N' in r1_seq or 'N' in r2_seq:
                continue
            if not all(c in 'ATGC' for c in r1_seq + r2_seq):
                continue
            
            # Generate high-quality scores
            r1_qual = generate_quality_string(read_length)
            r2_qual = generate_quality_string(read_length)
            
            # Create read IDs
            read_id = f"SYN_{reads_generated:08d}"
            r1_header = f"@{read_id}/1"
            r2_header = f"@{read_id}/2"
            
            # Write to FASTQ files
            r1_file.write(f"{r1_header}\n{r1_seq}\n+\n{r1_qual}\n")
            r2_file.write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")
            
            reads_generated += 1
            
            if reads_generated % 500000 == 0:
                print(f"  Generated {reads_generated:,} read pairs...")
    
    print(f"\nCompleted! Generated {reads_generated:,} paired-end reads")
    print(f"Output files:")
    print(f"  - {r1_filename}")
    print(f"  - {r2_filename}")
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic paired-end reads from genome FASTA files")
    parser.add_argument("genome_files", nargs="+", help="One or more FASTA files containing genome sequences")
    parser.add_argument("-n", "--num_reads", type=int, default=10000000, 
                       help="Number of read pairs to generate (default: 10M)")
    parser.add_argument("-o", "--output_prefix", default="synthetic_reads", 
                       help="Output file prefix (default: synthetic_reads)")
    parser.add_argument("-l", "--read_length", type=int, default=150, 
                       help="Length of each read (default: 150)")
    parser.add_argument("--insert_mean", type=int, default=350, 
                       help="Mean insert size (default: 350)")
    parser.add_argument("--insert_std", type=int, default=50, 
                       help="Insert size standard deviation (default: 50)")
    
    args = parser.parse_args()
    
    # Check input files exist
    for genome_file in args.genome_files:
        if not os.path.exists(genome_file):
            print(f"ERROR: File not found: {genome_file}")
            sys.exit(1)
    
    # Generate reads
    success = generate_paired_reads_from_genomes(
        args.genome_files,
        args.num_reads,
        args.output_prefix,
        args.read_length,
        args.insert_mean,
        args.insert_std
    )
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    # Can be run directly with default parameters or use command line arguments
    if len(sys.argv) == 1:
        # Default execution for testing
        genome_files = ["data/test_fastq/pneumo_oxytoca_Ecoli_aureus_serratia_burk.fasta"]
        generate_paired_reads_from_genomes(
            genome_files, 
            1000000,  # 1M reads for testing
            "synthetic_reads_clean",
            read_length=150,
            insert_size_mean=350,
            insert_size_std=50
        )
    else:
        main()