#!/usr/bin/env python3
"""
Generate test data for BioGPU microbial profiler
Creates synthetic microbial sequences and FASTQ files for testing
"""

import random
import sys
from typing import List, Tuple, Dict
import gzip

# Test organisms with characteristic sequences
TEST_ORGANISMS = {
    "Escherichia_coli": {
        "taxon_id": 100,
        "marker": "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
                 "CGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGT",
        "gc_content": 0.51
    },
    "Klebsiella_pneumoniae": {
        "taxon_id": 101,
        "marker": "ATGAACATTAAAGGTCTGGTTGTTGCCGCTGCTGCCTTAGGTGGTGGCGCAACTGTCGCAGCAGAT"
                 "ATTGCTCAGGATAAACTGGAAGAGAAACTGAACGCTGCGTTGCAGGCCTTTGATAAAAACAAAGAT",
        "gc_content": 0.57
    },
    "Clostridioides_difficile": {
        "taxon_id": 108,
        "marker": "ATGAAAGTAAAAGAATTTAAAGCAGAAGAAGAAATTGATGCTGAAAAAGCTAAAGCAGCTGTTGAA"
                 "GCTATGGCTGCTGAAGATGCTAAAGCTGATATTGATGCTATTGAAGCTGATGAAGCTAAAGAAGAA",
        "gc_content": 0.29
    },
    "Staphylococcus_aureus": {
        "taxon_id": 103,
        "marker": "ATGAATATCAAAGAGCAAATTAAAGAATTAAAAGCAGAAGATGAAAAAGATAAAAAGAAAGATAAT"
                 "GAACAAGCTCAAAAGAAACTAGGTGTTCTAGATATTGAAGTGACAAATAATGTAAAAGATAGTAAT",
        "gc_content": 0.33
    },
    "Pseudomonas_aeruginosa": {
        "taxon_id": 107,
        "marker": "ATGAGCAAATCCTTGTCGCGAATCGCCTTGCTGGCCGTCGCCCTGGCGGCCAACCTGGTGAGCAAA"
                 "GATATTGCCCAGGGCAAACTGGAAGAGAAGCTGAACGCCGCCCTGCAGGCCTTCGACAAGAACAAG",
        "gc_content": 0.66
    }
}

def generate_sequence_with_gc(length: int, gc_content: float) -> str:
    """Generate a random sequence with specified GC content"""
    gc_count = int(length * gc_content)
    at_count = length - gc_count
    
    bases = ['G', 'C'] * (gc_count // 2) + ['A', 'T'] * (at_count // 2)
    
    # Handle odd counts
    if gc_count % 2:
        bases.append(random.choice(['G', 'C']))
    if at_count % 2:
        bases.append(random.choice(['A', 'T']))
    
    random.shuffle(bases)
    return ''.join(bases[:length])

def mutate_sequence(seq: str, error_rate: float = 0.01) -> str:
    """Add random mutations to simulate sequencing errors"""
    bases = ['A', 'C', 'G', 'T']
    mutated = list(seq)
    
    for i in range(len(mutated)):
        if random.random() < error_rate:
            current = mutated[i]
            mutated[i] = random.choice([b for b in bases if b != current])
    
    return ''.join(mutated)

def generate_organism_genome(name: str, info: Dict, genome_size: int = 50000) -> str:
    """Generate a synthetic genome for an organism"""
    marker = info["marker"]
    gc_content = info["gc_content"]
    
    # Start with random sequence
    genome = generate_sequence_with_gc(genome_size, gc_content)
    
    # Insert marker sequences at random positions
    num_copies = 5
    for _ in range(num_copies):
        pos = random.randint(0, genome_size - len(marker))
        genome = genome[:pos] + marker + genome[pos + len(marker):]
    
    # Add some conserved regions (simulate housekeeping genes)
    conserved_regions = [
        "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGG",  # Mock GFP-like
        "ATGAGCAAAGGAGAAGAACTTTTTACAGGAGTTGTACCAATTCTGGTTGAATTGGATGG",    # Mock RFP-like
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGG"   # Mock YFP-like
    ]
    
    for region in conserved_regions:
        # Add with slight mutations
        mutated_region = mutate_sequence(region, 0.05)
        pos = random.randint(0, genome_size - len(region))
        genome = genome[:pos] + mutated_region + genome[pos + len(region):]
    
    return genome

def extract_reads(genome: str, num_reads: int, read_length: int = 150) -> List[str]:
    """Extract random reads from a genome"""
    reads = []
    
    for _ in range(num_reads):
        if len(genome) <= read_length:
            reads.append(genome)
        else:
            start = random.randint(0, len(genome) - read_length)
            read = genome[start:start + read_length]
            
            # Add sequencing errors
            read = mutate_sequence(read, error_rate=0.005)
            
            # Randomly reverse complement some reads
            if random.random() < 0.5:
                read = reverse_complement(read)
            
            reads.append(read)
    
    return reads

def reverse_complement(seq: str) -> str:
    """Get reverse complement of a sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def create_fastq_file(filename: str, organism_abundances: Dict[str, float], 
                     total_reads: int = 100000, read_length: int = 150):
    """Create a FASTQ file with reads from multiple organisms"""
    
    # Generate genomes
    print("Generating synthetic genomes...")
    genomes = {}
    for org_name, org_info in TEST_ORGANISMS.items():
        genomes[org_name] = generate_organism_genome(org_name, org_info)
        print(f"  Generated {org_name} genome ({len(genomes[org_name])} bp)")
    
    # Calculate reads per organism
    reads_per_org = {}
    for org, abundance in organism_abundances.items():
        reads_per_org[org] = int(total_reads * abundance)
    
    # Handle rounding errors
    total_assigned = sum(reads_per_org.values())
    if total_assigned < total_reads:
        # Add remaining reads to most abundant organism
        most_abundant = max(organism_abundances, key=organism_abundances.get)
        reads_per_org[most_abundant] += total_reads - total_assigned
    
    # Generate reads and write FASTQ
    print(f"\nGenerating {total_reads} reads...")
    
    use_gzip = filename.endswith('.gz')
    if use_gzip:
        f = gzip.open(filename, 'wt')
    else:
        f = open(filename, 'w')
    
    read_id = 0
    for org_name, num_reads in reads_per_org.items():
        if num_reads == 0:
            continue
            
        print(f"  Generating {num_reads} reads from {org_name}")
        reads = extract_reads(genomes[org_name], num_reads, read_length)
        
        for read in reads:
            # FASTQ format
            f.write(f"@read_{read_id}_{org_name}\n")
            f.write(f"{read}\n")
            f.write("+\n")
            f.write("I" * len(read) + "\n")  # Fake quality scores
            read_id += 1
    
    f.close()
    print(f"Created FASTQ file: {filename}")

def create_kmer_list(filename: str, k: int = 31):
    """Create a list of k-mers and their taxon IDs for database building"""
    print(f"\nGenerating k-mer list (k={k})...")
    
    with open(filename, 'w') as f:
        f.write("# kmer\ttaxon_id\n")
        
        for org_name, org_info in TEST_ORGANISMS.items():
            genome = generate_organism_genome(org_name, org_info)
            taxon_id = org_info["taxon_id"]
            
            # Extract all k-mers
            kmers_written = 0
            for i in range(len(genome) - k + 1):
                kmer = genome[i:i+k]
                
                # Skip k-mers with N
                if 'N' in kmer:
                    continue
                
                # Get canonical k-mer
                rc = reverse_complement(kmer)
                canonical = min(kmer, rc)
                
                # Write every 10th k-mer (to keep file size reasonable)
                if i % 10 == 0:
                    f.write(f"{canonical}\t{taxon_id}\n")
                    kmers_written += 1
            
            print(f"  {org_name}: {kmers_written} k-mers")

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_test_data.py <command> [options]")
        print("\nCommands:")
        print("  fastq    - Generate test FASTQ file")
        print("  kmers    - Generate k-mer list for database building")
        print("  both     - Generate both FASTQ and k-mer list")
        print("\nExamples:")
        print("  python generate_test_data.py fastq test.fastq")
        print("  python generate_test_data.py kmers kmers.txt")
        print("  python generate_test_data.py both test")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "fastq":
        if len(sys.argv) < 3:
            print("Usage: python generate_test_data.py fastq <output.fastq>")
            sys.exit(1)
        
        # Define community composition (relative abundances)
        abundances = {
            "Escherichia_coli": 0.3,
            "Klebsiella_pneumoniae": 0.2,
            "Clostridioides_difficile": 0.25,  # Your research focus!
            "Staphylococcus_aureus": 0.15,
            "Pseudomonas_aeruginosa": 0.1
        }
        
        create_fastq_file(sys.argv[2], abundances, total_reads=50000)
    
    elif command == "kmers":
        if len(sys.argv) < 3:
            print("Usage: python generate_test_data.py kmers <output.txt>")
            sys.exit(1)
        
        create_kmer_list(sys.argv[2])
    
    elif command == "both":
        if len(sys.argv] < 3:
            print("Usage: python generate_test_data.py both <prefix>")
            sys.exit(1)
        
        prefix = sys.argv[2]
        
        # Generate k-mer list
        create_kmer_list(f"{prefix}_kmers.txt")
        
        # Generate FASTQ with mixed community
        abundances = {
            "Escherichia_coli": 0.3,
            "Klebsiella_pneumoniae": 0.2,
            "Clostridioides_difficile": 0.25,
            "Staphylococcus_aureus": 0.15,
            "Pseudomonas_aeruginosa": 0.1
        }
        
        create_fastq_file(f"{prefix}.fastq", abundances)
        
        print(f"\nTest data generation complete!")
        print(f"K-mer list: {prefix}_kmers.txt")
        print(f"FASTQ file: {prefix}.fastq")
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

if __name__ == "__main__":
    main()