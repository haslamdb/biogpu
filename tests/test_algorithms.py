#!/usr/bin/env python3
"""
BioGPU Testing Framework
Tests core algorithms with CPU fallback when GPU is not available
"""

import numpy as np
import time
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import subprocess
import os

# Try to import GPU libraries
try:
    import cupy as cp
    GPU_AVAILABLE = cp.cuda.runtime.getDeviceCount() > 0
except:
    GPU_AVAILABLE = False
    cp = np  # Fallback to numpy

print(f"GPU Available: {GPU_AVAILABLE}")

@dataclass
class DNASequence:
    """Represents a DNA sequence with quality scores"""
    sequence: str
    quality: Optional[str] = None
    read_id: str = ""
    
    def encode(self) -> np.ndarray:
        """Convert DNA to numeric encoding"""
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        return np.array([encoding.get(c, 4) for c in self.sequence.upper()], dtype=np.uint8)

@dataclass
class ResistanceMutation:
    """Known resistance mutation"""
    gene: str
    position: int
    wild_type: str
    mutant: str
    resistance_level: float  # MIC fold change

class KmerCounter:
    """K-mer counting implementation"""
    
    def __init__(self, k: int = 31):
        self.k = k
        
    def count_kmers_cpu(self, sequences: List[DNASequence]) -> Dict[int, int]:
        """CPU implementation of k-mer counting"""
        kmer_counts = {}
        
        for seq in sequences:
            encoded = seq.encode()
            for i in range(len(encoded) - self.k + 1):
                kmer = encoded[i:i+self.k]
                # Simple hash function
                hash_val = 0
                for j, base in enumerate(kmer):
                    hash_val = (hash_val * 5 + base) % (2**32)
                
                kmer_counts[hash_val] = kmer_counts.get(hash_val, 0) + 1
                
        return kmer_counts
    
    def count_kmers_gpu(self, sequences: List[DNASequence]) -> Dict[int, int]:
        """GPU implementation of k-mer counting"""
        if not GPU_AVAILABLE:
            print("GPU not available, falling back to CPU")
            return self.count_kmers_cpu(sequences)
        
        # Convert sequences to GPU array
        max_len = max(len(seq.sequence) for seq in sequences)
        num_seqs = len(sequences)
        
        # Pad sequences and encode
        encoded_seqs = cp.zeros((num_seqs, max_len), dtype=cp.uint8)
        seq_lengths = cp.zeros(num_seqs, dtype=cp.int32)
        
        for i, seq in enumerate(sequences):
            encoded = seq.encode()
            encoded_seqs[i, :len(encoded)] = cp.asarray(encoded)
            seq_lengths[i] = len(encoded)
        
        # GPU kernel would be called here
        # For now, simulate with CuPy operations
        kmer_counts = {}
        
        # This is a simplified version - real implementation would use custom kernel
        for i in range(num_seqs):
            seq_len = int(seq_lengths[i])
            for j in range(seq_len - self.k + 1):
                kmer = encoded_seqs[i, j:j+self.k]
                hash_val = int(cp.sum(kmer * cp.arange(self.k)) % (2**32))
                kmer_counts[hash_val] = kmer_counts.get(hash_val, 0) + 1
                
        return kmer_counts

class ShortReadAligner:
    """Bowtie2-style short read aligner"""
    
    def __init__(self, reference_kmers: Dict[int, List[int]]):
        self.reference_kmers = reference_kmers
        self.seed_length = 22
        
    def align_cpu(self, reads: List[DNASequence]) -> List[Tuple[int, int, int]]:
        """CPU implementation of short read alignment"""
        alignments = []
        
        for read_idx, read in enumerate(reads):
            encoded = read.encode()
            
            # Extract seeds
            for i in range(0, len(encoded) - self.seed_length + 1, 11):
                seed = encoded[i:i+self.seed_length]
                seed_hash = hash(seed.tobytes()) % (2**32)
                
                # Check if seed exists in reference
                if seed_hash in self.reference_kmers:
                    # Simple alignment: just report first hit
                    ref_positions = self.reference_kmers[seed_hash]
                    if ref_positions:
                        alignments.append((read_idx, ref_positions[0], i))
                        break
                        
        return alignments

class ResistanceScanner:
    """Scan for antibiotic resistance mutations"""
    
    def __init__(self, known_mutations: List[ResistanceMutation]):
        self.mutations = known_mutations
        self.codon_table = self._build_codon_table()
        
    def _build_codon_table(self) -> Dict[str, str]:
        """Build DNA codon to amino acid table"""
        # Simplified codon table
        return {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
    def scan_mutations_cpu(self, gene_sequences: Dict[str, str]) -> List[Dict]:
        """Scan for known resistance mutations"""
        found_mutations = []
        
        for mutation in self.mutations:
            if mutation.gene in gene_sequences:
                gene_seq = gene_sequences[mutation.gene]
                
                # Extract codon at mutation position (1-based to 0-based)
                codon_start = (mutation.position - 1) * 3
                if codon_start + 3 <= len(gene_seq):
                    codon = gene_seq[codon_start:codon_start + 3]
                    
                    # Translate codon
                    if codon in self.codon_table:
                        aa = self.codon_table[codon]
                        
                        # Check if it's the mutant amino acid
                        if aa == mutation.mutant:
                            found_mutations.append({
                                'gene': mutation.gene,
                                'position': mutation.position,
                                'wild_type': mutation.wild_type,
                                'mutant': aa,
                                'codon': codon,
                                'resistance_level': mutation.resistance_level
                            })
                            
        return found_mutations

class MicrobialAbundanceCalculator:
    """Calculate microbial abundance from alignments"""
    
    def calculate_tpm(self, read_counts: Dict[str, int], 
                      gene_lengths: Dict[str, int]) -> Dict[str, float]:
        """Calculate Transcripts Per Million (TPM) normalization"""
        # Calculate RPK (reads per kilobase)
        rpk = {}
        for gene, count in read_counts.items():
            if gene in gene_lengths:
                rpk[gene] = count / (gene_lengths[gene] / 1000.0)
                
        # Calculate scaling factor
        scaling_factor = sum(rpk.values()) / 1e6
        
        # Calculate TPM
        tpm = {}
        for gene, rpk_value in rpk.items():
            tpm[gene] = rpk_value / scaling_factor if scaling_factor > 0 else 0
            
        return tpm

def test_kmer_counting():
    """Test k-mer counting functionality"""
    print("\n=== Testing K-mer Counting ===")
    
    # Generate test sequences
    test_sequences = [
        DNASequence("ATCGATCGATCGATCG", read_id="read1"),
        DNASequence("GCTAGCTAGCTAGCTA", read_id="read2"),
        DNASequence("ATCGATCGATCGATCG", read_id="read3"),  # Duplicate
    ]
    
    counter = KmerCounter(k=8)
    
    # Test CPU implementation
    start = time.time()
    cpu_counts = counter.count_kmers_cpu(test_sequences)
    cpu_time = time.time() - start
    
    print(f"CPU k-mer counting: {len(cpu_counts)} unique k-mers in {cpu_time:.4f}s")
    
    # Test GPU implementation
    start = time.time()
    gpu_counts = counter.count_kmers_gpu(test_sequences)
    gpu_time = time.time() - start
    
    print(f"GPU k-mer counting: {len(gpu_counts)} unique k-mers in {gpu_time:.4f}s")

def test_resistance_detection():
    """Test fluoroquinolone resistance detection"""
    print("\n=== Testing Resistance Detection ===")
    
    # Define known fluoroquinolone resistance mutations
    fq_mutations = [
        ResistanceMutation("gyrA", 83, "S", "L", 8.0),  # S83L
        ResistanceMutation("gyrA", 87, "D", "N", 4.0),  # D87N
        ResistanceMutation("parC", 80, "S", "I", 4.0),  # S80I
    ]
    
    # Test gene sequences with mutations
    test_genes = {
        "gyrA": "..." + "TTG" + "..." + "AAC" + "...",  # Has S83L (TTG=L) and D87N (AAC=N)
        "parC": "..." + "ATC" + "...",  # Has S80I (ATC=I)
    }
    
    # Properly construct sequences with correct codon positions
    # Position 83 means codon starts at (83-1)*3 = 246
    gyrA_seq = "A" * 246 + "TTG" + "A" * 9 + "AAC" + "A" * 100  # S83L and D87N
    parC_seq = "A" * 237 + "ATC" + "A" * 100  # S80I
    
    test_genes = {"gyrA": gyrA_seq, "parC": parC_seq}
    
    scanner = ResistanceScanner(fq_mutations)
    mutations_found = scanner.scan_mutations_cpu(test_genes)
    
    print(f"Found {len(mutations_found)} resistance mutations:")
    for mut in mutations_found:
        print(f"  - {mut['gene']} {mut['wild_type']}{mut['position']}{mut['mutant']} "
              f"(codon: {mut['codon']}, resistance: {mut['resistance_level']}x)")

def test_llvm_compilation():
    """Test LLVM IR generation"""
    print("\n=== Testing LLVM Compilation ===")
    
    # Compile and run the main executable
    result = subprocess.run(["./biogpu"], capture_output=True, text=True)
    
    if result.returncode == 0:
        print("✓ LLVM module generated successfully")
        print("Generated functions:")
        # Parse output to show generated functions
        for line in result.stdout.split('\n'):
            if 'define' in line and 'kernel' in line:
                func_name = line.split()[2].split('(')[0]
                print(f"  - {func_name}")
    else:
        print("✗ LLVM compilation failed")
        print(result.stderr)

def benchmark_algorithms():
    """Benchmark different algorithms"""
    print("\n=== Algorithm Benchmarks ===")
    
    # Generate larger test dataset
    num_reads = 10000
    read_length = 150
    
    print(f"Generating {num_reads} reads of length {read_length}...")
    
    # Random DNA sequences
    bases = ['A', 'C', 'G', 'T']
    sequences = []
    for i in range(num_reads):
        seq = ''.join(np.random.choice(bases, read_length))
        sequences.append(DNASequence(seq, read_id=f"read_{i}"))
    
    # Benchmark k-mer counting
    counter = KmerCounter(k=31)
    
    start = time.time()
    cpu_counts = counter.count_kmers_cpu(sequences[:1000])  # Subset for CPU
    cpu_time = time.time() - start
    
    print(f"\nK-mer counting (k=31):")
    print(f"  CPU: {len(cpu_counts)} k-mers from 1000 reads in {cpu_time:.2f}s")
    print(f"  Estimated for {num_reads} reads: {cpu_time * num_reads / 1000:.2f}s")
    
    if GPU_AVAILABLE:
        start = time.time()
        gpu_counts = counter.count_kmers_gpu(sequences)
        gpu_time = time.time() - start
        print(f"  GPU: {len(gpu_counts)} k-mers from {num_reads} reads in {gpu_time:.2f}s")
        print(f"  Speedup: {(cpu_time * num_reads / 1000) / gpu_time:.1f}x")

if __name__ == "__main__":
    print("BioGPU Testing Framework")
    print("=" * 50)
    
    # Run all tests
    test_kmer_counting()
    test_resistance_detection()
    
    # Only run if binary exists
    if os.path.exists("./biogpu"):
        test_llvm_compilation()
    else:
        print("\n[INFO] Skipping LLVM test - compile with 'make' first")
    
    benchmark_algorithms()
    
    print("\n" + "=" * 50)
    print("Testing complete!")
    
    if not GPU_AVAILABLE:
        print("\n⚠️  GPU not available - all tests ran on CPU")
        print("   Install GPU to see acceleration benefits")