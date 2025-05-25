#!/usr/bin/env python3
"""
BioGPU Database Builder
Prepares reference genomes and mutation databases for GPU processing
"""

import json
import struct
import numpy as np
from Bio import SeqIO
from typing import Dict, List, Tuple
import argparse
import os
import gzip
import urllib.request
from dataclasses import dataclass, asdict

@dataclass
class GeneInfo:
    """Information about a resistance gene"""
    gene_id: int
    gene_name: str
    organism: str
    sequence: str
    length: int
    qrdr_start: int = 0
    qrdr_end: int = 0

@dataclass
class MutationInfo:
    """Fluoroquinolone resistance mutation"""
    gene_name: str
    organism: str
    position: int
    wild_type: str
    mutant: str
    mic_change: float
    drugs: List[str]
    pubmed_ids: List[int]

class BiogpuDatabaseBuilder:
    def __init__(self, output_dir: str = "data/biogpu_db"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.genes = {}
        self.mutations = []
        self.organisms = {}
        
    def download_reference_genes(self):
        """Download key resistance genes from NCBI"""
        print("Downloading reference genes...")
        
        # Example genes for E. coli
        reference_genes = {
            # Gene: (GenBank ID, QRDR region)
            "gyrA": ("NC_000913.3:2336102-2338714", (67, 106)),
            "gyrB": ("NC_000913.3:3735943-3738348", (426, 447)),
            "parC": ("NC_000913.3:3168901-3171203", (64, 102)),
            "parE": ("NC_000913.3:3169898-3171898", (416, 420))
        }
        
        # In real implementation, fetch from NCBI
        # For now, create example sequences
        for gene_name, (accession, qrdr) in reference_genes.items():
            # Simulated gene sequence
            gene_seq = "ATG" + "GCT" * 300 + "TAA"  # Simplified
            
            gene_info = GeneInfo(
                gene_id=len(self.genes),
                gene_name=gene_name,
                organism="Escherichia coli",
                sequence=gene_seq,
                length=len(gene_seq),
                qrdr_start=qrdr[0],
                qrdr_end=qrdr[1]
            )
            
            self.genes[gene_name] = gene_info
            
        print(f"Loaded {len(self.genes)} reference genes")
        
    def load_mutation_database(self, json_file: str = None):
        """Load fluoroquinolone resistance mutations"""
        if json_file and os.path.exists(json_file):
            with open(json_file) as f:
                data = json.load(f)
                mutations = data.get("mutations", [])
        else:
            # Default mutations for E. coli
            mutations = [
                {
                    "gene": "gyrA",
                    "organism": "Escherichia coli",
                    "position": 83,
                    "wild_type": "S",
                    "mutations": ["L", "W"],
                    "mic_change": 8.0,
                    "drugs": ["ciprofloxacin", "levofloxacin"]
                },
                {
                    "gene": "gyrA",
                    "organism": "Escherichia coli",
                    "position": 87,
                    "wild_type": "D",
                    "mutations": ["N", "G", "Y"],
                    "mic_change": 4.0,
                    "drugs": ["ciprofloxacin", "levofloxacin"]
                },
                {
                    "gene": "parC",
                    "organism": "Escherichia coli",
                    "position": 80,
                    "wild_type": "S",
                    "mutations": ["I", "R"],
                    "mic_change": 4.0,
                    "drugs": ["ciprofloxacin", "levofloxacin"]
                }
            ]
        
        # Convert to MutationInfo objects
        for mut_data in mutations:
            for mutant in mut_data.get("mutations", []):
                mutation = MutationInfo(
                    gene_name=mut_data["gene"],
                    organism=mut_data["organism"],
                    position=mut_data["position"],
                    wild_type=mut_data["wild_type"],
                    mutant=mutant,
                    mic_change=mut_data["mic_change"],
                    drugs=mut_data["drugs"],
                    pubmed_ids=mut_data.get("pmid", [])
                )
                self.mutations.append(mutation)
                
        print(f"Loaded {len(self.mutations)} resistance mutations")
        
    def build_kmer_index(self, k: int = 31):
        """Build k-mer index for genes"""
        print(f"Building {k}-mer index...")
        
        kmer_index = {}
        
        for gene_name, gene_info in self.genes.items():
            seq = gene_info.sequence
            
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                # Simple hash function
                kmer_hash = hash(kmer) & 0xFFFFFFFF
                
                if kmer_hash not in kmer_index:
                    kmer_index[kmer_hash] = []
                    
                kmer_index[kmer_hash].append({
                    'gene_id': gene_info.gene_id,
                    'position': i
                })
                
        print(f"Created index with {len(kmer_index)} unique k-mers")
        return kmer_index
        
    def save_binary_database(self):
        """Save database in GPU-friendly binary format"""
        print("Saving binary database...")
        
        # Save genes
        with open(os.path.join(self.output_dir, "genes.bin"), "wb") as f:
            # Header: number of genes
            f.write(struct.pack("I", len(self.genes)))
            
            # Gene data
            for gene_name, gene_info in self.genes.items():
                # Gene ID (4 bytes)
                f.write(struct.pack("I", gene_info.gene_id))
                
                # Gene name (32 bytes, padded)
                name_bytes = gene_name.encode('utf-8')[:32]
                name_bytes += b'\0' * (32 - len(name_bytes))
                f.write(name_bytes)
                
                # Sequence length (4 bytes)
                f.write(struct.pack("I", gene_info.length))
                
                # QRDR region (8 bytes)
                f.write(struct.pack("II", gene_info.qrdr_start, gene_info.qrdr_end))
                
                # Sequence data (2-bit encoded)
                self._write_encoded_sequence(f, gene_info.sequence)
                
        # Save mutations
        with open(os.path.join(self.output_dir, "mutations.bin"), "wb") as f:
            # Header: number of mutations
            f.write(struct.pack("I", len(self.mutations)))
            
            # Mutation data
            for mutation in self.mutations:
                # Find gene ID
                gene_id = self.genes[mutation.gene_name].gene_id
                
                # Gene ID (4 bytes)
                f.write(struct.pack("I", gene_id))
                
                # Position (4 bytes)
                f.write(struct.pack("I", mutation.position))
                
                # Wild type and mutant (2 bytes)
                f.write(struct.pack("BB", ord(mutation.wild_type), ord(mutation.mutant)))
                
                # MIC change (4 bytes float)
                f.write(struct.pack("f", mutation.mic_change))
                
        print(f"Database saved to {self.output_dir}")
        
    def _write_encoded_sequence(self, f, sequence: str):
        """Write 2-bit encoded DNA sequence"""
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        # Pack 16 bases into each uint32
        packed = []
        for i in range(0, len(sequence), 16):
            chunk = sequence[i:i+16]
            value = 0
            for j, base in enumerate(chunk):
                value |= encoding.get(base, 0) << (j * 2)
            packed.append(value)
            
        # Write packed data
        for val in packed:
            f.write(struct.pack("I", val))
            
    def save_metadata(self):
        """Save metadata as JSON for reference"""
        metadata = {
            "version": "1.0",
            "genes": {name: asdict(info) for name, info in self.genes.items()},
            "mutations": [asdict(m) for m in self.mutations],
            "statistics": {
                "num_genes": len(self.genes),
                "num_mutations": len(self.mutations),
                "organisms": list(set(m.organism for m in self.mutations))
            }
        }
        
        with open(os.path.join(self.output_dir, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=2)
            
    def create_test_fastq(self, num_reads: int = 10000):
        """Create test FASTQ with known mutations"""
        print(f"Creating test FASTQ with {num_reads} reads...")
        
        output_file = os.path.join(self.output_dir, "test_reads.fastq")
        
        with open(output_file, "w") as f:
            # Add some reads with known mutations
            
            # 1. Wild type gyrA region
            wt_seq = "ATCGATCGTCTTCTGATCGATCG"  # Contains TCT (Ser83)
            for i in range(100):
                f.write(f"@read_wt_{i}\n")
                f.write(wt_seq + "\n")
                f.write("+\n")
                f.write("I" * len(wt_seq) + "\n")
                
            # 2. S83L mutant
            mut_seq = "ATCGATCGTTGTGATCGATCG"  # Contains TTG (Leu83)
            for i in range(50):
                f.write(f"@read_s83l_{i}\n")
                f.write(mut_seq + "\n")
                f.write("+\n")
                f.write("I" * len(mut_seq) + "\n")
                
            # 3. Random sequences
            import random
            bases = "ACGT"
            for i in range(num_reads - 150):
                seq = ''.join(random.choice(bases) for _ in range(150))
                f.write(f"@read_random_{i}\n")
                f.write(seq + "\n")
                f.write("+\n")
                f.write("I" * 150 + "\n")
                
        print(f"Test FASTQ saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Build BioGPU database")
    parser.add_argument("-o", "--output", default="data/biogpu_db",
                        help="Output directory")
    parser.add_argument("-m", "--mutations", 
                        help="Mutations JSON file")
    parser.add_argument("--test-data", action="store_true",
                        help="Create test FASTQ file")
    args = parser.parse_args()
    
    builder = BiogpuDatabaseBuilder(args.output)
    
    # Build database
    builder.download_reference_genes()
    builder.load_mutation_database(args.mutations)
    builder.build_kmer_index()
    builder.save_binary_database()
    builder.save_metadata()
    
    if args.test_data:
        builder.create_test_fastq()
        
    print("\nDatabase build complete!")
    print(f"Files created in: {args.output}")
    print("\nNext steps:")
    print("1. Use the binary files for GPU loading")
    print("2. Test with: biogpu run --db", args.output)

if __name__ == "__main__":
    main()