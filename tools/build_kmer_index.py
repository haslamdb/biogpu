#!/usr/bin/env python3
"""
BioGPU Genome Database Processor
Processes bacterial, fungal, and viral genomes to create a GPU-optimized k-mer index
"""

import os
import gzip
import csv
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
import struct
import pickle
import hashlib
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class GenomeInfo:
    """Information about a genome"""
    organism_name: str
    accession: str
    taxonomy_id: str
    genome_type: str  # bacterial, fungal, or viral
    file_path: str
    sequence_length: int = 0
    gc_content: float = 0.0

@dataclass
class KmerIndex:
    """GPU-optimized k-mer index structure"""
    kmer_size: int
    kmer_to_genomes: Dict[int, List[int]]  # kmer_hash -> list of genome_ids
    genome_info: List[GenomeInfo]
    kmer_positions: Dict[int, List[Tuple[int, int]]]  # kmer_hash -> [(genome_id, position)]
    
class GenomeProcessor:
    def __init__(self, data_dir: str, kmer_size: int = 31):
        self.data_dir = Path(data_dir)
        self.kmer_size = kmer_size
        self.genomes = []
        self.kmer_index = defaultdict(list)
        self.kmer_positions = defaultdict(list)
        
    def load_metadata(self, genome_type: str) -> Dict[str, GenomeInfo]:
        """Load genome metadata from CSV file"""
        metadata_path = self.data_dir / "genomes" / genome_type / "genome_metadata.csv"
        if not metadata_path.exists():
            logger.warning(f"Metadata file not found: {metadata_path}")
            return {}
        
        genome_map = {}
        with open(metadata_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Adjust column names based on actual CSV structure
                genome_info = GenomeInfo(
                    organism_name=row.get('organism_name', row.get('species', '')),
                    accession=row.get('accession', row.get('genome_accession', '')),
                    taxonomy_id=row.get('taxonomy_id', row.get('taxid', '')),
                    genome_type=genome_type,
                    file_path=row.get('file_path', row.get('filename', ''))
                )
                genome_map[genome_info.file_path] = genome_info
        
        logger.info(f"Loaded {len(genome_map)} {genome_type} genome metadata entries")
        return genome_map
    
    def find_genome_files(self, genome_type: str) -> List[Path]:
        """Find all .fna.gz files in the genome directory"""
        genome_dir = self.data_dir / "genomes" / genome_type
        fna_files = list(genome_dir.rglob("*.fna.gz"))
        logger.info(f"Found {len(fna_files)} {genome_type} genome files")
        return fna_files
    
    def encode_kmer(self, kmer: str) -> int:
        """Encode k-mer to 2-bit representation (A=00, C=01, G=10, T=11)"""
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        encoded = 0
        for base in kmer.upper():
            if base in encoding:
                encoded = (encoded << 2) | encoding[base]
            else:
                return -1  # Invalid k-mer with non-ACGT base
        return encoded
    
    def process_genome_file(self, file_path: Path, genome_info: GenomeInfo, genome_id: int) -> Tuple[GenomeInfo, Dict[int, List[int]]]:
        """Process a single genome file and extract k-mers"""
        local_kmer_counts = defaultdict(list)
        total_length = 0
        gc_count = 0
        
        try:
            with gzip.open(file_path, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    seq = str(record.seq).upper()
                    total_length += len(seq)
                    gc_count += seq.count('G') + seq.count('C')
                    
                    # Extract k-mers
                    for i in range(len(seq) - self.kmer_size + 1):
                        kmer = seq[i:i + self.kmer_size]
                        kmer_hash = self.encode_kmer(kmer)
                        if kmer_hash != -1:  # Valid k-mer
                            local_kmer_counts[kmer_hash].append(i)
            
            genome_info.sequence_length = total_length
            genome_info.gc_content = gc_count / total_length if total_length > 0 else 0
            
            # Convert to genome-specific format
            genome_kmers = {}
            for kmer_hash, positions in local_kmer_counts.items():
                genome_kmers[kmer_hash] = positions
                
            return genome_info, genome_kmers
            
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")
            return genome_info, {}
    
    def process_all_genomes(self, n_processes: int = None):
        """Process all genome files in parallel"""
        if n_processes is None:
            n_processes = mp.cpu_count()
        
        all_files = []
        genome_id_counter = 0
        
        # Collect all genome files
        for genome_type in ['bacteria', 'fungi', 'viral']:
            metadata = self.load_metadata(genome_type)
            files = self.find_genome_files(genome_type)
            
            for file_path in files:
                # Match file to metadata
                file_name = file_path.name
                genome_info = metadata.get(file_name, GenomeInfo(
                    organism_name=f"Unknown_{genome_type}_{file_name}",
                    accession=file_name.split('.')[0],
                    taxonomy_id="unknown",
                    genome_type=genome_type,
                    file_path=str(file_path)
                ))
                all_files.append((file_path, genome_info, genome_id_counter))
                genome_id_counter += 1
        
        logger.info(f"Processing {len(all_files)} genome files using {n_processes} processes")
        
        # Process files in parallel
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            futures = {executor.submit(self.process_genome_file, *args): args[2] 
                      for args in all_files}
            
            for future in as_completed(futures):
                genome_id = futures[future]
                try:
                    genome_info, genome_kmers = future.result()
                    self.genomes.append(genome_info)
                    
                    # Update global k-mer index
                    for kmer_hash, positions in genome_kmers.items():
                        self.kmer_index[kmer_hash].append(genome_id)
                        for pos in positions:
                            self.kmer_positions[kmer_hash].append((genome_id, pos))
                    
                    logger.info(f"Processed {genome_info.organism_name} ({genome_info.genome_type})")
                    
                except Exception as e:
                    logger.error(f"Failed to process genome {genome_id}: {e}")
    
    def create_gpu_optimized_index(self) -> 'GPUKmerIndex':
        """Create GPU-optimized data structures"""
        logger.info("Creating GPU-optimized index structures")
        
        # Convert to arrays for GPU efficiency
        num_kmers = len(self.kmer_index)
        max_genomes_per_kmer = max(len(genomes) for genomes in self.kmer_index.values())
        
        # Create fixed-size arrays for GPU
        kmer_hashes = np.zeros(num_kmers, dtype=np.uint64)
        genome_counts = np.zeros(num_kmers, dtype=np.uint32)
        genome_ids = np.full((num_kmers, max_genomes_per_kmer), -1, dtype=np.int32)
        
        # Fill arrays
        for i, (kmer_hash, genome_list) in enumerate(self.kmer_index.items()):
            kmer_hashes[i] = kmer_hash
            genome_counts[i] = len(genome_list)
            genome_ids[i, :len(genome_list)] = genome_list
        
        return GPUKmerIndex(
            kmer_size=self.kmer_size,
            kmer_hashes=kmer_hashes,
            genome_counts=genome_counts,
            genome_ids=genome_ids,
            genome_info=self.genomes,
            total_genomes=len(self.genomes)
        )
    
    def save_index(self, output_path: str):
        """Save the processed index to disk"""
        gpu_index = self.create_gpu_optimized_index()
        gpu_index.save(output_path)
        logger.info(f"Saved index to {output_path}")

@dataclass
class GPUKmerIndex:
    """GPU-optimized k-mer index with fixed-size arrays"""
    kmer_size: int
    kmer_hashes: np.ndarray  # Shape: (num_kmers,)
    genome_counts: np.ndarray  # Shape: (num_kmers,)
    genome_ids: np.ndarray  # Shape: (num_kmers, max_genomes_per_kmer)
    genome_info: List[GenomeInfo]
    total_genomes: int
    
    def save(self, filepath: str):
        """Save index in binary format for fast loading"""
        with open(filepath, 'wb') as f:
            # Write header
            f.write(b'BGPU')  # Magic number
            f.write(struct.pack('I', 1))  # Version
            f.write(struct.pack('I', self.kmer_size))
            f.write(struct.pack('I', self.total_genomes))
            f.write(struct.pack('Q', len(self.kmer_hashes)))
            
            # Write arrays
            self.kmer_hashes.tofile(f)
            self.genome_counts.tofile(f)
            self.genome_ids.tofile(f)
            
            # Write genome info
            pickle.dump(self.genome_info, f)
    
    @classmethod
    def load(cls, filepath: str) -> 'GPUKmerIndex':
        """Load index from binary format"""
        with open(filepath, 'rb') as f:
            # Read header
            magic = f.read(4)
            if magic != b'BGPU':
                raise ValueError("Invalid index file format")
            
            version = struct.unpack('I', f.read(4))[0]
            kmer_size = struct.unpack('I', f.read(4))[0]
            total_genomes = struct.unpack('I', f.read(4))[0]
            num_kmers = struct.unpack('Q', f.read(8))[0]
            
            # Read arrays
            kmer_hashes = np.fromfile(f, dtype=np.uint64, count=num_kmers)
            genome_counts = np.fromfile(f, dtype=np.uint32, count=num_kmers)
            
            # Calculate max_genomes_per_kmer from file
            remaining_ints = (os.path.getsize(filepath) - f.tell() - 1000000) // 4  # Rough estimate
            max_genomes = remaining_ints // num_kmers
            genome_ids = np.fromfile(f, dtype=np.int32, count=num_kmers * max_genomes)
            genome_ids = genome_ids.reshape(num_kmers, max_genomes)
            
            # Read genome info
            genome_info = pickle.load(f)
            
        return cls(
            kmer_size=kmer_size,
            kmer_hashes=kmer_hashes,
            genome_counts=genome_counts,
            genome_ids=genome_ids,
            genome_info=genome_info,
            total_genomes=total_genomes
        )
    
    def get_stats(self) -> Dict:
        """Get index statistics"""
        return {
            'total_kmers': len(self.kmer_hashes),
            'total_genomes': self.total_genomes,
            'kmer_size': self.kmer_size,
            'avg_genomes_per_kmer': np.mean(self.genome_counts),
            'max_genomes_per_kmer': np.max(self.genome_counts),
            'index_size_mb': (self.kmer_hashes.nbytes + self.genome_counts.nbytes + self.genome_ids.nbytes) / 1024 / 1024
        }

def main():
    """Main function to process genomes and build index"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Build GPU-optimized k-mer index for microbial genomes')
    parser.add_argument('--data-dir', type=str, default='data', help='Base data directory')
    parser.add_argument('--kmer-size', type=int, default=31, help='K-mer size (default: 31)')
    parser.add_argument('--output', type=str, default='data/indices/microbial_genomes.bgpu', help='Output index file')
    parser.add_argument('--processes', type=int, default=None, help='Number of parallel processes')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process genomes
    processor = GenomeProcessor(args.data_dir, args.kmer_size)
    processor.process_all_genomes(args.processes)
    processor.save_index(args.output)
    
    # Load and display stats
    index = GPUKmerIndex.load(args.output)
    stats = index.get_stats()
    logger.info(f"Index statistics: {stats}")

if __name__ == "__main__":
    main()