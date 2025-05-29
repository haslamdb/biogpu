#!/usr/bin/env python3
"""
BioGPU Genome Database Processor
Processes bacterial, fungal, and viral genomes to create a GPU-optimized k-mer index
with minimizer-based compression for memory efficiency
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
import json
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass, field
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import mmh3  # MurmurHash3 for fast hashing
from bitarray import bitarray
import h5py  # For hierarchical storage

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
class MinimizerConfig:
    """Configuration for minimizer selection"""
    kmer_size: int
    window_size: int
    hash_function: str = "murmurhash3"  # or "lexicographic"
    density: float = 1.0  # 1.0 = all minimizers, 0.5 = 50% subsampling
    
@dataclass
class MinimizerIndex:
    """Minimizer-based k-mer index structure"""
    config: MinimizerConfig
    minimizer_to_kmers: Dict[int, Set[int]]  # minimizer_hash -> set of k-mer hashes
    minimizer_to_genomes: Dict[int, List[int]]  # minimizer_hash -> list of genome_ids
    kmer_to_positions: Dict[int, List[Tuple[int, int]]]  # kmer_hash -> [(genome_id, position)]
    genome_info: List[GenomeInfo]
    
@dataclass
class HierarchicalIndex:
    """Hierarchical index supporting multiple memory footprints"""
    base_config: MinimizerConfig
    full_index: MinimizerIndex
    bloom_filters: Dict[str, 'BloomFilter']  # Different sizes/FP rates
    subsampled_indices: Dict[float, MinimizerIndex]  # density -> index

class BloomFilter:
    """Memory-efficient probabilistic data structure for k-mer membership testing"""
    def __init__(self, expected_items: int, fp_rate: float = 0.01):
        self.fp_rate = fp_rate
        self.size = self._optimal_size(expected_items, fp_rate)
        self.hash_count = self._optimal_hash_count(expected_items, self.size)
        self.bit_array = bitarray(self.size)
        self.bit_array.setall(0)
        
    def _optimal_size(self, n: int, p: float) -> int:
        """Calculate optimal bit array size"""
        return int(-n * np.log(p) / (np.log(2) ** 2))
    
    def _optimal_hash_count(self, n: int, m: int) -> int:
        """Calculate optimal number of hash functions"""
        return int((m / n) * np.log(2))
    
    def add(self, item: int):
        """Add an item to the bloom filter"""
        for i in range(self.hash_count):
            index = mmh3.hash(str(item), i) % self.size
            self.bit_array[index] = 1
    
    def contains(self, item: int) -> bool:
        """Check if an item might be in the set"""
        for i in range(self.hash_count):
            index = mmh3.hash(str(item), i) % self.size
            if not self.bit_array[index]:
                return False
        return True
    
    def save(self, filepath: str):
        """Save bloom filter to disk"""
        with open(filepath, 'wb') as f:
            pickle.dump({
                'fp_rate': self.fp_rate,
                'size': self.size,
                'hash_count': self.hash_count,
                'bit_array': self.bit_array
            }, f)
    
    @classmethod
    def load(cls, filepath: str) -> 'BloomFilter':
        """Load bloom filter from disk"""
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        bf = cls(1, 0.01)  # Dummy values
        bf.fp_rate = data['fp_rate']
        bf.size = data['size']
        bf.hash_count = data['hash_count']
        bf.bit_array = data['bit_array']
        return bf
    
class MinimizerSelector:
    """Efficient minimizer selection for sequences"""
    def __init__(self, config: MinimizerConfig):
        self.config = config
        
    def hash_kmer(self, kmer: str) -> int:
        """Hash k-mer using configured method"""
        if self.config.hash_function == "murmurhash3":
            return mmh3.hash(kmer, signed=False)
        else:  # lexicographic
            return hash(kmer)
    
    def select_minimizers(self, sequence: str) -> List[Tuple[int, int, str]]:
        """Select minimizers from sequence. Returns (position, hash, kmer)"""
        k = self.config.kmer_size
        w = self.config.window_size
        
        if len(sequence) < k:
            return []
        
        minimizers = []
        prev_minimizer = None
        
        for i in range(len(sequence) - k + 1):
            # Find minimizer in window starting at position i
            window_end = min(i + w, len(sequence) - k + 1)
            
            min_hash = float('inf')
            min_pos = -1
            min_kmer = ""
            
            for j in range(i, window_end):
                kmer = sequence[j:j + k]
                if 'N' in kmer:  # Skip k-mers with ambiguous bases
                    continue
                    
                kmer_hash = self.hash_kmer(kmer)
                if kmer_hash < min_hash:
                    min_hash = kmer_hash
                    min_pos = j
                    min_kmer = kmer
            
            # Add minimizer if it's new or different from previous
            if min_pos != -1 and (prev_minimizer is None or 
                                 min_pos != prev_minimizer[0]):
                minimizers.append((min_pos, min_hash, min_kmer))
                prev_minimizer = (min_pos, min_hash, min_kmer)
        
        # Apply density subsampling if needed
        if self.config.density < 1.0:
            n_keep = int(len(minimizers) * self.config.density)
            # Use hash to ensure consistent subsampling
            minimizers.sort(key=lambda x: x[1])
            minimizers = minimizers[:n_keep]
            minimizers.sort(key=lambda x: x[0])  # Restore position order
        
        return minimizers
    
class GenomeProcessor:
    def __init__(self, data_dir: str, config: MinimizerConfig):
        self.data_dir = Path(data_dir)
        self.config = config
        self.genomes = []
        self.minimizer_selector = MinimizerSelector(config)
        
        # Minimizer-based indexing
        self.minimizer_to_kmers = defaultdict(set)
        self.minimizer_to_genomes = defaultdict(list)
        self.kmer_to_positions = defaultdict(list)
        
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
    
    def process_genome_file(self, file_path: Path, genome_info: GenomeInfo, genome_id: int) -> Tuple[GenomeInfo, Dict[int, Dict[int, List[int]]]]:
        """Process a single genome file and extract minimizers with their k-mers"""
        # Structure: minimizer_hash -> {kmer_hash -> [positions]}
        local_minimizer_data = defaultdict(lambda: defaultdict(list))
        total_length = 0
        gc_count = 0
        
        try:
            with gzip.open(file_path, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    seq = str(record.seq).upper()
                    total_length += len(seq)
                    gc_count += seq.count('G') + seq.count('C')
                    
                    # Select minimizers
                    minimizers = self.minimizer_selector.select_minimizers(seq)
                    
                    # For each minimizer, collect associated k-mers
                    for min_pos, min_hash, min_kmer in minimizers:
                        # The minimizer itself is a k-mer
                        kmer_encoded = self.encode_kmer(min_kmer)
                        if kmer_encoded != -1:
                            local_minimizer_data[min_hash][kmer_encoded].append(min_pos)
                        
                        # Also collect k-mers in the window around this minimizer
                        window_start = max(0, min_pos - self.config.window_size + 1)
                        window_end = min(len(seq) - self.config.kmer_size + 1, 
                                       min_pos + self.config.window_size)
                        
                        for pos in range(window_start, window_end):
                            if pos != min_pos:  # Don't duplicate the minimizer
                                kmer = seq[pos:pos + self.config.kmer_size]
                                if 'N' not in kmer:
                                    kmer_encoded = self.encode_kmer(kmer)
                                    if kmer_encoded != -1:
                                        local_minimizer_data[min_hash][kmer_encoded].append(pos)
            
            genome_info.sequence_length = total_length
            genome_info.gc_content = gc_count / total_length if total_length > 0 else 0
            
            return genome_info, dict(local_minimizer_data)
            
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
                    genome_info, minimizer_data = future.result()
                    self.genomes.append(genome_info)
                    
                    # Update global minimizer index
                    for minimizer_hash, kmer_dict in minimizer_data.items():
                        # Track which genomes have this minimizer
                        if genome_id not in self.minimizer_to_genomes[minimizer_hash]:
                            self.minimizer_to_genomes[minimizer_hash].append(genome_id)
                        
                        # Track all k-mers associated with this minimizer
                        for kmer_hash, positions in kmer_dict.items():
                            self.minimizer_to_kmers[minimizer_hash].add(kmer_hash)
                            
                            # Track positions of k-mers
                            for pos in positions:
                                self.kmer_to_positions[kmer_hash].append((genome_id, pos))
                    
                    logger.info(f"Processed {genome_info.organism_name} ({genome_info.genome_type})")
                    
                except Exception as e:
                    logger.error(f"Failed to process genome {genome_id}: {e}")
    
    def create_minimizer_index(self) -> MinimizerIndex:
        """Create minimizer-based index"""
        logger.info("Creating minimizer index")
        
        return MinimizerIndex(
            config=self.config,
            minimizer_to_kmers=dict(self.minimizer_to_kmers),
            minimizer_to_genomes=dict(self.minimizer_to_genomes),
            kmer_to_positions=dict(self.kmer_to_positions),
            genome_info=self.genomes
        )
    
    def create_bloom_filters(self, base_index: MinimizerIndex) -> Dict[str, BloomFilter]:
        """Create bloom filters with different false positive rates"""
        logger.info("Creating bloom filters")
        
        bloom_filters = {}
        
        # Count total unique k-mers
        all_kmers = set()
        for kmer_set in base_index.minimizer_to_kmers.values():
            all_kmers.update(kmer_set)
        
        # Create filters with different FP rates for different memory constraints
        for name, fp_rate in [("high_accuracy", 0.001), ("balanced", 0.01), ("low_memory", 0.05)]:
            logger.info(f"Creating {name} bloom filter with FP rate {fp_rate}")
            bf = BloomFilter(len(all_kmers), fp_rate)
            
            for kmer_hash in all_kmers:
                bf.add(kmer_hash)
            
            bloom_filters[name] = bf
            size_mb = bf.size / 8 / 1024 / 1024
            logger.info(f"  {name} filter size: {size_mb:.2f} MB")
        
        return bloom_filters
    
    def create_subsampled_indices(self, base_index: MinimizerIndex) -> Dict[float, MinimizerIndex]:
        """Create subsampled indices for different memory footprints"""
        logger.info("Creating subsampled indices")
        
        subsampled = {}
        
        for density in [0.5, 0.25, 0.1]:
            logger.info(f"Creating {density*100}% density index")
            
            # Select subset of minimizers based on hash value for consistency
            minimizer_hashes = sorted(base_index.minimizer_to_kmers.keys())
            n_keep = int(len(minimizer_hashes) * density)
            
            # Use hash-based selection for deterministic subsampling
            selected_minimizers = sorted(minimizer_hashes)[:n_keep]
            selected_set = set(selected_minimizers)
            
            # Create new index with only selected minimizers
            sub_index = MinimizerIndex(
                config=MinimizerConfig(
                    kmer_size=base_index.config.kmer_size,
                    window_size=base_index.config.window_size,
                    hash_function=base_index.config.hash_function,
                    density=density
                ),
                minimizer_to_kmers={m: base_index.minimizer_to_kmers[m] 
                                  for m in selected_minimizers},
                minimizer_to_genomes={m: base_index.minimizer_to_genomes[m] 
                                    for m in selected_minimizers},
                kmer_to_positions=base_index.kmer_to_positions,  # Keep all k-mer positions
                genome_info=base_index.genome_info
            )
            
            subsampled[density] = sub_index
            
            # Report statistics
            n_kmers = sum(len(kmers) for kmers in sub_index.minimizer_to_kmers.values())
            logger.info(f"  {density*100}% index: {len(selected_minimizers)} minimizers, {n_kmers} k-mers")
        
        return subsampled
    
    def build_hierarchical_index(self) -> HierarchicalIndex:
        """Build complete hierarchical index with multiple access methods"""
        logger.info("Building hierarchical index")
        
        # Create base minimizer index
        base_index = self.create_minimizer_index()
        
        # Create bloom filters
        bloom_filters = self.create_bloom_filters(base_index)
        
        # Create subsampled indices
        subsampled_indices = self.create_subsampled_indices(base_index)
        
        return HierarchicalIndex(
            base_config=self.config,
            full_index=base_index,
            bloom_filters=bloom_filters,
            subsampled_indices=subsampled_indices
        )
    
    def save_hierarchical_index(self, output_dir: str):
        """Save the hierarchical index to disk"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Build the index
        h_index = self.build_hierarchical_index()
        
        # Save using HDF5 for efficient hierarchical storage
        logger.info(f"Saving hierarchical index to {output_path}")
        
        with h5py.File(output_path / "index.h5", 'w') as f:
            # Save configuration
            config_group = f.create_group('config')
            config_group.attrs['kmer_size'] = h_index.base_config.kmer_size
            config_group.attrs['window_size'] = h_index.base_config.window_size
            config_group.attrs['hash_function'] = h_index.base_config.hash_function
            config_group.attrs['density'] = h_index.base_config.density
            
            # Save genome info
            with open(output_path / "genome_info.pkl", 'wb') as gf:
                pickle.dump(h_index.full_index.genome_info, gf)
            
            # Save full index
            with open(output_path / "full_index.pkl", 'wb') as fi:
                pickle.dump({
                    'minimizer_to_kmers': h_index.full_index.minimizer_to_kmers,
                    'minimizer_to_genomes': h_index.full_index.minimizer_to_genomes,
                    'kmer_to_positions': h_index.full_index.kmer_to_positions
                }, fi)
            
            # Save bloom filters
            bf_group = f.create_group('bloom_filters')
            for name, bf in h_index.bloom_filters.items():
                bf.save(str(output_path / f"bloom_{name}.pkl"))
                bf_group.attrs[name] = str(output_path / f"bloom_{name}.pkl")
            
            # Save subsampled indices
            for density, sub_index in h_index.subsampled_indices.items():
                with open(output_path / f"index_density_{density}.pkl", 'wb') as sf:
                    pickle.dump({
                        'minimizer_to_kmers': sub_index.minimizer_to_kmers,
                        'minimizer_to_genomes': sub_index.minimizer_to_genomes,
                        'config': sub_index.config
                    }, sf)
        
        # Save index statistics
        stats = self.calculate_index_stats(h_index)
        with open(output_path / "index_stats.json", 'w') as f:
            json.dump(stats, f, indent=2)
        
        logger.info(f"Hierarchical index saved to {output_path}")
    
    def calculate_index_stats(self, h_index: HierarchicalIndex) -> Dict:
        """Calculate comprehensive index statistics"""
        stats = {
            'total_genomes': len(h_index.full_index.genome_info),
            'total_minimizers': len(h_index.full_index.minimizer_to_kmers),
            'total_unique_kmers': len(set().union(*h_index.full_index.minimizer_to_kmers.values())),
            'config': {
                'kmer_size': h_index.base_config.kmer_size,
                'window_size': h_index.base_config.window_size,
                'hash_function': h_index.base_config.hash_function
            },
            'bloom_filters': {},
            'subsampled_indices': {}
        }
        
        # Bloom filter stats
        for name, bf in h_index.bloom_filters.items():
            stats['bloom_filters'][name] = {
                'size_mb': bf.size / 8 / 1024 / 1024,
                'fp_rate': bf.fp_rate,
                'hash_functions': bf.hash_count
            }
        
        # Subsampled index stats
        for density, sub_index in h_index.subsampled_indices.items():
            n_kmers = sum(len(kmers) for kmers in sub_index.minimizer_to_kmers.values())
            stats['subsampled_indices'][f"{int(density*100)}%"] = {
                'minimizers': len(sub_index.minimizer_to_kmers),
                'kmers': n_kmers,
                'reduction_factor': len(h_index.full_index.minimizer_to_kmers) / len(sub_index.minimizer_to_kmers)
            }
        
        return stats

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
    """Main function to process genomes and build minimizer-based index"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Build GPU-optimized minimizer-based k-mer index for microbial genomes'
    )
    parser.add_argument('--data-dir', type=str, default='data', 
                       help='Base data directory')
    parser.add_argument('--kmer-size', type=int, default=31, 
                       help='K-mer size (default: 31)')
    parser.add_argument('--window-size', type=int, default=11, 
                       help='Minimizer window size (default: 11)')
    parser.add_argument('--hash-function', type=str, default='murmurhash3',
                       choices=['murmurhash3', 'lexicographic'],
                       help='Hash function for minimizer selection')
    parser.add_argument('--output-dir', type=str, default='data/indices/microbial_genomes',
                       help='Output directory for hierarchical index')
    parser.add_argument('--processes', type=int, default=None,
                       help='Number of parallel processes')
    parser.add_argument('--test-mode', action='store_true',
                       help='Process only a small subset for testing')
    
    args = parser.parse_args()
    
    # Create minimizer configuration
    config = MinimizerConfig(
        kmer_size=args.kmer_size,
        window_size=args.window_size,
        hash_function=args.hash_function,
        density=1.0  # Full density for base index
    )
    
    # Create output directory if needed
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Process genomes
    processor = GenomeProcessor(args.data_dir, config)
    
    if args.test_mode:
        logger.info("Running in test mode - processing limited data")
        # In test mode, process only a few files
        processor.process_all_genomes(args.processes)
    else:
        processor.process_all_genomes(args.processes)
    
    # Save hierarchical index
    processor.save_hierarchical_index(args.output_dir)
    
    # Display statistics
    stats_file = output_path / "index_stats.json"
    if stats_file.exists():
        with open(stats_file, 'r') as f:
            stats = json.load(f)
        
        logger.info("\n=== Index Statistics ===")
        logger.info(f"Total genomes: {stats['total_genomes']}")
        logger.info(f"Total minimizers: {stats['total_minimizers']}")
        logger.info(f"Total unique k-mers: {stats['total_unique_kmers']}")
        logger.info(f"K-mer size: {stats['config']['kmer_size']}")
        logger.info(f"Window size: {stats['config']['window_size']}")
        
        logger.info("\n=== Bloom Filter Sizes ===")
        for name, bf_stats in stats['bloom_filters'].items():
            logger.info(f"{name}: {bf_stats['size_mb']:.2f} MB (FP rate: {bf_stats['fp_rate']})")
        
        logger.info("\n=== Subsampled Index Sizes ===")
        for density, sub_stats in stats['subsampled_indices'].items():
            logger.info(f"{density}: {sub_stats['minimizers']} minimizers, "
                       f"{sub_stats['kmers']} k-mers (reduction: {sub_stats['reduction_factor']:.2f}x)")

if __name__ == "__main__":
    main()