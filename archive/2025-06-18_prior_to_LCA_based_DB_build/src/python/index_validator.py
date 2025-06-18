#!/usr/bin/env python3
"""
Index Validation and Testing Tool

This tool validates k-mer indices and tests them against synthetic and real data.
It helps debug the k-mer screening process before moving to GPU implementation.

Usage:
    python index_validator.py index_dir [--test-data data_dir] [--create-synthetic]
"""

import os
import struct
import json
import gzip
import random
import argparse
from collections import defaultdict, Counter
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class IndexValidator:
    def __init__(self, index_dir):
        self.index_dir = index_dir
        self.kmer_length = 15
        
        # Load index components
        self.metadata = None
        self.kmer_entries = []
        self.sequences = []
        self.load_index()
    
    def load_index(self):
        """Load the binary index files"""
        logger.info(f"Loading index from {self.index_dir}")
        
        # Load metadata
        metadata_path = os.path.join(self.index_dir, 'index_metadata.json')
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                self.metadata = json.load(f)
            logger.info(f"Metadata loaded: {self.metadata['num_kmer_entries']} k-mer entries")
            self.kmer_length = self.metadata['kmer_length']
        else:
            logger.warning("No metadata file found")
        
        # Load k-mer index
        kmer_path = os.path.join(self.index_dir, 'kmer_index.bin')
        if os.path.exists(kmer_path):
            with open(kmer_path, 'rb') as f:
                num_entries = struct.unpack('I', f.read(4))[0]
                kmer_length = struct.unpack('I', f.read(4))[0]
                
                logger.info(f"Loading {num_entries} k-mer entries (k={kmer_length})")
                
                for i in range(num_entries):
                    entry = {
                        'kmer': struct.unpack('Q', f.read(8))[0],
                        'gene_id': struct.unpack('I', f.read(4))[0],
                        'species_id': struct.unpack('I', f.read(4))[0],
                        'seq_id': struct.unpack('I', f.read(4))[0],
                        'position': struct.unpack('H', f.read(2))[0]
                    }
                    self.kmer_entries.append(entry)
                
                logger.info(f"Loaded {len(self.kmer_entries)} k-mer entries")
        else:
            logger.error("K-mer index file not found!")
        
        # Load sequence database
        seq_path = os.path.join(self.index_dir, 'sequences.bin')
        if os.path.exists(seq_path):
            with open(seq_path, 'rb') as f:
                num_sequences = struct.unpack('I', f.read(4))[0]
                logger.info(f"Loading {num_sequences} sequences")
                
                for i in range(num_sequences):
                    seq_len = struct.unpack('I', f.read(4))[0]
                    sequence = f.read(seq_len).decode('utf-8')
                    length = struct.unpack('I', f.read(4))[0]
                    organism_id = struct.unpack('I', f.read(4))[0]
                    gene_len = struct.unpack('I', f.read(4))[0]
                    gene_name = f.read(gene_len).decode('utf-8')
                    
                    self.sequences.append({
                        'id': i,
                        'sequence': sequence,
                        'length': length,
                        'organism_id': organism_id,
                        'gene_name': gene_name
                    })
                
                logger.info(f"Loaded {len(self.sequences)} sequences")
        else:
            logger.error("Sequence database file not found!")
    
    def decode_kmer(self, encoded_kmer, k):
        """Decode k-mer from integer back to string"""
        base_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        kmer = ""
        
        for i in range(k):
            base_code = encoded_kmer & 3  # Get last 2 bits
            kmer = base_map[base_code] + kmer
            encoded_kmer >>= 2
        
        return kmer
    
    def encode_kmer(self, kmer):
        """Encode k-mer string to integer"""
        encoded = 0
        base_map = {'A': 0, 'T': 3, 'C': 1, 'G': 2}
        
        for base in kmer:
            if base not in base_map:
                return None
            encoded = (encoded << 2) | base_map[base]
        
        return encoded
    
    def validate_index_integrity(self):
        """Validate internal consistency of the index"""
        logger.info("Validating index integrity...")
        
        issues = []
        
        # Check k-mer uniqueness and sorting
        seen_kmers = set()
        prev_kmer = 0
        
        for i, entry in enumerate(self.kmer_entries):
            kmer = entry['kmer']
            
            # Check if sorted
            if kmer < prev_kmer:
                issues.append(f"K-mer index not sorted at position {i}")
            prev_kmer = kmer
            
            # Check for duplicates (same k-mer, same sequence)
            key = (kmer, entry['seq_id'], entry['position'])
            if key in seen_kmers:
                issues.append(f"Duplicate k-mer entry at position {i}")
            seen_kmers.add(key)
        
        # Check sequence ID references
        max_seq_id = len(self.sequences) - 1
        for i, entry in enumerate(self.kmer_entries):
            if entry['seq_id'] > max_seq_id:
                issues.append(f"Invalid sequence ID {entry['seq_id']} at k-mer entry {i}")
        
        # Check k-mer encoding/decoding
        for i in range(min(100, len(self.kmer_entries))):  # Test first 100
            entry = self.kmer_entries[i]
            decoded = self.decode_kmer(entry['kmer'], self.kmer_length)
            reencoded = self.encode_kmer(decoded)
            
            if reencoded != entry['kmer']:
                issues.append(f"K-mer encoding/decoding mismatch at entry {i}")
        
        if issues:
            logger.error(f"Index integrity issues found:")
            for issue in issues:
                logger.error(f"  - {issue}")
            return False
        else:
            logger.info("Index integrity validation PASSED")
            return True
    
    def analyze_kmer_distribution(self):
        """Analyze k-mer distribution across genes and organisms"""
        logger.info("Analyzing k-mer distribution...")
        
        # Count k-mers per gene
        gene_counts = Counter()
        organism_counts = Counter()
        sequence_counts = Counter()
        
        for entry in self.kmer_entries:
            gene_counts[entry['gene_id']] += 1
            organism_counts[entry['species_id']] += 1
            sequence_counts[entry['seq_id']] += 1
        
        logger.info("K-mer distribution by gene ID:")
        for gene_id, count in gene_counts.most_common():
            logger.info(f"  Gene {gene_id}: {count} k-mers")
        
        logger.info("K-mer distribution by organism ID:")
        for org_id, count in organism_counts.most_common():
            logger.info(f"  Organism {org_id}: {count} k-mers")
        
        # Find sequences with most k-mers
        logger.info("Top sequences by k-mer count:")
        for seq_id, count in sequence_counts.most_common(10):
            if seq_id < len(self.sequences):
                seq_info = self.sequences[seq_id]
                logger.info(f"  Seq {seq_id} ({seq_info['gene_name']}): {count} k-mers, "
                           f"length {seq_info['length']}")
        
        return gene_counts, organism_counts, sequence_counts
    
    def test_kmer_lookup(self, test_sequences):
        """Test k-mer lookup functionality"""
        logger.info(f"Testing k-mer lookup with {len(test_sequences)} test sequences...")
        
        # Build lookup table for fast searching
        kmer_lookup = defaultdict(list)
        for i, entry in enumerate(self.kmer_entries):
            kmer_lookup[entry['kmer']].append(i)
        
        total_test_kmers = 0
        found_kmers = 0
        matches_by_gene = defaultdict(int)
        
        for seq_idx, test_seq in enumerate(test_sequences):
            sequence = test_seq['sequence']
            expected_gene = test_seq.get('expected_gene', 'unknown')
            
            # Generate k-mers from test sequence
            for i in range(len(sequence) - self.kmer_length + 1):
                kmer = sequence[i:i+self.kmer_length]
                
                # Only test valid k-mers
                if all(base in 'ATCG' for base in kmer):
                    total_test_kmers += 1
                    encoded_kmer = self.encode_kmer(kmer)
                    
                    if encoded_kmer in kmer_lookup:
                        found_kmers += 1
                        
                        # Check which genes this k-mer matches
                        for entry_idx in kmer_lookup[encoded_kmer]:
                            entry = self.kmer_entries[entry_idx]
                            seq_info = self.sequences[entry['seq_id']]
                            matches_by_gene[seq_info['gene_name']] += 1
        
        hit_rate = found_kmers / total_test_kmers if total_test_kmers > 0 else 0
        
        logger.info(f"K-mer lookup results:")
        logger.info(f"  Total test k-mers: {total_test_kmers}")
        logger.info(f"  Found in index: {found_kmers}")
        logger.info(f"  Hit rate: {hit_rate:.3f}")
        
        logger.info("Matches by gene:")
        for gene, count in sorted(matches_by_gene.items(), key=lambda x: x[1], reverse=True):
            logger.info(f"  {gene}: {count} k-mer matches")
        
        return hit_rate, matches_by_gene
    
    def create_synthetic_test_data(self, output_dir, num_reads=1000):
        """Create synthetic test data based on sequences in the index"""
        logger.info(f"Creating synthetic test data with {num_reads} reads...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Group sequences by gene
        sequences_by_gene = defaultdict(list)
        for seq in self.sequences:
            sequences_by_gene[seq['gene_name']].append(seq)
        
        synthetic_reads = []
        read_length = 150
        
        for read_id in range(num_reads):
            # Randomly select a gene
            gene_name = random.choice(list(sequences_by_gene.keys()))
            sequences = sequences_by_gene[gene_name]
            
            # Randomly select a sequence from that gene
            source_seq = random.choice(sequences)
            full_sequence = source_seq['sequence']
            
            if len(full_sequence) < read_length:
                # If sequence is shorter than read length, use the whole sequence
                read_sequence = full_sequence
            else:
                # Extract random substring
                start_pos = random.randint(0, len(full_sequence) - read_length)
                read_sequence = full_sequence[start_pos:start_pos + read_length]
            
            # Add some sequencing errors (5% error rate)
            if random.random() < 0.05:
                # Introduce single base substitution
                pos = random.randint(0, len(read_sequence) - 1)
                bases = 'ATCG'
                new_base = random.choice([b for b in bases if b != read_sequence[pos]])
                read_sequence = read_sequence[:pos] + new_base + read_sequence[pos+1:]
            
            synthetic_reads.append({
                'id': f"synthetic_read_{read_id}",
                'sequence': read_sequence,
                'quality': 'I' * len(read_sequence),  # High quality
                'expected_gene': gene_name,
                'source_organism': source_seq['organism_id']
            })
        
        # Write FASTQ files
        for read_type in ['R1', 'R2']:
            fastq_path = os.path.join(output_dir, f'synthetic_reads_{read_type}.fastq.gz')
            with gzip.open(fastq_path, 'wt') as f:
                for read in synthetic_reads:
                    f.write(f"@{read['id']}/{read_type}\n")
                    f.write(f"{read['sequence']}\n")
                    f.write(f"+\n")
                    f.write(f"{read['quality']}\n")
        
        # Write metadata
        metadata_path = os.path.join(output_dir, 'synthetic_metadata.json')
        with open(metadata_path, 'w') as f:
            metadata = {
                'num_reads': num_reads,
                'read_length': read_length,
                'error_rate': 0.05,
                'gene_distribution': {gene: len([r for r in synthetic_reads if r['expected_gene'] == gene])
                                    for gene in sequences_by_gene.keys()}
            }
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Synthetic test data created in {output_dir}")
        logger.info(f"Gene distribution: {metadata['gene_distribution']}")
        
        return synthetic_reads
    
    def test_with_synthetic_data(self, synthetic_reads):
        """Test the index with synthetic data"""
        logger.info("Testing with synthetic data...")
        
        hit_rate, matches_by_gene = self.test_kmer_lookup(synthetic_reads)
        
        # Analyze results
        gene_distribution = Counter(read['expected_gene'] for read in synthetic_reads)
        
        logger.info("Expected vs found gene distribution:")
        for gene in gene_distribution:
            expected = gene_distribution[gene]
            found = matches_by_gene.get(gene, 0)
            logger.info(f"  {gene}: expected {expected} reads, found {found} k-mer matches")
        
        return hit_rate > 0.1  # Expect at least 10% hit rate
    
    def load_real_fastq(self, fastq_path, max_reads=1000):
        """Load real FASTQ data for testing"""
        logger.info(f"Loading real FASTQ data from {fastq_path}")
        
        reads = []
        
        if fastq_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        
        try:
            with open_func(fastq_path, mode) as f:
                while len(reads) < max_reads:
                    header = f.readline().strip()
                    if not header:
                        break
                    
                    sequence = f.readline().strip()
                    plus = f.readline().strip()
                    quality = f.readline().strip()
                    
                    if header.startswith('@') and plus.startswith('+'):
                        reads.append({
                            'id': header[1:],
                            'sequence': sequence,
                            'quality': quality
                        })
            
            logger.info(f"Loaded {len(reads)} reads from {fastq_path}")
            return reads
            
        except Exception as e:
            logger.error(f"Error loading FASTQ file: {e}")
            return []
    
    def test_with_real_data(self, fastq_path):
        """Test the index with real metagenomic data"""
        logger.info(f"Testing with real data: {fastq_path}")
        
        reads = self.load_real_fastq(fastq_path, max_reads=10000)
        if not reads:
            logger.error("No reads loaded from real data")
            return False
        
        hit_rate, matches_by_gene = self.test_kmer_lookup(reads)
        
        logger.info(f"Real data test results:")
        logger.info(f"  Hit rate: {hit_rate:.4f}")
        logger.info(f"  Total gene matches: {sum(matches_by_gene.values())}")
        
        return hit_rate > 0.001  # Expect at least 0.1% hit rate for real data
    
    def run_comprehensive_validation(self, test_data_dir=None, create_synthetic=False):
        """Run all validation tests"""
        logger.info("=== Running Comprehensive Index Validation ===")
        
        all_tests_passed = True
        
        # Test 1: Index integrity
        if not self.validate_index_integrity():
            all_tests_passed = False
        
        # Test 2: Analyze distribution
        self.analyze_kmer_distribution()
        
        # Test 3: Synthetic data
        if create_synthetic:
            synthetic_dir = os.path.join(self.index_dir, 'synthetic_test_data')
            synthetic_reads = self.create_synthetic_test_data(synthetic_dir)
            if not self.test_with_synthetic_data(synthetic_reads):
                logger.error("Synthetic data test FAILED")
                all_tests_passed = False
        
        # Test 4: Real data
        if test_data_dir:
            for fastq_file in ['VRE12_R1.fastq.gz', 'VRE12_R2.fastq.gz']:
                fastq_path = os.path.join(test_data_dir, fastq_file)
                if os.path.exists(fastq_path):
                    if not self.test_with_real_data(fastq_path):
                        logger.warning(f"Real data test with {fastq_file} had low hit rate")
                        # Don't fail here - real data might have low hit rates
        
        if all_tests_passed:
            logger.info("=== All validation tests PASSED ===")
        else:
            logger.error("=== Some validation tests FAILED ===")
        
        return all_tests_passed

def main():
    parser = argparse.ArgumentParser(description='Index Validation and Testing Tool')
    parser.add_argument('index_dir', help='Directory containing index files')
    parser.add_argument('--test-data', help='Directory containing test FASTQ files')
    parser.add_argument('--create-synthetic', action='store_true', 
                       help='Create synthetic test data')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.index_dir):
        print(f"Error: Index directory not found: {args.index_dir}")
        return 1
    
    validator = IndexValidator(args.index_dir)
    success = validator.run_comprehensive_validation(
        test_data_dir=args.test_data,
        create_synthetic=args.create_synthetic
    )
    
    if success:
        print(f"\n✅ Index validation completed successfully!")
    else:
        print(f"\n❌ Index validation found issues!")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
