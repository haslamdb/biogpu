#!/usr/bin/env python3
"""
Enhanced K-mer Index Builder for FQ Resistance Detection

This script builds and validates k-mer indices from downloaded GenBank sequences.
Includes extensive debugging and validation to ensure proper index creation.

Usage:
    python enhanced_kmer_builder.py sequences_dir mutations.csv output_dir
"""

import json
import os
import pandas as pd
import numpy as np
from pathlib import Path
import struct
import argparse
from collections import defaultdict, Counter
import hashlib
import h5py
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EnhancedKmerIndexBuilder:
    def __init__(self, kmer_length=15):
        self.kmer_length = kmer_length
        
        # QRDR regions (amino acid ranges)
        self.qrdr_regions = {
            'gyrA': {'start': 67, 'end': 106},   # E. coli numbering
            'parC': {'start': 68, 'end': 108},   # E. coli numbering  
            'gyrB': {'start': 400, 'end': 470},  # E. coli numbering
            'parE': {'start': 400, 'end': 470},  # E. coli numbering
            'grlA': {'start': 68, 'end': 108},   # S. aureus equivalent of parC
            'grlB': {'start': 400, 'end': 470}   # S. aureus equivalent of parE
        }
        
        # Statistics for validation
        self.stats = {
            'organisms_processed': 0,
            'genes_processed': 0,
            'sequences_processed': 0,
            'kmers_generated': 0,
            'unique_kmers': 0,
            'qrdr_sequences_found': 0
        }
    
    def load_sequences(self, sequences_dir):
        """Load all sequence data from JSON files with validation"""
        logger.info(f"Loading sequences from {sequences_dir}")
        
        sequence_data = defaultdict(lambda: defaultdict(list))
        
        for organism_dir in Path(sequences_dir).iterdir():
            if not organism_dir.is_dir():
                continue
                
            organism = organism_dir.name.replace('_', ' ')
            logger.info(f"  Loading {organism}...")
            
            organism_gene_count = 0
            
            for gene_file in organism_dir.glob('*.json'):
                gene = gene_file.stem
                
                try:
                    with open(gene_file, 'r') as f:
                        gene_sequences = json.load(f)
                    
                    if gene_sequences:  # Only count non-empty
                        sequence_data[organism][gene] = gene_sequences
                        organism_gene_count += 1
                        logger.info(f"    {gene}: {len(gene_sequences)} sequences")
                        self.stats['sequences_processed'] += len(gene_sequences)
                    else:
                        logger.warning(f"    {gene}: No sequences found")
                        
                except Exception as e:
                    logger.error(f"    Error loading {gene_file}: {e}")
            
            if organism_gene_count > 0:
                self.stats['organisms_processed'] += 1
                self.stats['genes_processed'] += organism_gene_count
            else:
                logger.warning(f"  No valid genes found for {organism}")
        
        logger.info(f"Loaded {self.stats['organisms_processed']} organisms, "
                   f"{self.stats['genes_processed']} genes, "
                   f"{self.stats['sequences_processed']} sequences")
        
        return sequence_data
    
    def validate_nucleotide_sequence(self, sequence):
        """Validate that sequence contains only valid nucleotides"""
        valid_bases = set('ATCGUatcgu')
        sequence_bases = set(sequence)
        invalid_bases = sequence_bases - valid_bases
        
        if invalid_bases:
            return False, f"Invalid bases found: {invalid_bases}"
        
        return True, "Valid"
    
    def extract_qrdr_sequence(self, gene_feature, gene_name):
        """Extract QRDR region from gene feature with validation"""
        protein_seq = gene_feature.get('protein_sequence', '')
        nucleotide_seq = gene_feature.get('nucleotide_sequence', '')
        
        if not protein_seq or not nucleotide_seq:
            return None, "Missing sequence data"
        
        if gene_name not in self.qrdr_regions:
            return None, f"Unknown gene: {gene_name}"
        
        # Validate nucleotide sequence
        is_valid, validation_msg = self.validate_nucleotide_sequence(nucleotide_seq)
        if not is_valid:
            return None, f"Invalid nucleotide sequence: {validation_msg}"
        
        qrdr_range = self.qrdr_regions[gene_name]
        start_pos = max(0, qrdr_range['start'] - 1)  # Convert to 0-based
        end_pos = min(len(protein_seq), qrdr_range['end'])
        
        if end_pos <= start_pos:
            return None, f"Invalid QRDR range: {start_pos}-{end_pos}"
        
        qrdr_aa_seq = protein_seq[start_pos:end_pos]
        
        # Map back to nucleotide coordinates
        nt_start = start_pos * 3
        nt_end = end_pos * 3
        
        if nt_end > len(nucleotide_seq):
            return None, f"Nucleotide sequence too short: {len(nucleotide_seq)} < {nt_end}"
        
        qrdr_nt_seq = nucleotide_seq[nt_start:nt_end]
        
        qrdr_info = {
            'aa_sequence': qrdr_aa_seq,
            'nt_sequence': qrdr_nt_seq.upper(),  # Ensure uppercase
            'aa_start': start_pos + 1,  # 1-based for output
            'aa_end': end_pos,
            'nt_start': nt_start,
            'nt_end': nt_end,
            'length_aa': len(qrdr_aa_seq),
            'length_nt': len(qrdr_nt_seq),
            'full_gene_length': len(nucleotide_seq)
        }
        
        return qrdr_info, "Success"
    
    def generate_kmers(self, sequence, k=15):
        """Generate all k-mers from a DNA sequence with validation"""
        sequence = sequence.upper().replace('U', 'T')  # Convert RNA to DNA
        
        if len(sequence) < k:
            return [], f"Sequence too short: {len(sequence)} < {k}"
        
        valid_bases = set('ATCG')
        kmers = []
        invalid_count = 0
        
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_bases = set(kmer)
            
            # Only include k-mers with valid bases
            if kmer_bases.issubset(valid_bases):
                kmers.append(kmer)
            else:
                invalid_count += 1
        
        status = f"Generated {len(kmers)} k-mers, {invalid_count} invalid"
        return kmers, status
    
    def encode_kmer(self, kmer):
        """Encode k-mer as 64-bit integer (2 bits per base)"""
        if len(kmer) > 32:  # Can't fit in 64 bits
            return None, "K-mer too long for encoding"
        
        encoded = 0
        base_map = {'A': 0, 'T': 3, 'C': 1, 'G': 2}
        
        for base in kmer:
            if base not in base_map:
                return None, f"Invalid base: {base}"
            encoded = (encoded << 2) | base_map[base]
        
        return encoded, "Success"
    
    def build_comprehensive_index(self, sequence_data, mutation_map):
        """Build comprehensive k-mer index with full metadata"""
        logger.info("Building comprehensive k-mer index...")
        
        # Data structures
        organism_db = []
        kmer_to_sequences = defaultdict(list)  # kmer -> list of (org_id, gene, seq_id, pos)
        all_sequences = []  # Store all sequences for reference
        
        organism_id = 0
        global_seq_id = 0
        
        for organism, genes in sequence_data.items():
            logger.info(f"  Processing {organism}...")
            
            organism_entry = {
                'organism_id': organism_id,
                'species_name': organism,
                'genes': {},
                'total_sequences': 0,
                'total_kmers': 0
            }
            
            for gene_name, sequences in genes.items():
                if not sequences:
                    continue
                
                logger.info(f"    {gene_name}: processing {len(sequences)} sequences")
                
                gene_entry = {
                    'gene_name': gene_name,
                    'sequences': [],
                    'mutations': mutation_map.get(organism, {}).get(gene_name, []),
                    'qrdr_sequences': [],
                    'total_kmers': 0
                }
                
                seq_count = 0
                
                for seq_idx, seq_data in enumerate(sequences):
                    for feature_idx, feature in enumerate(seq_data.get('gene_features', [])):
                        # Extract QRDR if this is a resistance gene
                        qrdr_info, qrdr_status = self.extract_qrdr_sequence(feature, gene_name)
                        
                        if qrdr_info:
                            self.stats['qrdr_sequences_found'] += 1
                            gene_entry['qrdr_sequences'].append(qrdr_info['aa_sequence'])
                        
                        # Get full gene sequence for k-mer generation
                        full_sequence = feature.get('nucleotide_sequence', '')
                        if not full_sequence:
                            logger.warning(f"      Sequence {seq_idx}/{feature_idx}: No nucleotide sequence")
                            continue
                        
                        # Validate and generate k-mers
                        kmers, kmer_status = self.generate_kmers(full_sequence, self.kmer_length)
                        if not kmers:
                            logger.warning(f"      Sequence {seq_idx}/{feature_idx}: {kmer_status}")
                            continue
                        
                        # Store sequence metadata
                        sequence_entry = {
                            'global_seq_id': global_seq_id,
                            'local_seq_id': seq_idx,
                            'feature_id': feature_idx,
                            'accession': seq_data.get('accession', ''),
                            'description': seq_data.get('description', ''),
                            'organism_id': organism_id,
                            'gene_name': gene_name,
                            'sequence': full_sequence,
                            'length': len(full_sequence),
                            'qrdr_info': qrdr_info,
                            'num_kmers': len(kmers)
                        }
                        
                        all_sequences.append(sequence_entry)
                        gene_entry['sequences'].append(sequence_entry)
                        
                        # Process k-mers
                        for kmer_pos, kmer in enumerate(kmers):
                            encoded_kmer, encode_status = self.encode_kmer(kmer)
                            if encoded_kmer is not None:
                                kmer_info = {
                                    'organism_id': organism_id,
                                    'gene_name': gene_name,
                                    'global_seq_id': global_seq_id,
                                    'local_seq_id': seq_idx,
                                    'feature_id': feature_idx,
                                    'position': kmer_pos,
                                    'kmer_sequence': kmer,
                                    'accession': seq_data.get('accession', ''),
                                    'qrdr_info': qrdr_info is not None
                                }
                                kmer_to_sequences[encoded_kmer].append(kmer_info)
                                self.stats['kmers_generated'] += 1
                                gene_entry['total_kmers'] += 1
                        
                        global_seq_id += 1
                        seq_count += 1
                
                organism_entry['total_sequences'] += seq_count
                organism_entry['genes'][gene_name] = gene_entry
                logger.info(f"      Processed {seq_count} sequences, {gene_entry['total_kmers']} k-mers")
            
            organism_db.append(organism_entry)
            organism_entry['total_kmers'] = sum(g['total_kmers'] for g in organism_entry['genes'].values())
            logger.info(f"    Total for {organism}: {organism_entry['total_sequences']} sequences, "
                       f"{organism_entry['total_kmers']} k-mers")
            organism_id += 1
        
        self.stats['unique_kmers'] = len(kmer_to_sequences)
        
        logger.info(f"Index building complete:")
        logger.info(f"  {len(organism_db)} organisms")
        logger.info(f"  {global_seq_id} total sequences") 
        logger.info(f"  {self.stats['kmers_generated']} total k-mers")
        logger.info(f"  {self.stats['unique_kmers']} unique k-mers")
        logger.info(f"  {self.stats['qrdr_sequences_found']} QRDR sequences")
        
        return organism_db, kmer_to_sequences, all_sequences
    
    def create_debug_reports(self, organism_db, kmer_to_sequences, all_sequences, output_dir):
        """Create comprehensive debug reports"""
        logger.info("Creating debug reports...")
        
        debug_dir = os.path.join(output_dir, 'debug')
        os.makedirs(debug_dir, exist_ok=True)
        
        # 1. K-mer frequency distribution
        kmer_counts = [len(sequences) for sequences in kmer_to_sequences.values()]
        kmer_freq_dist = Counter(kmer_counts)
        
        with open(os.path.join(debug_dir, 'kmer_frequency_distribution.txt'), 'w') as f:
            f.write("K-mer Frequency Distribution\n")
            f.write("===========================\n")
            f.write("Frequency\tCount\n")
            for freq in sorted(kmer_freq_dist.keys()):
                f.write(f"{freq}\t{kmer_freq_dist[freq]}\n")
        
        # 2. Organism/gene summary
        with open(os.path.join(debug_dir, 'organism_gene_summary.txt'), 'w') as f:
            f.write("Organism/Gene Summary\n")
            f.write("====================\n")
            for org in organism_db:
                f.write(f"\nOrganism: {org['species_name']} (ID: {org['organism_id']})\n")
                f.write(f"  Total sequences: {org['total_sequences']}\n")
                f.write(f"  Total k-mers: {org['total_kmers']}\n")
                for gene_name, gene_data in org['genes'].items():
                    f.write(f"  Gene {gene_name}:\n")
                    f.write(f"    Sequences: {len(gene_data['sequences'])}\n")
                    f.write(f"    K-mers: {gene_data['total_kmers']}\n")
                    f.write(f"    QRDR sequences: {len(gene_data['qrdr_sequences'])}\n")
                    f.write(f"    Mutations: {len(gene_data['mutations'])}\n")
        
        # 3. Sample k-mers for manual inspection
        with open(os.path.join(debug_dir, 'sample_kmers.txt'), 'w') as f:
            f.write("Sample K-mers for Manual Inspection\n")
            f.write("===================================\n")
            f.write("Encoded K-mer\tK-mer Sequence\tOccurrences\tFirst Organism\tFirst Gene\n")
            
            # Show first 50 k-mers
            for i, (encoded_kmer, sequences) in enumerate(list(kmer_to_sequences.items())[:50]):
                first_seq = sequences[0]
                f.write(f"{encoded_kmer}\t{first_seq['kmer_sequence']}\t{len(sequences)}\t")
                f.write(f"{first_seq['organism_id']}\t{first_seq['gene_name']}\n")
        
        # 4. QRDR sequence analysis
        with open(os.path.join(debug_dir, 'qrdr_analysis.txt'), 'w') as f:
            f.write("QRDR Sequence Analysis\n")
            f.write("======================\n")
            
            for org in organism_db:
                for gene_name, gene_data in org['genes'].items():
                    if gene_data['qrdr_sequences']:
                        f.write(f"\n{org['species_name']} - {gene_name}:\n")
                        qrdr_counter = Counter(gene_data['qrdr_sequences'])
                        for qrdr_seq, count in qrdr_counter.most_common():
                            f.write(f"  {qrdr_seq} (count: {count})\n")
        
        # 5. Index validation summary
        with open(os.path.join(debug_dir, 'validation_summary.txt'), 'w') as f:
            f.write("Index Validation Summary\n")
            f.write("=======================\n")
            f.write(f"Total organisms: {len(organism_db)}\n")
            f.write(f"Total sequences: {len(all_sequences)}\n")
            f.write(f"Total k-mers generated: {self.stats['kmers_generated']}\n")
            f.write(f"Unique k-mers: {self.stats['unique_kmers']}\n")
            f.write(f"QRDR sequences found: {self.stats['qrdr_sequences_found']}\n")
            f.write(f"K-mer length: {self.kmer_length}\n")
            
            # Check for potential issues
            f.write("\nPotential Issues:\n")
            if self.stats['qrdr_sequences_found'] == 0:
                f.write("- WARNING: No QRDR sequences found!\n")
            if self.stats['unique_kmers'] < 1000:
                f.write("- WARNING: Very few unique k-mers generated\n")
            if len(organism_db) == 0:
                f.write("- ERROR: No organisms processed!\n")
        
        logger.info(f"Debug reports saved to {debug_dir}")
    
    def save_binary_index(self, organism_db, kmer_to_sequences, all_sequences, output_dir):
        """Save index in binary format for GPU loading"""
        logger.info("Saving binary index...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Sort k-mers for binary search
        sorted_kmers = sorted(kmer_to_sequences.keys())
        
        # Create k-mer lookup table
        kmer_entries = []
        for encoded_kmer in sorted_kmers:
            for seq_info in kmer_to_sequences[encoded_kmer]:
                kmer_entries.append({
                    'kmer': encoded_kmer,
                    'gene_id': seq_info['organism_id'],  # Using organism_id as gene_id for now
                    'species_id': seq_info['organism_id'],
                    'seq_id': seq_info['global_seq_id'],
                    'position': seq_info['position']
                })
        
        # Save binary k-mer index
        kmer_path = os.path.join(output_dir, 'kmer_index.bin')
        with open(kmer_path, 'wb') as f:
            # Header
            f.write(struct.pack('I', len(kmer_entries)))  # Number of entries
            f.write(struct.pack('I', self.kmer_length))   # K-mer length
            
            # K-mer entries
            for entry in kmer_entries:
                f.write(struct.pack('Q', entry['kmer']))      # 64-bit k-mer
                f.write(struct.pack('I', entry['gene_id']))   # Gene ID
                f.write(struct.pack('I', entry['species_id'])) # Species ID  
                f.write(struct.pack('I', entry['seq_id']))    # Sequence ID
                f.write(struct.pack('H', entry['position']))  # Position
        
        logger.info(f"Saved binary k-mer index: {kmer_path} ({len(kmer_entries)} entries)")
        
        # Save sequence database
        seq_path = os.path.join(output_dir, 'sequences.bin')
        with open(seq_path, 'wb') as f:
            # Header
            f.write(struct.pack('I', len(all_sequences)))
            
            # Sequences
            for seq in all_sequences:
                seq_bytes = seq['sequence'].encode('utf-8')
                f.write(struct.pack('I', len(seq_bytes)))
                f.write(seq_bytes)
                f.write(struct.pack('I', seq['length']))
                f.write(struct.pack('I', seq['organism_id']))
                
                # Gene name
                gene_bytes = seq['gene_name'].encode('utf-8')
                f.write(struct.pack('I', len(gene_bytes)))
                f.write(gene_bytes)
        
        logger.info(f"Saved sequence database: {seq_path} ({len(all_sequences)} sequences)")
        
        # Save JSON metadata for human readability
        metadata = {
            'creation_date': pd.Timestamp.now().isoformat(),
            'kmer_length': self.kmer_length,
            'num_organisms': len(organism_db),
            'num_sequences': len(all_sequences),
            'num_kmer_entries': len(kmer_entries),
            'num_unique_kmers': len(sorted_kmers),
            'statistics': self.stats,
            'organisms': [{'id': org['organism_id'], 'name': org['species_name']} 
                         for org in organism_db]
        }
        
        metadata_path = os.path.join(output_dir, 'index_metadata.json')
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Saved metadata: {metadata_path}")
        
        return kmer_path, seq_path, metadata_path
    
    def validate_index(self, kmer_path, seq_path, metadata_path):
        """Validate the created index by reading it back"""
        logger.info("Validating created index...")
        
        try:
            # Read metadata
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
            
            # Read k-mer index
            with open(kmer_path, 'rb') as f:
                num_entries = struct.unpack('I', f.read(4))[0]
                kmer_length = struct.unpack('I', f.read(4))[0]
                
                logger.info(f"K-mer index: {num_entries} entries, k={kmer_length}")
                
                # Read first few entries
                for i in range(min(5, num_entries)):
                    kmer = struct.unpack('Q', f.read(8))[0]
                    gene_id = struct.unpack('I', f.read(4))[0]
                    species_id = struct.unpack('I', f.read(4))[0]
                    seq_id = struct.unpack('I', f.read(4))[0]
                    position = struct.unpack('H', f.read(2))[0]
                    
                    logger.info(f"  Entry {i}: kmer={kmer}, gene={gene_id}, species={species_id}, "
                               f"seq={seq_id}, pos={position}")
            
            # Read sequence database
            with open(seq_path, 'rb') as f:
                num_sequences = struct.unpack('I', f.read(4))[0]
                logger.info(f"Sequence database: {num_sequences} sequences")
                
                # Read first sequence
                if num_sequences > 0:
                    seq_len = struct.unpack('I', f.read(4))[0]
                    sequence = f.read(seq_len).decode('utf-8')
                    length = struct.unpack('I', f.read(4))[0]
                    organism_id = struct.unpack('I', f.read(4))[0]
                    gene_len = struct.unpack('I', f.read(4))[0]
                    gene_name = f.read(gene_len).decode('utf-8')
                    
                    logger.info(f"  First sequence: organism={organism_id}, gene={gene_name}, "
                               f"length={length}, seq_start={sequence[:50]}...")
            
            logger.info("Index validation PASSED")
            return True
            
        except Exception as e:
            logger.error(f"Index validation FAILED: {e}")
            return False
    
    def build_index(self, sequences_dir, mutations_csv, output_dir):
        """Main function to build the complete index"""
        logger.info("=== Enhanced K-mer Index Builder ===")
        
        # Load input data
        sequence_data = self.load_sequences(sequences_dir)
        if not sequence_data:
            logger.error("No sequence data loaded!")
            return False
        
        mutations_df = pd.read_csv(mutations_csv)
        mutation_map = self.map_mutation_positions(mutations_df)
        
        # Build comprehensive index
        organism_db, kmer_to_sequences, all_sequences = self.build_comprehensive_index(
            sequence_data, mutation_map
        )
        
        # Create debug reports
        self.create_debug_reports(organism_db, kmer_to_sequences, all_sequences, output_dir)
        
        # Save binary index
        kmer_path, seq_path, metadata_path = self.save_binary_index(
            organism_db, kmer_to_sequences, all_sequences, output_dir
        )
        
        # Validate index
        validation_success = self.validate_index(kmer_path, seq_path, metadata_path)
        
        if validation_success:
            logger.info("=== Index Building COMPLETED SUCCESSFULLY ===")
            logger.info(f"Files created:")
            logger.info(f"  K-mer index: {kmer_path}")
            logger.info(f"  Sequences: {seq_path}")
            logger.info(f"  Metadata: {metadata_path}")
            logger.info(f"  Debug reports: {os.path.join(output_dir, 'debug')}")
        else:
            logger.error("=== Index Building FAILED ===")
        
        return validation_success
    
    def map_mutation_positions(self, mutations_df):
        """Create mapping of mutation positions for each organism/gene"""
        mutation_map = defaultdict(lambda: defaultdict(list))
        
        for _, row in mutations_df.iterrows():
            organism = row['Species']
            gene = row['Gene']
            position = int(row['Mutation Position'])
            wild_type = row['Wild-type Amino Acid']
            mutant = row['Mutant Amino Acid']
            
            mutation_info = {
                'position': position,
                'wild_type': wild_type,
                'mutant': mutant
            }
            
            mutation_map[organism][gene].append(mutation_info)
        
        return mutation_map

def main():
    parser = argparse.ArgumentParser(description='Enhanced K-mer Index Builder')
    parser.add_argument('sequences_dir', help='Directory containing downloaded sequences')
    parser.add_argument('mutations_csv', help='Path to mutations CSV file')
    parser.add_argument('output_dir', help='Output directory for index files')
    parser.add_argument('--kmer-length', type=int, default=15, help='K-mer length (default: 15)')
    
    args = parser.parse_args()
    
    builder = EnhancedKmerIndexBuilder(args.kmer_length)
    success = builder.build_index(args.sequences_dir, args.mutations_csv, args.output_dir)
    
    if success:
        print(f"\n✅ Index building completed successfully!")
        print(f"Output directory: {args.output_dir}")
        print(f"Check the debug/ subdirectory for validation reports")
    else:
        print(f"\n❌ Index building failed!")
        print(f"Check the logs above for details")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
