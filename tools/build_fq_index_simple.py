#!/usr/bin/env python3
"""
Simplified fluoroquinolone resistance index builder.
Creates a single HDF5 file with all necessary data for GPU-based detection.
"""

import h5py
import numpy as np
import pandas as pd
from pathlib import Path
import json
import logging
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

class SimpleFQIndexBuilder:
    def __init__(self, db_dir: str):
        self.db_dir = Path(db_dir)
        self.markers = []
        self.sequences = []
        self.kmers = {}
        self.kmer_length = 15
        
    def load_mutations(self):
        """Load mutation data from CSV files"""
        logger.info("Loading mutation data...")
        
        # Load mutation tables
        mutation_file = self.db_dir / "resistance_mutations.csv"
        if not mutation_file.exists():
            # Try the original filename
            mutation_file = Path("tools/Quinolone_resistance_mutation_table.csv")
        
        if mutation_file.exists():
            df = pd.read_csv(mutation_file)
            logger.info(f"Loaded {len(df)} mutations from {mutation_file}")
            
            parsed_count = 0
            skipped_count = 0
            
            # Process each mutation
            for idx, row in df.iterrows():
                # Extract gene info
                gene_name = row.get('gene', '')
                species = row.get('species', '')
                
                # Build mutation string from location, wt, and mut columns
                try:
                    location = row.get('location', 0)
                    # Handle float/NaN values from pandas
                    if pd.isna(location) or not location:
                        continue
                    position = int(float(location))
                    ref_aa = row.get('wt', '')
                    mut_aa = row.get('mut', '')
                    
                    # Skip if any required field is missing
                    if not position or not ref_aa or not mut_aa:
                        continue
                    if pd.isna(ref_aa) or pd.isna(mut_aa):
                        continue
                        
                    mutation_str = f"{ref_aa}{position}{mut_aa}"
                    
                    # Create marker entry
                    marker = {
                        'gene': gene_name,
                        'species': species,
                        'mutation': mutation_str,
                        'position': position,
                        'ref_aa': ref_aa,
                        'mut_aa': mut_aa,
                        'marker_id': len(self.markers)
                    }
                    self.markers.append(marker)
                    parsed_count += 1
                    
                except (ValueError, TypeError):
                    skipped_count += 1
                    continue
                    
            logger.info(f"Parsed {parsed_count} valid mutations, skipped {skipped_count}")
        else:
            logger.warning(f"Mutation file not found: {mutation_file}")
            
    def load_reference_sequences(self):
        """Load reference sequences for each gene"""
        logger.info("Loading reference sequences...")
        
        # Try both possible locations
        ref_file = self.db_dir / "reference_sequences.fasta"
        if not ref_file.exists():
            ref_file = self.db_dir / "reference_sequences" / "reference_sequences.fasta"
        
        if not ref_file.exists():
            logger.warning(f"Reference file not found: {ref_file}")
            return
            
        logger.info(f"Loading sequences from: {ref_file}")
            
        # Group markers by gene and species
        gene_species_markers = {}
        for marker in self.markers:
            key = (marker['gene'].lower(), marker['species'].lower())
            if key not in gene_species_markers:
                gene_species_markers[key] = []
            gene_species_markers[key].append(marker)
        
        # Load sequences
        seq_count = 0
        matched_count = 0
        for record in SeqIO.parse(ref_file, "fasta"):
            # The full header is in the description field
            # Example: "Escherichia coli_parE_EFM0471628.1 Escherichia coli parE"
            # We need to find the structured part with underscores
            desc = record.description
            
            # Parse the description which has format:
            # "Species name_gene_accession Species name gene"
            # The species name appears twice - once with underscores, once with spaces
            
            # First, let's get the species with spaces from the end
            parts = desc.split()
            gene = None
            species = None
            
            # Find gene name in the parts
            for i, part in enumerate(parts):
                if part.lower() in ['gyra', 'gyrb', 'parc', 'pare', 'grla']:
                    gene = part.lower()
                    # Species is everything before the gene at the end
                    if i > 0:
                        species = ' '.join(parts[1:i]).lower()  # Skip the first part (ID)
                    break
            
            if not gene or not species:
                continue
            
            seq_count += 1
                    
                    key = (gene, species)
                    if seq_count <= 3:  # Debug first few
                        logger.debug(f"Sequence {seq_count}: {species} {gene}")
                        logger.debug(f"  Looking for key: {key}")
                        logger.debug(f"  Available keys sample: {list(gene_species_markers.keys())[:3]}")
                    
                    if key in gene_species_markers:
                        # Extract QRDR region based on gene
                        seq_str = str(record.seq).upper()
                        
                        # Define QRDR regions for each gene (approximate)
                        qrdr_regions = {
                            'gyra': (200, 500),   # Around codons 83-84
                            'gyrb': (400, 700),   # Around codon 466
                            'parc': (150, 450),   # Around codons 80-84
                            'pare': (350, 650)    # Around codons 420-458
                        }
                        
                        if gene in qrdr_regions:
                            qrdr_start, qrdr_end = qrdr_regions[gene]
                            qrdr_start = max(0, min(qrdr_start, len(seq_str) - 100))
                            qrdr_end = min(len(seq_str), qrdr_end)
                        else:
                            # Default region
                            qrdr_start = 100
                            qrdr_end = min(len(seq_str), 600)
                        
                        qrdr_seq = seq_str[qrdr_start:qrdr_end]
                        
                        # Store sequence for each marker
                        for marker in gene_species_markers[key]:
                            self.sequences.append({
                                'marker_id': marker['marker_id'],
                                'sequence': qrdr_seq,
                                'offset': qrdr_start,
                                'accession': record.id
                            })
                            
                            # Generate k-mers around mutation position
                            self._generate_kmers(marker, qrdr_seq, qrdr_start)
                            
                        matched_count += 1
                    else:
                        if seq_count < 5:  # Debug first few
                            logger.debug(f"No markers for {species} {gene}")
                        
        logger.info(f"Processed {seq_count} sequences, matched {matched_count} to markers")
        logger.info(f"Total markers: {len(self.markers)}, sequences stored: {len(self.sequences)}")
    
    def _generate_kmers(self, marker: Dict, sequence: str, seq_offset: int):
        """Generate k-mers around mutation position"""
        k = self.kmer_length
        mutation_pos = marker['position'] - seq_offset
        
        # Generate k-mers in a window around the mutation
        window = 50  # bases around mutation to consider
        start = max(0, mutation_pos - window)
        end = min(len(sequence), mutation_pos + window + k)
        
        for i in range(start, end - k + 1):
            kmer = sequence[i:i + k]
            if 'N' not in kmer:  # Skip k-mers with unknown bases
                # Encode k-mer as integer
                kmer_encoded = self._encode_kmer(kmer)
                
                if kmer_encoded not in self.kmers:
                    self.kmers[kmer_encoded] = []
                    
                # Add marker ID to this k-mer's list
                if marker['marker_id'] not in self.kmers[kmer_encoded]:
                    self.kmers[kmer_encoded].append(marker['marker_id'])
    
    def _encode_kmer(self, kmer: str) -> int:
        """Encode k-mer string as integer"""
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        value = 0
        for base in kmer:
            value = (value << 2) | encoding.get(base, 0)
        return value
    
    def _decode_kmer(self, value: int, k: int = 15) -> str:
        """Decode integer back to k-mer string (for debugging)"""
        decoding = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        kmer = []
        for _ in range(k):
            kmer.append(decoding[value & 3])
            value >>= 2
        return ''.join(reversed(kmer))
    
    def save_index(self, output_path: str):
        """Save simplified index as HDF5"""
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Saving index to {output_file}")
        
        with h5py.File(output_file, 'w') as f:
            # Metadata
            meta = f.create_group('metadata')
            meta.attrs['version'] = '1.0'
            meta.attrs['num_markers'] = len(self.markers)
            meta.attrs['num_sequences'] = len(self.sequences)
            meta.attrs['num_kmers'] = len(self.kmers)
            meta.attrs['kmer_length'] = self.kmer_length
            
            # Markers table
            if self.markers:
                markers_grp = f.create_group('markers')
                markers_grp.create_dataset('marker_id', 
                    data=[m['marker_id'] for m in self.markers])
                markers_grp.create_dataset('gene_name', 
                    data=[m['gene'].encode() for m in self.markers])
                markers_grp.create_dataset('species', 
                    data=[m['species'].encode() for m in self.markers])
                markers_grp.create_dataset('mutation', 
                    data=[m['mutation'].encode() for m in self.markers])
                markers_grp.create_dataset('position', 
                    data=[m['position'] for m in self.markers])
                markers_grp.create_dataset('ref_aa', 
                    data=[m['ref_aa'].encode() for m in self.markers])
                markers_grp.create_dataset('mut_aa', 
                    data=[m['mut_aa'].encode() for m in self.markers])
            
            # Sequences
            if self.sequences:
                seq_grp = f.create_group('sequences')
                
                # Concatenate all sequences
                all_seqs = ''.join(s['sequence'] for s in self.sequences)
                seq_grp.create_dataset('sequence_data', data=all_seqs.encode())
                
                # Store offsets and lengths
                current_offset = 0
                marker_ids = []
                seq_starts = []
                seq_lengths = []
                
                for seq in self.sequences:
                    marker_ids.append(seq['marker_id'])
                    seq_starts.append(current_offset)
                    seq_lengths.append(len(seq['sequence']))
                    current_offset += len(seq['sequence'])
                
                seq_grp.create_dataset('marker_id', data=marker_ids)
                seq_grp.create_dataset('seq_start', data=seq_starts)
                seq_grp.create_dataset('seq_length', data=seq_lengths)
            
            # K-mers (sorted for binary search)
            if self.kmers:
                kmer_grp = f.create_group('kmers')
                
                # Sort k-mers
                sorted_kmers = sorted(self.kmers.items())
                
                # Store k-mer values
                kmer_values = [k for k, _ in sorted_kmers]
                kmer_grp.create_dataset('kmer_encoded', data=kmer_values, dtype='uint64')
                
                # Store marker IDs for each k-mer
                # Use ragged array format: store all IDs concatenated + offsets
                all_marker_ids = []
                kmer_offsets = [0]
                
                for _, marker_ids in sorted_kmers:
                    all_marker_ids.extend(marker_ids)
                    kmer_offsets.append(len(all_marker_ids))
                
                kmer_grp.create_dataset('marker_ids', data=all_marker_ids)
                kmer_grp.create_dataset('marker_id_offsets', data=kmer_offsets)
                
                # Debug: store some k-mers as strings for verification
                debug_grp = kmer_grp.create_group('debug')
                sample_size = min(10, len(sorted_kmers))
                debug_kmers = []
                for i in range(sample_size):
                    kmer_int = sorted_kmers[i][0]
                    kmer_str = self._decode_kmer(kmer_int, self.kmer_length)
                    debug_kmers.append(kmer_str)
                debug_grp.create_dataset('sample_kmers', 
                    data=[k.encode() for k in debug_kmers])
        
        # Save summary
        summary = {
            'num_markers': len(self.markers),
            'num_sequences': len(self.sequences),
            'num_kmers': len(self.kmers),
            'kmer_length': self.kmer_length,
            'genes': list(set(m['gene'] for m in self.markers)),
            'species': list(set(m['species'] for m in self.markers)),
            'mutations': [m['mutation'] for m in self.markers[:10]]  # First 10
        }
        
        summary_file = output_file.with_suffix('.summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Index saved successfully!")
        logger.info(f"  Markers: {len(self.markers)}")
        logger.info(f"  Sequences: {len(self.sequences)}")
        logger.info(f"  K-mers: {len(self.kmers)}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Build simplified FQ resistance index')
    parser.add_argument('--db-dir', type=str, default='data/resistance_db',
                       help='Path to resistance database directory')
    parser.add_argument('--output', type=str, 
                       default='data/indices/fq_mutations/fq_simple_index.h5',
                       help='Output HDF5 file')
    
    args = parser.parse_args()
    
    # Build index
    builder = SimpleFQIndexBuilder(args.db_dir)
    
    # Load data
    builder.load_mutations()
    builder.load_reference_sequences()
    
    # Generate some test data if no real data found
    if not builder.markers:
        logger.warning("No markers loaded, creating test data...")
        # Add some test markers
        test_genes = ['gyrA', 'gyrA', 'parC', 'parC']
        test_species = ['Escherichia coli', 'Staphylococcus aureus', 
                       'Escherichia coli', 'Klebsiella pneumoniae']
        test_mutations = ['S83L', 'S84L', 'S80I', 'E84K']
        
        for i, (gene, species, mutation) in enumerate(zip(test_genes, test_species, test_mutations)):
            marker = {
                'gene': gene,
                'species': species,
                'mutation': mutation,
                'position': int(mutation[1:-1]),
                'ref_aa': mutation[0],
                'mut_aa': mutation[-1],
                'marker_id': i
            }
            builder.markers.append(marker)
            
            # Add test sequence
            test_seq = 'A' * 50 + 'TCGATCGATCGATCG' * 10 + 'T' * 50
            builder.sequences.append({
                'marker_id': i,
                'sequence': test_seq,
                'offset': 0
            })
            
            # Generate test k-mers
            builder._generate_kmers(marker, test_seq, 0)
    
    # Save index
    builder.save_index(args.output)
    logger.info("Index building complete!")

if __name__ == "__main__":
    main()