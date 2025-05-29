#!/usr/bin/env python3
"""
Convert downloaded GenBank sequences to GPU-ready database format

This script processes the JSON files created by the NCBI downloader and:
1. Maps amino acid positions from CSV to nucleotide coordinates
2. Extracts QRDR regions for each organism
3. Builds k-mer indices for GPU screening
4. Creates binary database files for CUDA kernels

Usage:
    python genbank_to_gpu_db.py sequences_dir mutations.csv output_dir
"""

import json
import os
import pandas as pd
import numpy as np
from pathlib import Path
import struct
import argparse
from collections import defaultdict
import hashlib

class GenBankToGPUConverter:
    def __init__(self):
        self.kmer_length = 15
        
        # QRDR regions (approximate amino acid ranges)
        self.qrdr_regions = {
            'gyrA': {'start': 67, 'end': 106},   # E. coli numbering
            'parC': {'start': 68, 'end': 108},   # E. coli numbering  
            'gyrB': {'start': 400, 'end': 470},  # E. coli numbering
            'parE': {'start': 400, 'end': 470}   # E. coli numbering
        }
    
    def load_sequences(self, sequences_dir):
        """Load all sequence data from JSON files"""
        print(f"Loading sequences from {sequences_dir}")
        
        sequence_data = defaultdict(lambda: defaultdict(list))
        
        for organism_dir in Path(sequences_dir).iterdir():
            if not organism_dir.is_dir():
                continue
                
            organism = organism_dir.name.replace('_', ' ')
            print(f"  Loading {organism}...")
            
            for gene_file in organism_dir.glob('*.json'):
                gene = gene_file.stem
                
                with open(gene_file, 'r') as f:
                    gene_sequences = json.load(f)
                
                sequence_data[organism][gene] = gene_sequences
                print(f"    {gene}: {len(gene_sequences)} sequences")
        
        return sequence_data
    
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
    
    def extract_qrdr_sequence(self, gene_feature, gene_name):
        """Extract QRDR region from gene feature"""
        protein_seq = gene_feature.get('protein_sequence', '')
        
        if not protein_seq or gene_name not in self.qrdr_regions:
            return None
        
        qrdr_range = self.qrdr_regions[gene_name]
        start_pos = max(0, qrdr_range['start'] - 1)  # Convert to 0-based
        end_pos = min(len(protein_seq), qrdr_range['end'])
        
        if end_pos <= start_pos:
            return None
        
        qrdr_sequence = protein_seq[start_pos:end_pos]
        
        # Map back to nucleotide coordinates
        nt_start = start_pos * 3
        nt_end = end_pos * 3
        
        nucleotide_seq = gene_feature.get('nucleotide_sequence', '')
        qrdr_nt_seq = nucleotide_seq[nt_start:nt_end] if nt_end <= len(nucleotide_seq) else ''
        
        return {
            'aa_sequence': qrdr_sequence,
            'nt_sequence': qrdr_nt_seq,
            'aa_start': start_pos + 1,  # 1-based for output
            'aa_end': end_pos,
            'nt_start': nt_start,
            'nt_end': nt_end,
            'length_aa': len(qrdr_sequence),
            'length_nt': len(qrdr_nt_seq)
        }
    
    def generate_kmers(self, sequence, k=15):
        """Generate all k-mers from a DNA sequence"""
        if len(sequence) < k:
            return []
        
        kmers = []
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k].upper()
            # Only include k-mers with valid bases
            if all(base in 'ATCG' for base in kmer):
                kmers.append(kmer)
        
        return kmers
    
    def encode_kmer(self, kmer):
        """Encode k-mer as 64-bit integer (2 bits per base)"""
        if len(kmer) > 32:  # Can't fit in 64 bits
            return None
        
        encoded = 0
        base_map = {'A': 0, 'T': 3, 'C': 1, 'G': 2}
        
        for base in kmer:
            if base not in base_map:
                return None
            encoded = (encoded << 2) | base_map[base]
        
        return encoded
    
    def build_organism_database(self, sequence_data, mutation_map):
        """Build organism-specific database with QRDR sequences and mutations"""
        print("\nBuilding organism database...")
        
        organism_db = []
        kmer_index = defaultdict(list)  # kmer -> [(organism_id, gene, sequence_id)]
        
        organism_id = 0
        
        for organism, genes in sequence_data.items():
            print(f"  Processing {organism}...")
            
            organism_entry = {
                'organism_id': organism_id,
                'species_name': organism,
                'genes': {}
            }
            
            for gene_name, sequences in genes.items():
                if not sequences:
                    continue
                
                print(f"    {gene_name}: {len(sequences)} sequences")
                
                gene_entry = {
                    'gene_name': gene_name,
                    'sequences': [],
                    'consensus_qrdr': None,
                    'mutations': mutation_map.get(organism, {}).get(gene_name, [])
                }
                
                qrdr_sequences = []
                
                for seq_idx, seq_data in enumerate(sequences):
                    for feature in seq_data.get('gene_features', []):
                        # Extract QRDR if this is a resistance gene
                        qrdr_info = None
                        if gene_name in self.qrdr_regions:
                            qrdr_info = self.extract_qrdr_sequence(feature, gene_name)
                        
                        if qrdr_info:
                            qrdr_sequences.append(qrdr_info['aa_sequence'])
                        
                        sequence_entry = {
                            'sequence_id': seq_idx,
                            'accession': seq_data.get('accession', ''),
                            'description': seq_data.get('description', ''),
                            'gene_feature': feature,
                            'qrdr_info': qrdr_info
                        }
                        
                        gene_entry['sequences'].append(sequence_entry)
                        
                        # Generate k-mers for screening
                        nt_sequence = feature.get('nucleotide_sequence', '')
                        if nt_sequence:
                            kmers = self.generate_kmers(nt_sequence, self.kmer_length)
                            for kmer in kmers:
                                encoded_kmer = self.encode_kmer(kmer)
                                if encoded_kmer is not None:
                                    kmer_index[encoded_kmer].append({
                                        'organism_id': organism_id,
                                        'gene': gene_name,
                                        'sequence_id': seq_idx,
                                        'kmer_sequence': kmer
                                    })
                
                # Build consensus QRDR sequence
                if qrdr_sequences:
                    gene_entry['consensus_qrdr'] = self.build_consensus_sequence(qrdr_sequences)
                
                organism_entry['genes'][gene_name] = gene_entry
            
            organism_db.append(organism_entry)
            organism_id += 1
        
        print(f"Built database with {len(organism_db)} organisms")
        print(f"Generated {len(kmer_index)} unique k-mers")
        
        return organism_db, kmer_index
    
    def build_consensus_sequence(self, sequences):
        """Build consensus sequence from multiple aligned sequences"""
        if not sequences:
            return None
        
        # Find the most common sequence (simple approach)
        # In production, you'd want proper multiple sequence alignment
        seq_counts = defaultdict(int)
        for seq in sequences:
            seq_counts[seq] += 1
        
        consensus = max(seq_counts.keys(), key=lambda k: seq_counts[k])
        confidence = seq_counts[consensus] / len(sequences)
        
        return {
            'sequence': consensus,
            'confidence': confidence,
            'variant_count': len(seq_counts),
            'total_sequences': len(sequences)
        }
    
    def save_binary_database(self, organism_db, kmer_index, output_dir):
        """Save database in binary format for GPU loading"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save organism database
        org_db_path = os.path.join(output_dir, 'organisms.bin')
        with open(org_db_path, 'wb') as f:
            # Header
            f.write(struct.pack('I', len(organism_db)))  # Number of organisms
            
            for org in organism_db:
                # Organism ID and name
                f.write(struct.pack('I', org['organism_id']))
                name_bytes = org['species_name'].encode('utf-8')
                f.write(struct.pack('I', len(name_bytes)))
                f.write(name_bytes)
                
                # Genes
                f.write(struct.pack('I', len(org['genes'])))
                for gene_name, gene_data in org['genes'].items():
                    gene_bytes = gene_name.encode('utf-8')
                    f.write(struct.pack('I', len(gene_bytes)))
                    f.write(gene_bytes)
                    
                    # Consensus QRDR sequence
                    if gene_data['consensus_qrdr']:
                        qrdr_seq = gene_data['consensus_qrdr']['sequence']
                        qrdr_bytes = qrdr_seq.encode('utf-8')
                        f.write(struct.pack('I', len(qrdr_bytes)))
                        f.write(qrdr_bytes)
                        f.write(struct.pack('f', gene_data['consensus_qrdr']['confidence']))
                    else:
                        f.write(struct.pack('I', 0))  # No QRDR sequence
                        f.write(struct.pack('f', 0.0))
                    
                    # Mutations
                    mutations = gene_data['mutations']
                    f.write(struct.pack('I', len(mutations)))
                    for mut in mutations:
                        f.write(struct.pack('I', mut['position']))
                        f.write(struct.pack('c', mut['wild_type'].encode('utf-8')))
                        f.write(struct.pack('c', mut['mutant'].encode('utf-8')))
        
        print(f"Saved organism database: {org_db_path}")
        
        # Save k-mer index
        kmer_db_path = os.path.join(output_dir, 'kmers.bin')
        with open(kmer_db_path, 'wb') as f:
            # Header
            f.write(struct.pack('I', len(kmer_index)))  # Number of unique k-mers
            f.write(struct.pack('I', self.kmer_length))  # K-mer length
            
            for encoded_kmer, entries in kmer_index.items():
                f.write(struct.pack('Q', encoded_kmer))  # 64-bit encoded k-mer
                f.write(struct.pack('I', len(entries)))   # Number of entries for this k-mer
                
                for entry in entries:
                    f.write(struct.pack('I', entry['organism_id']))
                    gene_bytes = entry['gene'].encode('utf-8')
                    f.write(struct.pack('I', len(gene_bytes)))
                    f.write(gene_bytes)
                    f.write(struct.pack('I', entry['sequence_id']))
        
        print(f"Saved k-mer index: {kmer_db_path}")
        
        # Save human-readable summary
        summary_path = os.path.join(output_dir, 'database_summary.json')
        summary = {
            'creation_date': str(pd.Timestamp.now()),
            'num_organisms': len(organism_db),
            'num_unique_kmers': len(kmer_index),
            'kmer_length': self.kmer_length,
            'organisms': []
        }
        
        for org in organism_db:
            org_summary = {
                'organism_id': org['organism_id'],
                'species_name': org['species_name'],
                'genes': {}
            }
            
            for gene_name, gene_data in org['genes'].items():
                org_summary['genes'][gene_name] = {
                    'num_sequences': len(gene_data['sequences']),
                    'num_mutations': len(gene_data['mutations']),
                    'has_consensus_qrdr': gene_data['consensus_qrdr'] is not None
                }
            
            summary['organisms'].append(org_summary)
        
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Saved database summary: {summary_path}")
    
    def convert(self, sequences_dir, mutations_csv, output_dir):
        """Main conversion function"""
        print("=== GenBank to GPU Database Converter ===")
        
        # Load data
        sequence_data = self.load_sequences(sequences_dir)
        mutations_df = pd.read_csv(mutations_csv)
        mutation_map = self.map_mutation_positions(mutations_df)
        
        # Build database
        organism_db, kmer_index = self.build_organism_database(sequence_data, mutation_map)
        
        # Save binary database
        self.save_binary_database(organism_db, kmer_index, output_dir)
        
        print(f"\nGPU database ready at: {output_dir}")
        return organism_db, kmer_index

def main():
    parser = argparse.ArgumentParser(description='Convert GenBank sequences to GPU database')
    parser.add_argument('sequences_dir', help='Directory containing downloaded sequences')
    parser.add_argument('mutations_csv', help='Path to mutations CSV file')
    parser.add_argument('output_dir', help='Output directory for GPU database')
    
    args = parser.parse_args()
    
    converter = GenBankToGPUConverter()
    organism_db, kmer_index = converter.convert(
        args.sequences_dir,
        args.mutations_csv, 
        args.output_dir
    )
    
    print("Conversion complete!")

if __name__ == "__main__":
    main()