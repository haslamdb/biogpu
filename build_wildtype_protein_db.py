#!/usr/bin/env python3
"""
build_wildtype_protein_db.py

Build a wild-type protein reference database for mutation detection.
Only keeps canonical/reference sequences, removing variants to enable
direct mutation detection via Smith-Waterman alignment.
"""

import json
import os
import sys
import struct
import numpy as np
from collections import defaultdict
import pandas as pd
import argparse
from datetime import datetime
import hashlib

class WildTypeProteinDatabaseBuilder:
    def __init__(self, fq_genes_dir, output_dir):
        self.fq_genes_dir = fq_genes_dir
        self.output_dir = output_dir
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Database components
        self.proteins = []
        self.protein_kmers = defaultdict(list)  # k-mer -> list of (protein_id, position)
        self.species_map = {}
        self.gene_map = {}
        
        # K-mer parameters
        self.kmer_length = 5  # Match existing system
        
        # Statistics
        self.stats = {
            'total_sequences': 0,
            'unique_proteins': 0,
            'total_kmers': 0,
            'species_count': 0,
            'gene_count': 0
        }
        
        # Key fluoroquinolone resistance genes to prioritize
        self.priority_genes = {
            'gyrA', 'gyrB', 'parC', 'parE'
        }
    
    def build(self):
        """Main build process for wild-type database"""
        print(f"Building WILD-TYPE protein resistance database from {self.fq_genes_dir}")
        print(f"Output directory: {self.output_dir}")
        
        # Step 1: Load protein sequences and select wild-types
        self.load_wildtype_sequences()
        
        # Step 2: Build k-mer index
        self.build_kmer_index()
        
        # Step 3: Write binary database
        self.write_binary_database()
        
        # Step 4: Write metadata
        self.write_metadata()
        
        print("\nWild-type database build complete!")
        self.print_statistics()
    
    def load_wildtype_sequences(self):
        """Load protein sequences, selecting only wild-type/reference sequences"""
        print("\nLoading wild-type protein sequences...")
        
        species_id = 0
        gene_id = 0
        
        for species_dir in sorted(os.listdir(self.fq_genes_dir)):
            species_path = os.path.join(self.fq_genes_dir, species_dir)
            if not os.path.isdir(species_path):
                continue
            
            self.species_map[species_id] = species_dir
            gene_count_for_species = 0
            
            for gene_file in sorted(os.listdir(species_path)):
                if not gene_file.endswith('.json'):
                    continue
                
                gene_name = gene_file.replace('.json', '')
                
                # Skip if gene_id already exists (avoid duplicates)
                if gene_name not in [g for g in self.gene_map.values()]:
                    self.gene_map[gene_id] = gene_name
                else:
                    # Find existing gene_id
                    gene_id = next(k for k, v in self.gene_map.items() if v == gene_name)
                
                json_path = os.path.join(species_path, gene_file)
                
                try:
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                    
                    # Select the best wild-type sequence for this gene/species
                    wildtype_seq = self.select_wildtype_sequence(data, gene_name, species_dir)
                    
                    if wildtype_seq:
                        protein_entry = {
                            'id': len(self.proteins),
                            'species_id': species_id,
                            'gene_id': gene_id,
                            'sequence': wildtype_seq['sequence'],
                            'length': len(wildtype_seq['sequence']),
                            'gene_name': gene_name,
                            'species_name': species_dir,
                            'accession': wildtype_seq['accession'],
                            'is_wildtype': True,
                            'priority': gene_name in self.priority_genes
                        }
                        
                        self.proteins.append(protein_entry)
                        self.stats['total_sequences'] += 1
                        
                        print(f"  Added {species_dir} {gene_name}: {len(wildtype_seq['sequence'])} aa")
                    
                    gene_count_for_species += 1
                    
                except Exception as e:
                    print(f"  Warning: Could not process {json_path}: {e}")
                
                if gene_name not in [g for g in self.gene_map.values()]:
                    gene_id += 1
            
            if gene_count_for_species > 0:
                species_id += 1
        
        self.stats['species_count'] = len(self.species_map)
        self.stats['gene_count'] = len(self.gene_map)
        self.stats['unique_proteins'] = len(self.proteins)
        print(f"  Total: {self.stats['total_sequences']} wild-type protein sequences selected")
    
    def select_wildtype_sequence(self, data, gene_name, species_name):
        """Select the best wild-type sequence from multiple entries"""
        candidates = []
        
        for entry in data:
            if 'gene_features' not in entry:
                continue
            
            for feature in entry['gene_features']:
                if 'protein_sequence' not in feature:
                    continue
                
                protein_seq = feature['protein_sequence']
                
                # Quality filters
                if len(protein_seq) < 20 or 'X' in protein_seq or '*' in protein_seq[:-1]:
                    continue
                
                candidates.append({
                    'sequence': protein_seq.rstrip('*'),  # Remove terminal stop codon
                    'accession': entry['accession'],
                    'length': len(protein_seq),
                    'quality_score': self.score_sequence_quality(protein_seq, gene_name)
                })
        
        if not candidates:
            return None
        
        # Sort by quality score (higher is better)
        candidates.sort(key=lambda x: x['quality_score'], reverse=True)
        
        # For key resistance genes, prefer sequences closest to known lengths
        if gene_name in self.priority_genes:
            expected_lengths = self.get_expected_gene_lengths(gene_name)
            if expected_lengths:
                # Re-score based on length similarity
                for candidate in candidates:
                    length_score = max(1.0 - abs(candidate['length'] - exp_len) / exp_len 
                                     for exp_len in expected_lengths)
                    candidate['quality_score'] += length_score * 10
                
                candidates.sort(key=lambda x: x['quality_score'], reverse=True)
        
        return candidates[0]
    
    def score_sequence_quality(self, sequence, gene_name):
        """Score sequence quality for wild-type selection"""
        score = 0.0
        
        # Length score (reasonable protein lengths)
        if 50 <= len(sequence) <= 1000:
            score += 10.0
        elif 20 <= len(sequence) <= 50 or 1000 < len(sequence) <= 2000:
            score += 5.0
        
        # Composition score (avoid too many rare amino acids)
        rare_aas = sequence.count('W') + sequence.count('C') + sequence.count('M')
        if rare_aas / len(sequence) < 0.1:
            score += 5.0
        
        # Priority gene bonus
        if gene_name in self.priority_genes:
            score += 20.0
        
        # Penalize sequences with unusual features
        if sequence.count('*') > 1:  # Multiple stop codons
            score -= 10.0
        if 'XX' in sequence:  # Consecutive unknowns
            score -= 5.0
        
        return score
    
    def get_expected_gene_lengths(self, gene_name):
        """Get expected protein lengths for known genes (approximate)"""
        expected_lengths = {
            'gyrA': [875, 900],  # ~880 residues typical
            'gyrB': [800, 850],  # ~820 residues typical  
            'parC': [750, 800],  # ~780 residues typical
            'parE': [600, 650],  # ~630 residues typical
            'acrA': [350, 400],  # ~380 residues typical
            'acrB': [1000, 1100],  # ~1050 residues typical
            'tolC': [450, 500],  # ~470 residues typical
        }
        return expected_lengths.get(gene_name, [])
    
    def build_kmer_index(self, k=5):
        """Build k-mer index for protein sequences"""
        self.kmer_length = k
        print(f"\nBuilding {k}-mer index for wild-type sequences...")
        
        for protein_idx, protein in enumerate(self.proteins):
            seq = protein['sequence']
            
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                self.protein_kmers[kmer].append((protein_idx, i))
                self.stats['total_kmers'] += 1
        
        print(f"  Indexed {len(self.protein_kmers)} unique {k}-mers")
        print(f"  Total {k}-mer occurrences: {self.stats['total_kmers']}")
    
    def write_binary_database(self):
        """Write binary format compatible with existing GPU loader"""
        print("\nWriting wild-type binary database...")
        
        # Write protein sequences in simple format (number + concatenated sequences)
        with open(os.path.join(self.output_dir, 'proteins.bin'), 'wb') as f:
            # Header: number of proteins
            f.write(struct.pack('I', len(self.proteins)))
            
            # Write all sequences concatenated (this matches our fixed loader)
            for protein in self.proteins:
                f.write(protein['sequence'].encode('ascii'))
        
        # Write k-mer index (same format as original)
        with open(os.path.join(self.output_dir, 'protein_kmers.bin'), 'wb') as f:
            # Header: k-mer length and count
            f.write(struct.pack('I', self.kmer_length))
            f.write(struct.pack('I', len(self.protein_kmers)))
            
            # K-mer entries
            for kmer, positions in sorted(self.protein_kmers.items()):
                f.write(kmer.encode('ascii'))
                f.write(struct.pack('I', len(positions)))
                for protein_idx, pos in positions:
                    f.write(struct.pack('II', protein_idx, pos))
        
        print(f"  Wild-type binary database written to {self.output_dir}")
    
    def write_metadata(self):
        """Write metadata JSON file"""
        metadata = {
            'creation_date': datetime.now().isoformat(),
            'database_type': 'wild_type_references',
            'description': 'Wild-type protein sequences for mutation detection',
            'statistics': self.stats,
            'species_map': self.species_map,
            'gene_map': self.gene_map,
            'parameters': {
                'k': self.kmer_length,
                'source_dir': self.fq_genes_dir,
                'selection_method': 'quality_score_based'
            }
        }
        
        with open(os.path.join(self.output_dir, 'metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Write protein details for reference
        protein_details = []
        for protein in self.proteins:
            protein_details.append({
                'id': protein['id'],
                'gene': protein['gene_name'],
                'species': protein['species_name'],
                'length': protein['length'],
                'accession': protein['accession'],
                'sequence': protein['sequence'][:50] + '...' if len(protein['sequence']) > 50 else protein['sequence']
            })
        
        with open(os.path.join(self.output_dir, 'protein_details.json'), 'w') as f:
            json.dump(protein_details, f, indent=2)
    
    def print_statistics(self):
        """Print database statistics"""
        print("\n=== Wild-Type Database Statistics ===")
        print(f"Species: {self.stats['species_count']}")
        print(f"Genes: {self.stats['gene_count']}")
        print(f"Wild-type proteins: {self.stats['unique_proteins']}")
        print(f"Total k-mers: {self.stats['total_kmers']}")
        print(f"Unique k-mers: {len(self.protein_kmers)}")
        
        # Gene distribution
        gene_counts = defaultdict(int)
        for protein in self.proteins:
            gene_counts[protein['gene_name']] += 1
        
        print("\nGenes in wild-type database:")
        for gene, count in sorted(gene_counts.items()):
            priority = "⭐" if gene in self.priority_genes else ""
            print(f"  {gene}: {count} sequences {priority}")
        
        # Priority gene coverage
        priority_covered = sum(1 for gene in self.priority_genes if gene in gene_counts)
        print(f"\nFluoroquinolone resistance gene coverage: {priority_covered}/{len(self.priority_genes)}")
        print("Priority genes:", ", ".join(self.priority_genes))
        
        # Size information
        if os.path.exists(os.path.join(self.output_dir, 'proteins.bin')):
            proteins_size = os.path.getsize(os.path.join(self.output_dir, 'proteins.bin')) / 1024
            kmers_size = os.path.getsize(os.path.join(self.output_dir, 'protein_kmers.bin')) / 1024
            print(f"\nDatabase sizes:")
            print(f"  Proteins: {proteins_size:.1f} KB")
            print(f"  K-mer index: {kmers_size:.1f} KB")
            print(f"  Total: {proteins_size + kmers_size:.1f} KB")


def main():
    parser = argparse.ArgumentParser(description='Build wild-type protein database for mutation detection')
    parser.add_argument('fq_genes_dir', help='Directory containing species folders with gene JSON files')
    parser.add_argument('output_dir', help='Output directory for wild-type binary database')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isdir(args.fq_genes_dir):
        print(f"Error: {args.fq_genes_dir} is not a directory")
        sys.exit(1)
    
    # Build wild-type database
    builder = WildTypeProteinDatabaseBuilder(args.fq_genes_dir, args.output_dir)
    builder.build()


if __name__ == '__main__':
    main()