#!/usr/bin/env python3
"""
build_protein_resistance_db.py

Build a protein reference database for GPU-accelerated translated search
of fluoroquinolone resistance genes.
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

# BLOSUM62 matrix for scoring
BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1, '*': -4},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1, '*': -4},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -2, '*': -4},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1, '*': -4},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1, '*': -4},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1, '*': -4},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -4},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1, '*': -4},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1, '*': -4},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -2, '*': -4},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': 0, '*': -4},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': 0, '*': -4},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -2, '*': -4},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1, '*': -4},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1, '*': -4},
    'B': {'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4},
    'Z': {'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4},
    'X': {'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1, '*': -4},
    '*': {'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4, 'X': -4, '*': 1}
}

# Known resistance positions for key genes (from literature)
RESISTANCE_POSITIONS = {
    'gyrA': {
        'Escherichia_coli': {83: ['S', 'LFW'], 87: ['D', 'NGY']},
        'Klebsiella_pneumoniae': {83: ['S', 'LFI'], 87: ['D', 'NGY']},
        'Pseudomonas_aeruginosa': {83: ['T', 'I'], 87: ['D', 'NY']},
        'Staphylococcus_aureus': {84: ['S', 'LFW'], 88: ['E', 'KV']},
        'Clostridioides_difficile': {82: ['T', 'IV'], 86: ['D', 'N']}
    },
    'parC': {
        'Escherichia_coli': {80: ['S', 'IR'], 84: ['E', 'VGK']},
        'Klebsiella_pneumoniae': {80: ['S', 'IR'], 84: ['E', 'VG']},
        'Pseudomonas_aeruginosa': {87: ['S', 'L'], 91: ['E', 'G']},
        'Staphylococcus_aureus': {80: ['S', 'FY'], 84: ['E', 'KVG']}
    },
    'gyrB': {
        'Escherichia_coli': {426: ['D', 'N'], 447: ['K', 'E']},
        'Pseudomonas_aeruginosa': {464: ['E', 'K'], 468: ['S', 'Y']}
    },
    'parE': {
        'Escherichia_coli': {416: ['L', 'F'], 420: ['S', 'R']},
        'Staphylococcus_aureus': {432: ['D', 'N'], 436: ['E', 'K']}
    }
}

class ProteinDatabaseBuilder:
    def __init__(self, fq_genes_dir, mutations_csv, output_dir, kmer_length=8):
        self.fq_genes_dir = fq_genes_dir
        self.mutations_csv = mutations_csv
        self.output_dir = output_dir
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Database components
        self.proteins = []
        self.protein_kmers = defaultdict(list)  # k-mer -> list of (protein_id, position)
        self.species_map = {}
        self.gene_map = {}
        self.accession_map = {}
        
        # K-mer parameters
        self.kmer_length = kmer_length  # K-mer size from command line
        
        # Statistics
        self.stats = {
            'total_sequences': 0,
            'unique_proteins': 0,
            'total_kmers': 0,
            'species_count': 0,
            'gene_count': 0
        }
    
    def build(self):
        """Main build process"""
        print(f"Building protein resistance database from {self.fq_genes_dir}")
        print(f"Output directory: {self.output_dir}")
        
        # Step 1: Load all protein sequences
        self.load_protein_sequences()
        
        # Step 2: Remove redundancy
        self.cluster_proteins()
        
        # Step 3: Build k-mer index
        self.build_kmer_index(self.kmer_length)
        
        # Step 4: Map resistance mutations
        self.map_resistance_mutations()
        
        # Step 5: Write binary database
        self.write_binary_database()
        
        # Step 6: Write metadata
        self.write_metadata()
        
        print("\nDatabase build complete!")
        self.print_statistics()
    
    def load_protein_sequences(self):
        """Load all protein sequences from JSON files"""
        print("\nLoading protein sequences...")
        
        species_id = 0
        gene_id = 0
        accession_id = 0
        
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
                self.gene_map[gene_id] = gene_name
                
                json_path = os.path.join(species_path, gene_file)
                
                try:
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                    
                    for entry in data:
                        if 'gene_features' not in entry:
                            continue
                        
                        for feature in entry['gene_features']:
                            if 'protein_sequence' not in feature:
                                continue
                            
                            protein_seq = feature['protein_sequence']
                            
                            # Skip if too short or contains invalid characters
                            if len(protein_seq) < 20 or 'X' in protein_seq:
                                continue
                            
                            self.accession_map[accession_id] = entry['accession']
                            
                            protein_entry = {
                                'id': len(self.proteins),
                                'species_id': species_id,
                                'gene_id': gene_id,
                                'accession_id': accession_id,
                                'sequence': protein_seq,
                                'length': len(protein_seq),
                                'gene_name': gene_name,
                                'species_name': species_dir,
                                'accession': entry['accession'],
                                'resistance_positions': {}
                            }
                            
                            # Add known resistance positions if available
                            if gene_name in RESISTANCE_POSITIONS:
                                if species_dir in RESISTANCE_POSITIONS[gene_name]:
                                    protein_entry['resistance_positions'] = RESISTANCE_POSITIONS[gene_name][species_dir]
                            
                            self.proteins.append(protein_entry)
                            accession_id += 1
                            self.stats['total_sequences'] += 1
                    
                    gene_count_for_species += 1
                    
                except Exception as e:
                    print(f"  Warning: Could not process {json_path}: {e}")
                
                gene_id += 1
            
            if gene_count_for_species > 0:
                print(f"  Loaded {gene_count_for_species} genes from {species_dir}")
                species_id += 1
        
        self.stats['species_count'] = len(self.species_map)
        self.stats['gene_count'] = len(self.gene_map)
        print(f"  Total: {self.stats['total_sequences']} protein sequences loaded")
    
    def cluster_proteins(self):
        """Remove redundant sequences (simple version - exact match)"""
        print("\nClustering proteins to remove redundancy...")
        
        unique_seqs = {}
        clustered_proteins = []
        
        for protein in self.proteins:
            seq_hash = hashlib.md5(protein['sequence'].encode()).hexdigest()
            
            if seq_hash not in unique_seqs:
                unique_seqs[seq_hash] = len(clustered_proteins)
                clustered_proteins.append(protein)
            else:
                # Merge resistance positions if same sequence
                existing_idx = unique_seqs[seq_hash]
                existing_positions = clustered_proteins[existing_idx]['resistance_positions']
                new_positions = protein['resistance_positions']
                
                for pos, (wt, muts) in new_positions.items():
                    if pos not in existing_positions:
                        existing_positions[pos] = (wt, muts)
        
        self.proteins = clustered_proteins
        self.stats['unique_proteins'] = len(self.proteins)
        print(f"  Reduced from {self.stats['total_sequences']} to {self.stats['unique_proteins']} unique proteins")
    
    def build_kmer_index(self, k=5):
        """Build k-mer index for protein sequences"""
        self.kmer_length = k  # Store the actual k-mer length used
        print(f"\nBuilding {k}-mer index...")
        
        for protein_idx, protein in enumerate(self.proteins):
            seq = protein['sequence']
            
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                self.protein_kmers[kmer].append((protein_idx, i))
                self.stats['total_kmers'] += 1
        
        print(f"  Indexed {len(self.protein_kmers)} unique {k}-mers")
        print(f"  Total {k}-mer occurrences: {self.stats['total_kmers']}")
    
    def map_resistance_mutations(self):
        """Map known resistance mutations from CSV"""
        print("\nMapping resistance mutations from CSV...")
        
        if os.path.exists(self.mutations_csv):
            try:
                df = pd.read_csv(self.mutations_csv)
                # Process mutations CSV if you have specific format
                print(f"  Loaded {len(df)} mutation records")
            except Exception as e:
                print(f"  Warning: Could not load mutations CSV: {e}")
        else:
            print(f"  Using built-in resistance positions")
    
    def write_binary_database(self):
        """Write binary format for GPU"""
        print("\nWriting binary database...")
        
        # Write protein sequences
        with open(os.path.join(self.output_dir, 'proteins.bin'), 'wb') as f:
            # Header: number of proteins
            f.write(struct.pack('I', len(self.proteins)))
            
            # Protein entries
            for protein in self.proteins:
                # Write fixed-size metadata
                f.write(struct.pack('I', protein['id']))
                f.write(struct.pack('I', protein['species_id']))
                f.write(struct.pack('I', protein['gene_id']))
                f.write(struct.pack('I', protein['accession_id']))
                f.write(struct.pack('I', protein['length']))
                
                # Write sequence (null-terminated)
                f.write(protein['sequence'].encode('ascii') + b'\0')
                
                # Write resistance positions
                f.write(struct.pack('I', len(protein['resistance_positions'])))
                for pos, (wt, muts) in protein['resistance_positions'].items():
                    f.write(struct.pack('I', pos))
                    f.write(wt.encode('ascii'))
                    f.write(struct.pack('I', len(muts)))
                    f.write(muts.encode('ascii'))
        
        # Write k-mer index
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
        
        # Write BLOSUM matrix as binary
        self.write_blosum_matrix()
        
        print(f"  Binary database written to {self.output_dir}")
    
    def write_blosum_matrix(self):
        """Write BLOSUM62 matrix in GPU-friendly format"""
        aa_order = 'ARNDCQEGHILKMFPSTWYVBZX*'
        matrix_size = len(aa_order)
        
        # Create numpy array
        blosum_array = np.zeros((matrix_size, matrix_size), dtype=np.float32)
        
        for i, aa1 in enumerate(aa_order):
            for j, aa2 in enumerate(aa_order):
                if aa1 in BLOSUM62 and aa2 in BLOSUM62[aa1]:
                    blosum_array[i, j] = BLOSUM62[aa1][aa2]
        
        # Save as binary
        with open(os.path.join(self.output_dir, 'blosum62.bin'), 'wb') as f:
            f.write(struct.pack('I', matrix_size))
            f.write(aa_order.encode('ascii'))
            blosum_array.tofile(f)
        
        print("  BLOSUM62 matrix written")
    
    def write_metadata(self):
        """Write metadata JSON file"""
        metadata = {
            'creation_date': datetime.now().isoformat(),
            'statistics': self.stats,
            'species_map': self.species_map,
            'gene_map': self.gene_map,
            'parameters': {
                'k': self.kmer_length,
                'source_dir': self.fq_genes_dir
            }
        }
        
        with open(os.path.join(self.output_dir, 'metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Also write accession map separately (can be large)
        with open(os.path.join(self.output_dir, 'accession_map.json'), 'w') as f:
            json.dump(self.accession_map, f, indent=2)
    
    def print_statistics(self):
        """Print database statistics"""
        print("\n=== Database Statistics ===")
        print(f"Species: {self.stats['species_count']}")
        print(f"Genes: {self.stats['gene_count']}")
        print(f"Total sequences: {self.stats['total_sequences']}")
        print(f"Unique proteins: {self.stats['unique_proteins']}")
        print(f"Total k-mers: {self.stats['total_kmers']}")
        print(f"Unique k-mers: {len(self.protein_kmers)}")
        
        # Gene distribution
        gene_counts = defaultdict(int)
        for protein in self.proteins:
            gene_counts[protein['gene_name']] += 1
        
        print("\nTop genes by sequence count:")
        for gene, count in sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"  {gene}: {count}")
        
        # Size information
        proteins_size = os.path.getsize(os.path.join(self.output_dir, 'proteins.bin')) / (1024*1024)
        kmers_size = os.path.getsize(os.path.join(self.output_dir, 'protein_kmers.bin')) / (1024*1024)
        print(f"\nDatabase sizes:")
        print(f"  Proteins: {proteins_size:.2f} MB")
        print(f"  K-mer index: {kmers_size:.2f} MB")
        print(f"  Total: {proteins_size + kmers_size:.2f} MB")


def main():
    parser = argparse.ArgumentParser(description='Build protein resistance database for GPU processing')
    parser.add_argument('fq_genes_dir', help='Directory containing species folders with gene JSON files')
    parser.add_argument('mutations_csv', help='CSV file with known resistance mutations')
    parser.add_argument('output_dir', help='Output directory for binary database')
    parser.add_argument('--kmer-length', type=int, default=3, help='K-mer length for protein index (default: 3)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isdir(args.fq_genes_dir):
        print(f"Error: {args.fq_genes_dir} is not a directory")
        sys.exit(1)
    
    # Build database
    builder = ProteinDatabaseBuilder(args.fq_genes_dir, args.mutations_csv, args.output_dir, args.kmer_length)
    builder.build()


if __name__ == '__main__':
    main()