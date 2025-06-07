#!/usr/bin/env python3
"""
build_clean_dynamic_database.py
Build database components based ONLY on actual data - no imputed drug/resistance info
"""

import json
import os
import sys
import struct
import numpy as np
from collections import defaultdict
from pathlib import Path
from datetime import datetime
import hashlib
import argparse
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CleanDatabaseBuilder:
    def __init__(self, wildtype_dir, output_dir, mutations_csv=None):
        self.wildtype_dir = Path(wildtype_dir)
        self.output_dir = Path(output_dir)
        self.mutations_csv = mutations_csv
        
        # Create output structure
        self.nucleotide_dir = self.output_dir / "nucleotide"
        self.protein_dir = self.output_dir / "protein"
        self.resistance_dir = self.output_dir / "resistance"
        
        for dir in [self.nucleotide_dir, self.protein_dir, self.resistance_dir]:
            dir.mkdir(exist_ok=True, parents=True)
        
        # Dynamic mappings - populated from data
        self.gene_to_id = {}
        self.species_to_id = {}
        self.id_to_gene = {}
        self.id_to_species = {}
        
        # Protein data
        self.proteins = []
        self.protein_kmers = defaultdict(list)
        self.kmer_length = 5
        
        # Mutations from CSV (if provided)
        self.known_mutations = []
        
    def discover_sequences(self):
        """Discover all species and genes from FASTA files"""
        logger.info(f"Discovering sequences from {self.wildtype_dir}")
        
        discovered_species = set()
        discovered_genes = set()
        
        fasta_files = list(self.wildtype_dir.glob("*.fasta")) + list(self.wildtype_dir.glob("*.fa"))
        
        for fasta_file in fasta_files:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.strip().lstrip('>')
                        species, gene = self.parse_header(header)
                        
                        if species and gene:
                            discovered_species.add(species)
                            discovered_genes.add(gene)
        
        # Assign IDs based on sorted order for consistency
        for i, gene in enumerate(sorted(discovered_genes)):
            self.gene_to_id[gene] = i
            self.id_to_gene[i] = gene
            
        for i, species in enumerate(sorted(discovered_species)):
            self.species_to_id[species] = i
            self.id_to_species[i] = species
            
        logger.info(f"Discovered {len(discovered_genes)} genes: {sorted(discovered_genes)}")
        logger.info(f"Discovered {len(discovered_species)} species: {sorted(discovered_species)}")
        
        return True
        
    def parse_header(self, header):
        """Parse FASTA header to extract species and gene"""
        # Remove '>' if present
        header = header.strip().lstrip('>')
        
        # Try multiple parsing strategies
        gene = None
        species = None
        
        # Strategy 1: Look for gene names we know about
        header_lower = header.lower()
        for candidate_gene in ['gyrA', 'gyrB', 'parC', 'parE', 'grlA', 'grlB', 'acrA', 'acrB', 'tolC']:
            if candidate_gene.lower() in header_lower:
                gene = candidate_gene
                break
        
        # Strategy 2: Parse structured headers
        if '_' in header:
            parts = header.split('_')
            
            # Try to find gene in last few parts
            if not gene:
                for part in reversed(parts[-3:]):
                    if len(part) >= 3 and len(part) <= 6 and part.isalpha():
                        gene = part
                        break
            
            # Species is typically first two parts
            if len(parts) >= 2:
                species = f"{parts[0]}_{parts[1]}"
        
        # Strategy 3: Look for species in brackets
        if '[' in header and ']' in header:
            start = header.find('[')
            end = header.find(']', start)
            species_text = header[start+1:end]
            species = species_text.replace(' ', '_')
            
        return species, gene
        
    def load_mutations_from_csv(self):
        """Load mutations from CSV if provided - NO drug/resistance info"""
        if not self.mutations_csv or not Path(self.mutations_csv).exists():
            logger.info("No mutations CSV provided or file not found")
            return
            
        logger.info(f"Loading mutations from {self.mutations_csv}")
        
        try:
            df = pd.read_csv(self.mutations_csv)
            
            # We only care about: Gene, Position, Wildtype AA, and Mutant AA
            required_cols = ['Gene', 'Position']
            
            for col in required_cols:
                if col not in df.columns:
                    logger.warning(f"Missing required column: {col}")
                    return
                    
            for _, row in df.iterrows():
                gene = str(row.get('Gene', '')).strip()
                
                # Skip if gene not in our database
                if gene not in self.gene_to_id:
                    continue
                    
                try:
                    position = int(row.get('Position', 0))
                    if position <= 0:
                        continue
                        
                    # Extract wildtype and mutant amino acids
                    wildtype = None
                    mutants = []
                    
                    # Try different column names
                    if 'Wildtype' in row:
                        wildtype = str(row['Wildtype']).strip()
                    elif 'Wild-type Amino Acid' in row:
                        wildtype = str(row['Wild-type Amino Acid']).strip()
                    elif 'WT' in row:
                        wildtype = str(row['WT']).strip()
                        
                    if 'Mutant' in row:
                        mutant_str = str(row['Mutant']).strip()
                        if ',' in mutant_str:
                            mutants = [m.strip() for m in mutant_str.split(',')]
                        else:
                            mutants = [mutant_str]
                    elif 'Mutant Amino Acid' in row:
                        mutants = [str(row['Mutant Amino Acid']).strip()]
                        
                    # Also try to parse from mutation notation (e.g., "S83L")
                    if 'Mutation' in row:
                        mutation_str = str(row['Mutation']).strip()
                        if len(mutation_str) >= 3:
                            wt = mutation_str[0]
                            pos_str = ''
                            mut = mutation_str[-1]
                            
                            for char in mutation_str[1:-1]:
                                if char.isdigit():
                                    pos_str += char
                                    
                            if not wildtype and wt.isalpha():
                                wildtype = wt
                            if not mutants and mut.isalpha():
                                mutants = [mut]
                                
                    # Only add if we have the minimum required info
                    if wildtype and mutants and position > 0:
                        self.known_mutations.append({
                            'gene': gene,
                            'gene_id': self.gene_to_id[gene],
                            'position': position,
                            'wildtype_aa': wildtype,
                            'mutant_aas': mutants
                        })
                        
                except (ValueError, TypeError) as e:
                    logger.warning(f"Error parsing mutation row: {e}")
                    continue
                    
            logger.info(f"Loaded {len(self.known_mutations)} mutations")
            
        except Exception as e:
            logger.error(f"Error loading mutations CSV: {e}")
            
    def load_protein_sequences(self):
        """Load all protein sequences"""
        logger.info("Loading protein sequences...")
        
        protein_id = 0
        
        for fasta_file in sorted(self.wildtype_dir.glob("*.fasta")) + sorted(self.wildtype_dir.glob("*.fa")):
            with open(fasta_file, 'r') as f:
                current_header = None
                current_sequence = ""
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Process previous sequence
                        if current_header and current_sequence:
                            self.process_sequence(current_header, current_sequence, protein_id)
                            protein_id += 1
                            
                        current_header = line.lstrip('>')
                        current_sequence = ""
                    else:
                        current_sequence += line
                        
                # Process last sequence
                if current_header and current_sequence:
                    self.process_sequence(current_header, current_sequence, protein_id)
                    protein_id += 1
                    
        logger.info(f"Loaded {len(self.proteins)} protein sequences")
        
    def process_sequence(self, header, sequence, protein_id):
        """Process a single protein sequence"""
        species, gene = self.parse_header(header)
        
        if not species or not gene:
            logger.warning(f"Could not parse header: {header}")
            return
            
        # Clean sequence
        sequence = sequence.upper().replace('*', '').replace(' ', '').replace('\n', '')
        
        if len(sequence) < 20 or 'X' in sequence:
            return
            
        protein_entry = {
            'id': protein_id,
            'species_id': self.species_to_id.get(species, -1),
            'gene_id': self.gene_to_id.get(gene, -1),
            'sequence': sequence,
            'length': len(sequence),
            'gene_name': gene,
            'species_name': species,
            'header': header
        }
        
        self.proteins.append(protein_entry)
        
    def build_kmer_index(self):
        """Build k-mer index for proteins"""
        logger.info(f"Building {self.kmer_length}-mer index...")
        
        for protein_idx, protein in enumerate(self.proteins):
            seq = protein['sequence']
            
            for i in range(len(seq) - self.kmer_length + 1):
                kmer = seq[i:i + self.kmer_length]
                self.protein_kmers[kmer].append((protein_idx, i))
                
        logger.info(f"Indexed {len(self.protein_kmers)} unique {self.kmer_length}-mers")
        
    def write_protein_database(self):
        """Write protein database in format expected by GPU code"""
        logger.info("Writing protein database...")
        
        # Write proteins.bin
        with open(self.protein_dir / 'proteins.bin', 'wb') as f:
            # Header: number of proteins
            f.write(struct.pack('I', len(self.proteins)))
            
            # Write concatenated sequences
            for protein in self.proteins:
                f.write(protein['sequence'].encode('ascii'))
                
        # Write protein_kmers.bin
        with open(self.protein_dir / 'protein_kmers.bin', 'wb') as f:
            f.write(struct.pack('I', self.kmer_length))
            f.write(struct.pack('I', len(self.protein_kmers)))
            
            for kmer, positions in sorted(self.protein_kmers.items()):
                f.write(kmer.encode('ascii'))
                f.write(struct.pack('I', len(positions)))
                for protein_idx, pos in positions:
                    f.write(struct.pack('II', protein_idx, pos))
                    
        # Write metadata.json
        metadata = {
            'creation_date': datetime.now().isoformat(),
            'database_type': 'protein_wildtype',
            'num_proteins': len(self.proteins),
            'num_species': len(self.species_to_id),
            'num_genes': len(self.gene_to_id),
            'species_map': self.id_to_species,
            'gene_map': self.id_to_gene,
            'kmer_length': self.kmer_length,
            'source_directory': str(self.wildtype_dir)
        }
        
        with open(self.protein_dir / 'metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
            
        # Write detailed protein information
        protein_details = []
        for p in self.proteins:
            protein_details.append({
                'id': p['id'],
                'species': p['species_name'],
                'gene': p['gene_name'],
                'species_id': p['species_id'],
                'gene_id': p['gene_id'],
                'length': p['length'],
                'header': p['header']
            })
            
        with open(self.protein_dir / 'protein_details.json', 'w') as f:
            json.dump(protein_details, f, indent=2)
            
    def write_resistance_database(self):
        """Write minimal resistance database structure"""
        logger.info("Writing resistance database...")
        
        # Create minimal resistance_db.json
        resistance_db = {
            'metadata': {
                'version': '3.0',
                'description': 'Minimal resistance database - positions only',
                'creation_date': datetime.now().isoformat(),
                'note': 'No drug or resistance level information included'
            },
            'gene_map': self.id_to_gene,
            'species_map': self.id_to_species,
            'mutations': self.known_mutations
        }
        
        with open(self.resistance_dir / 'resistance_db.json', 'w') as f:
            json.dump(resistance_db, f, indent=2)
            
        # Also write a simple mutations list
        if self.known_mutations:
            with open(self.resistance_dir / 'mutations.tsv', 'w') as f:
                f.write("gene\tgene_id\tposition\twildtype\tmutants\n")
                for mut in self.known_mutations:
                    mutants_str = ','.join(mut['mutant_aas'])
                    f.write(f"{mut['gene']}\t{mut['gene_id']}\t{mut['position']}\t"
                           f"{mut['wildtype_aa']}\t{mutants_str}\n")
                           
    def write_mappings(self):
        """Write mapping files for reference"""
        logger.info("Writing mapping files...")
        
        # Combined mappings file
        mappings = {
            'gene_to_id': self.gene_to_id,
            'species_to_id': self.species_to_id,
            'id_to_gene': self.id_to_gene,
            'id_to_species': self.id_to_species
        }
        
        with open(self.output_dir / 'database_mappings.json', 'w') as f:
            json.dump(mappings, f, indent=2)
            
        # TSV files for easy viewing
        with open(self.output_dir / 'gene_mappings.tsv', 'w') as f:
            f.write("gene_name\tgene_id\n")
            for gene, gid in sorted(self.gene_to_id.items()):
                f.write(f"{gene}\t{gid}\n")
                
        with open(self.output_dir / 'species_mappings.tsv', 'w') as f:
            f.write("species_name\tspecies_id\n")
            for species, sid in sorted(self.species_to_id.items()):
                f.write(f"{species}\t{sid}\n")
                
    def build_all(self):
        """Build all database components"""
        logger.info("Building clean dynamic database...")
        
        # Step 1: Discover sequences
        if not self.discover_sequences():
            return False
            
        # Step 2: Load mutations if CSV provided
        if self.mutations_csv:
            self.load_mutations_from_csv()
            
        # Step 3: Load protein sequences
        self.load_protein_sequences()
        
        # Step 4: Build k-mer index
        self.build_kmer_index()
        
        # Step 5: Write databases
        self.write_protein_database()
        self.write_resistance_database()
        self.write_mappings()
        
        # Print summary
        self.print_summary()
        
        return True
        
    def print_summary(self):
        """Print build summary"""
        print("\n=== Clean Database Build Summary ===")
        print(f"Source directory: {self.wildtype_dir}")
        print(f"Output directory: {self.output_dir}")
        
        print(f"\nDiscovered {len(self.species_to_id)} species:")
        for i, species in sorted(self.id_to_species.items()):
            count = sum(1 for p in self.proteins if p['species_id'] == i)
            print(f"  ID {i}: {species} ({count} proteins)")
            
        print(f"\nDiscovered {len(self.gene_to_id)} genes:")
        for i, gene in sorted(self.id_to_gene.items()):
            count = sum(1 for p in self.proteins if p['gene_id'] == i)
            print(f"  ID {i}: {gene} ({count} proteins)")
            
        print(f"\nDatabase statistics:")
        print(f"  Total proteins: {len(self.proteins)}")
        print(f"  Unique {self.kmer_length}-mers: {len(self.protein_kmers)}")
        print(f"  Known mutations: {len(self.known_mutations)}")
        
        print(f"\nOutput files created:")
        print(f"  protein/proteins.bin")
        print(f"  protein/protein_kmers.bin")
        print(f"  protein/metadata.json")
        print(f"  protein/protein_details.json")
        print(f"  resistance/resistance_db.json")
        if self.known_mutations:
            print(f"  resistance/mutations.tsv")
        print(f"  database_mappings.json")
        print(f"  gene_mappings.tsv")
        print(f"  species_mappings.tsv")

def main():
    parser = argparse.ArgumentParser(description='Build clean dynamic resistance database')
    parser.add_argument('wildtype_dir', help='Directory containing wildtype protein FASTA files')
    parser.add_argument('output_dir', help='Output directory for database')
    parser.add_argument('--mutations-csv', help='Optional CSV file with known mutations (Gene, Position, etc.)')
    parser.add_argument('--kmer-length', type=int, default=5, help='K-mer length for protein index (default: 5)')
    
    args = parser.parse_args()
    
    # Create builder
    builder = CleanDatabaseBuilder(
        args.wildtype_dir,
        args.output_dir,
        args.mutations_csv
    )
    
    # Override k-mer length if specified
    if args.kmer_length:
        builder.kmer_length = args.kmer_length
    
    # Build database
    if builder.build_all():
        print("\nDatabase build completed successfully!")
        print(f"\nTo use this database with the pipeline:")
        print(f"  ./integrated_resistance_pipeline \\")
        print(f"    {args.output_dir}/nucleotide \\")
        print(f"    {args.output_dir}/protein \\")
        print(f"    {args.output_dir}/resistance \\")
        print(f"    reads_R1.fastq.gz reads_R2.fastq.gz output_prefix")
    else:
        print("\nDatabase build failed!")
        sys.exit(1)

if __name__ == '__main__':
    main()