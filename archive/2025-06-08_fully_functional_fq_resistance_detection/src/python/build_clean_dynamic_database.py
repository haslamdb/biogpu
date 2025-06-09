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
    def __init__(self, wildtype_dir, output_dir, mutations_csv=None, nucleotide_dir=None):
        self.wildtype_dir = Path(wildtype_dir)
        self.output_dir = Path(output_dir)
        self.mutations_csv = mutations_csv
        self.nucleotide_sequences_dir = Path(nucleotide_dir) if nucleotide_dir else None
        
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
        self.kmer_length = 8
        
        # Nucleotide k-mer data
        self.nucleotide_kmer_length = 15
        self.include_rc = True  # Include reverse complement k-mers
        
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
    
    def reverse_complement(self, seq):
        """Get reverse complement of sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
        seq = seq.upper().replace('U', 'T')  # Convert RNA to DNA
        return ''.join(complement.get(base, 'N') for base in seq[::-1])
    
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
    
    def generate_kmers_with_rc(self, sequence, k):
        """Generate k-mers from both forward and reverse complement strands"""
        sequence = sequence.upper().replace('U', 'T')  # Convert RNA to DNA
        
        if len(sequence) < k:
            return [], []
        
        valid_bases = set('ATCG')
        forward_kmers = []
        rc_kmers = []
        
        # Generate forward k-mers
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_bases = set(kmer)
            
            if kmer_bases.issubset(valid_bases):
                forward_kmers.append((kmer, i, False))  # (kmer, position, is_rc)
        
        # Generate reverse complement k-mers
        rc_sequence = self.reverse_complement(sequence)
        rc_length = len(rc_sequence)
        
        for i in range(len(rc_sequence) - k + 1):
            kmer = rc_sequence[i:i+k]
            kmer_bases = set(kmer)
            
            if kmer_bases.issubset(valid_bases):
                # Calculate original position (from end of sequence)
                original_pos = rc_length - i - k
                rc_kmers.append((kmer, original_pos, True))  # (kmer, position, is_rc)
        
        return forward_kmers, rc_kmers
        
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
            
            # We only care about: gene, location (position), wt, and mut
            # Check for both uppercase and lowercase column names
            gene_col = 'gene' if 'gene' in df.columns else 'Gene'
            pos_col = 'location' if 'location' in df.columns else 'Position'
            
            if gene_col not in df.columns or pos_col not in df.columns:
                logger.warning(f"Missing required columns. Found columns: {list(df.columns)}")
                return
                    
            for _, row in df.iterrows():
                gene = str(row.get(gene_col, '')).strip()
                
                # Skip if gene not in our database
                if gene not in self.gene_to_id:
                    continue
                    
                try:
                    position = int(row.get(pos_col, 0))
                    if position <= 0:
                        continue
                        
                    # Extract wildtype and mutant amino acids
                    wildtype = None
                    mutants = []
                    
                    # Try different column names
                    if 'wt' in row:
                        wildtype = str(row['wt']).strip()
                    elif 'Wildtype' in row:
                        wildtype = str(row['Wildtype']).strip()
                    elif 'Wild-type Amino Acid' in row:
                        wildtype = str(row['Wild-type Amino Acid']).strip()
                    elif 'WT' in row:
                        wildtype = str(row['WT']).strip()
                        
                    if 'mut' in row:
                        mutant_str = str(row['mut']).strip()
                        if ',' in mutant_str:
                            mutants = [m.strip() for m in mutant_str.split(',')]
                        else:
                            mutants = [mutant_str]
                    elif 'Mutant' in row:
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
    
    def build_nucleotide_index(self):
        """Build nucleotide k-mer index from sequence files"""
        if not self.nucleotide_sequences_dir or not self.nucleotide_sequences_dir.exists():
            logger.warning("No nucleotide sequences directory provided or directory doesn't exist")
            return False
            
        logger.info(f"Building nucleotide k-mer index from {self.nucleotide_sequences_dir}")
        
        # Data structures for k-mer index
        kmer_entries = []
        all_sequences = []
        sequence_id = 0
        
        # Create separate ID mappings for nucleotide sequences
        # (to include all species/genes, not just those in protein DB)
        nucleotide_species_to_id = {}
        nucleotide_gene_to_id = {}
        next_species_id = 0
        next_gene_id = 0
        
        # Process each species directory
        for species_dir in sorted(self.nucleotide_sequences_dir.iterdir()):
            if not species_dir.is_dir():
                continue
                
            species_name = species_dir.name
            
            # Assign species ID (create new if not in protein DB)
            if species_name not in nucleotide_species_to_id:
                if species_name in self.species_to_id:
                    nucleotide_species_to_id[species_name] = self.species_to_id[species_name]
                else:
                    nucleotide_species_to_id[species_name] = len(self.species_to_id) + next_species_id
                    next_species_id += 1
                    
            species_id = nucleotide_species_to_id[species_name]
            logger.info(f"Processing nucleotide sequences for {species_name} (ID: {species_id})")
            
            # Process JSON files in species directory
            for json_file in sorted(species_dir.glob("*.json")):
                gene_name = json_file.stem
                
                # Assign gene ID (create new if not in protein DB)
                if gene_name not in nucleotide_gene_to_id:
                    if gene_name in self.gene_to_id:
                        nucleotide_gene_to_id[gene_name] = self.gene_to_id[gene_name]
                    else:
                        nucleotide_gene_to_id[gene_name] = len(self.gene_to_id) + next_gene_id
                        next_gene_id += 1
                        
                gene_id = nucleotide_gene_to_id[gene_name]
                
                # Load sequence data
                with open(json_file, 'r') as f:
                    sequences = json.load(f)
                
                if not isinstance(sequences, list) or len(sequences) == 0:
                    logger.warning(f"No sequences in {json_file}")
                    continue
                
                # Process each sequence in the file
                for seq_idx, seq_data in enumerate(sequences):
                    if 'gene_features' not in seq_data or not seq_data['gene_features']:
                        continue
                        
                    # Get the first gene feature (usually there's only one)
                    gene_feature = seq_data['gene_features'][0]
                    
                    if 'nucleotide_sequence' not in gene_feature:
                        logger.warning(f"No nucleotide sequence in {json_file} entry {seq_idx}")
                        continue
                        
                    nucleotide_seq = gene_feature['nucleotide_sequence']
                    
                    # Store sequence
                    seq_entry = {
                        'sequence_id': sequence_id,
                        'species_id': species_id,
                        'gene_id': gene_id,
                        'species_name': species_name,
                        'gene_name': gene_name,
                        'sequence': nucleotide_seq,
                        'length': len(nucleotide_seq),
                        'accession': seq_data.get('accession', '')
                    }
                    all_sequences.append(seq_entry)
                    
                    # Generate k-mers
                    if self.include_rc:
                        forward_kmers, rc_kmers = self.generate_kmers_with_rc(nucleotide_seq, self.nucleotide_kmer_length)
                        all_kmers = forward_kmers + rc_kmers
                    else:
                        forward_kmers, rc_kmers = self.generate_kmers_with_rc(nucleotide_seq, self.nucleotide_kmer_length)
                        all_kmers = forward_kmers  # Only use forward k-mers
                    
                    # Add k-mer entries
                    for kmer_seq, position, is_rc in all_kmers:
                        encoded = self.encode_kmer(kmer_seq)
                        if encoded is not None:
                            kmer_entries.append({
                                'kmer': encoded,
                                'gene_id': gene_id,
                                'species_id': species_id,
                                'seq_id': sequence_id,
                                'position': position
                            })
                    
                    sequence_id += 1
                
        # Sort k-mer entries by encoded k-mer value for binary search
        kmer_entries.sort(key=lambda x: x['kmer'])
        
        # Write binary k-mer index
        kmer_index_path = self.nucleotide_dir / "kmer_index.bin"
        with open(kmer_index_path, 'wb') as f:
            # Header
            f.write(struct.pack('I', len(kmer_entries)))  # Number of entries
            f.write(struct.pack('I', self.nucleotide_kmer_length))  # K-mer length
            
            # K-mer entries
            for entry in kmer_entries:
                f.write(struct.pack('Q', entry['kmer']))      # 64-bit k-mer
                f.write(struct.pack('I', entry['gene_id']))   # Gene ID
                f.write(struct.pack('I', entry['species_id'])) # Species ID  
                f.write(struct.pack('I', entry['seq_id']))    # Sequence ID
                f.write(struct.pack('H', entry['position']))  # Position (16-bit)
        
        logger.info(f"Wrote {len(kmer_entries)} k-mer entries to {kmer_index_path}")
        
        # Write sequence database
        seq_db_path = self.nucleotide_dir / "sequences.bin"
        with open(seq_db_path, 'wb') as f:
            # Header
            f.write(struct.pack('I', len(all_sequences)))
            
            # Sequences
            for seq in all_sequences:
                seq_bytes = seq['sequence'].encode('utf-8')
                f.write(struct.pack('I', len(seq_bytes)))
                f.write(seq_bytes)
                f.write(struct.pack('I', seq['length']))
                f.write(struct.pack('I', seq['species_id']))
                f.write(struct.pack('I', seq['gene_id']))
                
                # Accession
                acc_bytes = seq['accession'].encode('utf-8')
                f.write(struct.pack('I', len(acc_bytes)))
                f.write(acc_bytes)
        
        logger.info(f"Wrote {len(all_sequences)} sequences to {seq_db_path}")
        
        # Write metadata
        metadata = {
            'kmer_length': self.nucleotide_kmer_length,
            'num_kmers': len(kmer_entries),
            'num_sequences': len(all_sequences),
            'include_rc': self.include_rc,
            'species_included': list(set(seq['species_name'] for seq in all_sequences)),
            'genes_included': list(set(seq['gene_name'] for seq in all_sequences)),
            'species_map': {str(v): k for k, v in nucleotide_species_to_id.items()},
            'gene_map': {str(v): k for k, v in nucleotide_gene_to_id.items()},
            'num_species': len(nucleotide_species_to_id),
            'num_genes': len(nucleotide_gene_to_id)
        }
        
        metadata_path = self.nucleotide_dir / "metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Also write nucleotide-specific ID mappings
        nucleotide_mappings = {
            'species_to_id': nucleotide_species_to_id,
            'gene_to_id': nucleotide_gene_to_id,
            'id_to_species': {v: k for k, v in nucleotide_species_to_id.items()},
            'id_to_gene': {v: k for k, v in nucleotide_gene_to_id.items()}
        }
        
        with open(self.nucleotide_dir / 'id_mappings.json', 'w') as f:
            json.dump(nucleotide_mappings, f, indent=2)
            
        logger.info(f"Nucleotide k-mer index built successfully")
        logger.info(f"  Total species: {len(nucleotide_species_to_id)} ({len(nucleotide_species_to_id) - len(self.species_to_id)} new)")
        logger.info(f"  Total genes: {len(nucleotide_gene_to_id)} ({len(nucleotide_gene_to_id) - len(self.gene_to_id)} new)")
        return True
                
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
        
        # Step 6: Build nucleotide index if sequences provided
        if self.nucleotide_sequences_dir:
            self.build_nucleotide_index()
        
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
        
        if self.nucleotide_sequences_dir and (self.nucleotide_dir / "kmer_index.bin").exists():
            print(f"\nNucleotide index files:")
            print(f"  nucleotide/kmer_index.bin")
            print(f"  nucleotide/sequences.bin")
            print(f"  nucleotide/metadata.json")
            print(f"  Nucleotide k-mer length: {self.nucleotide_kmer_length}")
            print(f"  Reverse complement k-mers: {'ENABLED' if self.include_rc else 'DISABLED'}")

def main():
    parser = argparse.ArgumentParser(description='Build clean dynamic resistance database')
    parser.add_argument('wildtype_dir', nargs='?', default='data/wildtype_protein_seqs', 
                        help='Directory containing wildtype protein FASTA files (default: data/wildtype_protein_seqs)')
    parser.add_argument('output_dir', nargs='?', default='data/integrated_clean_db',
                        help='Output directory for database (default: data/integrated_clean_db)')
    parser.add_argument('--mutations-csv', help='Optional CSV file with known mutations (Gene, Position, etc.)')
    parser.add_argument('--kmer-length', type=int, default=8, help='K-mer length for protein index (default: 8)')
    parser.add_argument('--nucleotide-kmer-length', type=int, default=15, help='K-mer length for nucleotide index (default: 15)')
    parser.add_argument('--nucleotide-sequences', default='data/fq_genes', 
                        help='Directory containing nucleotide sequences (default: data/fq_genes)')
    parser.add_argument('--no-rc', action='store_true', help='Disable reverse complement k-mers for nucleotide index')
    
    args = parser.parse_args()
    
    # Create builder
    builder = CleanDatabaseBuilder(
        args.wildtype_dir,
        args.output_dir,
        args.mutations_csv,
        args.nucleotide_sequences
    )
    
    # Set k-mer lengths from command line
    builder.kmer_length = args.kmer_length
    builder.nucleotide_kmer_length = args.nucleotide_kmer_length
    builder.include_rc = not args.no_rc
    
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