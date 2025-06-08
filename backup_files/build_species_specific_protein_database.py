#!/usr/bin/env python3
"""
build_species_specific_protein_database.py
Build GPU-ready protein database preserving species-specific sequences
WITHOUT imputing any resistance or drug information
"""

import json
import struct
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import hashlib

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ProteinSequence:
    """Represents a specific protein sequence from a specific species"""
    species_name: str
    gene_name: str
    sequence: str
    species_gene_id: int  # Unique ID for species-gene combination
    source_file: str  # Track where this came from
    
@dataclass 
class KnownMutation:
    """Known resistance mutation - only position and amino acids, no drug/MIC guessing"""
    species_gene_id: int
    position: int  # 1-based
    wildtype_aa: str
    mutant_aas: List[str]
    source: str  # Where this mutation info came from

class SpeciesSpecificProteinDatabase:
    """Build protein database maintaining species-specific sequences"""
    
    def __init__(self):
        self.protein_sequences: List[ProteinSequence] = []
        self.known_mutations: List[KnownMutation] = []
        
        # Mapping dictionaries
        self.species_gene_to_id: Dict[Tuple[str, str], int] = {}
        self.id_to_species_gene: Dict[int, Tuple[str, str]] = {}
        self.next_id = 0
        
        # Track what we've seen
        self.all_species: Set[str] = set()
        self.all_genes: Set[str] = set()
        self.sequence_hashes: Dict[str, Tuple[str, str]] = {}  # To detect duplicates
        
    def get_or_create_species_gene_id(self, species: str, gene: str) -> int:
        """Get ID for species-gene combination, creating if needed"""
        key = (species, gene)
        if key not in self.species_gene_to_id:
            self.species_gene_to_id[key] = self.next_id
            self.id_to_species_gene[self.next_id] = key
            self.next_id += 1
            
        self.all_species.add(species)
        self.all_genes.add(gene)
        
        return self.species_gene_to_id[key]
    
    def load_sequences_from_fasta(self, fasta_path: str):
        """Load protein sequences from a FASTA file"""
        fasta_file = Path(fasta_path)
        if not fasta_file.exists():
            logger.error(f"File not found: {fasta_path}")
            return
            
        logger.info(f"Loading sequences from {fasta_file.name}")
        sequences_loaded = 0
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Parse header to extract species and gene
            # Expected format: Species_name_gene or similar
            header = record.id
            
            # Try to intelligently parse the header
            parts = header.split('_')
            
            # Handle special cases
            if 'Pseudomonas_aeruginosa' in header:
                # Handle A/B variants
                if '_A_' in header:
                    species = 'Pseudomonas_aeruginosa_A'
                    gene = parts[-1]
                elif '_B_' in header:
                    species = 'Pseudomonas_aeruginosa_B'
                    gene = parts[-1]
                else:
                    species = 'Pseudomonas_aeruginosa'
                    gene = parts[-1]
            else:
                # Standard format: everything except last part is species
                gene = parts[-1]
                species = '_'.join(parts[:-1])
            
            # Validate we got reasonable values
            if not species or not gene:
                logger.warning(f"Could not parse header: {header}")
                continue
                
            # Get the actual sequence
            sequence = str(record.seq).upper()
            
            # Check for duplicates using sequence hash
            seq_hash = hashlib.md5(sequence.encode()).hexdigest()
            if seq_hash in self.sequence_hashes:
                existing_species, existing_gene = self.sequence_hashes[seq_hash]
                logger.info(f"Duplicate sequence found: {species}_{gene} matches {existing_species}_{existing_gene}")
                # Still add it as it's a different species-gene combination
            else:
                self.sequence_hashes[seq_hash] = (species, gene)
            
            # Create ID for this species-gene combination
            species_gene_id = self.get_or_create_species_gene_id(species, gene)
            
            # Store the protein sequence
            protein = ProteinSequence(
                species_name=species,
                gene_name=gene,
                sequence=sequence,
                species_gene_id=species_gene_id,
                source_file=fasta_file.name
            )
            
            self.protein_sequences.append(protein)
            sequences_loaded += 1
            
            logger.debug(f"Loaded {species} {gene}: {len(sequence)} aa")
            
        logger.info(f"Loaded {sequences_loaded} sequences from {fasta_file.name}")
    
    def load_sequences_from_directory(self, directory: str):
        """Load all FASTA files from a directory"""
        dir_path = Path(directory)
        if not dir_path.exists():
            logger.error(f"Directory not found: {directory}")
            return
            
        fasta_files = list(dir_path.glob("*.fasta")) + list(dir_path.glob("*.fa"))
        logger.info(f"Found {len(fasta_files)} FASTA files in {directory}")
        
        for fasta_file in sorted(fasta_files):
            self.load_sequences_from_fasta(str(fasta_file))
            
    def load_known_mutations_csv(self, csv_path: str):
        """Load known mutations from CSV - only positions, no drug/MIC guessing"""
        import pandas as pd
        
        logger.info(f"Loading mutations from {csv_path}")
        df = pd.read_csv(csv_path)
        
        mutations_loaded = 0
        
        for _, row in df.iterrows():
            # Required columns: Species, Gene, Position, Wildtype, Mutant
            if not all(col in row for col in ['Species', 'Gene', 'Position', 'Wildtype', 'Mutant']):
                logger.warning(f"Missing required columns in row: {row}")
                continue
                
            species = row['Species']
            gene = row['Gene']
            position = int(row['Position'])
            wildtype = row['Wildtype']
            mutants = row['Mutant'].split(',') if ',' in str(row['Mutant']) else [row['Mutant']]
            
            # Get species-gene ID
            if (species, gene) not in self.species_gene_to_id:
                logger.warning(f"No sequence found for {species} {gene}")
                continue
                
            species_gene_id = self.species_gene_to_id[(species, gene)]
            
            # Check if we already have this mutation position
            existing = None
            for mut in self.known_mutations:
                if mut.species_gene_id == species_gene_id and mut.position == position:
                    existing = mut
                    break
                    
            if existing:
                # Add new mutants to existing position
                for mutant in mutants:
                    if mutant not in existing.mutant_aas:
                        existing.mutant_aas.append(mutant)
            else:
                # Create new mutation entry
                mutation = KnownMutation(
                    species_gene_id=species_gene_id,
                    position=position,
                    wildtype_aa=wildtype,
                    mutant_aas=mutants,
                    source=csv_path
                )
                self.known_mutations.append(mutation)
                mutations_loaded += 1
                
        logger.info(f"Loaded {mutations_loaded} mutation positions from {csv_path}")
    
    def validate_mutations(self):
        """Validate that mutation positions match actual sequences"""
        validated = 0
        errors = 0
        
        for mutation in self.known_mutations:
            # Find the corresponding sequence
            seq_found = False
            for protein in self.protein_sequences:
                if protein.species_gene_id == mutation.species_gene_id:
                    seq_found = True
                    
                    # Check position validity
                    if mutation.position > len(protein.sequence):
                        logger.error(f"{protein.species_name} {protein.gene_name}: "
                                   f"position {mutation.position} exceeds sequence length {len(protein.sequence)}")
                        errors += 1
                        break
                        
                    # Check wildtype matches
                    actual_aa = protein.sequence[mutation.position - 1]
                    if actual_aa != mutation.wildtype_aa:
                        logger.warning(f"{protein.species_name} {protein.gene_name} position {mutation.position}: "
                                     f"expected {mutation.wildtype_aa}, found {actual_aa}")
                        errors += 1
                    else:
                        validated += 1
                    break
                    
            if not seq_found:
                species, gene = self.id_to_species_gene[mutation.species_gene_id]
                logger.error(f"No sequence found for mutation in {species} {gene}")
                errors += 1
                
        logger.info(f"Validated {validated} mutations, found {errors} errors")
        
    def export_protein_fasta(self, output_dir: str):
        """Export all protein sequences as FASTA files"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        # Write all wildtype sequences
        all_proteins = []
        for protein in self.protein_sequences:
            record = SeqRecord(
                Seq(protein.sequence),
                id=f"{protein.species_name}_{protein.gene_name}",
                description=f"Species-gene ID: {protein.species_gene_id}"
            )
            all_proteins.append(record)
            
        output_file = output_path / "wildtype_proteins.fasta"
        SeqIO.write(all_proteins, output_file, "fasta")
        logger.info(f"Wrote {len(all_proteins)} protein sequences to {output_file}")
        
        # Also write a file organized by gene
        for gene in sorted(self.all_genes):
            gene_proteins = []
            for protein in self.protein_sequences:
                if protein.gene_name == gene:
                    record = SeqRecord(
                        Seq(protein.sequence),
                        id=f"{protein.species_name}_{protein.gene_name}",
                        description=f"Species-gene ID: {protein.species_gene_id}"
                    )
                    gene_proteins.append(record)
                    
            if gene_proteins:
                gene_file = output_path / f"{gene}_all_species.fasta"
                SeqIO.write(gene_proteins, gene_file, "fasta")
                logger.info(f"Wrote {len(gene_proteins)} {gene} sequences to {gene_file}")
                
    def export_binary_database(self, output_file: str):
        """Export database in binary format for GPU"""
        with open(output_file, 'wb') as f:
            # Write header
            f.write(struct.pack('I', len(self.protein_sequences)))  # num_sequences
            f.write(struct.pack('I', len(self.known_mutations)))   # num_mutations
            f.write(struct.pack('I', len(self.all_species)))       # num_species
            f.write(struct.pack('I', len(self.all_genes)))         # num_genes
            
            # Write protein sequences
            for protein in self.protein_sequences:
                f.write(struct.pack('I', protein.species_gene_id))
                f.write(struct.pack('I', len(protein.sequence)))
                f.write(protein.sequence.encode('ascii'))
                
            # Write mutations (without drug/MIC info)
            for mut in self.known_mutations:
                f.write(struct.pack('I', mut.species_gene_id))
                f.write(struct.pack('H', mut.position))
                f.write(mut.wildtype_aa.encode('ascii'))
                
                # Write mutant AAs (pad to 10)
                for i in range(10):
                    if i < len(mut.mutant_aas):
                        f.write(mut.mutant_aas[i].encode('ascii'))
                    else:
                        f.write(b'\0')
                        
                f.write(struct.pack('B', len(mut.mutant_aas)))
                
        logger.info(f"Exported binary database to {output_file}")
        
    def export_json_database(self, output_file: str):
        """Export database as JSON for debugging and pipeline use"""
        
        # Create mappings
        species_gene_map = {}
        for (species, gene), id_val in self.species_gene_to_id.items():
            species_gene_map[f"{species}_{gene}"] = id_val
            
        # Organize sequences by species-gene
        sequences_dict = {}
        for protein in self.protein_sequences:
            key = f"{protein.species_name}_{protein.gene_name}"
            sequences_dict[key] = {
                'species': protein.species_name,
                'gene': protein.gene_name,
                'species_gene_id': protein.species_gene_id,
                'sequence': protein.sequence,
                'length': len(protein.sequence),
                'source': protein.source_file
            }
            
        # Organize mutations
        mutations_list = []
        for mut in self.known_mutations:
            species, gene = self.id_to_species_gene[mut.species_gene_id]
            mutations_list.append({
                'species': species,
                'gene': gene,
                'species_gene_id': mut.species_gene_id,
                'position': mut.position,
                'wildtype': mut.wildtype_aa,
                'mutants': mut.mutant_aas,
                'source': mut.source
            })
            
        db_dict = {
            'metadata': {
                'version': '2.0',
                'description': 'Species-specific protein sequences without imputed data',
                'num_sequences': len(self.protein_sequences),
                'num_mutations': len(self.known_mutations),
                'num_species': len(self.all_species),
                'num_genes': len(self.all_genes),
                'species_list': sorted(list(self.all_species)),
                'gene_list': sorted(list(self.all_genes))
            },
            'species_gene_id_map': species_gene_map,
            'sequences': sequences_dict,
            'known_mutations': mutations_list,
            'statistics': {
                'sequences_per_gene': {},
                'sequences_per_species': {},
                'average_length_per_gene': {}
            }
        }
        
        # Calculate statistics
        for gene in self.all_genes:
            gene_seqs = [p for p in self.protein_sequences if p.gene_name == gene]
            db_dict['statistics']['sequences_per_gene'][gene] = len(gene_seqs)
            if gene_seqs:
                avg_len = sum(len(p.sequence) for p in gene_seqs) / len(gene_seqs)
                db_dict['statistics']['average_length_per_gene'][gene] = avg_len
                
        for species in self.all_species:
            species_seqs = [p for p in self.protein_sequences if p.species_name == species]
            db_dict['statistics']['sequences_per_species'][species] = len(species_seqs)
            
        with open(output_file, 'w') as f:
            json.dump(db_dict, f, indent=2)
            
        logger.info(f"Exported JSON database to {output_file}")
        
    def export_id_mappings(self, output_dir: str):
        """Export ID mappings in simple format for pipeline use"""
        output_path = Path(output_dir)
        
        # Species-gene to ID mapping
        with open(output_path / 'species_gene_to_id.tsv', 'w') as f:
            f.write("species\tgene\tspecies_gene_id\n")
            for (species, gene), id_val in sorted(self.species_gene_to_id.items()):
                f.write(f"{species}\t{gene}\t{id_val}\n")
                
        # ID to species-gene mapping  
        with open(output_path / 'id_to_species_gene.tsv', 'w') as f:
            f.write("species_gene_id\tspecies\tgene\n")
            for id_val, (species, gene) in sorted(self.id_to_species_gene.items()):
                f.write(f"{id_val}\t{species}\t{gene}\n")
                
        logger.info(f"Exported ID mappings to {output_path}")
        
    def print_summary(self):
        """Print database summary"""
        print("\n=== Database Summary ===")
        print(f"Total sequences: {len(self.protein_sequences)}")
        print(f"Total species: {len(self.all_species)}")
        print(f"Total genes: {len(self.all_genes)}")
        print(f"Known mutations: {len(self.known_mutations)}")
        
        print("\nSequences per gene:")
        for gene in sorted(self.all_genes):
            count = sum(1 for p in self.protein_sequences if p.gene_name == gene)
            print(f"  {gene}: {count} sequences")
            
        print("\nSequences per species:")
        for species in sorted(self.all_species):
            count = sum(1 for p in self.protein_sequences if p.species_name == species)
            print(f"  {species}: {count} sequences")
            
        print("\nExample sequences:")
        for protein in self.protein_sequences[:5]:
            print(f"  {protein.species_name} {protein.gene_name} (ID: {protein.species_gene_id}): {len(protein.sequence)} aa")
            
def main():
    parser = argparse.ArgumentParser(description='Build species-specific protein database')
    parser.add_argument('--input-dir', required=True, help='Directory with protein FASTA files')
    parser.add_argument('--mutations-csv', help='Optional CSV file with known mutations')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Build database
    db = SpeciesSpecificProteinDatabase()
    
    # Load sequences
    db.load_sequences_from_directory(args.input_dir)
    
    # Load mutations if provided
    if args.mutations_csv:
        db.load_known_mutations_csv(args.mutations_csv)
        db.validate_mutations()
        
    # Export everything
    db.export_protein_fasta(args.output_dir)
    db.export_binary_database(output_path / 'protein_database.bin')
    db.export_json_database(output_path / 'protein_database.json')
    db.export_id_mappings(args.output_dir)
    
    # Print summary
    db.print_summary()
    
if __name__ == '__main__':
    main()