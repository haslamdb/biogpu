#!/usr/bin/env python3
"""
build_resistance_database.py
Build GPU-ready resistance database from various sources
"""

import json
import struct
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ResistanceMutation:
    """Single resistance mutation entry"""
    gene_name: str
    gene_id: int
    position: int  # 1-based
    wildtype_aa: str
    resistant_aas: List[str]
    resistance_level: float  # 0-1
    drugs_affected: List[str]
    pubmed_ids: List[str] = None
    notes: str = ""

@dataclass
class ResistanceDatabase:
    """Complete resistance database"""
    mutations: List[ResistanceMutation]
    gene_sequences: Dict[str, str]  # gene_name -> protein sequence
    gene_id_map: Dict[str, int]  # gene_name -> gene_id
    species_map: Dict[str, int]  # species_name -> species_id
    
class ResistanceDatabaseBuilder:
    """Build GPU-ready resistance database"""
    
    def __init__(self):
        self.mutations = []
        self.gene_sequences = {}
        self.gene_id_map = {}
        self.species_map = {}
        self.next_gene_id = 0
        self.next_species_id = 0
        
        # Initialize with known FQ resistance genes
        self._init_default_genes()
        
    def _init_default_genes(self):
        """Initialize default gene mappings"""
        default_genes = ['gyrA', 'gyrB', 'parC', 'parE', 'grlA', 'grlB']
        for gene in default_genes:
            self.gene_id_map[gene] = self.next_gene_id
            self.next_gene_id += 1
            
        # Default species
        default_species = [
            'Escherichia_coli',
            'Klebsiella_pneumoniae', 
            'Pseudomonas_aeruginosa',
            'Staphylococcus_aureus',
            'Enterococcus_faecium',
            'Enterococcus_faecalis'
        ]
        for species in default_species:
            self.species_map[species] = self.next_species_id
            self.next_species_id += 1
    
    def load_known_mutations_csv(self, csv_path: str):
        """Load mutations from CSV file"""
        import pandas as pd
        
        logger.info(f"Loading mutations from {csv_path}")
        df = pd.read_csv(csv_path)
        
        for _, row in df.iterrows():
            gene_name = row['Gene']
            
            # Parse position and mutation
            mutation_str = row['Mutation']
            wildtype = mutation_str[0]
            position = int(mutation_str[1:-1])
            mutant = mutation_str[-1]
            
            # Get or create gene ID
            if gene_name not in self.gene_id_map:
                self.gene_id_map[gene_name] = self.next_gene_id
                self.next_gene_id += 1
            
            # Parse drugs affected
            drugs = [d.strip() for d in row.get('Drugs_Affected', 'fluoroquinolones').split(',')]
            
            # Resistance level based on MIC increase
            mic_increase = row.get('MIC_Increase', '4-fold')
            if '64' in str(mic_increase) or 'high' in str(mic_increase).lower():
                resistance_level = 0.9
            elif '16' in str(mic_increase) or '32' in str(mic_increase):
                resistance_level = 0.7
            elif '8' in str(mic_increase):
                resistance_level = 0.5
            else:
                resistance_level = 0.3
            
            # Check if we already have this position
            existing = None
            for mut in self.mutations:
                if mut.gene_name == gene_name and mut.position == position:
                    existing = mut
                    break
            
            if existing:
                if mutant not in existing.resistant_aas:
                    existing.resistant_aas.append(mutant)
            else:
                self.mutations.append(ResistanceMutation(
                    gene_name=gene_name,
                    gene_id=self.gene_id_map[gene_name],
                    position=position,
                    wildtype_aa=wildtype,
                    resistant_aas=[mutant],
                    resistance_level=resistance_level,
                    drugs_affected=drugs
                ))
        
        logger.info(f"Loaded {len(self.mutations)} unique mutation positions")
    
    def load_gene_sequences(self, fasta_dir: str):
        """Load protein sequences from FASTA files"""
        fasta_path = Path(fasta_dir)
        
        for fasta_file in fasta_path.glob("*.fasta"):
            logger.info(f"Loading sequences from {fasta_file}")
            
            for record in SeqIO.parse(fasta_file, "fasta"):
                # Extract gene name from header
                gene_name = None
                if 'gyrA' in record.description:
                    gene_name = 'gyrA'
                elif 'gyrB' in record.description:
                    gene_name = 'gyrB'
                elif 'parC' in record.description:
                    gene_name = 'parC'
                elif 'parE' in record.description:
                    gene_name = 'parE'
                elif 'grlA' in record.description:
                    gene_name = 'grlA'
                
                if gene_name and gene_name not in self.gene_sequences:
                    self.gene_sequences[gene_name] = str(record.seq)
                    logger.info(f"Loaded {gene_name}: {len(record.seq)} aa")
    
    def add_manual_mutations(self):
        """Add well-characterized resistance mutations manually"""
        
        # E. coli gyrA mutations
        ecoli_gyra = [
            ResistanceMutation('gyrA', 0, 83, 'S', ['L', 'F', 'W'], 0.9, 
                             ['ciprofloxacin', 'levofloxacin', 'moxifloxacin']),
            ResistanceMutation('gyrA', 0, 87, 'D', ['N', 'G', 'Y', 'H'], 0.7,
                             ['ciprofloxacin', 'levofloxacin']),
        ]
        
        # E. coli parC mutations  
        ecoli_parc = [
            ResistanceMutation('parC', 2, 80, 'S', ['I', 'R'], 0.6,
                             ['ciprofloxacin', 'levofloxacin']),
            ResistanceMutation('parC', 2, 84, 'E', ['V', 'K', 'G'], 0.5,
                             ['ciprofloxacin']),
        ]
        
        # P. aeruginosa mutations
        pa_gyra = [
            ResistanceMutation('gyrA', 0, 83, 'T', ['I'], 0.9,
                             ['ciprofloxacin', 'levofloxacin']),
            ResistanceMutation('gyrA', 0, 87, 'D', ['N', 'Y'], 0.7,
                             ['ciprofloxacin']),
        ]
        
        # S. aureus mutations (grlA = parC homolog)
        sa_grla = [
            ResistanceMutation('grlA', 4, 80, 'S', ['F', 'Y'], 0.9,
                             ['ciprofloxacin', 'levofloxacin', 'moxifloxacin']),
            ResistanceMutation('grlA', 4, 84, 'E', ['K', 'V'], 0.6,
                             ['ciprofloxacin']),
        ]
        
        # Add all manual mutations
        for mut_list in [ecoli_gyra, ecoli_parc, pa_gyra, sa_grla]:
            for mut in mut_list:
                # Check if already exists
                exists = False
                for existing in self.mutations:
                    if (existing.gene_name == mut.gene_name and 
                        existing.position == mut.position):
                        # Merge resistant AAs
                        for aa in mut.resistant_aas:
                            if aa not in existing.resistant_aas:
                                existing.resistant_aas.append(aa)
                        exists = True
                        break
                
                if not exists:
                    self.mutations.append(mut)
        
        logger.info(f"Total mutations after manual additions: {len(self.mutations)}")
    
    def generate_resistant_variants(self, output_dir: str):
        """Generate protein sequences with resistance mutations"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        variants = []
        
        for mutation in self.mutations:
            if mutation.gene_name not in self.gene_sequences:
                logger.warning(f"No sequence found for {mutation.gene_name}")
                continue
                
            wt_seq = self.gene_sequences[mutation.gene_name]
            
            # Check position validity
            if mutation.position > len(wt_seq):
                logger.warning(f"Position {mutation.position} exceeds {mutation.gene_name} length")
                continue
            
            # Verify wildtype matches
            actual_wt = wt_seq[mutation.position - 1]
            if actual_wt != mutation.wildtype_aa:
                logger.warning(f"{mutation.gene_name} position {mutation.position}: "
                             f"expected {mutation.wildtype_aa}, found {actual_wt}")
            
            # Generate each variant
            for resistant_aa in mutation.resistant_aas:
                variant_seq = list(wt_seq)
                variant_seq[mutation.position - 1] = resistant_aa
                variant_seq = ''.join(variant_seq)
                
                variant_id = f"{mutation.gene_name}_{mutation.wildtype_aa}{mutation.position}{resistant_aa}"
                
                record = SeqRecord(
                    Seq(variant_seq),
                    id=variant_id,
                    description=f"Resistance: {','.join(mutation.drugs_affected)}"
                )
                variants.append(record)
        
        # Write variants
        variant_file = output_path / "resistant_variants.fasta"
        SeqIO.write(variants, variant_file, "fasta")
        logger.info(f"Generated {len(variants)} resistant variant sequences")
        
        # Also write wildtype sequences
        wt_records = []
        for gene_name, seq in self.gene_sequences.items():
            record = SeqRecord(
                Seq(seq),
                id=f"{gene_name}_WT",
                description="Wildtype sequence"
            )
            wt_records.append(record)
        
        wt_file = output_path / "wildtype_sequences.fasta"
        SeqIO.write(wt_records, wt_file, "fasta")
        logger.info(f"Wrote {len(wt_records)} wildtype sequences")
    
    def build_minimizer_index(self, k: int = 8, w: int = 4) -> Dict[int, List[Tuple[int, int]]]:
        """Build minimizer index for sequences"""
        minimizer_index = {}
        
        for gene_name, seq in self.gene_sequences.items():
            gene_id = self.gene_id_map[gene_name]
            
            for i in range(len(seq) - k - w + 1):
                # Find minimizer in window
                min_hash = float('inf')
                min_pos = -1
                
                for j in range(w):
                    kmer = seq[i + j:i + j + k]
                    # Simple hash function
                    kmer_hash = hash(kmer)
                    
                    if kmer_hash < min_hash:
                        min_hash = kmer_hash
                        min_pos = i + j
                
                if min_hash not in minimizer_index:
                    minimizer_index[min_hash] = []
                
                minimizer_index[min_hash].append((gene_id, min_pos))
        
        return minimizer_index
    
    def export_binary_database(self, output_file: str):
        """Export database in binary format for GPU"""
        
        with open(output_file, 'wb') as f:
            # Write header
            f.write(struct.pack('I', len(self.mutations)))  # num_mutations
            f.write(struct.pack('I', len(self.gene_sequences)))  # num_genes
            
            # Write mutations
            for mut in self.mutations:
                f.write(struct.pack('I', mut.gene_id))
                f.write(struct.pack('H', mut.position))
                f.write(mut.wildtype_aa.encode('ascii'))
                
                # Write resistant AAs (pad to 10)
                for i in range(10):
                    if i < len(mut.resistant_aas):
                        f.write(mut.resistant_aas[i].encode('ascii'))
                    else:
                        f.write(b'\0')
                
                f.write(struct.pack('B', len(mut.resistant_aas)))
                f.write(struct.pack('f', mut.resistance_level))
                
                # Drug name (pad to 32 chars)
                drug_name = mut.drugs_affected[0] if mut.drugs_affected else "fluoroquinolone"
                drug_bytes = drug_name.encode('ascii')[:32]
                f.write(drug_bytes.ljust(32, b'\0'))
            
            # Write gene sequences
            for gene_name in sorted(self.gene_sequences.keys()):
                seq = self.gene_sequences[gene_name]
                f.write(struct.pack('I', len(seq)))
                f.write(seq.encode('ascii'))
        
        logger.info(f"Exported binary database to {output_file}")
    
    def export_json_database(self, output_file: str):
        """Export database as JSON for debugging"""
        
        db_dict = {
            'metadata': {
                'version': '1.0',
                'num_mutations': len(self.mutations),
                'num_genes': len(self.gene_sequences),
                'genes': list(self.gene_id_map.keys()),
                'species': list(self.species_map.keys())
            },
            'gene_id_map': self.gene_id_map,
            'species_map': self.species_map,
            'mutations': [
                {
                    'gene': mut.gene_name,
                    'gene_id': mut.gene_id,
                    'position': mut.position,
                    'wildtype': mut.wildtype_aa,
                    'resistant': mut.resistant_aas,
                    'resistance_level': mut.resistance_level,
                    'drugs': mut.drugs_affected
                }
                for mut in self.mutations
            ],
            'gene_lengths': {
                gene: len(seq) for gene, seq in self.gene_sequences.items()
            }
        }
        
        with open(output_file, 'w') as f:
            json.dump(db_dict, f, indent=2)
        
        logger.info(f"Exported JSON database to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Build GPU-ready resistance database')
    parser.add_argument('--csv', help='Known mutations CSV file')
    parser.add_argument('--fasta-dir', help='Directory with gene FASTA files')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--add-manual', action='store_true', 
                       help='Add manually curated mutations')
    
    args = parser.parse_args()
    
    # Create output directory
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Build database
    builder = ResistanceDatabaseBuilder()
    
    if args.csv:
        builder.load_known_mutations_csv(args.csv)
    
    if args.fasta_dir:
        builder.load_gene_sequences(args.fasta_dir)
    
    if args.add_manual:
        builder.add_manual_mutations()
    
    # Generate outputs
    builder.generate_resistant_variants(args.output_dir)
    builder.export_binary_database(output_path / 'resistance_db.bin')
    builder.export_json_database(output_path / 'resistance_db.json')
    
    # Build and save minimizer index
    minimizer_index = builder.build_minimizer_index()
    logger.info(f"Built minimizer index with {len(minimizer_index)} unique minimizers")
    
    print("\nDatabase Summary:")
    print(f"Total mutations: {len(builder.mutations)}")
    print(f"Genes: {list(builder.gene_id_map.keys())}")
    print(f"Species: {list(builder.species_map.keys())}")
    
    # Print example mutations
    print("\nExample mutations:")
    for mut in builder.mutations[:5]:
        print(f"  {mut.gene_name} {mut.wildtype_aa}{mut.position}"
              f"{'/'.join(mut.resistant_aas)} - {mut.drugs_affected[0]}")

if __name__ == '__main__':
    main()