#!/usr/bin/env python3
"""
BioGPU Fluoroquinolone Resistance Database Builder
Processes resistance genes and mutations to create a comprehensive database
"""

import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Set, Optional
import json
import struct
import pickle
import logging
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import requests
import time

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set up Entrez
Entrez.email = "your.email@example.com"  # Change this to your email

@dataclass
class ResistanceMutation:
    """Represents a point mutation conferring fluoroquinolone resistance"""
    species: str
    gene: str
    position: int
    wild_type: str
    mutant: str
    protein_accession: str
    nucleotide_start: int
    nucleotide_stop: int
    resistance_level: str = "MODERATE"  # Can be LOW, MODERATE, HIGH
    drugs_affected: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        # Map common mutations to resistance levels
        if self.gene in ['gyrA', 'parC']:
            if self.position in [83, 87] and self.gene == 'gyrA':
                self.resistance_level = "HIGH"
            elif self.position in [80, 84] and self.gene == 'parC':
                self.resistance_level = "HIGH"
        
        # Set default drugs affected
        if not self.drugs_affected:
            self.drugs_affected = [
                "ciprofloxacin", "levofloxacin", "moxifloxacin", 
                "norfloxacin", "ofloxacin"
            ]

@dataclass
class ResistanceGene:
    """Represents a gene conferring fluoroquinolone resistance"""
    species: str
    gene_name: str
    gene_family: str
    protein_accession: str
    nucleotide_start: int
    nucleotide_stop: int
    mechanism: str = "unknown"
    
    def __post_init__(self):
        # Determine mechanism based on gene family
        if self.gene_family.startswith('qnr'):
            self.mechanism = "target_protection"
        elif self.gene_family == 'aac(6\')-Ib-cr':
            self.mechanism = "antibiotic_inactivation"
        elif self.gene_family in ['oqxA', 'oqxB', 'qepA']:
            self.mechanism = "efflux"

@dataclass
class QRDRRegion:
    """Quinolone Resistance-Determining Region"""
    gene: str
    species: str
    start_codon: int
    end_codon: int
    critical_positions: List[int]
    
class FluoroquinoloneResistanceDB:
    def __init__(self):
        self.mutations: List[ResistanceMutation] = []
        self.genes: List[ResistanceGene] = []
        self.qrdr_regions: List[QRDRRegion] = []
        self.species_set: Set[str] = set()
        
        # Define QRDR regions for common species
        self._define_qrdr_regions()
        
    def _define_qrdr_regions(self):
        """Define known QRDR regions for various species"""
        # E. coli and similar Enterobacteriaceae
        self.qrdr_regions.extend([
            QRDRRegion("gyrA", "Escherichia coli", 67, 106, [83, 87]),
            QRDRRegion("gyrB", "Escherichia coli", 426, 466, [426, 447]),
            QRDRRegion("parC", "Escherichia coli", 64, 102, [80, 84]),
            QRDRRegion("parE", "Escherichia coli", 410, 470, [416, 445]),
        ])
        
        # S. aureus
        self.qrdr_regions.extend([
            QRDRRegion("gyrA", "Staphylococcus aureus", 68, 108, [84, 88]),
            QRDRRegion("grlA", "Staphylococcus aureus", 64, 102, [80, 84]),  # parC equivalent
        ])
        
        # P. aeruginosa
        self.qrdr_regions.extend([
            QRDRRegion("gyrA", "Pseudomonas aeruginosa", 67, 106, [83, 87]),
            QRDRRegion("parC", "Pseudomonas aeruginosa", 69, 107, [87]),
        ])
        
    def load_mutations_from_csv(self, csv_path: str):
        """Load mutations from CSV file"""
        logger.info(f"Loading mutations from {csv_path}")
        df = pd.read_csv(csv_path)
        
        for _, row in df.iterrows():
            # Skip entries with missing mutation info
            if pd.isna(row['wt']) or pd.isna(row['mut']):
                continue
                
            mutation = ResistanceMutation(
                species=row['species'],
                gene=row['gene'],
                position=int(row['location']),
                wild_type=row['wt'],
                mutant=row['mut'],
                protein_accession=row['Protein'],
                nucleotide_start=int(row['start']),
                nucleotide_stop=int(row['stop'])
            )
            
            self.mutations.append(mutation)
            self.species_set.add(row['species'])
            
        logger.info(f"Loaded {len(self.mutations)} mutations from {len(self.species_set)} species")
        
    def load_genes_from_csv(self, csv_path: str):
        """Load resistance genes from CSV file"""
        logger.info(f"Loading resistance genes from {csv_path}")
        df = pd.read_csv(csv_path)
        
        for _, row in df.iterrows():
            # Extract gene family from gene name
            gene_family = row['gene'].split('_')[0] if '_' in row['gene'] else row['gene'][:3]
            
            gene = ResistanceGene(
                species=row['species'],
                gene_name=row['gene'],
                gene_family=gene_family,
                protein_accession=row['Protein'],
                nucleotide_start=int(row['start']),
                nucleotide_stop=int(row['stop'])
            )
            
            self.genes.append(gene)
            self.species_set.add(row['species'])
            
        logger.info(f"Loaded {len(self.genes)} resistance genes")
        
    def create_mutation_profiles(self) -> Dict[str, Dict]:
        """Create species-specific mutation profiles"""
        profiles = {}
        
        for species in self.species_set:
            species_mutations = [m for m in self.mutations if m.species == species]
            species_genes = [g for g in self.genes if g.species == species]
            
            profile = {
                'species': species,
                'total_mutations': len(species_mutations),
                'total_genes': len(species_genes),
                'mutations_by_gene': {},
                'resistance_genes': {}
            }
            
            # Group mutations by gene
            for mutation in species_mutations:
                if mutation.gene not in profile['mutations_by_gene']:
                    profile['mutations_by_gene'][mutation.gene] = []
                
                profile['mutations_by_gene'][mutation.gene].append({
                    'position': mutation.position,
                    'wild_type': mutation.wild_type,
                    'mutant': mutation.mutant,
                    'resistance_level': mutation.resistance_level
                })
            
            # Add resistance genes
            for gene in species_genes:
                profile['resistance_genes'][gene.gene_name] = {
                    'family': gene.gene_family,
                    'mechanism': gene.mechanism
                }
            
            profiles[species] = profile
            
        return profiles
        
    def create_gpu_optimized_structures(self):
        """Create GPU-optimized data structures for exact mutation detection"""
        # Create mutation lookup table
        num_mutations = len(self.mutations)
        
        # Arrays for GPU
        gene_ids = np.zeros(num_mutations, dtype=np.uint32)
        positions = np.zeros(num_mutations, dtype=np.uint32)
        wild_type_encoded = np.zeros(num_mutations, dtype=np.uint8)
        mutant_encoded = np.zeros(num_mutations, dtype=np.uint8)
        resistance_levels = np.zeros(num_mutations, dtype=np.uint8)
        
        # Create gene name to ID mapping
        unique_genes = list(set(m.gene for m in self.mutations))
        gene_to_id = {gene: i for i, gene in enumerate(unique_genes)}
        
        # Amino acid encoding
        aa_to_code = {
            'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7,
            'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
            'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19
        }
        
        resistance_to_level = {'LOW': 0, 'MODERATE': 1, 'HIGH': 2}
        
        # Fill arrays
        for i, mutation in enumerate(self.mutations):
            gene_ids[i] = gene_to_id[mutation.gene]
            positions[i] = mutation.position
            wild_type_encoded[i] = aa_to_code.get(mutation.wild_type, 255)
            mutant_encoded[i] = aa_to_code.get(mutation.mutant, 255)
            resistance_levels[i] = resistance_to_level[mutation.resistance_level]
            
        return {
            'gene_ids': gene_ids,
            'positions': positions,
            'wild_type': wild_type_encoded,
            'mutant': mutant_encoded,
            'resistance_levels': resistance_levels,
            'gene_mapping': gene_to_id,
            'aa_encoding': aa_to_code
        }
    
    def create_target_regions(self) -> Dict[str, Dict]:
        """Create target regions for exact alignment focusing on QRDRs"""
        target_regions = {}
        
        for qrdr in self.qrdr_regions:
            key = f"{qrdr.species}_{qrdr.gene}"
            
            # Calculate nucleotide positions (codon * 3)
            nt_start = (qrdr.start_codon - 1) * 3
            nt_end = qrdr.end_codon * 3
            
            # Add flanking regions for better alignment (50bp each side)
            flank_size = 50
            
            target_regions[key] = {
                'gene': qrdr.gene,
                'species': qrdr.species,
                'protein_start': qrdr.start_codon,
                'protein_end': qrdr.end_codon,
                'nucleotide_start': nt_start - flank_size,
                'nucleotide_end': nt_end + flank_size,
                'critical_codons': qrdr.critical_positions,
                'mutations': []
            }
            
            # Add known mutations for this region
            for mutation in self.mutations:
                if (mutation.gene == qrdr.gene and 
                    mutation.species == qrdr.species and
                    qrdr.start_codon <= mutation.position <= qrdr.end_codon):
                    
                    target_regions[key]['mutations'].append({
                        'position': mutation.position,
                        'wild_type': mutation.wild_type,
                        'mutant': mutation.mutant,
                        'resistance_level': mutation.resistance_level,
                        'codon_position': mutation.position - qrdr.start_codon
                    })
        
        return target_regions
    
    def create_alignment_templates(self) -> Dict[str, np.ndarray]:
        """Create pre-computed alignment scoring matrices for GPU"""
        templates = {}
        
        # DNA substitution matrix for exact matching
        # Higher penalties for mutations at critical positions
        dna_match = 2
        dna_mismatch = -3
        gap_open = -5
        gap_extend = -2
        
        # Create base scoring matrix
        base_matrix = np.full((4, 4), dna_mismatch, dtype=np.float32)
        np.fill_diagonal(base_matrix, dna_match)
        
        templates['base_scoring'] = base_matrix
        templates['gap_open'] = gap_open
        templates['gap_extend'] = gap_extend
        
        # Create position-specific scoring adjustments for critical codons
        for key, region in self.create_target_regions().items():
            if region['mutations']:
                # Create position weight matrix
                region_length = region['nucleotide_end'] - region['nucleotide_start']
                position_weights = np.ones(region_length, dtype=np.float32)
                
                # Increase weight at mutation positions
                for mutation in region['mutations']:
                    codon_nt_pos = (mutation['position'] - 1) * 3
                    rel_pos = codon_nt_pos - region['nucleotide_start'] + 50  # Account for flank
                    
                    # Weight the entire codon
                    for i in range(3):
                        if 0 <= rel_pos + i < region_length:
                            position_weights[rel_pos + i] = 2.0  # Double weight for critical positions
                
                templates[f"{key}_weights"] = position_weights
        
        return templates
        
    def save_database(self, output_dir: str):
        """Save the complete resistance database"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save mutations
        with open(output_path / 'mutations.pkl', 'wb') as f:
            pickle.dump(self.mutations, f)
            
        # Save genes
        with open(output_path / 'genes.pkl', 'wb') as f:
            pickle.dump(self.genes, f)
            
        # Save QRDR regions
        with open(output_path / 'qrdr_regions.pkl', 'wb') as f:
            pickle.dump(self.qrdr_regions, f)
            
        # Save species profiles as JSON for readability
        profiles = self.create_mutation_profiles()
        with open(output_path / 'species_profiles.json', 'w') as f:
            json.dump(profiles, f, indent=2)
            
        # Save GPU-optimized structures
        gpu_data = self.create_gpu_optimized_structures()
        np.savez(
            output_path / 'gpu_mutation_data.npz',
            **gpu_data
        )
        
        # Save target regions for exact alignment
        target_regions = self.create_target_regions()
        with open(output_path / 'target_regions.json', 'w') as f:
            json.dump(target_regions, f, indent=2)
            
        # Save alignment templates
        alignment_templates = self.create_alignment_templates()
        np.savez(
            output_path / 'alignment_templates.npz',
            **alignment_templates
        )
        
        # Save summary statistics
        stats = {
            'total_mutations': len(self.mutations),
            'total_genes': len(self.genes),
            'total_species': len(self.species_set),
            'total_target_regions': len(target_regions),
            'species_list': sorted(list(self.species_set)),
            'gene_families': sorted(list(set(g.gene_family for g in self.genes))),
            'mutation_genes': sorted(list(set(m.gene for m in self.mutations))),
            'alignment_params': {
                'match_score': 2,
                'mismatch_penalty': -3,
                'gap_open': -5,
                'gap_extend': -2,
                'critical_position_weight': 2.0
            }
        }
        
        with open(output_path / 'database_stats.json', 'w') as f:
            json.dump(stats, f, indent=2)
            
        logger.info(f"Database saved to {output_path}")
        
    def generate_pssm_for_qrdr(self, gene: str, species: str) -> np.ndarray:
        """Generate Position-Specific Scoring Matrix for QRDR regions"""
        # Find relevant mutations
        relevant_mutations = [
            m for m in self.mutations 
            if m.gene == gene and m.species == species
        ]
        
        if not relevant_mutations:
            return None
            
        # Find QRDR region
        qrdr = next(
            (q for q in self.qrdr_regions if q.gene == gene and q.species == species),
            None
        )
        
        if not qrdr:
            return None
            
        # Create PSSM (20 amino acids x region length)
        region_length = qrdr.end_codon - qrdr.start_codon + 1
        pssm = np.zeros((20, region_length), dtype=np.float32)
        
        # Add pseudocounts
        pssm += 0.1
        
        # Add mutation information
        for mutation in relevant_mutations:
            if qrdr.start_codon <= mutation.position <= qrdr.end_codon:
                pos_idx = mutation.position - qrdr.start_codon
                aa_idx = ord(mutation.mutant) - ord('A')  # Simplified encoding
                
                # Increase score for resistance mutations
                if mutation.resistance_level == "HIGH":
                    pssm[aa_idx, pos_idx] += 10.0
                elif mutation.resistance_level == "MODERATE":
                    pssm[aa_idx, pos_idx] += 5.0
                else:
                    pssm[aa_idx, pos_idx] += 2.0
                    
        # Normalize
        pssm = pssm / pssm.sum(axis=0, keepdims=True)
        
        return pssm

def main():
    """Main function to build the resistance database"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Build fluoroquinolone resistance database')
    parser.add_argument('--mutations', type=str, required=True, help='Path to mutations CSV')
    parser.add_argument('--genes', type=str, required=True, help='Path to resistance genes CSV')
    parser.add_argument('--output', type=str, default='data/resistance_db', help='Output directory')
    
    args = parser.parse_args()
    
    # Build database
    db = FluoroquinoloneResistanceDB()
    db.load_mutations_from_csv(args.mutations)
    db.load_genes_from_csv(args.genes)
    
    # Generate some example PSSMs
    logger.info("Generating example PSSMs for E. coli gyrA and parC")
    ecoli_gyra_pssm = db.generate_pssm_for_qrdr("gyrA", "Escherichia coli")
    ecoli_parc_pssm = db.generate_pssm_for_qrdr("parC", "Escherichia coli")
    
    # Save everything
    db.save_database(args.output)
    
    # Print summary
    profiles = db.create_mutation_profiles()
    print("\n=== Fluoroquinolone Resistance Database Summary ===")
    print(f"Total species: {len(db.species_set)}")
    print(f"Total mutations: {len(db.mutations)}")
    print(f"Total resistance genes: {len(db.genes)}")
    print("\nTop 5 species by mutation count:")
    
    species_counts = [(s, p['total_mutations']) for s, p in profiles.items()]
    species_counts.sort(key=lambda x: x[1], reverse=True)
    
    for species, count in species_counts[:5]:
        print(f"  {species}: {count} mutations")

if __name__ == "__main__":
    main()