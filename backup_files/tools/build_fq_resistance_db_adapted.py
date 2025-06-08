#!/usr/bin/env python3
"""
Quick adapter to build gpu_resistance_db with the existing CSV format
"""

import pandas as pd
import sys
import os
sys.path.append('backup_scripts/tools')

from build_fq_resistance_db import FluoroquinoloneResistanceDB, ResistanceMutation, ResistanceGene
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_mutations_adapted(db, csv_path):
    """Load mutations from CSV with adapted column names"""
    logger.info(f"Loading mutations from {csv_path}")
    df = pd.read_csv(csv_path)
    
    for _, row in df.iterrows():
        # Skip entries with missing mutation info
        if pd.isna(row['Wild-type Amino Acid']) or pd.isna(row['Mutant Amino Acid']):
            continue
            
        mutation = ResistanceMutation(
            species=row['Species'],
            gene=row['Gene'],
            position=int(row['Mutation Position']),
            wild_type=row['Wild-type Amino Acid'],
            mutant=row['Mutant Amino Acid'],
            protein_accession="Unknown",  # Not in your CSV
            nucleotide_start=0,  # Not in your CSV
            nucleotide_stop=0    # Not in your CSV
        )
        
        db.mutations.append(mutation)
        db.species_set.add(row['Species'])
        
    logger.info(f"Loaded {len(db.mutations)} mutations from {len(db.species_set)} species")

def load_genes_adapted(db, csv_path):
    """Load resistance genes from CSV with adapted column names"""
    logger.info(f"Loading resistance genes from {csv_path}")
    
    if not os.path.exists(csv_path):
        logger.warning(f"Genes file {csv_path} not found, creating minimal gene set")
        # Create minimal gene set based on mutations
        unique_genes = set((m.species, m.gene) for m in db.mutations)
        
        for i, (species, gene) in enumerate(unique_genes):
            gene_obj = ResistanceGene(
                species=species,
                gene_name=gene,
                gene_family=gene,
                protein_accession="Unknown",
                nucleotide_start=0,
                nucleotide_stop=1000
            )
            db.genes.append(gene_obj)
            db.species_set.add(species)
        
        logger.info(f"Created {len(db.genes)} gene entries from mutation data")
        return
    
    df = pd.read_csv(csv_path)
    
    for _, row in df.iterrows():
        gene_name = row['Efflux Gene']
        gene_family = gene_name.split('_')[0] if '_' in gene_name else gene_name[:3]
        
        gene = ResistanceGene(
            species=row['Species'],
            gene_name=gene_name,
            gene_family=gene_family,
            protein_accession="Unknown",
            nucleotide_start=0,
            nucleotide_stop=1000
        )
        
        db.genes.append(gene)
        db.species_set.add(row['Species'])
        
    logger.info(f"Loaded {len(db.genes)} resistance genes")

def main():
    # Build database
    db = FluoroquinoloneResistanceDB()
    
    # Load with adapted functions
    load_mutations_adapted(db, 'data/Known_Quinolone_Changes.csv')
    load_genes_adapted(db, 'data/Known_Efflux_Pump_Genes.csv')
    
    # Save everything
    db.save_database('data/gpu_resistance_db')
    
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