#!/usr/bin/env python3
"""
Test script to check if WP_ accession handling works
"""

import sys
sys.path.append('.')

from Bio import Entrez
from fq_reference_downloader import ResistanceGeneDownloader

# Set your email for NCBI Entrez
Entrez.email = "test@example.com"  # Change this!

def test_wp_accession():
    """Test the WP_ accession handling"""
    
    # Create a minimal downloader instance
    downloader = ResistanceGeneDownloader(
        mutations_csv="../tools/Quinolone_resistance_mutation_table.csv",
        genes_csv="../tools/Quinolone_resistance_gene_table.csv",
        flanking_size=100  # Smaller for testing
    )
    
    # Test cases
    test_cases = [
        ("WP_064768701.1", "Escherichia coli", "qepA7"),
        ("WP_014386436.1", "Salmonella enterica", "aac(6')-Ib-cr"),
        ("WP_000185469.1", "Klebsiella pneumoniae", "qnrB1")
    ]
    
    print("Testing WP_ accession handling...\n")
    
    for protein_acc, species, gene in test_cases:
        print(f"Testing {protein_acc} ({species} {gene})...")
        
        # Test the new method directly
        nuc_acc, start, end, is_complement = downloader.find_nucleotide_for_wp_accession(
            protein_acc, species, gene
        )
        
        if nuc_acc:
            print(f"  ✓ Found: {nuc_acc}:{start}..{end} (complement: {is_complement})")
        else:
            print(f"  ✗ Not found")
        
        # Also test the full process
        ref_gene = downloader.process_protein_accession(
            protein_acc, species, gene, 
            {'start': 99498, 'stop': 101033}  # Example coordinates from CSV
        )
        
        if ref_gene:
            print(f"  ✓ Full processing successful")
            print(f"    - Nucleotide: {ref_gene.nucleotide_accession}")
            print(f"    - Gene length: {ref_gene.gene_length} bp")
            print(f"    - Protein length: {len(ref_gene.protein_sequence)} aa")
        else:
            print(f"  ✗ Full processing failed")
        
        print()

if __name__ == "__main__":
    test_wp_accession()