#!/usr/bin/env python3
"""
Check protein alignments from CSV to see if mutations are present
"""

import csv
import sys
from collections import Counter

def load_reference_proteins(protein_db_path):
    """Load reference protein sequences from the database"""
    proteins = {}
    
    # Read protein details JSON to get metadata
    import json
    with open(f"{protein_db_path}/protein_details.json", 'r') as f:
        protein_details = json.load(f)
    
    # Read protein sequences
    with open(f"{protein_db_path}/proteins.bin", 'rb') as f:
        # Skip the count
        num_proteins_bytes = f.read(4)
        num_proteins = int.from_bytes(num_proteins_bytes, 'little')
        
        # Read all remaining sequences
        all_sequences = f.read().decode('ascii', errors='ignore')
    
    # Parse sequences based on metadata
    current_offset = 0
    for i, protein in enumerate(protein_details):
        length = protein['length']
        sequence = all_sequences[current_offset:current_offset+length]
        current_offset += length
        
        key = (protein['species'], protein['gene'])
        proteins[key] = {
            'sequence': sequence,
            'id': protein['id'],
            'gene_id': protein['gene_id'],
            'species_id': protein['species_id']
        }
    
    return proteins

def compare_sequences(query_peptide, ref_sequence, ref_start):
    """Compare query peptide to reference sequence and find mutations"""
    mutations = []
    
    for i, (q_aa, r_aa) in enumerate(zip(query_peptide, ref_sequence[ref_start:ref_start+len(query_peptide)])):
        if q_aa != r_aa and q_aa != 'X' and q_aa != '*':
            mutations.append({
                'position': ref_start + i + 1,  # 1-based position
                'ref_aa': r_aa,
                'query_aa': q_aa
            })
    
    return mutations

def main():
    if len(sys.argv) < 3:
        print("Usage: python check_protein_alignments.py <csv_file> <protein_db_path>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    protein_db_path = sys.argv[2]
    
    # Load reference proteins
    print("Loading reference proteins...")
    ref_proteins = load_reference_proteins(protein_db_path)
    print(f"Loaded {len(ref_proteins)} reference proteins")
    
    # Process CSV
    print(f"\nAnalyzing alignments from {csv_file}...")
    
    total_alignments = 0
    alignments_with_mutations = 0
    total_mutations = 0
    gene_counter = Counter()
    species_counter = Counter()
    mutation_positions = Counter()
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for i, row in enumerate(reader):
            if i >= 10000:  # Sample first 10k alignments
                break
                
            total_alignments += 1
            
            species = row['species_name']
            gene = row['gene_name']
            query_peptide = row['query_peptide']
            ref_start = int(row['ref_start'])
            
            gene_counter[gene] += 1
            species_counter[species] += 1
            
            # Get reference sequence
            key = (species, gene)
            if key in ref_proteins:
                ref_seq = ref_proteins[key]['sequence']
                
                # Compare sequences
                mutations = compare_sequences(query_peptide, ref_seq, ref_start)
                
                if mutations:
                    alignments_with_mutations += 1
                    total_mutations += len(mutations)
                    
                    if alignments_with_mutations <= 10:  # Show first 10 examples
                        print(f"\nAlignment {i}: {species} {gene}")
                        print(f"  Query: {query_peptide}")
                        print(f"  Ref:   {ref_seq[ref_start:ref_start+len(query_peptide)]}")
                        print(f"  Mutations: ", end="")
                        for mut in mutations:
                            print(f"{mut['ref_aa']}{mut['position']}{mut['query_aa']} ", end="")
                            mutation_positions[f"{gene}:{mut['position']}"] += 1
                        print()
    
    # Summary
    print(f"\n=== SUMMARY (first {total_alignments} alignments) ===")
    print(f"Total alignments: {total_alignments}")
    print(f"Alignments with mutations: {alignments_with_mutations} ({100*alignments_with_mutations/total_alignments:.1f}%)")
    print(f"Total mutations found: {total_mutations}")
    
    print(f"\nGene distribution:")
    for gene, count in gene_counter.most_common():
        print(f"  {gene}: {count}")
    
    print(f"\nSpecies distribution:")
    for species, count in species_counter.most_common():
        print(f"  {species}: {count}")
    
    if mutation_positions:
        print(f"\nTop mutation positions:")
        for pos, count in mutation_positions.most_common(10):
            print(f"  {pos}: {count} occurrences")

if __name__ == "__main__":
    main()