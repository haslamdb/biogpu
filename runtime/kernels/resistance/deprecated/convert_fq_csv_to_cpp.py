#!/usr/bin/env python3
"""
Convert FQ resistance CSV to C++ header file for hardcoding

This script converts a fluoroquinolone resistance mutations CSV file
into a C++ header file that can be compiled directly into the resistance pipeline.

Usage:
    python3 convert_fq_csv_to_cpp.py [csv_file] > fq_mutations_hardcoded.h

CSV Format:
    species,gene,location,wt,mut
    Escherichia_coli,gyrA,83,S,L
    ...

Author: BioGPU Pipeline
"""

import csv
import sys
import os

def convert_csv_to_cpp(csv_file):
    mutations = []
    line_num = 1  # Start at 1 for header
    
    if not os.path.exists(csv_file):
        print(f"// ERROR: File not found: {csv_file}", file=sys.stderr)
        sys.exit(1)
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            line_num += 1
            try:
                # Skip NA entries
                if row['location'] == 'NA' or row['wt'] == 'NA' or row['mut'] == 'NA':
                    continue
                
                # Validate amino acid codes
                if len(row['wt']) != 1 or len(row['mut']) != 1:
                    print(f"// WARNING: Line {line_num}: Invalid amino acid code length", file=sys.stderr)
                    continue
                
                mutations.append({
                    'species': row['species'],
                    'gene': row['gene'],
                    'position': int(row['location']),
                    'wt': row['wt'][0],  # Take first character
                    'mut': row['mut'][0]   # Take first character
                })
            except (KeyError, ValueError) as e:
                print(f"// WARNING: Line {line_num}: Error parsing row - {e}", file=sys.stderr)
                continue
    
    # Generate C++ header file
    print("#ifndef FQ_MUTATIONS_HARDCODED_H")
    print("#define FQ_MUTATIONS_HARDCODED_H")
    print()
    print("#include <vector>")
    print("#include <string>")
    print()
    print("// Forward declaration if not already included")
    print("struct FQResistanceMutation;")
    print()
    print("// Hardcoded FQ resistance mutations")
    print(f"// Generated from: {csv_file}")
    print(f"// Total mutations: {len(mutations)}")
    print("static const std::vector<FQResistanceMutation> HARDCODED_FQ_MUTATIONS = {")
    
    for i, mut in enumerate(mutations):
        comma = "," if i < len(mutations) - 1 else ""
        print(f'    {{"{mut["species"]}", "{mut["gene"]}", {mut["position"]}, '
              f'\'{mut["wt"]}\', \'{mut["mut"]}\'}}' + comma)
    
    print("};")
    print(f"\n// Total mutations: {len(mutations)}")
    
    # Also generate a summary by species and gene
    species_counts = {}
    gene_counts = {}
    
    for mut in mutations:
        species_counts[mut['species']] = species_counts.get(mut['species'], 0) + 1
        gene_counts[mut['gene']] = gene_counts.get(mut['gene'], 0) + 1
    
    print("\n// Summary by species:")
    for species, count in sorted(species_counts.items()):
        print(f"//   {species}: {count}")
    
    print("\n// Summary by gene:")
    for gene, count in sorted(gene_counts.items()):
        print(f"//   {gene}: {count}")
    
    print()
    print("#endif // FQ_MUTATIONS_HARDCODED_H")

def print_usage():
    print("""Usage: python3 convert_fq_csv_to_cpp.py [OPTIONS] [CSV_FILE]

Convert a fluoroquinolone resistance mutations CSV file to C++ code.

Options:
    -h, --help    Show this help message and exit

Arguments:
    CSV_FILE      Path to the mutations CSV file (default: 
                  /home/david/Documents/Code/biogpu/data/quinolone_resistance_mutation_table.csv)

Example:
    python3 convert_fq_csv_to_cpp.py my_mutations.csv > fq_mutations_hardcoded.h

CSV Format:
    The CSV file must have the following columns:
    - species: Species name with underscores (e.g., Escherichia_coli)
    - gene: Gene name (e.g., gyrA, parC)
    - location: Amino acid position (1-based)
    - wt: Wild-type amino acid
    - mut: Mutant amino acid

Output:
    Generates C++ header file to stdout. Redirect to save:
    > fq_mutations_hardcoded.h
""", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print_usage()
        sys.exit(0)
    
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = "/home/david/Documents/Code/biogpu/data/quinolone_resistance_mutation_table.csv"
    
    convert_csv_to_cpp(csv_file)