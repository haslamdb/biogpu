#!/usr/bin/env python3
"""
Process already downloaded genomes to extract k-mers
"""

import os
import glob
from Bio import SeqIO
import argparse

PRIORITY_PATHOGENS = {
    "Escherichia coli": {"taxid": "562"},
    "Klebsiella pneumoniae": {"taxid": "573"},
    "Proteus mirabilis": {"taxid": "584"},
    "Enterococcus faecalis": {"taxid": "1351"},
    "Pseudomonas aeruginosa": {"taxid": "287"},
    "Streptococcus pneumoniae": {"taxid": "1313"},
    "Haemophilus influenzae": {"taxid": "727"},
    "Moraxella catarrhalis": {"taxid": "480"},
    "Staphylococcus aureus": {"taxid": "1280"},
    "Clostridioides difficile": {"taxid": "1496"},
    "Salmonella enterica": {"taxid": "28901"},
    "Campylobacter jejuni": {"taxid": "197"},
    "Shigella flexneri": {"taxid": "623"},
    "Streptococcus agalactiae": {"taxid": "1311"},
    "Listeria monocytogenes": {"taxid": "1639"},
    "Neisseria meningitidis": {"taxid": "487"},
    "Enterobacter cloacae": {"taxid": "550"},
    "Citrobacter freundii": {"taxid": "546"},
    "Serratia marcescens": {"taxid": "615"},
    "Acinetobacter baumannii": {"taxid": "470"},
}

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def extract_kmers_from_genome(fasta_file, k=31, stride=1):
    kmers = set()
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            for i in range(0, len(seq) - k + 1, stride):
                kmer = seq[i:i+k]
                if all(base in 'ACGT' for base in kmer):
                    rc = reverse_complement(kmer)
                    canonical = min(kmer, rc)
                    kmers.add(canonical)
        return kmers
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
        return set()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--stride", type=int, default=5)
    args = parser.parse_args()
    
    genome_dir = os.path.join(args.output_dir, "genomes")
    
    # Map species names to taxids
    species_map = {}
    for species, info in PRIORITY_PATHOGENS.items():
        safe_name = species.replace(" ", "_")
        species_map[safe_name] = (species, info["taxid"])
    
    # Process all fasta files
    all_files = glob.glob(os.path.join(genome_dir, "*.fasta"))
    print(f"Found {len(all_files)} genome files")
    
    # Group by species
    species_files = {}
    for fasta_file in all_files:
        basename = os.path.basename(fasta_file)
        # Extract species from filename
        for safe_name, (species, taxid) in species_map.items():
            if basename.startswith(safe_name):
                if species not in species_files:
                    species_files[species] = []
                species_files[species].append(fasta_file)
                break
    
    print(f"\nSpecies summary:")
    for species, files in sorted(species_files.items()):
        print(f"  {species}: {len(files)} genomes")
    
    # Extract k-mers
    print(f"\nExtracting k-mers (k={args.k}, stride={args.stride})...")
    
    kmer_file = os.path.join(args.output_dir, "database_kmers.txt")
    species_summary = {}
    
    with open(kmer_file, 'w') as f:
        f.write("# K-mer database for BioGPU microbiome profiler\n")
        f.write(f"# K={args.k}, Stride={args.stride}\n")
        f.write("# Format: kmer<tab>taxon_id\n\n")
        
        for species, files in sorted(species_files.items()):
            taxid = PRIORITY_PATHOGENS[species]["taxid"]
            print(f"\nProcessing {species} ({len(files)} genomes)...")
            
            species_kmers = set()
            for i, fasta_file in enumerate(files):
                print(f"  [{i+1}/{len(files)}] {os.path.basename(fasta_file)}")
                kmers = extract_kmers_from_genome(fasta_file, k=args.k, stride=args.stride)
                species_kmers.update(kmers)
                
                # Write k-mers
                for kmer in kmers:
                    f.write(f"{kmer}\t{taxid}\n")
            
            species_summary[species] = {
                'taxid': taxid,
                'kmers': len(species_kmers),
                'genomes': len(files)
            }
    
    # Summary
    print("\n" + "="*70)
    print("K-mer Database Summary")
    print("="*70)
    
    total_kmers = 0
    for species, data in sorted(species_summary.items()):
        total_kmers += data['kmers']
        print(f"{species:<40} (taxid {data['taxid']:>6}): {data['kmers']:>10,} k-mers ({data['genomes']} genomes)")
    
    print(f"\nTotal unique k-mers: {total_kmers:,}")
    
    # Save summary
    summary_file = os.path.join(args.output_dir, "database_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("K-mer Database Summary\n")
        f.write("="*70 + "\n\n")
        
        for species, data in sorted(species_summary.items()):
            f.write(f"{species:<40} (taxid {data['taxid']:>6}): {data['kmers']:>10,} k-mers ({data['genomes']} genomes)\n")
        
        f.write(f"\nTotal unique k-mers: {total_kmers:,}\n")
    
    print(f"\nFiles written:")
    print(f"  - {kmer_file}")
    print(f"  - {summary_file}")

if __name__ == "__main__":
    main()