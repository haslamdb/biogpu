#!/usr/bin/env python3
"""Analyze HDF5 file to understand species distribution and E. coli presence"""

import h5py
import numpy as np
import json
from collections import defaultdict

def analyze_hdf5_alignments(h5_file):
    """Analyze HDF5 file for species distribution"""
    
    with h5py.File(h5_file, 'r') as f:
        print("=== HDF5 File Analysis ===\n")
        
        # Check metadata for species mapping if available
        if 'metadata' in f:
            metadata = f['metadata']
            print("Metadata groups:", list(metadata.keys()))
            
            # Check for species mapping
            if 'species_map' in metadata:
                species_map = {}
                for key in metadata['species_map'].keys():
                    species_map[int(key)] = metadata['species_map'][key][()]
                print(f"\nSpecies mapping found: {len(species_map)} species")
                for sid, name in sorted(species_map.items())[:10]:
                    print(f"  {sid}: {name}")
        
        # Analyze translated search results
        if 'translated_search' in f:
            ts = f['translated_search']
            
            # Get all relevant data
            species_ids = ts['species_id'][:]
            gene_ids = ts['gene_id'][:]
            is_resistance = ts['is_resistance_mutation'][:]
            num_mutations = ts['num_mutations'][:]
            
            print(f"\n=== Translated Search Analysis ===")
            print(f"Total alignments: {len(species_ids):,}")
            
            # Count by species
            species_counts = defaultdict(int)
            species_resistance = defaultdict(int)
            species_mutations = defaultdict(list)
            
            for i in range(len(species_ids)):
                sid = species_ids[i]
                species_counts[sid] += 1
                
                if is_resistance[i]:
                    species_resistance[sid] += 1
                
                if num_mutations[i] > 0:
                    species_mutations[sid].append(num_mutations[i])
            
            print(f"\nSpecies distribution:")
            for sid, count in sorted(species_counts.items(), key=lambda x: -x[1])[:20]:
                resistance_count = species_resistance.get(sid, 0)
                mutation_count = len(species_mutations.get(sid, []))
                print(f"  Species {sid}: {count:,} alignments, "
                      f"{resistance_count} resistance mutations, "
                      f"{mutation_count} alignments with mutations")
            
            # Check specifically for E. coli (species ID 562 or 511145)
            ecoli_ids = [562, 511145]
            for ecoli_id in ecoli_ids:
                if ecoli_id in species_counts:
                    print(f"\nE. coli (species {ecoli_id}) found:")
                    print(f"  Total alignments: {species_counts[ecoli_id]:,}")
                    print(f"  Resistance mutations: {species_resistance.get(ecoli_id, 0)}")
                    print(f"  Alignments with mutations: {len(species_mutations.get(ecoli_id, []))}")
                else:
                    print(f"\nE. coli (species {ecoli_id}) NOT FOUND in alignments")
            
            # Check what the actual species IDs might be
            unique_species = np.unique(species_ids)
            print(f"\nTotal unique species: {len(unique_species)}")
            print(f"Species ID range: {np.min(unique_species)} - {np.max(unique_species)}")
            
            # If species IDs look strange, try to decode them
            if np.max(unique_species) > 1000000:
                print("\nWARNING: Species IDs appear to be corrupted or encoded differently")
                print("Sample species IDs:", unique_species[:10])

# Also check the JSON report
def analyze_json_report(json_file):
    """Analyze the JSON report for species information"""
    print("\n\n=== JSON Report Analysis ===\n")
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Count species from alignments
    species_counts = defaultdict(int)
    gene_counts = defaultdict(int)
    
    if 'alignments' in data:
        for alignment in data['alignments']:
            species_id = alignment.get('species_id', -1)
            gene_id = alignment.get('gene_id', -1)
            species_counts[species_id] += 1
            gene_counts[gene_id] += 1
    
        print(f"Total alignments in JSON: {len(data['alignments'])}")
        print(f"\nSpecies distribution from JSON:")
        for sid, count in sorted(species_counts.items(), key=lambda x: -x[1])[:10]:
            print(f"  Species {sid}: {count} alignments")
        
        # Check for E. coli
        ecoli_count = species_counts.get(562, 0) + species_counts.get(511145, 0)
        print(f"\nE. coli alignments in JSON: {ecoli_count}")

if __name__ == "__main__":
    # Analyze HDF5 file
    analyze_hdf5_alignments("results/1M_reads_test.h5")
    
    # Analyze JSON report
    analyze_json_report("results/1M_reads_test.json")