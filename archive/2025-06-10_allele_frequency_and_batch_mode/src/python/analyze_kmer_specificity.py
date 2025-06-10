#!/usr/bin/env python3
"""
Analyze the specificity of different k-mer sizes for protein search
"""

import random
from collections import defaultdict

def calculate_kmer_statistics():
    """Calculate expected matches for different k-mer sizes"""
    
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    
    print("K-mer Size Analysis for Protein Search")
    print("=" * 60)
    
    # Calculate search space
    for k in range(3, 12):
        total_kmers = 20 ** k
        
        # Estimate false positive rate
        # Assuming average protein length of 500 AA
        kmers_per_protein = 500 - k + 1
        
        # With 12 proteins in database
        total_db_kmers = 12 * kmers_per_protein
        
        # Probability of random match
        prob_match = total_db_kmers / total_kmers
        
        # Expected matches for a random k-mer
        expected_matches = prob_match
        
        print(f"\nK-mer size: {k}")
        print(f"  Total possible k-mers: {total_kmers:,}")
        print(f"  K-mers in database (~12 proteins): ~{total_db_kmers:,}")
        print(f"  Probability of random match: {prob_match:.6f}")
        print(f"  Expected false matches per query: {expected_matches:.4f}")
        
        # For translated search (6-frame translation)
        # A 150bp read gives ~50 AA per frame = 300 AA total
        # That's ~295 5-mers per read in 6-frame translation
        if k == 5:
            kmers_per_read = (300 - k + 1) * 1  # Simplified
            expected_hits_per_read = kmers_per_read * prob_match
            print(f"  Expected hits per 150bp read (6-frame): ~{expected_hits_per_read:.1f}")

    print("\n" + "=" * 60)
    
    # Simulate actual k-mer matches
    print("\nSimulating k-mer matches in a small protein database...")
    
    # Generate some random protein sequences
    proteins = []
    for i in range(12):
        length = random.randint(200, 800)
        seq = ''.join(random.choice(amino_acids) for _ in range(length))
        proteins.append(seq)
    
    # Test different k-mer sizes
    for k in [5, 6, 7, 8, 9]:
        # Extract all k-mers from proteins
        db_kmers = defaultdict(list)
        for pid, protein in enumerate(proteins):
            for i in range(len(protein) - k + 1):
                kmer = protein[i:i+k]
                db_kmers[kmer].append((pid, i))
        
        # Count k-mers that appear multiple times
        multi_hits = sum(1 for kmer, locs in db_kmers.items() if len(locs) > 1)
        unique_kmers = len(db_kmers)
        
        print(f"\nK={k}: {unique_kmers} unique k-mers, {multi_hits} appear multiple times")
        
        # Simulate a query
        query_kmers = set()
        for _ in range(50):  # 50 random k-mers to simulate a read
            query_kmer = ''.join(random.choice(amino_acids) for _ in range(k))
            query_kmers.add(query_kmer)
        
        # Count hits
        hits = sum(1 for qk in query_kmers if qk in db_kmers)
        print(f"  Random query: {hits}/{len(query_kmers)} k-mers hit database")

if __name__ == "__main__":
    calculate_kmer_statistics()
    
    print("\n\nRECOMMENDATION:")
    print("-" * 60)
    print("For protein-based translated search, k-mer sizes of 7-9 are typically used:")
    print("- K=7: Good balance of sensitivity and specificity")
    print("- K=8: Better specificity, might miss some distant homologs") 
    print("- K=9: High specificity, best for close matches only")
    print("\nK=5 is too short and will generate many false positive seed matches,")
    print("overwhelming the Smith-Waterman extension step with spurious alignments.")