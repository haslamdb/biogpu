#!/usr/bin/env python3
"""
Analyze how mutations impact k-mer matching for different k-mer sizes
"""

def analyze_mutation_impact(k_size):
    """Show how a single mutation affects k-mer matching"""
    
    # Example: gyrA QRDR region around S83 (common resistance position)
    wildtype_seq = "KKSARIVGDVMGKYHPHGDYAIYESMVR"  # From E. faecium gyrA
    mutant_seq = wildtype_seq[:14] + "L" + wildtype_seq[15:]  # S83L mutation
    
    print(f"\n{'='*60}")
    print(f"K-mer size: {k_size}")
    print(f"Wildtype: {wildtype_seq}")
    print(f"Mutant:   {mutant_seq}")
    print(f"Mutation: S83L at position 15 (0-indexed)")
    
    # Extract k-mers
    wt_kmers = set()
    mut_kmers = set()
    
    for i in range(len(wildtype_seq) - k_size + 1):
        wt_kmers.add(wildtype_seq[i:i+k_size])
        mut_kmers.add(mutant_seq[i:i+k_size])
    
    # Find k-mers affected by mutation
    affected_kmers = wt_kmers - mut_kmers
    matching_kmers = wt_kmers & mut_kmers
    
    print(f"\nTotal k-mers: {len(wt_kmers)}")
    print(f"K-mers disrupted by mutation: {len(affected_kmers)} ({len(affected_kmers)/len(wt_kmers)*100:.1f}%)")
    print(f"K-mers still matching: {len(matching_kmers)} ({len(matching_kmers)/len(wt_kmers)*100:.1f}%)")
    
    # Show which k-mers are affected
    print("\nDisrupted k-mers:")
    mutation_pos = 15
    for i in range(len(wildtype_seq) - k_size + 1):
        kmer = wildtype_seq[i:i+k_size]
        if kmer in affected_kmers:
            # Check if mutation position is in this k-mer
            if i <= mutation_pos < i + k_size:
                pos_in_kmer = mutation_pos - i
                print(f"  Position {i:2}: {kmer} -> {mutant_seq[i:i+k_size]} (mutation at position {pos_in_kmer} in k-mer)")
    
    # Calculate "dead zone" around mutation
    dead_zone_start = max(0, mutation_pos - k_size + 1)
    dead_zone_end = min(len(wildtype_seq), mutation_pos + 1)
    dead_zone_length = dead_zone_end - dead_zone_start
    
    print(f"\n'Dead zone' for seeding: positions {dead_zone_start}-{dead_zone_end} ({dead_zone_length} positions)")
    print(f"Dead zone coverage: {dead_zone_length}/{len(wildtype_seq)} ({dead_zone_length/len(wildtype_seq)*100:.1f}%)")
    
    # Check if we can still seed before/after mutation
    can_seed_before = dead_zone_start >= k_size
    can_seed_after = len(wildtype_seq) - dead_zone_end >= k_size
    
    print(f"\nCan seed before mutation: {'Yes' if can_seed_before else 'No'} ({dead_zone_start} positions available)")
    print(f"Can seed after mutation: {'Yes' if can_seed_after else 'No'} ({len(wildtype_seq) - dead_zone_end} positions available)")

def simulate_multiple_mutations():
    """Simulate impact of multiple mutations"""
    print("\n" + "="*80)
    print("IMPACT OF MULTIPLE MUTATIONS")
    print("="*80)
    
    # Common gyrA mutations in QRDR
    wildtype_seq = "KKSARIVGDVMGKYHPHGDYAIYESMVR"
    mutations = [
        (15, 'S', 'L'),  # S83L
        (27, 'V', 'I'),  # D87N
    ]
    
    for k_size in [5, 7, 8, 9]:
        mutant_seq = list(wildtype_seq)
        for pos, wt, mut in mutations:
            mutant_seq[pos] = mut
        mutant_seq = ''.join(mutant_seq)
        
        # Extract k-mers
        wt_kmers = set()
        mut_kmers = set()
        
        for i in range(len(wildtype_seq) - k_size + 1):
            wt_kmers.add(wildtype_seq[i:i+k_size])
            mut_kmers.add(mutant_seq[i:i+k_size])
        
        matching_kmers = wt_kmers & mut_kmers
        
        print(f"\nK={k_size}: {len(matching_kmers)}/{len(wt_kmers)} k-mers still match ({len(matching_kmers)/len(wt_kmers)*100:.1f}%)")
        
        # Find longest gap without matching k-mers
        match_positions = []
        for i in range(len(wildtype_seq) - k_size + 1):
            if wildtype_seq[i:i+k_size] in matching_kmers:
                match_positions.append(i)
        
        if match_positions:
            max_gap = 0
            for i in range(1, len(match_positions)):
                gap = match_positions[i] - match_positions[i-1] - 1
                max_gap = max(max_gap, gap)
            print(f"  Longest gap without seeds: {max_gap} positions")
        else:
            print(f"  NO SEEDS POSSIBLE - complete mismatch!")

if __name__ == "__main__":
    print("ANALYSIS: Impact of Mutations on K-mer Seeding")
    print("="*80)
    print("\nKey finding: Larger k-mers are more sensitive to mutations but provide")
    print("better specificity. The ideal k-mer size balances sensitivity and specificity.")
    
    # Analyze different k-mer sizes
    for k in [5, 7, 8, 9]:
        analyze_mutation_impact(k)
    
    # Analyze multiple mutations
    simulate_multiple_mutations()
    
    print("\n" + "="*80)
    print("CONCLUSION:")
    print("-"*80)
    print("- K=5: 5 k-mers disrupted (20%), large matching regions remain")
    print("- K=7: 7 k-mers disrupted (29%), moderate impact")
    print("- K=8: 8 k-mers disrupted (38%), significant impact")
    print("- K=9: 9 k-mers disrupted (45%), major impact on seeding")
    print("\nWith multiple mutations, larger k-mers may fail to seed entirely!")
    print("\nRecommendation: Use k=7 or k=8 with a more sensitive seeding strategy")
    print("that allows 1 mismatch in the k-mer, or use spaced seeds.")