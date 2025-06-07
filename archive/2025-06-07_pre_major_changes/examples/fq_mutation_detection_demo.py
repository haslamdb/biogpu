#!/usr/bin/env python3
"""
Demonstration of two-stage fluoroquinolone mutation detection:
1. K-mer pre-filtering to find candidate reads
2. Exact alignment with position-specific weights around mutation sites
"""

import numpy as np
import h5py
import json
from pathlib import Path
from typing import List, Tuple, Dict

class FQMutationDetector:
    def __init__(self, index_dir: str):
        self.index_dir = Path(index_dir)
        self.load_index()
        
    def load_index(self):
        """Load the pre-built index"""
        # Load k-mer lookup
        with open(self.index_dir / "kmer_lookup.json", 'r') as f:
            kmer_data = json.load(f)
            self.kmer_lookup = kmer_data['kmer_index']
            self.kmer_to_idx = kmer_data['kmer_to_idx']
        
        # Load alignment parameters
        with open(self.index_dir / "alignment_params.json", 'r') as f:
            self.align_params = json.load(f)
            self.scoring_matrix = np.array(self.align_params['scoring_matrix'])
        
        # Load HDF5 data
        self.h5_file = h5py.File(self.index_dir / "fq_mutation_index.h5", 'r')
        
        print("Index loaded successfully")
        
    def stage1_kmer_filter(self, read_seq: str, k: int = 15) -> List[Dict]:
        """
        Stage 1: Fast k-mer pre-filtering
        Returns list of potential matches with gene/species info
        """
        candidates = set()
        
        # Extract k-mers from read
        for i in range(len(read_seq) - k + 1):
            kmer = read_seq[i:i + k]
            
            # Check if k-mer exists in our index
            if kmer in self.kmer_lookup:
                matches = self.kmer_lookup[kmer]
                for match in matches:
                    # Add unique (gene, species, seq_id) tuples
                    candidates.add((
                        match['gene'],
                        match['species'],
                        match['seq_id']
                    ))
        
        # Convert to list of dicts
        candidate_list = []
        for gene, species, seq_id in candidates:
            candidate_list.append({
                'gene': gene,
                'species': species,
                'seq_id': seq_id
            })
        
        return candidate_list
    
    def stage2_exact_alignment(self, read_seq: str, candidate: Dict) -> Dict:
        """
        Stage 2: Exact alignment with position-specific weights
        Uses Smith-Waterman with enhanced scoring at mutation sites
        """
        gene = candidate['gene']
        seq_id = candidate['seq_id']
        
        # Get reference data from HDF5
        gene_group = self.h5_file[gene]
        metadata = json.loads(gene_group.attrs['metadata'])
        
        # Find the sequence index
        seq_idx = metadata['seq_ids'].index(seq_id)
        
        # Get reference sequence and weights
        ref_seq_encoded = gene_group['sequences'][seq_idx]
        seq_length = gene_group['lengths'][seq_idx]
        position_weights = gene_group['position_weights'][seq_idx]
        mutation_mask = gene_group['mutation_masks'][seq_idx]
        
        # Convert reference to string for demo
        base_decode = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        ref_seq = ''.join([base_decode[b] for b in ref_seq_encoded[:seq_length]])
        
        # Get mutation info
        mutations = metadata['mutations'][seq_idx]
        
        # Simplified alignment scoring (full GPU version would use dynamic programming)
        alignment_score = 0
        matches = 0
        
        # For demo, we'll do a simple sliding window alignment
        best_score = -1
        best_pos = -1
        best_alignment = {}
        
        for start_pos in range(max(1, len(ref_seq) - len(read_seq) - 50), 
                              min(len(ref_seq) - len(read_seq) + 50, len(ref_seq))):
            if start_pos < 0:
                continue
                
            score = 0
            local_matches = 0
            
            for i, read_base in enumerate(read_seq):
                if start_pos + i >= len(ref_seq):
                    break
                    
                ref_base = ref_seq[start_pos + i]
                position_weight = position_weights[start_pos + i] if start_pos + i < len(position_weights) else 1.0
                
                if read_base == ref_base:
                    # Match - apply position weight
                    score += self.align_params['scoring_matrix'][0][0] * position_weight
                    local_matches += 1
                    
                    # Extra bonus if this is a mutation position and we match
                    if mutation_mask[start_pos + i] > 0:
                        score += self.align_params['mutation_match_bonus'] * position_weight
                else:
                    # Mismatch
                    score += self.align_params['scoring_matrix'][0][1] * position_weight
                    
                    # Extra penalty if this is a mutation position and we mismatch
                    if mutation_mask[start_pos + i] > 0:
                        score -= self.align_params['mutation_mismatch_penalty'] * position_weight
            
            # Track best alignment
            if score > best_score:
                best_score = score
                best_pos = start_pos
                best_alignment = {
                    'score': score,
                    'identity': local_matches / len(read_seq) if len(read_seq) > 0 else 0,
                    'start_pos': start_pos,
                    'matches': local_matches
                }
        
        # Check mutations covered
        mutations_detected = []
        if best_pos >= 0:
            for mut in mutations:
                mut_pos = mut['position']
                # Check if mutation position is covered by the read
                if best_pos <= mut_pos < best_pos + len(read_seq):
                    read_pos = mut_pos - best_pos
                    read_base = read_seq[read_pos]
                    ref_base = ref_seq[mut_pos]
                    
                    # Determine if read has mutant or wild-type
                    mutations_detected.append({
                        'position': mut['codon_position'],
                        'wild_type': mut['wild_type'],
                        'mutant': mut['mutant'],
                        'read_base': read_base,
                        'is_mutant': read_base == ref_base  # Since ref is mutant sequence
                    })
        
        result = {
            'gene': gene,
            'species': candidate['species'],
            'seq_id': seq_id,
            'alignment': best_alignment,
            'mutations_detected': mutations_detected,
            'passes_threshold': (best_score >= self.align_params['min_alignment_score'] and
                               best_alignment.get('identity', 0) >= self.align_params['min_identity'])
        }
        
        return result
    
    def detect_mutations(self, read_seq: str) -> List[Dict]:
        """
        Complete two-stage mutation detection pipeline
        """
        # Stage 1: K-mer filtering
        candidates = self.stage1_kmer_filter(read_seq)
        
        if not candidates:
            return []
        
        # Stage 2: Exact alignment for each candidate
        results = []
        for candidate in candidates:
            alignment_result = self.stage2_exact_alignment(read_seq, candidate)
            if alignment_result['passes_threshold']:
                results.append(alignment_result)
        
        return results
    
    def close(self):
        """Clean up resources"""
        self.h5_file.close()

def demonstrate():
    """Demonstrate the two-stage detection process"""
    
    # Initialize detector
    detector = FQMutationDetector("data/indices/fq_mutations")
    
    # Example 1: Simulated read with gyrA S83L mutation (TCG -> TTG)
    print("\n=== Example 1: Read with gyrA S83L mutation ===")
    # This would be around codon 83 region of gyrA
    mutant_read = "ATGGCGACGAATATCCCGCCGCACAACCTGACGGAAGTGATTAACGGCTGCCTGGCGTATATCGACAACGAAGACATCAGCATTGAAGGGCTGATGGAACA"
    
    results = detector.detect_mutations(mutant_read)
    
    for result in results:
        print(f"\nGene: {result['gene']}")
        print(f"Species: {result['species']}")
        print(f"Alignment Score: {result['alignment']['score']:.2f}")
        print(f"Identity: {result['alignment']['identity']:.2%}")
        
        if result['mutations_detected']:
            print("Mutations detected:")
            for mut in result['mutations_detected']:
                status = "MUTANT" if mut['is_mutant'] else "WILD-TYPE"
                print(f"  Position {mut['position']}: {mut['wild_type']} -> {mut['mutant']} [{status}]")
    
    # Example 2: Wild-type read (no mutations)
    print("\n=== Example 2: Wild-type read ===")
    wt_read = "ATGGCTATTATGTCGCGCCGGAACAGCCCGGATTGGGTCAGGAATTAAACGATGAGGTAGTGAAAGAATACCTGGCCTATGTGATTAAATAG"
    
    results = detector.detect_mutations(wt_read)
    
    if not results:
        print("No significant matches found (likely wild-type or non-target region)")
    else:
        for result in results:
            print(f"\nGene: {result['gene']}")
            print(f"Species: {result['species']}")
            if not result['mutations_detected']:
                print("No mutations detected in this region")
    
    # Show the position weighting effect
    print("\n=== Position Weight Demonstration ===")
    print("Position weights increase scoring sensitivity around mutation sites:")
    print("- At mutation site: 5.0x weight")
    print("- Within codon (±3bp): 3.0x weight")
    print("- Nearby (±15bp): 2.0x weight")
    print("- Flanking (±50bp): 1.5x weight")
    print("- Other positions: 1.0x weight")
    
    detector.close()

if __name__ == "__main__":
    demonstrate()