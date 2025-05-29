#!/usr/bin/env python3
"""
Build GPU-optimized index for fluoroquinolone resistance mutation detection
Uses existing mutation_positions.json and reference_sequences.fasta
"""

import json
import numpy as np
from pathlib import Path
from Bio import SeqIO
from typing import Dict, List, Tuple
import logging
import h5py

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class FQMutationIndexBuilder:
    def __init__(self, resistance_db_dir: str):
        self.db_dir = Path(resistance_db_dir)
        self.mutation_data = {}
        self.reference_sequences = {}
        self.qrdr_regions = {}
        
    def load_mutation_positions(self):
        """Load mutation positions from JSON"""
        with open(self.db_dir / "mutation_positions.json", 'r') as f:
            self.mutation_data = json.load(f)
        logger.info(f"Loaded {len(self.mutation_data)} mutation records")
        
    def load_reference_sequences(self):
        """Load reference sequences from FASTA"""
        fasta_path = self.db_dir / "reference_sequences.fasta"
        for record in SeqIO.parse(fasta_path, "fasta"):
            # The description contains the full ID we need
            # Format: "Species species_gene_accession Species species gene"
            # The ID we need is everything before the second occurrence of the species name
            desc_parts = record.description.split()
            
            # Find where the gene and accession info ends
            # This is typically at the third underscore-containing element
            id_parts = []
            for part in desc_parts:
                id_parts.append(part)
                if '_' in part and part.count('_') >= 2:
                    # This is the accession part, we're done
                    break
            
            seq_id = ' '.join(id_parts)
            self.reference_sequences[seq_id] = str(record.seq)
            
        logger.info(f"Loaded {len(self.reference_sequences)} reference sequences")
        
        # Debug: print first few IDs to verify matching
        ref_ids = list(self.reference_sequences.keys())[:3]
        mut_ids = list(self.mutation_data.keys())[:3]
        logger.info(f"Sample sequence IDs: {ref_ids}")
        logger.info(f"Sample mutation IDs: {mut_ids}")
        
        # Check for matches
        matches = sum(1 for k in self.mutation_data if k in self.reference_sequences)
        logger.info(f"Matching IDs between mutations and sequences: {matches}/{len(self.mutation_data)}")
        
    def extract_qrdr_regions(self):
        """Extract QRDR regions with mutations for GPU alignment"""
        
        # Define standard QRDR boundaries (in amino acids)
        qrdr_definitions = {
            'gyrA': {'start': 67, 'end': 106},
            'gyrB': {'start': 426, 'end': 466},
            'parC': {'start': 64, 'end': 102},
            'parE': {'start': 410, 'end': 470},
            'grlA': {'start': 64, 'end': 102}  # S. aureus parC equivalent
        }
        
        for seq_id, mut_info in self.mutation_data.items():
            gene = mut_info['gene']
            
            if gene not in qrdr_definitions:
                continue
                
            # Get QRDR boundaries
            qrdr = qrdr_definitions[gene]
            qrdr_nt_start = (qrdr['start'] - 1) * 3
            qrdr_nt_end = qrdr['end'] * 3
            
            # Add flanking regions for alignment
            flank = 150
            region_start = max(0, qrdr_nt_start - flank)
            region_end = qrdr_nt_end + flank
            
            # Extract sequence region
            if seq_id in self.reference_sequences:
                full_seq = self.reference_sequences[seq_id]
                qrdr_seq = full_seq[region_start:min(len(full_seq), region_end)]
                
                # Calculate mutation positions relative to QRDR region
                relative_mutations = []
                for mutation in mut_info['mutations']:
                    # Get mutation position relative to gene start
                    mut_nt_pos = mutation['nucleotide_position']
                    gene_start = mut_info['gene_start']
                    mut_pos_in_gene = mut_nt_pos - gene_start
                    
                    # Position relative to QRDR region
                    mut_pos_in_qrdr = mut_pos_in_gene - region_start
                    
                    if 0 <= mut_pos_in_qrdr < len(qrdr_seq):
                        relative_mutations.append({
                            'position': mut_pos_in_qrdr,
                            'codon_position': mutation['aa_position'],
                            'wild_type': mutation['wild_type_aa'],
                            'mutant': mutation['mutant_aa'],
                            'absolute_position': mut_nt_pos
                        })
                
                self.qrdr_regions[seq_id] = {
                    'species': mut_info['species'],
                    'gene': gene,
                    'sequence': qrdr_seq,
                    'region_start': region_start,
                    'region_end': region_end,
                    'mutations': relative_mutations,
                    'is_mutant': True  # These are mutant reference sequences
                }
                
        logger.info(f"Extracted {len(self.qrdr_regions)} QRDR regions")
        
    def create_gpu_structures(self):
        """Create GPU-optimized data structures"""
        
        # Encode sequences as 2-bit integers
        base_encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}
        
        # Group by gene for efficient batch processing
        genes = {}
        for seq_id, region in self.qrdr_regions.items():
            gene = region['gene']
            if gene not in genes:
                genes[gene] = []
            genes[gene].append((seq_id, region))
        
        gpu_data = {}
        
        for gene, regions in genes.items():
            # Find max sequence length for this gene
            max_len = max(len(r[1]['sequence']) for r in regions)
            
            # Create arrays
            n_seqs = len(regions)
            sequences = np.zeros((n_seqs, max_len), dtype=np.uint8)
            seq_lengths = np.zeros(n_seqs, dtype=np.int32)
            mutation_masks = np.zeros((n_seqs, max_len), dtype=np.uint8)
            position_weights = np.ones((n_seqs, max_len), dtype=np.float32)  # Position-specific weights
            
            # Metadata
            seq_ids = []
            species_list = []
            mutation_info = []
            
            for i, (seq_id, region) in enumerate(regions):
                # Encode sequence
                seq = region['sequence']
                seq_lengths[i] = len(seq)
                
                for j, base in enumerate(seq):
                    sequences[i, j] = base_encoding.get(base, 0)
                
                # Mark mutation positions and create position weights
                for mut in region['mutations']:
                    pos = mut['position']
                    if pos < max_len:
                        # Mark the entire codon
                        for k in range(3):
                            if pos + k < max_len:
                                mutation_masks[i, pos + k] = 1
                        
                        # Create position-specific weights
                        # Highest weight at mutation site, decreasing with distance
                        for j in range(max_len):
                            distance = abs(j - pos)
                            if distance == 0:
                                position_weights[i, j] = 5.0  # Maximum weight at mutation
                            elif distance <= 3:
                                position_weights[i, j] = 3.0  # High weight for codon
                            elif distance <= 15:
                                position_weights[i, j] = 2.0  # Medium weight nearby
                            elif distance <= 50:
                                position_weights[i, j] = 1.5  # Slight elevation in flanking
                            # Default weight is 1.0
                
                # Store metadata
                seq_ids.append(seq_id)
                species_list.append(region['species'])
                mutation_info.append(region['mutations'])
            
            gpu_data[gene] = {
                'sequences': sequences,
                'lengths': seq_lengths,
                'mutation_masks': mutation_masks,
                'position_weights': position_weights,
                'seq_ids': seq_ids,
                'species': species_list,
                'mutations': mutation_info,
                'max_length': max_len
            }
            
        return gpu_data
    
    def create_alignment_parameters(self):
        """Create parameters for GPU-accelerated Smith-Waterman alignment"""
        
        # DNA scoring matrix (4x4)
        # Match on diagonal, mismatch off-diagonal
        match_score = 2.0
        mismatch_penalty = -3.0
        
        scoring_matrix = np.full((4, 4), mismatch_penalty, dtype=np.float32)
        np.fill_diagonal(scoring_matrix, match_score)
        
        # Gap penalties
        gap_open = -5.0
        gap_extend = -2.0
        
        # Mutation-aware scoring adjustments
        # These multiply the base scores at mutation positions
        mutation_match_bonus = 2.0      # Bonus for matching at mutation site
        mutation_mismatch_penalty = 2.0  # Extra penalty for mismatch at mutation
        
        alignment_params = {
            'scoring_matrix': scoring_matrix,
            'gap_open': gap_open,
            'gap_extend': gap_extend,
            'mutation_match_bonus': mutation_match_bonus,
            'mutation_mismatch_penalty': mutation_mismatch_penalty,
            'min_alignment_score': 50.0,  # Minimum score to consider valid
            'min_identity': 0.85  # Minimum sequence identity
        }
        
        return alignment_params
    
    def create_kmer_index(self, k=15):
        """Create k-mer index for fast pre-filtering"""
        
        kmer_index = {}
        
        for seq_id, region in self.qrdr_regions.items():
            seq = region['sequence']
            gene = region['gene']
            
            # Extract k-mers around mutation sites
            for mut in region['mutations']:
                pos = mut['position']
                
                # Get k-mers covering the mutation
                for offset in range(-k+1, 3):  # Cover the codon
                    kmer_start = pos + offset
                    kmer_end = kmer_start + k
                    
                    if 0 <= kmer_start and kmer_end <= len(seq):
                        kmer = seq[kmer_start:kmer_end]
                        if 'N' not in kmer:
                            if kmer not in kmer_index:
                                kmer_index[kmer] = []
                            
                            kmer_index[kmer].append({
                                'seq_id': seq_id,
                                'gene': gene,
                                'mutation_pos': pos - kmer_start,
                                'species': region['species']
                            })
        
        # Convert to arrays for GPU
        unique_kmers = sorted(kmer_index.keys())
        kmer_array = np.array([list(k) for k in unique_kmers], dtype='S1')
        
        # Create lookup structures
        kmer_to_idx = {k: i for i, k in enumerate(unique_kmers)}
        
        return {
            'kmers': kmer_array,
            'kmer_to_idx': kmer_to_idx,
            'kmer_index': kmer_index,
            'k': k
        }
    
    def save_index(self, output_path: str):
        """Save the complete index"""
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create GPU structures
        gpu_data = self.create_gpu_structures()
        
        # Create k-mer index
        kmer_data = self.create_kmer_index()
        
        # Create alignment parameters
        alignment_params = self.create_alignment_parameters()
        
        # Save as HDF5 for efficient loading
        with h5py.File(output_dir / "fq_mutation_index.h5", 'w') as f:
            # Save GPU data by gene
            for gene, data in gpu_data.items():
                gene_group = f.create_group(gene)
                gene_group.create_dataset('sequences', data=data['sequences'], compression='gzip')
                gene_group.create_dataset('lengths', data=data['lengths'])
                gene_group.create_dataset('mutation_masks', data=data['mutation_masks'], compression='gzip')
                gene_group.create_dataset('position_weights', data=data['position_weights'], compression='gzip')
                gene_group.attrs['max_length'] = data['max_length']
                
                # Save metadata as JSON
                metadata = {
                    'seq_ids': data['seq_ids'],
                    'species': data['species'],
                    'mutations': data['mutations']
                }
                gene_group.attrs['metadata'] = json.dumps(metadata)
            
            # Save k-mer index
            kmer_group = f.create_group('kmer_index')
            kmer_group.create_dataset('kmers', data=kmer_data['kmers'])
            kmer_group.attrs['k'] = kmer_data['k']
            
            # Save alignment parameters
            align_group = f.create_group('alignment')
            align_group.create_dataset('scoring_matrix', data=alignment_params['scoring_matrix'])
            align_group.attrs['gap_open'] = alignment_params['gap_open']
            align_group.attrs['gap_extend'] = alignment_params['gap_extend']
            align_group.attrs['mutation_match_bonus'] = alignment_params['mutation_match_bonus']
            align_group.attrs['mutation_mismatch_penalty'] = alignment_params['mutation_mismatch_penalty']
            align_group.attrs['min_alignment_score'] = alignment_params['min_alignment_score']
            align_group.attrs['min_identity'] = alignment_params['min_identity']
            
        # Save additional metadata
        with open(output_dir / "kmer_lookup.json", 'w') as f:
            json.dump({
                'kmer_to_idx': kmer_data['kmer_to_idx'],
                'kmer_index': kmer_data['kmer_index']
            }, f)
        
        # Save alignment parameters separately for easy access
        with open(output_dir / "alignment_params.json", 'w') as f:
            # Convert numpy array to list for JSON serialization
            params_copy = alignment_params.copy()
            params_copy['scoring_matrix'] = alignment_params['scoring_matrix'].tolist()
            json.dump(params_copy, f, indent=2)
        
        # Save summary
        summary = {
            'total_sequences': len(self.qrdr_regions),
            'genes': list(gpu_data.keys()),
            'sequences_per_gene': {g: len(d['seq_ids']) for g, d in gpu_data.items()},
            'total_kmers': len(kmer_data['kmers']),
            'kmer_size': kmer_data['k'],
            'alignment_configured': True,
            'position_weights_included': True
        }
        
        with open(output_dir / "index_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Index saved to {output_dir}")
        logger.info(f"Total sequences: {summary['total_sequences']}")
        logger.info(f"Total k-mers indexed: {summary['total_kmers']}")
        logger.info("Position-specific weights and alignment parameters included")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Build FQ mutation detection index')
    parser.add_argument('--db-dir', type=str, default='data/resistance_db',
                       help='Path to resistance database directory')
    parser.add_argument('--output', type=str, default='data/indices/fq_mutations',
                       help='Output directory for index')
    
    args = parser.parse_args()
    
    # Build index
    builder = FQMutationIndexBuilder(args.db_dir)
    
    logger.info("Loading mutation positions...")
    builder.load_mutation_positions()
    
    logger.info("Loading reference sequences...")
    builder.load_reference_sequences()
    
    logger.info("Extracting QRDR regions...")
    builder.extract_qrdr_regions()
    
    logger.info("Building and saving index...")
    builder.save_index(args.output)
    
    logger.info("Index building complete!")

if __name__ == "__main__":
    main()