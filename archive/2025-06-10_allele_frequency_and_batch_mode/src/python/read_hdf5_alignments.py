#!/usr/bin/env python3
"""
Read and analyze HDF5 alignment output from BioGPU pipeline
"""

import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

class BioGPUAlignmentReader:
    def __init__(self, hdf5_path):
        self.h5_file = h5py.File(hdf5_path, 'r')
        self.alignments = None
        self.metadata = {}
        self._load_metadata()
        
    def _load_metadata(self):
        """Load metadata from HDF5 file"""
        meta_group = self.h5_file['metadata']
        for key in meta_group.attrs:
            self.metadata[key] = meta_group.attrs[key]
            
    def load_alignments(self):
        """Load alignment data into pandas DataFrame"""
        align_group = self.h5_file['alignments']
        
        # Load basic alignment data
        data = {
            'read_id': align_group['read_id'][:],
            'gene_id': align_group['gene_id'][:],
            'species_id': align_group['species_id'][:],
            'seq_id': align_group['seq_id'][:],
            'alignment_score': align_group['alignment_score'][:],
            'identity': align_group['identity'][:],
            'start_pos': align_group['start_pos'][:],
            'matches': align_group['matches'][:],
            'num_mutations': align_group['num_mutations'][:]
        }
        
        self.alignments = pd.DataFrame(data)
        return self.alignments
    
    def get_mutation_positions(self, alignment_idx):
        """Get mutation positions for a specific alignment"""
        offsets = self.h5_file['alignments/mutation_offsets'][:]
        positions = self.h5_file['alignments/mutation_positions'][:]
        
        start = offsets[alignment_idx]
        end = offsets[alignment_idx + 1] if alignment_idx + 1 < len(offsets) else len(positions)
        
        return positions[start:end]
    
    def summarize_by_gene(self):
        """Summarize alignments by gene"""
        if self.alignments is None:
            self.load_alignments()
            
        summary = self.alignments.groupby('gene_id').agg({
            'read_id': 'count',
            'alignment_score': ['mean', 'std'],
            'identity': ['mean', 'std'],
            'num_mutations': 'sum'
        })
        
        summary.columns = ['_'.join(col) for col in summary.columns]
        summary.rename(columns={'read_id_count': 'num_alignments'}, inplace=True)
        
        return summary
    
    def plot_alignment_distribution(self, output_dir=None):
        """Create visualization plots"""
        if self.alignments is None:
            self.load_alignments()
            
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Alignment score distribution
        axes[0, 0].hist(self.alignments['alignment_score'], bins=50, alpha=0.7)
        axes[0, 0].set_xlabel('Alignment Score')
        axes[0, 0].set_ylabel('Count')
        axes[0, 0].set_title('Distribution of Alignment Scores')
        
        # Plot 2: Identity distribution
        axes[0, 1].hist(self.alignments['identity'], bins=50, alpha=0.7, color='green')
        axes[0, 1].set_xlabel('Identity (%)')
        axes[0, 1].set_ylabel('Count')
        axes[0, 1].set_title('Distribution of Sequence Identity')
        
        # Plot 3: Alignments per gene
        gene_counts = self.alignments['gene_id'].value_counts().head(20)
        axes[1, 0].bar(range(len(gene_counts)), gene_counts.values)
        axes[1, 0].set_xlabel('Gene ID')
        axes[1, 0].set_ylabel('Number of Alignments')
        axes[1, 0].set_title('Top 20 Genes by Alignment Count')
        
        # Plot 4: Mutations detected
        mutation_counts = self.alignments['num_mutations'].value_counts().sort_index()
        axes[1, 1].bar(mutation_counts.index, mutation_counts.values, color='red')
        axes[1, 1].set_xlabel('Number of Mutations')
        axes[1, 1].set_ylabel('Number of Reads')
        axes[1, 1].set_title('Distribution of Mutations per Read')
        
        plt.tight_layout()
        
        if output_dir:
            output_path = Path(output_dir) / 'alignment_summary.png'
            plt.savefig(output_path, dpi=300)
            print(f"Saved plot to {output_path}")
        else:
            plt.show()
    
    def export_for_ml(self, output_path):
        """Export data in ML-ready format"""
        if self.alignments is None:
            self.load_alignments()
            
        # Create feature matrix
        features = self.alignments[['alignment_score', 'identity', 'matches', 'num_mutations']].values
        
        # Create labels (example: high confidence if identity > 0.95 and score > 50)
        labels = ((self.alignments['identity'] > 0.95) & 
                 (self.alignments['alignment_score'] > 50)).astype(int)
        
        # Save as numpy arrays
        np.savez(output_path,
                 features=features,
                 labels=labels,
                 gene_ids=self.alignments['gene_id'].values,
                 species_ids=self.alignments['species_id'].values)
        
        print(f"Exported ML data to {output_path}")
        print(f"Feature shape: {features.shape}")
        print(f"Positive labels: {labels.sum()} / {len(labels)}")

def main():
    parser = argparse.ArgumentParser(description='Read and analyze BioGPU HDF5 output')
    parser.add_argument('hdf5_file', help='Path to HDF5 alignment file')
    parser.add_argument('--summarize', action='store_true', help='Print summary statistics')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    parser.add_argument('--export-ml', help='Export ML-ready features to NPZ file')
    parser.add_argument('--output-dir', help='Directory for output files')
    
    args = parser.parse_args()
    
    reader = BioGPUAlignmentReader(args.hdf5_file)
    
    if args.summarize:
        print("\n=== Alignment Summary ===")
        print(f"Total alignments: {reader.metadata.get('total_alignments', 'N/A')}")
        print(f"Analysis timestamp: {reader.metadata.get('analysis_timestamp', 'N/A')}")
        print(f"Pipeline version: {reader.metadata.get('pipeline_version', 'N/A')}")
        
        alignments = reader.load_alignments()
        print(f"\nLoaded {len(alignments)} alignments")
        
        print("\n=== Summary by Gene ===")
        gene_summary = reader.summarize_by_gene()
        print(gene_summary.head(10))
    
    if args.plot:
        reader.plot_alignment_distribution(args.output_dir)
    
    if args.export_ml:
        reader.export_for_ml(args.export_ml)

if __name__ == '__main__':
    main()