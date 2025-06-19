#!/usr/bin/env python3
"""
Comprehensive Metagenomic Analysis Toolkit
==========================================

Python toolkit for comprehensive metagenomic profiling and analysis
Integrates with hybrid CPU-GPU profiler and provides rich analysis capabilities

Author: Metagenomic Analysis Pipeline
Usage: For comprehensive microbiome analysis and visualization
"""

import subprocess
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import sys
from typing import Dict, List, Tuple, Optional, Union
import logging
from scipy import stats
from scipy.spatial.distance import braycurtis, jaccard
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class ComprehensiveMetagenomicsAnalyzer:
    """
    Comprehensive metagenomic analysis and visualization toolkit
    Handles profiling, statistical analysis, and rich visualizations
    """
    
    def __init__(self, 
                 profiler_binary: str = "./build/hybrid_profiler",
                 database_path: str = None):
        self.profiler_binary = Path(profiler_binary)
        self.database_path = Path(database_path) if database_path else None
        
        # Analysis results storage
        self.abundance_data = None
        self.taxonomy_data = None
        self.coverage_data = None
        self.diversity_metrics = {}
        
        # Setup logging
        logging.basicConfig(level=logging.INFO,
                          format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
    
    def validate_inputs(self, fastq_file: str) -> bool:
        """Validate input files and system requirements"""
        
        if not self.profiler_binary.exists():
            self.logger.error(f"Profiler binary not found: {self.profiler_binary}")
            return False
        
        if not self.database_path or not self.database_path.exists():
            self.logger.error(f"Database not found: {self.database_path}")
            return False
        
        fastq_path = Path(fastq_file)
        if not fastq_path.exists():
            self.logger.error(f"FASTQ file not found: {fastq_file}")
            return False
        
        return True
    
    def run_comprehensive_profiling(self, 
                                  fastq_file: str, 
                                  output_prefix: str,
                                  **kwargs) -> Dict:
        """Run comprehensive metagenomic profiling"""
        
        self.logger.info(f"Starting comprehensive profiling of {fastq_file}")
        
        # Prepare command
        cmd = [
            str(self.profiler_binary),
            str(self.database_path),
            fastq_file,
            output_prefix
        ]
        
        # Add optional parameters
        if 'memory_limit' in kwargs:
            cmd.extend(['--memory', str(kwargs['memory_limit'])])
        
        # Run profiler
        try:
            self.logger.info("Running hybrid CPU-GPU profiler...")
            result = subprocess.run(cmd, 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=7200)  # 2 hour timeout
            
            if result.returncode != 0:
                self.logger.error(f"Profiler failed: {result.stderr}")
                raise RuntimeError(f"Profiling failed: {result.stderr}")
            
            self.logger.info("Profiling completed successfully")
            
            # Load and parse results
            return self.load_profiling_results(output_prefix)
            
        except subprocess.TimeoutExpired:
            self.logger.error("Profiling timed out after 2 hours")
            raise RuntimeError("Profiling timeout")
    
    def load_profiling_results(self, output_prefix: str) -> Dict:
        """Load and parse profiler output files"""
        
        results = {
            'abundance_table': None,
            'taxonomy_summary': None,
            'coverage_stats': None,
            'organism_report': None,
            'kraken_style': None
        }
        
        # Load abundance table
        abundance_file = f"{output_prefix}_abundance_table.tsv"
        if Path(abundance_file).exists():
            self.abundance_data = pd.read_csv(abundance_file, sep='\t')
            results['abundance_table'] = self.abundance_data
            self.logger.info(f"Loaded {len(self.abundance_data)} organism abundances")
        
        # Load taxonomy summary
        taxonomy_file = f"{output_prefix}_taxonomy_summary.tsv"
        if Path(taxonomy_file).exists():
            self.taxonomy_data = pd.read_csv(taxonomy_file, sep='\t')
            results['taxonomy_summary'] = self.taxonomy_data
        
        # Load coverage statistics
        coverage_file = f"{output_prefix}_coverage_stats.tsv"
        if Path(coverage_file).exists():
            self.coverage_data = pd.read_csv(coverage_file, sep='\t')
            results['coverage_stats'] = self.coverage_data
        
        # Load organism report
        report_file = f"{output_prefix}_organism_report.txt"
        if Path(report_file).exists():
            with open(report_file, 'r') as f:
                results['organism_report'] = f.read()
        
        # Load Kraken-style output
        kraken_file = f"{output_prefix}_kraken_style.txt"
        if Path(kraken_file).exists():
            results['kraken_style'] = pd.read_csv(kraken_file, 
                                                sep='\t', 
                                                names=['percentage', 'clade_reads', 'taxon_reads', 
                                                      'rank', 'taxid', 'name'])
        
        return results
    
    def calculate_diversity_metrics(self, abundance_column: str = 'relative_abundance') -> Dict:
        """Calculate comprehensive diversity metrics"""
        
        if self.abundance_data is None:
            self.logger.error("No abundance data loaded")
            return {}
        
        abundances = self.abundance_data[abundance_column].values
        abundances = abundances[abundances > 0]  # Remove zeros
        
        # Alpha diversity metrics
        metrics = {}
        
        # Species richness (number of observed species)
        metrics['richness'] = len(abundances)
        
        # Shannon diversity
        metrics['shannon'] = -np.sum(abundances * np.log(abundances))
        
        # Simpson diversity
        metrics['simpson'] = 1 - np.sum(abundances ** 2)
        
        # Pielou's evenness
        metrics['evenness'] = metrics['shannon'] / np.log(metrics['richness']) if metrics['richness'] > 1 else 0
        
        # Berger-Parker dominance
        metrics['dominance'] = np.max(abundances)
        
        # Chao1 richness estimator (simplified version)
        singletons = np.sum(abundances < 0.001)  # Very low abundance as proxy for singletons
        doubletons = np.sum((abundances >= 0.001) & (abundances < 0.002))
        if doubletons > 0:
            metrics['chao1'] = metrics['richness'] + (singletons ** 2) / (2 * doubletons)
        else:
            metrics['chao1'] = metrics['richness']
        
        self.diversity_metrics = metrics
        
        self.logger.info(f"Diversity metrics calculated:")
        for metric, value in metrics.items():
            self.logger.info(f"  {metric}: {value:.3f}")
        
        return metrics
    
    def analyze_taxonomic_composition(self, level: str = 'Genus') -> pd.DataFrame:
        """Analyze composition at specific taxonomic level"""
        
        if self.taxonomy_data is None:
            self.logger.error("No taxonomy data loaded")
            return pd.DataFrame()
        
        # Filter by taxonomic level
        level_data = self.taxonomy_data[self.taxonomy_data['level'] == level].copy()
        
        if level_data.empty:
            self.logger.warning(f"No data found for taxonomic level: {level}")
            return pd.DataFrame()
        
        # Sort by abundance
        level_data = level_data.sort_values('relative_abundance', ascending=False)
        
        # Calculate cumulative abundance
        level_data['cumulative_abundance'] = level_data['relative_abundance'].cumsum()
        
        self.logger.info(f"Found {len(level_data)} taxa at {level} level")
        self.logger.info(f"Top 5 {level} taxa:")
        for _, row in level_data.head().iterrows():
            self.logger.info(f"  {row['taxon']}: {row['relative_abundance']:.1%}")
        
        return level_data
    
    def create_abundance_plots(self, output_dir: str = "plots", top_n: int = 20):
        """Create comprehensive abundance visualization plots"""
        
        if self.abundance_data is None:
            self.logger.error("No abundance data available for plotting")
            return
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Filter to top N most abundant organisms
        top_organisms = self.abundance_data.nlargest(top_n, 'relative_abundance')
        
        # 1. Bar plot of top organisms
        plt.figure(figsize=(15, 8))
        
        # Create horizontal bar plot
        y_pos = np.arange(len(top_organisms))
        bars = plt.barh(y_pos, top_organisms['relative_abundance'] * 100)
        
        # Color bars by abundance
        colors = plt.cm.viridis(top_organisms['relative_abundance'] / top_organisms['relative_abundance'].max())
        for bar, color in zip(bars, colors):
            bar.set_color(color)
        
        plt.yticks(y_pos, [name[:50] + '...' if len(name) > 50 else name 
                          for name in top_organisms['organism_name']])
        plt.xlabel('Relative Abundance (%)')
        plt.title(f'Top {top_n} Most Abundant Organisms')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(output_path / 'abundance_barplot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Pie chart of top organisms
        plt.figure(figsize=(12, 10))
        
        # Group small abundances into "Others"
        plot_data = top_organisms.head(10).copy()
        others_abundance = top_organisms.iloc[10:]['relative_abundance'].sum()
        
        if others_abundance > 0:
            others_row = pd.DataFrame({
                'organism_name': ['Others'],
                'relative_abundance': [others_abundance]
            })
            plot_data = pd.concat([plot_data, others_row], ignore_index=True)
        
        # Create pie chart
        wedges, texts, autotexts = plt.pie(plot_data['relative_abundance'] * 100, 
                                          labels=[name[:30] + '...' if len(name) > 30 else name 
                                                 for name in plot_data['organism_name']],
                                          autopct='%1.1f%%',
                                          startangle=90)
        
        plt.title('Taxonomic Composition (Top 10 + Others)')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig(output_path / 'abundance_piechart.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Coverage vs Abundance scatter plot
        if self.coverage_data is not None:
            plt.figure(figsize=(12, 8))
            
            # Merge abundance and coverage data
            merged_data = pd.merge(self.abundance_data, self.coverage_data, on='organism_id')
            
            # Create scatter plot
            scatter = plt.scatter(merged_data['coverage_breadth'] * 100,
                                merged_data['relative_abundance'] * 100,
                                s=merged_data['unique_kmers'] / 100,
                                alpha=0.6,
                                c=merged_data['coverage_depth'],
                                cmap='viridis')
            
            plt.xlabel('Coverage Breadth (%)')
            plt.ylabel('Relative Abundance (%)')
            plt.title('Coverage vs Abundance (Size = Unique k-mers, Color = Depth)')
            plt.colorbar(scatter, label='Coverage Depth')
            plt.tight_layout()
            plt.savefig(output_path / 'coverage_vs_abundance.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        self.logger.info(f"Abundance plots saved to {output_path}")
    
    def create_taxonomy_plots(self, output_dir: str = "plots"):
        """Create taxonomic composition visualizations"""
        
        if self.taxonomy_data is None:
            self.logger.error("No taxonomy data available for plotting")
            return
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Plot composition at different taxonomic levels
        levels = self.taxonomy_data['level'].unique()
        
        for level in levels:
            level_data = self.taxonomy_data[self.taxonomy_data['level'] == level]
            level_data = level_data.nlargest(15, 'relative_abundance')
            
            plt.figure(figsize=(12, 8))
            
            # Create horizontal bar plot
            y_pos = np.arange(len(level_data))
            bars = plt.barh(y_pos, level_data['relative_abundance'] * 100)
            
            # Color by abundance
            colors = plt.cm.Set3(np.linspace(0, 1, len(level_data)))
            for bar, color in zip(bars, colors):
                bar.set_color(color)
            
            plt.yticks(y_pos, level_data['taxon'])
            plt.xlabel('Relative Abundance (%)')
            plt.title(f'Taxonomic Composition at {level} Level')
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.savefig(output_path / f'taxonomy_{level.lower()}.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        self.logger.info(f"Taxonomy plots saved to {output_path}")
    
    def create_diversity_summary(self, output_dir: str = "plots"):
        """Create diversity metrics summary visualization"""
        
        if not self.diversity_metrics:
            self.calculate_diversity_metrics()
        
        if not self.diversity_metrics:
            return
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Create diversity metrics radar chart
        fig = go.Figure()
        
        metrics = ['richness', 'shannon', 'simpson', 'evenness']
        values = [self.diversity_metrics[m] for m in metrics]
        
        # Normalize values for radar chart (0-1 scale)
        normalized_values = []
        for i, (metric, value) in enumerate(zip(metrics, values)):
            if metric == 'richness':
                # Normalize richness to 0-1 scale (assuming max 100 species)
                normalized_values.append(min(value / 100, 1.0))
            else:
                # Other metrics are already 0-1 scale or close to it
                normalized_values.append(min(value, 1.0))
        
        fig.add_trace(go.Scatterpolar(
            r=normalized_values,
            theta=[m.title() for m in metrics],
            fill='toself',
            name='Diversity Metrics'
        ))
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1]
                )),
            title="Alpha Diversity Metrics",
            showlegend=True
        )
        
        fig.write_html(str(output_path / 'diversity_radar.html'))
        
        # Create diversity bar chart
        plt.figure(figsize=(10, 6))
        
        metrics_to_plot = ['shannon', 'simpson', 'evenness', 'dominance']
        values_to_plot = [self.diversity_metrics[m] for m in metrics_to_plot]
        
        bars = plt.bar(metrics_to_plot, values_to_plot, color=['skyblue', 'lightgreen', 'orange', 'pink'])
        plt.ylabel('Diversity Index Value')
        plt.title('Alpha Diversity Metrics')
        plt.xticks(rotation=45)
        
        # Add value labels on bars
        for bar, value in zip(bars, values_to_plot):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{value:.3f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(output_path / 'diversity_metrics.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Diversity plots saved to {output_path}")
    
    def create_interactive_dashboard(self, output_file: str = "metagenome_dashboard.html"):
        """Create interactive HTML dashboard with all visualizations"""
        
        if self.abundance_data is None:
            self.logger.error("No data available for dashboard")
            return
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Top Organisms', 'Taxonomic Composition', 
                           'Coverage vs Abundance', 'Diversity Metrics'),
            specs=[[{"type": "bar"}, {"type": "pie"}],
                   [{"type": "scatter"}, {"type": "bar"}]]
        )
        
        # Top organisms bar plot
        top_organisms = self.abundance_data.nlargest(15, 'relative_abundance')
        fig.add_trace(
            go.Bar(x=top_organisms['relative_abundance'] * 100,
                   y=top_organisms['organism_name'],
                   orientation='h',
                   name='Abundance'),
            row=1, col=1
        )
        
        # Taxonomic pie chart
        if self.taxonomy_data is not None:
            genus_data = self.taxonomy_data[self.taxonomy_data['level'] == 'Genus']
            if not genus_data.empty:
                top_genera = genus_data.nlargest(10, 'relative_abundance')
                fig.add_trace(
                    go.Pie(labels=top_genera['taxon'],
                           values=top_genera['relative_abundance'] * 100,
                           name="Genera"),
                    row=1, col=2
                )
        
        # Coverage vs abundance scatter
        if self.coverage_data is not None:
            merged_data = pd.merge(self.abundance_data, self.coverage_data, on='organism_id')
            fig.add_trace(
                go.Scatter(x=merged_data['coverage_breadth'] * 100,
                          y=merged_data['relative_abundance'] * 100,
                          mode='markers',
                          marker=dict(
                              size=merged_data['unique_kmers'] / 500,
                              color=merged_data['coverage_depth'],
                              colorscale='Viridis',
                              showscale=True
                          ),
                          text=merged_data['organism_name'],
                          name='Organisms'),
                row=2, col=1
            )
        
        # Diversity metrics bar chart
        if self.diversity_metrics:
            metrics = ['shannon', 'simpson', 'evenness']
            values = [self.diversity_metrics[m] for m in metrics]
            fig.add_trace(
                go.Bar(x=metrics, y=values, name='Diversity'),
                row=2, col=2
            )
        
        # Update layout
        fig.update_layout(
            height=800,
            title_text="Comprehensive Metagenomic Analysis Dashboard",
            showlegend=False
        )
        
        # Update axis labels
        fig.update_xaxes(title_text="Relative Abundance (%)", row=1, col=1)
        fig.update_xaxes(title_text="Coverage Breadth (%)", row=2, col=1)
        fig.update_xaxes(title_text="Diversity Metric", row=2, col=2)
        fig.update_yaxes(title_text="Organism", row=1, col=1)
        fig.update_yaxes(title_text="Relative Abundance (%)", row=2, col=1)
        fig.update_yaxes(title_text="Index Value", row=2, col=2)
        
        fig.write_html(output_file)
        self.logger.info(f"Interactive dashboard saved to {output_file}")
    
    def export_results_for_downstream(self, output_prefix: str):
        """Export results in formats suitable for downstream analysis"""
        
        if self.abundance_data is None:
            self.logger.error("No data available for export")
            return
        
        # Export for phyloseq (R package)
        phyloseq_dir = Path(f"{output_prefix}_phyloseq")
        phyloseq_dir.mkdir(exist_ok=True)
        
        # OTU table (abundance matrix)
        otu_table = self.abundance_data[['organism_id', 'relative_abundance']].copy()
        otu_table.columns = ['OTU_ID', 'Sample1']
        otu_table.to_csv(phyloseq_dir / 'otu_table.csv', index=False)
        
        # Taxonomy table
        if self.abundance_data is not None:
            tax_table = self.abundance_data[['organism_id', 'organism_name', 'taxonomy_path']].copy()
            # Parse taxonomy path into separate columns
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            for i, level in enumerate(tax_levels):
                tax_table[level] = tax_table['taxonomy_path'].apply(
                    lambda x: x.split(';')[i] if len(x.split(';')) > i else 'Unknown'
                )
            tax_table.drop('taxonomy_path', axis=1, inplace=True)
            tax_table.to_csv(phyloseq_dir / 'taxonomy_table.csv', index=False)
        
        # Export for QIIME2
        qiime_dir = Path(f"{output_prefix}_qiime2")
        qiime_dir.mkdir(exist_ok=True)
        
        # Feature table
        feature_table = self.abundance_data[['organism_id', 'relative_abundance']].copy()
        feature_table.columns = ['#OTU ID', 'Sample1']
        feature_table.to_csv(qiime_dir / 'feature_table.tsv', sep='\t', index=False)
        
        # Export summary statistics
        summary_stats = {
            'total_organisms': len(self.abundance_data),
            'diversity_metrics': self.diversity_metrics,
            'top_organisms': self.abundance_data.nlargest(10, 'relative_abundance')[
                ['organism_name', 'relative_abundance']
            ].to_dict('records')
        }
        
        with open(f"{output_prefix}_summary_stats.json", 'w') as f:
            json.dump(summary_stats, f, indent=2, default=str)
        
        self.logger.info(f"Results exported for downstream analysis:")
        self.logger.info(f"  - phyloseq format: {phyloseq_dir}")
        self.logger.info(f"  - QIIME2 format: {qiime_dir}")
        self.logger.info(f"  - Summary stats: {output_prefix}_summary_stats.json")
    
    def run_complete_analysis(self, 
                            fastq_file: str, 
                            output_prefix: str,
                            create_plots: bool = True,
                            create_dashboard: bool = True) -> Dict:
        """Run complete metagenomic analysis pipeline"""
        
        self.logger.info("Starting complete metagenomic analysis pipeline")
        
        # Validate inputs
        if not self.validate_inputs(fastq_file):
            raise ValueError("Input validation failed")
        
        # Run profiling
        results = self.run_comprehensive_profiling(fastq_file, output_prefix)
        
        # Calculate diversity metrics
        self.calculate_diversity_metrics()
        
        # Create visualizations
        if create_plots:
            self.logger.info("Creating visualization plots...")
            plot_dir = f"{output_prefix}_plots"
            self.create_abundance_plots(plot_dir)
            if self.taxonomy_data is not None:
                self.create_taxonomy_plots(plot_dir)
            self.create_diversity_summary(plot_dir)
        
        # Create interactive dashboard
        if create_dashboard:
            self.logger.info("Creating interactive dashboard...")
            self.create_interactive_dashboard(f"{output_prefix}_dashboard.html")
        
        # Export for downstream analysis
        self.logger.info("Exporting results for downstream analysis...")
        self.export_results_for_downstream(output_prefix)
        
        # Generate final summary
        summary = {
            'sample': fastq_file,
            'total_organisms_detected': len(self.abundance_data) if self.abundance_data is not None else 0,
            'diversity_metrics': self.diversity_metrics,
            'output_files': {
                'abundance_table': f"{output_prefix}_abundance_table.tsv",
                'taxonomy_summary': f"{output_prefix}_taxonomy_summary.tsv",
                'coverage_stats': f"{output_prefix}_coverage_stats.tsv",
                'plots_directory': f"{output_prefix}_plots",
                'dashboard': f"{output_prefix}_dashboard.html",
                'phyloseq_export': f"{output_prefix}_phyloseq",
                'qiime2_export': f"{output_prefix}_qiime2"
            }
        }
        
        self.logger.info("Complete analysis finished successfully!")
        self.logger.info(f"Detected {summary['total_organisms_detected']} organisms")
        self.logger.info(f"Shannon diversity: {self.diversity_metrics.get('shannon', 0):.3f}")
        
        return summary

def main():
    """Command line interface for comprehensive metagenomic analysis"""
    
    parser = argparse.ArgumentParser(
        description="Comprehensive Metagenomic Analysis Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete analysis
  python metagenomics_analyzer.py --database /data/microbes.db --fastq sample.fastq --output analysis

  # Run profiling only (no plots)
  python metagenomics_analyzer.py --database /data/microbes.db --fastq sample.fastq --output analysis --no-plots
        """)
    
    parser.add_argument('--database', required=True,
                       help='Path to comprehensive microbial database')
    parser.add_argument('--fastq', required=True,
                       help='Input FASTQ file')
    parser.add_argument('--output', required=True,
                       help='Output prefix for analysis results')
    parser.add_argument('--profiler', default='./build/hybrid_profiler',
                       help='Path to profiler binary')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip plot generation')
    parser.add_argument('--no-dashboard', action='store_true',
                       help='Skip interactive dashboard creation')
    parser.add_argument('--memory-limit', type=int, default=10,
                       help='Memory limit in GB for profiler')
    
    args = parser.parse_args()
    
    try:
        # Initialize analyzer
        analyzer = ComprehensiveMetagenomicsAnalyzer(
            profiler_binary=args.profiler,
            database_path=args.database
        )
        
        # Run analysis
        results = analyzer.run_complete_analysis(
            fastq_file=args.fastq,
            output_prefix=args.output,
            create_plots=not args.no_plots,
            create_dashboard=not args.no_dashboard
        )
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print("="*60)
        print(f"Sample: {results['sample']}")
        print(f"Organisms detected: {results['total_organisms_detected']}")
        print(f"Shannon diversity: {results['diversity_metrics'].get('shannon', 0):.3f}")
        print(f"Simpson diversity: {results['diversity_metrics'].get('simpson', 0):.3f}")
        print("\nOutput files:")
        for desc, path in results['output_files'].items():
            print(f"  {desc}: {path}")
        print("="*60)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()