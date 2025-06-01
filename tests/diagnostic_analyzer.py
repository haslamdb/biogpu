#!/usr/bin/env python3
"""
diagnostic_analyzer.py
Analyze BioGPU diagnostic reports and provide troubleshooting insights
"""

import sys
import json
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class BioGPUDiagnosticAnalyzer:
    def __init__(self, diagnostic_file, json_file=None, hdf5_file=None):
        self.diagnostic_file = Path(diagnostic_file)
        self.json_file = Path(json_file) if json_file else None
        self.hdf5_file = Path(hdf5_file) if hdf5_file else None
        
        self.stats = {}
        self.alignments = []
        self.issues = []
        
    def parse_diagnostic_report(self):
        """Parse the text diagnostic report"""
        print(f"📋 Parsing diagnostic report: {self.diagnostic_file}")
        
        if not self.diagnostic_file.exists():
            print(f"❌ Diagnostic file not found: {self.diagnostic_file}")
            return False
            
        with open(self.diagnostic_file, 'r') as f:
            content = f.read()
            
        # Extract key statistics
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if "Total input reads:" in line:
                self.stats['total_reads'] = int(line.split()[-1])
            elif "After Bloom filter:" in line:
                parts = line.split()
                self.stats['reads_after_bloom'] = int(parts[-2])
                self.stats['bloom_retention'] = float(parts[-1].strip('()%'))
            elif "After k-mer enrichment:" in line:
                parts = line.split()
                self.stats['reads_after_kmer'] = int(parts[-2])
                self.stats['kmer_retention'] = float(parts[-1].strip('()%'))
            elif "Total protein alignments:" in line:
                self.stats['protein_alignments'] = int(line.split()[-1])
            elif "Total mutations found:" in line:
                self.stats['mutations_found'] = int(line.split()[-1])
            elif "QRDR mutations found:" in line:
                self.stats['qrdr_mutations'] = int(line.split()[-1])
                
        return True
    
    def parse_json_output(self):
        """Parse the JSON output file"""
        if not self.json_file or not self.json_file.exists():
            print("📄 No JSON file provided or file not found")
            return
            
        print(f"📄 Parsing JSON output: {self.json_file}")
        
        with open(self.json_file, 'r') as f:
            data = json.load(f)
            
        # Extract mutations
        if 'mutations' in data:
            self.alignments = data['mutations']
            
        # Extract enhanced summary if available
        if 'enhanced_summary' in data:
            summary = data['enhanced_summary']
            self.stats.update({
                'enhanced_total_reads': summary.get('total_reads', 0),
                'enhanced_protein_matches': summary.get('protein_matches', 0),
                'enhanced_resistance_mutations': summary.get('resistance_mutations_found', 0),
                'enhanced_qrdr_alignments': summary.get('qrdr_alignments', 0),
                'enhanced_high_confidence': summary.get('high_confidence_matches', 0)
            })
    
    def analyze_pipeline_performance(self):
        """Analyze pipeline performance and identify issues"""
        print("\n🔍 PIPELINE PERFORMANCE ANALYSIS")
        print("=" * 50)
        
        # Calculate filter efficiency
        total_reads = self.stats.get('total_reads', 0)
        bloom_reads = self.stats.get('reads_after_bloom', 0)
        kmer_reads = self.stats.get('reads_after_kmer', 0)
        protein_alignments = self.stats.get('protein_alignments', 0)
        mutations_found = self.stats.get('mutations_found', 0)
        
        if total_reads > 0:
            bloom_efficiency = (total_reads - bloom_reads) / total_reads * 100
            kmer_efficiency = (bloom_reads - kmer_reads) / bloom_reads * 100 if bloom_reads > 0 else 0
            protein_rate = protein_alignments / total_reads * 100 if total_reads > 0 else 0
            
            print(f"📊 Filter Performance:")
            print(f"   Bloom filter removed: {bloom_efficiency:.1f}% of reads")
            print(f"   K-mer filter removed: {kmer_efficiency:.1f}% of remaining reads")
            print(f"   Protein hit rate: {protein_rate:.3f}%")
            
            # Identify potential issues
            if bloom_efficiency > 98:
                self.issues.append("🚨 Bloom filter removed >98% of reads - input may not contain target sequences")
            elif bloom_efficiency < 50:
                self.issues.append("⚠️ Bloom filter removed <50% of reads - may need optimization")
                
            if protein_rate < 0.1:
                self.issues.append("🚨 Very low protein hit rate (<0.1%) - check database or thresholds")
            elif protein_rate > 20:
                self.issues.append("📈 High protein hit rate (>20%) - good target enrichment")
                
            if mutations_found == 0 and protein_alignments > 0:
                self.issues.append("❓ Protein alignments found but no mutations - likely using mutant reference database")
    
    def analyze_mutation_detection(self):
        """Analyze mutation detection results"""
        print(f"\n🧬 MUTATION DETECTION ANALYSIS")
        print("=" * 50)
        
        mutations_found = self.stats.get('mutations_found', 0)
        qrdr_mutations = self.stats.get('qrdr_mutations', 0)
        protein_alignments = self.stats.get('protein_alignments', 0)
        
        print(f"📈 Mutation Statistics:")
        print(f"   Total mutations: {mutations_found}")
        print(f"   QRDR mutations: {qrdr_mutations}")
        print(f"   Protein alignments: {protein_alignments}")
        
        if mutations_found == 0:
            if protein_alignments > 0:
                print("\n❌ MUTATION DETECTION ISSUE IDENTIFIED:")
                print("   • Protein alignments found but no mutations detected")
                print("   • This typically indicates a MUTANT REFERENCE DATABASE")
                print("   • Mutant reads align perfectly to mutant references")
                print("\n💡 SOLUTIONS:")
                print("   1. Use wild-type reference sequences for mutation calling")
                print("   2. Implement resistance pattern recognition instead of mutation calling")
                print("   3. Check if alignments map to known resistance positions")
            else:
                print("\n❌ NO MUTATIONS OR ALIGNMENTS FOUND:")
                print("   • May indicate reads don't contain resistance genes")
                print("   • Or pipeline parameters need adjustment")
        else:
            resistance_rate = qrdr_mutations / mutations_found * 100 if mutations_found > 0 else 0
            print(f"   QRDR mutation rate: {resistance_rate:.1f}%")
            
            if qrdr_mutations > 0:
                print("\n✅ RESISTANCE MUTATIONS DETECTED")
                print("   • Fluoroquinolone resistance likely present")
                self.issues.append("🦠 Fluoroquinolone resistance mutations detected")
    
    def analyze_alignments(self):
        """Analyze individual alignments from JSON"""
        if not self.alignments:
            print("\n📄 No alignment details available in JSON")
            return
            
        print(f"\n🎯 ALIGNMENT ANALYSIS")
        print("=" * 50)
        print(f"Found {len(self.alignments)} mutations in JSON output")
        
        # Analyze by gene
        gene_counts = {}
        for alignment in self.alignments:
            gene_id = alignment.get('gene_id', 'unknown')
            gene_counts[gene_id] = gene_counts.get(gene_id, 0) + 1
            
        print("\n📊 Mutations by gene:")
        gene_names = {0: 'gyrA', 1: 'parC', 2: 'gyrB', 3: 'parE'}
        for gene_id, count in gene_counts.items():
            gene_name = gene_names.get(gene_id, f'Gene_{gene_id}')
            print(f"   {gene_name}: {count} mutations")
            
        # Show top alignments
        print(f"\n🔝 Top 5 alignments by score:")
        sorted_alignments = sorted(self.alignments, 
                                 key=lambda x: x.get('alignment_score', 0), 
                                 reverse=True)
        
        for i, alignment in enumerate(sorted_alignments[:5]):
            gene_name = gene_names.get(alignment.get('gene_id', 'unknown'), 'Unknown')
            score = alignment.get('alignment_score', 0)
            identity = alignment.get('identity', 0)
            read_id = alignment.get('read_pair', alignment.get('read_id', 'unknown'))
            
            print(f"   {i+1}. Read {read_id}, {gene_name}, Score: {score:.1f}, Identity: {identity:.1%}")
    
    def generate_recommendations(self):
        """Generate specific recommendations based on analysis"""
        print(f"\n💡 RECOMMENDATIONS")
        print("=" * 50)
        
        if not self.issues:
            print("✅ No major issues detected - pipeline appears to be working well")
            return
            
        print("Issues identified:")
        for issue in self.issues:
            print(f"   {issue}")
            
        print(f"\n🔧 Specific recommendations:")
        
        mutations_found = self.stats.get('mutations_found', 0)
        protein_alignments = self.stats.get('protein_alignments', 0)
        
        if mutations_found == 0 and protein_alignments > 0:
            print("""
   DATABASE ISSUE - MUTANT REFERENCES DETECTED:
   
   1. IMMEDIATE FIX:
      • Check if your protein database contains mutant sequences
      • If yes, rebuild with wild-type sequences for mutation calling
      
   2. ALTERNATIVE APPROACH:
      • Use position-based resistance detection
      • Check if alignments map to known resistance positions (83, 87 for gyrA)
      • Infer resistance from alignment location rather than mutation calling
      
   3. VERIFICATION:
      • Examine the "Top Alignments" section in diagnostic report
      • Check if alignments cover QRDR regions (positions 80-90 for gyrA/parC)
      • Look for high-identity matches to resistance genes
   """)
            
        bloom_retention = self.stats.get('bloom_retention', 100)
        if bloom_retention < 5:
            print("""
   LOW BLOOM FILTER RETENTION:
   
   • Input reads may not contain target resistance genes
   • Verify input files are from resistance gene-containing samples
   • Consider using positive control samples
   """)
            
        protein_rate = protein_alignments / self.stats.get('total_reads', 1) * 100
        if protein_rate < 0.1:
            print("""
   LOW PROTEIN HIT RATE:
   
   • Consider lowering identity threshold (current: 90%)
   • Check protein database completeness
   • Verify 6-frame translation is working correctly
   """)
    
    def create_summary_plot(self, output_path=None):
        """Create summary visualization"""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
            fig.suptitle('BioGPU Pipeline Performance Summary', fontsize=16)
            
            # Filter efficiency
            total = self.stats.get('total_reads', 1)
            bloom = self.stats.get('reads_after_bloom', 0)
            kmer = self.stats.get('reads_after_kmer', 0)
            
            categories = ['Input', 'After Bloom', 'After K-mer']
            values = [total, bloom, kmer]
            colors = ['lightblue', 'orange', 'green']
            
            ax1.bar(categories, values, color=colors)
            ax1.set_ylabel('Number of Reads')
            ax1.set_title('Read Filtering Performance')
            ax1.tick_params(axis='x', rotation=45)
            
            # Mutation detection
            protein_aligns = self.stats.get('protein_alignments', 0)
            mutations = self.stats.get('mutations_found', 0)
            qrdr_mutations = self.stats.get('qrdr_mutations', 0)
            
            ax2.bar(['Protein\nAlignments', 'Total\nMutations', 'QRDR\nMutations'], 
                   [protein_aligns, mutations, qrdr_mutations],
                   color=['skyblue', 'lightcoral', 'red'])
            ax2.set_ylabel('Count')
            ax2.set_title('Mutation Detection Results')
            
            # Retention rates
            if total > 0:
                bloom_ret = bloom / total * 100
                kmer_ret = kmer / total * 100
                protein_ret = protein_aligns / total * 100
                
                ax3.bar(['Bloom\nRetention', 'K-mer\nRetention', 'Protein\nHit Rate'], 
                       [bloom_ret, kmer_ret, protein_ret],
                       color=['lightgreen', 'yellow', 'purple'])
                ax3.set_ylabel('Percentage (%)')
                ax3.set_title('Pipeline Efficiency')
            
            # Issue summary
            issue_types = []
            issue_counts = []
            
            for issue in self.issues:
                if "mutation" in issue.lower():
                    issue_types.append("Mutation\nDetection")
                elif "bloom" in issue.lower():
                    issue_types.append("Bloom\nFilter")
                elif "protein" in issue.lower():
                    issue_types.append("Protein\nSearch")
                else:
                    issue_types.append("Other")
                    
            if issue_types:
                issue_df = pd.Series(issue_types).value_counts()
                ax4.pie(issue_df.values, labels=issue_df.index, autopct='%1.0f%%')
                ax4.set_title('Issues Detected')
            else:
                ax4.text(0.5, 0.5, 'No Issues\nDetected', 
                        ha='center', va='center', fontsize=14, color='green')
                ax4.set_xlim(0, 1)
                ax4.set_ylim(0, 1)
                ax4.set_title('Pipeline Status')
            
            plt.tight_layout()
            
            if output_path:
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                print(f"📊 Summary plot saved: {output_path}")
            else:
                plt.show()
                
        except ImportError:
            print("📊 Matplotlib not available - skipping plot generation")
    
    def run_analysis(self):
        """Run complete analysis"""
        print("🔬 BioGPU Diagnostic Analysis")
        print("=" * 50)
        
        # Parse files
        if not self.parse_diagnostic_report():
            return False
            
        self.parse_json_output()
        
        # Run analyses
        self.analyze_pipeline_performance()
        self.analyze_mutation_detection()
        self.analyze_alignments()
        self.generate_recommendations()
        
        # Create summary
        print(f"\n📋 ANALYSIS SUMMARY")
        print("=" * 50)
        print(f"Total reads processed: {self.stats.get('total_reads', 'N/A')}")
        print(f"Protein alignments: {self.stats.get('protein_alignments', 'N/A')}")
        print(f"Mutations detected: {self.stats.get('mutations_found', 'N/A')}")
        print(f"Issues identified: {len(self.issues)}")
        
        return True

def main():
    parser = argparse.ArgumentParser(description='Analyze BioGPU diagnostic reports')
    parser.add_argument('diagnostic_file', help='Path to diagnostic report file')
    parser.add_argument('--json', help='Path to JSON output file')
    parser.add_argument('--hdf5', help='Path to HDF5 output file')
    parser.add_argument('--plot', help='Save summary plot to file')
    
    args = parser.parse_args()
    
    analyzer = BioGPUDiagnosticAnalyzer(args.diagnostic_file, args.json, args.hdf5)
    
    if analyzer.run_analysis():
        if args.plot:
            analyzer.create_summary_plot(args.plot)
        print("\n✅ Analysis complete!")
    else:
        print("\n❌ Analysis failed!")
        return 1
        
    return 0

if __name__ == '__main__':
    sys.exit(main())