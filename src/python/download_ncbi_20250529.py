#!/usr/bin/env python3
"""
NCBI Sequence Downloader for Fluoroquinolone Resistance Genes

Downloads GenBank sequences for gyrA, parC, parE, gyrB genes from organisms
listed in the resistance mutation CSV files. Downloads 100-500 isolates per
organism per gene to build comprehensive reference database.

Requires: biopython, pandas, requests
Install: pip install biopython pandas requests
"""

import pandas as pd
import requests
import time
import os
import sys
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import json
import argparse
from datetime import datetime

class NCBISequenceDownloader:
    def __init__(self, email, api_key=None):
        """
        Initialize NCBI downloader
        
        Args:
            email: Your email address (required by NCBI)
            api_key: NCBI API key (optional, increases rate limit)
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        # Rate limiting
        self.request_delay = 0.34 if api_key else 1.0  # seconds between requests
        
        # Gene mappings for different organisms
        self.gene_aliases = {
            'gyrA': ['gyrA', 'gyra', 'GyrA', 'DNA gyrase subunit A'],
            'gyrB': ['gyrB', 'gyrb', 'GyrB', 'DNA gyrase subunit B'],  
            'parC': ['parC', 'parc', 'ParC', 'topoisomerase IV subunit A'],
            'parE': ['parE', 'pare', 'ParE', 'topoisomerase IV subunit B'],
            'grlA': ['grlA', 'grla', 'GrlA'],  # S. aureus equivalent of parC
            'grlB': ['grlB', 'grlb', 'GrlB']   # S. aureus equivalent of parE
        }
        
        # Species name standardization
        self.species_standardization = {
            'Escherichia coli': 'Escherichia coli',
            'Salmonella enterica': 'Salmonella enterica',
            'Campylobacter jejuni': 'Campylobacter jejuni',
            'Clostridioides difficile': 'Clostridioides difficile',
            'Klebsiella pneumoniae': 'Klebsiella pneumoniae',
            'Pseudomonas aeruginosa': 'Pseudomonas aeruginosa',
            'Bacteroides fragilis': 'Bacteroides fragilis',
            'Enterococcus faecium': 'Enterococcus faecium',
            'Lactobacillus delbrueckii': 'Lactobacillus delbrueckii',
            'Lactobacillus acidophilus': 'Lactobacillus acidophilus',
            'Lactobacillus plantarum': 'Lactobacillus plantarum',
            'Acinetobacter baumannii': 'Acinetobacter baumannii',
            'Serratia marcescens': 'Serratia marcescens',
            'Enterobacter cloacae': 'Enterobacter cloacae',
            'Staphylococcus aureus': 'Staphylococcus aureus',
            'Streptococcus pneumoniae': 'Streptococcus pneumoniae'
        }
    
    def load_csv_files(self, mutations_csv, efflux_csv):
        """Load and parse the resistance CSV files"""
        print(f"Loading mutation data from {mutations_csv}")
        mutations_df = pd.read_csv(mutations_csv)
        
        print(f"Loading efflux data from {efflux_csv}")  
        efflux_df = pd.read_csv(efflux_csv)
        
        # Get unique organisms and genes
        organisms_from_mutations = set(mutations_df['Species'].unique())
        organisms_from_efflux = set(efflux_df['Species'].unique())
        
        all_organisms = organisms_from_mutations.union(organisms_from_efflux)
        genes_from_mutations = set(mutations_df['Gene'].unique())
        efflux_genes = set(efflux_df['Efflux Gene'].unique())
        
        print(f"\nFound {len(all_organisms)} unique organisms")
        print(f"Resistance genes: {genes_from_mutations}")
        print(f"Efflux pump genes: {len(efflux_genes)} total")
        
        return mutations_df, efflux_df, all_organisms, genes_from_mutations
    
    def search_sequences(self, organism, gene, max_sequences=500):
        """
        Search NCBI for sequences of a specific gene from an organism
        
        Args:
            organism: Species name (e.g., "Escherichia coli")
            gene: Gene name (e.g., "gyrA")
            max_sequences: Maximum number of sequences to retrieve
            
        Returns:
            List of GenBank accession numbers
        """
        print(f"  Searching for {gene} sequences from {organism}...")
        
        # Standardize organism name
        organism = self.species_standardization.get(organism, organism)
        
        # Build search query
        gene_terms = self.gene_aliases.get(gene, [gene])
        gene_query = " OR ".join([f'"{term}"[Gene Name]' for term in gene_terms])
        
        # More comprehensive search query
        search_query = f'({gene_query}) AND "{organism}"[Organism] AND ("complete cds"[All Fields] OR "partial cds"[All Fields])'
        
        try:
            # Search for sequences
            search_handle = Entrez.esearch(
                db="nuccore",
                term=search_query,
                retmax=max_sequences,
                retmode="xml"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            accession_list = search_results["IdList"]
            
            if len(accession_list) == 0:
                # Try alternative search with just gene name
                alt_query = f'{gene}[Gene Name] AND "{organism}"[Organism]'
                print(f"    No results with primary query, trying: {alt_query}")
                
                search_handle = Entrez.esearch(
                    db="nuccore", 
                    term=alt_query,
                    retmax=max_sequences,
                    retmode="xml"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                accession_list = search_results["IdList"]
            
            print(f"    Found {len(accession_list)} sequences")
            time.sleep(self.request_delay)
            
            return accession_list
            
        except Exception as e:
            print(f"    ERROR searching for {gene} in {organism}: {e}")
            return []
    
    def download_genbank_records(self, accession_list, batch_size=50):
        """
        Download GenBank records for a list of accessions
        
        Args:
            accession_list: List of GenBank accession numbers
            batch_size: Number of records to download per request
            
        Returns:
            List of SeqRecord objects
        """
        if not accession_list:
            return []
        
        records = []
        
        # Download in batches to avoid timeout
        for i in range(0, len(accession_list), batch_size):
            batch = accession_list[i:i + batch_size]
            print(f"    Downloading batch {i//batch_size + 1} ({len(batch)} records)...")
            
            try:
                fetch_handle = Entrez.efetch(
                    db="nuccore",
                    id=batch,
                    rettype="gbwithparts",
                    retmode="text"
                )
                
                batch_records = list(SeqIO.parse(fetch_handle, "genbank"))
                fetch_handle.close()
                
                records.extend(batch_records)
                print(f"      Successfully downloaded {len(batch_records)} records")
                
                time.sleep(self.request_delay)
                
            except Exception as e:
                print(f"      ERROR downloading batch: {e}")
                continue
        
        return records
    
    def extract_gene_features(self, record, gene_name):
        """
        Extract gene features from a GenBank record
        
        Args:
            record: SeqRecord object from GenBank
            gene_name: Name of gene to extract
            
        Returns:
            Dictionary with gene information
        """
        gene_info = {
            'accession': record.id,
            'organism': record.annotations.get('organism', 'Unknown'),
            'description': record.description,
            'sequence_length': len(record.seq),
            'gene_features': []
        }
        
        gene_aliases = self.gene_aliases.get(gene_name, [gene_name])
        
        for feature in record.features:
            if feature.type == "CDS":
                # Check if this CDS corresponds to our target gene
                gene_qual = feature.qualifiers.get('gene', [])
                product_qual = feature.qualifiers.get('product', [])
                
                # Check if gene name matches
                gene_match = any(alias.lower() in [g.lower() for g in gene_qual] 
                               for alias in gene_aliases)
                
                # Check product description
                product_match = any(alias.lower() in product.lower() 
                                  for product in product_qual 
                                  for alias in gene_aliases)
                
                if gene_match or product_match:
                    # Extract sequence coordinates
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    
                    # Extract nucleotide and protein sequences
                    gene_seq = feature.extract(record.seq)
                    
                    protein_seq = ""
                    if 'translation' in feature.qualifiers:
                        protein_seq = feature.qualifiers['translation'][0]
                    
                    gene_feature = {
                        'gene_name': gene_qual[0] if gene_qual else gene_name,
                        'product': product_qual[0] if product_qual else 'Unknown',
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'nucleotide_sequence': str(gene_seq),
                        'protein_sequence': protein_seq,
                        'length_nt': len(gene_seq),
                        'length_aa': len(protein_seq)
                    }
                    
                    gene_info['gene_features'].append(gene_feature)
        
        return gene_info
    
    def save_sequences(self, gene_data, output_dir):
        """Save downloaded sequences to files"""
        os.makedirs(output_dir, exist_ok=True)
        
        for organism in gene_data:
            org_dir = os.path.join(output_dir, organism.replace(' ', '_'))
            os.makedirs(org_dir, exist_ok=True)
            
            for gene in gene_data[organism]:
                gene_file = os.path.join(org_dir, f"{gene}.json")
                
                with open(gene_file, 'w') as f:
                    json.dump(gene_data[organism][gene], f, indent=2)
                
                print(f"  Saved {len(gene_data[organism][gene])} {gene} sequences for {organism}")
    
    def download_resistance_genes(self, mutations_csv, efflux_csv, output_dir, 
                                max_per_gene=300, target_organisms=None):
        """
        Main function to download all resistance gene sequences
        
        Args:
            mutations_csv: Path to mutations CSV file
            efflux_csv: Path to efflux genes CSV file  
            output_dir: Directory to save sequences
            max_per_gene: Maximum sequences per organism per gene
            target_organisms: List of specific organisms to download (optional)
        """
        print("=== NCBI Resistance Gene Sequence Downloader ===")
        print(f"Started at: {datetime.now()}")
        
        # Load CSV files
        mutations_df, efflux_df, all_organisms, resistance_genes = self.load_csv_files(
            mutations_csv, efflux_csv
        )
        
        # Filter organisms if specified
        if target_organisms:
            all_organisms = [org for org in all_organisms if org in target_organisms]
            print(f"Limiting to {len(all_organisms)} specified organisms")
        
        # Initialize data structure
        gene_data = defaultdict(lambda: defaultdict(list))
        download_stats = defaultdict(lambda: defaultdict(int))
        
        # Download resistance genes (gyrA, parC, etc.)
        print(f"\n=== Downloading Resistance Genes ===")
        for organism in sorted(all_organisms):
            if organism == "Various Gram-negatives":
                continue  # Skip generic entries
                
            print(f"\nProcessing {organism}...")
            
            # Get relevant genes for this organism from mutations CSV
            org_mutations = mutations_df[mutations_df['Species'] == organism]
            genes_for_org = set(org_mutations['Gene'].unique())
            
            for gene in sorted(genes_for_org):
                accessions = self.search_sequences(organism, gene, max_per_gene)
                
                if accessions:
                    records = self.download_genbank_records(accessions)
                    
                    for record in records:
                        gene_info = self.extract_gene_features(record, gene)
                        if gene_info['gene_features']:  # Only save if gene was found
                            gene_data[organism][gene].append(gene_info)
                            download_stats[organism][gene] += 1
        
        # Download efflux pump genes
        print(f"\n=== Downloading Efflux Pump Genes ===")
        for _, row in efflux_df.iterrows():
            organism = row['Species']
            efflux_gene = row['Efflux Gene']
            
            if organism == "Various Gram-negatives":
                continue
                
            if target_organisms and organism not in target_organisms:
                continue
            
            print(f"\nProcessing {efflux_gene} from {organism}...")
            
            accessions = self.search_sequences(organism, efflux_gene, max_per_gene)
            
            if accessions:
                records = self.download_genbank_records(accessions)
                
                for record in records:
                    gene_info = self.extract_gene_features(record, efflux_gene)
                    if gene_info['gene_features']:
                        gene_data[organism][efflux_gene].append(gene_info)
                        download_stats[organism][efflux_gene] += 1
        
        # Save all sequences
        print(f"\n=== Saving Sequences ===")
        self.save_sequences(gene_data, output_dir)
        
        # Print summary
        print(f"\n=== Download Summary ===")
        total_sequences = 0
        for organism in sorted(download_stats.keys()):
            print(f"\n{organism}:")
            for gene in sorted(download_stats[organism].keys()):
                count = download_stats[organism][gene]
                print(f"  {gene}: {count} sequences")
                total_sequences += count
        
        print(f"\nTotal sequences downloaded: {total_sequences}")
        print(f"Completed at: {datetime.now()}")
        
        return gene_data, download_stats

def main():
    parser = argparse.ArgumentParser(description='Download resistance gene sequences from NCBI')
    parser.add_argument('mutations_csv', help='Path to mutations CSV file')
    parser.add_argument('efflux_csv', help='Path to efflux genes CSV file')
    parser.add_argument('output_dir', help='Output directory for sequences')
    parser.add_argument('--email', required=True, help='Your email address (required by NCBI)')
    parser.add_argument('--api-key', help='NCBI API key (optional, increases rate limit)')
    parser.add_argument('--max-per-gene', type=int, default=300, 
                       help='Maximum sequences per organism per gene (default: 300)')
    parser.add_argument('--organisms', nargs='+', 
                       help='Specific organisms to download (optional)')
    
    args = parser.parse_args()
    
    # Create downloader
    downloader = NCBISequenceDownloader(args.email, args.api_key)
    
    # Download sequences
    gene_data, stats = downloader.download_resistance_genes(
        args.mutations_csv,
        args.efflux_csv, 
        args.output_dir,
        args.max_per_gene,
        args.organisms
    )
    
    print(f"\nSequences saved to: {args.output_dir}")
    print("Ready for GPU pipeline integration!")

if __name__ == "__main__":
    main()