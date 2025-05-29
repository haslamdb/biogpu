#!/usr/bin/env python3
"""
BioGPU Fluoroquinolone Resistance Gene Reference Downloader
Downloads reference sequences for all resistance genes from NCBI and builds a targeted database
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import time
import pickle
import json
import logging
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# IMPORTANT: Set your email for NCBI Entrez
Entrez.email = "your.email@example.com"  # Change this!

@dataclass
class ReferenceGene:
    """Reference sequence for a resistance gene"""
    gene_name: str
    species: str
    protein_accession: str
    nucleotide_accession: str
    nucleotide_sequence: str
    protein_sequence: str
    gene_start: int
    gene_end: int
    flanking_5prime: str = ""
    flanking_3prime: str = ""
    mutations: List[Dict] = field(default_factory=list)
    
    @property
    def full_sequence(self) -> str:
        """Return gene sequence with flanking regions"""
        return self.flanking_5prime + self.nucleotide_sequence + self.flanking_3prime
    
    @property
    def gene_length(self) -> int:
        return len(self.nucleotide_sequence)

class ResistanceGeneDownloader:
    def __init__(self, mutations_csv: str, genes_csv: str, flanking_size: int = 500):
        self.mutations_csv = mutations_csv
        self.genes_csv = genes_csv
        self.flanking_size = flanking_size
        self.reference_genes = {}
        self.failed_downloads = []
        
        # Target genes for fluoroquinolone resistance
        self.target_genes = {'gyrA', 'gyrB', 'parC', 'parE', 'grlA', 'grlB'}
        
        # Gene name variations (some species use different names)
        self.gene_aliases = {
            'grlA': ['parC'],  # In S. aureus, grlA is equivalent to parC
            'grlB': ['parE'],  # In S. aureus, grlB is equivalent to parE
        }
        
    def load_mutation_data(self) -> Dict[str, List[Dict]]:
        """Load mutation data from CSV"""
        logger.info("Loading mutation data...")
        df = pd.read_csv(self.mutations_csv)
        
        mutations_by_key = defaultdict(list)
        for _, row in df.iterrows():
            if pd.isna(row['wt']) or pd.isna(row['mut']):
                continue
                
            key = f"{row['species']}_{row['gene']}_{row['Protein']}"
            mutations_by_key[key].append({
                'position': int(row['location']),
                'wild_type': row['wt'],
                'mutant': row['mut'],
                'start': int(row['start']),
                'stop': int(row['stop'])
            })
            
        return mutations_by_key
    
    def load_resistance_genes(self) -> List[Dict]:
        """Load resistance gene data from CSV"""
        logger.info("Loading resistance gene data...")
        df = pd.read_csv(self.genes_csv)
        
        genes = []
        for _, row in df.iterrows():
            genes.append({
                'species': row['species'],
                'gene': row['gene'],
                'protein_accession': row['Protein'],
                'start': int(row['start']),
                'stop': int(row['stop'])
            })
            
        return genes
    
    def fetch_protein_record(self, accession: str, retries: int = 3) -> dict:
        """Fetch protein record from NCBI"""
        for attempt in range(retries):
            try:
                time.sleep(0.5)  # Be nice to NCBI
                handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                handle.close()
                return record
            except Exception as e:
                logger.warning(f"Attempt {attempt + 1} failed for {accession}: {e}")
                if attempt == retries - 1:
                    raise
                time.sleep(2 ** attempt)  # Exponential backoff
        return None
    
    def fetch_nucleotide_with_flanking(self, nucleotide_acc: str, start: int, end: int, 
                                     is_complement: bool = False) -> Tuple[str, str, str]:
        """Fetch nucleotide sequence with flanking regions"""
        try:
            # Calculate positions with flanking
            flank_start = max(1, start - self.flanking_size)
            flank_end = end + self.flanking_size
            
            time.sleep(0.5)
            
            # Try different approaches for fetching
            try:
                # First try: Direct fetch with seq_start and seq_stop
                handle = Entrez.efetch(
                    db="nucleotide", 
                    id=nucleotide_acc,
                    rettype="fasta",
                    seq_start=flank_start,
                    seq_stop=flank_end
                )
                record = SeqIO.read(handle, "fasta")
                handle.close()
            except:
                # Second try: Fetch full sequence and extract region
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=nucleotide_acc,
                    rettype="gb",
                    retmode="text"
                )
                record = SeqIO.read(handle, "genbank")
                handle.close()
                
                # Extract the region
                full_seq = str(record.seq)
                if flank_end > len(full_seq):
                    flank_end = len(full_seq)
                    
                extracted_seq = full_seq[flank_start-1:flank_end]
                record.seq = Seq(extracted_seq)
            
            full_seq = str(record.seq).upper()
            
            # If it's a complement, reverse complement the sequence
            if is_complement:
                full_seq = str(Seq(full_seq).reverse_complement())
            
            # Calculate the actual positions in the extracted sequence
            actual_start = start - flank_start + 1
            actual_end = end - flank_start + 1
            
            # Extract parts
            flanking_5prime = full_seq[:actual_start-1] if actual_start > 1 else ""
            gene_seq = full_seq[actual_start-1:actual_end]
            flanking_3prime = full_seq[actual_end:] if actual_end < len(full_seq) else ""
            
            return flanking_5prime, gene_seq, flanking_3prime
            
        except Exception as e:
            logger.error(f"Failed to fetch nucleotide {nucleotide_acc}: {e}")
            return "", "", ""
    
    def process_protein_accession(self, accession: str, species: str, gene_name: str, 
                                gene_info: Dict = None) -> ReferenceGene:
        """Process a single protein accession to get reference sequence"""
        try:
            # Fetch protein record
            protein_record = self.fetch_protein_record(accession)
            if not protein_record:
                raise ValueError(f"Could not fetch protein record {accession}")
            
            protein_seq = str(protein_record.seq)
            
            # Find nucleotide accession and location from protein record
            nucleotide_acc = None
            gene_start = None
            gene_end = None
            is_complement = False
            
            # Look for CDS feature with translation
            for feature in protein_record.features:
                if feature.type == "CDS":
                    # Get coded_by qualifier
                    if "coded_by" in feature.qualifiers:
                        coded_by = feature.qualifiers["coded_by"][0]
                        # Parse format like "NC_000913.3:2407763..2410355" or "complement(AATJPX010000002.1:65736..67628)"
                        
                        # Check if it's complement
                        if "complement(" in coded_by:
                            is_complement = True
                            coded_by = coded_by.replace("complement(", "").replace(")", "")
                        
                        if ":" in coded_by:
                            nuc_acc, positions = coded_by.split(":", 1)
                            nucleotide_acc = nuc_acc  # Keep the full accession with version
                            
                            # Parse positions
                            if ".." in positions:
                                start_str, end_str = positions.split("..")
                                # Clean up the position strings - remove all non-digit characters
                                start_str = ''.join(c for c in start_str if c.isdigit())
                                end_str = ''.join(c for c in end_str if c.isdigit())
                                
                                if start_str and end_str:
                                    gene_start = int(start_str)
                                    gene_end = int(end_str)
            
            # Use provided coordinates if available and we don't have them
            if gene_info and not gene_start:
                gene_start = gene_info.get('start')
                gene_end = gene_info.get('stop')
            
            if not nucleotide_acc:
                # Try alternative approach - look for db_xref
                for feature in protein_record.features:
                    if feature.type == "CDS" and "db_xref" in feature.qualifiers:
                        for xref in feature.qualifiers["db_xref"]:
                            if xref.startswith("GeneID:"):
                                # Could use this to look up nucleotide sequence
                                logger.info(f"Found GeneID: {xref} for {accession}")
            
            if not (nucleotide_acc and gene_start and gene_end):
                # Try using the protein accession itself to find linked nucleotide
                try:
                    # Search for nucleotide records linked to this protein
                    search_handle = Entrez.esearch(db="nucleotide", term=f"{accession}[Protein Accession]")
                    search_results = Entrez.read(search_handle)
                    search_handle.close()
                    
                    if search_results["IdList"]:
                        nucleotide_id = search_results["IdList"][0]
                        # Fetch the nucleotide record to get accession
                        fetch_handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="acc")
                        nucleotide_acc = fetch_handle.read().strip()
                        fetch_handle.close()
                        logger.info(f"Found linked nucleotide {nucleotide_acc} for protein {accession}")
                except:
                    pass
            
            if not (nucleotide_acc and gene_start and gene_end):
                raise ValueError(f"Could not determine nucleotide location for {accession}")
            
            # Fetch nucleotide sequence with flanking
            flank_5p, gene_seq, flank_3p = self.fetch_nucleotide_with_flanking(
                nucleotide_acc, gene_start, gene_end, is_complement
            )
            
            if not gene_seq:
                raise ValueError(f"Could not fetch nucleotide sequence for {accession}")
            
            # Create reference gene object
            ref_gene = ReferenceGene(
                gene_name=gene_name,
                species=species,
                protein_accession=accession,
                nucleotide_accession=nucleotide_acc,
                nucleotide_sequence=gene_seq,
                protein_sequence=protein_seq,
                gene_start=gene_start,
                gene_end=gene_end,
                flanking_5prime=flank_5p,
                flanking_3prime=flank_3p
            )
            
            return ref_gene
            
        except Exception as e:
            logger.error(f"Failed to process {accession} ({species} {gene_name}): {e}")
            self.failed_downloads.append({
                'accession': accession,
                'species': species,
                'gene': gene_name,
                'error': str(e)
            })
            return None
    
    def download_mutation_gene_references(self):
        """Download reference sequences for genes with known mutations"""
        logger.info("Downloading reference sequences for mutation genes...")
        
        mutations_by_key = self.load_mutation_data()
        processed_accessions = set()
        
        for key, mutations in mutations_by_key.items():
            species, gene, protein_acc = key.rsplit('_', 2)
            
            # Skip if not a target gene
            if gene not in self.target_genes:
                continue
                
            # Skip if already processed
            if protein_acc in processed_accessions:
                continue
                
            processed_accessions.add(protein_acc)
            
            logger.info(f"Processing {species} {gene} ({protein_acc})")
            
            # Get coordinate info from mutations
            gene_info = {
                'start': min(m['start'] for m in mutations),
                'stop': max(m['stop'] for m in mutations)
            }
            
            ref_gene = self.process_protein_accession(
                protein_acc, species, gene, gene_info
            )
            
            if ref_gene:
                ref_gene.mutations = mutations
                self.reference_genes[f"{species}_{gene}_{protein_acc}"] = ref_gene
    
    def download_resistance_gene_references(self):
        """Download reference sequences for plasmid-mediated resistance genes"""
        logger.info("Downloading reference sequences for resistance genes...")
        
        genes = self.load_resistance_genes()
        processed_accessions = set()
        
        for gene_info in genes:
            protein_acc = gene_info['protein_accession']
            
            # Skip if already processed
            if protein_acc in processed_accessions:
                continue
                
            processed_accessions.add(protein_acc)
            
            logger.info(f"Processing {gene_info['species']} {gene_info['gene']} ({protein_acc})")
            
            ref_gene = self.process_protein_accession(
                protein_acc,
                gene_info['species'],
                gene_info['gene'],
                gene_info
            )
            
            if ref_gene:
                key = f"{gene_info['species']}_{gene_info['gene']}_{protein_acc}"
                self.reference_genes[key] = ref_gene
    
    def create_gpu_index(self) -> Dict:
        """Create GPU-optimized index structures"""
        logger.info("Creating GPU-optimized index...")
        
        # Collect all sequences
        sequences = []
        sequence_info = []
        
        for key, ref_gene in self.reference_genes.items():
            # Add full sequence (with flanking)
            sequences.append(ref_gene.full_sequence)
            sequence_info.append({
                'key': key,
                'species': ref_gene.species,
                'gene': ref_gene.gene_name,
                'gene_start_in_seq': len(ref_gene.flanking_5prime),
                'gene_end_in_seq': len(ref_gene.flanking_5prime) + len(ref_gene.nucleotide_sequence),
                'mutations': ref_gene.mutations
            })
        
        # Create 2-bit encoded sequences for GPU
        encoded_sequences = []
        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        for seq in sequences:
            encoded = []
            for i in range(0, len(seq), 16):  # Pack 16 bases into 32 bits
                chunk = seq[i:i+16]
                packed = 0
                for j, base in enumerate(chunk):
                    if base in encoding:
                        packed |= (encoding[base] << (2 * j))
                encoded.append(packed)
            encoded_sequences.append(np.array(encoded, dtype=np.uint32))
        
        return {
            'sequences': sequences,
            'encoded_sequences': encoded_sequences,
            'sequence_info': sequence_info,
            'num_sequences': len(sequences)
        }
    
    def save_database(self, output_dir: str):
        """Save the reference database"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save reference genes
        with open(output_path / 'reference_genes.pkl', 'wb') as f:
            pickle.dump(self.reference_genes, f)
        
        # Save GPU index
        gpu_index = self.create_gpu_index()
        with open(output_path / 'gpu_index.pkl', 'wb') as f:
            pickle.dump(gpu_index, f)
        
        # Save sequences in FASTA format for verification
        with open(output_path / 'reference_sequences.fasta', 'w') as f:
            for key, ref_gene in self.reference_genes.items():
                f.write(f">{key} {ref_gene.species} {ref_gene.gene_name}\n")
                f.write(f"{ref_gene.full_sequence}\n")
        
        # Save failed downloads
        if self.failed_downloads:
            with open(output_path / 'failed_downloads.json', 'w') as f:
                json.dump(self.failed_downloads, f, indent=2)
        
        # Save summary statistics
        stats = {
            'total_references': len(self.reference_genes),
            'failed_downloads': len(self.failed_downloads),
            'species': len(set(r.species for r in self.reference_genes.values())),
            'genes_per_type': {}
        }
        
        for gene_type in self.target_genes:
            count = sum(1 for r in self.reference_genes.values() if r.gene_name == gene_type)
            stats['genes_per_type'][gene_type] = count
        
        with open(output_path / 'database_summary.json', 'w') as f:
            json.dump(stats, f, indent=2)
        
        logger.info(f"Database saved to {output_path}")
        logger.info(f"Total references: {stats['total_references']}")
        logger.info(f"Failed downloads: {stats['failed_downloads']}")

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Download fluoroquinolone resistance gene references')
    parser.add_argument('--mutations', type=str, required=True, 
                       help='Path to mutations CSV file')
    parser.add_argument('--genes', type=str, required=True,
                       help='Path to resistance genes CSV file')
    parser.add_argument('--output', type=str, default='data/fq_reference_db',
                       help='Output directory')
    parser.add_argument('--flanking', type=int, default=500,
                       help='Size of flanking regions to include')
    parser.add_argument('--email', type=str, required=True,
                       help='Your email for NCBI Entrez')
    
    args = parser.parse_args()
    
    # Set email for Entrez
    Entrez.email = args.email
    
    # Create downloader
    downloader = ResistanceGeneDownloader(
        args.mutations,
        args.genes,
        flanking_size=args.flanking
    )
    
    # Download references
    downloader.download_mutation_gene_references()
    downloader.download_resistance_gene_references()
    
    # Save database
    downloader.save_database(args.output)
    
    print("\n=== Download Summary ===")
    print(f"Successfully downloaded: {len(downloader.reference_genes)} references")
    print(f"Failed downloads: {len(downloader.failed_downloads)}")
    
    if downloader.failed_downloads:
        print("\nFailed downloads:")
        for fail in downloader.failed_downloads[:5]:  # Show first 5
            print(f"  - {fail['species']} {fail['gene']} ({fail['accession']}): {fail['error']}")
        if len(downloader.failed_downloads) > 5:
            print(f"  ... and {len(downloader.failed_downloads) - 5} more")

if __name__ == "__main__":
    main()