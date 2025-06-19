#!/usr/bin/env python3
"""
Microbiome Type Strain Genome Downloader
Downloads bacterial type strain genomes from NCBI and extracts taxonomic information
"""

import os
import time
import requests
import pandas as pd
from Bio import Entrez
from urllib.parse import urljoin
import xml.etree.ElementTree as ET
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MicrobiomeGenomeDownloader:
    def __init__(self, email, api_key=None, output_dir="microbiome_genomes"):
        """
        Initialize the genome downloader
        
        Args:
            email (str): Your email for NCBI Entrez
            api_key (str, optional): NCBI API key for increased rate limits
            output_dir (str): Directory to save downloaded genomes
        """
        self.email = email
        self.api_key = api_key
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Configure Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            
        # Common gut microbiome bacteria (type strains) - Updated with better availability
        self.target_organisms = [
            "Bacteroides fragilis NCTC 9343",
            "Bacteroides thetaiotaomicron VPI-5482", 
            "Escherichia coli str. K-12 substr. DH10B",  # Better availability than MG1655
            "Enterococcus faecalis V583",
            "Lactobacillus acidophilus NCFM",
            "Bifidobacterium longum NCC2705",
            "Clostridium difficile 630",
            "Akkermansia muciniphila ATCC BAA-835",
            "Prevotella melaninogenica ATCC 25845",
            "Faecalibacterium prausnitzii A2-165",
            "Roseburia intestinalis L1-82",
            "Eubacterium rectale ATCC 33656",
            "Blautia obeum ATCC 29174",
            "Ruminococcus bromii ATCC 27255",  # Better availability than L2-63
            "Coprococcus eutactus ATCC 27759",
            "Dorea formicigenerans ATCC 27755",
            "Lachnospira pectinoschiza ATCC 51485",  # Better availability than DSM 2713
            "Veillonella parvula DSM 2008",
            "Streptococcus thermophilus LMG 18311",
            "Enterobacter cloacae subsp. cloacae ATCC 13047",  # More specific designation
            "Klebsiella pneumoniae subsp. pneumoniae ATCC 13883",  # More specific designation
            "Citrobacter freundii ATCC 8090",
            "Proteus mirabilis ATCC 7002",
            "Pseudomonas aeruginosa PAO1",
            "Staphylococcus epidermidis ATCC 12228",
            "Enterococcus faecium DO",
            "Lactobacillus rhamnosus GG",
            "Bifidobacterium bifidum PRL2010",
            "Clostridium perfringens ATCC 13124",
            # Additional well-sequenced alternatives
            "Escherichia coli CFT073",  # Well-sequenced uropathogenic strain
            "Klebsiella pneumoniae HS11286",  # Well-sequenced clinical isolate
            "Enterobacter cloacae EcWSU1"  # Alternative strain with good assembly
        ]
        
        self.genome_data = []
        
    def search_genome_assembly(self, organism):
        """Search for genome assembly for a given organism"""
        try:
            # Search for assembly
            search_term = f'"{organism}"[Organism] AND "reference genome"[Filter]'
            logger.info(f"Searching for: {search_term}")
            
            handle = Entrez.esearch(db="assembly", term=search_term, retmax=5)
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results["IdList"]:
                # Try without reference genome filter
                search_term = f'"{organism}"[Organism] AND "complete genome"[Assembly Level]'
                handle = Entrez.esearch(db="assembly", term=search_term, retmax=5)
                search_results = Entrez.read(handle)
                handle.close()
                
            if not search_results["IdList"]:
                # Try even broader search
                search_term = f'"{organism.split()[0]} {organism.split()[1]}"[Organism]'
                handle = Entrez.esearch(db="assembly", term=search_term, retmax=10)
                search_results = Entrez.read(handle)
                handle.close()
            
            return search_results["IdList"]
            
        except Exception as e:
            logger.error(f"Error searching for {organism}: {e}")
            return []
    
    def get_assembly_info(self, assembly_id):
        """Get detailed information about an assembly"""
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            result = Entrez.read(handle)
            handle.close()
            
            # Handle the nested structure of NCBI API response
            if 'DocumentSummarySet' in result:
                doc_sum_set = result['DocumentSummarySet']
                if 'DocumentSummary' in doc_sum_set:
                    summaries = doc_sum_set['DocumentSummary']
                    if isinstance(summaries, list) and len(summaries) > 0:
                        summary = summaries[0]
                    else:
                        summary = summaries
                else:
                    logger.error(f"No DocumentSummary found for {assembly_id}")
                    return None
            else:
                # Fallback to old structure if it exists
                if isinstance(result, list) and len(result) > 0:
                    summary = result[0]
                else:
                    logger.error(f"Unexpected response structure for {assembly_id}")
                    return None
            
            return {
                'assembly_id': assembly_id,
                'accession': summary.get('AssemblyAccession', ''),
                'organism': summary.get('Organism', ''),
                'assembly_name': summary.get('AssemblyName', ''),
                'taxid': summary.get('Taxid', ''),
                'species_taxid': summary.get('SpeciesTaxid', ''),
                'ftp_path': summary.get('FtpPath_RefSeq', summary.get('FtpPath_GenBank', '')),
                'assembly_level': summary.get('AssemblyLevel', ''),
                'assembly_status': summary.get('AssemblyStatus', '')
            }
            
        except Exception as e:
            logger.error(f"Error getting assembly info for {assembly_id}: {e}")
            return None
    
    def download_genome(self, assembly_info):
        """Download genome FASTA file"""
        try:
            if not assembly_info['ftp_path']:
                logger.warning(f"No FTP path for {assembly_info['organism']}")
                return None
                
            # Construct download URL
            ftp_path = assembly_info['ftp_path']
            filename = os.path.basename(ftp_path) + "_genomic.fna.gz"
            download_url = urljoin(ftp_path + "/", filename)
            
            # Create safe filename
            safe_organism = "".join(c for c in assembly_info['organism'] if c.isalnum() or c in (' ', '-', '_')).rstrip()
            safe_organism = safe_organism.replace(' ', '_')
            output_filename = f"{assembly_info['accession']}_{safe_organism}_genomic.fna.gz"
            output_path = self.output_dir / output_filename
            
            logger.info(f"Downloading {assembly_info['organism']} from {download_url}")
            
            # Use urllib for FTP downloads
            import urllib.request
            urllib.request.urlretrieve(download_url, output_path)
                    
            logger.info(f"Downloaded: {output_filename}")
            
            # Add local path to assembly info
            assembly_info['local_path'] = str(output_path)
            return output_path
            
        except Exception as e:
            logger.error(f"Error downloading genome for {assembly_info['organism']}: {e}")
            return None
    
    def get_taxonomic_lineage(self, taxid):
        """Get full taxonomic lineage for a taxid"""
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            if records:
                lineage = records[0].get('LineageEx', [])
                lineage_dict = {}
                for item in lineage:
                    rank = item.get('Rank', '')
                    name = item.get('ScientificName', '')
                    if rank in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                        lineage_dict[rank] = name
                
                # Add the organism itself if it's at species level
                if records[0].get('Rank') == 'species':
                    lineage_dict['species'] = records[0].get('ScientificName', '')
                
                return lineage_dict
            
        except Exception as e:
            logger.error(f"Error getting lineage for taxid {taxid}: {e}")
            
        return {}
    
    def process_all_organisms(self):
        """Process all target organisms"""
        logger.info(f"Processing {len(self.target_organisms)} target organisms")
        
        for i, organism in enumerate(self.target_organisms, 1):
            logger.info(f"Processing {i}/{len(self.target_organisms)}: {organism}")
            
            # Search for assemblies
            assembly_ids = self.search_genome_assembly(organism)
            
            if not assembly_ids:
                logger.warning(f"No assemblies found for {organism}")
                continue
            
            # Get info for the first (best) assembly
            assembly_info = self.get_assembly_info(assembly_ids[0])
            
            if not assembly_info:
                continue
            
            # Download genome
            local_path = self.download_genome(assembly_info)
            
            # Get taxonomic information
            lineage = {}
            if assembly_info['taxid']:
                lineage = self.get_taxonomic_lineage(assembly_info['taxid'])
            
            # Combine all information
            genome_record = {
                **assembly_info,
                **lineage,
                'search_organism': organism,
                'download_success': local_path is not None
            }
            
            self.genome_data.append(genome_record)
            
            # Be nice to NCBI servers
            time.sleep(0.5)
        
        return self.genome_data
    
    def save_metadata(self):
        """Save metadata to CSV file"""
        if not self.genome_data:
            logger.warning("No genome data to save")
            return
        
        df = pd.DataFrame(self.genome_data)
        
        # Reorder columns for better readability
        columns_order = [
            'search_organism', 'organism', 'accession', 'assembly_name',
            'taxid', 'species_taxid', 'superkingdom', 'phylum', 'class', 
            'order', 'family', 'genus', 'species', 'assembly_level',
            'assembly_status', 'download_success', 'local_path', 'ftp_path'
        ]
        
        # Only include columns that exist
        existing_columns = [col for col in columns_order if col in df.columns]
        df = df[existing_columns]
        
        output_file = self.output_dir / "genome_metadata.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Metadata saved to: {output_file}")
        
        # Create a simple taxid list
        taxid_file = self.output_dir / "species_taxids.txt"
        with open(taxid_file, 'w') as f:
            f.write("# Species TaxIDs from downloaded genomes\n")
            for _, row in df.iterrows():
                if pd.notna(row.get('species_taxid')):
                    f.write(f"{row['species_taxid']}\t{row.get('species', row.get('organism', ''))}\n")
        
        logger.info(f"TaxID list saved to: {taxid_file}")
        
        return df

def main():
    # Configuration - REPLACE WITH YOUR EMAIL
    EMAIL = "dbhaslam@gmail.com"  # REQUIRED: Replace with your email
    API_KEY = None  # Optional: Add your NCBI API key for faster downloads
    
    if EMAIL == "your.email@institution.edu":
        print("Please replace EMAIL with your actual email address before running!")
        return
    
    # Initialize downloader
    downloader = MicrobiomeGenomeDownloader(
        email=EMAIL,
        api_key=API_KEY,
        output_dir="data/type_strain_reference_genomes"
    )
    
    # Process all organisms
    logger.info("Starting genome download process...")
    genome_data = downloader.process_all_organisms()
    
    # Save metadata
    df = downloader.save_metadata()
    
    # Summary
    successful_downloads = sum(1 for record in genome_data if record.get('download_success'))
    logger.info(f"Download complete: {successful_downloads}/{len(genome_data)} genomes downloaded successfully")
    
    if df is not None:
        print(f"\nSummary:")
        print(f"Total organisms processed: {len(df)}")
        print(f"Successful downloads: {successful_downloads}")
        print(f"Output directory: {downloader.output_dir}")
        print(f"Metadata file: {downloader.output_dir}/genome_metadata.csv")
        print(f"TaxID file: {downloader.output_dir}/species_taxids.txt")

if __name__ == "__main__":
    main()