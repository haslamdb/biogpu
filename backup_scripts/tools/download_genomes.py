#!/usr/bin/env python3
"""
Download reference genomes from NCBI RefSeq with proper organism mapping.
Downloads genomic.fna.gz files and maintains metadata for organism identification.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import json
import pandas as pd
from typing import List, Dict, Optional, Tuple
import ftplib
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlparse
import time

class NCBIGenomeDownloader:
    """Download genomes from NCBI RefSeq with metadata tracking."""
    
    def __init__(self, project_root: str = None, max_workers: int = 4):
        """
        Initialize the downloader.
        
        Args:
            project_root: Root directory of the project
            max_workers: Maximum number of concurrent downloads
        """
        self.project_root = Path(project_root) if project_root else Path.cwd()
        self.data_dir = self.project_root / "data" / "genomes"
        self.max_workers = max_workers
        
        # NCBI FTP settings
        self.ftp_host = "ftp.ncbi.nlm.nih.gov"
        self.libraries = {
            "fungi": "/genomes/refseq/fungi/",
            "bacteria": "/genomes/refseq/bacteria/",
            "viral": "/genomes/refseq/viral/", 
            "plasmid": "/genomes/refseq/plasmid/"
        }
        
        # Setup logging
        self.setup_logging()
        
        # Initialize FTP connection
        self.ftp = None
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_dir = self.project_root / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = log_dir / f"genome_download_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def connect_ftp(self):
        """Establish FTP connection."""
        try:
            self.ftp = ftplib.FTP(self.ftp_host)
            self.ftp.login()
            self.logger.info(f"Connected to {self.ftp_host}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to connect to FTP: {e}")
            return False
            
    def disconnect_ftp(self):
        """Close FTP connection."""
        if self.ftp:
            try:
                self.ftp.quit()
            except:
                pass
                
    def download_assembly_summary(self, library: str) -> Optional[pd.DataFrame]:
        """
        Download and parse the main assembly_summary.txt for a library.
        
        Returns:
            DataFrame with assembly information or None if failed
        """
        summary_path = self.libraries[library] + "assembly_summary.txt"
        local_file = self.data_dir / library / "assembly_summary.txt"
        local_file.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # Download using wget for reliability
            wget_cmd = [
                "wget",
                "-q",
                f"ftp://{self.ftp_host}{summary_path}",
                "-O", str(local_file)
            ]
            subprocess.run(wget_cmd, check=True)
            
            # Parse the assembly summary (skip comment lines)
            # First, count comment lines
            with open(local_file, 'r') as f:
                skip_lines = 0
                for line in f:
                    if line.startswith('#'):
                        skip_lines += 1
                    else:
                        break
            
            # Read the file, skipping comment lines (skip_lines already includes the header comment)
            df = pd.read_csv(local_file, sep='\t', skiprows=skip_lines-1)
            self.logger.info(f"Downloaded assembly summary for {library}: {len(df)} assemblies")
            
            # Remove '#' from first column name if present
            if df.columns[0].startswith('#'):
                df.columns = [df.columns[0].lstrip('#')] + list(df.columns[1:])
            
            self.logger.debug(f"DataFrame columns: {list(df.columns[:5])}...")  # Show first 5 columns
            
            # Add library column
            df['library'] = library
            
            return df
            
        except Exception as e:
            self.logger.error(f"Failed to download assembly summary for {library}: {e}")
            return None
            
    def get_genome_info(self, assembly_df: pd.DataFrame, assembly_level: List[str] = None, reference_only: bool = True) -> List[Dict]:
        """
        Extract genome download information from assembly summary.
        
        Args:
            assembly_df: DataFrame from assembly_summary.txt
            assembly_level: List of assembly levels to include 
                          (e.g., ['Complete Genome', 'Chromosome', 'Scaffold'])
            reference_only: Download only reference genomes (default: True)
                          
        Returns:
            List of dictionaries with genome information
        """
        genomes = []
        
        # Filter for reference genomes only if requested
        if reference_only:
            # Check if there are any reference genomes in this library
            reference_df = assembly_df[assembly_df['refseq_category'] == 'reference genome']
            if len(reference_df) > 0:
                working_df = reference_df
                self.logger.info(f"Found {len(working_df)} reference genomes out of {len(assembly_df)} total assemblies")
            else:
                # For libraries without reference genomes (like viral), use all genomes
                working_df = assembly_df
                self.logger.info(f"No reference genomes found in this library. Processing all {len(working_df)} assemblies")
        else:
            working_df = assembly_df
            self.logger.info(f"Processing all {len(working_df)} assemblies")
        
        # Filter by assembly level if specified
        if assembly_level:
            working_df = working_df[working_df['assembly_level'].isin(assembly_level)]
            self.logger.debug(f"After assembly level filter: {len(working_df)} rows")
            
        for _, row in working_df.iterrows():
            # Extract FTP path
            ftp_path = row.get('ftp_path', '')
            if not ftp_path or ftp_path == 'na':
                self.logger.debug(f"Skipping row with invalid FTP path: {ftp_path}")
                continue
                
            # Parse the path to get the genomic file
            # The ftp_path might be https:// or ftp:// - normalize it
            if ftp_path.startswith('https://'):
                path_parts = ftp_path.replace('https://ftp.ncbi.nlm.nih.gov', '')
            else:
                path_parts = ftp_path.replace('ftp://ftp.ncbi.nlm.nih.gov', '')
            assembly_accession = row['assembly_accession']  # Fixed column name
            
            genome_info = {
                'assembly_accession': assembly_accession,
                'organism_name': row.get('organism_name', 'Unknown'),
                'taxid': row.get('taxid', 'Unknown'),
                'species_taxid': row.get('species_taxid', 'Unknown'),
                'assembly_level': row.get('assembly_level', 'Unknown'),
                'genome_rep': row.get('genome_rep', 'Full'),
                'refseq_category': row.get('refseq_category', 'na'),
                'ftp_path': path_parts,
                'library': row['library'],
                'local_path': None,
                'download_status': 'pending'
            }
            
            # Construct the genomic file name
            if '_' in assembly_accession:
                genome_file = f"{assembly_accession}_{row.get('asm_name', 'unknown')}_genomic.fna.gz"
            else:
                genome_file = f"{assembly_accession}_genomic.fna.gz"
                
            genome_info['genome_file'] = genome_file
            # Keep the original URL format
            genome_info['full_ftp_path'] = f"{ftp_path}/{genome_file}"
            
            genomes.append(genome_info)
            
        return genomes
        
    def download_genome(self, genome_info: Dict) -> bool:
        """
        Download a single genome file.
        
        Args:
            genome_info: Dictionary with genome information
            
        Returns:
            True if successful, False otherwise
        """
        library = genome_info['library']
        organism_name = genome_info['organism_name'].replace(' ', '_').replace('/', '_')
        assembly_acc = genome_info['assembly_accession']
        
        # Create local directory structure
        local_dir = self.data_dir / library / organism_name / assembly_acc
        local_dir.mkdir(parents=True, exist_ok=True)
        
        # Local file path
        local_file = local_dir / genome_info['genome_file']
        genome_info['local_path'] = str(local_file)
        
        # Skip if already downloaded
        if local_file.exists():
            self.logger.info(f"Already downloaded: {organism_name} - {assembly_acc}")
            genome_info['download_status'] = 'exists'
            return True
            
        try:
            # Download using wget
            # Use the full_ftp_path as is (it already includes the protocol and domain)
            ftp_url = genome_info['full_ftp_path']
            wget_cmd = [
                "wget",
                "-q",
                "--show-progress",
                ftp_url,
                "-O", str(local_file)
            ]
            
            self.logger.info(f"Downloading: {organism_name} - {assembly_acc}")
            subprocess.run(wget_cmd, check=True)
            
            genome_info['download_status'] = 'completed'
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to download {organism_name}: {e}")
            genome_info['download_status'] = 'failed'
            
            # Remove partial file
            if local_file.exists():
                local_file.unlink()
                
            return False
            
    def download_library_genomes(self, library: str, assembly_level: List[str] = None,
                               limit: int = None, parallel: bool = True, reference_only: bool = True) -> pd.DataFrame:
        """
        Download all genomes for a library.
        
        Args:
            library: Library name (fungi, bacteria, viral, plasmid)
            assembly_level: Filter by assembly level
            limit: Maximum number of genomes to download
            parallel: Use parallel downloads
            reference_only: Download only reference genomes (default: True)
            
        Returns:
            DataFrame with download results
        """
        self.logger.info(f"\nProcessing {library} library...")
        
        # Get assembly summary
        assembly_df = self.download_assembly_summary(library)
        if assembly_df is None:
            return pd.DataFrame()
            
        # Get genome information
        genomes = self.get_genome_info(assembly_df, assembly_level, reference_only)
        self.logger.info(f"Found {len(genomes)} genomes to download")
        
        if limit:
            genomes = genomes[:limit]
            self.logger.info(f"Limited to {limit} genomes")
            
        # Download genomes
        if parallel and len(genomes) > 1:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {executor.submit(self.download_genome, genome): genome 
                          for genome in genomes}
                
                for future in as_completed(futures):
                    genome = futures[future]
                    try:
                        success = future.result()
                    except Exception as e:
                        self.logger.error(f"Download failed for {genome['organism_name']}: {e}")
        else:
            for genome in genomes:
                self.download_genome(genome)
                
        # Create results DataFrame
        results_df = pd.DataFrame(genomes)
        
        # Save metadata
        metadata_file = self.data_dir / library / "genome_metadata.csv"
        results_df.to_csv(metadata_file, index=False)
        self.logger.info(f"Saved metadata to {metadata_file}")
        
        # Summary statistics
        if len(results_df) > 0 and 'download_status' in results_df.columns:
            status_counts = results_df['download_status'].value_counts()
            self.logger.info(f"\nDownload summary for {library}:")
            for status, count in status_counts.items():
                self.logger.info(f"  {status}: {count}")
        else:
            self.logger.info(f"\nNo genomes processed for {library}")
            
        return results_df
        
    def create_contig_mapping(self, library: str) -> pd.DataFrame:
        """
        Create a mapping of contigs to organisms for a library.
        
        Returns:
            DataFrame with contig to organism mapping
        """
        metadata_file = self.data_dir / library / "genome_metadata.csv"
        if not metadata_file.exists():
            self.logger.error(f"No metadata found for {library}")
            return pd.DataFrame()
            
        # Load metadata
        metadata_df = pd.read_csv(metadata_file)
        
        contig_mapping = []
        
        for _, row in metadata_df.iterrows():
            if row['download_status'] != 'completed' and row['download_status'] != 'exists':
                continue
                
            genome_file = Path(row['local_path'])
            if not genome_file.exists():
                continue
                
            # Parse the genome file to get contig information
            try:
                with gzip.open(genome_file, 'rt') as f:
                    for line in f:
                        if line.startswith('>'):
                            contig_id = line.strip().split()[0][1:]  # Remove '>'
                            
                            contig_info = {
                                'contig_id': contig_id,
                                'assembly_accession': row['assembly_accession'],
                                'organism_name': row['organism_name'],
                                'taxid': row['taxid'],
                                'species_taxid': row['species_taxid'],
                                'assembly_level': row['assembly_level'],
                                'library': library,
                                'genome_file': str(genome_file)
                            }
                            contig_mapping.append(contig_info)
                            
            except Exception as e:
                self.logger.error(f"Failed to parse {genome_file}: {e}")
                
        # Create DataFrame
        contig_df = pd.DataFrame(contig_mapping)
        
        # Save contig mapping
        contig_file = self.data_dir / library / "contig_mapping.csv"
        contig_df.to_csv(contig_file, index=False)
        self.logger.info(f"Created contig mapping for {library}: {len(contig_df)} contigs")
        
        return contig_df
        
    def download_plasmid_files(self, limit: int = None) -> pd.DataFrame:
        """
        Download plasmid files directly from the plasmid directory.
        
        Args:
            limit: Maximum number of files to download
            
        Returns:
            DataFrame with download results
        """
        self.logger.info("\nProcessing plasmid library...")
        
        # Create local directory
        plasmid_dir = self.data_dir / "plasmid"
        plasmid_dir.mkdir(parents=True, exist_ok=True)
        
        # List files in the plasmid directory
        try:
            # Use wget to list directory contents
            list_cmd = [
                "wget", "-q", "-O-", 
                f"ftp://{self.ftp_host}/genomes/refseq/plasmid/"
            ]
            result = subprocess.run(list_cmd, capture_output=True, text=True, check=True)
            
            # Parse for genomic.fna.gz files
            import re
            pattern = r'href="[^"]*/(plasmid\.[0-9]+\.[0-9]+\.genomic\.fna\.gz)"'
            files = re.findall(pattern, result.stdout)
            
            self.logger.info(f"Found {len(files)} plasmid files to download")
            
            if limit and len(files) > limit:
                files = files[:limit]
                self.logger.info(f"Limited to {limit} files")
            
            # Download files
            plasmid_metadata = []
            for filename in files:
                local_file = plasmid_dir / filename
                
                # Skip if already downloaded
                if local_file.exists():
                    self.logger.info(f"Already downloaded: {filename}")
                    status = 'exists'
                else:
                    try:
                        # Download the file
                        wget_cmd = [
                            "wget", "-q", "--show-progress",
                            f"ftp://{self.ftp_host}/genomes/refseq/plasmid/{filename}",
                            "-O", str(local_file)
                        ]
                        self.logger.info(f"Downloading: {filename}")
                        subprocess.run(wget_cmd, check=True)
                        status = 'completed'
                    except Exception as e:
                        self.logger.error(f"Failed to download {filename}: {e}")
                        status = 'failed'
                        if local_file.exists():
                            local_file.unlink()
                
                # Add to metadata
                plasmid_metadata.append({
                    'filename': filename,
                    'library': 'plasmid',
                    'local_path': str(local_file),
                    'download_status': status,
                    'ftp_url': f"ftp://{self.ftp_host}/genomes/refseq/plasmid/{filename}"
                })
            
            # Create DataFrame and save metadata
            results_df = pd.DataFrame(plasmid_metadata)
            metadata_file = plasmid_dir / "plasmid_metadata.csv"
            results_df.to_csv(metadata_file, index=False)
            self.logger.info(f"Saved metadata to {metadata_file}")
            
            # Summary statistics
            if len(results_df) > 0:
                status_counts = results_df['download_status'].value_counts()
                self.logger.info(f"\nDownload summary for plasmids:")
                for status, count in status_counts.items():
                    self.logger.info(f"  {status}: {count}")
            
            return results_df
            
        except Exception as e:
            self.logger.error(f"Failed to list plasmid directory: {e}")
            return pd.DataFrame()
    
    def create_combined_metadata(self) -> pd.DataFrame:
        """
        Create a combined metadata file for all downloaded libraries.
        """
        all_metadata = []
        
        for library in self.libraries.keys():
            if library == 'plasmid':
                metadata_file = self.data_dir / library / "plasmid_metadata.csv"
            else:
                metadata_file = self.data_dir / library / "genome_metadata.csv"
                
            if metadata_file.exists():
                df = pd.read_csv(metadata_file)
                all_metadata.append(df)
                
        if all_metadata:
            combined_df = pd.concat(all_metadata, ignore_index=True)
            combined_file = self.data_dir / "all_genomes_metadata.csv"
            combined_df.to_csv(combined_file, index=False)
            self.logger.info(f"Created combined metadata: {len(combined_df)} genomes")
            return combined_df
        else:
            return pd.DataFrame()


def main():
    """Main function to handle command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download NCBI RefSeq genomes with metadata tracking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download reference bacterial genomes only (default behavior)
  python download_genomes.py --library bacteria --assembly-level "Complete Genome"
  
  # Download ALL bacterial genomes, not just reference
  python download_genomes.py --library bacteria --all-genomes
  
  # Download first 100 fungal reference genomes
  python download_genomes.py --library fungi --limit 100
  
  # Download all viral reference genomes
  python download_genomes.py --library viral
  
  # Download plasmid sequences (Note: plasmids are downloaded as bulk files)
  python download_genomes.py --library plasmid --limit 5
  
  # Create contig mapping for downloaded genomes
  python download_genomes.py --library bacteria --create-mapping
  
  # Download from multiple libraries (reference genomes only)
  python download_genomes.py --library bacteria fungi --assembly-level "Complete Genome" "Chromosome"
        """
    )
    
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="Project root directory (default: current directory)"
    )
    
    parser.add_argument(
        "--library",
        nargs="+",
        choices=["fungi", "bacteria", "viral", "plasmid"],
        required=True,
        help="Libraries to download from"
    )
    
    parser.add_argument(
        "--assembly-level",
        nargs="+",
        choices=["Complete Genome", "Chromosome", "Scaffold", "Contig"],
        help="Filter by assembly level (default: all levels)"
    )
    
    parser.add_argument(
        "--limit",
        type=int,
        help="Maximum number of genomes to download per library"
    )
    
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel download workers (default: 4)"
    )
    
    parser.add_argument(
        "--create-mapping",
        action="store_true",
        help="Create contig to organism mapping after download"
    )
    
    parser.add_argument(
        "--combine-metadata",
        action="store_true",
        help="Create combined metadata file for all libraries"
    )
    
    parser.add_argument(
        "--reference-only",
        action="store_true",
        default=True,
        help="Download only reference genomes (default: True)"
    )
    
    parser.add_argument(
        "--all-genomes",
        action="store_true",
        help="Download all genomes, not just reference genomes"
    )
    
    args = parser.parse_args()
    
    # Initialize downloader
    downloader = NCBIGenomeDownloader(
        project_root=args.project_root,
        max_workers=args.workers
    )
    
    # Determine if we should download only reference genomes
    reference_only = not args.all_genomes
    
    # Download genomes for each library
    for library in args.library:
        if library == 'plasmid':
            # Special handling for plasmids
            results = downloader.download_plasmid_files(limit=args.limit)
        else:
            results = downloader.download_library_genomes(
                library=library,
                assembly_level=args.assembly_level,
                limit=args.limit,
                parallel=args.workers > 1,
                reference_only=reference_only
            )
        
        # Create contig mapping if requested (skip for plasmids)
        if args.create_mapping and not results.empty and library != 'plasmid':
            downloader.create_contig_mapping(library)
            
    # Create combined metadata if requested
    if args.combine_metadata:
        downloader.create_combined_metadata()
        
    downloader.logger.info("\nDownload complete!")


if __name__ == "__main__":
    main()