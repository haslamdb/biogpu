#!/usr/bin/env python3
"""
Process compressed .fna.gz files for Kraken2 database construction.

This script:
1. Finds all .fna.gz files in nested directories
2. Extracts assembly_accession from filenames
3. Looks up taxid and metadata from assembly_summary.txt
4. Decompresses, modifies headers with assembly_accession and taxid
5. Saves processed files to output directory

Usage:
    python process_fna_db_build_compressed.py --library_dir /path/to/library --output_dir processed_fna
"""

import os
import re
import sys
import gzip
import shutil
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FNAProcessor:
    def __init__(self):
        """Initialize processor."""
        self.processed_count = 0
        self.error_count = 0
        self.missing_taxid_count = 0
        self.assembly_summary = {}
    
    def load_assembly_summary(self, summary_path: str) -> Dict[str, Dict]:
        """Load assembly summary file and create lookup dictionary."""
        logger.info(f"Loading assembly summary from {summary_path}")
        
        try:
            # Read the assembly summary file, skipping comment lines
            df = pd.read_csv(summary_path, sep='\t', comment='#', 
                           dtype={'taxid': str, 'species_taxid': str})
            
            # Create lookup dictionary keyed by assembly_accession
            lookup = {}
            for _, row in df.iterrows():
                accession = row['# assembly_accession'] if '# assembly_accession' in row else row.iloc[0]
                taxid = row['taxid'] if 'taxid' in row else row.iloc[5]
                
                lookup[accession] = {
                    'taxid': str(taxid),
                    'organism_name': row['organism_name'] if 'organism_name' in row else row.iloc[7]
                }
            
            logger.info(f"Loaded {len(lookup)} assembly records")
            return lookup
            
        except Exception as e:
            logger.error(f"Error loading assembly summary: {e}")
            return {}
    
    def extract_assembly_accession(self, filepath: str) -> Optional[str]:
        """Extract assembly accession from filepath."""
        # Pattern: GCF_XXXXXXXXX.X or GCA_XXXXXXXXX.X
        filename = os.path.basename(filepath)
        match = re.search(r'(GC[AF]_\d+\.\d+)', filename)
        if match:
            return match.group(1)
        return None
    
    def find_fna_files(self, input_dir: str) -> list:
        """Find all .fna.gz files in directory and subdirectories."""
        fna_files = []
        input_path = Path(input_dir)
        
        # Look for both .fna and .fna.gz files
        for pattern in ['*.fna', '*.fna.gz']:
            for fna_file in input_path.rglob(pattern):
                if fna_file.is_file() and 'genomic' in fna_file.name:
                    fna_files.append(fna_file)
        
        logger.info(f"Found {len(fna_files)} genomic .fna(.gz) files")
        return fna_files
    
    def process_fna_header(self, header_line: str, assembly_accession: str, 
                          taxid: str) -> str:
        """Modify header to include assembly_accession and taxid."""
        # Remove the '>' character for processing
        original_header = header_line.lstrip('>')
        
        # Create new header: >assembly_accession|taxid|original_header
        new_header = f">{assembly_accession}|{taxid}|{original_header}"
        
        return new_header
    
    def process_single_fna(self, fna_path: Path, output_dir: str) -> bool:
        """Process a single .fna or .fna.gz file."""
        assembly_accession = self.extract_assembly_accession(str(fna_path))
        
        if not assembly_accession:
            logger.warning(f"Could not extract assembly accession from {fna_path}")
            self.error_count += 1
            return False
        
        # Look up taxid
        if assembly_accession not in self.assembly_summary:
            logger.warning(f"Assembly accession {assembly_accession} not found in summary")
            self.missing_taxid_count += 1
            return False
        
        metadata = self.assembly_summary[assembly_accession]
        taxid = metadata['taxid']
        
        if not taxid or taxid == 'na':
            logger.warning(f"No valid taxid for {assembly_accession}")
            self.missing_taxid_count += 1
            return False
        
        # Prepare output filename (always .fna, not .fna.gz)
        output_filename = fna_path.name.replace('.gz', '') if fna_path.name.endswith('.gz') else fna_path.name
        output_path = Path(output_dir) / output_filename
        
        try:
            # Open input file (compressed or not)
            if fna_path.name.endswith('.gz'):
                infile = gzip.open(fna_path, 'rt')
            else:
                infile = open(fna_path, 'r')
            
            with infile:
                with open(output_path, 'w') as outfile:
                    for line in infile:
                        if line.startswith('>'):
                            # Modify header
                            new_header = self.process_fna_header(line.strip(), 
                                                               assembly_accession, taxid)
                            outfile.write(new_header + '\n')
                        else:
                            # Copy sequence lines as-is
                            outfile.write(line)
            
            self.processed_count += 1
            if self.processed_count % 100 == 0:
                logger.info(f"Processed {self.processed_count} files...")
            
            return True
            
        except Exception as e:
            logger.error(f"Error processing {fna_path}: {e}")
            self.error_count += 1
            return False
    
    def process_library_subdirectories(self, library_dir: str, output_dir: str):
        """Process all subdirectories in library directory (bacteria, fungi, etc.)."""
        # Create main output directory
        os.makedirs(output_dir, exist_ok=True)
        
        library_path = Path(library_dir)
        if not library_path.exists():
            logger.error(f"Library directory does not exist: {library_dir}")
            return
        
        # Find all subdirectories that contain assembly_summary.txt
        subdirs_to_process = []
        for subdir in library_path.iterdir():
            if subdir.is_dir():
                assembly_summary_path = subdir / "assembly_summary.txt"
                if assembly_summary_path.exists():
                    subdirs_to_process.append(subdir)
                    logger.info(f"Found subdirectory to process: {subdir.name}")
        
        if not subdirs_to_process:
            logger.error("No subdirectories with assembly_summary.txt found")
            return
        
        total_stats = {
            'processed': 0,
            'errors': 0,
            'missing_taxids': 0,
            'total_files': 0
        }
        
        # Process each subdirectory
        for subdir in subdirs_to_process:
            logger.info(f"\n=== Processing subdirectory: {subdir.name} ===")
            
            # Create subdirectory-specific output directory
            subdir_output = Path(output_dir) / subdir.name
            os.makedirs(subdir_output, exist_ok=True)
            
            # Reset counters for this subdirectory
            self.processed_count = 0
            self.error_count = 0
            self.missing_taxid_count = 0
            
            # Load assembly summary for this subdirectory
            assembly_summary_path = subdir / "assembly_summary.txt"
            self.assembly_summary = self.load_assembly_summary(str(assembly_summary_path))
            
            if not self.assembly_summary:
                logger.warning(f"Failed to load assembly summary for {subdir.name}")
                continue
            
            # Find and process .fna files in this subdirectory
            fna_files = self.find_fna_files(str(subdir))
            
            if not fna_files:
                logger.warning(f"No genomic .fna(.gz) files found in {subdir.name}")
                continue
            
            logger.info(f"Processing {len(fna_files)} .fna(.gz) files in {subdir.name}...")
            
            # Process each file
            for fna_path in fna_files:
                self.process_single_fna(fna_path, str(subdir_output))
            
            # Report results for this subdirectory
            logger.info(f"Subdirectory {subdir.name} complete:")
            logger.info(f"  Successfully processed: {self.processed_count}")
            logger.info(f"  Errors: {self.error_count}")
            logger.info(f"  Missing taxids: {self.missing_taxid_count}")
            logger.info(f"  Total files: {len(fna_files)}")
            
            # Add to total stats
            total_stats['processed'] += self.processed_count
            total_stats['errors'] += self.error_count
            total_stats['missing_taxids'] += self.missing_taxid_count
            total_stats['total_files'] += len(fna_files)
        
        # Report overall results
        logger.info("\n=== OVERALL PROCESSING COMPLETE ===")
        logger.info(f"Total successfully processed: {total_stats['processed']}")
        logger.info(f"Total errors: {total_stats['errors']}")
        logger.info(f"Total missing taxids: {total_stats['missing_taxids']}")
        logger.info(f"Total files: {total_stats['total_files']}")
        
        # Update instance variables for final stats
        self.processed_count = total_stats['processed']
        self.error_count = total_stats['errors']
        self.missing_taxid_count = total_stats['missing_taxids']
    
    def create_database_stats(self, output_dir: str) -> str:
        """Create a summary of processed files for database construction."""
        stats_file = os.path.join(output_dir, "database_stats.txt")
        
        with open(stats_file, 'w') as f:
            f.write("FNA Database Processing Summary\n")
            f.write("================================\n\n")
            f.write(f"Successfully processed files: {self.processed_count}\n")
            f.write(f"Files with errors: {self.error_count}\n")
            f.write(f"Files missing taxids: {self.missing_taxid_count}\n")
            f.write(f"\nProcessed subdirectories:\n")
            
            # List subdirectories in output
            output_path = Path(output_dir)
            for subdir in output_path.iterdir():
                if subdir.is_dir():
                    file_count = len(list(subdir.glob("*.fna")))
                    f.write(f"  {subdir.name}: {file_count} files\n")
            
            f.write(f"\nReady for Kraken2 database construction\n")
            f.write(f"Output directory: {output_dir}\n")
        
        return stats_file

def main():
    parser = argparse.ArgumentParser(description="Process .fna.gz files for Kraken2 database")
    parser.add_argument("--library_dir", required=True, 
                       help="Library directory containing subdirectories (bacteria, fungi, etc.) with assembly_summary.txt files")
    parser.add_argument("--output_dir", required=True,
                       help="Output directory for processed .fna files")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate inputs
    if not os.path.exists(args.library_dir):
        logger.error(f"Library directory does not exist: {args.library_dir}")
        sys.exit(1)
    
    # Process files
    processor = FNAProcessor()
    processor.process_library_subdirectories(args.library_dir, args.output_dir)
    
    # Create summary stats
    if processor.processed_count > 0:
        stats_file = processor.create_database_stats(args.output_dir)
        logger.info(f"Database stats written to: {stats_file}")

if __name__ == "__main__":
    main()