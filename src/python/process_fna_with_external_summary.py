#!/usr/bin/env python3
"""
Process compressed .fna.gz files with external assembly_summary.txt file.

This script:
1. Finds all .fna.gz files in the input directory
2. Uses an external assembly_summary.txt file for species_taxid lookup
3. Modifies headers with assembly_accession and species_taxid
4. Saves processed files to output directory

Usage:
    python process_fna_with_external_summary.py --library_dir /path/to/fna/files --assembly_summary /path/to/assembly_summary.txt --output_dir processed_fna
"""

import os
import re
import sys
import gzip
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, Optional
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FNAProcessor:
    def __init__(self, assembly_summary_path: str):
        """Initialize processor with external assembly summary."""
        self.processed_count = 0
        self.error_count = 0
        self.missing_taxid_count = 0
        self.assembly_summary = self.load_assembly_summary(assembly_summary_path)
    
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
                # Handle different column name formats
                if '# assembly_accession' in df.columns:
                    accession = row['# assembly_accession']
                elif 'assembly_accession' in df.columns:
                    accession = row['assembly_accession']
                else:
                    accession = row.iloc[0]  # First column
                
                # Get species_taxid (preferred) or regular taxid
                if 'species_taxid' in df.columns:
                    taxid = row['species_taxid']
                elif 'taxid' in df.columns:
                    taxid = row['taxid']
                else:
                    taxid = row.iloc[6]  # Usually 7th column for species_taxid
                
                # Get organism name
                if 'organism_name' in df.columns:
                    organism_name = row['organism_name']
                else:
                    organism_name = row.iloc[7] if len(row) > 7 else 'Unknown'
                
                lookup[accession] = {
                    'taxid': str(taxid),
                    'organism_name': organism_name
                }
            
            logger.info(f"Loaded {len(lookup)} assembly records (using species_taxid)")
            return lookup
            
        except Exception as e:
            logger.error(f"Error loading assembly summary: {e}")
            sys.exit(1)
    
    def extract_assembly_accession(self, filepath: str) -> Optional[str]:
        """Extract assembly accession from filepath."""
        # Pattern: GCF_XXXXXXXXX.X or GCA_XXXXXXXXX.X
        filename = os.path.basename(filepath)
        match = re.search(r'(GC[AF]_\d+\.\d+)', filename)
        if match:
            return match.group(1)
        return None
    
    def find_fna_files(self, input_dir: str) -> list:
        """Find all .fna.gz and .fna files in directory and subdirectories."""
        fna_files = []
        input_path = Path(input_dir)
        
        # Look for both .fna and .fna.gz files
        for pattern in ['*.fna', '*.fna.gz', '*.fasta', '*.fasta.gz', '*.fa', '*.fa.gz']:
            for fna_file in input_path.rglob(pattern):
                if fna_file.is_file():
                    fna_files.append(fna_file)
        
        logger.info(f"Found {len(fna_files)} FASTA files")
        return fna_files
    
    def process_fna_header(self, header_line: str, assembly_accession: str, 
                          taxid: str) -> str:
        """Modify header to include assembly_accession and taxid."""
        # Remove the '>' character for processing
        original_header = header_line.lstrip('>')
        
        # Create new header: >assembly_accession|species_taxid|original_header
        new_header = f">{assembly_accession}|{taxid}|{original_header}"
        
        return new_header
    
    def process_single_fna(self, fna_path: Path, output_dir: str) -> bool:
        """Process a single .fna or .fna.gz file."""
        assembly_accession = self.extract_assembly_accession(str(fna_path))
        
        if not assembly_accession:
            logger.warning(f"Could not extract assembly accession from {fna_path.name}")
            self.error_count += 1
            return False
        
        # Look up taxid
        if assembly_accession not in self.assembly_summary:
            logger.warning(f"Assembly accession {assembly_accession} not found in summary")
            self.missing_taxid_count += 1
            return False
        
        metadata = self.assembly_summary[assembly_accession]
        taxid = metadata['taxid']
        
        if not taxid or taxid == 'na' or taxid == 'NA':
            logger.warning(f"No valid taxid for {assembly_accession}")
            self.missing_taxid_count += 1
            return False
        
        # Preserve directory structure in output
        relative_path = fna_path.relative_to(Path(output_dir).parent.parent) if fna_path.is_absolute() else fna_path
        output_path = Path(output_dir) / relative_path.parent.name / fna_path.name
        
        # Remove .gz extension if present
        if output_path.name.endswith('.gz'):
            output_path = output_path.with_suffix('')
        
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
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
    
    def process_all_files(self, library_dir: str, output_dir: str):
        """Process all FASTA files in library directory."""
        # Create main output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Find all FASTA files
        fna_files = self.find_fna_files(library_dir)
        
        if not fna_files:
            logger.error(f"No FASTA files found in {library_dir}")
            return
        
        logger.info(f"Processing {len(fna_files)} files...")
        
        # Process each file
        for fna_path in fna_files:
            self.process_single_fna(fna_path, output_dir)
        
        # Report results
        logger.info("\n=== PROCESSING COMPLETE ===")
        logger.info(f"Successfully processed: {self.processed_count}")
        logger.info(f"Errors: {self.error_count}")
        logger.info(f"Missing taxids: {self.missing_taxid_count}")
        logger.info(f"Total files: {len(fna_files)}")
    
    def create_database_stats(self, output_dir: str) -> str:
        """Create a summary of processed files."""
        stats_file = os.path.join(output_dir, "processing_stats.txt")
        
        with open(stats_file, 'w') as f:
            f.write("FNA Processing Summary\n")
            f.write("=====================\n\n")
            f.write(f"Successfully processed files: {self.processed_count}\n")
            f.write(f"Files with errors: {self.error_count}\n")
            f.write(f"Files missing taxids: {self.missing_taxid_count}\n")
            f.write(f"\nOutput directory: {output_dir}\n")
        
        return stats_file

def main():
    parser = argparse.ArgumentParser(description="Process .fna.gz files with external assembly_summary.txt")
    parser.add_argument("--library_dir", required=True, 
                       help="Directory containing .fna.gz files")
    parser.add_argument("--assembly_summary", required=True,
                       help="Path to assembly_summary.txt file")
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
    
    if not os.path.exists(args.assembly_summary):
        logger.error(f"Assembly summary file does not exist: {args.assembly_summary}")
        sys.exit(1)
    
    # Process files
    processor = FNAProcessor(args.assembly_summary)
    processor.process_all_files(args.library_dir, args.output_dir)
    
    # Create summary stats
    if processor.processed_count > 0:
        stats_file = processor.create_database_stats(args.output_dir)
        logger.info(f"Processing stats written to: {stats_file}")

if __name__ == "__main__":
    main()