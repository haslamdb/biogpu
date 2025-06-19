#!/usr/bin/env python3
"""
Process large library.fna files to inject species_taxid into headers.

This script:
1. Reads large library.fna files containing multiple genome sequences
2. Extracts assembly accession from each sequence header
3. Looks up species_taxid from assembly_summary.txt
4. Modifies headers to include species_taxid
5. Writes processed sequences to output file

Usage:
    python process_library_fna_species_taxid.py --library_fna /path/to/library.fna --assembly_summary /path/to/assembly_summary.txt --output_file processed_library.fna
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

class LibraryFNAProcessor:
    def __init__(self, assembly_summary_path: str):
        """Initialize processor with external assembly summary."""
        self.processed_sequences = 0
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
    
    def extract_assembly_accession(self, header_line: str) -> Optional[str]:
        """Extract assembly accession from header line."""
        # First try to find GCF/GCA pattern
        match = re.search(r'(GC[AF]_\d+\.\d+)', header_line)
        if match:
            return match.group(1)
        
        # If not found, try to extract from sequence ID (e.g., NZ_CM000441.1)
        # These RefSeq IDs can be mapped if they're in assembly_summary
        match = re.search(r'(NZ_[A-Z]+\d+\.\d+)', header_line)
        if match:
            # For RefSeq IDs, we need to search through our assembly data
            # This is a limitation - we can't directly map these
            return None
        
        return None
    
    def process_header(self, header_line: str) -> str:
        """Process a header line to inject species_taxid."""
        # Remove the '>' character for processing
        original_header = header_line.lstrip('>')
        
        # Extract assembly accession
        assembly_accession = self.extract_assembly_accession(header_line)
        
        if not assembly_accession:
            logger.debug(f"Could not extract assembly accession from header: {header_line[:100]}...")
            self.error_count += 1
            return header_line  # Return original header if can't process
        
        # Look up taxid
        if assembly_accession not in self.assembly_summary:
            logger.debug(f"Assembly accession {assembly_accession} not found in summary")
            self.missing_taxid_count += 1
            return header_line  # Return original header if no taxid found
        
        metadata = self.assembly_summary[assembly_accession]
        taxid = metadata['taxid']
        
        if not taxid or taxid == 'na' or taxid == 'NA':
            logger.debug(f"No valid taxid for {assembly_accession}")
            self.missing_taxid_count += 1
            return header_line  # Return original header if no valid taxid
        
        # Check if header already has taxid format
        if '|' in original_header and original_header.startswith(assembly_accession):
            # Header might already be formatted, check if second field is a number
            parts = original_header.split('|')
            if len(parts) >= 2 and parts[1].isdigit():
                # Replace existing taxid with species_taxid
                parts[1] = taxid
                return '>' + '|'.join(parts)
        
        # Create new header: >assembly_accession|species_taxid|original_header
        new_header = f">{assembly_accession}|{taxid}|{original_header}"
        
        self.processed_sequences += 1
        if self.processed_sequences % 1000 == 0:
            logger.info(f"Processed {self.processed_sequences} sequences...")
        
        return new_header
    
    def process_library_file(self, input_file: str, output_file: str):
        """Process a large library.fna file."""
        logger.info(f"Processing library file: {input_file}")
        logger.info(f"Output will be written to: {output_file}")
        
        # Determine if input is compressed
        is_compressed = input_file.endswith('.gz')
        
        try:
            # Open input file
            if is_compressed:
                infile = gzip.open(input_file, 'rt')
            else:
                infile = open(input_file, 'r')
            
            # Process file line by line to handle large files efficiently
            with infile:
                with open(output_file, 'w') as outfile:
                    for line in infile:
                        if line.startswith('>'):
                            # Process header line
                            processed_header = self.process_header(line.strip())
                            outfile.write(processed_header + '\n')
                        else:
                            # Write sequence lines as-is
                            outfile.write(line)
            
            logger.info(f"\nProcessing complete!")
            logger.info(f"Total sequences processed: {self.processed_sequences}")
            logger.info(f"Sequences with errors: {self.error_count}")
            logger.info(f"Sequences with missing taxids: {self.missing_taxid_count}")
            
        except Exception as e:
            logger.error(f"Error processing library file: {e}")
            sys.exit(1)
    
    def create_processing_report(self, output_dir: str) -> str:
        """Create a report of processing statistics."""
        report_file = os.path.join(output_dir, "library_processing_report.txt")
        
        with open(report_file, 'w') as f:
            f.write("Library FNA Processing Report\n")
            f.write("============================\n\n")
            f.write(f"Total sequences processed: {self.processed_sequences}\n")
            f.write(f"Sequences with extraction errors: {self.error_count}\n")
            f.write(f"Sequences with missing taxids: {self.missing_taxid_count}\n")
            f.write(f"\nSuccess rate: {self.processed_sequences / (self.processed_sequences + self.error_count + self.missing_taxid_count) * 100:.2f}%\n")
        
        return report_file

def main():
    parser = argparse.ArgumentParser(description="Process library.fna files to inject species_taxid")
    parser.add_argument("--library_fna", required=True,
                       help="Path to library.fna file (can be .gz compressed)")
    parser.add_argument("--assembly_summary", required=True,
                       help="Path to assembly_summary.txt file")
    parser.add_argument("--output_file", required=True,
                       help="Output file path for processed library.fna")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate inputs
    if not os.path.exists(args.library_fna):
        logger.error(f"Library file does not exist: {args.library_fna}")
        sys.exit(1)
    
    if not os.path.exists(args.assembly_summary):
        logger.error(f"Assembly summary file does not exist: {args.assembly_summary}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Process the library file
    processor = LibraryFNAProcessor(args.assembly_summary)
    processor.process_library_file(args.library_fna, args.output_file)
    
    # Create report
    if output_dir:
        report_file = processor.create_processing_report(output_dir)
        logger.info(f"Processing report written to: {report_file}")

if __name__ == "__main__":
    main()