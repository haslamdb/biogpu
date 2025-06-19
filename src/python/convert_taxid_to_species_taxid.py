#!/usr/bin/env python3
"""
Convert kraken library.fna files from strain taxid to species taxid.

This script:
1. Reads library.fna files with headers like: >kraken:taxid|455631|...
2. Builds a mapping from strain taxid to species taxid using assembly_summary.txt
3. Replaces strain taxids with species taxids in headers
4. Writes the modified sequences to output file

Usage:
    python convert_taxid_to_species_taxid.py --library_fna /path/to/library.fna --assembly_summary /path/to/assembly_summary.txt --output_file library_species.fna
"""

import os
import re
import sys
import gzip
import argparse
import pandas as pd
from typing import Dict, Optional
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TaxidConverter:
    def __init__(self, assembly_summary_path: str):
        """Initialize converter with assembly summary."""
        self.processed_sequences = 0
        self.converted_count = 0
        self.unchanged_count = 0
        self.taxid_to_species_taxid = self.build_taxid_mapping(assembly_summary_path)
    
    def build_taxid_mapping(self, summary_path: str) -> Dict[str, str]:
        """Build mapping from strain taxid to species taxid."""
        logger.info(f"Building taxid to species_taxid mapping from {summary_path}")
        
        try:
            # Read the assembly summary file
            df = pd.read_csv(summary_path, sep='\t', comment='#', 
                           dtype={'taxid': str, 'species_taxid': str})
            
            # Build mapping from taxid to species_taxid
            mapping = {}
            for _, row in df.iterrows():
                # Get taxid and species_taxid
                if 'taxid' in df.columns and 'species_taxid' in df.columns:
                    taxid = str(row['taxid'])
                    species_taxid = str(row['species_taxid'])
                else:
                    # Use column positions if names not found
                    taxid = str(row.iloc[5])  # Usually 6th column
                    species_taxid = str(row.iloc[6])  # Usually 7th column
                
                # Only add to mapping if both are valid
                if taxid and taxid != 'na' and species_taxid and species_taxid != 'na':
                    mapping[taxid] = species_taxid
            
            # Add identity mappings for species that are their own species_taxid
            for species_id in set(mapping.values()):
                if species_id not in mapping:
                    mapping[species_id] = species_id
            
            logger.info(f"Built mapping for {len(mapping)} taxids")
            return mapping
            
        except Exception as e:
            logger.error(f"Error building taxid mapping: {e}")
            sys.exit(1)
    
    def process_kraken_header(self, header_line: str) -> str:
        """Process a kraken-formatted header to convert taxid to species_taxid."""
        # Pattern: >kraken:taxid|TAXID|rest_of_header
        match = re.match(r'^>kraken:taxid\|(\d+)\|(.+)$', header_line)
        
        if not match:
            # Not a kraken header format, return unchanged
            self.unchanged_count += 1
            return header_line
        
        strain_taxid = match.group(1)
        rest_of_header = match.group(2)
        
        # Look up species_taxid
        if strain_taxid in self.taxid_to_species_taxid:
            species_taxid = self.taxid_to_species_taxid[strain_taxid]
            
            # Only change if different
            if strain_taxid != species_taxid:
                self.converted_count += 1
                new_header = f">kraken:taxid|{species_taxid}|{rest_of_header}"
                logger.debug(f"Converted taxid {strain_taxid} to species {species_taxid}")
                return new_header
        else:
            logger.debug(f"No species mapping found for taxid {strain_taxid}")
        
        self.unchanged_count += 1
        return header_line
    
    def process_library_file(self, input_file: str, output_file: str):
        """Process a library.fna file to convert taxids."""
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
            
            # Process file line by line
            with infile:
                with open(output_file, 'w') as outfile:
                    for line in infile:
                        if line.startswith('>'):
                            # Process header line
                            self.processed_sequences += 1
                            processed_header = self.process_kraken_header(line.strip())
                            outfile.write(processed_header + '\n')
                            
                            if self.processed_sequences % 1000 == 0:
                                logger.info(f"Processed {self.processed_sequences} sequences "
                                          f"({self.converted_count} converted)...")
                        else:
                            # Write sequence lines as-is
                            outfile.write(line)
            
            logger.info(f"\nProcessing complete!")
            logger.info(f"Total sequences processed: {self.processed_sequences}")
            logger.info(f"Sequences converted to species taxid: {self.converted_count}")
            logger.info(f"Sequences unchanged: {self.unchanged_count}")
            
        except Exception as e:
            logger.error(f"Error processing library file: {e}")
            sys.exit(1)
    
    def create_conversion_report(self, output_dir: str) -> str:
        """Create a report of conversion statistics."""
        report_file = os.path.join(output_dir, "taxid_conversion_report.txt")
        
        with open(report_file, 'w') as f:
            f.write("Taxid to Species Taxid Conversion Report\n")
            f.write("=======================================\n\n")
            f.write(f"Total sequences processed: {self.processed_sequences}\n")
            f.write(f"Sequences converted: {self.converted_count}\n")
            f.write(f"Sequences unchanged: {self.unchanged_count}\n")
            f.write(f"\nConversion rate: {self.converted_count / self.processed_sequences * 100:.2f}%\n")
            f.write(f"Total unique taxid mappings: {len(self.taxid_to_species_taxid)}\n")
        
        return report_file

def main():
    parser = argparse.ArgumentParser(description="Convert library.fna from strain to species taxid")
    parser.add_argument("--library_fna", required=True,
                       help="Path to library.fna file with kraken headers")
    parser.add_argument("--assembly_summary", required=True,
                       help="Path to assembly_summary.txt file")
    parser.add_argument("--output_file", required=True,
                       help="Output file path for converted library.fna")
    parser.add_argument("--test_lines", type=int,
                       help="Process only first N lines for testing")
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
    
    # For testing, create a temporary sample file
    input_file = args.library_fna
    if args.test_lines:
        logger.info(f"Test mode: Processing first {args.test_lines} lines")
        test_file = "/tmp/library_test_sample.fna"
        with open(args.library_fna, 'r') as infile:
            with open(test_file, 'w') as outfile:
                for i, line in enumerate(infile):
                    if i >= args.test_lines:
                        break
                    outfile.write(line)
        input_file = test_file
    
    # Process the library file
    converter = TaxidConverter(args.assembly_summary)
    converter.process_library_file(input_file, args.output_file)
    
    # Create report
    if output_dir:
        report_file = converter.create_conversion_report(output_dir)
        logger.info(f"Conversion report written to: {report_file}")

if __name__ == "__main__":
    main()