#!/usr/bin/env python3
"""
Generate a CSV file for batch processing from a directory of FASTQ files.
Automatically pairs R1 and R2 files based on naming patterns.
"""

import os
import argparse
import re
import csv
from pathlib import Path
from collections import defaultdict

def find_fastq_files(directory, extensions=['.fastq.gz', '.fq.gz', '.fastq', '.fq']):
    """Find all FASTQ files in a directory."""
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if any(file.endswith(ext) for ext in extensions):
                fastq_files.append(os.path.join(root, file))
    return fastq_files

def extract_sample_name(filename, patterns):
    """Extract sample name from filename using regex patterns."""
    basename = os.path.basename(filename)
    
    for pattern in patterns:
        match = re.search(pattern, basename)
        if match:
            return match.group(1)
    
    # Fallback: remove common suffixes
    name = basename
    for suffix in ['_R1', '_R2', '_1', '_2', '.fastq.gz', '.fq.gz', '.fastq', '.fq']:
        name = name.replace(suffix, '')
    return name

def pair_reads(fastq_files, r1_patterns, r2_patterns):
    """Pair R1 and R2 files based on sample names."""
    r1_files = {}
    r2_files = {}
    
    # Classify files as R1 or R2
    for file in fastq_files:
        basename = os.path.basename(file)
        
        # Check if it's R1
        is_r1 = any(re.search(pattern, basename) for pattern in r1_patterns)
        # Check if it's R2
        is_r2 = any(re.search(pattern, basename) for pattern in r2_patterns)
        
        if is_r1 and not is_r2:
            sample_name = extract_sample_name(file, r1_patterns)
            r1_files[sample_name] = file
        elif is_r2 and not is_r1:
            sample_name = extract_sample_name(file, r2_patterns)
            r2_files[sample_name] = file
        else:
            print(f"Warning: Could not classify file as R1 or R2: {basename}")
    
    # Pair files
    paired_samples = []
    unpaired_r1 = []
    unpaired_r2 = []
    
    for sample_name, r1_path in r1_files.items():
        if sample_name in r2_files:
            paired_samples.append({
                'name': sample_name,
                'r1': r1_path,
                'r2': r2_files[sample_name]
            })
        else:
            unpaired_r1.append((sample_name, r1_path))
    
    for sample_name, r2_path in r2_files.items():
        if sample_name not in r1_files:
            unpaired_r2.append((sample_name, r2_path))
    
    return paired_samples, unpaired_r1, unpaired_r2

def write_csv(samples, output_file, use_absolute_paths=True):
    """Write samples to CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['SampleName', 'FilePath', 'R1 file', 'R2 file'])
        
        for sample in samples:
            # Extract directory and filenames
            r1_dir = os.path.dirname(sample['r1'])
            r1_filename = os.path.basename(sample['r1'])
            
            if sample['r2']:
                r2_dir = os.path.dirname(sample['r2'])
                r2_filename = os.path.basename(sample['r2'])
                
                # Verify both files are in the same directory
                if r1_dir != r2_dir:
                    print(f"Warning: R1 and R2 for sample {sample['name']} are in different directories!")
                    print(f"  R1: {r1_dir}")
                    print(f"  R2: {r2_dir}")
                    print(f"  Using R1 directory for both files")
            else:
                r2_filename = ''
            
            # Determine the file path
            file_path = r1_dir
            if use_absolute_paths:
                file_path = os.path.abspath(file_path)
            else:
                # Convert to relative paths with ~
                home = str(Path.home())
                if file_path.startswith(home):
                    file_path = '~' + file_path[len(home):]
            
            # Ensure path ends with separator
            if not file_path.endswith('/'):
                file_path += '/'
            
            writer.writerow([sample['name'], file_path, r1_filename, r2_filename])

def main():
    parser = argparse.ArgumentParser(
        description='Generate CSV file for batch processing from FASTQ directory'
    )
    parser.add_argument('directory', help='Directory containing FASTQ files')
    parser.add_argument('-o', '--output', default='samples.csv',
                        help='Output CSV file (default: samples.csv)')
    parser.add_argument('--r1-pattern', nargs='+', 
                        default=[r'(.+?)_R1[._]', r'(.+?)_1[._]', r'(.+?)\.R1\.'],
                        help='Regex patterns to identify R1 files and extract sample names')
    parser.add_argument('--r2-pattern', nargs='+',
                        default=[r'(.+?)_R2[._]', r'(.+?)_2[._]', r'(.+?)\.R2\.'],
                        help='Regex patterns to identify R2 files and extract sample names')
    parser.add_argument('--relative-paths', action='store_true',
                        help='Use relative paths with ~ for home directory')
    parser.add_argument('--include-single-end', action='store_true',
                        help='Include unpaired R1 files as single-end samples')
    parser.add_argument('--recursive', action='store_true',
                        help='Search recursively in subdirectories')
    
    args = parser.parse_args()
    
    # Find FASTQ files
    print(f"Searching for FASTQ files in: {args.directory}")
    if args.recursive:
        print("  (recursive search enabled)")
    
    fastq_files = find_fastq_files(args.directory)
    print(f"Found {len(fastq_files)} FASTQ files")
    
    if not fastq_files:
        print("No FASTQ files found!")
        return
    
    # Pair reads
    paired_samples, unpaired_r1, unpaired_r2 = pair_reads(
        fastq_files, args.r1_pattern, args.r2_pattern
    )
    
    print(f"\nPairing results:")
    print(f"  Paired samples: {len(paired_samples)}")
    print(f"  Unpaired R1 files: {len(unpaired_r1)}")
    print(f"  Unpaired R2 files: {len(unpaired_r2)}")
    
    # Show unpaired files
    if unpaired_r1:
        print("\nUnpaired R1 files:")
        for name, path in unpaired_r1[:5]:  # Show first 5
            print(f"  {name}: {os.path.basename(path)}")
        if len(unpaired_r1) > 5:
            print(f"  ... and {len(unpaired_r1) - 5} more")
    
    if unpaired_r2:
        print("\nUnpaired R2 files:")
        for name, path in unpaired_r2[:5]:  # Show first 5
            print(f"  {name}: {os.path.basename(path)}")
        if len(unpaired_r2) > 5:
            print(f"  ... and {len(unpaired_r2) - 5} more")
    
    # Add single-end samples if requested
    samples_to_write = paired_samples.copy()
    if args.include_single_end and unpaired_r1:
        print(f"\nAdding {len(unpaired_r1)} single-end samples")
        for name, r1_path in unpaired_r1:
            samples_to_write.append({
                'name': name,
                'r1': r1_path,
                'r2': ''  # Empty for single-end
            })
    
    # Write CSV
    if samples_to_write:
        write_csv(samples_to_write, args.output, 
                 use_absolute_paths=not args.relative_paths)
        print(f"\nWrote {len(samples_to_write)} samples to: {args.output}")
        
        # Show first few samples
        print("\nFirst few samples:")
        for i, sample in enumerate(samples_to_write[:3]):
            print(f"  {sample['name']}:")
            print(f"    R1: {os.path.basename(sample['r1'])}")
            if sample['r2']:
                print(f"    R2: {os.path.basename(sample['r2'])}")
    else:
        print("\nNo samples to write!")

if __name__ == '__main__':
    main()