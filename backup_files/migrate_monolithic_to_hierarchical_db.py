#!/usr/bin/env python3
"""
Migration script for BioGPU Hierarchical Database
Helps transition from monolithic to hierarchical database for memory-constrained systems
"""

import os
import sys
import subprocess
import argparse
import json
import time
from pathlib import Path

def get_gpu_memory():
    """Get available GPU memory in GB"""
    try:
        result = subprocess.run(['nvidia-smi', '--query-gpu=memory.total', '--format=csv,noheader,nounits'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            memory_mb = int(result.stdout.strip())
            return memory_mb / 1024  # Convert to GB
    except:
        pass
    return None

def estimate_database_size(kmer_file):
    """Estimate the size of the database from k-mer file"""
    try:
        # Count lines in k-mer file
        result = subprocess.run(['wc', '-l', kmer_file], capture_output=True, text=True)
        if result.returncode == 0:
            lines = int(result.stdout.split()[0])
            # Subtract header lines (rough estimate)
            kmers = max(0, lines - 10)
            # Each k-mer entry is ~16 bytes
            size_bytes = kmers * 16
            return size_bytes / (1024**3)  # Convert to GB
    except:
        pass
    return None

def recommend_strategy(gpu_memory_gb, db_size_gb):
    """Recommend the best strategy based on available memory and database size"""
    if gpu_memory_gb is None:
        return "hierarchical", 8, "Cannot detect GPU memory - using hierarchical as safe default"
    
    if db_size_gb is None:
        return "hierarchical", min(8, int(gpu_memory_gb * 0.8)), "Cannot estimate database size - using hierarchical"
    
    # If database fits comfortably in GPU memory (with 50% headroom)
    if db_size_gb * 1.5 < gpu_memory_gb:
        return "monolithic", int(gpu_memory_gb * 0.9), f"Database ({db_size_gb:.1f}GB) fits in GPU memory ({gpu_memory_gb:.1f}GB)"
    
    # If database is much larger than GPU memory
    elif db_size_gb > gpu_memory_gb * 3:
        tier_size = min(512, int(gpu_memory_gb * 0.8 * 1024 / 2))  # Half of available memory per tier
        return "hierarchical", tier_size, f"Large database ({db_size_gb:.1f}GB) vs GPU memory ({gpu_memory_gb:.1f}GB) - using streaming"
    
    # In between - hierarchical is safer
    else:
        tier_size = min(1024, int(gpu_memory_gb * 0.6 * 1024))  # 60% of GPU memory per tier
        return "hierarchical", tier_size, f"Moderate database ({db_size_gb:.1f}GB) - hierarchical for memory safety"

def build_hierarchical_database(kmer_file, output_prefix, tier_size_mb):
    """Build hierarchical database"""
    print(f"Building hierarchical database...")
    print(f"  Input: {kmer_file}")
    print(f"  Output: {output_prefix}")
    print(f"  Tier size: {tier_size_mb} MB")
    
    cmd = [
        './build_hierarchical_db',
        kmer_file,
        output_prefix,
        '--tier-size', str(tier_size_mb)
    ]
    
    print(f"Running: {' '.join(cmd)}")
    start_time = time.time()
    
    result = subprocess.run(cmd)
    
    build_time = time.time() - start_time
    print(f"Database build completed in {build_time:.1f} seconds")
    
    return result.returncode == 0

def build_monolithic_database(kmer_file, output_file):
    """Build traditional monolithic database"""
    print(f"Building monolithic database...")
    print(f"  Input: {kmer_file}")
    print(f"  Output: {output_file}")
    
    cmd = ['./build_db_from_kmers', kmer_file, output_file]
    
    print(f"Running: {' '.join(cmd)}")
    start_time = time.time()
    
    result = subprocess.run(cmd)
    
    build_time = time.time() - start_time
    print(f"Database build completed in {build_time:.1f} seconds")
    
    return result.returncode == 0

def create_subset_database(kmer_file, output_file, target_size):
    """Create a stratified subset for testing"""
    print(f"Creating stratified subset...")
    print(f"  Input: {kmer_file}")
    print(f"  Output: {output_file}")
    print(f"  Target size: {target_size} k-mers")
    
    # First create the subset tool if it doesn't exist
    subset_script = """#!/usr/bin/env python3
import random
import sys

def create_stratified_subset(input_file, output_file, target_size):
    # Read all k-mers
    kmers = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                kmers.append(line.strip())
    
    print(f"Loaded {len(kmers)} k-mers")
    
    if len(kmers) <= target_size:
        print(f"Input size ({len(kmers)}) <= target size ({target_size}), copying all")
        with open(output_file, 'w') as f:
            f.write("# Stratified subset (all k-mers included)\\n")
            for kmer in kmers:
                f.write(kmer + "\\n")
        return
    
    # Stratified sampling - ensure we get representatives from different hash ranges
    kmers.sort()  # Sort by k-mer (roughly by hash)
    
    # Take every Nth k-mer to get uniform distribution
    step = len(kmers) // target_size
    subset = []
    for i in range(0, len(kmers), step):
        subset.append(kmers[i])
        if len(subset) >= target_size:
            break
    
    print(f"Created subset with {len(subset)} k-mers")
    
    with open(output_file, 'w') as f:
        f.write(f"# Stratified subset of {input_file}\\n")
        f.write(f"# Target size: {target_size}, actual size: {len(subset)}\\n")
        for kmer in subset:
            f.write(kmer + "\\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: create_subset.py <input_kmers> <output_kmers> <target_size>")
        sys.exit(1)
    
    create_stratified_subset(sys.argv[1], sys.argv[2], int(sys.argv[3]))
"""
    
    # Write temporary script
    script_path = "temp_create_subset.py"
    with open(script_path, 'w') as f:
        f.write(subset_script)
    
    try:
        cmd = ['python3', script_path, kmer_file, output_file, str(target_size)]
        result = subprocess.run(cmd)
        return result.returncode == 0
    finally:
        if os.path.exists(script_path):
            os.remove(script_path)

def test_profiler(db_path, fastq_file, strategy, memory_limit, is_hierarchical):
    """Test the profiler with the built database"""
    print(f"\nTesting profiler...")
    print(f"  Database: {db_path}")
    print(f"  FASTQ: {fastq_file}")
    print(f"  Strategy: {strategy}")
    
    if is_hierarchical:
        cmd = [
            './hierarchical_profiler_pipeline',
            db_path,
            fastq_file,
            '--memory', str(memory_limit),
            '--output-prefix', 'migration_test'
        ]
    else:
        cmd = [
            './gpu_profiler_pipeline',
            db_path,
            fastq_file,
            '--output-prefix', 'migration_test'
        ]
    
    print(f"Running: {' '.join(cmd)}")
    start_time = time.time()
    
    result = subprocess.run(cmd)
    
    test_time = time.time() - start_time
    print(f"Profiler test completed in {test_time:.1f} seconds")
    
    return result.returncode == 0

def main():
    parser = argparse.ArgumentParser(description="BioGPU Hierarchical Database Migration Tool")
    parser.add_argument("kmer_file", help="Input k-mer file")
    parser.add_argument("output_prefix", help="Output database prefix/path")
    parser.add_argument("--fastq", help="Test FASTQ file (optional)")
    parser.add_argument("--force-strategy", choices=["monolithic", "hierarchical", "subset"], 
                       help="Force a specific strategy instead of auto-detection")
    parser.add_argument("--memory", type=int, help="Override GPU memory limit in GB")
    parser.add_argument("--tier-size", type=int, help="Override tier size in MB")
    parser.add_argument("--subset-size", type=int, default=5000000, 
                       help="Size of subset database (default: 5M k-mers)")
    parser.add_argument("--dry-run", action="store_true", help="Show recommendations without building")
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.kmer_file):
        print(f"Error: K-mer file not found: {args.kmer_file}")
        return 1
    
    if args.fastq and not os.path.exists(args.fastq):
        print(f"Error: FASTQ file not found: {args.fastq}")
        return 1
    
    print("BioGPU Hierarchical Database Migration Tool")
    print("=" * 50)
    
    # Analyze system and database
    gpu_memory = args.memory or get_gpu_memory()
    db_size = estimate_database_size(args.kmer_file)
    
    print(f"\nSystem Analysis:")
    print(f"  GPU Memory: {gpu_memory:.1f} GB" if gpu_memory else "  GPU Memory: Not detected")
    print(f"  Database Size: {db_size:.1f} GB" if db_size else "  Database Size: Cannot estimate")
    
    # Get recommendation
    if args.force_strategy:
        strategy = args.force_strategy
        if strategy == "hierarchical":
            param = args.tier_size or (int(gpu_memory * 0.6 * 1024) if gpu_memory else 512)
        elif strategy == "subset":
            param = args.subset_size
        else:
            param = int(gpu_memory * 0.9) if gpu_memory else 8
        reason = f"User forced strategy: {strategy}"
    else:
        strategy, param, reason = recommend_strategy(gpu_memory, db_size)
    
    print(f"\nRecommendation:")
    print(f"  Strategy: {strategy}")
    if strategy == "hierarchical":
        print(f"  Tier size: {param} MB")
    elif strategy == "subset":
        print(f"  Subset size: {param} k-mers")
    else:
        print(f"  Memory limit: {param} GB")
    print(f"  Reason: {reason}")
    
    if args.dry_run:
        print("\nDry run mode - not building database")
        return 0
    
    print(f"\nBuilding database...")
    
    # Build database based on strategy
    success = False
    
    if strategy == "subset":
        # Create subset first
        subset_file = args.output_prefix + "_subset_kmers.txt"
        if create_subset_database(args.kmer_file, subset_file, param):
            # Then build hierarchical database from subset
            success = build_hierarchical_database(subset_file, args.output_prefix, 256)
        is_hierarchical = True
        memory_limit = min(4, gpu_memory or 4)
        
    elif strategy == "hierarchical":
        success = build_hierarchical_database(args.kmer_file, args.output_prefix, param)
        is_hierarchical = True
        memory_limit = param // 1024 + 2  # Convert MB to GB with some headroom
        
    else:  # monolithic
        output_file = args.output_prefix + ".bin"
        success = build_monolithic_database(args.kmer_file, output_file)
        is_hierarchical = False
        memory_limit = param
        args.output_prefix = output_file  # For testing
    
    if not success:
        print("Database build failed!")
        return 1
    
    print("Database build successful!")
    
    # Test with FASTQ if provided
    if args.fastq:
        success = test_profiler(args.output_prefix, args.fastq, strategy, memory_limit, is_hierarchical)
        if success:
            print("\nProfiler test successful!")
            
            # Print usage recommendations
            print(f"\nUsage for your setup:")
            if is_hierarchical:
                print(f"  ./hierarchical_profiler_pipeline {args.output_prefix} your_reads.fastq --memory {memory_limit}")
            else:
                print(f"  ./gpu_profiler_pipeline {args.output_prefix} your_reads.fastq")
        else:
            print("\nProfiler test failed!")
            return 1
    else:
        print(f"\nDatabase ready for use:")
        if is_hierarchical:
            print(f"  ./hierarchical_profiler_pipeline {args.output_prefix} your_reads.fastq --memory {memory_limit}")
        else:
            print(f"  ./gpu_profiler_pipeline {args.output_prefix} your_reads.fastq")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())