#!/usr/bin/env python3
"""
run_pipeline_comparison.py
Script to compile and run the pipeline comparison test
"""

import os
import sys
import subprocess
import argparse
import json
from datetime import datetime

def build_comparison_test(build_dir="build_comparison"):
    """Build the pipeline comparison test executable"""
    print("Building pipeline comparison test...")
    
    # Create build directory
    os.makedirs(build_dir, exist_ok=True)
    
    # Create CMakeLists.txt for the comparison test
    cmake_content = """
cmake_minimum_required(VERSION 3.20)
project(PipelineComparison LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Find packages
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(ZLIB REQUIRED)

# Include directories from main project
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(../include)
include_directories(../runtime/kernels/resistance)

# CUDA include directories
include_directories(/usr/local/cuda/include)
if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# CUDA settings
set(CMAKE_CUDA_ARCHITECTURES "61;70;75;80;86")
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

# Add the comparison test executable
add_executable(pipeline_comparison_test
    ../pipeline_comparison_test.cpp
    ../runtime/kernels/resistance/fq_mutation_detector.cu
    ../runtime/kernels/resistance/bloom_filter.cu
    ../runtime/kernels/resistance/kmer_screening.cu
    ../runtime/kernels/resistance/translated_search.cu
    ../runtime/kernels/resistance/enhanced_mutation_detection.cu
    ../runtime/kernels/resistance/diagnostic_report.cpp
    ../runtime/kernels/resistance/hdf5_alignment_writer.cpp
)

set_target_properties(pipeline_comparison_test PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Link libraries
target_link_libraries(pipeline_comparison_test
    ${ZLIB_LIBRARIES}
    cudart
    cuda
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)
"""
    
    with open(os.path.join(build_dir, "CMakeLists.txt"), "w") as f:
        f.write(cmake_content)
    
    # Run cmake and make
    try:
        subprocess.run(["cmake", "."], cwd=build_dir, check=True)
        subprocess.run(["make", "-j4", "pipeline_comparison_test"], cwd=build_dir, check=True)
        print("Build successful!")
        return os.path.join(build_dir, "pipeline_comparison_test")
    except subprocess.CalledProcessError as e:
        print(f"Build failed: {e}")
        return None

def run_comparison_test(executable, index_path, r1_path, r2_path, protein_db_path, 
                       batch_size=10000, enable_sw=False, output_file=None):
    """Run the pipeline comparison test"""
    cmd = [
        executable,
        index_path,
        r1_path,
        r2_path,
        protein_db_path,
        str(batch_size)
    ]
    
    if enable_sw:
        cmd.append("--smith-waterman")
    
    print(f"\nRunning: {' '.join(cmd)}")
    print("This may take several minutes depending on file size...\n")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Print output
        print(result.stdout)
        
        # Save to file if requested
        if output_file:
            with open(output_file, "w") as f:
                f.write(f"Pipeline Comparison Test Results\n")
                f.write(f"Generated: {datetime.now()}\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write("=" * 80 + "\n\n")
                f.write(result.stdout)
            print(f"\nResults saved to: {output_file}")
        
        return parse_results(result.stdout)
        
    except subprocess.CalledProcessError as e:
        print(f"Test failed: {e}")
        print(f"Error output: {e.stderr}")
        return None

def parse_results(output):
    """Parse the test output to extract key metrics"""
    results = {}
    
    # Extract timing and performance data
    lines = output.split('\n')
    current_mode = None
    
    for line in lines:
        if "Running" in line and "..." in line:
            # Extract pipeline mode being tested
            mode = line.split("Running ")[1].split("...")[0]
            current_mode = mode
            results[mode] = {}
        
        # Parse the results table
        if any(mode in line for mode in ["Full Pipeline", "Skip Bloom", "Skip K-mer", "Translation Only"]):
            parts = line.split()
            if len(parts) >= 6:
                mode = " ".join(parts[:-5])
                if mode in results:
                    try:
                        results[mode]["total_time_ms"] = float(parts[-5].replace("ms", ""))
                        results[mode]["reads_per_sec"] = int(parts[-4])
                        results[mode]["proteins_found"] = int(parts[-3])
                        results[mode]["mutations_found"] = int(parts[-2])
                        results[mode]["hit_rate"] = float(parts[-1].replace("%", ""))
                    except (ValueError, IndexError):
                        pass
    
    return results

def generate_recommendation(results):
    """Generate recommendations based on results"""
    if not results:
        return "No results to analyze"
    
    recommendations = []
    
    # Check if we have all pipeline modes
    full = results.get("Full Pipeline", {})
    skip_bloom = results.get("Skip Bloom", {})
    skip_kmer = results.get("Skip K-mer", {})
    trans_only = results.get("Translation Only", {})
    
    if full and skip_kmer:
        # Compare protein detection rates
        if "proteins_found" in full and "proteins_found" in skip_kmer:
            protein_ratio = skip_kmer["proteins_found"] / full["proteins_found"] if full["proteins_found"] > 0 else 0
            
            if protein_ratio >= 0.95:
                time_saved = full.get("total_time_ms", 0) - skip_kmer.get("total_time_ms", 0)
                speedup = full.get("total_time_ms", 1) / skip_kmer.get("total_time_ms", 1)
                
                recommendations.append(
                    f"✓ K-mer filtering appears redundant. Skipping it maintains {protein_ratio*100:.1f}% "
                    f"protein detection with {speedup:.2f}x speedup ({time_saved:.0f}ms saved)."
                )
            else:
                recommendations.append(
                    f"✗ K-mer filtering is necessary. Removing it reduces protein detection to {protein_ratio*100:.1f}%."
                )
    
    # Find optimal configuration
    best_mode = None
    best_time = float('inf')
    
    for mode, stats in results.items():
        if mode == "Full Pipeline":
            continue
        
        time = stats.get("total_time_ms", float('inf'))
        proteins = stats.get("proteins_found", 0)
        full_proteins = full.get("proteins_found", 1)
        
        if proteins / full_proteins >= 0.95 and time < best_time:
            best_time = time
            best_mode = mode
    
    if best_mode:
        recommendations.append(f"\n✓ Recommended configuration: {best_mode}")
        recommendations.append(f"  - Maintains >95% protein detection")
        recommendations.append(f"  - {(full.get('total_time_ms', 0) / best_time):.2f}x faster than full pipeline")
    
    return "\n".join(recommendations)

def main():
    parser = argparse.ArgumentParser(description="Compare FQ resistance detection pipeline configurations")
    parser.add_argument("index_path", help="Path to k-mer index")
    parser.add_argument("r1_fastq", help="R1 FASTQ file (can be gzipped)")
    parser.add_argument("r2_fastq", help="R2 FASTQ file (can be gzipped)")
    parser.add_argument("protein_db", help="Path to protein database")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size for processing")
    parser.add_argument("--smith-waterman", action="store_true", help="Enable Smith-Waterman alignment")
    parser.add_argument("--build-dir", default="build_comparison", help="Build directory")
    parser.add_argument("--output", help="Save results to file")
    parser.add_argument("--skip-build", action="store_true", help="Skip building if executable exists")
    
    args = parser.parse_args()
    
    # Check if executable exists or build it
    executable = os.path.join(args.build_dir, "pipeline_comparison_test")
    
    if not os.path.exists(executable) or not args.skip_build:
        executable = build_comparison_test(args.build_dir)
        if not executable:
            print("Failed to build comparison test")
            return 1
    else:
        print(f"Using existing executable: {executable}")
    
    # Run the comparison test
    results = run_comparison_test(
        executable,
        args.index_path,
        args.r1_fastq,
        args.r2_fastq,
        args.protein_db,
        args.batch_size,
        args.smith_waterman,
        args.output
    )
    
    # Generate and print recommendations
    if results:
        print("\n" + "=" * 80)
        print("ANALYSIS AND RECOMMENDATIONS")
        print("=" * 80)
        print(generate_recommendation(results))
        
        # Save results as JSON for further analysis
        json_file = args.output.replace(".txt", ".json") if args.output else "pipeline_comparison_results.json"
        with open(json_file, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nDetailed results saved to: {json_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())