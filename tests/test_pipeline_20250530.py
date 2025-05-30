#!/usr/bin/env python3
"""
Complete Pipeline Test Script

This script orchestrates the complete testing workflow:
1. Build enhanced k-mer index
2. Validate index integrity
3. Test with synthetic data
4. Test with real metagenomic data
5. Debug any issues found

Usage:
    python test_pipeline.py --sequences-dir data/sequences --mutations-csv data/mutations.csv --test-data data/test_fastq
"""

import os
import sys
import subprocess
import argparse
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PipelineTester:
    def __init__(self, sequences_dir, mutations_csv, test_data_dir, output_dir):
        self.sequences_dir = sequences_dir
        self.mutations_csv = mutations_csv
        self.test_data_dir = test_data_dir
        self.output_dir = output_dir
        self.index_dir = os.path.join(output_dir, 'index')
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(self.index_dir, exist_ok=True)
    
    def check_prerequisites(self):
        """Check that all required files and directories exist"""
        logger.info("Checking prerequisites...")
        
        missing = []
        
        if not os.path.exists(self.sequences_dir):
            missing.append(f"Sequences directory: {self.sequences_dir}")
        
        if not os.path.exists(self.mutations_csv):
            missing.append(f"Mutations CSV: {self.mutations_csv}")
        
        if self.test_data_dir and not os.path.exists(self.test_data_dir):
            missing.append(f"Test data directory: {self.test_data_dir}")
        
        # Check for sequence files
        seq_count = 0
        if os.path.exists(self.sequences_dir):
            for org_dir in Path(self.sequences_dir).iterdir():
                if org_dir.is_dir():
                    for gene_file in org_dir.glob('*.json'):
                        seq_count += 1
        
        if seq_count == 0:
            missing.append("No sequence JSON files found in sequences directory")
        
        if missing:
            logger.error("Missing prerequisites:")
            for item in missing:
                logger.error(f"  - {item}")
            return False
        
        logger.info(f"Prerequisites check passed. Found {seq_count} sequence files.")
        return True
    
    def build_index(self):
        """Build the k-mer index using enhanced builder"""
        logger.info("Building enhanced k-mer index...")
        
        # Create the enhanced k-mer builder script if needed
        builder_script = os.path.join(os.path.dirname(__file__), '..', 'src', 'python', 'enhanced_kmer_builder.py')
        
        cmd = [
            'python3', builder_script,
            self.sequences_dir,
            self.mutations_csv,
            self.index_dir,
            '--kmer-length', '15'
        ]
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(__file__))
            
            if result.returncode != 0:
                logger.error("Index building failed!")
                logger.error(f"STDOUT: {result.stdout}")
                logger.error(f"STDERR: {result.stderr}")
                return False
            
            logger.info("Index building completed successfully")
            logger.info(result.stdout)
            return True
            
        except Exception as e:
            logger.error(f"Error running index builder: {e}")
            return False
    
    def validate_index(self):
        """Validate the created index"""
        logger.info("Validating index...")
        
        validator_script = os.path.join(os.path.dirname(__file__), '..', 'src', 'python', 'index_validator.py')
        
        cmd = [
            'python3', validator_script,
            self.index_dir,
            '--create-synthetic'
        ]
        
        if self.test_data_dir:
            cmd.extend(['--test-data', self.test_data_dir])
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(__file__))
            
            if result.returncode != 0:
                logger.error("Index validation failed!")
                logger.error(f"STDOUT: {result.stdout}")
                logger.error(f"STDERR: {result.stderr}")
                return False
            
            logger.info("Index validation completed successfully")
            logger.info(result.stdout)
            return True
            
        except Exception as e:
            logger.error(f"Error running index validator: {e}")
            return False
    
    def test_cuda_screening(self):
        """Test CUDA k-mer screening if binary is available"""
        logger.info("Testing CUDA k-mer screening...")
        
        # Look for the CUDA test binary
        cuda_binary = None
        possible_paths = [
            './build/fq_pipeline_gpu',
            './fq_pipeline_gpu',
            '../build/fq_pipeline_gpu'
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                cuda_binary = path
                break
        
        if not cuda_binary:
            logger.warning("CUDA binary not found, skipping CUDA test")
            return True
        
        # Test with synthetic data first
        synthetic_dir = os.path.join(self.index_dir, 'synthetic_test_data')
        if os.path.exists(synthetic_dir):
            r1_path = os.path.join(synthetic_dir, 'synthetic_reads_R1.fastq.gz')
            r2_path = os.path.join(synthetic_dir, 'synthetic_reads_R2.fastq.gz')
            output_path = os.path.join(self.output_dir, 'synthetic_cuda_results.json')
            
            if os.path.exists(r1_path) and os.path.exists(r2_path):
                # Use binary index for CUDA test
                index_file = os.path.join(self.index_dir, 'kmer_index.bin')
                
                cmd = [cuda_binary, index_file, r1_path, r2_path, output_path]
                
                logger.info(f"Running CUDA test: {' '.join(cmd)}")
                
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                    
                    if result.returncode != 0:
                        logger.error("CUDA test failed!")
                        logger.error(f"STDOUT: {result.stdout}")
                        logger.error(f"STDERR: {result.stderr}")
                        return False
                    
                    logger.info("CUDA test completed successfully")
                    logger.info(result.stdout)
                    
                    # Analyze results
                    if os.path.exists(output_path):
                        self.analyze_cuda_results(output_path, "synthetic")
                    
                    return True
                    
                except subprocess.TimeoutExpired:
                    logger.error("CUDA test timed out after 5 minutes")
                    return False
                except Exception as e:
                    logger.error(f"Error running CUDA test: {e}")
                    return False
        
        logger.warning("Synthetic test data not found, skipping CUDA test")
        return True
    
    def test_real_data(self):
        """Test with real metagenomic data"""
        if not self.test_data_dir:
            logger.info("No test data directory specified, skipping real data test")
            return True
        
        logger.info("Testing with real metagenomic data...")
        
        # Look for VRE12 files
        r1_path = os.path.join(self.test_data_dir, 'VRE12_R1.fastq.gz')
        r2_path = os.path.join(self.test_data_dir, 'VRE12_R2.fastq.gz')
        
        if not (os.path.exists(r1_path) and os.path.exists(r2_path)):
            logger.warning("VRE12 test files not found, skipping real data test")
            return True
        
        # Test with CUDA binary if available
        cuda_binary = None
        possible_paths = [
            './build/fq_pipeline_gpu',
            './fq_pipeline_gpu',
            '../build/fq_pipeline_gpu'
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                cuda_binary = path
                break
        
        if cuda_binary:
            index_file = os.path.join(self.index_dir, 'kmer_index.bin')
            output_path = os.path.join(self.output_dir, 'real_data_cuda_results.json')
            
            cmd = [cuda_binary, index_file, r1_path, r2_path, output_path]
            
            logger.info(f"Running real data test: {' '.join(cmd)}")
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                
                if result.returncode != 0:
                    logger.error("Real data test failed!")
                    logger.error(f"STDOUT: {result.stdout}")
                    logger.error(f"STDERR: {result.stderr}")
                    return False
                
                logger.info("Real data test completed successfully")
                logger.info(result.stdout)
                
                # Analyze results
                if os.path.exists(output_path):
                    self.analyze_cuda_results(output_path, "real")
                
                return True
                
            except subprocess.TimeoutExpired:
                logger.error("Real data test timed out after 10 minutes")
                return False
            except Exception as e:
                logger.error(f"Error running real data test: {e}")
                return False
        
        logger.warning("CUDA binary not available for real data test")
        return True
    
    def analyze_cuda_results(self, results_path, data_type):
        """Analyze CUDA screening results"""
        logger.info(f"Analyzing {data_type} data results...")
        
        try:
            with open(results_path, 'r') as f:
                results = json.load(f)
            
            summary = results.get('summary', {})
            mutations = results.get('mutations', [])
            
            logger.info(f"Results for {data_type} data:")
            logger.info(f"  Total reads processed: {summary.get('total_reads', 'unknown')}")
            logger.info(f"  Mutations found: {len(mutations)}")
            
            if mutations:
                logger.info("Sample mutations found:")
                for i, mut in enumerate(mutations[:5]):  # Show first 5
                    logger.info(f"    Mutation {i+1}: read_pair={mut.get('read_pair')}, "
                               f"gene_id={mut.get('gene_id')}, species_id={mut.get('species_id')}")
            else:
                if data_type == "synthetic":
                    logger.warning("WARNING: No mutations found in synthetic data - this may indicate a problem")
                else:
                    logger.info("No mutations found in real data (this may be expected)")
            
        except Exception as e:
            logger.error(f"Error analyzing results: {e}")
    
    def generate_report(self):
        """Generate a comprehensive test report"""
        logger.info("Generating test report...")
        
        report_path = os.path.join(self.output_dir, 'test_report.md')
        
        with open(report_path, 'w') as f:
            f.write("# BioGPU Pipeline Test Report\n\n")
            f.write(f"Generated on: {pd.Timestamp.now().isoformat()}\n\n")
            
            f.write("## Test Configuration\n")
            f.write(f"- Sequences directory: {self.sequences_dir}\n")
            f.write(f"- Mutations CSV: {self.mutations_csv}\n")
            f.write(f"- Test data directory: {self.test_data_dir}\n")
            f.write(f"- Output directory: {self.output_dir}\n\n")
            
            f.write("## Index Files Created\n")
            index_files = [
                'kmer_index.bin',
                'sequences.bin', 
                'index_metadata.json',
                'debug/kmer_frequency_distribution.txt',
                'debug/organism_gene_summary.txt',
                'debug/validation_summary.txt'
            ]
            
            for file_name in index_files:
                file_path = os.path.join(self.index_dir, file_name)
                status = "✅ Created" if os.path.exists(file_path) else "❌ Missing"
                f.write(f"- {file_name}: {status}\n")
            
            f.write("\n## Test Results\n")
            
            # Check for various result files
            result_files = [
                ('synthetic_cuda_results.json', 'Synthetic data CUDA test'),
                ('real_data_cuda_results.json', 'Real data CUDA test')
            ]
            
            for file_name, description in result_files:
                file_path = os.path.join(self.output_dir, file_name)
                if os.path.exists(file_path):
                    f.write(f"### {description}\n")
                    try:
                        with open(file_path, 'r') as rf:
                            results = json.load(rf)
                        summary = results.get('summary', {})
                        mutations = results.get('mutations', [])
                        f.write(f"- Total reads: {summary.get('total_reads', 'unknown')}\n")
                        f.write(f"- Mutations found: {len(mutations)}\n")
                    except:
                        f.write("- Error reading results file\n")
                    f.write("\n")
            
            f.write("## Debug Information\n")
            f.write("Check the following files for detailed debugging information:\n")
            f.write(f"- Index debug reports: {os.path.join(self.index_dir, 'debug')}\n")
            f.write(f"- Synthetic test data: {os.path.join(self.index_dir, 'synthetic_test_data')}\n")
        
        logger.info(f"Test report generated: {report_path}")
    
    def run_complete_test(self):
        """Run the complete test pipeline"""
        logger.info("=== Starting Complete Pipeline Test ===")
        
        # Check prerequisites
        if not self.check_prerequisites():
            logger.error("Prerequisites check failed, aborting")
            return False
        
        # Step 1: Build index
        if not self.build_index():
            logger.error("Index building failed, aborting")
            return False
        
        # Step 2: Validate index
        if not self.validate_index():
            logger.error("Index validation failed, aborting")
            return False
        
        # Step 3: Test CUDA screening
        if not self.test_cuda_screening():
            logger.error("CUDA screening test failed")
            # Don't abort here - continue with other tests
        
        # Step 4: Test with real data
        if not self.test_real_data():
            logger.error("Real data test failed")
            # Don't abort here - this is less critical
        
        # Step 5: Generate report
        self.generate_report()
        
        logger.info("=== Complete Pipeline Test Finished ===")
        logger.info(f"Results available in: {self.output_dir}")
        
        return True

def main():
    parser = argparse.ArgumentParser(description='Complete Pipeline Test Script')
    parser.add_argument('--sequences-dir', required=True,
                       help='Directory containing downloaded sequences')
    parser.add_argument('--mutations-csv', required=True,
                       help='Path to mutations CSV file')
    parser.add_argument('--test-data', 
                       help='Directory containing test FASTQ files (optional)')
    parser.add_argument('--output-dir', default='pipeline_test_output',
                       help='Output directory for test results')
    
    args = parser.parse_args()
    
    # Import pandas for timestamps
    try:
        import pandas as pd
    except ImportError:
        logger.error("pandas is required for this script")
        return 1
    
    tester = PipelineTester(
        args.sequences_dir,
        args.mutations_csv, 
        args.test_data,
        args.output_dir
    )
    
    success = tester.run_complete_test()
    
    if success:
        print(f"\n✅ Pipeline test completed successfully!")
        print(f"Check results in: {args.output_dir}")
        return 0
    else:
        print(f"\n❌ Pipeline test failed!")
        print(f"Check logs and results in: {args.output_dir}")
        return 1

if __name__ == "__main__":
    exit(main())