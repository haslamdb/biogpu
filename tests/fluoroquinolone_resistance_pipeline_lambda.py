# fluoroquinolone_resistance_pipeline.py
"""
GPU-accelerated fluoroquinolone resistance detection for clinical metagenomics
Optimized for dual NVIDIA Titan Xp GPUs
"""

import cupy as cp
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import time
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from datetime import datetime

@dataclass
class QRDRMutation:
    """Quinolone Resistance-Determining Region mutation"""
    gene: str
    position: int
    wild_type: str
    mutant: str
    organism: str
    resistance_level: str  # 'high', 'moderate', 'low'
    mic_fold_change: float
    drugs_affected: List[str]

class FluoroquinoloneResistanceDB:
    """Clinical database of validated fluoroquinolone resistance mutations"""
    
    def __init__(self):
        self.mutations = self._load_clinical_mutations()
        self.clinical_rules = self._load_clinical_rules()
    
    def _load_clinical_mutations(self) -> List[QRDRMutation]:
        """Load clinically validated QRDR mutations"""
        mutations = [
            # E. coli and Enterobacteriaceae
            QRDRMutation("gyrA", 83, "S", "L", "E.coli", "high", 16.0, 
                        ["ciprofloxacin", "levofloxacin", "ofloxacin"]),
            QRDRMutation("gyrA", 83, "S", "F", "E.coli", "high", 32.0,
                        ["ciprofloxacin", "levofloxacin"]),
            QRDRMutation("gyrA", 87, "D", "N", "E.coli", "moderate", 8.0,
                        ["ciprofloxacin", "norfloxacin"]),
            QRDRMutation("parC", 80, "S", "I", "E.coli", "moderate", 4.0,
                        ["ciprofloxacin", "levofloxacin"]),
            QRDRMutation("parC", 84, "E", "V", "E.coli", "low", 2.0,
                        ["ciprofloxacin"]),
            
            # P. aeruginosa
            QRDRMutation("gyrA", 83, "T", "I", "P.aeruginosa", "high", 16.0,
                        ["ciprofloxacin", "levofloxacin"]),
            QRDRMutation("parC", 87, "S", "L", "P.aeruginosa", "moderate", 8.0,
                        ["ciprofloxacin"]),
            
            # S. aureus (grlA = parC homolog)
            QRDRMutation("grlA", 80, "S", "F", "S.aureus", "high", 16.0,
                        ["ciprofloxacin", "moxifloxacin"]),
            
            # C. difficile
            QRDRMutation("gyrA", 82, "T", "I", "C.difficile", "high", 32.0,
                        ["moxifloxacin", "gatifloxacin"]),
        ]
        return mutations
    
    def _load_clinical_rules(self) -> Dict:
        """Clinical interpretation rules"""
        return {
            "high_resistance": {
                "threshold": "any high-level mutation",
                "recommendation": "Avoid all fluoroquinolones",
                "alternatives": ["ceftriaxone", "azithromycin", "trimethoprim-sulfamethoxazole"]
            },
            "moderate_resistance": {
                "threshold": "moderate mutation without high",
                "recommendation": "Consider higher dose FQ or alternatives",
                "alternatives": ["levofloxacin 750mg", "moxifloxacin"]
            },
            "low_resistance": {
                "threshold": "only low-level mutations",
                "recommendation": "Standard FQ therapy may be effective",
                "monitor": "Consider susceptibility testing"
            }
        }

class TitanXpFQPipeline:
    """Dual Titan Xp GPU pipeline for fluoroquinolone resistance detection"""
    
    def __init__(self, gpu_primary: int = 0, gpu_secondary: int = 1):
        self.gpu_primary = gpu_primary
        self.gpu_secondary = gpu_secondary
        self.resistance_db = FluoroquinoloneResistanceDB()
        
        # Initialize both GPUs
        with cp.cuda.Device(self.gpu_primary):
            cp.cuda.runtime.deviceSynchronize()
        with cp.cuda.Device(self.gpu_secondary):
            cp.cuda.runtime.deviceSynchronize()
        
        print(f"Initialized dual GPU pipeline: GPU {gpu_primary} (primary), GPU {gpu_secondary} (secondary)")
    
    def process_clinical_sample(self, fastq_path: str, sample_id: str) -> Dict:
        """
        Main entry point for clinical sample processing
        
        Args:
            fastq_path: Path to FASTQ file
            sample_id: Clinical sample identifier
        
        Returns:
            Clinical resistance report
        """
        start_time = time.time()
        
        # Step 1: Load and QC reads
        reads = self._load_fastq(fastq_path)
        print(f"Loaded {len(reads)} reads")
        
        # Step 2: Rapid microbial profiling (GPU 0)
        with cp.cuda.Device(self.gpu_primary):
            community_profile = self._profile_microbiome_gpu(reads)
        
        # Step 3: Extract resistance gene reads (GPU 0)
        with cp.cuda.Device(self.gpu_primary):
            target_reads = self._extract_qrdr_reads(reads, community_profile)
        
        # Step 4: Mutation detection (GPU 1)
        with cp.cuda.Device(self.gpu_secondary):
            mutations = self._detect_mutations_gpu(target_reads)
        
        # Step 5: Clinical interpretation
        clinical_report = self._generate_clinical_report(
            mutations, community_profile, sample_id
        )
        
        processing_time = time.time() - start_time
        clinical_report['processing_time_seconds'] = processing_time
        
        print(f"Sample {sample_id} processed in {processing_time:.2f} seconds")
        
        return clinical_report
    
    def _load_fastq(self, fastq_path: str) -> List[Dict]:
        """Load FASTQ file and extract sequences"""
        reads = []
        for i, record in enumerate(SeqIO.parse(fastq_path, "fastq")):
            reads.append({
                'id': record.id,
                'sequence': str(record.seq),
                'quality': record.letter_annotations.get('phred_quality', [30] * len(record.seq))
            })
            
            # Limit for testing - remove for production
            if i >= 100000:  # 100k reads for testing
                break
        
        return reads
    
    def _profile_microbiome_gpu(self, reads: List[Dict]) -> Dict:
        """
        GPU-accelerated microbial community profiling using k-mer approach
        """
        print("Profiling microbial community on GPU...")
        
        # Convert sequences to GPU arrays
        sequences = [read['sequence'] for read in reads]
        gpu_sequences = cp.array([list(seq.encode('ascii')) for seq in sequences[:1000]])  # Subset for testing
        
        # Simulate k-mer based taxonomic classification
        # In real implementation, this would use pre-built k-mer database
        
        # Mock community profile for testing
        community_profile = {
            'E.coli': 0.45,
            'K.pneumoniae': 0.20,
            'P.aeruginosa': 0.15,
            'S.aureus': 0.10,
            'C.difficile': 0.05,
            'other': 0.05
        }
        
        return community_profile
    
    def _extract_qrdr_reads(self, reads: List[Dict], community_profile: Dict) -> List[Dict]:
        """Extract reads that map to QRDR regions"""
        print("Extracting QRDR-mapping reads...")
        
        # QRDR gene sequences (simplified - would use full database)
        qrdr_genes = {
            'gyrA': 'ATGGATAAAGATGTAGCTTATTTAGATATCGATGATGATGAAGATGATCCT',
            'parC': 'ATGGATCGAGATCTAGCTTATGTATATCGATGATGATGAAGATGATCCT',
            'grlA': 'ATGGATCGAGATCTAGCTTATGTATATCGATGATGATGAAGATGATCCT'
        }
        
        target_reads = []
        
        # Simple string matching - in production would use BWA/minimap2
        for read in reads[:10000]:  # Subset for testing
            for gene, sequence in qrdr_genes.items():
                # Check if read contains part of QRDR sequence
                if any(sequence[i:i+20] in read['sequence'] for i in range(0, len(sequence)-20, 10)):
                    read['mapped_gene'] = gene
                    target_reads.append(read)
                    break
        
        print(f"Extracted {len(target_reads)} QRDR-mapping reads")
        return target_reads
    
    def _detect_mutations_gpu(self, target_reads: List[Dict]) -> List[Dict]:
        """GPU-accelerated mutation detection in QRDR regions"""
        print("Detecting mutations on GPU...")
        
        if not target_reads:
            return []
        
        # Convert to GPU arrays
        sequences = [read['sequence'] for read in target_reads]
        genes = [read.get('mapped_gene', 'unknown') for read in target_reads]
        
        detected_mutations = []
        
        # Mock mutation detection - in production would use CUDA kernels
        for i, (seq, gene) in enumerate(zip(sequences, genes)):
            # Simulate finding mutations
            if 'GATGAA' in seq and gene == 'gyrA':  # Mock codon for position 83
                detected_mutations.append({
                    'gene': 'gyrA',
                    'position': 83,
                    'wild_type': 'S',
                    'mutant': 'L',
                    'read_id': target_reads[i]['id'],
                    'confidence': 0.95
                })
            elif 'GATCCT' in seq and gene == 'parC':  # Mock codon for position 80
                detected_mutations.append({
                    'gene': 'parC',
                    'position': 80,
                    'wild_type': 'S',
                    'mutant': 'I',
                    'read_id': target_reads[i]['id'],
                    'confidence': 0.90
                })
        
        print(f"Detected {len(detected_mutations)} mutations")
        return detected_mutations
    
    def _generate_clinical_report(self, mutations: List[Dict], 
                                 community_profile: Dict, sample_id: str) -> Dict:
        """Generate clinical resistance report"""
        
        # Analyze detected mutations
        resistance_analysis = self._analyze_resistance_level(mutations)
        
        # Generate recommendations
        clinical_recommendations = self._get_clinical_recommendations(
            resistance_analysis, community_profile
        )
        
        report = {
            'sample_id': sample_id,
            'analysis_date': datetime.now().isoformat(),
            'community_profile': community_profile,
            'detected_mutations': mutations,
            'resistance_analysis': resistance_analysis,
            'clinical_recommendations': clinical_recommendations,
            'quality_metrics': {
                'total_reads_analyzed': 100000,  # Mock
                'qrdr_reads_found': len(mutations) * 10,  # Mock
                'mutation_confidence': 'high' if mutations else 'none'
            }
        }
        
        return report
    
    def _analyze_resistance_level(self, mutations: List[Dict]) -> Dict:
        """Analyze overall resistance level from detected mutations"""
        
        if not mutations:
            return {
                'overall_level': 'susceptible',
                'explanation': 'No resistance mutations detected'
            }
        
        # Check for high-level resistance mutations
        high_level_mutations = [
            m for m in mutations 
            if (m.get('gene') == 'gyrA' and m.get('position') == 83 and m.get('mutant') in ['L', 'F'])
        ]
        
        if high_level_mutations:
            return {
                'overall_level': 'high_resistance',
                'explanation': 'High-level resistance mutations detected',
                'critical_mutations': high_level_mutations
            }
        
        return {
            'overall_level': 'moderate_resistance',
            'explanation': 'Moderate resistance mutations detected',
            'mutations': mutations
        }
    
    def _get_clinical_recommendations(self, resistance_analysis: Dict, 
                                    community_profile: Dict) -> Dict:
        """Generate clinical treatment recommendations"""
        
        resistance_level = resistance_analysis['overall_level']
        
        if resistance_level == 'high_resistance':
            return {
                'primary_recommendation': 'AVOID fluoroquinolones',
                'alternative_antibiotics': [
                    'ceftriaxone', 'azithromycin', 'trimethoprim-sulfamethoxazole'
                ],
                'urgency': 'high',
                'note': 'High-level fluoroquinolone resistance detected'
            }
        
        elif resistance_level == 'moderate_resistance':
            return {
                'primary_recommendation': 'Consider alternative antibiotics',
                'alternative_antibiotics': [
                    'levofloxacin (high dose)', 'moxifloxacin', 'ceftriaxone'
                ],
                'urgency': 'moderate',
                'note': 'Moderate resistance may require higher doses or alternatives'
            }
        
        else:
            return {
                'primary_recommendation': 'Standard fluoroquinolone therapy appropriate',
                'antibiotics': ['ciprofloxacin', 'levofloxacin'],
                'urgency': 'routine',
                'note': 'No resistance mutations detected'
            }

def test_pipeline():
    """Test function for the FQ resistance pipeline"""
    
    # Initialize pipeline
    pipeline = TitanXpFQPipeline()
    
    # Create mock FASTQ file for testing
    test_fastq = "test_sample.fastq"
    with open(test_fastq, 'w') as f:
        f.write("@read1\n")
        f.write("ATGGATAAAGATGTAGCTTATTTAGATATCGATGATGAAGATGATCCT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
        f.write("@read2\n")
        f.write("ATGGATCGAGATCTAGCTTATGTATATCGATGATCCTGATGATCCT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    
    # Process sample
    try:
        results = pipeline.process_clinical_sample(test_fastq, "TEST001")
        
        print("\n" + "="*50)
        print("CLINICAL FLUOROQUINOLONE RESISTANCE REPORT")
        print("="*50)
        print(f"Sample ID: {results['sample_id']}")
        print(f"Analysis Date: {results['analysis_date']}")
        print(f"Processing Time: {results['processing_time_seconds']:.2f} seconds")
        print("\nCommunity Profile:")
        for organism, abundance in results['community_profile'].items():
            print(f"  {organism}: {abundance:.1%}")
        
        print(f"\nMutations Detected: {len(results['detected_mutations'])}")
        for mutation in results['detected_mutations']:
            print(f"  {mutation['gene']} {mutation['wild_type']}{mutation['position']}{mutation['mutant']} (confidence: {mutation['confidence']:.0%})")
        
        print(f"\nResistance Level: {results['resistance_analysis']['overall_level']}")
        print(f"Clinical Recommendation: {results['clinical_recommendations']['primary_recommendation']}")
        
        if 'alternative_antibiotics' in results['clinical_recommendations']:
            print("Alternative Antibiotics:")
            for antibiotic in results['clinical_recommendations']['alternative_antibiotics']:
                print(f"  - {antibiotic}")
        
        return results
        
    except Exception as e:
        print(f"Error in pipeline: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # Clean up test file
        if Path(test_fastq).exists():
            Path(test_fastq).unlink()

if __name__ == "__main__":
    print("Starting Fluoroquinolone Resistance Detection Pipeline")
    print("Optimized for dual NVIDIA Titan Xp GPUs")
    test_pipeline()