#!/usr/bin/env python3
"""
simple_fq_diagnostic.py - Diagnose FQ resistance detection issues
Run this to understand what your pipeline should be finding
"""

import os
import sys
import json
from datetime import datetime

# Genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# FQ resistance positions (1-based amino acid positions)
FQ_RESISTANCE_SITES = {
    'gyrA': {
        83: {'wild_type': 'S', 'resistant': ['L', 'F', 'W'], 'codon_wt': 'TCG'},
        87: {'wild_type': 'D', 'resistant': ['N', 'G', 'Y'], 'codon_wt': 'GAC'}
    },
    'parC': {
        80: {'wild_type': 'S', 'resistant': ['I', 'R'], 'codon_wt': 'AGC'},
        84: {'wild_type': 'E', 'resistant': ['V', 'K', 'G'], 'codon_wt': 'GAA'}
    }
}

def translate_dna(dna_seq):
    """Translate DNA to protein"""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)

def create_test_sequences():
    """Create test sequences with known mutations"""
    
    # E. coli gyrA QRDR region (includes codons 67-106)
    # Position 83 (S) is at nucleotides 247-249 from ATG
    gyrA_wt_dna = (
        "ATGAGCGACCTTGCGAGAGAAATTACACCGGTAACATTGAGGATGGGTCGCAGCCGCTG"  # 1-60
        "GAGCAACTCGACGCGACCCTGCGTAAACTCGACGATCTCTTTGCCGATGAGCGTTTTAT"  # 61-120
        "GACGCGTTATGGTGATCAGCCAGGATCGGGCGGTAAAGCCGGAAACCATCTCGATGGCG"  # 121-180
        "GAGCAGGGGCTGTACGGTTTCTCGGTGTACGACCTCGGTTCGATTGGCCCGACCAACGG"  # 181-240
        "GATCTCGCACGGCGATCCTGCGCGACGGCAGCTATGAGCGTCTGATCAAAGATATGCGT"  # 241-300
        "ATGCGCTCAGCGTATGCGATTTCCGTTGGGCAATGAACAGTGGAACAACCTGTTCCGCT"  # 301-360
    )
    
    # Create S83L mutant (TCG -> TTG at position 247-249)
    gyrA_s83l = gyrA_wt_dna[:247] + "TTG" + gyrA_wt_dna[250:]
    
    # Create D87N mutant (GAC -> AAC at position 259-261)
    gyrA_d87n = gyrA_wt_dna[:259] + "AAC" + gyrA_wt_dna[262:]
    
    # parC wild type (S80 at 238-240, E84 at 250-252)
    parC_wt_dna = (
        "ATGAGCGAAATCGCGCCGTGCGTGACCGCCTGGAAGCGGTAGATCAGGAACGCTCGATG"  # 1-60
        "GCGGTACAGCAGGTTATCCACCTTCGACGCGGCGTTTACGCCATGGAAGAAGTTTGGCG"  # 61-120
        "GCGGTGATTGCCAGTGTAGAGGGGAAGAGTCGGTGGAAGATGCGCGTCGCCGTTTGCCA"  # 121-180
        "GAAGCGGCGGTGAAAGCCGTATCTGACATTTATAAAGAGCTGTATCTGACCAATATGAG"  # 181-240
        "CGACCATATTCAGCGGGAAGTTCTGCAAGCAGAACCGGCTAAAGCGATCAAAACCTTGC"  # 241-300
    )
    
    # Create S80I mutant (AGC -> ATC at 238-240)
    parC_s80i = parC_wt_dna[:238] + "ATC" + parC_wt_dna[241:]
    
    return [
        ("gyrA_WT", gyrA_wt_dna, "gyrA"),
        ("gyrA_S83L", gyrA_s83l, "gyrA"),
        ("gyrA_D87N", gyrA_d87n, "gyrA"),
        ("parC_WT", parC_wt_dna, "parC"),
        ("parC_S80I", parC_s80i, "parC")
    ]

def analyze_sequence(name, dna_seq, gene_name):
    """Analyze a sequence for FQ resistance mutations"""
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"Gene: {gene_name}")
    print(f"Length: {len(dna_seq)} bp")
    
    # Translate to protein
    protein = translate_dna(dna_seq)
    print(f"Protein length: {len(protein)} aa")
    
    # Check each resistance position
    resistance_found = []
    
    for codon_pos, site_info in FQ_RESISTANCE_SITES.get(gene_name, {}).items():
        # Convert to 0-based index
        aa_index = codon_pos - 1
        
        if aa_index < len(protein):
            observed_aa = protein[aa_index]
            dna_start = aa_index * 3
            observed_codon = dna_seq[dna_start:dna_start+3]
            
            print(f"\nPosition {codon_pos}:")
            print(f"  Expected WT: {site_info['wild_type']} ({site_info['codon_wt']})")
            print(f"  Observed: {observed_aa} ({observed_codon})")
            
            if observed_aa == site_info['wild_type']:
                print(f"  → Wild-type")
            elif observed_aa in site_info['resistant']:
                print(f"  → RESISTANT MUTATION: {site_info['wild_type']}{codon_pos}{observed_aa}")
                resistance_found.append(f"{site_info['wild_type']}{codon_pos}{observed_aa}")
            else:
                print(f"  → Novel mutation: {site_info['wild_type']}{codon_pos}{observed_aa}")
        else:
            print(f"\nPosition {codon_pos}: Not covered (sequence too short)")
    
    return resistance_found

def create_fastq_files():
    """Create test FASTQ files with known mutations"""
    
    # Create R1 file
    with open("test_fq_resistance_R1.fastq", "w") as f:
        sequences = create_test_sequences()
        
        for i, (name, seq, gene) in enumerate(sequences):
            # Write FASTQ record
            f.write(f"@read_{i}_{name}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write("I" * len(seq) + "\n")
    
    # Create R2 file (reverse complement)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    with open("test_fq_resistance_R2.fastq", "w") as f:
        sequences = create_test_sequences()
        
        for i, (name, seq, gene) in enumerate(sequences):
            # Reverse complement
            rc_seq = ''.join(complement[base] for base in seq[::-1])
            
            f.write(f"@read_{i}_{name}_RC\n")
            f.write(f"{rc_seq}\n")
            f.write("+\n")
            f.write("I" * len(rc_seq) + "\n")
    
    print("\nCreated test files:")
    print("  - test_fq_resistance_R1.fastq")
    print("  - test_fq_resistance_R2.fastq")

def diagnose_alignment_issue():
    """Diagnose common alignment issues"""
    
    print("\n" + "="*60)
    print("DIAGNOSTIC SUMMARY")
    print("="*60)
    
    print("""
Common issues in FQ resistance detection:

1. PROTEIN vs DNA COORDINATES:
   - Your pipeline uses 6-frame translation
   - Protein positions are 1-based (S83 means position 83)
   - The code checks positions 82 (0-based) for S83
   
2. GENE ID MAPPING:
   - Verify that gene_id=0 is gyrA
   - Verify that gene_id=1 is parC
   - Check your protein database metadata

3. ALIGNMENT COVERAGE:
   - Ensure reads are long enough to cover QRDR
   - Check that ref_start + match_length covers positions 82-86 (gyrA)
   
4. DEBUGGING STEPS:
   a) Enable DEBUG_TRANS in translated_search.cu
   b) Add the debug code provided earlier
   c) Look for these outputs:
      - "Triggering Smith-Waterman"
      - "MUTATION DEBUG"
      - "QRDR MUTATION"
""")

def main():
    print("FQ Resistance Detection Diagnostic Tool")
    print("="*60)
    
    # Create test sequences and analyze them
    sequences = create_test_sequences()
    
    all_mutations = []
    for name, seq, gene in sequences:
        mutations = analyze_sequence(name, seq, gene)
        if mutations:
            all_mutations.extend(mutations)
    
    # Create test FASTQ files
    create_fastq_files()
    
    # Show diagnostic summary
    diagnose_alignment_issue()
    
    # Create a summary report
    report = {
        "timestamp": datetime.now().isoformat(),
        "test_sequences": len(sequences),
        "mutations_found": all_mutations,
        "expected_mutations": ["S83L", "D87N", "S80I"],
        "recommendations": [
            "Enable DEBUG_TRANS=1 in translated_search.cu",
            "Add debug printf statements in enhanced_protein_kmer_match_kernel",
            "Verify gene_id mapping matches your database",
            "Check SW_SCORE_THRESHOLD (try lowering to 50.0f)"
        ]
    }
    
    with open("diagnostic_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    print(f"\n✓ Diagnostic complete. Check diagnostic_report.json")
    print(f"✓ Test FASTQ files created for pipeline testing")

if __name__ == "__main__":
    main()