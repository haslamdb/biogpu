#!/usr/bin/env python3
"""
Test CuPy installation and basic fluoroquinolone resistance detection
"""

import numpy as np
try:
    import cupy as cp
    print("✓ CuPy successfully imported!")
    print(f"  CUDA version: {cp.cuda.runtime.runtimeGetVersion()}")
    print(f"  CuPy version: {cp.__version__}")
    print(f"  Available GPU memory: {cp.cuda.MemoryPool().free_bytes() / 1024**3:.2f} GB")
except ImportError as e:
    print(f"✗ CuPy import failed: {e}")
    print("  Please install with: pip install cupy-cuda12x")
    exit(1)

# Test data: Simulated gyrA sequences with and without S83L mutation
# Wild type: ...TCT... (Serine)
# Mutant:    ...TTG... (Leucine)

def create_test_sequences():
    """Create test sequences with known mutations"""
    # gyrA QRDR region (simplified)
    wild_type = "ATGGCAGATTCCGTTGATGGTGAATCTTCCGTTGATGGT"  # Contains TCT (Ser83)
    s83l_mutant = "ATGGCAGATTCCGTTGATGGTGAATTGTTCCGTTGATGGT"  # TCT -> TTG
    
    # Create batch of sequences
    sequences = []
    labels = []
    
    # Add 50 wild type
    for _ in range(50):
        sequences.append(wild_type)
        labels.append(0)
    
    # Add 30 S83L mutants
    for _ in range(30):
        sequences.append(s83l_mutant)
        labels.append(1)
    
    # Add 20 with random mutations (noise)
    import random
    bases = ['A', 'T', 'C', 'G']
    for _ in range(20):
        seq = list(wild_type)
        # Random mutation outside QRDR
        pos = random.choice([0, 1, 2, 35, 36, 37])
        seq[pos] = random.choice(bases)
        sequences.append(''.join(seq))
        labels.append(0)
    
    return sequences, labels

def sequence_to_gpu_array(sequences):
    """Convert DNA sequences to GPU array"""
    # Encode as integers: A=0, C=1, G=2, T=3
    encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    max_len = max(len(s) for s in sequences)
    encoded = np.zeros((len(sequences), max_len), dtype=np.int8)
    
    for i, seq in enumerate(sequences):
        for j, base in enumerate(seq):
            encoded[i, j] = encoding.get(base, -1)
    
    return cp.asarray(encoded)

def detect_s83l_mutation_gpu(encoded_sequences):
    """Detect S83L mutation using GPU"""
    # Position 26-28 in our test sequence corresponds to codon 83
    codon_pos = 26
    
    # Extract codon at position 83
    codon = encoded_sequences[:, codon_pos:codon_pos+3]
    
    # Check for S83L mutation
    # TCT (Ser) = [3,1,3], TTG (Leu) = [3,3,2]
    ser_codon = cp.array([3, 1, 3])
    leu_codon = cp.array([3, 3, 2])
    
    # Check each sequence
    is_wild_type = cp.all(codon == ser_codon, axis=1)
    is_s83l = cp.all(codon == leu_codon, axis=1)
    
    return is_s83l

def run_resistance_detection():
    """Run the test pipeline"""
    print("\n=== Fluoroquinolone Resistance Detection Test ===")
    
    # Create test data
    print("Creating test sequences...")
    sequences, true_labels = create_test_sequences()
    print(f"  Generated {len(sequences)} sequences")
    print(f"  True S83L mutants: {sum(true_labels)}")
    
    # Convert to GPU arrays
    print("\nConverting to GPU arrays...")
    gpu_sequences = sequence_to_gpu_array(sequences)
    print(f"  GPU array shape: {gpu_sequences.shape}")
    print(f"  GPU memory used: {gpu_sequences.nbytes / 1024:.2f} KB")
    
    # Detect mutations
    print("\nDetecting S83L mutations on GPU...")
    start = cp.cuda.Event()
    end = cp.cuda.Event()
    
    start.record()
    s83l_detected = detect_s83l_mutation_gpu(gpu_sequences)
    end.record()
    end.synchronize()
    
    gpu_time = cp.cuda.get_elapsed_time(start, end)
    
    # Results
    detected_count = int(cp.sum(s83l_detected))
    print(f"\nResults:")
    print(f"  S83L mutations detected: {detected_count}")
    print(f"  Detection accuracy: {(detected_count == sum(true_labels)) * 100:.1f}%")
    print(f"  GPU processing time: {gpu_time:.3f} ms")
    print(f"  Sequences processed per second: {len(sequences) / (gpu_time/1000):.0f}")
    
    # Show sample results
    print("\nSample results (first 10):")
    results = cp.asnumpy(s83l_detected)
    for i in range(10):
        status = "S83L MUTANT" if results[i] else "Wild Type"
        print(f"  Sequence {i+1}: {status} (True: {'S83L' if true_labels[i] else 'WT'})")

if __name__ == "__main__":
    run_resistance_detection()
    
    print("\n✓ CuPy is working correctly for bioinformatics applications!")
    print("\nNext steps:")
    print("1. Scale up to process millions of reads")
    print("2. Add detection for other QRDR mutations (D87N, parC S80I, etc.)")
    print("3. Implement parallel mapping to reference genomes")
    print("4. Build the full BioGPU pipeline!")
