#!/usr/bin/env python3
"""
Check peptide matches and mutations for QRDR region alignments
"""

# Example QRDR alignment from diagnostic output:
# Position: 74-123 (0-based, length 50)
# Peptide: GKYHPHGDAAVYDTIVRMAQPFSLRYMLV

# Known QRDR positions for E. coli:
# gyrA: 83, 87 (1-based)
# parC: 80, 84 (1-based) 
# parE: 416, 420 (1-based)

alignment_start = 74  # 0-based from alignment
alignment_end = 123   # 0-based from alignment
alignment_length = 50

# Convert known positions to 0-based
gyrA_positions = [82, 86]  # 83, 87 in 1-based
parC_positions = [79, 83]  # 80, 84 in 1-based
parE_positions = [415, 419]  # 416, 420 in 1-based

print("QRDR Position Analysis")
print("======================")
print(f"Alignment covers protein positions: {alignment_start}-{alignment_end} (0-based)")
print(f"Alignment length: {alignment_length}")
print()

print("gyrA QRDR positions (0-based):")
for pos in gyrA_positions:
    if alignment_start <= pos <= alignment_end:
        relative_pos = pos - alignment_start
        print(f"  Position {pos} (1-based: {pos+1}) is COVERED at relative position {relative_pos}")
    else:
        print(f"  Position {pos} (1-based: {pos+1}) is NOT covered")

print()
print("Mutation position mapping:")
print("If mutation_positions[] contains relative positions within the alignment:")
print("  - Relative position 8 = global position 74+8 = 82 (1-based: 83) ✓")
print("  - Relative position 12 = global position 74+12 = 86 (1-based: 87) ✓")

print()
print("But the code is checking:")
print("  if (pm.mutation_positions[m] == 83 || pm.mutation_positions[m] == 87)")
print("This is incorrect! It should check relative positions 8 and 12, OR")
print("convert to global positions first: pm.ref_start + pm.mutation_positions[m]")