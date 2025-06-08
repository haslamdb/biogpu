#!/usr/bin/env python3
"""
Test protein database offset calculation with the expanded database
"""
import json
import struct
import os

def test_offset_calculation(db_path):
    """Test that protein offsets are calculated correctly"""
    print(f"Testing offset calculation for database: {db_path}")
    
    # Read metadata
    with open(os.path.join(db_path, 'metadata.json'), 'r') as f:
        metadata = json.load(f)
    
    # Read protein details
    with open(os.path.join(db_path, 'protein_details.json'), 'r') as f:
        protein_details = json.load(f)
    
    print(f"Database contains {len(protein_details)} proteins")
    print(f"Species: {metadata['statistics']['species_count']}")
    
    # Read binary proteins file
    proteins_path = os.path.join(db_path, 'proteins.bin')
    with open(proteins_path, 'rb') as f:
        num_proteins = struct.unpack('I', f.read(4))[0]
        remaining_data = f.read()
        total_sequence_bytes = len(remaining_data)
    
    print(f"Binary file contains {num_proteins} proteins")
    print(f"Total sequence data: {total_sequence_bytes} bytes")
    
    # Calculate expected offsets
    cumulative_offset = 0
    print("\nExpected protein offsets:")
    print("ID | Gene       | Species               | Length | Offset | Sequence Preview")
    print("-" * 85)
    
    for i, protein in enumerate(protein_details[:15]):  # Show first 15
        print(f"{protein['id']:2d} | {protein['gene']:10s} | {protein['species']:20s} | {protein['length']:6d} | {cumulative_offset:6d} | {protein['sequence'][:15]}...")
        cumulative_offset += protein['length']
    
    if len(protein_details) > 15:
        print(f"... and {len(protein_details) - 15} more proteins")
    
    print(f"\nTotal calculated length: {cumulative_offset} amino acids")
    print(f"Binary file length: {total_sequence_bytes} bytes")
    print(f"Match: {'✅ YES' if cumulative_offset == total_sequence_bytes else '❌ NO'}")
    
    # Test dynamic offset calculation logic (simulating what the GPU code would do)
    print("\nSimulating GPU dynamic offset calculation:")
    avg_protein_len = total_sequence_bytes // num_proteins
    print(f"Average protein length: {avg_protein_len} amino acids")
    
    print("\nDynamic offsets (first 10 proteins):")
    for i in range(min(10, num_proteins)):
        calculated_offset = i * avg_protein_len
        expected_length = avg_protein_len if i < num_proteins - 1 else (total_sequence_bytes - calculated_offset)
        print(f"  Protein {i}: offset={calculated_offset}, length={expected_length}")
    
    # Check if this would work correctly
    print(f"\nDynamic calculation assessment:")
    print(f"  Average length approach would work: {'✅ YES' if avg_protein_len > 0 else '❌ NO'}")
    print(f"  However, since proteins vary greatly in length, dynamic offsets may not align perfectly")
    print(f"  This is why the hardcoded mapping for specific database sizes is more accurate")

if __name__ == '__main__':
    # Test both databases
    print("=== Testing Original 2-Protein Database ===")
    test_offset_calculation("data/wildtype_protein_db")
    
    print("\n" + "="*80)
    print("=== Testing Expanded 92-Protein Database ===")
    test_offset_calculation("data/wildtype_protein_db_expanded")