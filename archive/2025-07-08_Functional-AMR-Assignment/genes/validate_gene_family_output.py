#!/usr/bin/env python3
"""
Validation script to verify gene_id and gene_family propagation
"""

import json
import csv
import sys
import os

# Expected mappings
EXPECTED_MAPPINGS = {
    1: {"gene_name": "blaKPC-2", "gene_family": "blaKPC", "drug_class": "Beta-lactamase"},
    2: {"gene_name": "qnrA1", "gene_family": "qnrA1", "drug_class": "Quinolone resistance"},
    3: {"gene_name": "vanA", "gene_family": "vanA", "drug_class": "Glycopeptide resistance"},
    4: {"gene_name": "mecA", "gene_family": "mecA", "drug_class": "Methicillin resistance"},
    5: {"gene_name": "aac6-Ib", "gene_family": "aac6", "drug_class": "Aminoglycoside resistance"}
}

def validate_json_output(json_file):
    """Validate JSON report output"""
    print(f"\n=== Validating JSON output: {json_file} ===")
    
    if not os.path.exists(json_file):
        print(f"ERROR: JSON file not found: {json_file}")
        return False
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    success = True
    gene_summaries = data.get('gene_summaries', {})
    
    print(f"Found {len(gene_summaries)} genes in JSON report")
    
    for gene_name, summary in gene_summaries.items():
        print(f"\nGene: {gene_name}")
        print(f"  Gene family: {summary.get('gene_family', 'NOT FOUND')}")
        print(f"  Drug class: {summary.get('drug_class', 'NOT FOUND')}")
        print(f"  Read count: {summary.get('read_count', 0)}")
        print(f"  Coverage: {summary.get('percent_coverage', 0):.1f}%")
        
        # Check if gene_family is present
        if 'gene_family' not in summary:
            print("  ERROR: gene_family field missing!")
            success = False
    
    return success

def validate_tsv_output(tsv_file):
    """Validate TSV abundance output"""
    print(f"\n=== Validating TSV output: {tsv_file} ===")
    
    if not os.path.exists(tsv_file):
        print(f"ERROR: TSV file not found: {tsv_file}")
        return False
    
    success = True
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # Check if gene_family column exists
        if 'gene_family' not in reader.fieldnames:
            print("ERROR: gene_family column not found in TSV!")
            print(f"Available columns: {reader.fieldnames}")
            return False
        
        print(f"TSV columns: {', '.join(reader.fieldnames)}")
        print("\nGene entries:")
        
        for row in reader:
            gene_name = row['gene_name']
            gene_family = row.get('gene_family', 'NOT FOUND')
            drug_class = row.get('drug_class', 'NOT FOUND')
            
            print(f"\n  {gene_name}:")
            print(f"    Family: {gene_family}")
            print(f"    Class: {drug_class}")
            print(f"    Reads: {row.get('read_count', 0)}")
            print(f"    Coverage: {row.get('percent_coverage', 0)}%")
            
            # Validate against expected values
            for gene_id, expected in EXPECTED_MAPPINGS.items():
                if gene_name == expected['gene_name']:
                    if gene_family != expected['gene_family']:
                        print(f"    ERROR: Expected family '{expected['gene_family']}', got '{gene_family}'")
                        success = False
                    else:
                        print(f"    ✓ Gene family correct")
                    break
    
    return success

def validate_database_metadata(json_file):
    """Validate protein database metadata"""
    print(f"\n=== Validating database metadata: {json_file} ===")
    
    if not os.path.exists(json_file):
        print(f"ERROR: Database metadata file not found: {json_file}")
        return False
    
    with open(json_file, 'r') as f:
        proteins = json.load(f)
    
    print(f"Found {len(proteins)} proteins in database")
    
    success = True
    for protein in proteins[:5]:  # Check first 5
        protein_id = protein.get('id', -1)
        gene_id = protein.get('gene_id', -1)
        gene_name = protein.get('gene_name', 'NOT FOUND')
        gene_family = protein.get('gene_family', 'NOT FOUND')
        
        print(f"\nProtein {protein_id}:")
        print(f"  Gene ID: {gene_id}")
        print(f"  Gene name: {gene_name}")
        print(f"  Gene family: {gene_family}")
        
        if 'gene_family' not in protein:
            print("  ERROR: gene_family field missing!")
            success = False
        
        # Validate against expected
        if gene_id in EXPECTED_MAPPINGS:
            expected = EXPECTED_MAPPINGS[gene_id]
            if gene_name != expected['gene_name']:
                print(f"  ERROR: Expected name '{expected['gene_name']}', got '{gene_name}'")
                success = False
            if gene_family != expected['gene_family']:
                print(f"  ERROR: Expected family '{expected['gene_family']}', got '{gene_family}'")
                success = False
    
    return success

def main():
    """Main validation function"""
    if len(sys.argv) < 2:
        print("Usage: python3 validate_gene_family_output.py <output_directory>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    
    print("=== Gene Family Validation Tool ===")
    print(f"Validating outputs in: {output_dir}")
    
    all_success = True
    
    # Validate database metadata
    db_json = os.path.join(output_dir, "../test_amr_db/protein_details.json")
    if os.path.exists(db_json):
        all_success &= validate_database_metadata(db_json)
    
    # Validate TSV output
    tsv_file = os.path.join(output_dir, "test_results_amr_abundance.tsv")
    if os.path.exists(tsv_file):
        all_success &= validate_tsv_output(tsv_file)
    
    # Validate JSON report
    json_file = os.path.join(output_dir, "test_results_clinical_amr_report.json")
    if os.path.exists(json_file):
        all_success &= validate_json_output(json_file)
    
    print("\n=== VALIDATION SUMMARY ===")
    if all_success:
        print("✓ All validations PASSED!")
        print("Gene families are correctly propagated through the pipeline.")
    else:
        print("✗ Some validations FAILED!")
        print("Check the errors above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()