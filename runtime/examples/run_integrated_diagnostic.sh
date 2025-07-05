#!/bin/bash

# Script to demonstrate integrated diagnostic analysis for clinical samples
# This shows how all three pipelines work together to guide treatment decisions

echo "=== BioGPU Integrated Diagnostic Workflow ==="
echo "Purpose: Comprehensive antimicrobial resistance profiling"
echo ""

# Check if sample file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <sample.fastq.gz> [patient_id]"
    echo ""
    echo "Example workflows:"
    echo "  UTI:      $0 urine_sample.fastq.gz UTI_2024_001"
    echo "  Sepsis:   $0 blood_culture.fastq.gz SEPSIS_2024_002" 
    echo "  Pneumonia: $0 bal_sample.fastq.gz HAP_2024_003"
    exit 1
fi

SAMPLE_FILE=$1
PATIENT_ID=${2:-"SAMPLE_$(date +%Y%m%d_%H%M%S)"}
OUTPUT_DIR="diagnostic_results/${PATIENT_ID}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create integrated configuration
cat > "$OUTPUT_DIR/integrated_config.json" << EOF
{
    "pipeline": {
        "enable_resistance": true,
        "enable_genes": true,
        "enable_profiler": true,
        "batch_size": 100000,
        "output_dir": "$OUTPUT_DIR"
    },
    "resistance": {
        "min_identity": 0.95,
        "min_coverage": 0.80,
        "calculate_allele_frequencies": true,
        "generate_clinical_report": true,
        "target_genes": ["gyrA", "gyrB", "parC", "parE", "grlA", "grlB"]
    },
    "genes": {
        "min_gene_coverage": 0.90,
        "min_gene_identity": 0.90,
        "amr_cds_path": "AMR_CDS.fa",
        "amr_prot_path": "AMRProt.fa"
    },
    "profiler": {
        "min_abundance": 0.01,
        "confidence_threshold": 0.80,
        "report_unclassified": false
    }
}
EOF

echo "Patient ID: $PATIENT_ID"
echo "Sample: $SAMPLE_FILE"
echo "Output: $OUTPUT_DIR"
echo ""

# Step 1: Quick pathogen screen
echo "Step 1: Rapid pathogen identification..."
./biogpu_diagnostic \
    --sample "$SAMPLE_FILE" \
    --analysis profiler \
    --output "$OUTPUT_DIR/pathogen_screen.txt" \
    --quick-mode &

PID1=$!

# Step 2: Resistance mutation detection
echo "Step 2: Detecting resistance mutations..."
./biogpu_diagnostic \
    --sample "$SAMPLE_FILE" \
    --analysis resistance \
    --output "$OUTPUT_DIR/resistance_mutations.txt" \
    --detect-qrdr &

PID2=$!

# Step 3: Resistance gene quantification
echo "Step 3: Quantifying resistance genes..."
./biogpu_diagnostic \
    --sample "$SAMPLE_FILE" \
    --analysis genes \
    --output "$OUTPUT_DIR/resistance_genes.txt" \
    --amr-panel &

PID3=$!

# Wait for all analyses to complete
wait $PID1 $PID2 $PID3

echo ""
echo "Step 4: Integrating results..."

# Run integrated analysis
./example_integrated_analysis \
    "$SAMPLE_FILE" \
    "$OUTPUT_DIR/integrated_config.json" \
    "$OUTPUT_DIR"

# Generate clinical summary
echo ""
echo "=== Clinical Summary for $PATIENT_ID ==="
echo ""

# Parse key findings (simplified example)
if [ -f "$OUTPUT_DIR/integrated_clinical_report.html" ]; then
    echo "✓ Integrated analysis complete"
    echo ""
    
    # Extract key information (would use proper parsing in practice)
    echo "Detected Pathogens:"
    grep -A2 "Primary pathogen" "$OUTPUT_DIR/integrated_clinical_report.html" | \
        sed 's/<[^>]*>//g' | grep -v "^$"
    
    echo ""
    echo "Resistance Status:"
    grep -E "resistance detected|No significant" "$OUTPUT_DIR/integrated_clinical_report.html" | \
        sed 's/<[^>]*>//g' | sed 's/⚠️/WARNING:/' | sed 's/✓/OK:/'
    
    echo ""
    echo "Key Findings:"
    # Check for specific resistance patterns
    if grep -q "gyrA_S83L" "$OUTPUT_DIR/resistance_mutations.txt" 2>/dev/null; then
        echo "• Fluoroquinolone resistance detected (gyrA mutation)"
        echo "  → Avoid ciprofloxacin/levofloxacin"
    fi
    
    if grep -q "CTX-M" "$OUTPUT_DIR/resistance_genes.txt" 2>/dev/null; then
        echo "• ESBL producer detected"
        echo "  → Use carbapenems for serious infections"
    fi
    
    if grep -q "mecA" "$OUTPUT_DIR/resistance_genes.txt" 2>/dev/null; then
        echo "• MRSA detected"
        echo "  → Use vancomycin or alternative anti-MRSA agents"
    fi
fi

echo ""
echo "=== Output Files ==="
echo "Clinical Report: $OUTPUT_DIR/integrated_clinical_report.html"
echo "Detailed Results: $OUTPUT_DIR/"
echo "  - Pathogen profile: pathogen_screen.txt"
echo "  - Resistance mutations: resistance_mutations.txt"
echo "  - Resistance genes: resistance_genes.txt"
echo ""

# Generate one-page clinical summary
cat > "$OUTPUT_DIR/clinical_summary.txt" << EOF
PATIENT ID: $PATIENT_ID
DATE: $(date)
SAMPLE TYPE: Metagenomic sequencing

RAPID DIAGNOSTIC SUMMARY
========================

1. PATHOGENS IDENTIFIED:
$(head -10 "$OUTPUT_DIR/pathogen_screen.txt" 2>/dev/null | grep -E "coli|aureus|pneumoniae" || echo "See detailed report")

2. RESISTANCE DETECTED:
$(grep -E "DETECTED|Resistant" "$OUTPUT_DIR/resistance_mutations.txt" 2>/dev/null | head -5 || echo "No high-confidence resistance mutations")

3. TREATMENT RECOMMENDATIONS:
$(grep -A3 "Treatment Recommendations" "$OUTPUT_DIR/integrated_clinical_report.html" 2>/dev/null | sed 's/<[^>]*>//g' | grep -v "^$" | head -10 || echo "See clinical report")

4. CLINICAL ACTION:
- Review integrated clinical report for detailed findings
- Consider local antibiogram data
- Adjust therapy based on clinical response

Generated by BioGPU Diagnostic Pipeline v2.0
This is a diagnostic aid - clinical correlation required
EOF

echo "Clinical summary saved to: $OUTPUT_DIR/clinical_summary.txt"
echo ""
echo "Analysis complete in $(($SECONDS / 60)) minutes $(($SECONDS % 60)) seconds"