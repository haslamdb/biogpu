#!/bin/bash
# run_fq_diagnostics.sh - Complete diagnostic workflow for FQ resistance detection

echo "========================================"
echo "FQ Resistance Detection Diagnostic Tool"
echo "========================================"
echo ""

# Step 1: Create diagnostic directory
DIAG_DIR="fq_diagnostic_$(date +%Y%m%d_%H%M%S)"
mkdir -p $DIAG_DIR
cd $DIAG_DIR

echo "Created diagnostic directory: $DIAG_DIR"

# Step 2: Copy diagnostic files
echo "Setting up diagnostic files..."
cat > diagnostic_fq_detection.cu << 'EOF'
// diagnostic_fq_detection.cu
// Diagnostic version to troubleshoot why FQ resistance mutations aren't being detected

#include <cuda_runtime.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>

// Genetic code for translation
__device__ const char CODON_TABLE[64] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
};

// Base to index mapping
__device__ int base_to_index(char base) {
    switch(base) {
        case 'T': case 't': return 0;
        case 'C': case 'c': return 1;
        case 'A': case 'a': return 2;
        case 'G': case 'g': return 3;
        default: return -1;
    }
}

// Translate codon to amino acid
__device__ char translate_codon(const char* codon) {
    int idx1 = base_to_index(codon[0]);
    int idx2 = base_to_index(codon[1]);
    int idx3 = base_to_index(codon[2]);
    
    if (idx1 < 0 || idx2 < 0 || idx3 < 0) return 'X'; // Unknown
    
    int codon_index = idx1 * 16 + idx2 * 4 + idx3;
    return CODON_TABLE[codon_index];
}

// Structure for QRDR positions we care about
struct QRDRPosition {
    const char* gene;
    int codon_position;      // 1-based codon position
    int gene_start;          // 0-based start position in reference
    char wild_type_aa;
    const char* resistant_aas;
};

// Key positions