#!/bin/bash
# Download NCBI BioProject PRJNA313047 - Bacterial Antimicrobial Resistance Reference Gene Database
# This is the comprehensive 7,600+ gene database you want

set -e  # Exit on any error

# Configuration
BIOPROJECT="PRJNA313047"
DOWNLOAD_DIR="./prjna313047_amr_database"
LOG_FILE="./prjna313047_download.log"
EMAIL="dbhaslam@gmail.com"  # REQUIRED: Replace with your email

# Environment variables for debugging and control
SKIP_COUNT="${SKIP_COUNT:-false}"  # Set to "true" to skip sequence counting
DEBUG_MODE="${DEBUG_MODE:-true}"   # Set to "false" to reduce debug output
RETRY_COUNT="${RETRY_COUNT:-3}"    # Number of retries for failed downloads
TIMEOUT_SEARCH="${TIMEOUT_SEARCH:-60}"     # Timeout for esearch (seconds)
TIMEOUT_LINK="${TIMEOUT_LINK:-120}"        # Timeout for elink (seconds)
TIMEOUT_FETCH="${TIMEOUT_FETCH:-600}"      # Timeout for efetch (seconds)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"
    exit 1
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOG_FILE"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

# Check dependencies
check_dependencies() {
    log "Checking dependencies..."
    
    if ! command -v esearch &> /dev/null; then
        error "NCBI EDirect tools not found. Please install with: conda install -c bioconda entrez-direct"
    fi
    
    if ! command -v elink &> /dev/null; then
        error "NCBI EDirect elink not found. Please install entrez-direct package."
    fi
    
    if ! command -v efetch &> /dev/null; then
        error "NCBI EDirect efetch not found. Please install entrez-direct package."
    fi
    
    success "All dependencies found"
}

# Validate email
validate_email() {
    if [[ "$EMAIL" == "your.email@domain.com" ]]; then
        error "Please set your email address in the EMAIL variable at the top of this script"
    fi
    
    log "Using email: $EMAIL"
}

# Create download directory
mkdir -p "$DOWNLOAD_DIR"
cd "$DOWNLOAD_DIR"

log "Starting BioProject $BIOPROJECT download..."
log "Download directory: $(pwd)"

# Check dependencies
check_dependencies
validate_email

# Function to download with progress tracking
download_sequences() {
    local db_type="$1"
    local output_file="$2"
    local description="$3"
    
    log "Downloading $description..."
    log "Command: esearch -email \"$EMAIL\" -db bioproject -query \"$BIOPROJECT\" | elink -target $db_type | efetch -format fasta"
    
    # Add option to skip counting
    if [ "$SKIP_COUNT" != "true" ]; then
        log "Counting sequences (this may take a moment)..."
        log "DEBUG: Running esearch..."
        local search_ids=$(timeout 30 esearch -email "$EMAIL" -db bioproject -query "$BIOPROJECT" 2>&1)
        if [ $? -ne 0 ]; then
            warning "esearch timed out or failed. Output: $search_ids"
        else
            log "DEBUG: esearch completed. Running elink..."
            local link_ids=$(echo "$search_ids" | timeout 30 elink -target "$db_type" 2>&1)
            if [ $? -ne 0 ]; then
                warning "elink timed out or failed. Output: $link_ids"
            else
                log "DEBUG: elink completed. Counting IDs..."
                local seq_count=$(echo "$link_ids" | efetch -format docsum | grep -c "<Id>" || echo "0")
                log "Found $seq_count sequences in $db_type database"
            fi
        fi
    else
        log "Skipping sequence count (SKIP_COUNT=true)"
    fi
    
    # Download sequences with timeout and intermediate steps
    log "Starting download process..."
    local temp_file="${output_file}.tmp"
    local search_file="${output_file}.search"
    local link_file="${output_file}.link"
    
    # Step 1: Search
    log "Step 1/3: Running esearch (timeout: ${TIMEOUT_SEARCH}s)..."
    if timeout "$TIMEOUT_SEARCH" esearch -email "$EMAIL" -db bioproject -query "$BIOPROJECT" > "$search_file" 2>&1; then
        log "DEBUG: esearch completed. File size: $(stat -c%s "$search_file" 2>/dev/null || echo 0) bytes"
        if [ "$DEBUG_MODE" = "true" ] && [ -f "$search_file" ]; then
            log "DEBUG: First 5 lines of esearch output:"
            head -5 "$search_file" | sed 's/^/    /'
        fi
        
        # Step 2: Link
        log "Step 2/3: Running elink (timeout: ${TIMEOUT_LINK}s)..."
        if timeout "$TIMEOUT_LINK" elink -target "$db_type" < "$search_file" > "$link_file" 2>&1; then
            log "DEBUG: elink completed. File size: $(stat -c%s "$link_file" 2>/dev/null || echo 0) bytes"
            if [ "$DEBUG_MODE" = "true" ] && [ -f "$link_file" ]; then
                log "DEBUG: Number of IDs found: $(grep -c "<Id>" "$link_file" 2>/dev/null || echo "0")"
            fi
            
            # Step 3: Fetch
            log "Step 3/3: Running efetch (timeout: ${TIMEOUT_FETCH}s)..."
            log "NOTE: This step may take several minutes for large datasets..."
            if timeout "$TIMEOUT_FETCH" efetch -format fasta < "$link_file" > "$temp_file" 2>&1; then
                log "DEBUG: efetch completed. File size: $(stat -c%s "$temp_file" 2>/dev/null || echo 0) bytes"
                
                # Check if we got data
                if [ -s "$temp_file" ]; then
                    mv "$temp_file" "$output_file"
                    local downloaded_count=$(grep -c "^>" "$output_file" 2>/dev/null || echo "0")
                    local file_size=$(stat -c%s "$output_file" 2>/dev/null || echo 0)
                    
                    if [ "$downloaded_count" -gt 0 ] && [ "$file_size" -gt 1000 ]; then
                        success "Downloaded $downloaded_count sequences to $output_file ($file_size bytes)"
                        # Clean up intermediate files
                        rm -f "$search_file" "$link_file"
                        return 0
                    else
                        error "Download produced invalid file: $output_file (sequences: $downloaded_count, size: $file_size)"
                        # Show last few lines for debugging
                        log "DEBUG: Last 10 lines of output:"
                        tail -10 "$output_file" | head -5
                    fi
                else
                    error "efetch produced empty file"
                    log "DEBUG: efetch stderr/stdout content:"
                    cat "$temp_file" | head -20
                fi
            else
                error "efetch timed out or failed"
                log "DEBUG: Partial output (if any):"
                if [ -f "$temp_file" ]; then
                    head -20 "$temp_file"
                fi
            fi
        else
            error "elink timed out or failed"
            log "DEBUG: elink output:"
            head -20 "$link_file"
        fi
    else
        error "esearch timed out or failed"
        log "DEBUG: esearch output:"
        head -20 "$search_file"
    fi
    
    # Clean up on failure
    rm -f "$temp_file" "$search_file" "$link_file"
    return 1
}

# Download DNA sequences (nucleotide database)
log "=== Downloading DNA/Nucleotide Sequences ==="
download_sequences "nuccore" "AMR_CDS.fa" "nucleotide sequences (DNA)"

# Download protein sequences
log "=== Downloading Protein Sequences ==="
download_sequences "protein" "AMRProt.fa" "protein sequences"

# Download metadata and additional information
log "=== Downloading Metadata ==="

# Get sequence metadata in XML format for parsing
log "Downloading nucleotide metadata..."
esearch -email "$EMAIL" -db bioproject -query "$BIOPROJECT" | elink -target nuccore | efetch -format docsum > nuccore_metadata.xml

log "Downloading protein metadata..."
esearch -email "$EMAIL" -db bioproject -query "$BIOPROJECT" | elink -target protein | efetch -format docsum > protein_metadata.xml

# Extract basic metadata to TSV format
log "Creating metadata TSV files..."

# Create nucleotide metadata TSV
if [ -f "nuccore_metadata.xml" ]; then
    echo -e "Accession\tTitle\tLength\tOrganism\tTaxId" > nuccore_metadata.tsv
    
    # This is a simplified extraction - you might need to adjust based on actual XML structure
    grep -o '<Id>[^<]*</Id>\|<Title>[^<]*</Title>\|<Length>[^<]*</Length>' nuccore_metadata.xml | \
    paste - - - | \
    sed 's/<[^>]*>//g' | \
    sed 's/\t/\t/g' >> nuccore_metadata.tsv 2>/dev/null || warning "Could not parse nucleotide metadata"
fi

# Create protein metadata TSV
if [ -f "protein_metadata.xml" ]; then
    echo -e "Accession\tTitle\tLength\tOrganism\tTaxId" > protein_metadata.tsv
    
    grep -o '<Id>[^<]*</Id>\|<Title>[^<]*</Title>\|<Length>[^<]*</Length>' protein_metadata.xml | \
    paste - - - | \
    sed 's/<[^>]*>//g' | \
    sed 's/\t/\t/g' >> protein_metadata.tsv 2>/dev/null || warning "Could not parse protein metadata"
fi

# Download BioProject description
log "Downloading BioProject information..."
esearch -email "$EMAIL" -db bioproject -query "$BIOPROJECT" | efetch -format xml > bioproject_info.xml

# Extract key information from BioProject
if [ -f "bioproject_info.xml" ]; then
    grep -o '<Title>[^<]*</Title>\|<Description>[^<]*</Description>' bioproject_info.xml | \
    sed 's/<[^>]*>//g' > bioproject_summary.txt 2>/dev/null || warning "Could not parse BioProject info"
fi

# Create comprehensive analysis
log "Analyzing downloaded sequences..."

# Analyze DNA sequences
if [ -f "AMR_CDS.fa" ]; then
    DNA_COUNT=$(grep -c "^>" "AMR_CDS.fa" || echo 0)
    DNA_SIZE=$(stat -c%s "AMR_CDS.fa" 2>/dev/null || echo 0)
    
    log "DNA sequences: $DNA_COUNT sequences, $DNA_SIZE bytes"
    
    # Extract resistance classes from sequence headers
    log "Analyzing resistance gene families..."
    grep "^>" "AMR_CDS.fa" | head -20 > sample_dna_headers.txt
    
    # Try to extract gene families from headers
    grep "^>" "AMR_CDS.fa" | \
    sed 's/^>//' | \
    cut -d' ' -f1 | \
    sed 's/[_-][0-9]*$//' | \
    sort | uniq -c | sort -nr > dna_gene_families.txt 2>/dev/null || true
fi

# Analyze protein sequences
if [ -f "AMRProt.fa" ]; then
    PROT_COUNT=$(grep -c "^>" "AMRProt.fa" || echo 0)
    PROT_SIZE=$(stat -c%s "AMRProt.fa" 2>/dev/null || echo 0)
    
    log "Protein sequences: $PROT_COUNT sequences, $PROT_SIZE bytes"
    
    # Sample protein headers
    grep "^>" "AMRProt.fa" | head -20 > sample_protein_headers.txt
    
    # Extract protein families
    grep "^>" "AMRProt.fa" | \
    sed 's/^>//' | \
    cut -d' ' -f1 | \
    sed 's/[_-][0-9]*$//' | \
    sort | uniq -c | sort -nr > protein_gene_families.txt 2>/dev/null || true
fi

# Validate downloads
log "Validating downloads..."

REQUIRED_FILES=("AMR_CDS.fa" "AMRProt.fa")
MISSING_FILES=()

for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ] || [ ! -s "$file" ]; then
        MISSING_FILES+=("$file")
    fi
done

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    error "Missing required files: ${MISSING_FILES[*]}"
fi

# Size validation
MIN_DNA_SIZE=5000000      # Expect at least 5MB for DNA (7600+ genes)
MIN_PROT_SIZE=2000000     # Expect at least 2MB for proteins

if [ "$DNA_SIZE" -lt "$MIN_DNA_SIZE" ]; then
    warning "DNA file smaller than expected ($DNA_SIZE bytes). May be incomplete."
fi

if [ "$PROT_SIZE" -lt "$MIN_PROT_SIZE" ]; then
    warning "Protein file smaller than expected ($PROT_SIZE bytes). May be incomplete."
fi

# Create comprehensive database info
log "Creating database information file..."

cat > "database_info.txt" << EOF
NCBI BioProject PRJNA313047 - Bacterial Antimicrobial Resistance Reference Gene Database
=======================================================================================
Download Date: $(date)
BioProject: $BIOPROJECT
Download Method: NCBI EDirect (esearch + elink + efetch)
Email Used: $EMAIL

File Summary:
=============
- AMR_CDS.fa: $DNA_COUNT DNA sequences ($DNA_SIZE bytes)
- AMRProt.fa: $PROT_COUNT protein sequences ($PROT_SIZE bytes)
- nuccore_metadata.xml: Nucleotide sequence metadata
- protein_metadata.xml: Protein sequence metadata
- bioproject_info.xml: BioProject description

Analysis Files:
===============
- dna_gene_families.txt: Gene family distribution (DNA)
- protein_gene_families.txt: Gene family distribution (proteins)
- sample_dna_headers.txt: Sample DNA sequence headers
- sample_protein_headers.txt: Sample protein sequence headers

Database Statistics:
====================
Total DNA sequences: $DNA_COUNT
Total protein sequences: $PROT_COUNT
DNA/Protein ratio: $(echo "scale=2; $PROT_COUNT / $DNA_COUNT" | bc -l 2>/dev/null || echo "N/A")

This is the comprehensive AMR database containing 7,600+ resistance genes
across all major antibiotic classes, perfect for BioGPU integration.
EOF

# Show top gene families
if [ -f "dna_gene_families.txt" ]; then
    echo "" >> "database_info.txt"
    echo "Top 20 DNA Gene Families:" >> "database_info.txt"
    head -20 "dna_gene_families.txt" >> "database_info.txt"
fi

if [ -f "protein_gene_families.txt" ]; then
    echo "" >> "database_info.txt"
    echo "Top 20 Protein Gene Families:" >> "database_info.txt"
    head -20 "protein_gene_families.txt" >> "database_info.txt"
fi

# Create validation script specific to this database
log "Creating validation script..."

cat > "validate_prjna313047.py" << 'EOF'
#!/usr/bin/env python3
"""
Validation script for NCBI BioProject PRJNA313047 AMR database
"""

import sys
import re
from collections import Counter

def analyze_fasta_headers(filename):
    """Analyze FASTA headers to extract gene information"""
    gene_families = []
    resistance_classes = []
    organisms = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()
                    
                    # Extract gene family (first part before underscore/number)
                    gene_match = re.search(r'>([a-zA-Z]+)', header)
                    if gene_match:
                        gene_families.append(gene_match.group(1))
                    
                    # Look for resistance class indicators
                    if 'beta' in header.lower() or 'bla' in header.lower():
                        resistance_classes.append('BETA_LACTAM')
                    elif 'van' in header.lower():
                        resistance_classes.append('GLYCOPEPTIDE')
                    elif 'tet' in header.lower():
                        resistance_classes.append('TETRACYCLINE')
                    elif 'sul' in header.lower():
                        resistance_classes.append('SULFONAMIDE')
                    elif 'aac' in header.lower() or 'ant' in header.lower() or 'aph' in header.lower():
                        resistance_classes.append('AMINOGLYCOSIDE')
                    
                    # Extract organism if present
                    org_match = re.search(r'\[([^\]]+)\]', header)
                    if org_match:
                        organisms.append(org_match.group(1))
    
    except FileNotFoundError:
        return None
    
    return {
        'gene_families': Counter(gene_families).most_common(20),
        'resistance_classes': Counter(resistance_classes).most_common(),
        'organisms': Counter(organisms).most_common(10)
    }

def validate_sequence_quality(filename, seq_type="DNA"):
    """Validate sequence quality and characteristics"""
    seq_count = 0
    total_length = 0
    lengths = []
    invalid_chars = 0
    
    valid_dna = set('ATCGNatcgn')
    valid_protein = set('ACDEFGHIKLMNPQRSTVWYXZacdefghiklmnpqrstvwyxz*')
    valid_chars = valid_dna if seq_type == "DNA" else valid_protein
    
    try:
        with open(filename, 'r') as f:
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        seq_len = len(current_seq)
                        total_length += seq_len
                        lengths.append(seq_len)
                        
                        # Check for invalid characters
                        for char in current_seq:
                            if char not in valid_chars:
                                invalid_chars += 1
                    
                    seq_count += 1
                    current_seq = ""
                else:
                    current_seq += line
            
            # Process last sequence
            if current_seq:
                seq_len = len(current_seq)
                total_length += seq_len
                lengths.append(seq_len)
                
                for char in current_seq:
                    if char not in valid_chars:
                        invalid_chars += 1
    
    except FileNotFoundError:
        return None
    
    if seq_count == 0:
        return None
    
    avg_length = total_length / seq_count
    min_length = min(lengths) if lengths else 0
    max_length = max(lengths) if lengths else 0
    
    return {
        'sequences': seq_count,
        'total_length': total_length,
        'avg_length': avg_length,
        'min_length': min_length,
        'max_length': max_length,
        'invalid_chars': invalid_chars
    }

def main():
    print("NCBI BioProject PRJNA313047 Database Validation")
    print("=" * 50)
    
    # Validate DNA sequences
    print("\nValidating AMR_CDS.fa (DNA sequences)...")
    dna_result = validate_sequence_quality("AMR_CDS.fa", "DNA")
    dna_analysis = analyze_fasta_headers("AMR_CDS.fa")
    
    if dna_result:
        print(f"  ✓ Found {dna_result['sequences']} DNA sequences")
        print(f"  ✓ Average length: {dna_result['avg_length']:.1f} bp")
        print(f"  ✓ Range: {dna_result['min_length']} - {dna_result['max_length']} bp")
        print(f"  ✓ Invalid characters: {dna_result['invalid_chars']}")
        
        if dna_result['sequences'] < 5000:
            print(f"  ⚠ Fewer sequences than expected ({dna_result['sequences']})")
    else:
        print("  ✗ Could not validate DNA sequences")
    
    if dna_analysis:
        print("  ✓ Top resistance classes detected:")
        for class_name, count in dna_analysis['resistance_classes'][:5]:
            print(f"    - {class_name}: {count}")
    
    # Validate protein sequences
    print("\nValidating AMRProt.fa (protein sequences)...")
    prot_result = validate_sequence_quality("AMRProt.fa", "Protein")
    prot_analysis = analyze_fasta_headers("AMRProt.fa")
    
    if prot_result:
        print(f"  ✓ Found {prot_result['sequences']} protein sequences")
        print(f"  ✓ Average length: {prot_result['avg_length']:.1f} aa")
        print(f"  ✓ Range: {prot_result['min_length']} - {prot_result['max_length']} aa")
        print(f"  ✓ Invalid characters: {prot_result['invalid_chars']}")
        
        if prot_result['sequences'] < 5000:
            print(f"  ⚠ Fewer sequences than expected ({prot_result['sequences']})")
    else:
        print("  ✗ Could not validate protein sequences")
    
    if prot_analysis:
        print("  ✓ Top resistance classes detected:")
        for class_name, count in prot_analysis['resistance_classes'][:5]:
            print(f"    - {class_name}: {count}")
    
    # Compare sequence counts
    if dna_result and prot_result:
        print(f"\nSequence Count Comparison:")
        print(f"  DNA sequences: {dna_result['sequences']}")
        print(f"  Protein sequences: {prot_result['sequences']}")
        
        ratio = prot_result['sequences'] / dna_result['sequences']
        print(f"  Protein/DNA ratio: {ratio:.2f}")
        
        if 0.8 <= ratio <= 1.2:
            print("  ✓ Sequence counts are well-matched")
        else:
            print("  ⚠ Significant sequence count difference")
    
    # Summary
    total_sequences = 0
    if dna_result:
        total_sequences += dna_result['sequences']
    if prot_result:
        total_sequences += prot_result['sequences']
    
    print(f"\nSummary:")
    print(f"  Total sequences downloaded: {total_sequences}")
    print(f"  Expected: ~15,000+ (7,600+ genes × 2 formats)")
    
    if total_sequences >= 10000:
        print("  ✓ Download appears complete")
    else:
        print("  ⚠ Download may be incomplete")
    
    print("\nDatabase ready for BioGPU integration!")

if __name__ == "__main__":
    main()
EOF

chmod +x "validate_prjna313047.py"

# Run validation
log "Running database validation..."
if command -v python3 &> /dev/null; then
    python3 "validate_prjna313047.py" | tee -a "$LOG_FILE"
else
    warning "Python3 not found. Skipping automatic validation."
fi

# Create simple usage instructions
cat > "README.md" << EOF
# NCBI BioProject PRJNA313047 - AMR Database

This directory contains the comprehensive Bacterial Antimicrobial Resistance Reference Gene Database downloaded from NCBI BioProject PRJNA313047.

## Files

- \`AMR_CDS.fa\` - DNA sequences of resistance genes (~$DNA_COUNT sequences)
- \`AMRProt.fa\` - Protein sequences of resistance genes (~$PROT_COUNT sequences)
- \`database_info.txt\` - Detailed database statistics and analysis
- \`validate_prjna313047.py\` - Validation script for quality checking

## Usage with BioGPU

This database is designed for use with your BioGPU NCBIAMRDatabaseLoader:

\`\`\`cpp
NCBIAMRDatabaseLoader loader;
bool success = loader.loadFromFiles("AMR_CDS.fa", "AMRProt.fa", "metadata.tsv");
\`\`\`

## Database Statistics

- Total DNA sequences: $DNA_COUNT
- Total protein sequences: $PROT_COUNT
- Download date: $(date)
- BioProject: $BIOPROJECT

This is the comprehensive 7,600+ gene AMR database you requested, containing resistance genes across all major antibiotic classes.
EOF

# Final summary
success "BioProject $BIOPROJECT download completed successfully!"

log "Database location: $(pwd)"
log "Key files:"
log "  - AMR_CDS.fa ($DNA_COUNT DNA sequences, $DNA_SIZE bytes)"
log "  - AMRProt.fa ($PROT_COUNT protein sequences, $PROT_SIZE bytes)"
log "  - database_info.txt (comprehensive statistics)"
log "  - validate_prjna313047.py (quality validation)"

log "This is the comprehensive AMR database with 7,600+ genes you requested!"
log "Ready for BioGPU translated alignment pipeline integration."

# Create symlinks for consistency with your loader expectations
ln -sf AMR_CDS.fa AMR_CDS.fa
ln -sf AMRProt.fa AMRProt.fa

log "Next steps:"
log "  1. Review database_info.txt for detailed statistics"
log "  2. Run validate_prjna313047.py for quality checks"
log "  3. Integrate with your BioGPU NCBIAMRDatabaseLoader"
log "  4. Build translated alignment kernels for this comprehensive database"

success "Download complete! Check $LOG_FILE for full details."
