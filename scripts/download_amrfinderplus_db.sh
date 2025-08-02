#!/bin/bash
# Download AMRFinderPlus database files from NCBI
# These files are used for antimicrobial resistance gene detection

set -e  # Exit on error

# Configuration
BASE_URL="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
OUTPUT_DIR="${1:-./}"  # Default to current directory if no argument provided

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "================================================"
echo "AMRFinderPlus Database Downloader"
echo "================================================"
echo "Downloading to: $OUTPUT_DIR"
echo ""

# Function to download a file with progress and error handling
download_file() {
    local filename=$1
    local url="${BASE_URL}/${filename}"
    local output_path="${OUTPUT_DIR}/${filename}"
    
    echo "Downloading ${filename}..."
    
    if wget --progress=bar:force:noscroll "${url}" -O "${output_path}"; then
        echo "✓ Successfully downloaded ${filename}"
        echo "  Size: $(ls -lh "${output_path}" | awk '{print $5}')"
    else
        echo "✗ Failed to download ${filename}"
        return 1
    fi
    echo ""
}

# Download all required files
echo "Starting downloads..."
echo ""

# Main database files
download_file "AMRProt.fa"
download_file "AMR_CDS.fa"
download_file "ReferenceGeneCatalog.txt"
download_file "AMRProt-mutation.tsv"

# Optional: Download version info if available
echo "Checking for version information..."
if wget -q --spider "${BASE_URL}/version.txt" 2>/dev/null; then
    download_file "version.txt"
else
    echo "No version.txt file available"
fi

echo ""
echo "================================================"
echo "Download Summary"
echo "================================================"
echo "Location: $OUTPUT_DIR"
echo ""
echo "Files downloaded:"
ls -lh "$OUTPUT_DIR"/*.fa "$OUTPUT_DIR"/*.txt "$OUTPUT_DIR"/*.tsv 2>/dev/null | grep -E "(AMRProt|AMR_CDS|ReferenceGeneCatalog|version)" || true
echo ""
echo "Total size: $(du -sh "$OUTPUT_DIR" | cut -f1)"
echo ""

# Copy to data directory if we're in the genes directory
if [[ "$OUTPUT_DIR" == *"/runtime/kernels/genes"* ]] && [[ -d "../../../data" ]]; then
    echo "Copying AMRProt.fa to data directory as AMR_protein.fa..."
    cp "${OUTPUT_DIR}/AMRProt.fa" "../../../data/AMR_protein.fa"
    echo "✓ Copied to data/AMR_protein.fa"
fi

echo "================================================"
echo "Download complete!"
echo "================================================"