#!/bin/bash
# Install script for BioGPU Unified Pipeline

set -e

EXECUTABLE="build_unified/bio_gpu_pipeline"
INSTALL_PATH="/usr/local/bin/bio_gpu_pipeline"

if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please run ./build_unified_minimal.sh first"
    exit 1
fi

echo "Installing BioGPU Unified Pipeline..."
echo "Source: $EXECUTABLE"
echo "Destination: $INSTALL_PATH"
echo ""
echo "This script needs sudo access to install to /usr/local/bin"
echo "You will be prompted for your password."
echo ""

# Copy and set permissions
sudo cp "$EXECUTABLE" "$INSTALL_PATH"
sudo chmod +x "$INSTALL_PATH"

echo "Installation complete!"
echo ""
echo "Test the installation with:"
echo "  bio_gpu_pipeline --help"
echo ""
echo "Example usage:"
echo "  bio_gpu_pipeline --use-multi-gpu \\"
echo "    --r1 /data/sample_R1.fastq.gz \\"
echo "    --r2 /data/sample_R2.fastq.gz \\"
echo "    --output-dir /data/results \\"
echo "    --reference-db /data/microbial_ref_db \\"
echo "    --resistance-db /data/quinolone_resistance_db \\"
echo "    --sample-id SAMPLE001 \\"
echo "    --progress-json"