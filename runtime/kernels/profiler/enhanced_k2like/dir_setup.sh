#!/bin/bash
# Step 1: Create the new modular directory structure

echo "Creating modular directory structure for GPU Kraken Database Builder..."

# Navigate to the enhanced_k2like directory
cd runtime/kernels/profiler/enhanced_k2like

# Create the new modular directories
mkdir -p core
mkdir -p memory
mkdir -p processing
mkdir -p gpu
mkdir -p taxonomy
mkdir -p output
mkdir -p utils

echo "Directory structure created:"
echo "enhanced_k2like/"
echo "├── core/           # Main orchestration class"
echo "├── memory/         # GPU memory management"
echo "├── processing/     # File and FNA processing"
echo "├── gpu/            # CUDA kernels and device functions"
echo "├── taxonomy/       # NCBI taxonomy and phylogenetics"
echo "├── output/         # Database serialization"
echo "└── utils/          # Statistics and validation"

# Backup existing files
echo ""
echo "Creating backups of existing files..."
cp gpu_kraken_database_builder.cu gpu_kraken_database_builder.cu.backup
cp gpu_kraken_database_builder.h gpu_kraken_database_builder.h.backup

echo "Backup files created:"
echo "- gpu_kraken_database_builder.cu.backup"
echo "- gpu_kraken_database_builder.h.backup"

echo ""
echo "✓ Directory setup complete!"
echo "Next: Extract CUDA kernels to gpu/ directory"
