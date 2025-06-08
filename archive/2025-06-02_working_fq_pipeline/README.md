# BioGPU Working FQ Pipeline - 2025-06-02 Snapshot

## Version
v0.4.1 - Enhanced with diagnostic reporting and mutation detection

## Build Instructions

```bash
# Create build directory
mkdir -p build
cd build

# Configure and build
cmake ..
make -j$(nproc)
```

## Working Command

```bash
./build/fq_pipeline_gpu data/fq_resistance_index \
    569_A_038_R1.fastq.gz \
    569_A_038_R2.fastq.gz \
    --enable-translated-search \
    --enable-smith-waterman \
    --protein-db data/wildtype_protein_db \
    --enable-diagnostic-reporting
```

## Key Features in This Version
- Enhanced mutation detection with wild-type comparison
- Comprehensive diagnostic reporting
- Alignment visualization and quality metrics
- Automatic troubleshooting recommendations
- Post-processing analysis tools
- Bloom filter integration for efficient k-mer screening
- Translated search capability
- Smith-Waterman alignment support

## Output Files
- Main output: JSON report with resistance predictions
- Diagnostic output: Detailed diagnostic text file (when --enable-diagnostic-reporting is used)

## Dependencies
- CUDA (tested with architectures 61, 70, 75, 80, 86)
- HDF5 (version 1.10.10)
- C++17 compiler
- CMake 3.18+

## Notes
This snapshot represents a fully working version tested with:
- Input files: 569_A_038_R1.fastq.gz and 569_A_038_R2.fastq.gz
- All diagnostic and alignment features enabled
- Successful resistance detection and reporting