# BioGPU Integrated Resistance Pipeline - Complete Guide

## Overview

The integrated resistance pipeline combines three complementary approaches for comprehensive fluoroquinolone resistance detection:

1. **Nucleotide k-mer screening** - Fast initial filtering
2. **6-frame translated search** - Protein-level mutation detection
3. **Pileup-based variant calling** - Population-level resistance frequencies

## Building the Integrated Database

### Step 1: Prepare Source Data

Ensure you have the following files:
- `data/fq_genes/` - Directory containing FASTA files of resistance genes
- `data/Known_Quinolone_Changes.csv` - Known resistance mutations

### Step 2: Copy the Database Builder Script

The `build_integrated_resistance_db.py` script from the documents needs to be placed in the correct location:

```bash
# Create the directory
mkdir -p src/python

# Copy the script content from the documents to:
# src/python/build_integrated_resistance_db.py

# Make it executable
chmod +x src/python/build_integrated_resistance_db.py
```

### Step 3: Run the Database Build Script

```bash
# Make the build script executable
chmod +x build_integrated_database.sh

# Run the build script
./build_integrated_database.sh
```

This will create the complete database structure in `data/integrated_resistance_db/`.

## Building the Pipeline Executable

### Step 1: Create Build Directory

```bash
mkdir -p build
cd build
```

### Step 2: Configure with CMake

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
```

### Step 3: Build the Integrated Pipeline

```bash
# Build just the integrated pipeline
make integrated_resistance_pipeline

# Or build everything
make -j$(nproc)
```

## Running the Integrated Pipeline

### Basic Usage

```bash
./integrated_resistance_pipeline \
    data/integrated_resistance_db/nucleotide \
    data/integrated_resistance_db/protein \
    data/integrated_resistance_db \
    reads_R1.fastq.gz \
    reads_R2.fastq.gz \
    output_prefix
```

### Input Parameters

1. **Nucleotide index path** - Directory containing k-mer index
2. **Protein database path** - Directory with protein sequences and index
3. **Resistance database path** - Root directory of integrated database
4. **R1 reads** - Forward reads (gzipped FASTQ)
5. **R2 reads** - Reverse reads (gzipped FASTQ)
6. **Output prefix** - Base name for output files

### Output Files

The pipeline generates:
- `output_prefix.h5` - HDF5 file with all alignment data
- `output_prefix.json` - Human-readable resistance calls
- `output_prefix_summary.json` - Pipeline statistics
- `output_prefix_diagnostic.txt` - Diagnostic report (if enabled)

## Configuration Options

The pipeline can be configured by modifying the source code:

```cpp
// In integrate_resistance_pipeline.cpp
pipeline.setConfig("min_alignment_score", 50.0f);
pipeline.setConfig("min_allele_frequency", 0.1f);
pipeline.setConfig("min_depth", 5.0f);
pipeline.setConfig("min_confidence", 0.8f);
```

## Interpreting Results

### JSON Output Structure

```json
{
  "sample": "sample_name",
  "total_reads": 1000000,
  "resistance_calls": [
    {
      "gene": "gyrA",
      "position": 83,
      "wildtype": "S",
      "mutant": "L",
      "frequency": 0.95,
      "depth": 100,
      "confidence": 0.99,
      "drug": "ciprofloxacin",
      "resistance_level": "high"
    }
  ]
}
```

### HDF5 Data Structure

The HDF5 file contains:
- `/alignments/` - Nucleotide alignments
- `/translated_search/` - Protein alignments
- `/resistance_calls/` - Variant calls
- `/metadata/` - Pipeline parameters

## Troubleshooting

### Common Issues

1. **Missing Python Scripts**
   - Ensure all scripts from the documents are copied to `src/python/`
   - Check file permissions (should be executable)

2. **CUDA Errors**
   - Verify GPU compute capability (needs 6.1 or higher)
   - Check available GPU memory
   - Update CUDA drivers if needed

3. **Database Build Failures**
   - Check input file formats
   - Ensure sufficient disk space
   - Verify Python dependencies: `pip install biopython pandas numpy`

4. **Low Detection Rates**
   - Adjust sensitivity parameters
   - Check read quality
   - Verify database completeness

### Performance Optimization

- **Batch Size**: Adjust in source code (default 10,000 reads)
- **GPU Memory**: Monitor with `nvidia-smi`
- **Enable/Disable Components**: Modify config flags for specific use cases

## Advanced Usage

### Running with Custom Parameters

Create a configuration file:

```bash
cat > pipeline_config.json << EOF
{
  "enable_nucleotide_search": true,
  "enable_protein_search": true,
  "enable_pileup_calling": true,
  "enable_smith_waterman": true,
  "min_alignment_score": 40.0,
  "min_allele_frequency": 0.05,
  "min_depth": 3,
  "batch_size": 20000
}
EOF
```

### Analyzing Results

```bash
# Generate HTML report
python3 src/python/resistance_profile_analyzer.py \
    output_prefix.h5 \
    output_prefix.json \
    --generate-html \
    --output resistance_report.html

# Create clinical interpretation
python3 src/python/clinical_interpreter.py \
    output_prefix.json \
    --guidelines CLSI \
    --output clinical_report.pdf
```

## Extending the Pipeline

### Adding New Resistance Genes

1. Add gene sequences to `data/fq_genes/`
2. Update `Known_Quinolone_Changes.csv` with mutations
3. Rebuild the database

### Supporting New Drug Classes

1. Modify `ResistanceMutation` structure in source
2. Update database builder to include new mutations
3. Add clinical interpretation rules

## Citation

If you use this pipeline, please cite:
- BioGPU: GPU-accelerated bioinformatics pipelines
- Method papers for each component approach

## Support

For issues or questions:
1. Check the diagnostic report for troubleshooting hints
2. Review GPU memory usage and errors
3. Validate database integrity
4. Contact the development team