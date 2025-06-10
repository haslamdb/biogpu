# BioGPU Archive - June 10, 2025
## Allele Frequency Analysis and Batch Processing Mode

This archive captures the state of the BioGPU project after implementing comprehensive allele frequency analysis and batch processing capabilities.

### Archive Date: June 10, 2025
### Version: 0.6.3

## Key Features in This Snapshot

### 1. Allele Frequency Analysis (v0.6.2)
- Comprehensive allele frequency reporting for all detected positions
- Configurable depth filtering with `--min-allele-depth` and `--min-report-depth`
- Tracks wildtype vs mutant amino acid distributions
- Essential for monitoring mutation prevalence in metagenomic samples
- Outputs detailed CSV with per-position coverage and mutation summaries

### 2. Batch Processing Support (v0.6.3)
- **NEW: CSV-based batch mode** for processing multiple samples
- Added `--csv` option to specify sample manifest file
- CSV parser library (`sample_csv_parser.cpp/h`) for flexible input handling
- Supports multiple path formats (absolute, relative, home-relative)
- `--dry-run` option for validation without processing
- `--stop-on-error` option for batch error handling
- Includes `generate_sample_csv.py` utility for creating CSV from directories

### 3. Automatic Sample Management (v0.6.3)
- Automatic sample name extraction from filenames
- Default output structure: `results/<sample_name>/`
- Simplified command line interface
- Configurable output directory with `--output-dir`

## Major Components

### Core Pipeline
- `runtime/kernels/resistance/clean_resistance_pipeline_main.cpp` - Main pipeline with batch processing
- `runtime/kernels/resistance/sample_csv_parser.cpp/h` - CSV parsing library
- `runtime/kernels/resistance/fq_mutation_reporter.cpp` - Enhanced mutation reporting
- `runtime/kernels/resistance/clinical_fq_report_generator.cpp` - Clinical report generation

### Python Utilities
- `runtime/kernels/resistance/generate_sample_csv.py` - Generate CSV from directory of FASTQ files
- `src/python/build_clean_dynamic_database.py` - Unified database builder

### Documentation
- `docs/workflow_documentation.md` - Comprehensive pipeline documentation
- `docs/bloom_sw_configuration_guide.md` - Performance optimization guide

## Command Examples

### Single Sample Mode
```bash
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --no-bloom \
    --min-allele-depth 10
```

### Batch Mode
```bash
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --no-bloom \
    --output-dir /results/batch_analysis
```

## Performance Characteristics
- Throughput: ~16,667-21,739 reads/second on NVIDIA TITAN Xp
- Bloom filter disabled by default (--no-bloom) for 6% speed improvement
- Smith-Waterman enabled for sensitivity (30% more matches)

## Known Issues
- Bloom filter shows 0% pass rate (non-critical - pipeline works without it)
- CSV parsing error for qnrB19 entries with "NA" position (non-critical)

## Build Instructions
```bash
mkdir build && cd build
cmake ..
make clean_resistance_pipeline
```

## Testing
- Example CSV provided: `runtime/kernels/resistance/example_csv.csv`
- Dry run mode available with `--dry-run` flag
- Validation errors reported before processing

---

This archive represents a major milestone in making BioGPU practical for large-scale clinical and research applications with batch processing and detailed allele frequency analysis.