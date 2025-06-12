# BioGPU Workflow Documentation

## Overview
This document tracks the complete BioGPU pipeline for GPU-accelerated fluoroquinolone resistance detection from metagenomic data.

**Last Updated**: June 11 2025  
**Pipeline Version**: 0.7.0  
**Status**: Fully functional FQ resistance detection with clinical reporting and GPU-accelerated microbiome profiler with paired-end support

---

## 📁 Project Structure

```
biogpu/
├── data/
│   ├── Known_Quinolone_Changes.csv      # Source resistance mutations from NCBI
│   ├── fq_genes/                        # Downloaded resistance gene sequences by species
│   │   ├── Escherichia_coli/           
│   │   ├── Enterococcus_faecium/       
│   │   └── [other species]/            
│   ├── integrated_clean_db/             # ✅ ACTIVELY USED integrated database
│   │   ├── nucleotide/                  # Nucleotide k-mer index
│   │   │   ├── kmer_index.bin          # Binary k-mer hash table (15-mers)
│   │   │   ├── sequences.bin           # Reference sequences
│   │   │   └── metadata.json           # Index metadata
│   │   ├── protein/                     # Protein database
│   │   │   ├── proteins.bin            # Binary protein sequences
│   │   │   ├── protein_kmers.bin       # 5-mer protein k-mer index
│   │   │   ├── metadata.json           # Protein database metadata
│   │   │   └── species_gene_map.json   # Species-gene ID mappings
│   │   └── resistance_db.json          # Resistance mutation positions
│   └── clean_resistance_db/             # Alternative simplified DB format
├── runtime/                             # ✅ PRODUCTION CODE
│   └── kernels/
│       ├── resistance/
│       │   ├── clean_resistance_pipeline_main.cpp  # Main clean pipeline executable
│       │   ├── enhanced_mutation_detection.cu      # GPU mutation detection with wildtype comparison
│       │   ├── fix_pipeline_resistance_loader.cpp  # Clean resistance DB loader
│       │   ├── bloom_filter.cu                    # Bloom filter pre-screening
│       │   ├── kmer_screening.cu                  # K-mer filtering
│       │   ├── translated_search_revised.cu       # 6-frame translated search
│       │   ├── hdf5_alignment_writer.cpp          # HDF5 output formatting
│       │   ├── clinical_fq_report_generator.cpp   # ✅ NEW: Clinical report generation
│       │   ├── global_fq_resistance_mapper.cpp    # ✅ FQ resistance database interface
│       │   ├── fq_resistance_positions.h          # ✅ QRDR position definitions
│       │   ├── fq_mutations_hardcoded.h           # ✅ NEW: Hardcoded FQ mutations (v0.6.4+)
│       │   ├── convert_fq_csv_to_cpp.py           # ✅ NEW: Script to convert CSV to C++ header
│       │   └── sample_csv_parser.cpp              # ✅ NEW: CSV parser for batch processing
│       └── profiler/                          # ✅ ENHANCED: GPU-accelerated microbiome profiler
│           ├── CMakeLists.txt                 # Build configuration
│           ├── minimizer_extraction.cu        # GPU minimizer extraction kernel
│           ├── minimizer_extractor.h          # Minimizer extractor interface
│           ├── minimizer_common.h             # Common data structures
│           ├── fastq_processing.h/cpp         # FASTQ I/O and pipeline (supports gzip)
│           ├── gpu_kmer_database.cu/h         # ✅ NEW: GPU k-mer database for taxonomic classification
│           ├── gpu_profiler_pipeline.cu       # ✅ NEW: Main GPU profiler with paired-end support
│           ├── minimizer_io.cpp/h             # ✅ NEW: I/O utilities for minimizers
│           ├── build_test_database.cpp        # ✅ NEW: Build test microbial database
│           ├── build_db_from_kmers.cpp        # ✅ NEW: Build database from k-mer list
│           ├── generate_test_data.py          # ✅ NEW: Generate synthetic test data
│           ├── test_minimizer.cpp             # Minimizer testing utility
│           ├── debug_minimizer.cpp            # Debug tool for minimizer extraction
│           └── hybrid_profiler_pipeline.cu    # Hybrid CPU-GPU profiler (Kraken2 integration)
├── src/                                 
│   └── python/                         
│       ├── build_integrated_resistance_db.py  # ✅ Builds integrated database
│       ├── build_clean_dynamic_database.py    # ✅ Builds both protein & nucleotide indices
│       ├── enhanced_kmer_builder.py           # ⚠️ DEPRECATED - use build_clean_dynamic_database.py
│       └── download_ncbi_20250529.py          # ✅ Downloads sequences from NCBI
├── CMakeLists.txt                       # ✅ Build system
└── build/                               
    ├── clean_resistance_pipeline        # ✅ Main executable
    └── integrated_resistance_pipeline   # Alternative full pipeline
```

---

## 🔄 Current Pipeline Workflow

### Stage 0: Data Acquisition and Database Building

#### 1. Download Reference Sequences
```bash
python src/python/download_ncbi_20250529.py \
    data/quinolone_resistance_mutation_table.csv \
    data/fq_genes \
    --email your_email@example.com \
    --max-per-gene 300
```

This downloads resistance gene sequences (gyrA, gyrB, parC, parE) from multiple species based on the Known_Quinolone_Changes.csv file.

#### 2. Build Integrated Clean Database (Protein and Nucleotide)

The database is generated by:
```bash
# Minimal command using all defaults
python3 src/python/build_clean_dynamic_database.py \
    --mutations-csv data/quinolone_resistance_mutation_table.csv
```

This builds both protein and nucleotide k-mer indices with the following defaults:
- **Protein sequences**: `data/wildtype_protein_seqs` (positional arg 1 to override)
- **Output directory**: `data/integrated_clean_db` (positional arg 2 to override)
- **Nucleotide sequences**: `data/fq_genes` (use `--nucleotide-sequences` to override)
- **Protein k-mer length**: 8 (use `--kmer-length` to override)
- **Nucleotide k-mer length**: 15 (use `--nucleotide-kmer-length` to override)
- **Reverse complement**: Enabled (use `--no-rc` to disable)

Example with all options specified:
```bash
python3 src/python/build_clean_dynamic_database.py \
    /custom/protein/dir \                    # Override protein input directory
    /custom/output/dir \                      # Override output directory
    --mutations-csv /path/to/mutations.csv \  # Required: mutations CSV file
    --nucleotide-sequences /custom/nucl/dir \ # Override nucleotide sequences directory
    --kmer-length 10 \                        # Override protein k-mer length
    --nucleotide-kmer-length 20 \             # Override nucleotide k-mer length
    --no-rc                                   # Disable reverse complement k-mers
```

The database builder (`build_clean_dynamic_database.py`):
- Reads wildtype protein sequences from `data/wildtype_protein_seqs/`
- Reads nucleotide sequences from `data/fq_genes/` (JSON format from NCBI download)
- Discovers species and genes by parsing headers
- Creates both protein and nucleotide k-mer indices in one run
- Assigns dynamic IDs to discovered species and genes
- Optionally loads known mutations from CSV

**Note**: The `enhanced_kmer_builder.py` script is no longer required as both nucleotide and protein k-mer building functionality is now integrated into `build_clean_dynamic_database.py`.

This creates:
- **Nucleotide index**: 15-mer k-mer index from downloaded sequences
  - Binary format for fast GPU searching
  - Includes reverse complement k-mers by default
  - Located in `nucleotide/` subdirectory
- **Protein database**: Translated protein sequences with 8-mer index
  - From wildtype protein sequences
  - Located in `protein/` subdirectory
- **Resistance mappings**: Gene-species-position mappings for known mutations

**Key Features of Integrated Database**:
- Species-aware gene tracking (e.g., Enterococcus_faecium_gyrA vs E.coli_gyrA)
- Position-only mutation information (no drug/resistance level data in clean version)
- Both nucleotide and protein search capabilities
- Consistent ID mappings across all components

### Stage 1: GPU Pipeline Execution

#### Running the Clean Resistance Pipeline

**New in v0.6.3**: Batch processing support and automatic output management
- Process single samples or multiple samples from a CSV file
- Sample names automatically extracted from filenames
- Default output location: `results/<sample_name>/`
- Use `--output-dir` to specify a custom output directory

**Single Sample Mode (Recommended Command)**:
```bash
# Default behavior - outputs to results/569_A_038/
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    /home/david/Documents/Code/biogpu/data/integrated_clean_db/nucleotide \
    /home/david/Documents/Code/biogpu/data/integrated_clean_db/protein \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    --no-bloom

# Custom output directory - outputs to /home/david/fq_analysis_results/569_A_038/
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    --no-bloom \
    --min-allele-depth 10 \
    --min-report-depth 20 \
    --output-dir /home/david/fq_analysis_results
```

**Batch Mode (NEW in v0.6.3)**:
```bash
# Process multiple samples from CSV file
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --no-bloom

# Validate CSV without processing (dry run to check sample names and paths)
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --dry-run

# Process with custom options and stop on first error
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --no-bloom \
    --min-allele-depth 10 \
    --output-dir /custom/output/path \
    --stop-on-error
```

**CSV File Format**:
The CSV file should have the following columns:
- **SampleName**: Name for the sample (used for output directory)
- **FilePath**: Directory containing the FASTQ files
- **R1 file**: Filename of forward reads (R1)
- **R2 file**: Filename of reverse reads (R2) - optional for single-end

Example CSV (`samples.csv`):
```csv
SampleName,FilePath,R1 file,R2 file
Sample_001,~/data/sequencing/,Sample_001_R1.fastq.gz,Sample_001_R2.fastq.gz
Sample_002,/absolute/path/to/data/,Sample_002_R1.fastq.gz,Sample_002_R2.fastq.gz
Sample_003,../relative/path/,Sample_003_R1.fastq.gz,Sample_003_R2.fastq.gz
```

Path specifications in CSV:
- Absolute paths: `/full/path/to/data/`
- Home-relative paths: `~/path/to/data/` (expands to user's home)
- Working directory relative paths: `../relative/path/` (relative to where command is run)

**Rationale**: Based on performance testing (June 8, 2025):
- `--no-bloom` provides ~6% performance improvement (212,813 vs 200,295 reads/sec)
- Bloom filter has minimal impact on results (same 1,326 protein matches)
- Smith-Waterman should remain ENABLED for sensitivity (30% more matches: 1,326 vs 922)
- This configuration balances speed and accuracy

**Pipeline Stages**:

1. **Bloom Filter Pre-screening** (Optional, default ON - **recommend OFF**)
   - Fast k-mer presence check
   - Minimal performance benefit (~6% improvement when disabled)
   - No impact on sensitivity when disabled
   - Recommendation: Use `--no-bloom` for better performance

2. **Nucleotide K-mer Matching** (Optional, default ON)
   - Binary search in sorted k-mer index
   - Species and gene-aware tracking
   - Identifies candidate resistance genes

3. **6-Frame Translated Search with banded Smith-Waterman**
   - Translates reads in all 6 frames
   - 8-mer protein k-mer seeding for efficiency
   - Banded Smith-Waterman alignment for very rapid identification of top-scoring matches
   - Enhanced mutation detection against known wildtype sequences

### Stage 2: Enhanced Mutation Detection

The new `enhanced_mutation_detection.cu` module compares observed sequences against known wildtype amino acids at resistance positions:

```cpp
// Known resistance positions are checked against wildtype
struct ResistancePosition {
    uint32_t gene_id;
    uint16_t position;      // 1-based amino acid position
    char wildtype_aa;       // Expected wildtype
    char resistant_aas[8];  // Known resistant variants
    float resistance_score;
};
```

**Key Improvements**:
- Wildtype-aware mutation calling (not dependent on reference being wildtype)
- Species-specific resistance pattern detection
- QRDR region coverage tracking
- Identity calculation based on wildtype comparison

### Stage 3: Output Generation

The pipeline generates multiple output files in the `results/<sample_name>/` directory (or custom directory if specified with `--output-dir`):

1. **HDF5 File** (`<sample_name>.h5`)
   - All alignment data
   - Nucleotide and protein results
   - Structured for downstream analysis

2. **JSON Results** (`<sample_name>.json`)
   - Human-readable resistance calls
   - Mutation details with positions
   - Pipeline statistics

3. **CSV Report** (`<sample_name>_protein_matches.csv`)
   - All protein alignments with mutation details
   - QRDR coverage information (`is_qrdr_alignment` flag)
   - Species and gene identification

4. **Clinical Reports** (NEW in v0.6.0)
   - `<sample_name>_clinical_fq_report.html` - Web-viewable clinical report
   - `<sample_name>_clinical_fq_report.json` - Machine-readable clinical data
   - `<sample_name>_clinical_fq_report.txt` - Text summary
   - Includes confidence scoring and clinical interpretation

5. **Allele Frequency Report** (NEW in v0.6.1, Enhanced in v0.6.2)
   - `<sample_name>_allele_frequencies.csv` - Comprehensive allele frequency data
   - Reports mutation frequencies at all detected positions
   - Includes wildtype and mutant amino acid counts and percentages
   - Tracks resistance mutations with position-specific depth information
   - Essential for monitoring mutation prevalence in metagenomic samples
   - **v0.6.2 Enhancement**: Configurable depth filtering (see Configuration Options)

## 🎉 Key Changes in Version 0.6.3 (June 10, 2025)

### Batch Processing Support via CSV Input
- **NEW: Batch mode** for processing multiple samples with a single command
- Added `--csv` option to specify a CSV file containing sample information
- CSV format supports flexible path specifications (absolute, relative, or home-relative)
- Added `--dry-run` option to validate CSV and preview what will be processed
- Added `--stop-on-error` option for batch error handling
- Maintains full backward compatibility with single-sample mode

### Automatic Sample Name Extraction and Output Directory Management
- **Automatic sample name extraction** from input filename (splits on `_R1.fastq.gz`)
- **Output directory is no longer a command line argument**
- Default output location: `results/<sample_name>/`
- Added `--output-dir` option for custom output directory specification
- Simplified command line interface - no need to specify output prefix
- Maintains backward compatibility with explicit output directory option

## 🎉 Key Changes in Version 0.6.2 (June 10, 2025)

### Configurable Depth Filtering for Allele Frequencies
- Added `--min-allele-depth N` parameter to control minimum depth for allele frequency analysis
  - Default: 5 (only positions with ≥5 reads are included in frequency calculations)
  - Set to 0 to include all positions regardless of coverage
- Added `--min-report-depth N` parameter to control minimum depth for reporting polymorphisms
  - Default: 0 (no filtering - all calculated frequencies are reported)
  - Set to higher values to only report high-confidence positions in CSV output
- These parameters work independently:
  - `--min-allele-depth` filters which positions are analyzed
  - `--min-report-depth` filters which analyzed positions are written to CSV
- Example: `--min-allele-depth 1 --min-report-depth 10` analyzes all positions with any coverage but only reports those with ≥10 reads

## 🎉 Key Changes in Version 0.6.1 (June 9, 2025)

### Allele Frequency Calculation
- Added comprehensive allele frequency reporting for all detected positions
- Calculates mutation frequencies with depth-aware statistics
- Tracks wildtype vs mutant amino acid distributions
- Particularly useful for monitoring E. coli gyrA D87N/G and other key resistance mutations
- Outputs detailed CSV with per-position coverage and mutation summaries

### Microbiome Profiler Development

1. **New Minimizer Extraction Module**
   - Added GPU-accelerated minimizer extraction (`runtime/kernels/profiler/minimizer_extraction.cu`)
   - Achieves ~87 Mbases/second throughput on TITAN Xp
   - K-mer size: 31, Window size: 15
   - Extracts ~17 minimizers per 150bp read

2. **FASTQ Processing Pipeline**
   - Modularized FASTQ reading and processing (`fastq_processing.h/cpp`)
   - Supports gzipped FASTQ files
   - Multi-threaded GPU processing with configurable batch sizes
   - Namespace organization: `biogpu::FastqReader`, `biogpu::GPUMinimizerPipeline`

3. **Build System Updates**
   - Separated profiler components into `runtime/kernels/profiler/`
   - Fixed CUDA device linking issues with static libraries
   - Added debug tools: `fastq_pipeline_debug`, `debug_minimizer`

4. **Performance Optimization Attempts**
   - Attempted pinned memory optimizations (reverted due to issues)
   - Current implementation uses standard CUDA memory transfers
   - Processing 1M reads takes ~1 second
   - File reading: 77-84ms for 100k reads
   - GPU processing: 10-16ms for 10k sequences

5. **Known Issues and Solutions**
   - Fixed namespace confusion: `ReadBatch` and `Minimizer` are NOT in biogpu namespace
   - Fixed CUDA linking issues by enabling `CUDA_RESOLVE_DEVICE_SYMBOLS`
   - Statistics collection hanging issue still under investigation
   - Low minimizer count (0.71 per read) in main pipeline needs debugging

6. **Files to Keep/Remove**
   - KEEP: `minimizer_extraction.cu` (working version)
   - REMOVE: `minimizer_extraction_optimized_attempt2.cu`
   - REMOVE: `minimizer_extraction_optimized.cu.backup`

## 🎉 Key Changes in Version 0.6.0 (June 8, 2025)

### Major Fixes and Improvements

1. **Fixed QRDR Detection**
   - `is_qrdr_alignment` now properly set on host side after protein matches
   - Checks if alignments or mutations fall within QRDR regions
   - Uses FQ resistance database for accurate QRDR position identification

2. **Fixed Species/Gene ID Mapping**
   - Resolved issue where all genes showed as "unknown"
   - Fixed metadata loading that was blocked by disabled fq_mutation_reporter
   - Now correctly identifies species (e.g., Escherichia_coli) and genes (e.g., gyrA)

3. **New Clinical Report Generator**
   - Comprehensive clinical reports in HTML, JSON, and text formats
   - Confidence scoring (0-95%) for resistance detection
   - Clinical interpretation with actionable recommendations
   - Distinguishes between known FQ resistance and QRDR mutations

4. **Performance Metrics**
   - Added reads/second calculation to pipeline output
   - Typical performance: 16,667-21,739 reads/second on NVIDIA TITAN Xp

### Example Clinical Report Output
```
CLINICAL INTERPRETATION
----------------------
HIGH CONFIDENCE: Fluoroquinolone resistance detected
Confidence: 95%

SPECIES BREAKDOWN
-----------------
Escherichia_coli: RESISTANT
  FQ resistance mutations: 658
  QRDR mutations: 1726
  Genes affected: gyrA, parE

KNOWN FQ RESISTANCE MUTATIONS:
  Escherichia_coli gyrA D87G - High-level FQ resistance
```

## 🚧 Key Changes in Version 0.5.0

### Database Structure Changes

1. **Integrated Database Design**
   - Single source of truth for nucleotide and protein data
   - Consistent species-gene ID mappings
   - Clean separation of sequence data and resistance annotations

2. **Species-Gene Aware Tracking**
   - Each gene is tracked with its species (e.g., gene_id maps to "Enterococcus_faecium_gyrA")
   - Allows detection of species-specific resistance patterns
   - Prevents cross-species false positives

3. **Enhanced Protein Database**
   - Contains both mutant and wildtype sequences
   - 5-mer k-mer indexing for fast seeding
   - Metadata includes accession mappings

### Algorithm Improvements

1. **Two-Stage Alignment**
   - Fast k-mer seeding (5-mers for proteins)
   - Precise Smith-Waterman for candidate regions
   - Banded alignment for efficiency

2. **Wildtype-Aware Mutation Detection**
   - Compares against known wildtype amino acids
   - Not dependent on reference sequence being wildtype
   - Handles both known and novel mutations

3. **Multi-Level Filtering**
   - Optional Bloom filter (memory efficient)
   - Optional nucleotide k-mer matching
   - Mandatory protein search for resistance detection

## 📋 Configuration Options

### Pipeline Flags (v0.6.0+)
- `--no-bloom`: Disable Bloom filter pre-screening (RECOMMENDED for 6% speed improvement)
- `--no-sw`: Disable Smith-Waterman alignment (NOT recommended - reduces sensitivity by 30%)
- `--min-allele-depth N`: Minimum read depth for allele frequency analysis (default: 5) (v0.6.2+)
- `--min-report-depth N`: Minimum read depth for reporting in CSV output (default: 0) (v0.6.2+)
- `--output-dir PATH`: Custom output directory (default: results/<sample_name>) (v0.6.3+)
- `--csv FILE`: Process multiple samples from CSV file (batch mode) (v0.6.3+)
- `--dry-run`: Validate CSV and show what would be processed without running (v0.6.3+)
- `--stop-on-error`: Stop batch processing on first error (batch mode only) (v0.6.3+)
- `--fq-csv FILE`: Path to fluoroquinolone resistance mutations CSV (optional, uses hardcoded data if not provided) (v0.6.4+)

### Fluoroquinolone Resistance Mutations Database (v0.6.4+)

The pipeline includes a comprehensive database of known fluoroquinolone resistance mutations. By default, it uses 255 hardcoded mutations across 33 species without requiring any external files.

**Default Behavior**: The pipeline automatically uses built-in resistance mutations data:
```bash
# No --fq-csv needed, uses hardcoded mutations
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz
```

**Custom Mutations Database**: You can override the hardcoded mutations with your own CSV file:
```bash
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --fq-csv /path/to/custom_mutations.csv
```

**Adding New Mutations**: To extend the database with additional mutations:

1. **Start with the existing database**:
   ```bash
   cp data/quinolone_resistance_mutation_table.csv my_mutations.csv
   ```

2. **Add your mutations** to the CSV file:
   ```csv
   species,gene,location,wt,mut
   Klebsiella_pneumoniae,gyrA,91,T,I
   Mycobacterium_tuberculosis,gyrA,94,G,S
   ```

3. **Convert to C++ header** using the helper script:
   ```bash
   python3 runtime/kernels/resistance/convert_fq_csv_to_cpp.py my_mutations.csv > \
       runtime/kernels/resistance/fq_mutations_hardcoded.h
   ```

4. **Rebuild the pipeline**:
   ```bash
   cd runtime/kernels/resistance/build
   make clean_resistance_pipeline
   ```

**CSV Format Requirements**:
- **species**: Species name with underscores (e.g., `Escherichia_coli`)
- **gene**: Gene name (e.g., `gyrA`, `parC`, `gyrB`, `parE`)
- **location**: Amino acid position (1-based)
- **wt**: Wild-type amino acid (single letter)
- **mut**: Mutant amino acid that confers resistance (single letter)

**Why Add Custom Mutations?**:
- Include newly discovered resistance mutations from recent literature
- Add species-specific mutations not in the default database
- Include mutations specific to your local epidemiology
- Support novel or emerging resistance patterns

The `convert_fq_csv_to_cpp.py` script:
- Validates CSV format and amino acid codes
- Generates proper C++ header with include guards
- Provides summary statistics by species and gene
- Skips invalid entries (with warnings)

### Performance Parameters
- Batch size: 10,000 reads (configurable)
- Bloom filter threshold: 3 k-mers minimum
- Smith-Waterman identity threshold: 90%
- K-mer sizes: 15 (nucleotide), 5 (protein)

## 🎯 Usage Examples

### Standard Analysis
```bash
# Build database (one time)
make build_integrated_resistance_db

# Run pipeline (outputs to results/sample/)
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz

# Run pipeline with custom output directory
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --output-dir /custom/output/path
```

### High-Sensitivity Mode (Disable Pre-filtering, Report all alignments and alleles)
```bash
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --disable-bloom-filter \
    --disable-kmer-match \
    --min-allele-depth 1 \
    --min-report-depth 1 
```

### Validation with Synthetic Data
```bash
# Generate synthetic reads with known mutations
python src/python/generate_synthetic_reads.py \
    --num-reads 1000000 \
    --mutation-rate 0.01 \
    --output-prefix 1M_synthetic_reads

# Run pipeline (outputs to results/1M_synthetic_reads/)
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    1M_synthetic_reads_R1.fastq.gz \
    1M_synthetic_reads_R2.fastq.gz
```

### Batch Processing Examples (v0.6.3+)
```bash
# Create a CSV file listing your samples
cat > samples.csv << EOF
SampleName,FilePath,R1 file,R2 file
Patient_001,~/sequencing/batch1/,Patient_001_R1.fastq.gz,Patient_001_R2.fastq.gz
Patient_002,~/sequencing/batch1/,Patient_002_R1.fastq.gz,Patient_002_R2.fastq.gz
Patient_003,~/sequencing/batch2/,Patient_003_R1.fastq.gz,Patient_003_R2.fastq.gz
EOF

# Validate the CSV and check file paths
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --dry-run

# Process all samples with optimized settings
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --no-bloom \
    --min-allele-depth 10 \
    --output-dir /results/batch_analysis

# Generate the CSV automatically from a directory
python runtime/kernels/resistance/generate_sample_csv.py \
    ~/sequencing/batch1/ \
    -o batch1_samples.csv \
    --recursive
```

### Depth Filtering Examples (v0.6.2+)
```bash
# Include all positions in analysis, report only high-confidence (≥20 reads)
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --min-allele-depth 1 \
    --min-report-depth 20

# Standard analysis (≥5 reads) but report everything
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --min-allele-depth 5 \
    --min-report-depth 0

# Ultra-sensitive: analyze and report all positions
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --min-allele-depth 0 \
    --min-report-depth 0
```

## 🔧 Building from Source

```bash
# Clone repository
git clone [repo] && cd biogpu

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build clean pipeline
make clean_resistance_pipeline

# Build database
make build_integrated_resistance_db
```

## 📊 Output Interpretation

### Resistance Calls in JSON
```json
{
  "read_id": 12345,
  "gene_id": 1,
  "species_id": 0,
  "frame": 1,
  "alignment_score": 85.2,
  "identity": 0.96,
  "mutations": [
    {
      "position": 83,
      "ref_aa": "S",
      "query_aa": "L"
    }
  ],
  "peptide": "MKSRISTGDMLR..."
}
```

### Key Metrics
- **Identity**: Percentage match to wildtype sequence
- **Alignment Score**: Smith-Waterman score (higher = better match)
- **Frame**: Translation frame (1-3 forward, -1 to -3 reverse)
- **Position**: 1-based amino acid position in gene

## 🚀 Performance Characteristics

- **Throughput**: 
  - ~16,667-21,739 reads/second on NVIDIA TITAN Xp (v0.6.0)
  - ~50,000-100,000 reads/second on NVIDIA A100 (estimated)
- **Memory Usage**: ~4GB GPU memory for standard run
- **Accuracy**: >99% for known resistance mutations at 10x coverage
- **Sensitivity**: Detects mutations at 10% allele frequency

## 🤝 Version History

### v0.7.0 (June 11 2025) - GPU MICROBIOME PROFILER WITH PAIRED-END SUPPORT
- **NEW: GPU-accelerated microbiome profiler** for taxonomic classification
- **Added paired-end read support** to GPU profiler pipeline
- **GPU k-mer database** with hash table implementation
- **Test database builder** for rapid development and testing
- **Binary database format** for fast loading and GPU transfer
- Minimizer extraction achieving ~210 Mbases/second
- Classification performance: ~500k reads/sec (single-end), ~250k pairs/sec (paired-end)
- Fixed database loading issue (separated upload from rebuild)
- Supports gzipped FASTQ files (single and paired)
- Command: `./gpu_profiler_pipeline database.db reads_R1.fastq.gz [reads_R2.fastq.gz]`

### v0.6.3 (June 10 2025) - BATCH PROCESSING & AUTOMATIC SAMPLE NAMING
- **NEW: Batch processing mode** - process multiple samples from CSV file
- **Added `--csv` option** for specifying sample manifest
- **Added `--dry-run` option** to validate CSV without processing
- **Added `--stop-on-error` option** for batch error handling
- **Automatic sample name extraction** from input filenames
- **Simplified command line** - removed output prefix argument
- **Default output structure**: `results/<sample_name>/`
- **Added `--output-dir` option** for custom output locations
- Improved user experience with sensible defaults
- Includes `generate_sample_csv.py` utility to create CSV from directory

### v0.6.2 (June 10 2025) - CONFIGURABLE DEPTH FILTERING
- **Added configurable minimum depth parameters** for allele frequency analysis
- `--min-allele-depth`: Controls which positions are included in analysis
- `--min-report-depth`: Controls which positions are written to output
- Allows flexible filtering based on coverage requirements
- Default behavior unchanged (min 5 reads for analysis, no filtering for output)

### v0.6.0 (June 8 2025) - FULLY FUNCTIONAL FQ RESISTANCE DETECTION
- **Fixed QRDR detection** - now properly identifies QRDR alignments
- **Fixed species/gene mapping** - no more "unknown" identifications
- **Added clinical report generation** with confidence scoring
- **Added performance metrics** (reads/second)
- Successfully detecting FQ resistance mutations (e.g., E. coli gyrA D87G)
- Clinical interpretation with 95% confidence for high-level resistance

### v0.5.1 (June 7 2025 - Evening Update)
- **Implemented banded Smith-Waterman alignment** for improved sensitivity
- **Expanded protein matches** - now detecting ~2.3M matches per 1M reads
- Successfully detecting mutations relative to reference sequences
- Integrated global FQ resistance database from `backup_scripts/tools/Quinolone_resistance_mutation_table.csv`
- **KNOWN ISSUE**: FQ resistance database mapping not working properly yet

### v0.5.0 (June 7 2025)
- Integrated clean resistance pipeline
- Species-aware mutation detection
- Enhanced protein search with wildtype comparison
- Optional filtering stages

### v0.4.0 (June 2025)
- Added Bloom filter pre-screening
- 6-frame translated search
- Smith-Waterman alignment

### v0.3.0
- Basic k-mer matching pipeline
- Simple mutation detection

---

## 🐛 Known Issues and TODOs

### TODO for Next Version
1. **Apply bloom_and_sw_flags.patch**
   - Add command-line flags to enable/disable Bloom filter
   - Add command-line flags to enable/disable Smith-Waterman alignment
   - This will allow users to control pipeline stages for performance tuning
   - See `docs/bloom_sw_configuration_guide.md` for implementation details

### Current Known Issues
1. **Bloom filter showing 0% pass rate**
   - All reads are failing Bloom filter check
   - Pipeline still works because protein search is independent
   - Needs investigation - may be related to k-mer size or filter parameters

2. **No nucleotide matches reported**
   - Possibly related to Bloom filter issue
   - Nucleotide k-mer matching stage may be bypassed

3. **CSV parsing error**
   - Error parsing qnrB19 entries with "NA" position
   - Non-critical - other mutations still processed correctly

### Optimization Opportunities
- Performance optimization for large-scale datasets
- Memory usage optimization for protein database
- Parallel processing of paired-end reads

---

## 🧬 GPU-Accelerated Microbiome Profiler (NEW in v0.7.0)

### Overview

The BioGPU microbiome profiler is a high-performance, GPU-accelerated tool for taxonomic classification of metagenomic reads. It uses minimizer-based k-mer matching for rapid species identification.

**Key Features**:
- GPU-accelerated minimizer extraction (~210 Mbases/second)
- GPU hash table lookups for taxonomic classification
- Support for both single-end and paired-end reads
- Gzipped FASTQ file support
- Efficient memory usage with streaming processing

### Building the Profiler Components

```bash
# Navigate to profiler directory
cd runtime/kernels/profiler
mkdir build && cd build
cmake ..
make

# Available targets:
# - gpu_profiler_pipeline    : Main GPU profiler (with paired-end support)
# - build_test_database      : Build test microbial database
# - build_db_from_kmers      : Build database from k-mer list
# - test_minimizer           : Test minimizer extraction
# - hybrid_profiler_pipeline : Hybrid CPU-GPU profiler (requires Kraken2 DB)
```

### Creating a K-mer Database

#### Option 1: Build Test Database (for development/testing)
```bash
# Create a small test database with 9 common species
./build_test_database test_microbiome.db

# This creates:
# - test_microbiome.db (4MB) with 15,789 k-mers from 9 species
# - test_microbiome.db.test.fastq with synthetic test reads
```

**Test Database Specifications**:
- Contains k-mers from 9 common bacterial species:
  - Escherichia coli (taxon ID: 100)
  - Klebsiella pneumoniae (taxon ID: 101)
  - Enterococcus faecalis (taxon ID: 102)
  - Staphylococcus aureus (taxon ID: 103)
  - Streptococcus pneumoniae (taxon ID: 104)
  - Haemophilus influenzae (taxon ID: 105)
  - Proteus mirabilis (taxon ID: 106)
  - Pseudomonas aeruginosa (taxon ID: 107)
  - Clostridioides difficile (taxon ID: 108)
- ~1,750 unique k-mers per species
- Uses MurmurHash3 with proper finalizer for k-mer hashing
- Applies canonical k-mer selection (min of forward/reverse)
- Generates matching test reads with known ground truth
- Perfect for verifying GPU pipeline functionality

#### Option 2: Build from K-mer List (for production)
```bash
# Format: one k-mer per line with taxon ID
# kmer<tab>taxon_id
./build_db_from_kmers kmers.txt microbiome.db

# K-mer file format example:
# ATCGATCGATCGATCGATCGATCGATCGATC    562
# GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT    1280
```

#### Option 3: Build from Reference Genomes (production workflow)
```bash
# 1. Download reference genomes from NCBI
python download_genomes.py --species-list species.txt --output-dir genomes/

# 2. Extract k-mers with taxonomic labels
python extract_kmers.py genomes/ --k 31 --output kmers.txt

# 3. Build GPU database
./build_db_from_kmers kmers.txt microbiome_production.db
```

### Important Database Building Considerations

**K-mer Selection**:
- K=31 is standard for species-level classification
- Larger k provides more specificity but less sensitivity
- Smaller k increases sensitivity but may reduce specificity

**Database Size Considerations**:
- Each k-mer entry requires 16 bytes (8 bytes hash + 4 bytes taxon ID + 4 bytes padding)
- Hash table size is automatically set to next power of 2 above 2x the number of k-mers
- Example: 1M k-mers → 2M buckets → ~32MB GPU memory

**Taxonomic ID Requirements**:
- Must use NCBI taxonomy IDs for compatibility
- Species-level IDs recommended (strain-level can cause over-classification)
- Database builder validates IDs against NCBI taxonomy

**Minimizer Parameters**:
- Window size (w=15) is hardcoded to match database k-mer size (k=31)
- This ensures compatible minimizer selection between database and query
- Changing k requires recompiling with matching window size

### Running the GPU Profiler - Complete Examples

#### Quick Test with Synthetic Data
```bash
# 1. Build test database and generate test reads
./build_test_database test_microbiome.db

# 2. Run profiler on test data - should get 100% classification
./gpu_profiler_pipeline test_microbiome.db test_microbiome.db.test.fastq

# Expected output:
# Running in single-end mode
# Loaded database from test_microbiome.db
#   Table size: 32768 buckets
#   Total k-mers: 15789
# Database uploaded to GPU (4 MB)
# Batch timing (μs): extract=1537 transfer=252 lookup=96 classify=46 return=27 total=1958 minimizers=19204
# 
# === Final Results (Single-End) ===
# Total reads: 1400
# Classified reads: 1400 (100%)
# Processing time: 0.003 seconds
# Reads per second: 466667
```

#### Real Data - Single-End Mode
```bash
# Profile a single FASTQ file
./gpu_profiler_pipeline microbiome.db sample.fastq.gz

# Profile with timing details (redirect stderr to see batch timings)
./gpu_profiler_pipeline microbiome.db sample.fastq.gz 2>&1 | tee profile_output.txt

# Extract only summary statistics
./gpu_profiler_pipeline microbiome.db sample.fastq.gz | grep -A 5 "Final Results"
```

#### Real Data - Paired-End Mode
```bash
# Profile paired-end reads (R1 and R2)
./gpu_profiler_pipeline microbiome.db sample_R1.fastq.gz sample_R2.fastq.gz

# Real example with E. coli data
./gpu_profiler_pipeline test_microbiome.db \
    /home/david/Documents/Code/biogpu/data/569_A_038_R1.fastq.gz \
    /home/david/Documents/Code/biogpu/data/569_A_038_R2.fastq.gz

# Example output:
# Running in paired-end mode
# Loaded database from test_microbiome.db
#   Table size: 32768 buckets  
#   Total k-mers: 15789
# Database uploaded to GPU (4 MB)
# Paired batch timing (μs): extract=15439 transfer=693 lookup=100 classify=75 return=86 total=16393 minimizers=284490
# Processed 100000 read pairs, 0 classified (0%)
# ...
# === Final Results (Paired-End) ===
# Total read pairs: 2802881
# Classified pairs: 0 (0%)  
# Processing time: 31.449 seconds
# Pairs per second: 89135.8
```

#### Performance Testing
```bash
# Test minimizer extraction only
./test_minimizer sample.fastq.gz

# Expected output:
# Reading sequences from sample.fastq.gz...
# Read 10000 sequences
# Testing minimizer extraction (k=31, m=15)...
# Results:
# - Processed 10000 sequences
# - Total minimizers: 142178
# - Average minimizers per sequence: 14.2178
# - Processing time: 16 ms
# - Throughput: 625000 sequences/second

# Debug mode with detailed timing
./debug_minimizer sample.fastq.gz | head -20

# Benchmark with large file
time ./gpu_profiler_pipeline large_db.bin large_sample.fastq.gz | tail -10
```

#### Batch Processing Script
```bash
#!/bin/bash
# process_samples.sh - Process multiple samples

DB="microbiome_production.db"
OUTPUT_DIR="profiling_results"
mkdir -p $OUTPUT_DIR

for sample in samples/*.fastq.gz; do
    basename=$(basename $sample .fastq.gz)
    echo "Processing $basename..."
    
    ./gpu_profiler_pipeline $DB $sample > $OUTPUT_DIR/${basename}_profile.txt
    
    # Extract summary
    grep -A 5 "Final Results" $OUTPUT_DIR/${basename}_profile.txt > $OUTPUT_DIR/${basename}_summary.txt
done

# Combine all summaries
cat $OUTPUT_DIR/*_summary.txt > $OUTPUT_DIR/all_samples_summary.txt
```

#### Integration with Downstream Analysis
```bash
# 1. Profile to identify abundant species
./gpu_profiler_pipeline microbiome.db sample_R1.fastq.gz sample_R2.fastq.gz \
    > sample_profile.txt

# 2. Parse results to get species counts (example parser)
grep "Classified" sample_profile.txt | awk '{print $2, $4}' > species_counts.txt

# 3. Filter reads by species (requires custom script)
python filter_reads_by_taxon.py \
    --profile sample_profile.txt \
    --taxon 100 \  # E. coli taxon ID
    --input sample_R1.fastq.gz sample_R2.fastq.gz \
    --output ecoli_R1.fastq.gz ecoli_R2.fastq.gz

# 4. Run targeted analysis on filtered reads
./clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    ecoli_R1.fastq.gz ecoli_R2.fastq.gz
```

#### Troubleshooting Examples
```bash
# Check database integrity
./gpu_profiler_pipeline test_microbiome.db test_microbiome.db.test.fastq 2>&1 | \
    grep -E "(Total k-mers|Table size|Classified)"

# Profile a subset of reads
zcat sample.fastq.gz | head -40000 | \
    ./gpu_profiler_pipeline microbiome.db /dev/stdin

# Compare single vs paired-end classification
echo "Single-end results:"
./gpu_profiler_pipeline db.bin sample_R1.fastq.gz | grep "Classified reads"

echo "Paired-end results:"  
./gpu_profiler_pipeline db.bin sample_R1.fastq.gz sample_R2.fastq.gz | grep "Classified pairs"

# Monitor GPU usage during profiling
nvidia-smi dmon -i 0 -s u -d 1 &
./gpu_profiler_pipeline large_db.bin large_sample.fastq.gz
kill %1
```

### Performance Characteristics

**Minimizer Extraction**:
- K-mer size: 31
- Window size: 15
- Throughput: ~210 Mbases/second on NVIDIA TITAN Xp
- ~13-17 minimizers extracted per 150bp read

**Classification Performance**:
- Single-end: ~500,000 reads/second
- Paired-end: ~250,000 pairs/second
- GPU memory usage: ~100MB + database size
- Batch size: 10,000 reads

**Database Specifications**:
- Hash table with power-of-2 buckets
- Linear probing for collision resolution
- ~20% collision rate at 50% load factor
- Binary format for fast loading

### Technical Implementation

**GPU Kernels**:
1. **Minimizer Extraction** (`minimizer_extraction.cu`)
   - Parallel k-mer generation
   - Window-based minimizer selection
   - MurmurHash3 for k-mer hashing

2. **K-mer Lookup** (`gpu_kmer_database.cu`)
   - GPU hash table with linear probing
   - Batch lookups for all minimizers
   - Taxon ID retrieval

3. **Read Classification** (`classify_reads_kernel`)
   - Majority vote classification
   - Confidence scoring based on k-mer hits
   - Species-level assignment

**Paired-End Processing**:
- Minimizers extracted from both reads independently
- Combined for classification (improves accuracy)
- Single taxonomic assignment per read pair
- Synchronized batch processing

### Database Format

The binary database format includes:
```cpp
struct HashTableParams {
    uint64_t table_size;    // Number of buckets (power of 2)
    uint64_t bucket_size;   // Entries per bucket (default: 4)
    uint64_t total_entries; // Total k-mers in database
};

struct KmerEntry {
    uint64_t hash;          // MurmurHash3 of k-mer
    uint32_t taxon_id;      // NCBI taxonomic ID
    uint32_t reserved;      // Padding for alignment
};
```

### Integration with Existing Pipeline

The profiler can be used as a pre-screening step before resistance detection:
```bash
# 1. Profile the sample to identify species
./gpu_profiler_pipeline microbiome.db sample_R1.fastq.gz sample_R2.fastq.gz > profile.txt

# 2. Extract species of interest (e.g., E. coli reads)
python filter_by_species.py sample_R1.fastq.gz sample_R2.fastq.gz \
    --species "Escherichia coli" --profile profile.txt \
    --output ecoli_R1.fastq.gz ecoli_R2.fastq.gz

# 3. Run resistance detection on filtered reads
./clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    ecoli_R1.fastq.gz ecoli_R2.fastq.gz
```

### Future Enhancements

1. **Abundance Estimation**
   - Read counting per taxon
   - Coverage-based abundance calculation
   - Relative abundance reporting

2. **Hierarchical Classification**
   - Genus/family level classification
   - LCA (Lowest Common Ancestor) algorithm
   - Taxonomic tree integration

3. **Output Formats**
   - Kraken-style reports
   - BIOM format export
   - Interactive visualization

4. **Database Tools**
   - Reference genome downloader
   - K-mer database merger
   - Database update utilities

---

*This is a living document. Update with each significant change.*
