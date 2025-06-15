# BioGPU Workflow Documentation

## Overview
This document tracks the complete BioGPU pipeline for GPU-accelerated fluoroquinolone resistance detection from metagenomic data and taxonomic classification.

**Last Updated**: June 15 2025  
**Pipeline Version**: 0.11.0  
**Status**: Fully functional FQ resistance detection with clinical reporting, GPU-accelerated Kraken2-style taxonomic classifier with streaming support, paired-end processing, and batch CSV processing

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
│           ├── hybrid_profiler_pipeline.cu    # Hybrid CPU-GPU profiler (Kraken2 integration)
│           ├── hybrid_gpu_cpu_profiler.cu     # ✅ NEW: Hybrid GPU/CPU profiler main implementation
│           ├── hybrid_gpu_cpu_profiler.py     # ✅ NEW: Python analysis toolkit for hybrid profiler
│           ├── create_hybrid_gpu_cpu_db.cpp   # ✅ NEW: Hybrid database builder
│           ├── hybrid_gpu_cpu_build_script.sh # ✅ NEW: Build script for hybrid components
│           ├── hierarchical_profiler_pipeline.cu   # ✅ NEW: Hierarchical GPU database profiler
│           ├── build_hierarchical_db.cpp      # ✅ NEW: Build hierarchical database
│           └── minimized/                     # ✅ NEW: Adaptive profiler components
│               ├── adaptive_paired_end_profiler.cu  # Adaptive paired-end profiler main
│               ├── build_minimizer_db.cpp     # Build minimizer-based database
│               ├── gpu_community_profiler.cu  # GPU community profiling kernels
│               ├── streaming_gpu_profiler.cu  # Streaming mode implementation
│               └── analyze_database.cu        # Database analysis utilities
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

### v0.10.0 (June 14 2025) - ADAPTIVE PAIRED-END PROFILER WITH AUTOMATIC OPTIMIZATION
- **NEW: Adaptive paired-end profiler** that automatically selects optimal processing strategy
- **Automatic database analysis** determines if direct GPU or streaming mode is needed
- **No manual configuration required** - adapts to available GPU memory
- **Minimizer-based approach** using Kraken2-inspired algorithm (k=35, m=31)
- **Incremental database building** supports adding new genomes without full rebuild
- **Confidence scoring** based on coverage breadth and unique minimizers
- **Paired-end optimization** improves classification accuracy
- Performance: ~400k reads/sec (direct mode), ~200k-300k reads/sec (streaming mode)
- Executables: `adaptive_paired_end_profiler`, `build_minimizer_db`

### v0.9.0 (June 13 2025) - HYBRID GPU/CPU PROFILER WITH COMPREHENSIVE ANALYSIS
- **NEW: Hybrid GPU/CPU profiler** for memory-efficient large database processing
- **Two-stage processing**: CPU screening followed by GPU precision matching
- **Memory-mapped database** for instant loading and efficient access
- **Comprehensive output formats**: 
  - Human-readable organism reports with confidence scores
  - Detailed abundance tables with coverage metrics
  - Taxonomic summaries at all hierarchical levels
  - Kraken-compatible output format
- **Python analysis toolkit** with diversity metrics and visualizations
- **Support for strain-level resolution** (e.g., 190 strains → 19 species)
- Performance: ~300k reads/second with minimal GPU memory usage
- Executable: `hybrid_gpu_cpu_profiler`

### v0.8.0 (June 12 2025) - HIERARCHICAL GPU DATABASE FOR LARGE-SCALE PROFILING
- **NEW: Hierarchical GPU database** for databases exceeding GPU memory
- **Dynamic tier loading** with LRU cache management
- **Memory-configurable profiling** with --memory parameter
- **Tested with 143M k-mers** (1.2GB database, 81M unique k-mers)
- Supports databases 10x larger than available GPU memory
- **Critical fixes**: 
  - Implemented missing `load_database()` function for proper tier metadata loading
  - Fixed tier file loading to skip 64-byte headers
  - Added canonical k-mer support to match minimizer extraction
- New executables: hierarchical_profiler_pipeline, build_hierarchical_db
- Transparent API - same usage as standard profiler

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

## 🚀 Hierarchical GPU Database for Large-Scale Pathogen Profiling (NEW in v0.8.0)

### Overview

The hierarchical GPU database implementation addresses the challenge of profiling large pathogen databases that exceed GPU memory limits. Instead of loading the entire database into GPU memory, it uses a tiered approach with dynamic loading and caching.

**Key Features**:
- **Memory-efficient**: Works with databases larger than GPU memory
- **Dynamic tier loading**: Loads database portions on-demand
- **LRU cache management**: Keeps frequently accessed tiers in GPU memory
- **Transparent to users**: Same API as standard database

### Building with Hierarchical Support

```bash
cd runtime/kernels/profiler
mkdir build && cd build
cmake ..
make

# New hierarchical targets:
# - hierarchical_profiler_pipeline : GPU profiler with hierarchical database
# - build_hierarchical_db         : Build hierarchical database from k-mers
# - test_hierarchical_db          : Test hierarchical functionality
```

### Creating a Hierarchical Database

```bash
# Build hierarchical database from k-mer list
./build_hierarchical_db database_kmers.txt pathogen_hierarchical_db [options]

# Options:
#   --tier-size <MB>      Target size per tier in MB (default: 512)
#   --sort-hash           Sort by hash instead of frequency
#   --no-manifest         Don't create manifest file
#   --progress <N>        Progress update interval (default: 1000000)

# Example:
./build_hierarchical_db data/pathogen_db/database_kmers.txt \
    data/pathogen_hierarchical.db \
    --tier-size 512 \
    --progress 10000000

# The builder will:
# 1. Read k-mers and compute canonical forms (min of forward/reverse complement)
# 2. Apply MurmurHash3 finalizer for k-mer hashing
# 3. Sort k-mers by frequency (most frequent in tier 0)
# 4. Partition into memory-efficient tiers (default: 512MB each)
# 5. Create manifest.json with tier metadata
# 6. Output structure:
#    pathogen_hierarchical.db/
#    ├── manifest.json     # Database metadata
#    ├── tier_0.bin       # Most frequent k-mers
#    ├── tier_1.bin       # Next tier
#    └── tier_N.bin       # Least frequent k-mers
```

### Important: Canonical K-mer Requirement

The hierarchical database MUST use canonical k-mers (minimum of forward and reverse complement) to match the minimizer extraction pipeline. The database builder automatically:

1. Computes both forward and reverse complement for each k-mer
2. Selects the lexicographically smaller one (canonical)
3. Applies MurmurHash3 to the canonical k-mer

This ensures that k-mers extracted from reads will match those in the database regardless of strand orientation.

### Using the Hierarchical Profiler

```bash
# Run with hierarchical database (same command syntax)
./hierarchical_profiler_pipeline pathogen_hierarchical.db \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    --memory 8  # Limit GPU memory to 8GB

# Additional options:
# --memory <GB>      : Maximum GPU memory to use (default: 10)
# --optimize-cache   : Preload frequently accessed tiers
# --output-prefix    : Output file prefix
```

### How It Works

1. **Database Structure**:
   - Database split into fixed-size tiers (default: 512MB)
   - Each tier contains k-mers within a hash range
   - Metadata tracks tier boundaries for fast lookup

2. **Query Processing**:
   - Queries grouped by tier based on hash value
   - Required tiers loaded on-demand
   - LRU eviction when memory limit reached
   - Batch processing minimizes tier switching

3. **Performance Optimization**:
   - Frequently accessed tiers kept in cache
   - Async loading overlaps with computation
   - Hash-based partitioning ensures good locality

### Example: Processing Large Pathogen Database

```bash
# 1. Download comprehensive pathogen database (e.g., from NCBI)
python src/python/download_microbial_genomes.py \
    data/pathogen_db \
    --email your_email@example.com \
    --genomes-per-species 10

# 2. Process genomes to extract k-mers
python src/python/process_existing_genomes.py \
    data/pathogen_db \
    --k 31 --stride 5

# 3. Build hierarchical database (handles 143M k-mers → 33GB)
./build_hierarchical_db \
    data/pathogen_db/database_kmers.txt \
    data/pathogen_hierarchical.db

# 4. Run profiler with memory limit
./hierarchical_profiler_pipeline \
    data/pathogen_hierarchical.db \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    --memory 10 \
    --output-prefix results/569_A_038_hierarchical

# Output includes:
# - results/569_A_038_hierarchical_abundance.tsv
# - results/569_A_038_hierarchical_detailed_report.txt
# - results/569_A_038_hierarchical_performance.json
```

### Performance Characteristics

**Memory Usage**:
- Configurable GPU memory limit (default: 10GB)
- Each tier uses ~512MB when loaded
- Metadata overhead: <100MB

**Performance Impact**:
- ~10-20% slower than fully-loaded database
- Cache hit rate typically >90% after warmup
- Tier loading: ~100ms per 512MB tier

**Scalability**:
- Tested with 143M k-mers (33GB database)
- Can handle databases 10x larger than GPU memory
- Linear scaling with database size

### Configuration Options

```cpp
struct HierarchicalDBConfig {
    size_t max_gpu_memory_gb = 10;    // Maximum GPU memory to use
    size_t tier_size_mb = 512;        // Size of each tier
    size_t cache_tiers = 3;           // Number of tiers to keep cached
    bool use_lru_eviction = true;     // Use LRU for tier eviction
    bool preload_frequent_tiers = true; // Preload most accessed tiers
};
```

### Monitoring and Debugging

The hierarchical profiler provides detailed statistics:

```
=== Hierarchical Database Statistics ===
Total lookups: 10000000
Cache hit rate: 92.5%
Tier loads: 45
Tier evictions: 12
Average lookup time: 2.3 μs
GPU memory usage: 8192.0 MB
Active tiers: 16/64
```

### When to Use Hierarchical vs Standard Database

**Use Hierarchical Database when**:
- Database size exceeds GPU memory
- Processing diverse samples with varying species
- Memory constraints on shared GPU systems
- Need to scale to larger databases over time

**Use Standard Database when**:
- Database fits comfortably in GPU memory
- Processing focused on specific species
- Maximum performance is critical
- Consistent species distribution across samples

### Troubleshooting Hierarchical Database Issues

**Problem: 0 classified reads despite database loading**

Common causes and solutions:

1. **Missing canonical k-mer support** (fixed in v0.8.0):
   - Ensure database was built with `build_hierarchical_db` that uses canonical k-mers
   - Rebuild database if it was created before the canonical k-mer fix

2. **Database loading issues**:
   - Check that `manifest.json` exists and is readable
   - Verify tier files are present and have correct permissions
   - Ensure tier file headers are being skipped (64-byte header)

3. **Hash mismatch**:
   - Verify both database builder and query pipeline use MurmurHash3 finalizer
   - Check that k-mer encoding matches (2-bit encoding, same base mapping)

4. **Memory limits**:
   - Increase `--memory` parameter if tiers are being evicted too frequently
   - Monitor cache hit rate in output statistics

**Debugging commands**:
```bash
# Test with small database first
./test_hierarchical_db data/test_hierarchical.db

# Check database integrity
cat data/pathogen_hierarchical.db/manifest.json | jq .

# Monitor tier loading
./hierarchical_profiler_pipeline db_path reads.fq 2>&1 | grep "Loading tier"
```

### Future Improvements

1. **Multi-level Hierarchy**:
   - Two-level tier system for better cache efficiency
   - Species-aware tier organization

2. **Predictive Loading**:
   - Anticipate tier access patterns
   - Background prefetching

3. **Compressed Tiers**:
   - On-GPU decompression
   - Further memory savings

4. **Distributed Processing**:
   - Multi-GPU tier distribution
   - Network-attached tier storage

---

## 🧬 Minimizer-Based Adaptive Paired-End Profiler (NEW in v0.10.0)

### Overview

The minimizer-based adaptive paired-end profiler provides an intelligent, memory-aware approach to metagenomic profiling. It automatically analyzes database characteristics and selects the optimal processing strategy—either direct GPU loading for smaller databases or streaming for larger ones.

**Key Features**:
- **Automatic strategy selection**: Analyzes database size and automatically chooses between direct GPU or streaming mode
- **Paired-end support**: Native support for paired-end reads with improved classification accuracy
- **Minimizer-based matching**: Uses Kraken2-inspired minimizer approach (k=35, m=31 by default)
- **Memory-aware processing**: Adapts to available GPU memory without manual configuration
- **Confidence scoring**: Provides confidence scores based on coverage and unique minimizers

### Building the Minimizer Database

#### Build Minimizer Database Tool

The `build_minimizer_db` tool creates a minimizer-based database from FASTA files:

```bash
# Build a new database from genome directory
./build_minimizer_db <input_fasta_directory> <output_database> [options]

# Options:
#   --incremental      Update existing database instead of rebuilding
#   --k-size N         K-mer size (default: 35)
#   --minimizer-size N Minimizer size (default: 31)

# Example: Build database from pathogen genomes
./build_minimizer_db data/pathogen_genomes minimized_pathogens_db

# Example: Incrementally add new genomes to existing database
./build_minimizer_db data/new_genomes minimized_pathogens_db --incremental
```

**Database Building Process**:
1. Scans directory for FASTA files (.fasta, .fa, .fna extensions)
2. Extracts organism names from filenames or headers
3. Generates taxonomy IDs using consistent hashing
4. Extracts minimizers using sliding window approach
5. Calculates uniqueness scores for each minimizer
6. Creates binary database with efficient lookup structure

**Output**:
- Binary database file (e.g., `minimized_pathogens_db`)
- Summary file with statistics (`minimized_pathogens_db.summary`)

### Running the Adaptive Paired-End Profiler

```bash
# Basic usage with paired-end reads
./adaptive_paired_end_profiler <database> <R1.fastq> <R2.fastq> [output_prefix]

# Single-end or interleaved reads
./adaptive_paired_end_profiler <database> <reads.fastq> [output_prefix]

# Examples:
./adaptive_paired_end_profiler minimized_pathogens_db \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_results

# With custom output prefix
./adaptive_paired_end_profiler minimized_pathogens_db \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    minimized_results
```

### How Adaptive Mode Works

The profiler analyzes the database before processing:

1. **Database Analysis Phase**:
   - Calculates estimated GPU memory requirements
   - Counts organisms and total minimizer entries
   - Determines if streaming is needed

2. **Strategy Selection**:
   - **Direct GPU Mode** (⚡): Selected when database fits in GPU memory
     - Entire database loaded to GPU
     - Maximum performance
     - Best for databases <30GB
   
   - **Streaming Mode** (🔄): Selected for larger databases
     - Database processed in chunks
     - Memory-efficient
     - Slightly slower but handles any size

3. **Processing**:
   - Extracts minimizers from reads
   - Matches against database using selected strategy
   - Calculates organism abundances and confidence scores

### Output Files

The profiler generates the following outputs:

1. **Abundance Table** (`<prefix>_adaptive_abundance.tsv`):
   ```
   taxonomy_id  organism_name  taxonomy_path  relative_abundance  coverage_breadth  unique_minimizers  total_hits  paired_hits  confidence_score  processing_mode
   ```

2. **Console Output**:
   - Database analysis results
   - Processing statistics
   - Top organisms detected
   - Performance metrics

### Complete Workflow Example

```bash
# Step 1: Download reference genomes (using existing tools)
python src/python/download_microbial_genomes.py \
    data/pathogen_genomes \
    --email your@email.com \
    --genomes-per-species 10

# Step 2: Build minimizer database
cd runtime/kernels/profiler
./build_minimizer_db ../../data/pathogen_genomes pathogen_minimized_db

# Step 3: Check database statistics
cat pathogen_minimized_db.summary

# Step 4: Profile clinical samples
./adaptive_paired_end_profiler pathogen_minimized_db \
    ../../data/clinical_R1.fastq.gz \
    ../../data/clinical_R2.fastq.gz \
    clinical_profile

# Step 5: Analyze results
head -20 clinical_profile_adaptive_abundance.tsv
```

### Batch Processing Support (NEW in v0.10.1)

The adaptive paired-end profiler now supports batch processing of multiple samples, maintaining consistency with the resistance detection pipeline's CSV format.

#### Batch Processing Modes

```bash
# Usage
./adaptive_paired_profiler <database> <mode> [options]

# Modes:
# - single: Process a single sample (backwards compatible)
# - batch: Process multiple samples from a list file

# Examples:
# Single sample mode
./adaptive_paired_profiler microbes.db single reads_R1.fq.gz reads_R2.fq.gz

# Batch mode with CSV file (recommended)
./adaptive_paired_profiler microbes.db batch samples.csv batch_results/

# Batch mode with text file (legacy format)
./adaptive_paired_profiler microbes.db batch samples.txt batch_results/
```

#### CSV Format (Consistent with Resistance Pipeline)

The profiler supports the same CSV format used by the resistance detection pipeline:

```csv
SampleName,FilePath,R1 file,R2 file
Sample1,~/data/sequencing/,sample1_R1.fq.gz,sample1_R2.fq.gz
Sample2,/absolute/path/,sample2_R1.fastq.gz,sample2_R2.fastq.gz
Sample3,../relative/path/,sample3_reads.fastq
```

**CSV Features**:
- **Flexible column naming**: Recognizes variants like "Sample Name", "SampleName", "sample_name"
- **Path handling**: Supports absolute, relative, and home-relative paths (~)
- **Auto-detection**: Automatically detects delimiter (comma, tab, or semicolon)
- **Single-end support**: Leave R2 file column empty for single-end reads
- **File validation**: Checks file existence before processing

#### Legacy Text Format

For backwards compatibility, a simple tab-delimited format is also supported:

```
# Text file format (tab-separated)
sample_name    R1_file    R2_file
sample1    /path/to/sample1_R1.fq.gz    /path/to/sample1_R2.fq.gz
sample2    /path/to/sample2_R1.fq.gz    /path/to/sample2_R2.fq.gz
```

#### Batch Processing Workflow

```bash
# 1. Create a CSV file listing your samples
cat > samples.csv << EOF
SampleName,FilePath,R1 file,R2 file
Patient_001,~/sequencing/batch1/,Patient_001_R1.fastq.gz,Patient_001_R2.fastq.gz
Patient_002,~/sequencing/batch1/,Patient_002_R1.fastq.gz,Patient_002_R2.fastq.gz
Patient_003,~/sequencing/batch2/,Patient_003_R1.fastq.gz,Patient_003_R2.fastq.gz
Control_001,~/sequencing/controls/,Control_001_R1.fastq.gz,Control_001_R2.fastq.gz
EOF

# 2. Run batch profiling
./adaptive_paired_profiler pathogen_minimized_db batch samples.csv profiling_results/

# 3. View batch summary
cat profiling_results/batch_summary.tsv
```

#### Output Structure

Batch processing creates an organized output structure:

```
profiling_results/
├── batch_summary.tsv              # Summary of all samples
├── Patient_001/
│   └── Patient_001_abundance.tsv  # Individual sample results
├── Patient_002/
│   └── Patient_002_abundance.tsv
├── Patient_003/
│   └── Patient_003_abundance.tsv
└── Control_001/
    └── Control_001_abundance.tsv
```

**Batch Summary Format**:
```
sample_name    total_organisms_detected    top_organism    top_abundance
Patient_001    15                         E. coli         0.4523
Patient_002    23                         K. pneumoniae   0.3892
```

#### Key Advantages

1. **Database loaded once**: Significant performance improvement for multiple samples
2. **Consistent format**: Same CSV format as resistance pipeline enables sequential analysis
3. **Progress tracking**: Real-time updates on sample processing
4. **Error handling**: Continues processing if individual samples fail
5. **Organized output**: Each sample gets its own directory

#### Performance Characteristics

- **Batch efficiency**: ~50% faster than processing samples individually
- **Memory usage**: Database kept in GPU memory throughout batch
- **Typical throughput**: 
  - First sample: Includes database loading time
  - Subsequent samples: 200-500k reads/second depending on mode

#### Integration with Resistance Pipeline

Process the same samples through both pipelines:

```bash
# 1. Profile microbiome composition
./adaptive_paired_profiler microbes.db batch samples.csv profiling_results/

# 2. Detect resistance in same samples
./clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    --csv samples.csv \
    --output-dir resistance_results/

# 3. Combine results for comprehensive analysis
python combine_results.py \
    --profiling profiling_results/ \
    --resistance resistance_results/ \
    --output combined_analysis/
```

### Performance Characteristics

**Database Building**:
- Processing speed: ~100 MB/s of FASTA input
- Memory usage: Proportional to unique minimizers
- Incremental updates supported

**Adaptive Profiling**:
- Direct GPU mode: ~400k-500k reads/second
- Streaming mode: ~200k-300k reads/second
- Automatic mode selection optimizes for available hardware
- Paired-end processing adds ~20% overhead vs single-end

### Advantages of Adaptive Approach

1. **No Manual Configuration**: Automatically selects optimal strategy
2. **Memory Efficient**: Handles databases of any size
3. **Performance Optimized**: Uses direct GPU when possible
4. **Future-Proof**: Adapts to different GPU memory sizes
5. **Simplified Workflow**: Single tool for all database sizes

### When to Use This Profiler

**Best for**:
- General metagenomic profiling
- Mixed database sizes
- Paired-end sequencing data
- Users who want automatic optimization
- Situations where GPU memory varies

**Consider alternatives when**:
- You need hierarchical taxonomic assignment
- Database is extremely large (>100GB)
- You require Kraken2 compatibility
- Single-species focused analysis

## 🧬 Monolithic vs Hybrid GPU/CPU Database Profiling Workflows

### Overview of Database Approaches

BioGPU now supports multiple database architectures for metagenomic profiling:

1. **Monolithic GPU Database**: Entire database loaded into GPU memory
   - Best for databases that fit in GPU memory (<10GB)
   - Maximum performance with all data on GPU
   - Single file binary format

2. **Hybrid GPU/CPU Database**: Database stored in CPU memory, streamed to GPU
   - Handles databases larger than GPU memory (tested up to 700MB)
   - CPU memory mapped file for efficient access
   - GPU processes batches of queries
   - Two-stage processing: CPU screening followed by GPU precision matching
   - Comprehensive taxonomic tracking from Kingdom to Strain level

### Complete Workflow Commands

#### Step 1: Download Sequences from NCBI

Both database types use the same source data:

```bash
# Download microbial genomes from NCBI
python src/python/download_microbial_genomes.py \
    data/pathogen_profiler_db \
    --email dbhaslam@gmail.com \
    --genomes-per-species 10 \
    --k 31 \
    --stride 5

# This creates:
# - data/pathogen_profiler_db/ with downloaded genomes
# - Organized by species subdirectories
# - K-mer extraction parameters saved
```

#### Step 2A: Build Monolithic GPU Database

```bash
# Extract k-mers from downloaded genomes
cd data/pathogen_profiler_db
ls */*.fasta | while read genome; do
    # Extract k-mers (example - adjust based on your extraction tool)
    extract_kmers.py $genome --k 31 --stride 5 >> all_kmers.txt
done

# Build monolithic database
cd ../../runtime/kernels/profiler
./build_db_from_kmers data/pathogen_profiler_db/all_kmers.txt \
    monolithic_pathogen.db

# Expected format for all_kmers.txt:
# ACGTACGTACGTACGTACGTACGTACGTACG<tab>562
# (31-mer sequence<tab>taxonomy_id)
```

#### Step 2B: Build Hybrid GPU/CPU Database

```bash
# Build hybrid database directly from FASTA directory
cd runtime/kernels/profiler
./create_hybrid_gpu_cpu_db \
    data/pathogen_profiler_db \
    hybrid_pathogen.db

# This creates:
# - hybrid_pathogen.db (binary database file)
# - Includes all sequences and taxonomy mappings
# - Memory-mapped format for CPU access
```

#### Step 3A: Profile Reads with Monolithic Database

```bash
# Single-end reads
./gpu_profiler_pipeline monolithic_pathogen.db \
    sample.fastq.gz \
    --output-prefix monolithic_results

# Paired-end reads
./gpu_profiler_pipeline monolithic_pathogen.db \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    --output-prefix monolithic_results

# Output files:
# - monolithic_results_summary.txt
# - monolithic_results_abundance.tsv
# - monolithic_results_kraken_style.txt
```

#### Step 3B: Profile Reads with Hybrid Database

```bash
# Single-end reads
./hybrid_gpu_cpu_profiler hybrid_pathogen.db \
    sample.fastq.gz \
    hybrid_results

# Paired-end reads (not yet supported in hybrid mode)
# Use monolithic pipeline for paired-end analysis

# Output files:
# - hybrid_results_organism_report.txt     # Human-readable comprehensive report
# - hybrid_results_abundance_table.tsv     # Detailed abundance with confidence scores
# - hybrid_results_taxonomy_summary.tsv    # Taxonomic composition at all levels
# - hybrid_results_coverage_stats.tsv      # Coverage breadth and depth metrics
# - hybrid_results_kraken_style.txt        # Kraken-compatible output format
```

### Hybrid Profiler Configuration

The hybrid profiler uses these default parameters:

```cpp
struct HybridConfig {
    int kmer_size = 31;              // K-mer length for matching
    int kmer_step = 50;              // Sample every N bp from reads
    float min_abundance = 1e-9;      // Minimum abundance threshold for reporting
    int max_gpu_organisms = 500;     // Maximum organisms loaded to GPU at once
    float coverage_threshold = 0.0001; // Minimum coverage for organism reporting
};
```

These can be adjusted by modifying the source code and recompiling.

### Performance Comparison

**Monolithic GPU Database**:
- Load time: ~1 second for 100MB database
- Classification: ~500k reads/second (single-end)
- Memory: Entire database in GPU memory
- Best for: Small to medium databases (<10GB)

**Hybrid GPU/CPU Database**:
- Load time: Instant (memory mapped)
- Classification: ~300k reads/second (single-end)
- Memory: Only active batches in GPU memory
- Best for: Large databases or limited GPU memory

### Database Format Details

**Monolithic Database Binary Format**:
```cpp
struct Header {
    uint64_t table_size;      // Hash table buckets
    uint64_t total_entries;   // Total k-mers
};
struct Entry {
    uint64_t hash;            // K-mer hash
    uint32_t taxon_id;        // Taxonomy ID
    uint32_t reserved;        // Padding
};
```

**Hybrid Database Binary Format**:
```cpp
struct Header {
    char magic[8];            // "HYBRIDDB"
    uint32_t version;         // Format version
    uint32_t num_organisms;   // Total organisms
    uint64_t total_size;      // Database size
};
struct Organism {
    uint32_t tax_id;          // Taxonomy ID
    uint32_t name_offset;     // Name string offset
    uint32_t seq_offset;      // Sequence offset
    uint32_t seq_length;      // Sequence length
    // Additional metadata...
};
```

### Choosing Between Monolithic and Hybrid

**Use Monolithic When**:
- Database fits in GPU memory
- Maximum performance needed
- Processing many samples with same database
- Need paired-end support

**Use Hybrid When**:
- Database exceeds GPU memory
- Running on shared GPU systems
- Need quick startup (no load time)
- Processing diverse databases

### Example: Complete E. coli Profiling

```bash
# 1. Download E. coli and related genomes
python src/python/download_microbial_genomes.py \
    data/ecoli_db \
    --email your@email.com \
    --species "Escherichia coli,Klebsiella pneumoniae,Salmonella enterica" \
    --genomes-per-species 20 \
    --k 31 --stride 5

# 2. Build both database types
# Monolithic
./build_db_from_kmers data/ecoli_db/kmers.txt ecoli_mono.db

# Hybrid  
./create_hybrid_gpu_cpu_db data/ecoli_db ecoli_hybrid.db

# 3. Profile clinical sample
# Monolithic (with paired-end)
./gpu_profiler_pipeline ecoli_mono.db \
    clinical_R1.fastq.gz clinical_R2.fastq.gz \
    --output-prefix mono_clinical

# Hybrid (single-end only currently)
./hybrid_gpu_cpu_profiler ecoli_hybrid.db \
    clinical_R1.fastq.gz \
    hybrid_clinical

# 4. Compare results
diff mono_clinical_abundance.tsv hybrid_clinical_abundance_table.tsv
```

### Output Interpretation for Hybrid Profiler

The hybrid profiler generates comprehensive output files with detailed metrics:

1. **Organism Report** (`*_organism_report.txt`):
   ```
   Rank: 1
   Organism: Escherichia coli
   Taxonomy: Bacteria > Proteobacteria > Gammaproteobacteria > Enterobacterales > Enterobacteriaceae > Escherichia > Escherichia coli
   Relative Abundance: 20.69%
   Coverage Breadth: 88.47%
   Coverage Depth: 5.63x
   Confidence Score: 95.79%
   ```

2. **Abundance Table** (`*_abundance_table.tsv`):
   - Tab-separated format for easy parsing
   - Includes organism ID, name, full taxonomy path
   - Metrics: relative abundance, coverage breadth/depth, unique k-mers, confidence score
   - Confidence based on coverage breadth and unique k-mer count

3. **Taxonomy Summary** (`*_taxonomy_summary.tsv`):
   - Aggregated abundances at each taxonomic level
   - Useful for understanding community structure
   - Format: level, taxon, abundance, sample_count

4. **Coverage Statistics** (`*_coverage_stats.tsv`):
   - Detailed per-organism coverage metrics
   - Includes genome size estimates
   - Useful for assessing detection reliability

### Python Analysis Toolkit

The hybrid profiler includes a comprehensive Python analysis script (`hybrid_gpu_cpu_profiler.py`) that provides:

```bash
# Basic usage
python hybrid_gpu_cpu_profiler.py \
    --database hybrid_pathogen.db \
    --fastq sample.fastq.gz \
    --output analysis_results

# With all options
python hybrid_gpu_cpu_profiler.py \
    --database hybrid_pathogen.db \
    --fastq sample_R1.fastq.gz \
    --output analysis_results \
    --min-abundance 0.001 \
    --top-n 50 \
    --export-formats phyloseq,qiime2,biom \
    --visualize \
    --html-report
```

**Features**:
- **Diversity Metrics**: Shannon, Simpson, Chao1, and more
- **Interactive Visualizations**: 
  - Abundance bar plots
  - Rarefaction curves
  - PCoA plots (for multiple samples)
- **Export Formats**: 
  - phyloseq R objects
  - QIIME2 artifacts
  - BIOM format
- **HTML Dashboard**: Interactive report with all results

### Two-Stage Processing Algorithm

The hybrid profiler uses an innovative two-stage approach:

**Stage 1: CPU K-mer Screening**
- Rapid k-mer matching using memory-mapped database
- Identifies candidate organisms based on k-mer hits
- Filters out organisms with insufficient evidence

**Stage 2: GPU Precision Matching**
- Top candidates loaded to GPU for detailed analysis
- Precise abundance calculation
- Coverage breadth and depth computation
- Confidence scoring based on multiple metrics

This approach balances speed and memory efficiency while maintaining accuracy.

### Building the Hybrid Profiler

```bash
# Quick build using provided script
cd runtime/kernels/profiler
./hybrid_gpu_cpu_build_script.sh

# Manual build with custom options
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make hybrid_gpu_cpu_profiler

# Build all profiler tools
make all
# Creates:
# - hybrid_gpu_cpu_profiler
# - create_hybrid_gpu_cpu_db
# - test_hybrid_functionality
```

### Database Creation Details

The hybrid database builder supports various input formats:

```bash
# From directory of FASTA files
./create_hybrid_gpu_cpu_db \
    /path/to/fasta/directory \
    output_database.db

# The builder:
# - Automatically detects taxonomy from headers
# - Supports multiple strains per species
# - Creates canonical k-mer index
# - Calculates uniqueness scores
# - Generates memory-mapped binary format
```

**Database Structure**:
- Header with magic number and version
- Organism metadata (taxonomy, genome stats)
- K-mer index with hash values
- Sequential storage for optimal memory access

### Troubleshooting

**Common Issues**:

1. **Out of GPU memory with monolithic**:
   - Switch to hybrid database
   - Reduce database size
   - Use hierarchical database (see v0.8.0)

2. **Slow performance with hybrid**:
   - Increase batch size in code
   - Ensure database file is on fast storage (SSD)
   - Check CPU-GPU transfer bottlenecks
   - Adjust `kmer_step` parameter for sparser sampling

3. **Different results between approaches**:
   - Verify same k-mer parameters (k=31)
   - Check taxonomy ID mappings
   - Ensure canonical k-mer handling matches
   - Compare confidence scores, not just abundances

4. **Low confidence scores**:
   - May indicate insufficient sequencing depth
   - Check coverage breadth metrics
   - Consider adjusting min_abundance threshold

## 🚀 GPU-Accelerated Kraken2-Style Taxonomic Classifier (NEW - December 15, 2025)

### Overview

Added a GPU-accelerated implementation of the Kraken2 algorithm for ultra-fast taxonomic classification. This implementation builds Kraken2-style databases from reference genomes and performs GPU-accelerated classification of metagenomic reads.

**Key Features**:
- GPU-accelerated minimizer extraction and database building
- Kraken2-compatible algorithm with spaced seeds (k=35, m=31, spaces=7)
- Memory-aware batch processing to handle large genome collections
- LCA (Lowest Common Ancestor) computation for accurate taxonomic assignment
- Configurable batch sizes to prevent out-of-memory errors
- Support for NCBI taxonomy integration
- **NEW**: Support for gzip-compressed FASTQ input files
- **NEW**: Paired-end read classification with concordance scoring
- **NEW**: Integrated taxonomy loading from database file

### Building the Kraken2-Style Database

The GPU Kraken pipeline builds a database in two stages:

```bash
# Navigate to the build directory
cd /home/david/Documents/Code/biogpu
cd build/runtime/kernels/profiler/k2like_from_scratch

# Build the database from a directory of genome files
./gpu_kraken_pipeline build \
    --genome-dir /path/to/genomes \
    --output ./kraken_db \
    --taxonomy /path/to/taxonomy \
    --gpu-batch 1  # Process 1 genome at a time to avoid OOM

# Example with actual paths:
./gpu_kraken_pipeline build \
    --genome-dir data/pathogen_profiler_db/genomes \
    --output ./test_kraken_db \
    --taxonomy data/ \
    --gpu-batch 1
```

**Database Building Options**:
- `--k <int>`: k-mer length (default: 35)
- `--minimizer-len <int>`: Minimizer length (default: 31)
- `--spaces <int>`: Spaced seed spacing (default: 7)
- `--gpu-batch <int>`: Number of genomes to process at once (default: 1000)
- `--taxonomy <dir>`: NCBI taxonomy directory containing nodes.dmp and names.dmp

### Memory Management Improvements

To handle out-of-memory errors when processing large genome collections, we implemented several improvements:

1. **Reduced Default Batch Sizes**:
   - `MAX_SEQUENCE_BATCH`: Reduced from 10,000 to 1,000 sequences
   - `MAX_MINIMIZERS_PER_BATCH`: Reduced from 5M to 500K minimizers

2. **Memory Pre-allocation Check**:
   ```cpp
   // Check available GPU memory before allocation
   size_t free_memory, total_memory_gpu;
   cudaMemGetInfo(&free_memory, &total_memory_gpu);
   
   // Only allocate if enough memory available
   if (total_required > free_memory * 0.9) {
       std::cerr << "ERROR: Not enough GPU memory!" << std::endl;
       return false;
   }
   ```

3. **Dynamic Batch Size Control**:
   - Added `--gpu-batch` parameter to control processing batch size
   - Allows processing as few as 1 genome at a time for memory-constrained situations
   - Memory allocation scales with batch size

4. **Optimized GPU Kernel**:
   - Changed from one-thread-per-genome to grid-stride loops
   - Better parallelism for processing large genome sequences
   - Prevents kernel timeouts on large genomes

### Example Workflow

```bash
# 1. Download reference genomes (if needed)
python src/python/download_microbial_genomes.py \
    data/kraken_genomes \
    --email your@email.com \
    --genomes-per-species 10

# 2. Download NCBI taxonomy (if needed)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C data/

# 3. Build the Kraken database with memory-safe settings
./gpu_kraken_pipeline build \
    --genome-dir data/kraken_genomes \
    --output ./my_kraken_db \
    --taxonomy data/ \
    --gpu-batch 1 \
    --k 35 \
    --minimizer-len 31

# 4. Classify reads
./gpu_kraken_pipeline classify \
    --database ./my_kraken_db \
    --reads sample.fastq.gz \
    --output results.txt \
    --confidence 0.1 \
    --batch-size 10000

# 5. Complete pipeline (build + classify)
./gpu_kraken_pipeline pipeline \
    --genome-dir data/kraken_genomes \
    --reads sample.fastq.gz \
    --output ./results \
    --taxonomy data/ \
    --gpu-batch 1
```

### Performance Characteristics

**Database Building**:
- Successfully processes 190 genomes (725MB) in 12 seconds
- Extracts ~5.1M minimizers at 56 MB/second
- Memory usage: ~11MB GPU memory for batch size 1
- Scales linearly with genome count

**Example Output**:
```
=== DATABASE BUILD COMPLETE ===
Total build time: 12 seconds

=== BUILD STATISTICS ===
Total sequences processed: 190
Total bases processed: 725519735
Valid minimizers extracted: 5112927
Sequence processing time: 6.99s
Minimizer extraction time: 5.94s
Processing rate: 5.61e+07 bases/second
```

### Troubleshooting Out-of-Memory Errors

If you encounter "out of memory" errors:

1. **Reduce GPU batch size**:
   ```bash
   # Process one genome at a time
   --gpu-batch 1
   ```

2. **Check available GPU memory**:
   ```bash
   nvidia-smi
   ```

3. **Monitor memory during execution**:
   ```bash
   watch -n 1 nvidia-smi
   ```

4. **Clear GPU memory if needed**:
   ```bash
   # Reset GPU
   sudo nvidia-smi --gpu-reset
   ```

### Implementation Details

**GPU Database Builder** (`gpu_kraken_database_builder.cu`):
- Processes genomes in configurable batches
- Extracts minimizers using GPU kernels
- Computes LCA assignments for taxonomic classification
- Saves database in binary format for fast loading

**Pipeline Main** (`kraken_pipeline_main.cu`):
- Unified executable for building and classification
- Command-line interface similar to Kraken2
- Support for single-end and paired-end reads
- Generates classification reports and statistics

**Key Data Structures**:
```cpp
struct GPUGenomeInfo {
    uint32_t taxon_id;
    uint32_t sequence_offset;
    uint32_t sequence_length;
    uint32_t genome_id;
};

struct GPUMinimizerHit {
    uint64_t minimizer_hash;
    uint32_t taxon_id;
    uint32_t position;
    uint32_t genome_id;
};
```

### New Classification Features (December 2025)

#### Gzip-Compressed Input Support

The classifier now supports reading gzip-compressed FASTQ files directly:

```bash
# Single-end classification with gzipped input
./gpu_kraken_pipeline classify \
    --database full_kraken_db \
    --reads data/569_A_038_R1.fastq.gz \
    --output test_classify

# Paired-end classification with gzipped inputs
./gpu_kraken_pipeline classify \
    --database full_kraken_db \
    --reads data/569_A_038_R1.fastq.gz \
    --reads2 data/569_A_038_R2.fastq.gz \
    --output test_classify_paired
```

The implementation automatically detects `.gz` file extensions and uses zlib for decompression, maintaining compatibility with uncompressed FASTQ files.

#### Enhanced Taxonomy Loading

Added `load_taxonomy_tree_from_db()` function that reads taxonomy information directly from the minimizer database:

```cpp
bool PairedEndGPUKrakenClassifier::load_taxonomy_tree_from_db(const std::string& db_file) {
    // Reads organism metadata from database header
    // Extracts taxonomy IDs, names, and hierarchical relationships
    // Builds taxonomy tree structure on GPU for classification
}
```

This eliminates the need for separate taxonomy files during classification, as all necessary taxonomic information is embedded in the database.

## GPU-Accelerated Kraken2-Style Taxonomic Classifier

### Overview

The `gpu_kraken_pipeline` is a complete GPU-accelerated implementation of Kraken2-style taxonomic classification with enhanced features:

- **Streaming processing** for large datasets
- **Paired-end support** with concordance scoring
- **Batch processing** from CSV files
- **Compressed input** support (.gz files)
- **Memory-efficient** design with automatic scaling

### Key Commands

```bash
# Build database from genomes
./gpu_kraken_pipeline build --genome-dir ./genomes --output kraken_db

# Single-end classification
./gpu_kraken_pipeline classify --database kraken_db --reads sample.fastq.gz --output results.txt

# Paired-end classification
./gpu_kraken_pipeline classify --database kraken_db --reads R1.fastq.gz --reads2 R2.fastq.gz --output results.txt

# Batch processing from CSV
./gpu_kraken_pipeline batch --csv samples.csv --database kraken_db --batch-output batch_results
```

#### Classification Output Format

The classifier produces Kraken2-compatible output with the following columns:
1. **Classification status**: `C` (classified) or `U` (unclassified)
2. **Read ID**: Sequential read identifier
3. **Taxon ID**: NCBI taxonomy ID of assigned taxon (0 for unclassified)
4. **Read length**: Length of the read in base pairs
5. **Confidence score**: Classification confidence (0.0-1.0)
6. **K-mer matches**: Number of k-mers matching the assigned taxon
7. **Total k-mers**: Total k-mers in the read

Example output:
```
C	read_1	680797	151	0.026	3	117
C	read_2	986813	149	0.053	6	114
U	read_3	0	151	0.000	0	117
```

#### Performance Metrics

With the new optimizations:
- **Classification speed**: ~787,402 reads/second on NVIDIA TITAN Xp
- **Database loading**: ~11 seconds for 2.6M taxonomy nodes
- **Memory usage**: Efficient batch processing with 10,000 reads per batch
- **Gzip overhead**: Minimal impact on performance

### New Features (v0.11.0)

#### 1. Streaming Input for Large Datasets

The pipeline now supports streaming processing to handle large FASTQ files without loading all reads into memory:

```cpp
struct StreamingConfig {
    size_t read_batch_size = 100000;          // Process 100K read pairs per batch
    size_t max_gpu_memory_for_reads_mb = 8192; // Reserve 8GB for read processing
    bool enable_streaming = true;
    bool show_batch_progress = true;
};
```

**Features:**
- Automatic streaming for files > 1GB
- Configurable batch size (default: 100,000 reads)
- Memory-efficient processing
- Real-time progress reporting

**Usage:**
```bash
# Force streaming mode with custom batch size
./gpu_kraken_pipeline classify \
    --database kraken_db \
    --reads large_R1.fastq.gz \
    --reads2 large_R2.fastq.gz \
    --output results.txt \
    --streaming-batch 50000 \
    --force-streaming
```

#### 2. Paired-End Support with Concordance

Full support for paired-end reads with concordance scoring:

```bash
# Paired-end classification
./gpu_kraken_pipeline classify \
    --database kraken_db \
    --reads sample_R1.fastq.gz \
    --reads2 sample_R2.fastq.gz \
    --output paired_results.txt \
    --confidence 0.1
```

**Output format for paired-end:**
```
C	pair_1	680797	151|149	0.667	2|3	117|115
U	pair_2	0	150|151	0.000	0|0	116|117
```

#### 3. Batch Processing from CSV

Process multiple samples efficiently with a single database load:

```bash
# Batch processing
./gpu_kraken_pipeline batch \
    --csv samples.csv \
    --database kraken_db \
    --batch-output batch_results \
    --streaming-batch 100000
```

**CSV Format:**
```csv
Sample Name,File Path,R1 file,R2 file
Sample001,/data/fastq/,sample001_R1.fastq.gz,sample001_R2.fastq.gz
Sample002,/data/fastq/,sample002_R1.fastq.gz,sample002_R2.fastq.gz
Sample003,/data/fastq/,sample003_R1.fastq.gz,
```

**Features:**
- Database loaded once for all samples
- Automatic sample directory creation
- Path validation with error reporting
- Batch summary report generation
- Support for mixed single/paired-end samples

**Options:**
- `--no-sample-dirs`: Output all results in one directory
- `--stop-on-error`: Stop batch on first error
- `--no-validate-paths`: Skip file existence checks

#### 4. Performance Optimizations

- **Streaming classification**: ~843,032 reads/second
- **Batch processing**: Database kept in GPU memory across samples
- **Compressed input**: Native gzip support with minimal overhead
- **Memory management**: Automatic memory scaling based on GPU capacity

### TODO: Future Enhancements

1. **HTML Report Generation**: 
   - Parse classification output to generate comprehensive HTML reports
   - Include taxonomic distribution charts
   - Sample comparison matrices for batch runs
   - Integration with clinical reporting system

2. **Output Parsing and Analysis**:
   - Aggregate statistics across batches
   - Species abundance calculations
   - Diversity metrics computation
   - Contamination detection

3. **Enhanced Batch Processing**:
   - Parallel sample processing on multiple GPUs
   - Real-time monitoring dashboard
   - Integration with LIMS systems

2. **Summary Statistics Output**:
   - Generate TSV file with aggregated classification results
   - Include species abundance calculations
   - Provide read count summaries by taxonomic level
   - Create Kraken2-style report format

3. **Paired-End Enhancements**:
   - Implement concordance-based classification for paired reads
   - Add pair-aware confidence scoring
   - Support for interleaved FASTQ files

4. **Performance Optimizations**:
   - Multi-GPU support for large-scale processing
   - Streaming mode for ultra-large databases
   - Memory-mapped database loading

5. **Integration Features**:
   - Bracken-style abundance estimation
   - Export to standard microbiome analysis formats
   - Integration with downstream analysis pipelines

### Example Complete Workflow

```bash
# 1. Build Kraken2-style database
./gpu_kraken_pipeline build \
    --genome-dir ./refseq_genomes \
    --output ./full_kraken_db \
    --taxonomy ./taxonomy \
    --gpu-batch 1

# 2. Classify single-end reads
./gpu_kraken_pipeline classify \
    --database ./full_kraken_db \
    --reads sample.fastq.gz \
    --output sample_classification.txt

# 3. Analyze results (future enhancement)
# ./summarize_kraken_output sample_classification.txt > sample_summary.tsv
```

### Future Enhancements

1. **Streaming Mode**: Process genomes without loading entire sequences
2. **Multi-GPU Support**: Distribute processing across multiple GPUs
3. **Incremental Updates**: Add new genomes without rebuilding entire database
4. **Compression**: Reduce database size with GPU-friendly compression
5. **Bracken Integration**: Abundance estimation at species level

*This is a living document. Update with each significant change.*
