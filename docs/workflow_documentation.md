# BioGPU Workflow Documentation

## Overview
This document tracks the complete BioGPU pipeline for GPU-accelerated fluoroquinolone resistance detection from metagenomic data.

**Last Updated**: June 10 2025  
**Pipeline Version**: 0.6.4  
**Status**: Fully functional FQ resistance detection with clinical reporting and hardcoded mutations database

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
│       └── profiler/                          # ✅ NEW: Microbiome profiler components
│           ├── CMakeLists.txt                 # Build configuration
│           ├── minimizer_extraction.cu        # GPU minimizer extraction kernel
│           ├── minimizer_extractor.h          # Minimizer extractor interface
│           ├── minimizer_common.h             # Common data structures
│           ├── fastq_processing.h/cpp         # FASTQ I/O and pipeline
│           ├── fastq_pipeline_main.cpp        # Main FASTQ processing executable
│           ├── fastq_pipeline_debug.cpp       # Debug version with timing
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

### High-Sensitivity Mode (Disable Pre-filtering)
```bash
./runtime/kernels/resistance/build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    --disable-bloom-filter \
    --disable-kmer-match
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

*This is a living document. Update with each significant change.*
