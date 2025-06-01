# BioGPU Workflow Documentation

## Overview
This document tracks the complete BioGPU pipeline for GPU-accelerated fluoroquinolone resistance detection from metagenomic data.

**Last Updated**: May 32, 2025  
**Pipeline Version**: 0.2.0  
**Status**: Working prototype with k-mer screening + simplified alignment

---

## 📁 Project Structure

```
biogpu/
├── data/
│   ├── Known_Quinolone_Changes.csv      # Source: (https://www.ncbi.nlm.nih.gov/pathogens/microbigge//#Quinolone%20AND%20Point)
│   ├── Known_Efflux_Pump_Genes.csv      # built using /src/R/parse_quinolone_mutations.R
│   ├── genomes/                         # Downloaded from NCBI (bacteria/, fungi/, plasmids/, viral/)
│   │   ├── bacteria/                    # Bacterial reference genomes
│   │   ├── fungi/                       # Fungal genomes
│   │   ├── plasmids/                    # Plasmid sequences
│   │   └── viral/                       # Viral genomes
│   ├── gpu_resistance_db/               # ❌ Created but NOT USED (mutations.pkl, genes.pkl, etc.)
│   ├── fq_resistance_index/             # ✅ ACTIVELY USED k-mer index
│   │   ├── kmer_index.bin              # Binary k-mer hash table (15-mers)
│   │   ├── sequences.bin               # Reference sequence database
│   │   ├── index_metadata.json         # Index metadata and statistics
│   │   └── debug/                      # Validation and analysis files
│   └── resistance_db/                   # Alternative resistance database format
├── output/
│   └── GeneFiles/                       # JSON files per species/gene (input for kmer builder)
│       ├── Escherichia_coli/           # E. coli resistance genes
│       ├── Pseudomonas_aeruginosa/     # P. aeruginosa resistance genes
│       └── [other species]/            # Additional organism gene files
├── runtime/                             # ✅ PRODUCTION CODE (used by CMake)
│   └── kernels/resistance/
│       ├── kmer_screening.cu           # Stage 1: K-mer filtering (GPU)
│       ├── fq_mutation_detector.cu     # Stage 2: Alignment/mutation detection (GPU)
│       ├── fq_mutation_detector.cuh    # CUDA header definitions
│       └── fq_pipeline_host.cpp        # Main pipeline orchestrator (CPU)
├── src/                                 # Development/experimental versions
│   ├── kernels/resistance/             # Development CUDA kernels (not used in build)
│   ├── R/                              # R tools
│   │   ├── parse_quinolone_mutations.R # ✅parse files downloaded from MicroBIGG-E
│   ├── python/                         # Python tools and builders
│   │   ├── enhanced_kmer_builder.py    # ✅ K-mer index builder (called by CMake)
│   │   ├── download_ncbi_20250529.py   # NCBI sequence downloader
│   │   ├── generate_synthetic_reads.py # Test data generator
│   │   └── index_validator.py          # Index validation tool
│   └── compiler/                       # BioGPU language compiler components
├── backup_scripts/
│   └── tools/
│       ├── build_fq_resistance_db_adapted.py  # Mutation DB builder (adapted for CSV)
│       ├── parse_quinolone_mutations.R        # ✅ ORIGINAL mutation CSV creator
│       ├── build_kmer_index.py               # Alternative k-mer builders
│       └── [other database builders]         # Additional build tools
├── include/biogpu/                      # Header files for C++/CUDA compilation
├── tests/                               # Test scripts and validation tools
│   ├── debug_gene_species_ids.cpp      # Debug tool for ID assignments
│   └── [test files]                    # Validation and integration tests
├── CMakeLists.txt                       # ✅ Build system (uses /runtime/ kernels)
└── build/                               # Build artifacts (created by CMake)
    ├── fq_pipeline_gpu                 # Main executable
    └── debug_ids                       # Debug utilities

Key Build Relationships:

✅ Active Pipeline (CMake builds):
- Source: /runtime/kernels/resistance/
- Executable: build/fq_pipeline_gpu
- Index Builder: src/python/enhanced_kmer_builder.py
- Data Input: data/fq_resistance_index/kmer_index.bin

🔧 Development/Backup:
- Development kernels: /src/kernels/resistance/ (not compiled)
- Database builders: /backup_scripts/tools/ (Python/R scripts)
- Unused database: data/gpu_resistance_db/ (rich mutation data, not integrated)
```

---

## 🔄 Current Pipeline Workflow

### Stage 0: Data Acquisition
**Purpose**: Download reference sequences for resistance genes from multiple species

The `output/GeneFiles` directory serves as the reference sequence database for the fluoroquinolone resistance detection pipeline.

**How it's created**:
```bash
# Download sequences from NCBI
python src/python/download_ncbi_20250529.py \
    data/Known_Quinolone_Changes.csv \
    data/Known_Efflux_Pump_Genes.csv \
    output/GeneFiles \
    --email your.email@example.com \
    --max-per-gene 300
```

**Role in workflow**:
1. Downloads GenBank sequences for resistance genes (gyrA, gyrB, parC, parE, efflux pumps, etc.) across multiple species
2. Saves JSON files with sequences, annotations, and metadata
3. Feeds into the k-mer index builder (`enhanced_kmer_builder.py`)
4. Creates the binary index used by the GPU resistance detection pipeline

**Output**: 
- `output/GeneFiles/[Species]/[gene].json` files
- Contains: GenBank sequences with CDS features for resistance genes
- Each JSON file includes sequence data, gene annotations, and species metadata

**Source Rationale**: Need comprehensive sequence database across multiple species to detect resistance genes in metagenomic samples

---

### Stage 1: Build K-mer Index
**Purpose**: Create GPU-optimized k-mer index for rapid sequence screening

```bash
python src/python/enhanced_kmer_builder.py \
    output/GeneFiles \
    data/Known_Quinolone_Changes.csv \
    data/fq_resistance_index \
    --kmer-length 15
```

**Output**: `data/fq_resistance_index/`
- `kmer_index.bin` - Binary k-mer lookup table
- `sequences.bin` - Reference sequences
- `index_metadata.json` - Index statistics
- `debug/` - Validation reports

**Key Features**:
- 15-mer index of all resistance genes
- Species and gene ID tracking
- QRDR region extraction
- ~230K unique k-mers indexed

**Status**: ✅ ACTIVELY USED by GPU pipeline

---

### Stage 2: Build Mutation Database (Currently Unused)
**Purpose**: Create comprehensive resistance mutation database

```bash
python backup_scripts/tools/build_fq_resistance_db_adapted.py \
    --mutations data/Known_Quinolone_Changes.csv \
    --genes data/Known_Efflux_Pump_Genes.csv \
    --output data/gpu_resistance_db
```

**Output**: `data/gpu_resistance_db/`
- `mutations.pkl` - 59 resistance mutations across 19 species
- `genes.pkl` - 117 resistance genes
- `gpu_mutation_data.npz` - GPU-optimized arrays
- `alignment_templates.npz` - Position-specific scoring matrices
- `target_regions.json` - QRDR coordinate definitions
- `species_profiles.json` - Per-species resistance summaries

**Status**: ❌ NOT CURRENTLY USED
- Created by Python but no C++ loaders implemented
- Contains rich data for precise mutation detection
- TODO: Integrate into GPU pipeline

---

### Stage 3: GPU Pipeline Execution
**Purpose**: Process metagenomic reads to detect resistance

```bash
# Build the pipeline
mkdir build && cd build
cmake ..
make -j8

# Run on sample data
./fq_pipeline_gpu \
    ../data/fq_resistance_index \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    results.json
```

**Current Implementation**:

#### Stage 3.1: K-mer Screening (`kmer_screening.cu`)
- Function: `enhanced_kmer_filter_kernel`
- Input: Raw FASTQ reads
- Process: 
  1. Extract 15-mers from each read
  2. Binary search in k-mer index
  3. Track hits by gene/species
- Output: Candidate reads with k-mer hits

#### Stage 3.2: Alignment (`fq_mutation_detector.cu`)
- Function: `simple_alignment_kernel`
- Input: Candidate reads from Stage 3.1
- Process:
  1. Score based on k-mer hit density
  2. Calculate simple identity metric
  3. Flag high-scoring alignments
- Output: Alignment results with gene/species attribution

**Current Limitations**:
- No actual sequence alignment (just k-mer counting)
- No codon-aware mutation detection
- No use of mutation database
- No species disambiguation

---

## 🚧 The Gap: What's Not Connected

### Available but Unused Data:
1. **Mutation Database** (`gpu_resistance_db/`)
   - Known resistance positions (e.g., gyrA S83L)
   - Species-specific mutation patterns
   - Position weight matrices

2. **QRDR Regions** (`target_regions.json`)
   - Exact coordinates of resistance-determining regions
   - Could focus alignment efforts

3. **Species Profiles** (`species_profiles.json`)
   - Species-specific resistance patterns
   - Could improve species attribution

### Missing Implementation:
```cpp
// Needed in C++/CUDA:
class MutationDatabaseLoader {
    void loadGPUMutationData(const char* npz_path);
    void loadTargetRegions(const char* json_path);
    void loadScoringMatrices(const char* npz_path);
};

// Enhanced mutation detection kernel
__global__ void targeted_mutation_detection_kernel(
    const char* reads,
    const MutationDatabase* db,
    const TargetRegions* qrdr,
    MutationResults* results
);
```

---

## 📋 Workflow Versioning Strategy

### Version Tracking
Each major change should update:
1. This document
2. `CMakeLists.txt` version number
3. Git tag (e.g., `v0.2.0-kmer-only`)

### Change Log Format
```markdown
## Version 0.3.0 (Date)
### Added
- Integrated mutation database loading
- Codon-aware alignment kernel
### Changed
- Enhanced k-mer screening with Bloom filter
### Fixed
- Species ID assignment bug
```

### Recreating from Scratch
```bash
# 1. Clone repository
git clone [repo] && cd biogpu

# 2. Download reference data
python src/python/download_ncbi_20250529.py ...

# 3. Build indices
python src/python/enhanced_kmer_builder.py ...
python backup_scripts/tools/build_fq_resistance_db_adapted.py ...

# 4. Compile GPU code
mkdir build && cd build
cmake .. && make -j8

# 5. Run pipeline
./fq_pipeline_gpu ...
```

---

## 🎯 Next Steps Priority

### Immediate (Version 0.3.0)
1. [ ✅] Add Bloom filter pre-screening : 5/31/25:  bloom_filter.cu, bloom_filter_integration.cpp
2. [ ] Implement species tracking through pipeline

### Short-term (Version 0.4.0)
4. [ ] Add C++ loaders for mutation database
5. [ ] Implement translated (6-frame) search
6. [ ] Add proper alignment kernel (banded SW or k-mer extension)

### Medium-term (Version 0.5.0)
7. [ ] Multi-species disambiguation
8. [ ] ML feature extraction
9. [ ] Clinical report generation

---

## 🔧 Development Guidelines

### Adding New Features
1. Document the purpose and rationale
2. Update this workflow document
3. Add integration tests
4. Tag version after testing

### File Naming Conventions
- Python scripts: `snake_case.py`
- CUDA kernels: `snake_case.cu`
- C++ headers: `snake_case.h` or `.cuh` for CUDA

### Testing Checklist
- [ ] Synthetic data validation
- [ ] Real sample testing (VRE12)
- [ ] Performance benchmarking
- [ ] Memory leak checking

---


### Technical Documentation
- CUDA Programming Guide
- BioGPU DSL Specification (internal)

---

## 🤝 Contributors
- David Haslam - Pipeline architecture, clinical integration
-

---

*This is a living document. Update with each significant change.*
