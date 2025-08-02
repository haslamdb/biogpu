# Unified BioGPU Pipeline Integration Summary
**Date**: August 2, 2025

## Overview
This document summarizes the complete integration of the BioGPU resistance mutation detection and AMR gene detection pipelines into a single unified executable that processes FASTQ files once and runs both detection algorithms.

## Key Accomplishments

### 1. Fixed AMR Pipeline CSV Format Issue
**Problem**: The unified pipeline was generating CSV files with incorrect format for the AMR pipeline.
- Original format: `sample_name,fastq_path` with both R1 and R2 in one field
- Required format: `SampleName,FilePath,R1 file,R2 file` with separate columns

**Solution**: Modified `unified_pipeline_process_spawning.cpp` (lines 363-384) to:
- Extract directory path and filenames separately
- Create proper column headers based on paired/single-end reads
- Format CSV correctly for the AMR pipeline parser

### 2. Rebuilt Protein Database with Correct K-mer Size
**Problem**: K-mer size mismatch between protein database (k=5) and resistance pipeline expectation (k=8).
- Error: `[PROTEIN DB ERROR] K-mer size mismatch: expected 8, got 5`

**Solution**: 
1. Modified the protein database build script to accept k-mer size parameter
2. Created `build_protein_resistance_db_k8.py` with proper k-mer length handling
3. Rebuilt the database:
```bash
cd /home/david/projects/biogpu/runtime/kernels/resistance
python3 build_protein_resistance_db_k8.py \
    /home/david/projects/biogpu/data/fq_genes \
    /home/david/projects/biogpu/data/quinolone_resistance_mutation_table.csv \
    /home/david/projects/biogpu/data/protein_resistance_db_new \
    --kmer-length 8
```
4. Replaced the old database:
```bash
cd /home/david/projects/biogpu/data
mv protein_resistance_db protein_resistance_db_old_k5
mv protein_resistance_db_new protein_resistance_db
```

### 3. Fixed Incorrect Protein Database Usage (August 2, 2025)
**Problem**: The resistance pipeline was matching proteins against the wrong database, resulting in:
- Matches against non-FQ genes (e.g., adeB efflux pump) instead of FQ target genes
- 0 FQ resistance mutations detected despite finding protein matches
- The unified pipeline was using `/home/david/projects/biogpu/data/protein_resistance_db` which contains ALL 92 genes from the fq_genes directory

**Root Cause**: The protein database contained all genes including efflux pumps, not just fluoroquinolone target genes (gyrA, parC, gyrB, parE, grlA, grlB).

**Solution**: Updated the unified pipeline to use the correct databases:
```cpp
// OLD (incorrect):
std::string resistance_index = "/home/david/projects/biogpu/data/fq_resistance_index";
std::string protein_db = "/home/david/projects/biogpu/data/protein_resistance_db";

// NEW (correct):
std::string resistance_index = "/home/david/projects/biogpu/data/integrated_clean_db/nucleotide";
std::string protein_db = "/home/david/projects/biogpu/data/integrated_clean_db/protein";
```

The integrated_clean_db protein database:
- Contains only 4 FQ target genes: grlA, gyrA, parC, parE
- Already has the correct k-mer size of 8
- Results in proper FQ resistance mutation detection

### 4. Added EM Algorithm Support to Unified Pipeline
**Enhancement**: Added support for the Expectation-Maximization (EM) algorithm for better AMR gene abundance estimation.

**Changes to `unified_pipeline_process_spawning.cpp`**:
1. Added to PipelineOptions struct (lines 47-48):
```cpp
bool use_em = true;  // Default to enabled
int em_iterations = 10;
```

2. Added command-line options (lines 205-206):
```cpp
{"no-em", no_argument, 0, 0},
{"em-iterations", required_argument, 0, 0},
```

3. Added parsing logic (lines 236-237):
```cpp
else if (opt_name == "no-em") options.use_em = false;
else if (opt_name == "em-iterations") options.em_iterations = std::stoi(optarg);
```

4. Added EM flags to AMR pipeline call (lines 408-412):
```cpp
if (options.use_em) {
    amr_proc.addArg("--em");
    amr_proc.addArg("--em-iterations");
    amr_proc.addArg(std::to_string(options.em_iterations));
}
```

## Final Build Process

### Build Script Used: `build_final_integration.sh`

Here's the complete build script that creates the unified pipeline:

```bash
#!/bin/bash

echo -e "\033[0;32m========================================\033[0m"
echo -e "\033[0;32mBioGPU Final Integration Build\033[0m"
echo -e "\033[0;32m========================================\033[0m"

# Build the unified pipeline executable
echo -e "\033[1;33mBuilding process-spawning unified pipeline...\033[0m"
# Note: After August 2, 2025 fix, CUDA libraries are required for linking
g++ -std=c++17 -O3 unified_pipeline_process_spawning.cpp -o bio_gpu_pipeline_integrated -lstdc++fs -pthread -lcudart -L/usr/local/cuda/lib64

if [ $? -eq 0 ]; then
    echo -e "\033[0;32mBuild successful!\033[0m"
else
    echo -e "\033[0;31mBuild failed!\033[0m"
    exit 1
fi

# Check for pipeline executables
echo -e "\033[1;33mChecking for pipeline executables...\033[0m"

# Check AMR pipeline
if [ -f "../genes/amr_detection" ]; then
    echo -e "\033[0;32m✓ Found AMR detection executable\033[0m"
    AMR_FOUND=true
else
    echo -e "\033[0;31m✗ AMR detection executable not found at ../genes/amr_detection\033[0m"
    AMR_FOUND=false
fi

# Check resistance pipeline  
if [ -f "../resistance/clean_resistance_pipeline" ]; then
    echo -e "\033[0;32m✓ Found resistance pipeline: ../resistance/clean_resistance_pipeline\033[0m"
    RESISTANCE_FOUND=true
else
    echo -e "\033[0;31m✗ Resistance pipeline not found\033[0m"
    RESISTANCE_FOUND=false
fi

# Create test script
echo -e "\033[1;33mCreating test script...\033[0m"
cat > test_integrated_pipeline.sh << 'EOF'
#!/bin/bash
./bio_gpu_pipeline_integrated \
    --r1 ../../../data/test_R1.fastq \
    --r2 ../../../data/test_R2.fastq \
    --output-dir test_output \
    --reference-db ../../../data/AMR_CDS.fa,../../../data/AMR_protein.fa \
    --resistance-db ../../../data/quinolone_resistance_mutation_table.csv \
    --sample-id test_sample \
    --gpu-device 0 \
    --progress-json
EOF
chmod +x test_integrated_pipeline.sh

# Create installation package
echo -e "\033[1;33mCreating installation package...\033[0m"
mkdir -p install_final
cp bio_gpu_pipeline_integrated install_final/
cat > install_final/install.sh << 'EOF'
#!/bin/bash
echo "Installing BioGPU unified pipeline..."
sudo cp bio_gpu_pipeline_integrated /usr/local/bin/bio_gpu_pipeline
sudo chmod +x /usr/local/bin/bio_gpu_pipeline
echo "Installation complete!"
echo "Usage: bio_gpu_pipeline --r1 <reads1> --r2 <reads2> --output-dir <dir> ..."
EOF
chmod +x install_final/install.sh

echo -e "\033[0;32m========================================\033[0m"
echo -e "\033[0;32mBuild Complete!\033[0m"
echo -e "\033[0;32m========================================\033[0m"
echo ""
echo "Executable: $(pwd)/bio_gpu_pipeline_integrated"
echo "Size: $(du -h bio_gpu_pipeline_integrated | cut -f1)"
echo ""
echo "This version uses process spawning to run existing pipelines:"
echo "- Searches for pipeline executables in multiple locations"
echo "- Runs pipelines as separate processes with proper GPU assignment"
echo "- Captures and reports progress in real-time"
echo "- Supports all original command-line arguments"
echo ""
echo "Test with: ./test_integrated_pipeline.sh"
echo "Install with: cd install_final && ./install.sh"
echo ""

if [ "$AMR_FOUND" = true ] && [ "$RESISTANCE_FOUND" = true ]; then
    echo -e "\033[0;32mBoth pipeline executables found - ready for full integration!\033[0m"
else
    echo -e "\033[0;33mWarning: Some pipeline executables are missing\033[0m"
fi
```

## Key Components Modified

### 1. AMR K-mer Database Build
**Location**: `/home/david/projects/biogpu/runtime/kernels/genes/`
**Built with**: `build_amr_protein_db` executable
```bash
cd /home/david/projects/biogpu/runtime/kernels/genes
./build_amr_protein_db AMRProt.fa amr_protein_db
```
**Note**: Created symlink in build_complete directory for runtime access

### 2. Integrated Clean Database Build (Correct FQ Database)
**Purpose**: Build the correct resistance database containing only FQ target genes (gyrA, parC, parE, grlA)
**Location**: `/home/david/projects/biogpu/data/integrated_clean_db/`
**Build command**:
```bash
cd /home/david/projects/biogpu
python3 src/python/build_clean_dynamic_database.py \
    data/wildtype_protein_seqs \  # input dir with wildtype sequences
    data/integrated_clean_db \     # output dir
    --mutations-csv data/quinolone_resistance_mutation_table.csv \
    --kmer-length 8  # k-mer length must be 8 to match PROTEIN_KMER_SIZE in clean_resistance_pipeline
```
**Important notes**:
- K-mers are built from **wildtype sequences only** (not mutated sequences)
- The resistance mutation information is stored as metadata for lookup during detection
- The pipeline identifies mutations by finding mismatches against wildtype k-mers
- K-mer length of 8 is hardcoded in `translated_search_revised.cu`, so this must match
- **Current limitation**: Using a limited set of wildtype sequences from `data/wildtype_protein_seqs/`
  - Only includes 6 species and 4 genes (grlA, gyrA, parC, parE)
  - Should be expanded to include more species and sequence variants for better coverage

### 3. Incorrect Protein Database (Old - DO NOT USE)
**Location**: `/home/david/projects/biogpu/data/protein_resistance_db/`
**Issue**: Contains ALL 92 genes from fq_genes directory including efflux pumps (adeA/B/C, mexA-Y, etc.)
**K-mer size**: Initially k=5, rebuilt to k=8, but still contains wrong genes

### 4. AMR Pipeline Build
**Location**: `/home/david/projects/biogpu/runtime/kernels/genes/`
**Build command**:
```bash
cd /home/david/projects/biogpu/runtime/kernels/genes
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```
**Executable**: `amr_detection`
**Note**: Requires AMR k-mer database built with `build_amr_protein_db`

### 5. Resistance Pipeline Build
**Location**: `/home/david/projects/biogpu/runtime/kernels/resistance/`
**Build command**:
```bash
cd /home/david/projects/biogpu/runtime/kernels/resistance
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```
**Executable**: `clean_resistance_pipeline`
**Note**: The resistance pipeline expects the integrated_clean_db databases as of August 2, 2025

### 6. Unified Pipeline Executable
**Location**: `/home/david/projects/biogpu/runtime/kernels/build_complete/bio_gpu_pipeline_integrated`
**Build command**:
```bash
cd /home/david/projects/biogpu/runtime/kernels/build_complete
g++ -std=c++17 -O3 unified_pipeline_process_spawning.cpp -o bio_gpu_pipeline_integrated \
    -lstdc++fs -pthread -lcudart -L/usr/local/cuda/lib64
```
**Key features**:
- Process spawning architecture
- JSON progress reporting
- Multi-GPU support (assigns different GPUs to each pipeline)
- EM algorithm enabled by default
- Command-line compatibility with original pipelines
- Uses integrated_clean_db databases for correct FQ resistance detection

## Test Results

### Small Test (3 reads)
- Both pipelines ran successfully
- No mutations/genes found (expected for synthetic data)
- Validated CSV format and database loading

### Full Metagenome Test (3.5M reads)
**Before fix (incorrect database)**:
- **Resistance Pipeline**: 8 seconds, 28 protein matches, 81 mutations detected, **0 FQ resistance mutations**
- Incorrectly matching adeB (efflux pump) instead of FQ target genes

**After fix (correct integrated_clean_db)**:
- **Resistance Pipeline**: 9 seconds, 1,326 protein matches, 7,251 mutations detected, **7 FQ resistance mutations**
- Correctly identifying mutations in gyrA, parC genes from K. pneumoniae, S. aureus, etc.
- **AMR Pipeline**: ~107 seconds, 110 AMR genes detected, EM algorithm resolved 14,088 hits
- **Total time**: 118 seconds for complete analysis
- **Performance**: >375K reads/second for resistance, full EM resolution for AMR

## Command-Line Usage

```bash
./bio_gpu_pipeline_integrated \
    --r1 <forward_reads.fastq.gz> \
    --r2 <reverse_reads.fastq.gz> \
    --output-dir <output_directory> \
    --reference-db <AMR_CDS.fa,AMR_protein.fa> \
    --resistance-db <resistance_mutations.csv> \
    --sample-id <sample_name> \
    --gpu-device <gpu_id> \
    --progress-json \
    --threads <num_threads> \
    [--no-em] \
    [--em-iterations <num>] \
    [--disable-resistance] \
    [--disable-amr] \
    [--use-bloom-filter]
```

## Future Modifications Notes

When making modifications:
1. The unified pipeline uses process spawning - modify `unified_pipeline_process_spawning.cpp`
2. Individual pipelines remain separate executables in `../genes/` and `../resistance/`
3. Database paths are currently hardcoded in some places - may want to make fully configurable
4. EM is enabled by default - use `--no-em` to disable
5. Progress reporting uses JSON format when `--progress-json` is specified

## Important File Locations

- Unified pipeline source: `/runtime/kernels/build_complete/unified_pipeline_process_spawning.cpp`
- Build script: `/runtime/kernels/build_complete/build_final_integration.sh`
- AMR pipeline: `/runtime/kernels/genes/amr_detection`
- Resistance pipeline: `/runtime/kernels/resistance/clean_resistance_pipeline`
- AMR k-mer database: `/runtime/kernels/genes/amr_protein_db/`
- **Correct FQ databases (as of August 2, 2025)**:
  - Nucleotide: `/data/integrated_clean_db/nucleotide/`
  - Protein: `/data/integrated_clean_db/protein/` (contains only grlA, gyrA, parC, parE)
- Incorrect protein database (contains all genes): `/data/protein_resistance_db/`
- Test data: `/data/569_A_038_R1.fastq.gz` and `/data/569_A_038_R2.fastq.gz`

## TODO: Phase 5 - Database Improvements

### 5.1 Expand Wildtype Protein Sequences
- [ ] Current database only includes 6 species and 4 genes from `data/wildtype_protein_seqs/`
- [ ] Add more bacterial species that are clinically relevant for FQ resistance:
  - [ ] Neisseria gonorrhoeae
  - [ ] Haemophilus influenzae
  - [ ] Campylobacter species
  - [ ] Salmonella enterica serovars
  - [ ] Additional Enterobacteriaceae
- [ ] Add sequence variants within species to capture natural variation
- [ ] Include gyrB and parE genes where applicable
- [ ] Validate new sequences against known resistance mutations in literature
- [ ] Rebuild integrated_clean_db with expanded wildtype sequences

## TODO: Phase 6 - Integration Testing

### 6.1 Test with Streaming Service
- [ ] Start infrastructure services:
  ```bash
  docker-compose up -d redis kafka zookeeper
  ```
- [ ] Run streaming service:
  ```bash
  source venv/bin/activate
  python3 run_service.py
  ```
- [ ] Submit test job via manifest
- [ ] Verify WebSocket updates are received
- [ ] Check output files are generated correctly

### 6.2 Validation Tests
- [ ] Compare results with individual pipelines
- [ ] Verify both resistance mutations and AMR genes are detected
- [ ] Test with various FASTQ file sizes
- [ ] Benchmark performance improvements

### 6.3 Error Handling
- [ ] Test with corrupted FASTQ files
- [ ] Test with missing reference databases
- [ ] Verify error messages are properly streamed
- [ ] Check retry logic works correctly

## TODO: Phase 7 - Production Deployment

### 7.1 Documentation
- [ ] Update README with unified pipeline usage
- [ ] Document all command-line arguments
- [ ] Create example output format documentation

### 7.2 Performance Tuning
- [ ] Profile GPU utilization
- [ ] Optimize batch sizes based on GPU memory
- [ ] Set appropriate thread counts

### 7.3 Final Verification
- [ ] Run full test suite
- [ ] Process production-scale datasets
- [ ] Monitor resource usage
- [ ] Verify streaming updates work at scale

## Success Criteria

The unified pipeline is complete when:
1. ✓ Single executable handles both resistance and gene detection
2. ✓ Executable is installed at `/usr/local/bin/bio_gpu_pipeline`
3. ✓ Streaming service successfully executes pipeline via Kafka jobs
4. ✓ Real-time progress updates stream via WebSocket
5. ✓ Results match individual pipeline outputs
6. ✓ Performance is improved (single FASTQ read, shared GPU memory)
7. ✓ All tests pass

## Quick Commands Reference

```bash
# Build unified pipeline
cd /home/david/projects/biogpu/runtime/kernels/build_complete
./build_final_integration.sh

# Install
cd install_final
sudo ./install.sh
# This installs to: /usr/local/bin/bio_gpu_pipeline

# Test with small synthetic data
/home/david/projects/biogpu/runtime/kernels/build_complete/bio_gpu_pipeline_integrated \
  --r1 /home/david/projects/biogpu/data/test_R1.fastq \
  --r2 /home/david/projects/biogpu/data/test_R2.fastq \
  --output-dir /tmp/test_output \
  --reference-db "/home/david/projects/biogpu/runtime/kernels/genes/AMR_CDS.fa,/home/david/projects/biogpu/runtime/kernels/genes/AMRProt.fa" \
  --resistance-db /home/david/projects/biogpu/data/quinolone_resistance_mutation_table.csv \
  --sample-id test_sample \
  --gpu-device 0 \
  --progress-json \
  --threads 8

# Test with real metagenome data (3.5M reads)
/home/david/projects/biogpu/runtime/kernels/build_complete/bio_gpu_pipeline_integrated \
  --r1 /home/david/projects/biogpu/data/569_A_038_R1.fastq.gz \
  --r2 /home/david/projects/biogpu/data/569_A_038_R2.fastq.gz \
  --output-dir /tmp/metagenome_test_output \
  --reference-db "/home/david/projects/biogpu/runtime/kernels/genes/AMR_CDS.fa,/home/david/projects/biogpu/runtime/kernels/genes/AMRProt.fa" \
  --resistance-db /home/david/projects/biogpu/data/quinolone_resistance_mutation_table.csv \
  --sample-id 569_A_038 \
  --gpu-device 0 \
  --progress-json \
  --threads 8

# Test with specific options
/home/david/projects/biogpu/runtime/kernels/build_complete/bio_gpu_pipeline_integrated \
  --r1 /home/david/projects/biogpu/data/569_A_038_R1.fastq.gz \
  --r2 /home/david/projects/biogpu/data/569_A_038_R2.fastq.gz \
  --output-dir /tmp/custom_test_output \
  --reference-db "/home/david/projects/biogpu/runtime/kernels/genes/AMR_CDS.fa,/home/david/projects/biogpu/runtime/kernels/genes/AMRProt.fa" \
  --resistance-db /home/david/projects/biogpu/data/quinolone_resistance_mutation_table.csv \
  --sample-id 569_A_038 \
  --gpu-device 0 \
  --progress-json \
  --threads 8 \
  --no-em \
  --disable-resistance  # Run only AMR detection without EM

# Start streaming service
cd /home/david/projects/biogpu/runtime/common/io
source venv/bin/activate
python3 run_service.py
```

## Notes

- The pipeline outputs JSON progress updates to stdout when `--progress-json` is used
- All actual results are written to files in the output directory
- The executable returns 0 on success, non-zero on failure
- GPU device selection is respected for multi-GPU systems (use `--gpu-device` to specify)
- EM algorithm is enabled by default; use `--no-em` to disable
- Both pipelines run by default; use `--disable-resistance` or `--disable-amr` to run only one
- The AMR database path must include both DNA and protein FASTA files separated by comma