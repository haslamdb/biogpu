# BioGPU Archive: 2025-06-08 - Fully Functional FQ Resistance Detection

This archive represents a milestone in the BioGPU project with fully functional fluoroquinolone (FQ) resistance detection capabilities.

## Key Features Working in This Version

1. **Complete FQ Resistance Detection Pipeline**
   - Accurate species and gene identification
   - QRDR (Quinolone Resistance-Determining Region) detection
   - Known FQ resistance mutation identification
   - Clinical confidence scoring

2. **Clinical Report Generation**
   - HTML, JSON, and text format reports
   - Clinical interpretation with confidence levels
   - Species-level resistance summaries
   - Detailed mutation analysis

3. **Performance Metrics**
   - ~16,667-21,739 reads/second on NVIDIA TITAN Xp
   - Efficient GPU-accelerated protein search
   - Bloom filter pre-screening (though currently showing 0% pass rate)

## Major Components

### Core Pipeline Files
- `runtime/kernels/resistance/clean_resistance_pipeline_main.cpp` - Main pipeline orchestrator
- `runtime/kernels/resistance/clinical_fq_report_generator.cpp` - Clinical report generation
- `runtime/kernels/resistance/global_fq_resistance_mapper.h/cpp` - FQ resistance database interface
- `runtime/kernels/resistance/fq_resistance_positions.h` - QRDR position definitions

### Key Improvements in This Version
- Fixed gene/species ID mapping that was broken when disabling fq_mutation_reporter
- Proper QRDR detection on host side after protein matches
- Clinical report generator that doesn't filter mutations
- Performance reporting (reads/second)

## Example Commands

### Build the Pipeline
```bash
mkdir -p build && cd build
cmake ..
make clean_resistance_pipeline -j8
cd ..
```

### Run FQ Resistance Detection on Synthetic Data
```bash
./build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    1M_synthetic_reads_R1.fastq.gz \
    1M_synthetic_reads_R2.fastq.gz \
    results/1M_reads_test \
    data/quinolone_resistance_mutation_table.csv
```

### Run on Real Data
```bash
./build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    /path/to/your_reads_R1.fastq.gz \
    /path/to/your_reads_R2.fastq.gz \
    results/your_sample_name \
    data/quinolone_resistance_mutation_table.csv
```

### Build Integrated Database (if needed)
```bash
# Build nucleotide database
python src/python/build_clean_dynamic_database.py \
    --species_list "Escherichia_coli,Klebsiella_pneumoniae,Staphylococcus_aureus" \
    --gene_list "gyrA,gyrB,parC,parE,grlA,grlB" \
    --output_dir data/integrated_clean_db \
    --kmer_length 15

# Build protein database
python src/python/build_wildtype_protein_db.py \
    --input_dir data/wildtype_protein_seqs \
    --output_dir data/integrated_clean_db/protein \
    --kmer_length 8
```

### Generate Synthetic Test Data
```bash
python src/python/generate_synthetic_reads.py \
    --reference data/reference_genomes/ecoli.fasta \
    --num_reads 1000000 \
    --read_length 150 \
    --mutation_rate 0.001 \
    --output_prefix synthetic_reads
```

## Output Files

The pipeline generates several output files:

1. **CSV File** (`*_protein_matches.csv`)
   - All protein matches with mutation details
   - Includes is_qrdr_alignment flag
   - For validation and debugging

2. **Clinical Reports**
   - `*_clinical_fq_report.html` - Web-viewable report for clinicians
   - `*_clinical_fq_report.json` - Machine-readable results
   - `*_clinical_fq_report.txt` - Text summary

3. **HDF5 File** (`*.h5`)
   - Detailed alignment data
   - For downstream analysis

4. **JSON Summary** (`*.json`)
   - Pipeline statistics
   - Overall results summary

## Example Output Interpretation

### High Confidence Resistance (from clinical report):
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

### Performance Metrics:
```
=== PROCESSING COMPLETE ===
Total reads: 1000000
Protein matches: 2876829
FQ resistance mutations: 658
Processing time: 60 seconds
Performance: 16667 reads/second
```

## Known Issues in This Version

1. Bloom filter showing 0% pass rate - needs investigation
2. CSV parsing error for qnrB19 gene entries with "NA" position
3. No nucleotide matches being reported (possibly related to bloom filter issue)

## Database Requirements

The pipeline requires:
1. Nucleotide k-mer index in `data/integrated_clean_db/nucleotide`
2. Protein database in `data/integrated_clean_db/protein`
3. FQ resistance mutation CSV (e.g., `data/quinolone_resistance_mutation_table.csv`)

## System Requirements

- CUDA-capable GPU (tested on NVIDIA TITAN Xp)
- CUDA toolkit installed
- HDF5 libraries
- zlib for compressed FASTQ support

## Version Info
- BioGPU Version: 0.5.0
- Archive Date: 2025-06-08
- Status: Fully functional FQ resistance detection