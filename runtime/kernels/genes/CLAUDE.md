# AMR Detection Pipeline - Next Steps

## TODO for Next Session

### 1. Unify Sample Processing
- **Issue**: Currently using a copy of `sample_csv_parser.cpp` from `../resistance`
- **Action**: Reference the original file in `../resistance` directly instead of copying
- **Benefit**: Ensures consistency across pipelines

### 2. Integrate with Existing Read Processing
- **Current State**: Three separate workflows that each read FASTQ files independently:
  1. FQ resistance detection (`../resistance`)
  2. AMR gene detection (this directory)
  3. Taxonomic assignment
- **Goal**: Single pass over sequence data that feeds all three pipelines
- **Implementation Strategy**:
  - Import reads once into GPU memory
  - Generate k-mers/minimizers once
  - Share these data structures across all three detection pipelines

### 3. Code Reuse Opportunities
- Look for existing k-mer/minimizer generation code in other kernels
- Identify common sequence processing patterns
- Create shared utilities for:
  - FASTQ reading and batching
  - K-mer generation
  - Minimizer extraction
  - Bloom filter construction

### 4. Unified Pipeline Architecture
```
FASTQ Input → Read Processing → K-mer/Minimizer Generation
                                           ↓
                    ┌──────────────────────┼──────────────────────┐
                    ↓                      ↓                      ↓
            FQ Resistance           AMR Gene Detection      Taxonomic Assignment
                    ↓                      ↓                      ↓
                    └──────────────────────┼──────────────────────┘
                                           ↓
                                   Integrated Report
```

### 5. Current Working State
- AMR detection pipeline is functional and tested
- Successfully processes FASTQ files
- Generates bloom filters and performs translated search
- Outputs TSV files and clinical reports

### 6. File Locations for Reference
- Sample processor: `../resistance/sample_csv_parser.cpp`
- FQ pipeline: `../resistance/clean_resistance_pipeline_main.cpp`
- AMR pipeline: `./amr_detection_pipeline.cpp`
- Kernel functions: `./amr_detection_kernels.cu`

### 7. Integration Points to Consider
- Shared GPU memory allocations
- Common data structures (ReadBatch, Minimizer, etc.)
- Unified configuration management
- Combined output reporting