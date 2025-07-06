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
- Bloom filter is now DISABLED by default (use `--use-bloom-filter` flag to enable)
- Fixed paired-end read ID matching for new Illumina format
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

## Current Issues (2025-01-06)

### GPU Memory Error in Protein Database Loading
**Error**: "CUDA error before allocation: invalid argument" when loading protein database

**Symptoms**:
- Occurs on every batch during `performTranslatedAlignment()`
- Error happens BEFORE the protein k-mer allocation, suggesting an earlier operation failed
- The error persists across all batches (not a memory exhaustion issue)

**What We've Done**:
1. Added `--use-bloom-filter` flag (bloom filter disabled by default)
2. Fixed paired-end read ID matching for new Illumina format (e.g., "1:N:0" vs "2:N:0")
3. Added error checking and synchronization around GPU operations
4. Added null pointer checks for GPU memory

**Next Debugging Steps**:
1. **Check earlier GPU operations**: The "invalid argument" error suggests a previous CUDA call failed
   - Look at minimizer generation kernel calls
   - Check read copying to GPU
   - Verify all GPU memory allocations in `allocateGPUMemory()`

2. **Investigate the gene entries copy**: 
   - The error occurs around `cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(), ...)`
   - Check if `d_gene_entries` is properly allocated in `NCBIAMRDatabaseLoader`
   - Verify the size calculations are correct

3. **Check for uninitialized pointers**:
   - Some GPU pointers might not be initialized when bloom filter is disabled
   - Look for any conditional allocations that might leave pointers null

4. **Memory alignment issues**:
   - Check if data structures have proper alignment for GPU access
   - Verify sizeof(AMRGeneEntry) matches between CPU and GPU code

5. **Consider running with cuda-memcheck**:
   ```bash
   cuda-memcheck ./amr_detection AMR_CDS.fa,AMRProt.fa test_samples.csv amr_results/
   ```

**Temporary Workaround**: 
- The pipeline continues to process despite the error, but no AMR genes are detected
- Each batch independently fails to load the protein database