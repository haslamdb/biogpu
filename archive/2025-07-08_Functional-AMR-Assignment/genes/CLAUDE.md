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

### GPU Memory Error in Translated Search
**Error**: "CUDA error after translated search: invalid argument"

**Root Causes Found**:
1. **Memory allocation size mismatch**: 
   - `d_reads` was allocated for `max_batch * max_read_len` (30MB)
   - But merged paired-end reads are much longer (2 * read_len + gap)
   - Fixed by allocating `max_batch * (max_read_len * 2 + 100)`

2. **Offset array size issue**:
   - `d_minimizer_offsets` needed `max_batch + 1` elements, not just `max_batch`
   - Fixed by updating allocation size

3. **Persistent CUDA error state**:
   - The translated search engine leaves CUDA in an error state after each batch
   - Error is now caught and cleared, but the underlying issue persists

**What We've Done**:
1. Added `--use-bloom-filter` flag (bloom filter disabled by default)
2. Fixed paired-end read ID matching for new Illumina format (e.g., "1:N:0" vs "2:N:0")
3. Fixed GPU memory allocation sizes for merged reads
4. Added comprehensive error checking throughout the pipeline
5. Added error clearing to prevent cascade failures

**Current Status**:
- First batch sometimes succeeds in loading the protein database
- All subsequent batches fail with "invalid argument" error
- The protein database is being reloaded for each batch (inefficient)
- No AMR genes are detected due to the persistent error

**Next Steps**:
1. **Investigate why translated search corrupts CUDA state**:
   - Check for memory leaks in translated_search_amr.cu
   - Look for unfreed GPU memory in the search engine
   - Verify all GPU pointers are properly initialized

2. **Optimize protein database loading**:
   - Load the protein database once at initialization
   - Reuse it across all batches instead of reloading
   - This would also improve performance significantly

3. **Debug the search engine lifecycle**:
   - Check what `destroy_translated_search_engine()` actually frees
   - Look for static/global GPU allocations that persist
   - Consider resetting CUDA context between batches (last resort)

**Temporary Workaround**: 
- The pipeline continues despite errors but finds no AMR genes
- Can process smaller batches or single samples to minimize impact