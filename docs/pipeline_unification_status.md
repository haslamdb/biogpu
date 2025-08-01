# Pipeline Unification Status

## Overview

This document tracks the progress of unifying the bioinformatics pipeline to quantitate antibiotic resistance genes. The goal is to merge the `runtime/kernels/resistance` and `runtime/kernels/genes` pipelines to use shared resources instead of duplicate effort.

## Current Architecture

### Existing Pipelines
1. **Resistance Pipeline** (`runtime/kernels/resistance/`)
   - Detects specific resistance mutations
   - Uses bloom filters for k-mer screening
   - Processes FASTQ files independently

2. **Genes Pipeline** (`runtime/kernels/genes/`)
   - Detects AMR (Antimicrobial Resistance) genes
   - Uses translated search against protein database
   - Processes FASTQ files independently

### Unification Progress Indicator
The file `kernels/pipeline_tester_with_bloom.cpp` shows the beginning of the unified approach:
- Implements shared FASTQ reading with batch processing
- Integrates bloom filter pre-screening capability
- Demonstrates GPU memory management patterns
- Uses shared data structures (`ReadBatch`, `FastqRecord`)

## Identified Duplications

### 1. Sample CSV Parser
**Status**: Complete duplication
- `/runtime/kernels/resistance/sample_csv_parser.{h,cpp}`
- `/runtime/kernels/genes/sample_csv_parser.{h,cpp}`
- Both implementations are identical

### 2. Bloom Filter Implementation
**Status**: Complete duplication
- `/runtime/kernels/resistance/bloom_filter.cu`
- `/runtime/kernels/genes/bloom_filter.cu`
- Identical CUDA implementations

### 3. FASTQ Processing
**Status**: Independent implementations
- Resistance: `fq_pipeline_host.cpp`
- Genes: `amr_detection_pipeline.cpp`
- Each reads FASTQ files separately, causing redundant I/O

### 4. Common Data Structures
**Status**: Partially shared
- `FastqRecord` and `ReadBatch` structures exist in `fastq_reader.h`
- K-mer and minimizer generation code is duplicated

## Recommended Unification Architecture

```
Unified Pipeline Entry Point
├── Sample CSV Parsing (shared)
├── FASTQ Reading (single pass)
│   ├── Batch loading
│   └── GPU memory transfer
├── K-mer/Minimizer Generation (once)
├── Bloom Filter Screening (optional, shared)
└── Parallel Processing
    ├── Resistance Detection Module
    │   └── Mutation-specific analysis
    ├── AMR Gene Detection Module
    │   └── Translated protein search
    └── [Future] Taxonomic Assignment Module
```

## Implementation Steps

### Phase 1: Create Shared Infrastructure
1. **Create shared directory structure**:
   ```bash
   mkdir -p runtime/kernels/shared
   ```

2. **Move common components**:
   - `sample_csv_parser.{h,cpp}` → `shared/`
   - `bloom_filter.cu` → `shared/`
   - `fastq_reader.{h,cpp}` → `shared/`

3. **Update CMakeLists.txt files** in both pipelines to reference shared components

### Phase 2: Unify Data Processing
1. **Create unified pipeline entry point**: `unified_amr_pipeline.cpp`
   - Single FASTQ reading pass
   - Shared k-mer/minimizer generation
   - Common GPU memory management

2. **Refactor existing pipelines as modules**:
   - Extract resistance detection logic into callable module
   - Extract AMR gene detection logic into callable module
   - Define clear interfaces for data passing

### Phase 3: Optimize Resource Usage
1. **Implement shared GPU memory pools**:
   - Single allocation for read data
   - Reusable k-mer/minimizer buffers
   - Shared bloom filter on GPU

2. **Batch processing optimization**:
   - Process larger batches through all modules
   - Minimize GPU memory transfers

## Current State Summary

### What's Working
- `pipeline_tester_with_bloom.cpp` demonstrates the unified approach
- Both pipelines function independently
- Shared data structures are defined

### What Needs Work
- Duplicate code removal
- Single-pass FASTQ processing
- Unified reporting system
- Integration of both detection methods

### Known Issues (from genes/CLAUDE.md)
- GPU memory errors in translated search need resolution
- Protein database is reloaded for each batch (inefficient)
- Bloom filter integration needs testing with unified pipeline

## Next Actions

1. **Immediate** (can be done now):
   - Create `shared/` directory
   - Move duplicate files
   - Update build configurations

2. **Short-term** (1-2 days):
   - Implement unified pipeline entry point
   - Test integrated processing
   - Benchmark performance improvements

3. **Medium-term** (1 week):
   - Optimize GPU memory usage
   - Add taxonomic assignment module
   - Create unified reporting system

## Benefits of Unification

1. **Performance**: 
   - Single FASTQ read pass (50% I/O reduction)
   - Shared k-mer generation (30-40% compute reduction)
   - Better GPU memory utilization

2. **Maintainability**:
   - Single source of truth for common components
   - Easier debugging and updates
   - Reduced code duplication

3. **Extensibility**:
   - Easy to add new analysis modules
   - Shared infrastructure for future pipelines
   - Consistent data formats

## Testing Strategy

1. **Validation Tests**:
   - Ensure unified pipeline produces same results as individual pipelines
   - Compare performance metrics
   - Verify memory usage improvements

2. **Integration Tests**:
   - Test with various FASTQ file sizes
   - Validate paired-end and single-end processing
   - Check bloom filter effectiveness

3. **Regression Tests**:
   - Maintain test cases from both pipelines
   - Ensure no functionality is lost
   - Monitor for performance regressions

## References

- Pipeline tester implementation: `pipeline_tester_with_bloom.cpp`
- AMR gene pipeline notes: `genes/CLAUDE.md`
- Resistance pipeline: `resistance/clean_resistance_pipeline_main.cpp`
- AMR detection pipeline: `genes/amr_detection_pipeline.cpp`