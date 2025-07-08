# Archive: Functional AMR Assignment
**Date**: July 8, 2025

## Purpose
This archive captures the functional state of the AMR (Antimicrobial Resistance) detection pipeline after achieving working AMR gene assignment functionality.

## Archive Contents

### resistance/
Contains the fluoroquinolone resistance detection pipeline code:
- Clean resistance pipeline implementation
- Clinical report generation
- Mutation detection kernels
- HDF5 output writers
- Bloom filter implementation
- Sample CSV parser

### genes/
Contains the AMR gene detection pipeline code:
- AMR detection pipeline and kernels
- Translated protein search implementation
- Clinical AMR report generator
- NCBI AMR database loader
- HDF5 AMR writer
- Argument database structures

## Key Achievements at Archive Time
1. **AMR Gene Detection**: Successfully implemented functional AMR gene detection from FASTQ reads
2. **Protein Database Integration**: Working protein database loading and search
3. **Clinical Reports**: Generating both HTML and JSON clinical reports for AMR findings
4. **Pipeline Integration**: Both resistance and genes pipelines are operational

## Known Issues (as of archive date)
- GPU memory error in translated search persists but pipeline continues
- Bloom filter disabled by default (use --use-bloom-filter flag to enable)
- Protein database reloaded for each batch (inefficient)

## Build Status
Both pipelines compile and run successfully with the standard build process:
- resistance: Uses CMakeLists.txt
- genes: Uses Makefile

## Next Steps (documented in CLAUDE.md)
- Unify sample processing between pipelines
- Integrate with existing read processing for single-pass analysis
- Optimize protein database loading (load once, reuse across batches)
- Fix GPU memory corruption in translated search engine