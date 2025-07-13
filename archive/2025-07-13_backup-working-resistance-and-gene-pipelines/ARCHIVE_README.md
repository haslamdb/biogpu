# Archive: Working Resistance and Gene Pipelines
**Date**: July 13, 2025

## Purpose
This archive captures the current working state of both the fluoroquinolone resistance detection pipeline and the AMR gene detection pipeline before attempting to unify them with a single main file.

## Archive Contents

### resistance/
Contains the fluoroquinolone resistance detection pipeline code:
- Clean resistance pipeline implementation
- Clinical report generation  
- Enhanced mutation detection kernels
- HDF5 output writers
- Bloom filter implementation
- K-mer screening functionality
- Sample CSV parser
- Global FQ resistance mapper
- Diagnostic report generation

### genes/
Contains the AMR gene detection pipeline code:
- AMR detection pipeline and kernels
- Translated protein search implementation
- Clinical AMR report generator
- NCBI AMR database loader
- HDF5 AMR writer
- ARG (Antimicrobial Resistance Gene) database structures
- Bloom filter implementation
- Sample CSV parser

## Current Status
Both pipelines are fully functional and operational as separate entities:

1. **Resistance Pipeline**: 
   - Detects fluoroquinolone resistance mutations
   - Generates clinical reports in HTML/JSON formats
   - Uses enhanced mutation detection with GPU acceleration
   - Processes FASTQ files through k-mer screening and alignment

2. **Gene Pipeline**:
   - Detects AMR genes from FASTQ reads
   - Performs translated protein search
   - Generates clinical AMR reports
   - Integrates with NCBI AMR database

## Build Status
- **resistance**: Uses CMakeLists.txt for building
- **genes**: Uses Makefile for building

Both pipelines compile and run successfully with their respective build systems.

## Reason for Archive
This archive was created to preserve the working state before attempting to:
1. Create a unified main file that can run both pipelines
2. Potentially integrate the pipelines for single-pass analysis
3. Share common components between the pipelines

## Next Steps
- Create unified main file to run both pipelines
- Identify shared components for potential consolidation
- Optimize for single-pass read processing