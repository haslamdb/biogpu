# BioGPU Runtime Kernels Cleanup Summary
**Date**: August 2, 2025

## Overview
This directory contains files that were moved from the active runtime/kernels directories during cleanup. These files are no longer used by the current resistance and AMR gene detection pipelines.

## Resistance Pipeline Cleanup

### Moved Files:
1. **Deprecated/Old Versions:**
   - `resistance_pipeline_v2.cpp/.h` - Old version of pipeline
   - `fq_pipeline_host_backup.cpp` - Backup of old host code
   - `integrate_resistance_pipeline.cpp` - Old integration attempt
   - `pipeline_integration_example.cpp` - Example code
   - `gpu_diagnostic_adapters.cu/.h` - Deprecated diagnostic code

2. **Unused Headers/Wrappers:**
   - `clean_resistance_pipeline.h` - Not referenced
   - `CleanResistancePipeline.h` - Not referenced
   - `CleanResistancePipelineWrapper.cpp` - Not referenced
   - `resistance_detector.cpp/.h` - Not used in current pipeline

3. **Test Files:**
   - `test_bloom_filter.cu` - Commented out in CMake
   - `test_bloom_filter.cu.backup` - Backup of test

4. **Build Files:**
   - `CMakeLists_backup_20250705.txt` - Old CMake backup
   - `CMakeLists_working.txt` - Another CMake backup
   - `build_resistance_fixed.sh` - Alternative build script
   - `build_resistance_minimal.sh` - Alternative build script
   - `build_resistance_stub.sh` - Alternative build script
   - `build_resistance_with_deps.sh` - Alternative build script

5. **Build Directories:**
   - `build_cmake/` - Old build directory
   - `build_fixed/` - Alternative build directory
   - `build_minimal/` - Alternative build directory

### Active Files Kept:
- Core pipeline: `clean_resistance_pipeline_main.cpp`, `fq_mutation_detector.cu`, `fq_mutation_reporter.cpp`
- Libraries: `translated_search_revised.cu`, `kmer_screening.cu`, `resistance_detection_gpu.cu`
- Support: `hdf5_alignment_writer.cpp`, `global_fq_resistance_mapper.cpp`, etc.

## AMR Gene Pipeline Cleanup

### Moved Files:
1. **Duplicate Files (using shared/ versions):**
   - `bloom_filter.cu` - Using shared/bloom_filter.cu
   - `sample_csv_parser.cpp/.h` - Using shared/sample_csv_parser.*

2. **Unused Components:**
   - `amr_gene_detector.cpp/.h` - Not referenced in CMake
   - `arg_database.h` - Not used
   - `build_arg_database.cpp` - Not referenced
   - `ncbi_amr_database_loader.h` - Not used
   - `translated_search_amr_fixed.h` - Appears to be a backup

### Active Files Kept:
- Main pipeline: `amr_detection_main.cpp`, `amr_detection_pipeline.cpp`
- Kernels: `amr_detection_kernels.cu`, `amr_detection_kernels_wrapper.cu`
- Support: `translated_search_amr.cu`, `hdf5_amr_writer.cpp`, etc.

## Shared Directory
All files in shared/ are actively used and were kept:
- `bloom_filter.cu` - Used by both pipelines
- `sample_csv_parser.cpp/.h` - Used by both pipelines
- Various interface headers

## Notes
- The `deprecated/` and `experimental/` directories in resistance were left as-is
- The `build/` directory in resistance contains the active build and was kept
- Test output files (*.json, *.h5, etc.) can be cleaned up separately as needed
- Database files and reference sequences were not touched