# Changes Summary - June 10, 2025
## Since Last Archive (June 8, 2025)

### Major New Features

#### 1. Batch Processing Support
- Added CSV-based batch processing mode with `--csv` option
- Implemented `sample_csv_parser.cpp/h` library for flexible CSV parsing
- Created `BatchProcessor` class for managing multi-sample workflows
- Added `generate_sample_csv.py` utility to create CSV files from directories
- Supports validation with `--dry-run` before processing
- Error handling with `--stop-on-error` option

#### 2. Enhanced Allele Frequency Analysis
- Added configurable depth filtering:
  - `--min-allele-depth`: Minimum depth for analysis (default: 5)
  - `--min-report-depth`: Minimum depth for reporting (default: 0)
- Improved allele frequency calculations with depth-aware statistics
- Better handling of low-coverage positions

#### 3. Automatic Sample Management
- Automatic sample name extraction from R1 filename
- Removed output prefix as required argument
- Default output structure: `results/<sample_name>/`
- Simplified command-line interface

### Code Changes

#### Modified Files
1. **clean_resistance_pipeline_main.cpp**
   - Added `processSingleSample()` function for modular processing
   - Implemented command-line parsing for batch mode
   - Added CSV file detection and processing logic
   - Integrated `BatchProcessor` for multi-sample handling

2. **CMakeLists.txt**
   - Added `sample_csv_parser` library
   - Linked CSV parser to clean_resistance_pipeline
   - Updated build dependencies

3. **docs/workflow_documentation.md**
   - Added batch processing documentation
   - Updated command examples for both modes
   - Added CSV format specification
   - Updated version history for v0.6.3

#### New Files
1. **sample_csv_parser.cpp/h**
   - Flexible CSV parsing with column name detection
   - Path expansion (home directory, relative paths)
   - File existence validation
   - Support for paired-end and single-end samples

2. **generate_sample_csv.py**
   - Automatically creates CSV from directory of FASTQ files
   - Pairs R1/R2 files based on naming patterns
   - Supports recursive directory scanning
   - Flexible output options

3. **example_csv.csv**
   - Sample CSV file for testing
   - Demonstrates different path formats

### Bug Fixes
- Fixed const-correctness issues in CSV parser
- Added missing header includes (<functional>, <chrono>)
- Improved error handling for file paths

### Performance Improvements
- Batch mode reduces overhead for processing multiple samples
- Parallel processing preparation (infrastructure in place)
- Maintained single-sample performance characteristics

### Documentation Updates
- Comprehensive batch processing examples
- CSV format documentation with examples
- Updated configuration options
- Added troubleshooting for path issues

### Backward Compatibility
- Single-sample mode fully preserved
- Optional output prefix still supported
- All existing flags continue to work
- Smooth migration path for users

---

This update makes BioGPU significantly more practical for clinical and research laboratories processing multiple samples.