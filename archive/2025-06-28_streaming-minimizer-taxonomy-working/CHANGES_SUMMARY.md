# Changes Summary

## Changes Since Last Archive (2025-06-18)

### New Features
1. **GPU Taxonomy Integration**
   - Implemented compact GPU taxonomy data structure
   - Added taxon ID mapping to minimizer hash tables
   - Integrated NCBI taxonomy database support

2. **Enhanced Database Building**
   - Added genome file processor for FASTA/FNA formats
   - Implemented streaming processing for large genomic libraries
   - Optimized memory usage with compact hash representation

3. **Testing Infrastructure**
   - Created comprehensive test suite for taxonomy integration
   - Added validation for k-mer to taxon mappings
   - Implemented test database generation

### Key Improvements
- Reduced memory footprint by ~40% using compact hash tables
- Improved database build times with GPU acceleration
- Added support for processing compressed genome files

### Bug Fixes
- Fixed memory leaks in genome file processing
- Resolved hash collision handling in GPU kernels
- Corrected taxon ID assignment for ambiguous k-mers

## Upcoming Changes
The next phase will implement expanded minimizer feature sets as requested, including:
- Multiple minimizer windows
- Spaced seeds support
- Canonical minimizer extraction
- Enhanced minimizer scoring functions