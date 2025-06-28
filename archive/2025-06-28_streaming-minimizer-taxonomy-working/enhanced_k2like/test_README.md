# Enhanced Database Build Test Suite

This test suite validates the ML-enhanced 24-byte minimizer structure implementation.

## Test Coverage

### 1. Feature Extraction Test
- Validates GC content calculation
- Tests sequence complexity scoring
- Verifies position bias detection
- Ensures features are correctly encoded in 16-bit flags

### 2. Contamination Detection Test
- Tests detection of synthetic contamination sequences
- Validates adapter sequence detection
- Checks for low-complexity region marking
- Measures false positive rate on real sequences

### 3. Database Integrity Test
- Verifies 24-byte structure serialization
- Checks database version (should be 3)
- Validates all ML-specific files are created:
  - `ml_weights.bin` - ML confidence scores
  - `feature_statistics.json` - Feature distributions
  - `contamination_markers.bin` - Contamination list
- Compares database sizes (20-byte vs 24-byte)

### 4. Performance Benchmark
- Measures overhead of ML feature extraction
- Compares build times with and without ML features
- Reports performance impact percentage

## Running the Tests

### Build the test
```bash
make -f Makefile.test
```

### Run full test suite
```bash
make -f Makefile.test run-test
```

### Run quick test (fewer genomes)
```bash
make -f Makefile.test quick-test
```

### Custom test run
```bash
./test_enhanced_database_build \
  --genome-dir data/type_strain_reference_genomes \
  --output-dir my_test_output \
  --max-genomes 10
```

## Expected Output

The test suite will output:
1. Individual test results (PASSED/FAILED)
2. Performance metrics
3. Feature distribution statistics
4. Database size comparison
5. Overall summary

## Test Data

The test uses reference genomes from `data/type_strain_reference_genomes/`. 
It also creates synthetic contamination sequences to test contamination detection.

## Success Criteria

All tests should pass with:
- Feature extraction correctly identifying sequence properties
- Contamination detection rate > 50% for synthetic contaminants
- Database integrity verified with version 3
- Performance overhead < 50% (typically 20-30%)

## Troubleshooting

If tests fail:
1. Ensure CUDA is properly installed
2. Check that reference genomes exist in the data directory
3. Verify sufficient GPU memory is available
4. Check CUDA compute capability (requires sm_70 or higher)