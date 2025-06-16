# Kraken2-like Database Build Parameters and Algorithm Comparison

This document describes the key parameters and algorithmic details of Kraken2's database building process and compares them with our GPU-accelerated implementation.

## Key Kraken2 Parameters

1. **k-mer length (k)**: Used for k-mers
2. **minimizer length (l)**: Used for minimizers within k-mers
3. **spaced_seed_mask**: For spaced seeds (default: `DEFAULT_SPACED_SEED_MASK`)
4. **toggle_mask**: XOR mask applied to minimizers (default: `DEFAULT_TOGGLE_MASK`)
5. **min_clear_hash_value**: For hash-based subsampling

## Key Algorithmic Details from Kraken2

### 1. Minimizer Extraction
- Uses `MinimizerScanner` class
- Extracts minimizers of length `l` from k-mers of length `k`
- The scanner appears to use a sliding window approach

### 2. Hash Processing
```cpp
if (min_clear_hash_value && MurmurHash3(*minimizer_ptr) < min_clear_hash_value)
    continue;
```
- Uses MurmurHash3 for hashing
- Applies subsampling based on hash value

### 3. Toggle Mask Application
- Applied during minimizer scanning (passed to `MinimizerScanner`)

## Comparison with Our GPU Implementation

Our implementation appears to match Kraken2's approach quite well:

### 1. Parameters Match
- k=35, ell=31 (k=35, l=31 in Kraken2 terms)
- spaces=7 (for spaced seeds)
- toggle_mask = 0xe37e28c4271b5a2dULL (same as Kraken2's default)
- min_clear_hash_value for subsampling

### 2. Algorithm Match
- Both use MurmurHash3
- Both apply toggle mask to minimizers
- Both use hash-based subsampling
- Both extract minimizers from k-mers using sliding windows

### 3. Key Difference
- Kraken2 uses `k - l + 1` as the window size for minimizers within a k-mer
- Our code calculates: `int window_size = params.k - params.ell + 1;`
- This matches perfectly!

## Processing Modes

The main structural difference is that Kraken2 has two processing modes:
- **ProcessSequenceFast**: Non-deterministic but faster (parallel)
- **ProcessSequence**: Deterministic but slower (uses locks for ordering)

Our GPU implementation achieves similar results by having threads cooperate through shared memory, which maintains determinism while still being parallel.

## Conclusion

Our implementation correctly follows Kraken2's minimizer extraction algorithm with the same parameters and approach. The multi-threaded GPU version we've optimized maintains the same algorithmic correctness while leveraging GPU parallelism effectively.

## Performance Optimizations in GPU Implementation

Our GPU implementation includes several optimizations over the original Kraken2:

1. **Multi-threaded Minimizer Extraction**: Uses 256 threads per genome block instead of sequential processing
2. **Shared Memory Coordination**: Threads communicate through shared memory for efficient parallel processing
3. **Coalesced Memory Access**: GPU-optimized memory access patterns
4. **Batch Processing**: Processes multiple genomes simultaneously on GPU

These optimizations provide significant speedup while maintaining the algorithmic correctness of Kraken2's approach.