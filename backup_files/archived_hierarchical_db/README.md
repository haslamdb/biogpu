# Archived Hierarchical Database Implementation

This directory contains the archived hierarchical database implementation for the BioGPU profiler.

## Overview

The hierarchical database was designed to handle very large k-mer databases (100s of GB to TBs) by:
- Dividing the database into tiers that could be loaded/unloaded from GPU memory
- Using LRU caching to keep frequently accessed tiers in GPU memory
- Supporting both frequency-sorted and hash-sorted storage

## Files

- `hierarchical_gpu_database.h` - Header-only implementation of the hierarchical database
- `hierarchical_profiler_pipeline.cu` - Profiler pipeline using hierarchical database
- `build_hierarchical_db.cpp` - Tool to build hierarchical databases from k-mer lists
- `test_hierarchical_db.cpp` - Test program for hierarchical database functionality
- `hierarchical_gpu_database_backup.cu` - Earlier implementation attempt

## Status

This implementation was archived because:
1. For databases that fit in CPU RAM (< 100GB), streaming from CPU memory is simpler and faster
2. The disk I/O overhead of loading tiers was significant
3. Complexity of tier management vs. benefits for moderate-sized databases

## Known Issues

1. **Zero hits problem**: When database is sorted by frequency (default), binary search doesn't work. Linear search is too slow.
2. **Floating point exception**: Crash during initialization, possibly in statistics calculation
3. **Performance**: Tier loading from disk creates significant latency

## Future Improvements

If revisiting this approach:
1. Implement CPU RAM as L2 cache between disk and GPU
2. Add async tier loading to hide I/O latency
3. Support both sorted-by-hash and sorted-by-frequency lookups efficiently
4. Implement better tier partitioning strategies

## Usage (when it was working)

```bash
# Build hierarchical database
./build_hierarchical_db kmers.txt output_db --tier-size 512 --sort-hash

# Run profiler
./hierarchical_profiler_pipeline output_db reads.fastq --memory 10
```

## Database Format

The hierarchical database consists of:
- `manifest.json` - Metadata about tiers
- `tier_N.bin` - Binary files containing k-mer entries
- Each tier can be loaded independently based on memory constraints