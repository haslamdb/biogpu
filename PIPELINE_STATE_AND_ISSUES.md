# Enhanced Kraken Pipeline - Current State and Issues

## Current Error
```
malloc(): invalid size (unsorted)
Aborted (core dumped)
```

This error occurs immediately after:
- "Initializing GPU Kraken database builder (High Capacity)..."
- "Parameters: k=35, ell=31, spaces=7"
- "Initial safe settings: Minimizers per batch: 1000000, Sequences per batch: 10"
- "Enabled auto memory scaling (using 80% of GPU memory)"

## Recent Changes Made
1. Implemented auto memory scaling feature (--auto-memory flag)
2. Added high capacity mode with larger batch sizes
3. Modified GPUKrakenDatabaseBuilder constructor to support dynamic memory allocation
4. Added memory profiling and optimization logic

## Key Issues to Fix

### 1. Malloc Error in Constructor
- The malloc error happens during GPUKrakenDatabaseBuilder initialization
- Likely related to incorrect size calculation or uninitialized variable
- Could be in the auto memory scaling logic

### 2. Memory Calculation Issues
- Auto memory scaling tries to use 80% of GPU memory
- May be calculating incorrect sizes for allocations
- Need to verify all size_t calculations are valid

### 3. Potential Problem Areas
- Constructor parameter initialization order
- Size calculations for minimizer/sequence batches
- Memory alignment requirements not being met
- Possible integer overflow in size calculations

## Current Implementation Status
- Basic pipeline structure is complete
- GPU kernels for k-mer processing implemented
- Taxonomy loading and management working
- Database building framework in place
- Auto memory optimization partially implemented

## What Needs Investigation
1. Check all malloc/cudaMalloc calls in the constructor
2. Verify size_t calculations don't overflow
3. Ensure all variables are initialized before use
4. Check if auto memory scaling calculations are correct
5. Verify batch size calculations are reasonable

## System Info
- GPU: NVIDIA TITAN Xp (11 GB)
- Running on Linux
- Using CUDA for GPU acceleration