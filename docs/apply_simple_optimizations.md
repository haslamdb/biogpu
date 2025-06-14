# Simple Optimizations for minimizer_extraction.cu

## Step 1: Add pinned memory allocation
- Pre-allocate pinned memory in constructor for typical batch size
- Use cudaMallocHost for h_sequences, h_offsets, h_lengths, h_minimizers, h_counts

## Step 2: Use async memory transfers
- Replace cudaMemcpy with cudaMemcpyAsync
- Add CUDA streams for overlapping transfers

## Step 3: Keep the same kernel launch parameters
- Don't change block size from 256 to 128 yet
- Keep same kernel code

## Key: Test after each step to ensure it still works!

The issue with our first attempt was likely:
1. Changed too many things at once
2. Possible issue with reallocation tracking
3. Might have broken the device memory allocation logic