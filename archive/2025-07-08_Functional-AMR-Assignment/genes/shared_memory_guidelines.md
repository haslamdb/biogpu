# GPU Shared Memory Calculation Guidelines for BioGPU

Based on analysis of the existing BioGPU kernel implementations, here are the standard patterns and guidelines for shared memory calculations:

## 1. Static Shared Memory Allocation

For fixed-size shared memory allocations, declare directly in the kernel:

```cuda
__shared__ int shared_stats[4];  // From bloom_filter.cu
```

## 2. Dynamic Shared Memory Allocation

### Pattern A: Simple Dynamic Allocation
```cuda
extern __shared__ int shared_mem[];
int* scores = shared_mem;  // Dynamic shared memory
```

### Pattern B: Multi-purpose Dynamic Allocation with Offset Calculation
From `amr_detection_kernels.cu`:
```cuda
extern __shared__ char shared_translations[];
const int max_protein_len = (length / 3) + 1;
char* translations = &shared_translations[threadIdx.x * max_protein_len * 6];
```

## 3. Shared Memory Size Calculation Patterns

### For Translated Search (from amr_detection_kernels_wrapper.cu):
```cuda
size_t sharedMemSize = blockSize.x * 6 * 67 * sizeof(char); 
// Formula: threads * frames * max_protein_length * element_size
// 6 frames for 6-frame translation
// 67 is max protein length per read segment
```

### For Local Alignment (from resistance_detection_gpu.cu):
```cuda
extern __shared__ float dp_matrix[];
// Size would be: (query_len + 1) * (target_len + 1) * sizeof(float)
```

## 4. Kernel Launch with Dynamic Shared Memory

Standard pattern for launching kernels with dynamic shared memory:
```cuda
kernel_name<<<gridSize, blockSize, sharedMemSize>>>(parameters...);
```

## 5. Best Practices Found in BioGPU

1. **Calculate shared memory per thread block**, not per thread
2. **Account for alignment requirements** - CUDA requires proper alignment
3. **Consider bank conflicts** - organize data to minimize conflicts
4. **Use shared memory for**:
   - Frequently accessed read-only data
   - Inter-thread communication within a block
   - Temporary working buffers (DP matrices, translations)
   
## 6. Common Shared Memory Uses in BioGPU

1. **Statistics Collection**: Small fixed arrays for block-level statistics
2. **Sequence Processing**: Temporary buffers for translations, reverse complements
3. **Dynamic Programming**: Matrices for local/global alignment
4. **K-mer Processing**: Temporary storage for k-mer generation and hashing

## 7. Memory Size Limits

- Maximum shared memory per block: 48KB (compute capability 3.x-6.x) or 96KB (7.x+)
- Configure with `cudaFuncSetAttribute` if needed:
```cuda
cudaFuncSetAttribute(kernel_name, cudaFuncAttributeMaxDynamicSharedMemorySize, size);
```

## 8. Example Calculation Template

```cuda
// Calculate shared memory requirements
size_t sharedMemPerThread = /* data per thread */;
size_t sharedMemPerBlock = blockSize.x * sharedMemPerThread;

// Add padding for alignment if needed
sharedMemPerBlock = ((sharedMemPerBlock + 15) / 16) * 16; // 16-byte alignment

// Launch kernel
kernel<<<gridSize, blockSize, sharedMemPerBlock>>>(...);
```