# Uniqueness Score Integration Guide

## Overview
This guide helps integrate the uniqueness scoring feature into the existing GPU Kraken database builder to reduce false positives.

## Files to Add
1. `features/uniqueness_score_implementation.cu` - Core uniqueness calculation
2. `features/enhanced_minimizer_flags.h` - Extended feature flag definitions  
3. `features/uniqueness_integration_patch.cu` - Integration examples

## Integration Steps

### Step 1: Update Headers
Add to `core/gpu_database_builder_core.h`:
```cpp
#include "../features/enhanced_minimizer_flags.h"
```

### Step 2: Update Configuration
In `DatabaseBuildConfig` struct (gpu_database_builder_core.h), add:
```cpp
// NEW: Uniqueness scoring configuration
bool enable_uniqueness_scoring = true;
bool enable_uniqueness_filtering = false;     
float uniqueness_threshold = 0.3f;
bool filter_extremely_common = true;
```

### Step 3: Update Statistics  
In `EnhancedBuildStats` struct (gpu_kraken_types.h), add:
```cpp
// NEW: Uniqueness statistics
size_t minimizers_with_uniqueness_scores = 0;
size_t minimizers_filtered_by_uniqueness = 0;
size_t unique_minimizers_count = 0;
size_t rare_minimizers_count = 0;
size_t reliable_minimizers_count = 0;
double average_uniqueness_score = 0.0;
```

### Step 4: Update Core Processing
In `GPUKrakenDatabaseBuilder::process_accumulated_sequences()`:

**Find this existing code:**
```cpp
if (feature_extractor_) {
    if (!feature_extractor_->process_first_pass(...)) {
        // error handling
    }
}
```

**Replace with:**
```cpp
if (feature_extractor_) {
    // Include uniqueness calculation
    if (!integrate_uniqueness_with_feature_extractor(
            feature_extractor_.get(),
            d_minimizer_hits, 
            total_hits_extracted, 
            genome_taxon_ids)) {
        std::cerr << "Uniqueness integration failed" << std::endl;
        return false;
    }
}
```

### Step 5: Add Function Declarations
In `core/gpu_database_builder_core.h`, add:
```cpp
// Forward declarations
bool integrate_uniqueness_with_feature_extractor(
    MinimizerFeatureExtractor* feature_extractor,
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids);

bool compute_and_encode_uniqueness_scores(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& genome_taxon_ids,
    float total_genomes_processed);
```

### Step 6: Include Implementation
At the top of `core/gpu_database_builder_core.cu`, add:
```cpp
#include "../features/uniqueness_score_implementation.cu"
```

## Testing the Integration

### Compile Test
```bash
cd runtime/kernels/profiler/enhanced_k2like/
make clean && make
```

### Expected Output
When running database build, you should see:
```
Computing uniqueness scores for X minimizer hits...
Found Y unique minimizers
Uniqueness Score Distribution:
  Category 0 (Extremely Common): 1234 minimizers
  Category 7 (Singleton-like): 5678 minimizers
✓ Uniqueness scores computed and encoded successfully
```

### Verification
Check that `feature_flags` now contain uniqueness information:
- Bits 8-10 should contain uniqueness categories (0-7)
- Bit 11 should be set for highly unique minimizers
- Bit 13 should be set for reliable minimizers

## Configuration Examples

### High Precision (Fewer False Positives)
```cpp
config.enable_uniqueness_scoring = true;
config.enable_uniqueness_filtering = true;
config.uniqueness_threshold = 0.7f;  // Only keep highly unique
config.filter_extremely_common = true;
```

### Balanced (Default)
```cpp
config.enable_uniqueness_scoring = true;
config.enable_uniqueness_filtering = false;  // Score but don't filter
config.uniqueness_threshold = 0.3f;
```

## Troubleshooting

### Common Issues
1. **Compilation errors**: Check include paths and CUDA compute capability
2. **No uniqueness output**: Verify `enable_uniqueness_scoring = true`
3. **GPU memory errors**: May need to increase memory allocation limits

### Debug Mode
Enable debug output:
```cpp
config.enable_debug_mode = true;
```

This will print detailed uniqueness statistics and feature flag analysis.

## Next Features to Add
1. Co-occurrence scoring (bits 14-16)
2. Contamination detection integration (bit 7)
3. Position clustering detection (bit 17)
