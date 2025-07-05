# AMR Pipeline Migration Guide

This guide explains how to migrate from the original AMR detection pipeline to the refactored version that uses shared components.

## Key Changes

### 1. **Inheritance from PipelineBase**
The new `AMRDetectionPipeline` inherits from `PipelineBase`, providing:
- Standardized initialization and cleanup
- Built-in profiling support
- Unified configuration management
- Shared GPU resource management

### 2. **Shared I/O Components**
- **Old**: Custom FASTQ reading in the pipeline
- **New**: Uses `StreamingFastqReader` from common components
- Benefits: Async reading, gzip support, paired-end handling

### 3. **GPU Memory Management**
- **Old**: Manual GPU memory allocation/deallocation
- **New**: Uses `GPUSequenceBuffer` with pinned memory
- Benefits: Faster transfers, automatic resizing, double-buffering support

### 4. **Configuration System**
- **Old**: Hardcoded `AMRDetectionConfig` struct
- **New**: Uses `GenesConfig` from unified configuration
- Benefits: JSON/YAML config files, validation, runtime overrides

## Migration Steps

### Step 1: Update Headers

Replace:
```cpp
#include "amr_detection_pipeline.h"
```

With:
```cpp
#include "amr_detection_pipeline_v2.h"
#include "../../common/config/unified_config.h"
```

### Step 2: Update Configuration

Old code:
```cpp
AMRDetectionConfig config;
config.bloom_filter_size = 1ULL << 30;
config.minimizer_k = 15;
config.min_identity = 0.85f;

AMRDetectionPipeline pipeline(config);
```

New code:
```cpp
// Load from file
auto config = ConfigLoader::loadFromFile("config.json");

// Or use defaults and modify
UnifiedConfig config = ConfigLoader::getDefault();
config.genes_config.minimizer_k = 15;
config.genes_config.min_gene_identity = 0.85f;

AMRDetectionPipeline pipeline(config.genes_config);
```

### Step 3: Update Initialization

Old code:
```cpp
pipeline.initialize(amr_db_path);
```

New code:
```cpp
// Database paths are now in configuration
pipeline.initialize();  // Uses paths from config
```

### Step 4: Update Batch Processing

Old code:
```cpp
std::vector<std::string> reads = readFastq(filename);
std::vector<std::string> read_ids = extractIds(reads);
pipeline.processBatch(reads, read_ids);
```

New code:
```cpp
// Create reader and GPU buffer
StreamingFastqReader reader(batch_size);
reader.open(filename);

GPUSequenceBuffer gpu_buffer(max_sequences, max_bases);
gpu_buffer.allocate();

// Process batches
while (reader.hasNext()) {
    auto batch = reader.getNextBatch();
    gpu_buffer.transferBatch(*batch);
    pipeline.processBatch(&gpu_buffer);
}
```

### Step 5: Update Results Handling

Old code:
```cpp
std::vector<AMRHit> hits = pipeline.getAMRHits();
pipeline.writeResults(output_prefix);
```

New code:
```cpp
auto results = pipeline.getResults();
auto* amr_results = dynamic_cast<AMRResults*>(results.get());

// Write reports
amr_results->writeReport("output.tsv");
amr_results->writeClinicalReport("clinical_report.txt");

// Access high-confidence hits
for (const auto& hit : amr_results->high_confidence_hits) {
    // Process hit
}
```

## Configuration File Example

Create a `config.json` file:

```json
{
    "genes": {
        "amr_cds_path": "AMR_CDS.fa",
        "amr_prot_path": "AMRProt.fa",
        "minimizer_k": 15,
        "minimizer_w": 10,
        "min_gene_coverage": 0.95,
        "min_gene_identity": 0.95,
        "enable_profiling": true
    },
    "pipeline": {
        "batch_size": 100000,
        "output_dir": "output"
    }
}
```

## Performance Considerations

1. **Batch Size**: The new pipeline automatically handles batching through `StreamingFastqReader`
2. **GPU Streams**: Multiple CUDA streams are managed by `PipelineBase`
3. **Memory Usage**: `GPUSequenceBuffer` manages pinned memory for optimal transfers
4. **Profiling**: Enable with `config.genes_config.enable_profiling = true`

## Backward Compatibility

To maintain backward compatibility during migration:

1. Keep the old pipeline files until migration is complete
2. Create wrapper functions that translate old config to new format
3. Test thoroughly with your existing test datasets

## Benefits of Migration

1. **Unified Architecture**: Share components with resistance and profiler pipelines
2. **Better Performance**: Async I/O and optimized GPU transfers
3. **Easier Maintenance**: Less code duplication
4. **Configuration Management**: JSON-based configs with validation
5. **Built-in Profiling**: Performance metrics without custom code
6. **Future-proof**: Easy to add new features to shared components

## Troubleshooting

### Common Issues:

1. **Missing includes**: Make sure to add include paths for common directories
2. **Link errors**: Update your Makefile to include common object files
3. **Config validation fails**: Check database paths and parameter ranges
4. **GPU memory errors**: Adjust batch sizes in configuration

### Debug Tips:

1. Enable verbose mode: `config.genes_config.verbose = true`
2. Check profiling output for performance bottlenecks
3. Validate configuration before running: `ConfigLoader::validate()`

## Example Migration

See `test_amr_pipeline_v2.cpp` for a complete example of using the new pipeline.