# Profiler Pipeline Status

## Overview

The profiler pipeline is a GPU-accelerated metagenomic classification system similar to Kraken2. The project has evolved through multiple iterations, with the current focus on building an enhanced ML-enabled database that attaches extensive metadata to k-mers/minimizers.

## Current Architecture

### Directory Structure
```
runtime/kernels/profiler/
├── enhanced_k2like/          # Current active development
│   ├── core/                 # Core database builder
│   ├── features/             # Feature extraction and ML components
│   ├── gpu/                  # GPU kernels
│   ├── memory/               # Memory management
│   ├── processing/           # Genome processing
│   └── taxonomy/             # Taxonomy handling
├── k2like_from_scratch/      # Alternative implementation
├── deprecated/               # Previous iterations
├── tools/                    # Utility tools
└── minimized/               # Simplified versions
```

## Development History

### Phase 1: Basic Kraken2-like Implementation
- Successfully implemented GPU-accelerated minimizer extraction
- Basic database building and classification working
- Paired-end read support

### Phase 2: Enhanced ML Features (Current)
- Added ML confidence scoring
- Implemented feature extraction (GC content, complexity, position bias)
- Enhanced minimizer structure from 20 to 24 bytes
- Added contamination detection capabilities

## Current Status

### What's Working
1. **Enhanced Minimizer Structure** (24 bytes):
   ```cpp
   struct GPUMinimizerHit {
       uint64_t minimizer_hash;  // 8 bytes
       uint32_t genome_id;       // 4 bytes
       uint32_t position;        // 4 bytes
       uint16_t strand;          // 2 bytes (includes classification flags)
       uint16_t taxon_id;        // 2 bytes
       uint16_t ml_weight;       // 2 bytes (ML confidence 0.0-1.0)
       uint32_t feature_flags;   // 4 bytes (GC, complexity, bias, contamination)
   };
   ```

2. **Feature Encoding**:
   - GC content categories (8 levels)
   - Sequence complexity scores (8 levels)
   - Position bias detection
   - Contamination risk flagging

3. **ML Weight System**:
   - Successfully encodes confidence scores (0.0-1.0) in 16-bit values
   - Tested and validated in TEST_RESULTS.md

### Memory Issues and Blockers

#### 1. **Metadata Explosion Problem**
The addition of extensive metadata to each minimizer causes memory issues:
- **Original structure**: 20 bytes per minimizer
- **Enhanced structure**: 24 bytes per minimizer (20% increase)
- **MinimizerStatistics tracking**: Additional ~200+ bytes per unique minimizer

```cpp
struct MinimizerStatistics {
    uint32_t occurrence_count;
    std::vector<uint32_t> taxon_occurrences;         // Variable size
    std::unordered_set<uint64_t> neighbor_minimizers; // Variable size
    std::unordered_map<uint64_t, uint32_t> cooccurring_minimizers;
    // ... more fields
};
```

#### 2. **Memory Calculation Overflows**
The deprecated code shows attempts to handle memory overflows:
- Sequence memory calculations could overflow with large databases
- Minimizer count estimates were exceeding size_t limits
- Total memory requirements exceeding GPU capacity

#### 3. **Statistical Tracking Overhead (PRIMARY MEMORY BOTTLENECK)**
The `minimizer_stats_` map in the database builder is the most memory-intensive component:
- Stores detailed statistics for EVERY unique minimizer
- Can consume hundreds of gigabytes of RAM for large databases
- No efficient pruning mechanism implemented

**Memory Multiplication Effect**:
- Microbial databases have **billions of unique minimizers**
- Each minimizer stores ~200+ bytes of statistics (unbounded due to vectors/maps)
- Total memory: billions × 200+ bytes = **hundreds of gigabytes**
- The statistics tracking dwarfs the actual minimizer data by 10x or more

**Why MinimizerStatistics is the Primary Culprit**:
```cpp
std::unordered_map<uint64_t, MinimizerStatistics> minimizer_stats_;
// For each unique minimizer hash, stores:
// - Variable-length vector of taxon occurrences
// - Variable-length set of ALL neighbor minimizers ever seen
// - Variable-length map of co-occurring minimizers with counts
// - Additional metadata fields
```

The co-occurrence tracking is particularly problematic:
- For each minimizer, tracks ALL its neighbors across ALL genomes
- Creates O(n²) memory growth in worst case
- Example: If a minimizer appears 1000 times with different neighbors, it stores 1000+ entries

## Technical Challenges

### 1. **Scale Mismatch**
- Microbial databases contain billions of minimizers
- Each minimizer now carries 24 bytes + statistical tracking
- Memory requirements grow exponentially with database size

### 2. **GPU Memory Constraints**
- Most GPUs have 8-48GB memory
- Database + working memory must fit within this limit
- Current design requires too much metadata per minimizer

### 3. **Co-occurrence Tracking**
- Tracking minimizer neighbors for co-occurrence scoring
- Each minimizer potentially tracks hundreds of neighbors
- Memory usage becomes O(n²) in worst case

## Attempted Solutions

### 1. **Memory Pool Implementation**
- Created GPUMemoryPool for efficient allocation
- Helps with fragmentation but not total size

### 2. **Auto-scaling Configuration**
- Attempts to dynamically adjust batch sizes
- Still hits hard limits with metadata overhead

### 3. **Safe Memory Checking**
- Added overflow protection in calculations
- Prevents crashes but doesn't solve fundamental issue

## Recommendations

### Short-term Fixes
1. **Reduce Metadata Granularity**:
   - Track statistics only for subset of minimizers
   - Use sampling instead of exhaustive tracking
   - Implement pruning for rare minimizers

2. **Hierarchical Storage**:
   - Keep only essential data in GPU memory
   - Move statistics to disk-based storage
   - Load on-demand for analysis

3. **Feature Quantization**:
   - Reduce feature_flags from 32 to 16 bits
   - Use fewer GC/complexity categories
   - Pack more efficiently

### Long-term Solutions
1. **Streaming Architecture**:
   - Process database in chunks
   - Don't load entire database into memory
   - Stream minimizers through classification

2. **Tiered Metadata System**:
   - Level 1: Essential classification data (in GPU)
   - Level 2: Common statistics (in RAM)
   - Level 3: Detailed analytics (on disk)

3. **Approximate Methods**:
   - Use sketching for co-occurrence tracking
   - Probabilistic data structures for statistics
   - Trade accuracy for memory efficiency

## Current Blockers

1. **MinimizerStatistics Memory Usage**
   - Unbounded growth with database size
   - No compression or approximation
   - Needs fundamental redesign

2. **Co-occurrence Matrix Storage**
   - Attempting to store all pairwise relationships
   - Quadratic memory growth
   - Should use sparse representations

3. **Feature Richness vs. Scalability**
   - Current design prioritizes features over scale
   - Need to identify essential vs. nice-to-have metadata
   - Consider feature importance analysis

## Path Forward

### Phase 1: Stabilization (1-2 weeks)
1. Implement minimizer sampling (e.g., track stats for 1% of minimizers)
2. Add memory usage estimation before building
3. Create configuration presets for different memory constraints

### Phase 2: Optimization (2-4 weeks)
1. Redesign MinimizerStatistics to use fixed-size representation
2. Implement sparse co-occurrence tracking
3. Add streaming support for large databases

### Phase 3: Production Ready (1-2 months)
1. Benchmark against standard databases
2. Validate classification accuracy with reduced metadata
3. Create memory-efficient production configuration

## Testing Strategy

1. **Memory Scaling Tests**:
   - Start with small test genomes
   - Gradually increase to full databases
   - Monitor memory usage at each scale

2. **Feature Ablation**:
   - Test classification accuracy with different feature subsets
   - Identify minimal feature set for good performance

3. **Performance Benchmarks**:
   - Compare to original Kraken2
   - Measure speed/memory/accuracy tradeoffs

## Conclusion

The profiler pipeline has made significant progress in implementing ML-enhanced features, but the addition of extensive metadata has created memory scalability issues. The core challenge is balancing feature richness with practical memory constraints. The recommended path forward involves reducing metadata granularity, implementing hierarchical storage, and moving to a streaming architecture for large-scale databases.