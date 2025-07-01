# Taxonomy Architecture Overview

## Components

### Compact Taxonomy (tools/compact_gpu_taxonomy.h)
This is a **GPU-optimized data structure** for taxonomy:
- **Purpose**: Efficient taxonomy storage and lookup on GPU during classification
- **Focus**: Memory efficiency and GPU performance
- **Usage**: Real-time classification operations
- **Key features**:
  - Compact memory representation
  - GPU-friendly data structures
  - Fast parallel lookups
  - Optimized for read-only access during classification

### Taxonomy Module (enhanced_k2like/taxonomy/)
This is the **host-side taxonomy processing** module:
- **Purpose**: Database building, taxonomy loading, and complex operations
- **Focus**: Flexibility and full-featured taxonomy operations
- **Usage**: Database construction phase
- **Key features**:
  - Load taxonomy from various sources (NCBI files, TSV, etc.)
  - Complex LCA computations
  - Phylogenetic analysis
  - Species tracking
  - Both simple and enhanced processors

## How They Work Together

```
Database Building Phase:
1. taxonomy/ module loads NCBI taxonomy files
2. Processes and validates taxonomy data
3. Computes LCAs and phylogenetic information
4. Passes data to compact_gpu_taxonomy for GPU optimization
5. compact_gpu_taxonomy creates efficient GPU representation

Classification Phase:
1. compact_gpu_taxonomy loads pre-built compact format
2. Provides fast GPU lookups during classification
3. Minimal memory footprint on GPU
```

## Key Architectural Decisions

### Separation of Concerns
- **Building complexity vs runtime efficiency**: The taxonomy module handles all the complex processing during database construction, while compact_gpu_taxonomy focuses solely on efficient runtime access
- **Two-phase approach**: Rich processing capabilities during build phase, lean and fast structure during classification phase
- **Flexibility**: Can use SimpleTaxonomyProcessor without GPU requirements, or EnhancedNCBITaxonomyProcessor with full GPU integration

### Design Analogy
Think of it like:
- `taxonomy/` = Full-featured kitchen for preparing the meal (all tools and ingredients available)
- `compact_gpu_taxonomy` = Optimized lunchbox for eating on the go (only what you need, packaged efficiently)

### Benefits
1. **Performance**: Classification doesn't pay for features it doesn't need
2. **Modularity**: Can update/improve either component independently
3. **Memory efficiency**: GPU memory is precious, compact format maximizes capacity
4. **Development flexibility**: Complex algorithms stay on host side where debugging is easier

## Usage Patterns

### During Database Building
```cpp
// Use the full-featured taxonomy module
EnhancedNCBITaxonomyProcessor taxonomy;
taxonomy.load_ncbi_taxonomy(nodes_file, names_file);

// Process and analyze
uint32_t lca = taxonomy.compute_lca_of_species(species_list);
uint8_t spread = taxonomy.calculate_phylogenetic_spread(species_list, lca);

// Convert to compact format for GPU
auto compact = taxonomy.get_compact_taxonomy();
compact->save_compact_format(output_file);
```

### During Classification
```cpp
// Load only the compact GPU taxonomy
CompactGPUTaxonomy gpu_taxonomy(true);  // enable caching
gpu_taxonomy.load_compact_taxonomy(compact_file);

// Fast GPU lookups during classification
// Minimal memory footprint, maximum performance
```

## File Organization
- `tools/compact_gpu_taxonomy.h` - GPU-optimized taxonomy structure
- `enhanced_k2like/taxonomy/taxonomy_processor.h` - Full taxonomy processing
- `enhanced_k2like/taxonomy/taxonomy_processor.cu` - Implementation with all algorithms

This separation allows for complex taxonomy operations during database building while maintaining high performance during classification.