# BioGPU Taxonomy System Documentation

## Overview

The BioGPU project uses a **two-phase taxonomy architecture** optimized for different stages of the pipeline:

1. **Build Phase**: `EnhancedNCBITaxonomyProcessor` - Full-featured taxonomy processing during database construction
2. **Runtime Phase**: `CompactGPUTaxonomy` - GPU-optimized compact format for classification

## Architecture

```
Database Building Phase:
┌─────────────────────────────────┐
│ EnhancedNCBITaxonomyProcessor   │
├─────────────────────────────────┤
│ • Load NCBI taxonomy files      │
│ • Compute complex LCAs          │
│ • Calculate phylogenetic spread │
│ • Species tracking              │
│ • Generate enhanced candidates  │
└─────────────────────────────────┘
                ↓
         Builds and exports
                ↓
┌─────────────────────────────────┐
│    CompactGPUTaxonomy           │
├─────────────────────────────────┤
│ • GPU-optimized hash table      │
│ • Minimal memory footprint      │
│ • Fast parallel lookups         │
│ • Cached distance computations  │
└─────────────────────────────────┘
                ↓
         Used during
                ↓
    Runtime Classification (GPU)
```

## Component Details

### 1. EnhancedNCBITaxonomyProcessor (Build Time)

Located in: `enhanced_k2like/taxonomy/`

**Features:**
- Full NCBI taxonomy loading from nodes.dmp/names.dmp
- Complex LCA computations with multiple species
- Phylogenetic distance and spread calculations
- Species diversity metrics
- Integration with PhylogeneticLCACandidate structures

**Usage during database building:**
```cpp
EnhancedNCBITaxonomyProcessor taxonomy;
taxonomy.load_ncbi_taxonomy("nodes.dmp", "names.dmp");

// Complex operations for database construction
auto lca = taxonomy.compute_lca_of_species(species_list);
auto spread = taxonomy.calculate_phylogenetic_spread(species_list, lca);
```

### 2. CompactGPUTaxonomy (Runtime)

Located in: `tools/`

**Features:**
- Compact binary format (~70% smaller than raw taxonomy)
- GPU-accelerated hash table lookups
- Batch operations for classification
- Optional distance caching
- Minimal CPU-GPU memory transfers

**Usage during classification:**
```cpp
CompactGPUTaxonomy gpu_taxonomy(true);  // Enable caching
gpu_taxonomy.load_compact_taxonomy("taxonomy.bin");

// Fast GPU lookups during read classification
gpu_taxonomy.lookup_taxons_gpu(query_taxons, parents, depths, ranks);
gpu_taxonomy.compute_lca_batch_gpu(taxon_pairs, lca_results);
```

## Building the Taxonomy

### Step 1: Download NCBI Taxonomy

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
```

### Step 2: Build Compact Format

```bash
./build_compact_taxonomy \
    --nodes nodes.dmp \
    --names names.dmp \
    --output compact_taxonomy.bin \
    --validate
```

### Step 3: Use in Database Building

```cpp
// In your database builder
EnhancedNCBITaxonomyProcessor taxonomy;
taxonomy.load_ncbi_taxonomy("nodes.dmp", "names.dmp");

// ... perform database building operations ...

// Export compact format for runtime
auto* compact = taxonomy.get_compact_taxonomy();
compact->save_compact_taxonomy("database/taxonomy.bin");
```

### Step 4: Use in Classification

```cpp
// In your classifier
CompactGPUTaxonomy gpu_taxonomy(true);
gpu_taxonomy.load_compact_taxonomy("database/taxonomy.bin");

// Ready for GPU-accelerated classification
```

## Performance Characteristics

### Build Time (EnhancedNCBITaxonomyProcessor)
- **Memory**: ~2GB for full NCBI taxonomy
- **Load time**: 10-30 seconds
- **LCA computation**: ~1μs per pair (CPU)
- **Use case**: One-time database construction

### Runtime (CompactGPUTaxonomy)
- **Memory**: ~500MB for compact format
- **Load time**: <1 second
- **GPU lookups**: ~100M lookups/second
- **GPU LCA**: ~50M LCAs/second
- **Use case**: High-throughput classification

## File Formats

### Compact Taxonomy Binary Format (.bin)

```
Header (16 bytes):
  - Magic number: 0x54415843 ("CTAX")
  - Version: uint32_t
  - Max taxon ID: uint32_t
  - Number of nodes: uint32_t

Node data:
  - Array of CompactTaxonNode structures

Name data:
  - Count: uint32_t
  - Array of (taxon_id, name_length, name_string)
```

## API Reference

### CompactGPUTaxonomy Key Methods

```cpp
// Build from NCBI files
bool build_from_ncbi_files(const string& nodes_dmp, const string& names_dmp);

// Save/load compact format
bool save_compact_taxonomy(const string& output_path);
bool load_compact_taxonomy(const string& input_path);

// GPU operations
bool lookup_taxons_gpu(
    const vector<uint32_t>& query_taxons,
    vector<uint32_t>& parent_results,
    vector<uint8_t>& depth_results,
    vector<uint8_t>& rank_results
);

bool compute_lca_batch_gpu(
    const vector<pair<uint32_t, uint32_t>>& taxon_pairs,
    vector<uint32_t>& lca_results
);
```

### EnhancedNCBITaxonomyProcessor Key Methods

```cpp
// Load taxonomy
bool load_ncbi_taxonomy(const string& nodes_dmp, const string& names_dmp);

// Compute LCAs
uint32_t compute_lca_of_species(const vector<uint32_t>& species_list);

// Phylogenetic metrics
uint8_t calculate_phylogenetic_spread(
    const vector<uint32_t>& species_list, 
    uint32_t lca
);

// Get compact taxonomy for export
CompactGPUTaxonomy* get_compact_taxonomy();
```

## Integration with Fluoroquinolone Resistance Pipeline

The taxonomy system integrates with your resistance detection pipeline:

1. **During Database Building**:
   ```cpp
   // Associate resistance genes with taxa
   EnhancedNCBITaxonomyProcessor taxonomy;
   for (const auto& genome : genomes) {
       uint32_t taxon_id = genome.taxon_id;
       // Extract gyrA, gyrB, parC, parE from genome
       // Store associations with taxonomy
   }
   ```

2. **During Classification**:
   ```cpp
   // Fast lookup of taxon for resistance interpretation
   CompactGPUTaxonomy gpu_taxonomy;
   // When resistance mutation found in read
   uint32_t source_taxon = gpu_taxonomy.lookup_taxon(hit.genome_id);
   ```

## Troubleshooting

### Out of GPU Memory
- Reduce batch sizes in classification
- Disable distance caching: `CompactGPUTaxonomy(false)`
- Use streaming mode for large databases

### Taxonomy Not Found
- Ensure taxon IDs match between database and taxonomy
- Check that all parent relationships are valid
- Verify root node (ID=1) exists

### Performance Issues
- Enable distance caching for repeated LCA queries
- Increase batch sizes for GPU operations
- Profile with `nvprof` to identify bottlenecks

## Building and Testing

```bash
# Build all components
make all

# Run tests
make test

# Build and validate taxonomy
./build_compact_taxonomy --nodes nodes.dmp --names names.dmp \
    --output test.bin --validate

# Run integration example
./taxonomy_integration_example nodes.dmp names.dmp
```

## Future Enhancements

1. **Incremental Updates**: Add new taxa without rebuilding
2. **Custom Taxonomies**: Support non-NCBI taxonomic systems
3. **Distributed Processing**: Multi-GPU taxonomy operations
4. **Compression**: Further reduce memory footprint
5. **ML Integration**: Taxonomic distance features for resistance prediction