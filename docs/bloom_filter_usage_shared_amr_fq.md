# Combined Minimizer Bloom Filter for AMR and FQ Resistance Detection

## Overview

We have developed a unified approach for building a **combined minimizer-based bloom filter** that includes both fluoroquinolone (FQ) resistance genes and comprehensive antimicrobial resistance (AMR) genes. This approach uses **minimizers** instead of all k-mers, providing significant space savings while maintaining sensitivity.

## Key Tool: `build_minimizer_bloom_filter`

This tool builds a minimizer-based bloom filter from nucleotide sequences, combining multiple FASTA sources into a single efficient pre-screening filter.

### Basic Usage

```bash
./build_minimizer_bloom_filter <output_bloom_filter.bin> [input1] [input2] ...
```

### Default Behavior (Recommended)

When run with just an output filename, the tool automatically uses both default databases:

```bash
./build_minimizer_bloom_filter combined_bloom_filter.bin
```

This command automatically processes:
1. **FQ Resistance Genes**: `../../data/resistance_db/reference_sequences.fasta`
2. **AMR Genes**: `../../data/AMR_CDS.fa`

### Custom Input Files

You can specify custom FASTA files or directories:

```bash
# Single custom source
./build_minimizer_bloom_filter custom_bloom.bin /path/to/sequences.fasta

# Multiple custom sources
./build_minimizer_bloom_filter multi_bloom.bin /path/to/fq_genes.fa /path/to/amr_genes.fa /path/to/more_genes.fa

# Mix of files and directories
./build_minimizer_bloom_filter mixed_bloom.bin /path/to/genes_dir/ /path/to/additional.fasta
```

## Default Database Sources

### 1. FQ Resistance Genes
- **Location**: `data/resistance_db/reference_sequences.fasta`
- **Content**: Reference sequences for fluoroquinolone resistance genes (gyrA, gyrB, parC, parE) from multiple bacterial species
- **Source**: Downloaded from NCBI based on known quinolone resistance mutations
- **Format**: Standard FASTA format with headers containing species and gene information

### 2. AMR Genes Database
- **Location**: `data/AMR_CDS.fa`
- **Content**: Comprehensive antimicrobial resistance gene sequences from the NCBI AMR database
- **Source**: NCBI's curated AMRFinderPlus database
- **Format**: FASTA format with 9,257+ resistance gene sequences covering multiple drug classes

## Technical Configuration

- **K-mer length**: 15 (optimal for bacterial sequences)
- **Window size**: 10 (for minimizer selection)
- **Reverse complements**: Included automatically
- **Minimizer advantage**: ~79% reduction in elements vs all k-mers
- **Bloom filter size**: 64MB (2^26 bits)
- **Hash functions**: 3 (MurmurHash3)
- **False positive rate**: Typically < 0.01%

## Building Process Example

```bash
# Build combined bloom filter with default sources
./build_minimizer_bloom_filter combined_bloom_filter.bin

# Output:
# Processing FQ resistance genes: ../../data/resistance_db/reference_sequences.fasta
# Processed 171 sequences from ../../data/resistance_db/reference_sequences.fasta
# 
# Processing AMR genes: ../../data/AMR_CDS.fa
# Processed 9257 sequences from ../../data/AMR_CDS.fa
# 
# Extracted 1058426 unique minimizers (k=15, w=10, including RC)
# 
# Building bloom filter...
# Bloom Filter Statistics:
#   Size: 8 MB
#   Bits set: 3101886 / 67108864
#   Fill ratio: 4.62217%
#   Estimated FP rate: 0.00987502%
#   Contains RC k-mers: YES
# Successfully saved bloom filter to: combined_bloom_filter.bin
```

## Updating the Bloom Filter

When either database is updated:

1. **Update FQ resistance genes**:
   ```bash
   # Re-download sequences if needed
   python src/python/download_ncbi_20250529.py ...
   
   # Rebuild bloom filter (will use new sequences automatically)
   ./build_minimizer_bloom_filter combined_bloom_filter.bin
   ```

2. **Update AMR genes**:
   ```bash
   # Download new AMR_CDS.fa from NCBI
   # Then rebuild bloom filter
   ./build_minimizer_bloom_filter combined_bloom_filter.bin
   ```

## Integration with Pipelines

The generated bloom filter can be used by both resistance detection pipelines:

1. **Copy to appropriate location**:
   ```bash
   cp combined_bloom_filter.bin data/integrated_clean_db/nucleotide/
   ```

2. **Pipeline will automatically load if present**:
   - FQ resistance pipeline checks for `bloom_filter.bin` in the index directory
   - AMR gene pipeline can be configured to use the same bloom filter

## Advantages of Minimizer-Based Approach

1. **Space Efficiency**: ~1M minimizers vs ~5M k-mers (79% reduction)
2. **Lower False Positive Rate**: 0.01% vs 0.7% with all k-mers
3. **Unified Pre-screening**: Single bloom filter for both pipelines
4. **Fast Updates**: Rebuilding takes seconds, not minutes
5. **Maintains Sensitivity**: Minimizers preserve important sequence features

## Legacy Tools (Still Available)

### `build_bloom_filter`
Builds bloom filter from all k-mers (not minimizers):
```bash
./build_bloom_filter <input.fasta> <output_bloom_filter.bin>
```

### `build_bloom_filter_simple`
Builds from existing k-mer index files:
```bash
./build_bloom_filter_simple <kmer_index.bin> <output_bloom_filter.bin>
```

## Pipeline Integration Tool: `pipeline_tester_bloom`

We have developed an integrated testing tool that combines sample CSV parsing, FASTQ.gz reading, and bloom filter screening in a single executable.

### Usage

```bash
# Without bloom filter (default)
./pipeline_tester_bloom samples.csv

# With bloom filter screening
./pipeline_tester_bloom samples.csv --bloom-filter

# With custom bloom filter file
./pipeline_tester_bloom samples.csv --bloom-filter custom_bloom.bin
```

### Features

- **Batch Processing**: Processes reads in batches of 1000 for efficiency
- **Minimizer Extraction**: Extracts 15-mer minimizers with window size 10
- **GPU Acceleration**: Uses CUDA for bloom filter screening
- **Detailed Statistics**: Reports comprehensive filtering metrics
- **Paired-End Support**: Either R1 or R2 passing causes the read pair to pass

### Performance Analysis

Based on testing with 3.4 million reads:

```
=== Bloom Filter Statistics ===
Total reads processed: 3405007
Reads passed filter: 1169822 (34.3559%)
Total minimizers extracted: 868663516
Minimizers found in bloom filter: 2502102 (0.28804%)
Processing time: 25351.8 ms
==============================
```

### Performance Considerations

**Current Status**: The bloom filter screening takes ~25 seconds for 3.4M reads, with 34% of reads passing the filter.

**Recommendation**: The bloom filter may actually be **slower** than direct processing in some cases:
- 34% pass rate means we're still processing 1/3 of all reads
- 25 seconds overhead may exceed the time saved by filtering
- The next pipeline steps may process all reads faster than bloom filter pre-screening

**When to Use Bloom Filter**:
1. When you expect very low hit rates (<5% of reads containing resistance genes)
2. When downstream processing is computationally expensive (>10 seconds per million reads)
3. When working with extremely large datasets where I/O is the bottleneck

**When to Skip Bloom Filter**:
1. For typical clinical samples where resistance genes are common
2. When the full pipeline processes reads in <30 seconds
3. When you need maximum sensitivity and don't want any false negatives

### Making the Bloom Filter More Stringent

To reduce the pass rate, you can:

1. **Increase min_kmer_threshold** (currently 1):
   ```cpp
   int min_kmers_threshold = 5;  // Require at least 5 minimizers to pass
   ```

2. **Use a smaller bloom filter** to increase false positive rate intentionally
3. **Build bloom filter from a more curated gene set** (fewer reference sequences)

## Best Practices

1. **Benchmark first**: Test your full pipeline with and without bloom filter to determine if it provides a speed benefit
2. **Consider your data**: Clinical samples with enriched resistance genes may not benefit from bloom filtering
3. **Monitor performance**: If bloom filter screening takes longer than downstream processing, disable it
4. **Rebuild after database updates** to ensure coverage of new resistance genes
5. **Use consistent k-mer length** (15) across all tools for compatibility