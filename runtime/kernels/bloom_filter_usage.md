# Bloom Filter Generation Guide

## Overview
We now have C++ tools to generate 15-mer nucleotide bloom filters for both the FQ resistance and AMR gene pipelines. Both use the same k-mer length (15) and include forward and reverse complement k-mers.

## Tools Built

### 1. `build_bloom_filter`
Extracts k-mers directly from FASTA files and builds a bloom filter.

**Usage:**
```bash
./build_bloom_filter <input> <output_bloom_filter.bin>
```

**Examples:**
```bash
# For AMR genes database (from AMR_CDS.fa)
./build_bloom_filter ../../data/AMR_CDS.fa amr_bloom_filter.bin

# For FQ resistance genes (if converted to FASTA)
./build_bloom_filter fq_resistance_genes.fasta fq_bloom_filter.bin
```

### 2. `build_bloom_filter_simple`
Uses existing k-mer index files (kmer_index.bin) to build a bloom filter.

**Usage:**
```bash
./build_bloom_filter_simple <kmer_index.bin> <output_bloom_filter.bin>
```

**Examples:**
```bash
# For FQ resistance genes (using existing index)
./build_bloom_filter_simple ../../data/integrated_clean_db/nucleotide/kmer_index.bin fq_bloom_filter.bin

# For AMR genes (if k-mer index exists)
./build_bloom_filter_simple ../../data/amr_nucleotide_db/kmer_index.bin amr_bloom_filter.bin
```

## Recommended Workflow

### For FQ Resistance Genes:
Since the FQ genes are already indexed during database building:
```bash
./build_bloom_filter_simple ../../data/integrated_clean_db/nucleotide/kmer_index.bin fq_bloom_filter.bin
```

### For AMR Genes:
Since AMR_CDS.fa is a standard FASTA file:
```bash
./build_bloom_filter ../../data/AMR_CDS.fa amr_bloom_filter.bin
```

## Updating Bloom Filters
When databases are updated, simply re-run the appropriate command:

1. **After updating FQ resistance genes:**
   - Re-run the Python database builder to update kmer_index.bin
   - Then: `./build_bloom_filter_simple <new_index> fq_bloom_filter.bin`

2. **After updating AMR_CDS.fa:**
   - Simply run: `./build_bloom_filter ../../data/AMR_CDS.fa amr_bloom_filter.bin`

## Technical Details
- K-mer length: 15 (fixed for both pipelines)
- Includes reverse complement k-mers
- Bloom filter size: 64MB (2^26 bits)
- Hash functions: 3 (MurmurHash3)
- Compatible with existing bloom_filter.cu implementation

## Integration with Pipeline Tester
The generated bloom filters can be used with the pipeline_tester by updating it to load the pre-built bloom filters instead of building them on-the-fly.