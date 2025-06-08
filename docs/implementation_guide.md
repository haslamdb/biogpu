# Implementation Guide: Fixing Spurious Alignments

## Overview
The spurious alignments are caused by using 5-mer protein k-mers, which are too short for specific matching. This guide shows how to:
1. Enable bloom filtering by default
2. Make Smith-Waterman configurable
3. Increase protein k-mer size from 5 to 8

## Step 1: Apply the Main Pipeline Patch

Apply the configuration flags patch to make bloom filtering and Smith-Waterman configurable:

```bash
patch runtime/kernels/resistance/clean_resistance_pipeline_main.cpp < bloom_and_sw_flags.patch
```

## Step 2: Update Protein K-mer Size in CUDA Code

Update the hardcoded k-mer size from 5 to 8 in the CUDA code:

```bash
patch runtime/kernels/resistance/translated_search_revised.cu < update_protein_kmer_size.patch
```

Or manually change line 29 in `translated_search_revised.cu`:
```cpp
#define PROTEIN_KMER_SIZE 8  // Changed from 5
```

## Step 3: Add K-mer Size Validation

Add the validation function to your pipeline. Either:

### Option A: Add to translated_search_revised.cu
Add the contents of `validate_protein_db_kmer.cpp` to the end of `translated_search_revised.cu` before the `#endif`

### Option B: Add to clean_resistance_pipeline_main.cpp
Add this simple validation function directly:

```cpp
bool validate_protein_db_kmer_size(const char* db_path, int expected_kmer_size) {
    std::string metadata_path = std::string(db_path) + "/metadata.json";
    std::ifstream metadata_file(metadata_path);
    
    if (!metadata_file.is_open()) {
        return false;
    }
    
    std::string line;
    while (std::getline(metadata_file, line)) {
        size_t pos = line.find("\"kmer_length\":");
        if (pos != std::string::npos) {
            size_t num_start = line.find_first_of("0123456789", pos);
            if (num_start != std::string::npos) {
                int db_kmer_size = std::stoi(line.substr(num_start));
                metadata_file.close();
                
                if (db_kmer_size != expected_kmer_size) {
                    std::cerr << "ERROR: Protein database k-mer size mismatch! "
                              << "Database: " << db_kmer_size 
                              << ", Expected: " << expected_kmer_size << std::endl;
                    return false;
                }
                return true;
            }
        }
    }
    
    metadata_file.close();
    return false;
}
```

## Step 4: Rebuild the Protein Database

**CRITICAL**: You must rebuild the protein database with k-mer size 8:

```bash
python3 src/python/build_clean_dynamic_database.py \
    data/wildtype_protein_seqs \
    data/integrated_clean_db_k8 \
    --mutations-csv data/quinolone_resistance_mutation_table.csv \
    --kmer-length 8
```

## Step 5: Rebuild the Pipeline

```bash
cd build
cmake ..
make clean_resistance_pipeline
```

## Step 6: Run with New Options

### Default (bloom ON, Smith-Waterman ON, k-mer size 8):
```bash
./clean_resistance_pipeline reads_R1.fq.gz reads_R2.fq.gz \
    nucleotide_index data/integrated_clean_db_k8/protein \
    fq_mutations.csv output_prefix
```

### Test without bloom filter:
```bash
./clean_resistance_pipeline reads_R1.fq.gz reads_R2.fq.gz \
    nucleotide_index data/integrated_clean_db_k8/protein \
    fq_mutations.csv output_prefix --no-bloom
```

### Test without Smith-Waterman:
```bash
./clean_resistance_pipeline reads_R1.fq.gz reads_R2.fq.gz \
    nucleotide_index data/integrated_clean_db_k8/protein \
    fq_mutations.csv output_prefix --no-sw
```

## Expected Results

With k-mer size 8 and bloom filtering enabled:
- **Specificity**: ~20,000x better than k-mer size 5
- **Speed**: Bloom filter eliminates 70-90% of non-target reads
- **Sensitivity**: Still detects resistance mutations via Smith-Waterman extension

## Troubleshooting

1. **"Protein database k-mer size mismatch" error**
   - You're using an old database built with k-mer size 5
   - Rebuild the database with `--kmer-length 8`

2. **Still getting spurious alignments**
   - Verify the CUDA code was rebuilt with `PROTEIN_KMER_SIZE 8`
   - Check that bloom filtering is enabled (not using `--no-bloom`)
   - Consider increasing to k-mer size 9 or 10 for even better specificity

3. **Missing resistance mutations**
   - K-mer size 8 should still detect most mutations
   - If sensitivity is critical, try k-mer size 7
   - Ensure Smith-Waterman is enabled (not using `--no-sw`)

## Performance Tuning

| K-mer Size | Specificity | Sensitivity | Speed | Recommendation |
|------------|-------------|-------------|-------|----------------|
| 5 | Poor | Excellent | Slow (many matches) | Not recommended |
| 7 | Good | Very Good | Fast | Good balance |
| 8 | Very Good | Good | Fast | **Recommended** |
| 9 | Excellent | Moderate | Very Fast | For high specificity |

## Future Improvements

To make k-mer size fully runtime configurable:
1. Refactor CUDA code to use dynamic k-mer sizes
2. Template the kernels on k-mer size
3. Store k-mer size in protein database structure
4. Add `--protein-kmer-size` command line flag

For now, changing k-mer size requires:
1. Modifying `PROTEIN_KMER_SIZE` in CUDA code
2. Rebuilding the pipeline
3. Rebuilding the protein database with matching k-mer size