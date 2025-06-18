# Bloom Filter and Smith-Waterman Test Results

## Test Configuration
- **Input files**: data/569_A_038_R1.fastq.gz, data/569_A_038_R2.fastq.gz
- **Nucleotide index**: data/integrated_clean_db/nucleotide
- **Protein database**: data/integrated_clean_db/protein
- **Mutations CSV**: data/quinolone_resistance_mutation_table.csv

## Test Commands and Results

### 1. Default (Bloom: ON, SW: ON)
```bash
./build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    results/569_A_default_bloom_sw \
    data/quinolone_resistance_mutation_table.csv
```
- **Start time**: Sun Jun  8 09:29:25 PM EDT 2025
- **End time**: Sun Jun  8 09:29:43 PM EDT 2025
- **Runtime**: 18.807s (17s processing time)
- **Reads processed**: 3,405,007
- **Reads/second**: 200,295
- **Protein matches found**: 1,326

### 2. No Bloom (Bloom: OFF, SW: ON)
```bash
./build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    results/569_A_no_bloom \
    data/quinolone_resistance_mutation_table.csv \
    --no-bloom
```
- **Start time**: Sun Jun  8 09:29:59 PM EDT 2025
- **End time**: Sun Jun  8 09:30:17 PM EDT 2025
- **Runtime**: 17.993s (16s processing time)
- **Reads processed**: 3,405,007
- **Reads/second**: 212,813
- **Protein matches found**: 1,326

### 3. No SW (Bloom: ON, SW: OFF)
```bash
./build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    results/569_A_no_sw \
    data/quinolone_resistance_mutation_table.csv \
    --no-sw
```
- **Start time**: Sun Jun  8 09:30:32 PM EDT 2025
- **End time**: Sun Jun  8 09:30:50 PM EDT 2025
- **Runtime**: 18.432s (16s processing time)
- **Reads processed**: 3,405,007
- **Reads/second**: 212,813
- **Protein matches found**: 922

### 4. No Bloom, No SW (Bloom: OFF, SW: OFF)
```bash
./build/clean_resistance_pipeline \
    data/integrated_clean_db/nucleotide \
    data/integrated_clean_db/protein \
    data/569_A_038_R1.fastq.gz \
    data/569_A_038_R2.fastq.gz \
    results/569_A_no_bloom_no_sw \
    data/quinolone_resistance_mutation_table.csv \
    --no-bloom --no-sw
```
- **Start time**: Sun Jun  8 09:31:05 PM EDT 2025
- **End time**: Sun Jun  8 09:31:22 PM EDT 2025
- **Runtime**: 17.608s (15s processing time)
- **Reads processed**: 3,405,007
- **Reads/second**: 227,000
- **Protein matches found**: 922

## Summary

### Performance Comparison
| Configuration | Runtime | Reads/sec | Protein Matches |
|--------------|---------|-----------|----------------|
| Bloom ON, SW ON (default) | 18.8s | 200,295 | 1,326 |
| Bloom OFF, SW ON | 18.0s | 212,813 | 1,326 |
| Bloom ON, SW OFF | 18.4s | 212,813 | 922 |
| Bloom OFF, SW OFF | 17.6s | 227,000 | 922 |

### Key Findings
1. **Bloom filter impact**: Minimal performance difference (~6% improvement when disabled)
2. **Smith-Waterman impact**: Reduces protein matches by ~30% when disabled (1,326 → 922)
3. **Fastest configuration**: Both disabled (227,000 reads/sec)
4. **Most sensitive configuration**: Default with both enabled (1,326 matches)

### Working Commands
All four command variations worked successfully:
- Default: `./build/clean_resistance_pipeline <args>`
- No bloom: `./build/clean_resistance_pipeline <args> --no-bloom`
- No SW: `./build/clean_resistance_pipeline <args> --no-sw`
- Both disabled: `./build/clean_resistance_pipeline <args> --no-bloom --no-sw`