# Fluoroquinolone Resistance Mutations

## Overview

The resistance detection pipeline includes a database of known fluoroquinolone resistance mutations. By default, the pipeline uses hardcoded mutations (255 mutations across 33 species) so no external CSV file is required.

## Default Behavior

When you run the pipeline without specifying `--fq-csv`, it uses the built-in hardcoded mutations:

```bash
./build/clean_resistance_pipeline \
    /path/to/nucleotide_index \
    /path/to/protein_db \
    --csv samples.csv
```

## Using a Custom Mutations CSV

You can override the hardcoded mutations by providing your own CSV file:

```bash
./build/clean_resistance_pipeline \
    /path/to/nucleotide_index \
    /path/to/protein_db \
    --csv samples.csv \
    --fq-csv /path/to/custom_mutations.csv
```

## CSV Format

The FQ resistance mutations CSV must have the following format:

```csv
species,gene,location,wt,mut
Escherichia_coli,gyrA,83,S,L
Escherichia_coli,gyrA,87,D,N
Staphylococcus_aureus,grlA,80,S,F
```

### Column Descriptions:
- **species**: Species name with underscores instead of spaces (e.g., `Escherichia_coli`)
- **gene**: Gene name (e.g., `gyrA`, `gyrB`, `parC`, `parE`, `grlA`)
- **location**: Amino acid position (1-based)
- **wt**: Wild-type (expected) amino acid at this position
- **mut**: Mutant amino acid that confers resistance

## Adding Custom Mutations

To add your own mutations to the hardcoded database:

1. **Start with the existing CSV** (recommended):
   ```bash
   cp /home/david/Documents/Code/biogpu/data/quinolone_resistance_mutation_table.csv my_mutations.csv
   ```

2. **Add your mutations** to the CSV file using any text editor or spreadsheet program

3. **Convert to C++ header**:
   ```bash
   python3 convert_fq_csv_to_cpp.py my_mutations.csv > fq_mutations_hardcoded.h
   ```

4. **Rebuild the pipeline**:
   ```bash
   cd build
   make clean_resistance_pipeline
   ```

## Helper Script Usage

The `convert_fq_csv_to_cpp.py` script converts a mutations CSV file to C++ code:

```bash
# Convert specific CSV file
python3 convert_fq_csv_to_cpp.py /path/to/mutations.csv > fq_mutations_hardcoded.h

# Use default path (if no argument provided)
python3 convert_fq_csv_to_cpp.py > fq_mutations_hardcoded.h
```

The script will:
- Parse the CSV file
- Generate a C++ static array with all mutations
- Include summary statistics (mutations per species and gene)
- Skip any entries with "NA" values

## Example: Adding New Mutations

Let's say you want to add a new resistance mutation for *Klebsiella pneumoniae*:

1. Edit your CSV file and add:
   ```csv
   Klebsiella_pneumoniae,gyrA,91,T,I
   ```

2. Regenerate the C++ header:
   ```bash
   python3 convert_fq_csv_to_cpp.py my_mutations.csv > fq_mutations_hardcoded.h
   ```

3. Rebuild:
   ```bash
   cd build && make clean_resistance_pipeline
   ```

## Supported Species and Genes

The default database includes mutations for:

### Species (33 total):
- Acinetobacter baumannii
- Burkholderia species
- Campylobacter species
- Clostridioides difficile
- Enterobacteriaceae (E. coli, Klebsiella, Salmonella, etc.)
- Enterococcus species
- Haemophilus influenzae
- Neisseria species
- Pseudomonas aeruginosa
- Staphylococcus species
- Streptococcus pneumoniae
- Vibrio species

### Genes (7 total):
- **gyrA** - DNA gyrase subunit A (121 mutations)
- **gyrB** - DNA gyrase subunit B (19 mutations)
- **parC** - Topoisomerase IV subunit C (70 mutations)
- **parE** - Topoisomerase IV subunit E (36 mutations)
- **grlA** - Topoisomerase IV subunit A (S. aureus) (6 mutations)
- **emrR** - Multidrug resistance regulator (2 mutations)
- **nfxB** - Transcriptional regulator (1 mutation)

## Notes

- Species names must use underscores, not spaces
- Gene names are case-sensitive (use lowercase)
- Position numbering is 1-based (amino acid positions)
- The pipeline will warn but continue if the FQ CSV file cannot be loaded
- Invalid entries (with "NA" values) are automatically skipped