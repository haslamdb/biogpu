# Batch Processing with CSV Input

The Clean Resistance Pipeline now supports batch processing of multiple samples using a CSV input file. This allows you to process many samples without manual command line entry for each one.

## CSV File Format

The CSV file should contain the following columns (in this order):

1. **Sample Name Column** (required): Can be named any of:
   - `Sample Name`, `SampleName`, `Sample ID`, `SampleID`
   - `sample_name`, `samplename`, `sample_id`, `sampleid`
   - `Sample`, `sample`, `ID`, `id`, `Name`, `name`
  
2. **File Path Column** (required): The directory containing the FASTQ files. Can be named any of:
   - `FilePath`, `File Path`, `filepath`, `file_path`
   - `Path`, `path`, `Directory`, `directory`, `Dir`, `dir`
   - `BaseDir`, `basedir`, `Base Dir`, `base_dir`
  
3. **R1 File Column** (required): The filename for Read 1 (without path). Can be named any of:
   - `R1 file`, `R1file`, `R1_file`, `R1`, `r1`
   - `Read1`, `read1`, `Forward`, `forward`
   - `File1`, `file1`, `FASTQ1`, `fastq1`
  
4. **R2 File Column** (optional): The filename for Read 2 (without path). Can be named any of:
   - `R2 file`, `R2file`, `R2_file`, `R2`, `r2`
   - `Read2`, `read2`, `Reverse`, `reverse`
   - `File2`, `file2`, `FASTQ2`, `fastq2`

### Example CSV Files

**Standard format (comma-separated):**
```csv
SampleName,FilePath,R1 file,R2 file
569_06_B,~/MSSData/,569_06_B_R1.fastq.gz,569_06_B_R2.fastq.gz
570_07_A,~/MSSData/,570_07_A_R1.fastq.gz,570_07_A_R2.fastq.gz
Patient_001,/data/sequencing/,Patient_001_S1_L001_R1_001.fastq.gz,Patient_001_S1_L001_R2_001.fastq.gz
```

**Tab-separated format:**
```tsv
Sample ID	FilePath	R1 file	R2 file
Patient_001	/data/sequencing/	Patient_001_R1.fastq.gz	Patient_001_R2.fastq.gz
Patient_002	/data/sequencing/	Patient_002_R1.fastq.gz	Patient_002_R2.fastq.gz
```

**With additional metadata columns:**
```csv
SampleName,FilePath,R1 file,R2 file,PatientID,CollectionDate,SampleType
569_06_B,~/MSSData/,569_06_B_R1.fastq.gz,569_06_B_R2.fastq.gz,P123,2024-01-15,Urine
570_07_A,~/MSSData/,570_07_A_R1.fastq.gz,570_07_A_R2.fastq.gz,P124,2024-01-16,Blood
```

### Important Notes:
- The `FilePath` should end with a directory separator (`/`), but the parser will add it if missing
- If R1 and R2 filenames contain path information, it will be stripped and only the filename will be used
- The full path to each file is constructed as: `FilePath + R1/R2 filename`
- Tilde (`~`) in paths is automatically expanded to the home directory

## Command Line Usage

### Batch Mode
```bash
./clean_resistance_pipeline <nucleotide_index> <protein_db> --csv <samples.csv> [options]
```

### Example Commands

**Basic batch processing:**
```bash
./clean_resistance_pipeline \
    /data/indices/biogpu_nucleotide \
    /data/databases/biogpu_protein \
    --csv samples.csv \
    --output-dir results
```

**With all options:**
```bash
./clean_resistance_pipeline \
    /data/indices/biogpu_nucleotide \
    /data/databases/biogpu_protein \
    --csv samples.csv \
    --output-dir /data/results \
    --fq-csv /data/quinolone_resistance_mutations.csv \
    --min-allele-depth 10 \
    --min-report-depth 5 \
    --stop-on-error
```

**Dry run to validate CSV:**
```bash
./clean_resistance_pipeline \
    /data/indices/biogpu_nucleotide \
    /data/databases/biogpu_protein \
    --csv samples.csv \
    --dry-run
```

## Options

- `--csv <file>`: CSV file containing sample information
- `--output-dir <dir>`: Base output directory (default: results)
- `--fq-csv <path>`: Path to FQ resistance mutations CSV
- `--no-bloom`: Disable bloom filter pre-screening
- `--no-sw`: Disable Smith-Waterman alignment
- `--min-allele-depth <N>`: Minimum depth for allele frequency analysis (default: 5)
- `--min-report-depth <N>`: Minimum depth for reporting polymorphisms (default: 0)
- `--stop-on-error`: Stop batch processing on first error
- `--dry-run`: Validate CSV and show what would be processed without running

## Output Structure

When using batch mode, the pipeline creates a structured output directory:

```
results/
├── Sample1/
│   ├── Sample1.h5
│   ├── Sample1.json
│   ├── Sample1_protein_matches.csv
│   ├── Sample1_allele_frequencies.csv
│   └── Sample1_clinical_report.html
├── Sample2/
│   ├── Sample2.h5
│   ├── Sample2.json
│   ├── Sample2_protein_matches.csv
│   ├── Sample2_allele_frequencies.csv
│   └── Sample2_clinical_report.html
└── batch_summary.log
```

## Features

1. **Automatic delimiter detection**: The parser automatically detects whether your CSV uses commas, tabs, or semicolons
2. **Path expansion**: Tilde (~) in paths is automatically expanded to home directory
3. **File validation**: The parser checks if input files exist before processing
4. **Flexible column naming**: Recognizes many variations of column names
5. **Progress tracking**: Shows progress for each sample in the batch
6. **Error handling**: Can continue processing even if some samples fail

## Tips

1. **Use absolute paths** in your CSV when possible to avoid path resolution issues
2. **Test with dry-run** first to validate your CSV format
3. **Use meaningful sample names** as they will be used for output directories
4. **Check file permissions** - ensure read access to all FASTQ files
5. **Monitor GPU memory** when processing many samples consecutively

## Integration with Other Tools

The CSV parser is implemented as a reusable module that can be used in other BioGPU pipelines:

```cpp
#include "sample_csv_parser.h"

BioGPU::SampleCSVParser parser;
if (parser.parseFile("samples.csv")) {
    for (size_t i = 0; i < parser.getSampleCount(); i++) {
        const BioGPU::SampleInfo* sample = parser.getSample(i);
        // Process sample...
    }
}
```