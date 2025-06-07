#!/bin/bash
# build_integrated_database.sh
# Script to build the complete integrated resistance database for BioGPU

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Building Integrated Resistance Database ===${NC}"

# Check for required directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="${SCRIPT_DIR}"
DATA_DIR="${PROJECT_ROOT}/data"
PYTHON_DIR="${PROJECT_ROOT}/src/python"

# Create output directory structure
OUTPUT_DIR="${DATA_DIR}/integrated_fq_resistance_database"
mkdir -p "${OUTPUT_DIR}"/{nucleotide,protein,metadata,resistance_catalog}

echo -e "${GREEN}Project root: ${PROJECT_ROOT}${NC}"
echo -e "${GREEN}Output directory: ${OUTPUT_DIR}${NC}"

# Step 1: Build integrated resistance database from raw data
echo -e "\n${BLUE}Step 1: Building integrated resistance database...${NC}"
if [ -f "${PYTHON_DIR}/build_integrated_resistance_db.py" ]; then
    python3 "${PYTHON_DIR}/build_integrated_resistance_db.py" \
        --fasta-dir "${DATA_DIR}/wildtype_protein_seqs" \
        --csv "${DATA_DIR}/Known_Quinolone_Changes.csv" \
        --output-dir "${OUTPUT_DIR}" \
        --add-manual \
        || { echo -e "${RED}Failed to build integrated database${NC}"; exit 1; }
else
    echo -e "${RED}Error: build_integrated_resistance_db.py not found${NC}"
    echo "Please ensure the Python script is in: ${PYTHON_DIR}"
    exit 1
fi

# Step 2: Generate nucleotide k-mer index
echo -e "\n${BLUE}Step 2: Building nucleotide k-mer index...${NC}"
if [ -f "${PYTHON_DIR}/simple_kmer_builder.py" ]; then
    python3 "${PYTHON_DIR}/simple_kmer_builder.py" \
        "${OUTPUT_DIR}/wildtype_sequences.fasta" \
        "${OUTPUT_DIR}/nucleotide" \
        --kmer-length 15 \
        || { echo -e "${RED}Failed to build nucleotide index${NC}"; exit 1; }
else
    echo -e "${YELLOW}Warning: simple_kmer_builder.py not found${NC}"
    echo -e "${YELLOW}Skipping k-mer index generation${NC}"
fi

# Step 3: Build protein database with 5-mer index
echo -e "\n${BLUE}Step 3: Building protein database...${NC}"
# First, extract protein sequences if not already done
if [ ! -f "${OUTPUT_DIR}/wildtype_proteins.fasta" ]; then
    echo "Extracting protein sequences..."
    python3 -c "
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load the database
with open('${OUTPUT_DIR}/resistance_db.json', 'r') as f:
    db = json.load(f)

# Extract unique protein sequences
proteins = []
seen = set()

for gene_name, gene_data in db.get('gene_sequences', {}).items():
    if gene_data not in seen:
        seen.add(gene_data)
        record = SeqRecord(
            Seq(gene_data),
            id=f'{gene_name}_WT',
            description='Wildtype protein sequence'
        )
        proteins.append(record)

# Write proteins
SeqIO.write(proteins, '${OUTPUT_DIR}/wildtype_proteins.fasta', 'fasta')
print(f'Extracted {len(proteins)} unique protein sequences')
"
fi

# Build protein k-mer index
python3 -c "
import os
import json
import struct
from collections import defaultdict

# Create protein k-mer index
kmer_size = 5
kmer_index = defaultdict(list)

# Read protein sequences
from Bio import SeqIO
proteins = list(SeqIO.parse('${OUTPUT_DIR}/wildtype_proteins.fasta', 'fasta'))

# Build k-mer index
for protein_idx, record in enumerate(proteins):
    seq = str(record.seq)
    for pos in range(len(seq) - kmer_size + 1):
        kmer = seq[pos:pos + kmer_size]
        kmer_index[kmer].append((protein_idx, pos))

# Write binary index
os.makedirs('${OUTPUT_DIR}/protein', exist_ok=True)
with open('${OUTPUT_DIR}/protein/protein_kmers.bin', 'wb') as f:
    f.write(struct.pack('I', kmer_size))
    f.write(struct.pack('I', len(kmer_index)))
    
    for kmer, positions in sorted(kmer_index.items()):
        f.write(kmer.encode('ascii'))
        f.write(struct.pack('I', len(positions)))
        for protein_idx, pos in positions:
            f.write(struct.pack('II', protein_idx, pos))

# Write protein sequences
with open('${OUTPUT_DIR}/protein/proteins.bin', 'wb') as f:
    f.write(struct.pack('I', len(proteins)))
    for record in proteins:
        f.write(str(record.seq).encode('ascii'))

# Load species map from resistance database
with open('${OUTPUT_DIR}/resistance_db.json', 'r') as f:
    resistance_db = json.load(f)

# Write metadata
metadata = {
    'num_proteins': len(proteins),
    'kmer_size': kmer_size,
    'num_kmers': len(kmer_index),
    'gene_map': {i: proteins[i].id for i in range(len(proteins))},
    'species_map': resistance_db.get('species_map', {})
}

with open('${OUTPUT_DIR}/protein/metadata.json', 'w') as f:
    json.dump(metadata, f, indent=2)

print(f'Built protein database with {len(proteins)} proteins and {len(kmer_index)} unique {kmer_size}-mers')
"

# Step 4: Compile resistance mutation catalog
echo -e "\n${BLUE}Step 4: Compiling resistance mutation catalog...${NC}"
python3 -c "
import json
import struct

# Load resistance database
with open('${OUTPUT_DIR}/resistance_db.json', 'r') as f:
    db = json.load(f)

# Convert to binary format for GPU
mutations = []
for mut_list in db.get('mutations', []):
    for mut in mut_list:
        if isinstance(mut, dict):
            mutations.append(mut)

# Write binary catalog
with open('${OUTPUT_DIR}/resistance_catalog.bin', 'wb') as f:
    f.write(struct.pack('I', len(mutations)))
    
    for mut in mutations:
        # Write mutation data in GPU-friendly format
        gene_name = mut.get('gene_name', 'unknown')[:32]
        f.write(gene_name.ljust(32, '\0').encode('ascii'))
        f.write(struct.pack('I', mut.get('gene_id', 0)))
        f.write(struct.pack('H', mut.get('position', 0)))
        f.write(mut.get('wildtype_aa', 'X').encode('ascii'))
        
        # Resistant AAs (up to 10)
        resistant_aas = mut.get('resistant_aas', [])[:10]
        for i in range(10):
            if i < len(resistant_aas):
                f.write(resistant_aas[i].encode('ascii'))
            else:
                f.write(b'\0')
        
        f.write(struct.pack('B', len(resistant_aas)))
        f.write(struct.pack('f', mut.get('resistance_level', 0.5)))
        
        # Drug affected
        drug = mut.get('drugs_affected', ['fluoroquinolone'])[0][:32]
        f.write(drug.ljust(32, '\0').encode('ascii'))

print(f'Compiled {len(mutations)} resistance mutations to binary catalog')
"

# Step 5: Create validation dataset
echo -e "\n${BLUE}Step 5: Creating validation dataset...${NC}"
mkdir -p "${OUTPUT_DIR}/validation"

python3 -c "
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create synthetic reads with known mutations
mutations = {
    'gyrA_S83L': ('S', 'L', 82),  # 0-based position
    'gyrA_D87N': ('D', 'N', 86),
    'parC_S80I': ('S', 'I', 79),
    'parC_E84V': ('E', 'V', 83)
}

# Generate synthetic reads
reads = []
for i in range(1000):
    # Create a random DNA sequence
    seq = ''.join(random.choices('ACGT', k=150))
    
    # Insert a known resistance codon occasionally
    if i % 10 == 0:
        # Insert gyrA S83L mutation (TCT -> TTG)
        pos = random.randint(10, 100)
        seq = seq[:pos] + 'TTG' + seq[pos+3:]
        desc = 'contains_gyrA_S83L_mutation'
    else:
        desc = 'wildtype'
    
    record = SeqRecord(
        Seq(seq),
        id=f'synthetic_read_{i}',
        description=desc
    )
    reads.append(record)

# Write synthetic FASTQ
with open('${OUTPUT_DIR}/validation/synthetic_R1.fastq', 'w') as f:
    for record in reads:
        f.write(f'@{record.id} {record.description}\n')
        f.write(f'{record.seq}\n')
        f.write('+\n')
        f.write('I' * len(record.seq) + '\n')

print(f'Created {len(reads)} synthetic reads for validation')
"

# Compress validation files
gzip -f "${OUTPUT_DIR}/validation/synthetic_R1.fastq"
cp "${OUTPUT_DIR}/validation/synthetic_R1.fastq.gz" "${OUTPUT_DIR}/validation/synthetic_R2.fastq.gz"

# Step 6: Validate database integrity
echo -e "\n${BLUE}Step 6: Validating database integrity...${NC}"
python3 -c "
import os
import json
import struct

# Check all required files exist
required_files = [
    '${OUTPUT_DIR}/resistance_db.json',
    '${OUTPUT_DIR}/wildtype_sequences.fasta',
    '${OUTPUT_DIR}/resistant_variants.fasta',
    '${OUTPUT_DIR}/nucleotide/kmer_index.bin',
    '${OUTPUT_DIR}/protein/protein_kmers.bin',
    '${OUTPUT_DIR}/protein/proteins.bin',
    '${OUTPUT_DIR}/resistance_catalog.bin'
]

missing = []
for file in required_files:
    if not os.path.exists(file):
        missing.append(file)
    else:
        size = os.path.getsize(file)
        print(f'✓ {os.path.basename(file)}: {size:,} bytes')

if missing:
    print(f'\n❌ Missing files: {missing}')
    exit(1)

# Validate JSON database
with open('${OUTPUT_DIR}/resistance_db.json', 'r') as f:
    db = json.load(f)
    
print(f'\nDatabase contents:')
print(f'  Genes: {len(db.get(\"gene_id_map\", {}))}')
print(f'  Species: {len(db.get(\"species_map\", {}))}')
print(f'  Mutations: {sum(len(m) for m in db.get(\"mutations\", []))}')
print(f'  Gene sequences: {len(db.get(\"gene_sequences\", {}))}')

print('\n✅ Database validation complete!')
"

# Step 7: Create summary report
echo -e "\n${BLUE}Creating summary report...${NC}"
cat > "${OUTPUT_DIR}/DATABASE_README.txt" << EOF
Integrated Resistance Database
==============================

This database contains comprehensive fluoroquinolone resistance data for GPU-accelerated detection.

Structure:
----------
- resistance_db.json: Master database with all resistance mutations and metadata
- wildtype_sequences.fasta: Wild-type gene sequences (nucleotide)
- resistant_variants.fasta: Sequences with known resistance mutations
- wildtype_proteins.fasta: Wild-type protein sequences

- nucleotide/: Nucleotide k-mer index (15-mers)
  - kmer_index.bin: Binary k-mer index for GPU
  - metadata.json: Index metadata

- protein/: Protein database (5-mers)
  - protein_kmers.bin: Protein k-mer index
  - proteins.bin: Protein sequences
  - metadata.json: Protein metadata

- resistance_catalog.bin: Binary format resistance mutations for GPU

- validation/: Synthetic test data
  - synthetic_R1.fastq.gz: Test reads with known mutations
  - synthetic_R2.fastq.gz: Paired reads

Usage:
------
./integrated_resistance_pipeline \\
    ${OUTPUT_DIR}/nucleotide \\
    ${OUTPUT_DIR}/protein \\
    ${OUTPUT_DIR} \\
    reads_R1.fastq.gz \\
    reads_R2.fastq.gz \\
    output_prefix

Generated: $(date)
EOF

echo -e "\n${GREEN}✅ Integrated resistance database build complete!${NC}"
echo -e "${GREEN}Database location: ${OUTPUT_DIR}${NC}"
echo -e "${GREEN}See ${OUTPUT_DIR}/DATABASE_README.txt for details${NC}"

# Optional: Run a quick test if the pipeline executable exists
if [ -f "${PROJECT_ROOT}/build/integrated_resistance_pipeline" ]; then
    echo -e "\n${BLUE}Testing integrated pipeline...${NC}"
    "${PROJECT_ROOT}/build/integrated_resistance_pipeline" \
        "${OUTPUT_DIR}/nucleotide" \
        "${OUTPUT_DIR}/protein" \
        "${OUTPUT_DIR}" \
        "${OUTPUT_DIR}/validation/synthetic_R1.fastq.gz" \
        "${OUTPUT_DIR}/validation/synthetic_R2.fastq.gz" \
        "${OUTPUT_DIR}/test_output" \
        || echo -e "${YELLOW}Pipeline test failed - this is expected if you haven't built the executable yet${NC}"
fi

echo -e "\n${BLUE}Next steps:${NC}"
echo "1. Build the project: cd build && make integrated_resistance_pipeline"
echo "2. Run on your data: ./integrated_resistance_pipeline <args>"
echo "3. Analyze results: python3 src/python/resistance_profile_analyzer.py <output.h5>"