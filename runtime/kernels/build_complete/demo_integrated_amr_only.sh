#!/bin/bash
# Demo script showing integrated pipeline with AMR detection only

echo "BioGPU Integrated Pipeline Demo"
echo "==============================="
echo ""
echo "This demo shows the unified pipeline running AMR detection."
echo "The resistance pipeline will be skipped since it's not built."
echo ""

# Create demo directories
mkdir -p /tmp/biogpu_demo/{input,output,databases}

echo "Creating minimal test data..."
# Create minimal test FASTQ files
cat > /tmp/biogpu_demo/input/demo_R1.fastq << 'EOF'
@SEQ_ID_001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID_002
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > /tmp/biogpu_demo/input/demo_R2.fastq << 'EOF'
@SEQ_ID_001
CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID_002
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Create placeholder database files
echo ">demo_ref" > /tmp/biogpu_demo/databases/ref.fasta
echo "ACGTACGTACGT" >> /tmp/biogpu_demo/databases/ref.fasta

echo "species,gene,position,mutation" > /tmp/biogpu_demo/databases/resistance.csv
echo "E.coli,gyrA,83,S83L" >> /tmp/biogpu_demo/databases/resistance.csv

echo ""
echo "Running integrated pipeline (AMR only)..."
echo ""

# Run the pipeline with AMR detection only
./bio_gpu_pipeline_integrated \
    --r1 /tmp/biogpu_demo/input/demo_R1.fastq \
    --r2 /tmp/biogpu_demo/input/demo_R2.fastq \
    --output-dir /tmp/biogpu_demo/output \
    --reference-db /tmp/biogpu_demo/databases/ref.fasta \
    --resistance-db /tmp/biogpu_demo/databases/resistance.csv \
    --sample-id DEMO001 \
    --disable-resistance \
    --progress-json

echo ""
echo "Demo complete!"
echo ""
echo "Output files:"
ls -la /tmp/biogpu_demo/output/
echo ""
echo "This demonstrates that the unified pipeline framework is working"
echo "and can integrate with the existing AMR detection pipeline."