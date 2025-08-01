#!/bin/bash
# Final integration build script using process spawning approach

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}BioGPU Final Integration Build${NC}"
echo -e "${GREEN}========================================${NC}"

# Step 1: Build the process-spawning unified pipeline
echo -e "${YELLOW}Building process-spawning unified pipeline...${NC}"

nvcc -o bio_gpu_pipeline_integrated unified_pipeline_process_spawning.cpp \
     -std=c++17 \
     -lcudart \
     -lpthread \
     -O3 \
     -Wno-deprecated-gpu-targets

if [ ! -f "bio_gpu_pipeline_integrated" ]; then
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi

echo -e "${GREEN}Build successful!${NC}"

# Step 2: Check for existing pipeline executables
echo -e "${YELLOW}Checking for pipeline executables...${NC}"

AMR_FOUND=false
RES_FOUND=false

if [ -f "../genes/amr_detection" ]; then
    echo -e "${GREEN}✓ Found AMR detection executable${NC}"
    AMR_FOUND=true
else
    echo -e "${YELLOW}✗ AMR detection executable not found${NC}"
fi

# Look for resistance pipeline in multiple locations
for exe in "clean_resistance_pipeline" "resistance_pipeline"; do
    for dir in "../resistance" "../resistance/build" "./resistance"; do
        if [ -f "$dir/$exe" ]; then
            echo -e "${GREEN}✓ Found resistance pipeline: $dir/$exe${NC}"
            RES_FOUND=true
            break 2
        fi
    done
done

if [ "$RES_FOUND" = false ]; then
    echo -e "${YELLOW}✗ Resistance pipeline executable not found${NC}"
fi

# Step 3: Create test script
echo -e "${YELLOW}Creating test script...${NC}"

cat > test_integrated_pipeline.sh << 'EOF'
#!/bin/bash
# Test script for integrated pipeline

echo "BioGPU Integrated Pipeline Test"
echo "==============================="

# Create test directories
mkdir -p /tmp/biogpu_test/{input,output}

# Create dummy test files (empty for validation test)
touch /tmp/biogpu_test/input/test_R1.fastq.gz
touch /tmp/biogpu_test/input/test_R2.fastq.gz

# Test help
echo -e "\n1. Testing help..."
./bio_gpu_pipeline_integrated --help

# Test validation (will fail on empty files, but tests argument parsing)
echo -e "\n2. Testing argument validation..."
./bio_gpu_pipeline_integrated \
    --r1 /tmp/biogpu_test/input/test_R1.fastq.gz \
    --r2 /tmp/biogpu_test/input/test_R2.fastq.gz \
    --output-dir /tmp/biogpu_test/output \
    --reference-db /tmp/ref_db \
    --resistance-db /tmp/res_db \
    --sample-id TEST001 \
    --progress-json 2>&1 | head -5

echo -e "\nTest complete!"
EOF

chmod +x test_integrated_pipeline.sh

# Step 4: Create installation package
echo -e "${YELLOW}Creating installation package...${NC}"

mkdir -p install_final
cp bio_gpu_pipeline_integrated install_final/bio_gpu_pipeline

cat > install_final/install.sh << 'INSTALL_EOF'
#!/bin/bash
echo "Installing BioGPU Unified Pipeline (Process Integration)..."

# Check for required executables
echo "Checking dependencies..."
AMR_OK=false
RES_OK=false

if [ -f "/usr/local/bin/amr_detection" ] || find .. -name "amr_detection" -type f -executable 2>/dev/null | head -1; then
    AMR_OK=true
    echo "✓ AMR detection available"
else
    echo "⚠ AMR detection not found - install separately"
fi

if [ -f "/usr/local/bin/clean_resistance_pipeline" ] || find .. -name "*resistance*pipeline*" -type f -executable 2>/dev/null | head -1; then
    RES_OK=true
    echo "✓ Resistance pipeline available"
else
    echo "⚠ Resistance pipeline not found - install separately"
fi

# Install unified pipeline
sudo cp bio_gpu_pipeline /usr/local/bin/
sudo chmod +x /usr/local/bin/bio_gpu_pipeline

echo ""
echo "Installation complete!"
echo "Run: bio_gpu_pipeline --help"
echo ""
if [ "$AMR_OK" = false ] || [ "$RES_OK" = false ]; then
    echo "Note: Some pipeline components are missing. Install them separately."
fi
INSTALL_EOF

chmod +x install_final/install.sh

# Step 5: Summary
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Build Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Executable: $(pwd)/bio_gpu_pipeline_integrated"
echo "Size: $(du -h bio_gpu_pipeline_integrated | cut -f1)"
echo ""
echo "This version uses process spawning to run existing pipelines:"
echo "- Searches for pipeline executables in multiple locations"
echo "- Runs pipelines as separate processes with proper GPU assignment"
echo "- Captures and reports progress in real-time"
echo "- Supports all original command-line arguments"
echo ""
echo "Test with: ./test_integrated_pipeline.sh"
echo "Install with: cd install_final && ./install.sh"
echo ""
if [ "$AMR_FOUND" = true ] && [ "$RES_FOUND" = true ]; then
    echo -e "${GREEN}Both pipeline executables found - ready for full integration!${NC}"
else
    echo -e "${YELLOW}Some pipeline executables not found - build them separately${NC}"
fi