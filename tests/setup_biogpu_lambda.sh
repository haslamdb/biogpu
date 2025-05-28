#!/bin/bash
# setup_biogpu_env.sh
# Setup script for BioGPU development on Titan Xp server

echo "Setting up BioGPU development environment..."
echo "Optimized for dual NVIDIA Titan Xp GPUs"

# # Create project directory
# mkdir -p ~/biogpu-dev
# cd ~/biogpu-dev

# # Check GPU status
# echo "Checking GPU status..."
# nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv

# Create conda environment
echo "Creating conda environment..."
conda create -n biogpu python=3.12 -y
source ~/miniconda3/etc/profile.d/conda.sh  # or ~/anaconda3/etc/profile.d/conda.sh
conda activate biogpu

# Install Python packages
echo "Installing Python packages..."
# pip install cupy-cuda11x  # Use cuda12x if you have CUDA 12
pip install biopython pysam pandas numpy matplotlib seaborn
pip install jupyterlab ipython

# Install bioinformatics tools via conda
echo "Installing bioinformatics tools..."
conda install -c bioconda minimap2 bwa samtools -y

# Create project structure
echo "Creating project structure..."
# mkdir -p {src/{python,cuda,cpp},tests,data/{raw,processed},results,scripts}

# Create CMakeLists.txt for CUDA compilation
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.18)
project(BioGPU LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

find_package(CUDA REQUIRED)

# Set CUDA architecture for Titan Xp (compute capability 6.1)
set(CMAKE_CUDA_ARCHITECTURES 61)

# CUDA kernel library
add_library(fq_resistance_kernels SHARED
    src/cuda/fq_resistance_kernels.cu
)

target_link_libraries(fq_resistance_kernels
    ${CUDA_LIBRARIES}
)

# Python extension (optional)
find_package(pybind11 QUIET)
if(pybind11_FOUND)
    pybind11_add_module(biogpu_cuda src/python/cuda_wrapper.cpp)
    target_link_libraries(biogpu_cuda PRIVATE fq_resistance_kernels)
endif()
EOF

# Create test data generation script
cat > scripts/generate_test_data.py << 'EOF'
#!/usr/bin/env python3
"""Generate test FASTQ data with known fluoroquinolone resistance mutations"""

import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def generate_test_fastq():
    """Generate synthetic FASTQ with known mutations"""
    
    # E. coli gyrA sequence around position 83 (simplified)
    gyr_a_wt = "ATGGATAAAGATGTAGCTTATTTAGATATCGATTCTGATGAAGATGATCCT"
    gyr_a_mut = "ATGGATAAAGATGTAGCTTATTTAGATATCGATTTGGATGAAGATGATCCT"  # S83L
    
    # P. aeruginosa parC sequence around position 87 (simplified)  
    par_c_wt = "ATGGATCGAGATCTAGCTTATGTATATCGATTCTGATGATGATGATCCT"
    par_c_mut = "ATGGATCGAGATCTAGCTTATGTATATCGATTTGGATGATGATGATCCT"  # S87L
    
    records = []
    
    # Generate reads with wild-type sequences
    for i in range(5000):
        seq = random.choice([gyr_a_wt, par_c_wt])
        record = SeqRecord(
            Seq(seq),
            id=f"read_wt_{i}",
            description="wild_type"
        )
        record.letter_annotations["phred_quality"] = [35] * len(seq)
        records.append(record)
    
    # Generate reads with resistance mutations
    for i in range(1000):
        seq = random.choice([gyr_a_mut, par_c_mut])
        record = SeqRecord(
            Seq(seq),
            id=f"read_mut_{i}",
            description="resistant_mutation"
        )
        record.letter_annotations["phred_quality"] = [35] * len(seq)
        records.append(record)
    
    # Write FASTQ
    with open("data/raw/test_sample.fastq", "w") as f:
        SeqIO.write(records, f, "fastq")
    
    print(f"Generated {len(records)} test reads")
    print("Wild-type reads: 5000")
    print("Resistant reads: 1000")

if __name__ == "__main__":
    generate_test_fastq()
EOF

# Create quick test script
cat > scripts/test_pipeline.py << 'EOF'
#!/usr/bin/env python3
"""Quick test of the fluoroquinolone resistance pipeline"""

import sys
import os
sys.path.append('src/python')

try:
    from fluoroquinolone_resistance_pipeline import TitanXpFQPipeline, test_pipeline
    
    print("Testing BioGPU Fluoroquinolone Resistance Pipeline")
    print("=" * 50)
    
    # Run the test
    test_pipeline()
    
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure you've activated the biogpu conda environment")
except Exception as e:
    print(f"Test failed: {e}")
    import traceback
    traceback.print_exc()
EOF

# Create build script for CUDA kernels
cat > scripts/build_cuda.sh << 'EOF'
#!/bin/bash
# Build CUDA kernels for Titan Xp

echo "Building CUDA kernels for compute capability 6.1..."

# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

echo "CUDA build complete!"
echo "Library location: build/libfq_resistance_kernels.so"
EOF

chmod +x scripts/*.py scripts/*.sh

# Test CUDA compilation
echo "Testing CUDA compilation..."
nvcc --version

# Create a simple CUDA test
cat > src/cuda/test_cuda.cu << 'EOF'
#include <cuda_runtime.h>
#include <stdio.h>

__global__ void hello_gpu() {
    printf("Hello from GPU thread %d in block %d!\n", threadIdx.x, blockIdx.x);
}

int main() {
    printf("Testing CUDA on Titan Xp...\n");
    
    // Query GPU properties
    int device_count;
    cudaGetDeviceCount(&device_count);
    
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("GPU %d: %s\n", i, prop.name);
        printf("  Compute capability: %d.%d\n", prop.major, prop.minor);
        printf("  Global memory: %.2f GB\n", prop.totalGlobalMem / 1e9);
    }
    
    // Launch simple kernel
    hello_gpu<<<2, 4>>>();
    cudaDeviceSynchronize();
    
    return 0;
}
EOF

# Test CUDA compilation
echo "Compiling CUDA test..."
cd src/cuda
nvcc -o test_cuda test_cuda.cu
if [ $? -eq 0 ]; then
    echo "CUDA compilation successful!"
    echo "Running test..."
    ./test_cuda
else
    echo "CUDA compilation failed. Check your CUDA installation."
fi

cd ~/biogpu-dev

echo ""
echo "Setup complete!"
echo "===================="
echo "To get started:"
echo "1. conda activate biogpu"
echo "2. cd ~/biogpu-dev"
echo "3. python scripts/generate_test_data.py"
echo "4. python scripts/test_pipeline.py"
echo ""
echo "Project structure:"
echo "~/biogpu-dev/"
echo "├── src/"
echo "│   ├── python/          # Python pipeline code"
echo "│   ├── cuda/            # CUDA kernels"
echo "│   └── cpp/             # C++ utilities"
echo "├── tests/               # Unit tests"
echo "├── data/"
echo "│   ├── raw/             # Input FASTQ files"
echo "│   └── processed/       # Analysis results"
echo "├── scripts/             # Utility scripts"
echo "└── results/             # Final reports"
echo ""
echo "Next steps:"
echo "- Copy the Python pipeline code to src/python/"
echo "- Copy the CUDA kernels to src/cuda/"
echo "- Build and test with real data!"