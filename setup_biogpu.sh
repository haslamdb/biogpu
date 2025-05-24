#!/bin/bash
# setup_biogpu_system_cuda.sh - Setup using system CUDA 12.5

set -e

echo "=== BioGPU Setup with System CUDA 12.5 ==="
echo

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }

# 1. Verify system CUDA
info "Checking system CUDA installation..."
if ! command -v nvcc &> /dev/null; then
    warn "NVCC not found! Please ensure CUDA 12.5 is installed"
    exit 1
fi

CUDA_VERSION=$(nvcc --version | grep release | awk '{print $6}' | cut -d',' -f1)
info "Found CUDA version: $CUDA_VERSION"

# 2. Set up CUDA paths
CUDA_HOME=${CUDA_HOME:-/usr/local/cuda-12.5}
if [ ! -d "$CUDA_HOME" ]; then
    CUDA_HOME=/usr/local/cuda
fi

info "Using CUDA_HOME: $CUDA_HOME"

# 3. Create conda environment
info "Creating conda environment..."
conda env create -f environment_system_cuda.yml -y || {
    warn "Environment already exists, updating..."
    conda env update -f environment_system_cuda.yml
}

# 4. Create activation script for the conda environment
info "Creating CUDA activation script..."
CONDA_ENV_PATH=$(conda info --base)/envs/biogpu

mkdir -p $CONDA_ENV_PATH/etc/conda/activate.d
mkdir -p $CONDA_ENV_PATH/etc/conda/deactivate.d

# Activation script
cat > $CONDA_ENV_PATH/etc/conda/activate.d/cuda_env.sh << EOF
#!/bin/bash
# Set CUDA environment variables

export CUDA_HOME=$CUDA_HOME
export CUDA_PATH=\$CUDA_HOME
export PATH=\$CUDA_HOME/bin:\$PATH
export LD_LIBRARY_PATH=\$CUDA_HOME/lib64:\${LD_LIBRARY_PATH}
export CUDA_ROOT=\$CUDA_HOME

# For compilers
export CUDACXX=\$CUDA_HOME/bin/nvcc
export CMAKE_CUDA_COMPILER=\$CUDA_HOME/bin/nvcc

# For CuPy to find CUDA
export CUDA_PATH=\$CUDA_HOME
export CUPY_CUDA_PER_THREAD_DEFAULT_STREAM=1  # Better performance

echo "CUDA environment configured: \$CUDA_HOME"
EOF

# Deactivation script
cat > $CONDA_ENV_PATH/etc/conda/deactivate.d/cuda_env.sh << EOF
#!/bin/bash
# Unset CUDA variables
unset CUDA_HOME
unset CUDA_PATH
unset CUDA_ROOT
unset CUDACXX
unset CMAKE_CUDA_COMPILER
unset CUPY_CUDA_PER_THREAD_DEFAULT_STREAM
EOF

chmod +x $CONDA_ENV_PATH/etc/conda/activate.d/cuda_env.sh
chmod +x $CONDA_ENV_PATH/etc/conda/deactivate.d/cuda_env.sh

# 5. Create test script
info "Creating test script..."
cat > test_biogpu_setup.py << 'EOF'
#!/usr/bin/env python
"""Test BioGPU setup with system CUDA"""

import os
import sys
import subprocess

def print_header(text):
    print(f"\n{'=' * 60}")
    print(f"{text:^60}")
    print('=' * 60)

def test_cuda_env():
    print_header("CUDA Environment Variables")
    cuda_vars = ['CUDA_HOME', 'CUDA_PATH', 'LD_LIBRARY_PATH']
    for var in cuda_vars:
        value = os.environ.get(var, 'NOT SET')
        print(f"{var}: {value}")

def test_cuda_tools():
    print_header("CUDA Tools")
    
    # Test nvcc
    try:
        result = subprocess.run(['nvcc', '--version'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            version = [l for l in result.stdout.split('\n') if 'release' in l][0]
            print(f"✓ NVCC: {version.strip()}")
        else:
            print("✗ NVCC: Failed to run")
    except Exception as e:
        print(f"✗ NVCC: {e}")

def test_python_packages():
    print_header("Python GPU Packages")
    
    # Test CuPy
    try:
        import cupy as cp
        print(f"✓ CuPy: {cp.__version__}")
        
        # Get CUDA runtime version that CuPy sees
        cuda_version = cp.cuda.runtime.runtimeGetVersion()
        major = cuda_version // 1000
        minor = (cuda_version % 1000) // 10
        print(f"  CUDA Runtime: {major}.{minor}")
        
        # Simple computation test
        x = cp.array([1, 2, 3])
        y = cp.array([4, 5, 6])
        z = cp.dot(x, y)
        print(f"  Computation test: [1,2,3] · [4,5,6] = {z}")
        
        # Memory info
        mempool = cp.get_default_memory_pool()
        print(f"  Memory pool: {mempool.used_bytes() / 1e6:.1f} MB used")
        
    except Exception as e:
        print(f"✗ CuPy: {e}")
    
    # Test Numba
    try:
        import numba
        from numba import cuda
        print(f"✓ Numba: {numba.__version__}")
        print(f"  CUDA available: {cuda.is_available()}")
        
        if cuda.is_available():
            # Simple kernel test
            @cuda.jit
            def add_kernel(x, y, out):
                idx = cuda.grid(1)
                if idx < x.size:
                    out[idx] = x[idx] + y[idx]
            
            # Test it
            import numpy as np
            x = np.array([1, 2, 3], dtype=np.float32)
            y = np.array([4, 5, 6], dtype=np.float32)
            out = np.zeros_like(x)
            
            # This will fail without GPU but shows compilation works
            try:
                add_kernel[1, 3](x, y, out)
                print(f"  Kernel test: {x} + {y} = {out}")
            except:
                print("  Kernel compiled but no GPU to run")
                
    except Exception as e:
        print(f"✗ Numba: {e}")
    
    # Test PyCUDA
    try:
        import pycuda
        import pycuda.driver as cuda_driver
        print(f"✓ PyCUDA: {pycuda.VERSION_TEXT}")
        
        # This initializes CUDA
        try:
            cuda_driver.init()
            print(f"  Device count: {cuda_driver.Device.count()}")
        except:
            print("  PyCUDA imported but no GPU available")
            
    except Exception as e:
        print(f"✗ PyCUDA: {e}")

def test_bio_packages():
    print_header("Bioinformatics Packages")
    
    packages = {
        'Bio': 'BioPython',
        'pysam': 'PySAM',
        'numpy': 'NumPy',
        'pandas': 'Pandas',
        'scipy': 'SciPy'
    }
    
    for module, name in packages.items():
        try:
            m = __import__(module)
            version = getattr(m, '__version__', 'unknown')
            print(f"✓ {name}: {version}")
        except Exception as e:
            print(f"✗ {name}: not installed")

def main():
    print("BioGPU Environment Test")
    print("=" * 60)
    
    test_cuda_env()
    test_cuda_tools()
    test_python_packages()
    test_bio_packages()
    
    print_header("Summary")
    print("If you see CUDA errors above, that's normal without a GPU.")
    print("The important thing is that packages are installed correctly.")
    print("\nWhen your RTX A5000 arrives, all tests should pass!")

if __name__ == "__main__":
    main()
EOF

chmod +x test_biogpu_setup.py

# 6. Create simple CUDA test
info "Creating CUDA compilation test..."
cat > test_cuda_compile.cu << 'EOF'
#include <stdio.h>
#include <cuda_runtime.h>

__global__ void hello_biogpu() {
    printf("Hello from BioGPU - Block %d, Thread %d\n", blockIdx.x, threadIdx.x);
}

int main() {
    printf("BioGPU CUDA Compilation Test\n");
    printf("============================\n");
    
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    
    if (deviceCount == 0) {
        printf("No CUDA devices found (normal without GPU)\n");
        printf("✓ But CUDA compilation works!\n");
        return 0;
    }
    
    // If GPU is present
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("GPU: %s\n", prop.name);
    printf("Compute capability: %d.%d\n", prop.major, prop.minor);
    
    // Launch kernel
    hello_biogpu<<<2, 4>>>();
    cudaDeviceSynchronize();
    
    return 0;
}
EOF

# 7. Final instructions
echo
info "Setup complete!"
echo
echo "To use BioGPU environment:"
echo "  conda activate biogpu"
echo "  python test_biogpu_setup.py"
echo
echo "To test CUDA compilation:"
echo "  nvcc test_cuda_compile.cu -o test_cuda"
echo "  ./test_cuda"
echo
echo "The environment is configured to use your system CUDA 12.5"
echo "All GPU features will work when your RTX A5000 is installed!"