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
