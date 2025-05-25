# BioGPU Development Environment Setup

This guide will help you set up a complete development environment for BioGPU.

## Prerequisites

### 1. CUDA Toolkit Installation

BioGPU requires CUDA 12.0 or later for GPU acceleration.

#### Ubuntu/Debian
```bash
# Add NVIDIA package repositories
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update

# Install CUDA Toolkit
sudo apt-get -y install cuda-toolkit-12-3

# Add to PATH (add to ~/.bashrc for persistence)
export PATH=/usr/local/cuda-12.3/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-12.3/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

#### macOS
CUDA support on macOS has been discontinued. Consider using a Linux VM or cloud instance.

#### Windows
Download and run the installer from [NVIDIA CUDA Downloads](https://developer.nvidia.com/cuda-downloads)

#### Verify Installation
```bash
nvcc --version
nvidia-smi
```

### 2. LLVM Development Environment

BioGPU uses LLVM for its optimization pipeline.

#### Ubuntu/Debian
```bash
# Install LLVM 15
wget https://apt.llvm.org/llvm.sh
chmod +x llvm.sh
sudo ./llvm.sh 15

# Install development packages
sudo apt-get install llvm-15-dev libclang-15-dev clang-15
```

#### macOS (using Homebrew)
```bash
brew install llvm@15
export PATH="/usr/local/opt/llvm@15/bin:$PATH"
```

#### Windows
Use the pre-built binaries from [LLVM Releases](https://github.com/llvm/llvm-project/releases)

### 3. Build Tools

```bash
# Ubuntu/Debian
sudo apt-get install cmake build-essential git python3-dev

# macOS
brew install cmake

# Rust (for compiler implementation)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### 4. Python Environment (for prototype)

```bash
# Create virtual environment
python3 -m venv biogpu-env
source biogpu-env/bin/activate  # On Windows: biogpu-env\Scripts\activate

# Install dependencies
pip install --upgrade pip
pip install pycuda numpy biopython pytest black mypy
```

### 5. Additional Dependencies

```bash
# ANTLR4 (for parser generation)
cd /usr/local/lib
sudo wget https://www.antlr.org/download/antlr-4.13.1-complete.jar
export CLASSPATH=".:/usr/local/lib/antlr-4.13.1-complete.jar:$CLASSPATH"
alias antlr4='java -jar /usr/local/lib/antlr-4.13.1-complete.jar'

# Bio-specific tools (optional, for testing)
sudo apt-get install samtools bcftools
```

## Verification Script

Create a file `check_setup.py` to verify your environment:

```python
#!/usr/bin/env python3
import subprocess
import sys

def check_command(cmd, name):
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ {name} found")
            return True
        else:
            print(f"✗ {name} not found")
            return False
    except:
        print(f"✗ {name} not found")
        return False

def check_python_package(package):
    try:
        __import__(package)
        print(f"✓ Python package '{package}' found")
        return True
    except ImportError:
        print(f"✗ Python package '{package}' not found")
        return False

print("Checking BioGPU development environment...\n")

# Check CUDA
check_command("nvcc --version", "CUDA Compiler (nvcc)")
check_command("nvidia-smi", "NVIDIA GPU Driver")

# Check LLVM
check_command("llvm-config --version", "LLVM")
check_command("clang --version", "Clang")

# Check build tools
check_command("cmake --version", "CMake")
check_command("cargo --version", "Rust/Cargo")
check_command("antlr4", "ANTLR4")

# Check Python packages
print("\nPython packages:")
check_python_package("pycuda")
check_python_package("numpy")
check_python_package("Bio")

print("\nSetup verification complete!")
```

## IDE Setup

### VS Code (Recommended)
Install these extensions:
- C/C++ (Microsoft)
- CUDA C++ (NVIDIA)
- Rust Analyzer
- Python
- LLVM IR

### Settings.json for VS Code
```json
{
    "files.associations": {
        "*.cu": "cuda-cpp",
        "*.cuh": "cuda-cpp",
        "*.biogpu": "python"
    },
    "C_Cpp.default.includePath": [
        "/usr/local/cuda/include",
        "${workspaceFolder}/include"
    ]
}
```

## Testing Your Setup

Once everything is installed, test with a simple CUDA program:

```cuda
// test_cuda.cu
#include <stdio.h>

__global__ void hello_cuda() {
    printf("Hello from GPU thread %d!\n", threadIdx.x);
}

int main() {
    hello_cuda<<<1, 10>>>();
    cudaDeviceSynchronize();
    return 0;
}
```

Compile and run:
```bash
nvcc test_cuda.cu -o test_cuda
./test_cuda
```

## Troubleshooting

### CUDA Issues
- **No CUDA-capable device**: Ensure you have an NVIDIA GPU
- **Version mismatch**: Check that CUDA toolkit version matches your driver

### LLVM Issues
- **llvm-config not found**: Add LLVM bin directory to PATH
- **Missing headers**: Install llvm-dev package

### Python Issues
- **PyCUDA installation fails**: Ensure CUDA is properly installed first
- **Import errors**: Activate the virtual environment

## Next Steps

1. Clone the repository: `git clone https://github.com/haslamdb/biogpu.git`
2. Run the verification script
3. Start with the examples in the `examples/` directory
4. Check out the [Getting Started Guide](docs/getting-started.md)

## Questions?

Open an issue on [GitHub](https://github.com/haslamdb/biogpu/issues) if you encounter any problems!