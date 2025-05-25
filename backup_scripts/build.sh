#!/bin/bash
# BioGPU Development Build Script

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Project directories
PROJECT_ROOT=$(pwd)
BUILD_DIR="${PROJECT_ROOT}/build"
DATA_DIR="${PROJECT_ROOT}/data"
TESTS_DIR="${PROJECT_ROOT}/tests"

echo -e "${BLUE}BioGPU Build System${NC}"
echo "===================="

# Function to check dependencies
check_dependencies() {
    echo -e "\n${YELLOW}Checking dependencies...${NC}"
    
    # Check CUDA
    if command -v nvcc &> /dev/null; then
        CUDA_VERSION=$(nvcc --version | grep "release" | awk '{print $6}' | cut -c2-)
        echo -e "${GREEN}✓ CUDA found: ${CUDA_VERSION}${NC}"
    else
        echo -e "${RED}✗ CUDA not found${NC}"
        echo "  Please install CUDA toolkit: https://developer.nvidia.com/cuda-downloads"
    fi
    
    # Check LLVM
    if command -v llvm-config &> /dev/null; then
        LLVM_VERSION=$(llvm-config --version)
        echo -e "${GREEN}✓ LLVM found: ${LLVM_VERSION}${NC}"
    else
        echo -e "${RED}✗ LLVM not found${NC}"
        echo "  Please install LLVM development packages"
    fi
    
    # Check CMake
    if command -v cmake &> /dev/null; then
        CMAKE_VERSION=$(cmake --version | head -n1 | awk '{print $3}')
        echo -e "${GREEN}✓ CMake found: ${CMAKE_VERSION}${NC}"
    else
        echo -e "${RED}✗ CMake not found${NC}"
    fi
    
    # Check Python environment
    if command -v python3 &> /dev/null; then
        PYTHON_VERSION=$(python3 --version | awk '{print $2}')
        echo -e "${GREEN}✓ Python found: ${PYTHON_VERSION}${NC}"
        
        # Check for CuPy
        if python3 -c "import cupy" 2>/dev/null; then
            echo -e "${GREEN}✓ CuPy installed${NC}"
        else
            echo -e "${YELLOW}! CuPy not installed (optional for GPU)${NC}"
        fi
    else
        echo -e "${RED}✗ Python 3 not found${NC}"
    fi
}

# Function to build the project
build_project() {
    echo -e "\n${YELLOW}Building BioGPU...${NC}"
    
    # Create build directory
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    # Configure with CMake
    echo "Configuring with CMake..."
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCUDA_ARCH_LIST="60;70;75;80;86" \
             -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    
    # Build
    echo "Building..."
    make -j$(nproc)
    
    echo -e "${GREEN}✓ Build complete${NC}"
    cd "$PROJECT_ROOT"
}

# Function to run tests
run_tests() {
    echo -e "\n${YELLOW}Running tests...${NC}"
    
    # Run C++ tests
    if [ -f "$BUILD_DIR/biogpu" ]; then
        echo "Running LLVM IR generation test..."
        "$BUILD_DIR/biogpu"
    fi
    
    # Run Python tests
    echo -e "\nRunning Python tests..."
    python3 tests/test_biogpu_setup.py
    
    # Run algorithm tests
    echo -e "\nRunning algorithm tests..."
    python3 tests/test_algorithms.py
    
    # Check for GPU and run CUDA tests
    if nvidia-smi &> /dev/null; then
        echo -e "\nRunning GPU tests..."
        if [ -f "$BUILD_DIR/tests/qrdr_detector" ]; then
            "$BUILD_DIR/tests/qrdr_detector"
        fi
    else
        echo -e "${YELLOW}! No GPU detected, skipping GPU tests${NC}"
    fi
}

# Function to download reference data
download_data() {
    echo -e "\n${YELLOW}Setting up reference data...${NC}"
    
    mkdir -p "$DATA_DIR"/{genomes,resistance,test}
    
    # Create example resistance mutation database
    cat > "$DATA_DIR/resistance/fluoroquinolone_mutations.json" << 'EOF'
{
  "mutations": [
    {
      "gene": "gyrA",
      "organism": "Escherichia coli",
      "position": 83,
      "wild_type": "S",
      "mutations": ["L", "W"],
      "mic_change": 8.0,
      "drugs": ["ciprofloxacin", "levofloxacin"]
    },
    {
      "gene": "gyrA",
      "organism": "Escherichia coli",
      "position": 87,
      "wild_type": "D",
      "mutations": ["N", "G", "Y"],
      "mic_change": 4.0,
      "drugs": ["ciprofloxacin", "levofloxacin"]
    },
    {
      "gene": "parC",
      "organism": "Escherichia coli",
      "position": 80,
      "wild_type": "S",
      "mutations": ["I", "R"],
      "mic_change": 4.0,
      "drugs": ["ciprofloxacin", "levofloxacin"]
    },
    {
      "gene": "gyrA",
      "organism": "Klebsiella pneumoniae",
      "position": 83,
      "wild_type": "S",
      "mutations": ["L", "F"],
      "mic_change": 16.0,
      "drugs": ["ciprofloxacin", "levofloxacin", "moxifloxacin"]
    }
  ]
}
EOF
    
    echo -e "${GREEN}✓ Created example mutation database${NC}"
    
    # Download test data if needed
    if [ ! -f "$DATA_DIR/test/test_reads.fastq" ]; then
        echo "Creating test FASTQ file..."
        python3 -c "
import random
bases = 'ACGT'
with open('$DATA_DIR/test/test_reads.fastq', 'w') as f:
    for i in range(1000):
        seq = ''.join(random.choice(bases) for _ in range(150))
        qual = ''.join('I' for _ in range(150))
        f.write(f'@read_{i}\n{seq}\n+\n{qual}\n')
"
        echo -e "${GREEN}✓ Created test FASTQ file${NC}"
    fi
}

# Function to generate documentation
generate_docs() {
    echo -e "\n${YELLOW}Generating documentation...${NC}"
    
    # Generate LLVM IR documentation
    if [ -f "$BUILD_DIR/biogpu" ]; then
        "$BUILD_DIR/biogpu" > "$PROJECT_ROOT/docs/generated_ir.ll"
        echo -e "${GREEN}✓ Generated LLVM IR documentation${NC}"
    fi
    
    # TODO: Add Doxygen/Sphinx generation
}

# Function to setup development environment
setup_dev_env() {
    echo -e "\n${YELLOW}Setting up development environment...${NC}"
    
    # Create Python virtual environment if it doesn't exist
    if [ ! -d "venv" ]; then
        echo "Creating Python virtual environment..."
        python3 -m venv venv
        source venv/bin/activate
        pip install --upgrade pip
        pip install numpy biopython pytest black mypy
        
        # Try to install GPU packages
        pip install cupy-cuda12x 2>/dev/null || echo "CuPy installation failed (normal without GPU)"
        pip install pycuda 2>/dev/null || echo "PyCUDA installation failed (normal without GPU)"
        
        echo -e "${GREEN}✓ Virtual environment created${NC}"
    else
        echo "Virtual environment already exists"
    fi
}

# Function to clean build artifacts
clean_build() {
    echo -e "\n${YELLOW}Cleaning build artifacts...${NC}"
    rm -rf "$BUILD_DIR"
    find . -name "*.pyc" -delete
    find . -name "__pycache__" -delete
    echo -e "${GREEN}✓ Clean complete${NC}"
}

# Main menu
show_menu() {
    echo -e "\n${BLUE}BioGPU Development Menu${NC}"
    echo "1) Check dependencies"
    echo "2) Build project"
    echo "3) Run tests"
    echo "4) Download/setup data"
    echo "5) Generate documentation"
    echo "6) Setup development environment"
    echo "7) Clean build"
    echo "8) Full build (all of the above)"
    echo "9) Exit"
}

# Process command line arguments
if [ $# -eq 0 ]; then
    # Interactive mode
    while true; do
        show_menu
        read -p "Select option: " choice
        
        case $choice in
            1) check_dependencies ;;
            2) build_project ;;
            3) run_tests ;;
            4) download_data ;;
            5) generate_docs ;;
            6) setup_dev_env ;;
            7) clean_build ;;
            8) 
                check_dependencies
                setup_dev_env
                build_project
                download_data
                run_tests
                generate_docs
                ;;
            9) exit 0 ;;
            *) echo "Invalid option" ;;
        esac
    done
else
    # Command line mode
    case $1 in
        check) check_dependencies ;;
        build) build_project ;;
        test) run_tests ;;
        data) download_data ;;
        docs) generate_docs ;;
        setup) setup_dev_env ;;
        clean) clean_build ;;
        all)
            check_dependencies
            setup_dev_env
            build_project
            download_data
            run_tests
            generate_docs
            ;;
        *)
            echo "Usage: $0 [check|build|test|data|docs|setup|clean|all]"
            echo "  or run without arguments for interactive mode"
            exit 1
            ;;
    esac
fi

echo -e "\n${GREEN}Done!${NC}"