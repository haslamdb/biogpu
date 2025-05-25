#!/bin/bash
# scripts/build_grammar.sh
# Build ANTLR grammar for BioGPU

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"
PARSER_DIR=~/Documents/Code/biogpu/src/parser

echo "Building BioGPU grammar..."

# Check if ANTLR is installed
if ! command -v antlr4 &> /dev/null; then
    echo "Error: antlr4 not found. Install it with:"
    echo "  Ubuntu: sudo apt-get install antlr4"
    echo "  Mac: brew install antlr"
    echo "  Or download from: https://www.antlr.org/download.html"
    exit 1
fi

# Navigate to parser directory
cd ~/Documents/Code/biogpu/src/parser

# Clean old generated files
echo "Cleaning old generated files..."
rm -f BioGPULexer.* BioGPUParser.* BioGPUVisitor.* BioGPUBaseVisitor.*

# Generate C++ parser
echo "Generating C++ parser from BioGPU.g4..."
antlr4 -Dlanguage=Cpp -visitor -no-listener -o . BioGPU.g4

# Check if generation succeeded
if [ $? -eq 0 ]; then
    echo "✓ Grammar build successful!"
    echo "Generated files:"
    ls -la BioGPU*.{h,cpp} 2>/dev/null || true
else
    echo "✗ Grammar build failed!"
    exit 1
fi

# Create a simple test
cat > test_parser.cpp << 'EOF'
#include <iostream>
#include <fstream>
#include "antlr4-runtime.h"
#include "BioGPULexer.h"
#include "BioGPUParser.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.biogpu>" << std::endl;
        return 1;
    }
    
    std::ifstream stream(argv[1]);
    antlr4::ANTLRInputStream input(stream);
    BioGPULexer lexer(&input);
    antlr4::CommonTokenStream tokens(&lexer);
    BioGPUParser parser(&tokens);
    
    auto tree = parser.program();
    std::cout << "Parse tree: " << tree->toStringTree(&parser) << std::endl;
    
    return 0;
}
EOF

echo ""
echo "To test the parser:"
echo "  g++ -o test_parser test_parser.cpp BioGPU*.cpp -lantlr4-runtime"
echo "  ./test_parser ../../examples/simple_fq_resistance.biogpu"