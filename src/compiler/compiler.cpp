// src/compiler/compiler.cpp
#include "biogpu/parser/parser.h"
#include "biogpu/codegen/llvm_codegen.h"
#include "biogpu/codegen/ast_codegen.h"

#include <iostream>
#include <fstream>
#include <memory>

namespace biogpu {

class Compiler {
public:
    Compiler() : parser(), llvmGen(), astGen(llvmGen) {}
    
    bool compileFile(const std::string& inputFile, const std::string& outputFile) {
        std::cout << "Compiling " << inputFile << "..." << std::endl;
        
        // Step 1: Parse
        std::cout << "  Parsing..." << std::endl;
        auto ast = parser.parseFile(inputFile);
        
        if (!ast) {
            std::cerr << "Parse errors:" << std::endl;
            for (const auto& error : parser.getErrors()) {
                std::cerr << "  Line " << error.line << ":" << error.column 
                         << " - " << error.message << std::endl;
            }
            return false;
        }
        
        // Step 2: Generate LLVM IR
        std::cout << "  Generating LLVM IR..." << std::endl;
        
        // First generate the built-in kernels
        llvmGen.generateBioKernels();
        
        // Then generate code from the AST
        astGen.generatePipeline(ast.get());
        
        // Step 3: Verify
        if (!llvmGen.verify()) {
            std::cerr << "LLVM module verification failed!" << std::endl;
            return false;
        }
        
        // Step 4: Write output
        std::cout << "  Writing output to " << outputFile << "..." << std::endl;
        std::error_code EC;
        llvm::raw_fd_ostream fileStream(outputFile, EC);
        if (EC) {
            std::cerr << "Error opening output file: " << EC.message() << std::endl;
            return false;
        }
        llvmGen.getModule()->print(fileStream, nullptr);
        
        std::cout << "✓ Compilation successful!" << std::endl;
        return true;
    }
    
    bool run(const std::string& inputFile) {
        // Compile to temporary file
        std::string irFile = inputFile + ".ll";
        if (!compileFile(inputFile, irFile)) {
            return false;
        }
        
        // In the future, this would:
        // 1. Compile IR to PTX
        // 2. Load and execute on GPU
        // 3. Return results
        
        std::cout << "Execution not yet implemented" << std::endl;
        return true;
    }
    
private:
    BioGPUParser parser;
    LLVMCodeGenerator llvmGen;
    ASTCodeGenerator astGen;
};

} // namespace biogpu

// Main entry point
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <command> [options]" << std::endl;
        std::cerr << "Commands:" << std::endl;
        std::cerr << "  compile <input.biogpu> [output.ll]  - Compile to LLVM IR" << std::endl;
        std::cerr << "  run <input.biogpu>                   - Compile and run" << std::endl;
        std::cerr << "  test                                 - Run built-in tests" << std::endl;
        return 1;
    }
    
    biogpu::Compiler compiler;
    std::string command = argv[1];
    
    if (command == "compile") {
        if (argc < 3) {
            std::cerr << "Error: Missing input file" << std::endl;
            return 1;
        }
        
        std::string inputFile = argv[2];
        std::string outputFile = (argc > 3) ? argv[3] : inputFile + ".ll";
        
        return compiler.compileFile(inputFile, outputFile) ? 0 : 1;
        
    } else if (command == "run") {
        if (argc < 3) {
            std::cerr << "Error: Missing input file" << std::endl;
            return 1;
        }
        
        return compiler.run(argv[2]) ? 0 : 1;
        
    } else if (command == "test") {
        // Run built-in test
        std::cout << "Running built-in tests..." << std::endl;
        
        // Test parsing
        const char* testProgram = R"(
            pipeline TestPipeline {
                input: fastq_file reads
                output: json results
                
                @gpu_kernel
                stage detect_resistance {
                    mutations = scan_mutations(reads) {
                        database: "fq_mutations",
                        min_quality: 30
                    }
                    emit: mutations
                }
            }
        )";
        
        biogpu::BioGPUParser parser;
        auto ast = parser.parseString(testProgram);
        
        if (ast) {
            std::cout << "✓ Parse test passed" << std::endl;
        } else {
            std::cout << "✗ Parse test failed" << std::endl;
        }
        
        return 0;
        
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        return 1;
    }
}