#include <iostream>
#include <fstream>
#include <system_error>
#include <llvm/Support/raw_ostream.h>
#include <llvm/IR/Module.h>
#include "biogpu/codegen/llvm_codegen.h"  // Use header instead of .cpp

int main(int argc, char* argv[]) {
    std::cout << "BioGPU Compiler v0.1.0" << std::endl;
    std::cout << "======================\n" << std::endl;
    
    // Test LLVM code generation
    biogpu::LLVMCodeGenerator codegen;
    
    std::cout << "Generating bioinformatics kernels..." << std::endl;
    
    // Generate all kernels
    codegen.generateBioKernels();
    
    // Verify module
    if (codegen.verify()) {
        std::cout << "✓ Module verified successfully!\n" << std::endl;
        
        // Print generated IR
        if (argc > 1 && std::string(argv[1]) == "--print-ir") {
            std::cout << "Generated LLVM IR:" << std::endl;
            std::cout << "==================" << std::endl;
            codegen.print();
        } else {
            // Save to file
            std::string filename = "biogpu_generated.ll";
            std::error_code EC;
            llvm::raw_fd_ostream fileStream(filename, EC);
            
            if (!EC) {
                codegen.getModule()->print(fileStream, nullptr);
                std::cout << "Generated IR saved to: " << filename << std::endl;
            } else {
                std::cerr << "Error opening file: " << EC.message() << std::endl;
            }
            
            // Print summary
            std::cout << "\nGenerated functions:" << std::endl;
            for (auto& F : codegen.getModule()->functions()) {
                if (!F.isDeclaration()) {
                    std::cout << "  - " << F.getName().str();
                    if (F.hasFnAttribute("nvvm.annotations")) {
                        std::cout << " [GPU kernel]";
                    }
                    std::cout << std::endl;
                }
            }
        }
        
        // Generate PTX if requested
        if (argc > 1 && std::string(argv[1]) == "--gen-ptx") {
            std::cout << "\nGenerating PTX assembly..." << std::endl;
            std::string ptx = codegen.generatePTX();
            
            std::ofstream ptxFile("biogpu_kernels.ptx");
            ptxFile << ptx;
            ptxFile.close();
            
            std::cout << "PTX saved to: biogpu_kernels.ptx" << std::endl;
        }
        
    } else {
        std::cerr << "✗ Module verification failed!" << std::endl;
        return 1;
    }
    
    std::cout << "\nNext steps:" << std::endl;
    std::cout << "1. Compile LLVM IR to PTX: llc -mcpu=sm_75 biogpu_generated.ll -o biogpu.ptx" << std::endl;
    std::cout << "2. Run tests: ./build.sh test" << std::endl;
    std::cout << "3. Try examples: biogpu run examples/simple_fq_resistance.biogpu" << std::endl;
    
    return 0;
}