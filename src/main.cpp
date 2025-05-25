#include <iostream>
#include "llvm_codegen.cpp"  // In practice, use proper headers

int main() {
    std::cout << "BioGPU Compiler v0.1.0" << std::endl;
    
    // Test LLVM code generation
    biogpu::LLVMCodeGenerator codegen;
    
    // Generate some functions
    codegen.generateSmithWaterman();
    codegen.generateKernelLauncher("map_sequences", 256, 1024);
    
    // Verify and print
    if (codegen.verify()) {
        std::cout << "Module verified successfully!" << std::endl;
        codegen.print();
    } else {
        std::cerr << "Module verification failed!" << std::endl;
        return 1;
    }
    
    return 0;
}
