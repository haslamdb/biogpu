#ifndef BIOGPU_LLVM_CODEGEN_H
#define BIOGPU_LLVM_CODEGEN_H

#include <memory>
#include <string>
#include <vector>
#include <llvm/IR/IRBuilder.h>

// Forward declarations for LLVM types
namespace llvm {
    class LLVMContext;
    class Module;
    class StructType;
    class Function;
    class Value;
    class Type;
}

namespace biogpu {

/**
 * @brief LLVM Code Generator for BioGPU
 * 
 * This class generates LLVM IR for bioinformatics kernels that can be
 * compiled to CUDA PTX for GPU execution.
 */
class LLVMCodeGenerator {
private:
    std::unique_ptr<llvm::LLVMContext> context;
    std::unique_ptr<llvm::Module> module;
    std::unique_ptr<llvm::IRBuilder<>> builder;
    
    // Custom types for bioinformatics
    llvm::StructType* dnaSequenceType;
    llvm::StructType* genomeDBType;
    llvm::StructType* kmerIndexType;
    llvm::StructType* alignmentResultType;
    llvm::StructType* resistanceMutationType;

    // Private helper methods
    void initializeTypes();
    void generateIntrinsics();
    llvm::Value* getGlobalThreadId();

public:
    /**
     * @brief Construct a new LLVMCodeGenerator
     */
    LLVMCodeGenerator();
    
    /**
     * @brief Destructor
     */
    ~LLVMCodeGenerator();
    
    // Prevent copying
    LLVMCodeGenerator(const LLVMCodeGenerator&) = delete;
    LLVMCodeGenerator& operator=(const LLVMCodeGenerator&) = delete;
    
    // Allow moving
    LLVMCodeGenerator(LLVMCodeGenerator&&) = default;
    LLVMCodeGenerator& operator=(LLVMCodeGenerator&&) = default;

    /**
     * @brief Generate CUDA kernel launch function
     * @param kernelName Name of the kernel
     * @param blockSize CUDA block size
     * @param gridSize CUDA grid size
     * @return Generated LLVM Function
     */
    llvm::Function* generateKernelLauncher(
        const std::string& kernelName,
        int blockSize,
        int gridSize
    );
    
    /**
     * @brief Generate k-mer counting kernel
     * @return Generated LLVM Function
     */
    llvm::Function* generateKmerCounter();
    
    /**
     * @brief Generate short read alignment kernel (Bowtie2-style)
     * @return Generated LLVM Function
     */
    llvm::Function* generateShortReadAligner();
    
    /**
     * @brief Generate resistance mutation scanning kernel
     * @return Generated LLVM Function
     */
    llvm::Function* generateResistanceScanner();
    
    /**
     * @brief Generate microbial abundance calculation kernel
     * @return Generated LLVM Function
     */
    llvm::Function* generateAbundanceCalculator();
    
    /**
     * @brief Generate optimized Smith-Waterman kernel for GPU
     * @return Generated LLVM Function
     */
    llvm::Function* generateSmithWatermanGPU();
    
    /**
     * @brief Generate basic Smith-Waterman implementation
     * @return Generated LLVM Function
     */
    llvm::Function* generateSmithWaterman();
    
    /**
     * @brief Generate all bioinformatics kernels
     */
    void generateBioKernels();
    
    /**
     * @brief Print the generated LLVM IR to stdout
     */
    void print();
    
    /**
     * @brief Verify the generated module
     * @return true if module is valid, false otherwise
     */
    bool verify();
    
    /**
     * @brief Generate PTX assembly from LLVM IR
     * @return PTX assembly as string
     */
    std::string generatePTX();
    
    /**
     * @brief Get the LLVM Module
     * @return Pointer to the module (not owned)
     */
    llvm::Module* getModule() { return module.get(); }
};

} // namespace biogpu

#endif // BIOGPU_LLVM_CODEGEN_H
