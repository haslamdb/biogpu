#!/bin/bash
# Script to reorganize BioGPU project structure

echo "Reorganizing BioGPU project structure..."

# Create include directory structure
mkdir -p include/biogpu/{codegen,kernels,runtime,types}

# Create header file from the artifact (you'll need to copy the content)
cat > include/biogpu/codegen/llvm_codegen.h << 'EOF'
#ifndef BIOGPU_LLVM_CODEGEN_H
#define BIOGPU_LLVM_CODEGEN_H

#include <memory>
#include <string>
#include <vector>

// Forward declarations for LLVM types
namespace llvm {
    class LLVMContext;
    class Module;
    template<typename T> class IRBuilder;
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
EOF

# Create additional header files
cat > include/biogpu/types/bio_types.h << 'EOF'
#ifndef BIOGPU_BIO_TYPES_H
#define BIOGPU_BIO_TYPES_H

#include <cstdint>
#include <string>

namespace biogpu {

// Forward declarations
struct DNASequence;
struct KmerIndex;
struct AlignmentResult;
struct ResistanceMutation;

// DNA encoding
enum class DNABase : uint8_t {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4  // Unknown
};

// Quality score type
using QualityScore = uint8_t;

// Sequence ID type
using SequenceID = uint32_t;

} // namespace biogpu

#endif // BIOGPU_BIO_TYPES_H
EOF

# Create a kernel header as an example
cat > include/biogpu/kernels/kmer.h << 'EOF'
#ifndef BIOGPU_KERNELS_KMER_H
#define BIOGPU_KERNELS_KMER_H

namespace biogpu {
namespace kernels {

// K-mer related kernel declarations would go here

} // namespace kernels
} // namespace biogpu

#endif // BIOGPU_KERNELS_KMER_H
EOF

# Now you need to update the implementation file
# Remove the old full implementation and create a clean one
# This should be done manually to preserve your specific implementation

echo "Directory structure created!"
echo ""
echo "Next steps:"
echo "1. Copy the implementation from the artifact to src/codegen/llvm_codegen.cpp"
echo "2. Update src/main.cpp to use the header file"
echo "3. Run: mkdir build && cd build && cmake .. && make"
echo ""
echo "The separation allows:"
echo "- Faster compilation (headers change less often)"
echo "- Better organization"
echo "- Easier testing (can mock interfaces)"
echo "- Clear API documentation"