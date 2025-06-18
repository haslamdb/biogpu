// include/biogpu/codegen/ast_codegen.h
#ifndef BIOGPU_AST_CODEGEN_H
#define BIOGPU_AST_CODEGEN_H

#include "biogpu/parser/ast.h"
#include <memory>

namespace llvm {
    class Module;
    class Value;
}

namespace biogpu {

class LLVMCodeGenerator;

/**
 * Translates BioGPU AST to LLVM IR
 */
class ASTCodeGenerator {
public:
    ASTCodeGenerator(LLVMCodeGenerator& codegen);
    ~ASTCodeGenerator();
    
    /**
     * Generate LLVM IR from a pipeline AST
     */
    void generatePipeline(const Pipeline* pipeline);
    
    /**
     * Generate IR for a stage
     */
    void generateStage(const Stage* stage);
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

} // namespace biogpu

#endif // BIOGPU_AST_CODEGEN_H