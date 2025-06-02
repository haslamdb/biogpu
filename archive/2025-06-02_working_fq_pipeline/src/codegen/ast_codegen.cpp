// src/codegen/ast_codegen.cpp
#include "biogpu/codegen/ast_codegen.h"
#include "biogpu/codegen/llvm_codegen.h"
#include "biogpu/parser/ast.h"

#include <llvm/IR/Module.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/IRBuilder.h>
#include <unordered_map>

namespace biogpu {

class ASTCodeGenerator::Impl {
public:
    LLVMCodeGenerator& codegen;
    std::unordered_map<std::string, llvm::Value*> symbolTable;
    
    Impl(LLVMCodeGenerator& cg) : codegen(cg) {}
    
    void generatePipeline(const Pipeline* pipeline) {
        // Generate a function for the entire pipeline
        auto pipelineName = "pipeline_" + pipeline->name;
        
        // For each stage in the pipeline
        for (const auto& stage : pipeline->stages) {
            generateStage(stage.get());
        }
        
        // Generate pipeline orchestration function
        generatePipelineOrchestrator(pipeline);
    }
    
    void generateStage(const Stage* stage) {
        if (!stage) return;
        
        // Check if this is a GPU kernel stage
        bool isGPUKernel = false;
        for (const auto& decorator : stage->decorators) {
            if (decorator == "@gpu_kernel") {
                isGPUKernel = true;
                break;
            }
        }
        
        if (isGPUKernel) {
            generateGPUStage(stage);
        } else {
            generateCPUStage(stage);
        }
    }
    
private:
    void generateGPUStage(const Stage* stage) {
        // Map stage operations to our predefined kernels
        for (const auto& stmt : stage->statements) {
            if (auto* call = dynamic_cast<const FunctionCall*>(stmt.get())) {
                generateKernelCall(call);
            }
        }
    }
    
    void generateCPUStage(const Stage* stage) {
        // Generate regular LLVM IR for CPU execution
        // This would handle report generation, I/O, etc.
    }
    
    void generateKernelCall(const FunctionCall* call) {
        const auto& funcName = call->functionName;
        
        // Map BioGPU functions to our LLVM kernels
        if (funcName == "filter_reads") {
            // Generate call to filter kernel
            generateFilterReadsCall(call);
        } else if (funcName == "scan_mutations") {
            // Use our resistance scanner kernel
            codegen.generateResistanceScanner();
        } else if (funcName == "parallel_map") {
            // Use our alignment kernel
            codegen.generateShortReadAligner();
        } else if (funcName == "calculate_abundance") {
            // Use our abundance calculator
            codegen.generateAbundanceCalculator();
        }
        // Add more mappings as needed
    }
    
    void generateFilterReadsCall(const FunctionCall* call) {
        // Extract parameters from the call
        // Generate appropriate LLVM IR
        
        // For now, generate a simple function
        // In reality, this would create the actual filtering kernel
    }
    
    void generatePipelineOrchestrator(const Pipeline* pipeline) {
        // Generate main function that coordinates all stages
        // Handle data flow between stages
        // Manage GPU memory allocation/deallocation
    }
};

ASTCodeGenerator::ASTCodeGenerator(LLVMCodeGenerator& codegen) 
    : pImpl(std::make_unique<Impl>(codegen)) {}

ASTCodeGenerator::~ASTCodeGenerator() = default;

void ASTCodeGenerator::generatePipeline(const Pipeline* pipeline) {
    pImpl->generatePipeline(pipeline);
}

void ASTCodeGenerator::generateStage(const Stage* stage) {
    pImpl->generateStage(stage);
}

} // namespace biogpu