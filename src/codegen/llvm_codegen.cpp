#include "biogpu/codegen/llvm_codegen.h"

#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Verifier.h>
#include <llvm/IR/DerivedTypes.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/BasicBlock.h>
#include <llvm/IR/Constants.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Target/TargetMachine.h>

namespace biogpu {

LLVMCodeGenerator::LLVMCodeGenerator() {
    // Initialize LLVM
    llvm::InitializeAllTargetInfos();
    llvm::InitializeAllTargets();
    llvm::InitializeAllTargetMCs();
    llvm::InitializeAllAsmPrinters();
    
    context = std::make_unique<llvm::LLVMContext>();
    module = std::make_unique<llvm::Module>("biogpu_module", *context);
    builder = std::make_unique<llvm::IRBuilder<>>(*context);
    
    // Set target triple for NVPTX
    module->setTargetTriple("nvptx64-nvidia-cuda");
    
    // Initialize custom types
    initializeTypes();
    
    // Generate intrinsic declarations
    generateIntrinsics();
}

LLVMCodeGenerator::~LLVMCodeGenerator() = default;

void LLVMCodeGenerator::initializeTypes() {
    // DNA Sequence type: { i8*, i32, i8*, i32 }
    // Fields: data, length, quality, read_id
    dnaSequenceType = llvm::StructType::create(*context, "DNASequence");
    std::vector<llvm::Type*> dnaFields = {
        llvm::PointerType::getUnqual(llvm::Type::getInt8Ty(*context)),  // sequence data
        llvm::Type::getInt32Ty(*context),    // length
        llvm::PointerType::getUnqual(llvm::Type::getInt8Ty(*context)),  // quality scores
        llvm::Type::getInt32Ty(*context)     // read_id
    };
    dnaSequenceType->setBody(dnaFields);
    
    // K-mer Index type: { i64*, i32*, i32, i32 }
    // Fields: kmer_hashes, positions, num_kmers, k
    kmerIndexType = llvm::StructType::create(*context, "KmerIndex");
    std::vector<llvm::Type*> kmerFields = {
        llvm::PointerType::getUnqual(llvm::Type::getInt64Ty(*context)), // kmer hashes
        llvm::PointerType::getUnqual(llvm::Type::getInt32Ty(*context)), // positions
        llvm::Type::getInt32Ty(*context),    // num_kmers
        llvm::Type::getInt32Ty(*context)     // k value
    };
    kmerIndexType->setBody(kmerFields);
    
    // Alignment Result type
    alignmentResultType = llvm::StructType::create(*context, "AlignmentResult");
    std::vector<llvm::Type*> alignFields = {
        llvm::Type::getInt32Ty(*context),    // reference_id
        llvm::Type::getInt32Ty(*context),    // query_start
        llvm::Type::getInt32Ty(*context),    // reference_start
        llvm::Type::getInt32Ty(*context),    // alignment_length
        llvm::Type::getFloatTy(*context),    // score
        llvm::Type::getInt8Ty(*context)      // strand (0=forward, 1=reverse)
    };
    alignmentResultType->setBody(alignFields);
    
    // Resistance Mutation type
    resistanceMutationType = llvm::StructType::create(*context, "ResistanceMutation");
    std::vector<llvm::Type*> mutationFields = {
        llvm::Type::getInt32Ty(*context),    // gene_id
        llvm::Type::getInt32Ty(*context),    // position
        llvm::Type::getInt8Ty(*context),     // wild_type_aa
        llvm::Type::getInt8Ty(*context),     // mutant_aa
        llvm::Type::getFloatTy(*context)     // resistance_level
    };
    resistanceMutationType->setBody(mutationFields);
    
    // Genome Database type
    genomeDBType = llvm::StructType::create(*context, "GenomeDB");
    std::vector<llvm::Type*> genomeFields = {
        llvm::PointerType::get(dnaSequenceType, 0),  // sequences array
        llvm::Type::getInt32Ty(*context),            // num sequences
        llvm::PointerType::get(kmerIndexType, 0),    // k-mer index
        llvm::PointerType::getUnqual(llvm::Type::getInt8Ty(*context))  // metadata
    };
    genomeDBType->setBody(genomeFields);
}

void LLVMCodeGenerator::generateIntrinsics() {
    // Declare CUDA intrinsics we'll use
    // threadIdx.x
    llvm::FunctionType* tidxType = llvm::FunctionType::get(
        llvm::Type::getInt32Ty(*context), false);
    llvm::Function::Create(tidxType, llvm::Function::ExternalLinkage,
        "llvm.nvvm.read.ptx.sreg.tid.x", module.get());
    
    // blockIdx.x
    llvm::FunctionType* bidxType = llvm::FunctionType::get(
        llvm::Type::getInt32Ty(*context), false);
    llvm::Function::Create(bidxType, llvm::Function::ExternalLinkage,
        "llvm.nvvm.read.ptx.sreg.ctaid.x", module.get());
    
    // blockDim.x
    llvm::FunctionType* bdimType = llvm::FunctionType::get(
        llvm::Type::getInt32Ty(*context), false);
    llvm::Function::Create(bdimType, llvm::Function::ExternalLinkage,
        "llvm.nvvm.read.ptx.sreg.ntid.x", module.get());
}

llvm::Value* LLVMCodeGenerator::getGlobalThreadId() {
    auto tid_x = builder->CreateCall(
        module->getFunction("llvm.nvvm.read.ptx.sreg.tid.x"));
    auto bid_x = builder->CreateCall(
        module->getFunction("llvm.nvvm.read.ptx.sreg.ctaid.x"));
    auto bdim_x = builder->CreateCall(
        module->getFunction("llvm.nvvm.read.ptx.sreg.ntid.x"));
    
    return builder->CreateAdd(tid_x, builder->CreateMul(bid_x, bdim_x));
}

llvm::Function* LLVMCodeGenerator::generateKernelLauncher(
    const std::string& kernelName,
    int blockSize,
    int gridSize
) {
    // Function type: void launch_kernel(DNASequence*, GenomeDB*, float*)
    std::vector<llvm::Type*> paramTypes = {
        llvm::PointerType::get(dnaSequenceType, 0),
        llvm::PointerType::get(genomeDBType, 0),
        llvm::PointerType::getUnqual(llvm::Type::getFloatTy(*context))
    };
    
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getVoidTy(*context),
        paramTypes,
        false
    );
    
    auto func = llvm::Function::Create(
        funcType,
        llvm::Function::ExternalLinkage,
        "launch_" + kernelName,
        module.get()
    );
    
    // Create entry block
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    
    // TODO: Generate CUDA kernel launch code
    // This would include:
    // 1. Allocate GPU memory
    // 2. Copy data to GPU
    // 3. Launch kernel
    // 4. Copy results back
    
    builder->CreateRetVoid();
    
    return func;
}

llvm::Function* LLVMCodeGenerator::generateSmithWaterman() {
    // Function signature: i32 smith_waterman(i8*, i32, i8*, i32, i32, i32, i32)
    std::vector<llvm::Type*> swParamTypes = {
        llvm::PointerType::getUnqual(llvm::Type::getInt8Ty(*context)),  // seq1
        llvm::Type::getInt32Ty(*context),    // len1
        llvm::PointerType::getUnqual(llvm::Type::getInt8Ty(*context)),  // seq2
        llvm::Type::getInt32Ty(*context),    // len2
        llvm::Type::getInt32Ty(*context),    // match_score
        llvm::Type::getInt32Ty(*context),    // mismatch_score
        llvm::Type::getInt32Ty(*context)     // gap_penalty
    };
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getInt32Ty(*context),
        swParamTypes,
        false
    );
    
    auto func = llvm::Function::Create(
        funcType,
        llvm::Function::ExternalLinkage,
        "smith_waterman",
        module.get()
    );
    
    // Implementation would go here
    // For now, just return 0
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    builder->CreateRet(llvm::ConstantInt::get(*context, llvm::APInt(32, 0)));
    
    return func;
}

llvm::Function* LLVMCodeGenerator::generateKmerCounter() {
    // __global__ void count_kmers(DNASequence* seqs, int num_seqs, 
    //                            int k, uint64_t* kmer_counts)
    std::vector<llvm::Type*> paramTypes = {
        llvm::PointerType::get(dnaSequenceType, 1),     // sequences (global)
        llvm::Type::getInt32Ty(*context),               // num_sequences
        llvm::Type::getInt32Ty(*context),               // k
        llvm::PointerType::get(llvm::Type::getInt64Ty(*context), 1)  // output counts
    };
    
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getVoidTy(*context), paramTypes, false);
    
    auto func = llvm::Function::Create(
        funcType, llvm::Function::ExternalLinkage,
        "count_kmers_kernel", module.get());
    
    // Set NVVM annotations for kernel
    func->addFnAttr("nvvm.annotations", "kernel");
    
    // Create entry block
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    
    // Get thread/block indices
    auto tid_x = builder->CreateCall(
        module->getFunction("llvm.nvvm.read.ptx.sreg.tid.x"));
    auto bid_x = builder->CreateCall(
        module->getFunction("llvm.nvvm.read.ptx.sreg.ctaid.x"));
    auto bdim_x = builder->CreateCall(
        module->getFunction("llvm.nvvm.read.ptx.sreg.ntid.x"));
    
    // Calculate global thread ID
    auto globalId = builder->CreateAdd(
        tid_x, 
        builder->CreateMul(bid_x, bdim_x)
    );
    
    // Get function arguments
    auto args = func->args().begin();
    llvm::Value* sequences = &*args++;
    llvm::Value* numSeqs = &*args++;
    llvm::Value* k = &*args++;
    llvm::Value* kmerCounts = &*args++;
    
    // Check bounds
    auto cmp = builder->CreateICmpSLT(globalId, numSeqs);
    auto thenBB = llvm::BasicBlock::Create(*context, "then", func);
    auto endBB = llvm::BasicBlock::Create(*context, "end", func);
    
    builder->CreateCondBr(cmp, thenBB, endBB);
    
    // Process sequence
    builder->SetInsertPoint(thenBB);
    
    // Get sequence pointer
    auto seqPtr = builder->CreateGEP(dnaSequenceType, sequences, globalId);
    
    // TODO: Implement k-mer extraction and hashing
    // For now, just placeholder
    
    builder->CreateBr(endBB);
    
    // End block
    builder->SetInsertPoint(endBB);
    builder->CreateRetVoid();
    
    return func;
}

llvm::Function* LLVMCodeGenerator::generateShortReadAligner() {
    // Simplified FM-index based alignment
    std::vector<llvm::Type*> paramTypes = {
        llvm::PointerType::get(dnaSequenceType, 1),          // reads
        llvm::Type::getInt32Ty(*context),                    // num_reads
        llvm::PointerType::get(genomeDBType, 1),             // reference
        llvm::PointerType::get(alignmentResultType, 1),      // results
        llvm::Type::getInt32Ty(*context)                     // min_score
    };
    
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getVoidTy(*context), paramTypes, false);
    
    auto func = llvm::Function::Create(
        funcType, llvm::Function::ExternalLinkage,
        "align_short_reads_kernel", module.get());
    
    func->addFnAttr("nvvm.annotations", "kernel");
    
    // Implementation outline
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    
    // Get thread ID
    auto globalId = getGlobalThreadId();
    
    // TODO: Implement FM-index search
    // 1. Extract seeds from read
    // 2. Search seeds in FM-index
    // 3. Extend alignments
    // 4. Score and filter
    
    builder->CreateRetVoid();
    
    return func;
}

llvm::Function* LLVMCodeGenerator::generateResistanceScanner() {
    std::vector<llvm::Type*> paramTypes = {
        llvm::PointerType::get(dnaSequenceType, 1),              // reads
        llvm::Type::getInt32Ty(*context),                        // num_reads
        llvm::PointerType::get(resistanceMutationType, 1),       // known_mutations
        llvm::Type::getInt32Ty(*context),                        // num_mutations
        llvm::PointerType::get(llvm::Type::getInt8Ty(*context), 1), // results
    };
    
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getVoidTy(*context), paramTypes, false);
    
    auto func = llvm::Function::Create(
        funcType, llvm::Function::ExternalLinkage,
        "scan_resistance_mutations_kernel", module.get());
    
    func->addFnAttr("nvvm.annotations", "kernel");
    
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    
    // Implementation for QRDR mutation detection
    auto globalId = getGlobalThreadId();
    
    // TODO: Implement mutation scanning
    // 1. Map read to gene
    // 2. Check specific codon positions
    // 3. Translate to amino acid
    // 4. Compare with known mutations
    
    builder->CreateRetVoid();
    
    return func;
}

llvm::Function* LLVMCodeGenerator::generateAbundanceCalculator() {
    std::vector<llvm::Type*> paramTypes = {
        llvm::PointerType::get(alignmentResultType, 1),     // alignments
        llvm::Type::getInt32Ty(*context),                   // num_alignments
        llvm::Type::getInt32Ty(*context),                   // num_species
        llvm::PointerType::get(llvm::Type::getFloatTy(*context), 1)  // abundances
    };
    
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getVoidTy(*context), paramTypes, false);
    
    auto func = llvm::Function::Create(
        funcType, llvm::Function::ExternalLinkage,
        "calculate_abundance_kernel", module.get());
    
    func->addFnAttr("nvvm.annotations", "kernel");
    
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    
    // Calculate TPM normalization
    // TODO: Implement proper abundance calculation
    
    builder->CreateRetVoid();
    
    return func;
}

llvm::Function* LLVMCodeGenerator::generateSmithWatermanGPU() {
    std::vector<llvm::Type*> paramTypes = {
        llvm::PointerType::get(llvm::Type::getInt8Ty(*context), 1),  // seq1
        llvm::Type::getInt32Ty(*context),    // len1
        llvm::PointerType::get(llvm::Type::getInt8Ty(*context), 1),  // seq2
        llvm::Type::getInt32Ty(*context),    // len2
        llvm::PointerType::get(llvm::Type::getInt32Ty(*context), 3), // scoring matrix (shared)
        llvm::Type::getInt32Ty(*context),    // gap_open
        llvm::Type::getInt32Ty(*context),    // gap_extend
        llvm::PointerType::get(llvm::Type::getInt32Ty(*context), 1)  // result
    };
    
    auto funcType = llvm::FunctionType::get(
        llvm::Type::getVoidTy(*context), paramTypes, false);
    
    auto func = llvm::Function::Create(
        funcType, llvm::Function::ExternalLinkage,
        "smith_waterman_gpu_kernel", module.get());
    
    func->addFnAttr("nvvm.annotations", "kernel");
    
    // Mark shared memory parameter
    func->getArg(4)->addAttr(llvm::Attribute::get(*context, "nvvm.shared"));
    
    auto entryBB = llvm::BasicBlock::Create(*context, "entry", func);
    builder->SetInsertPoint(entryBB);
    
    // Implement tiled dynamic programming
    // TODO: Full implementation with shared memory optimization
    
    builder->CreateRetVoid();
    
    return func;
}

void LLVMCodeGenerator::generateBioKernels() {
    generateKmerCounter();
    generateShortReadAligner();
    generateResistanceScanner();
    generateAbundanceCalculator();
    generateSmithWatermanGPU();
}

void LLVMCodeGenerator::print() {
    module->print(llvm::outs(), nullptr);
}

bool LLVMCodeGenerator::verify() {
    return !llvm::verifyModule(*module, &llvm::errs());
}

std::string LLVMCodeGenerator::generatePTX() {
    // This would use LLVM's NVPTX backend
    // For now, return placeholder
    return "; PTX assembly would be generated here\n";
}

} // namespace biogpu