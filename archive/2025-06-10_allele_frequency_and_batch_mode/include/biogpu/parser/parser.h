// include/biogpu/parser/parser.h
#ifndef BIOGPU_PARSER_H
#define BIOGPU_PARSER_H

#include <string>
#include <memory>
#include <vector>

namespace biogpu {

// Forward declarations
class Pipeline;
class CompilerError;

class BioGPUParser {
public:
    BioGPUParser();
    ~BioGPUParser();
    
    /**
     * Parse a BioGPU source file
     * @param filename Path to .biogpu file
     * @return Parsed pipeline AST
     */
    std::unique_ptr<Pipeline> parseFile(const std::string& filename);
    
    /**
     * Parse BioGPU source code
     * @param source Source code string
     * @return Parsed pipeline AST
     */
    std::unique_ptr<Pipeline> parseString(const std::string& source);
    
    /**
     * Get parse errors
     */
    const std::vector<CompilerError>& getErrors() const { return errors; }
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    std::vector<CompilerError> errors;
};

// Forward declarations - actual definitions in ast.h
class ASTNode;
class Pipeline;

class CompilerError {
public:
    std::string message;
    int line;
    int column;
    
    CompilerError(const std::string& msg, int l, int c) 
        : message(msg), line(l), column(c) {}
};

} // namespace biogpu

#endif // BIOGPU_PARSER_H