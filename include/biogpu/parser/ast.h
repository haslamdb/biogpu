#ifndef BIOGPU_AST_H
#define BIOGPU_AST_H

#include <string>
#include <vector>
#include <memory>
#include <variant>
#include <map>

namespace biogpu {

// Forward declarations
class ASTVisitor;

// Base class for all AST nodes
class ASTNode {
public:
    virtual ~ASTNode() = default;
    virtual void accept(ASTVisitor* visitor) = 0;
    
    // Source location info
    int line = 0;
    int column = 0;
};

// Type definitions
enum class DataType {
    FASTQ_FILE,
    FASTA_FILE,
    CSV_FILE,
    JSON,
    PDF,
    STRING,
    INT,
    FLOAT,
    BOOLEAN,
    ARRAY,
    MAP,
    GENOME_DATABASE,
    MUTATION_DATABASE,
    GENE_DATABASE,
    UNKNOWN
};

// Literals
class Literal : public ASTNode {
public:
    virtual ~Literal() = default;
};

class NumberLiteral : public Literal {
public:
    double value;
    bool isInteger;
    
    explicit NumberLiteral(double val, bool isInt = false) 
        : value(val), isInteger(isInt) {}
    void accept(ASTVisitor* visitor) override;
};

class StringLiteral : public Literal {
public:
    std::string value;
    
    explicit StringLiteral(const std::string& val) : value(val) {}
    void accept(ASTVisitor* visitor) override;
};

class BooleanLiteral : public Literal {
public:
    bool value;
    
    explicit BooleanLiteral(bool val) : value(val) {}
    void accept(ASTVisitor* visitor) override;
};

class DNASequenceLiteral : public Literal {
public:
    std::string sequence;
    
    explicit DNASequenceLiteral(const std::string& seq) : sequence(seq) {}
    void accept(ASTVisitor* visitor) override;
};

// Expressions
class Expression : public ASTNode {
public:
    virtual ~Expression() = default;
};

class Identifier : public Expression {
public:
    std::string name;
    
    explicit Identifier(const std::string& n) : name(n) {}
    void accept(ASTVisitor* visitor) override;
};

class ArrayExpression : public Expression {
public:
    std::vector<std::unique_ptr<Expression>> elements;
    
    void accept(ASTVisitor* visitor) override;
};

class MemberAccess : public Expression {
public:
    std::unique_ptr<Expression> object;
    std::string member;
    
    MemberAccess(std::unique_ptr<Expression> obj, const std::string& mem)
        : object(std::move(obj)), member(mem) {}
    void accept(ASTVisitor* visitor) override;
};

class IndexAccess : public Expression {
public:
    std::unique_ptr<Expression> array;
    std::unique_ptr<Expression> index;
    
    IndexAccess(std::unique_ptr<Expression> arr, std::unique_ptr<Expression> idx)
        : array(std::move(arr)), index(std::move(idx)) {}
    void accept(ASTVisitor* visitor) override;
};

// Function/method calls
class FunctionCall : public Expression {
public:
    std::string functionName;
    std::vector<std::unique_ptr<Expression>> arguments;
    std::map<std::string, std::unique_ptr<Expression>> config;
    
    explicit FunctionCall(const std::string& name) : functionName(name) {}
    void accept(ASTVisitor* visitor) override;
};

class ParallelMap : public Expression {
public:
    std::unique_ptr<Expression> input;
    std::unique_ptr<Expression> database;
    std::map<std::string, std::unique_ptr<Expression>> config;
    
    void accept(ASTVisitor* visitor) override;
};

// Statements
class Statement : public ASTNode {
public:
    virtual ~Statement() = default;
};

class Assignment : public Statement {
public:
    std::string variable;
    std::unique_ptr<Expression> value;
    
    Assignment(const std::string& var, std::unique_ptr<Expression> val)
        : variable(var), value(std::move(val)) {}
    void accept(ASTVisitor* visitor) override;
};

class EmitStatement : public Statement {
public:
    std::vector<std::string> variables;
    
    void accept(ASTVisitor* visitor) override;
};

class PrintStatement : public Statement {
public:
    std::string message;
    std::vector<std::string> interpolations;  // Variables to interpolate
    
    explicit PrintStatement(const std::string& msg) : message(msg) {}
    void accept(ASTVisitor* visitor) override;
};

class IfStatement : public Statement {
public:
    std::unique_ptr<Expression> condition;
    std::vector<std::unique_ptr<Statement>> thenBody;
    std::vector<std::unique_ptr<Statement>> elseBody;
    
    void accept(ASTVisitor* visitor) override;
};

class ForStatement : public Statement {
public:
    std::string variable;
    std::unique_ptr<Expression> iterable;
    std::vector<std::unique_ptr<Statement>> body;
    
    void accept(ASTVisitor* visitor) override;
};

// Report block statements
class ReportBlock : public Statement {
public:
    std::vector<std::unique_ptr<Statement>> statements;
    
    void accept(ASTVisitor* visitor) override;
};

class AlertStatement : public Statement {
public:
    std::string message;
    
    explicit AlertStatement(const std::string& msg) : message(msg) {}
    void accept(ASTVisitor* visitor) override;
};

class RecommendationStatement : public Statement {
public:
    std::string message;
    
    explicit RecommendationStatement(const std::string& msg) : message(msg) {}
    void accept(ASTVisitor* visitor) override;
};

// Pipeline components
class InputParameter {
public:
    std::string name;
    DataType type;
    bool optional = false;
    
    InputParameter(const std::string& n, DataType t, bool opt = false)
        : name(n), type(t), optional(opt) {}
};

class OutputParameter {
public:
    std::string name;
    DataType type;
    
    OutputParameter(const std::string& n, DataType t) : name(n), type(t) {}
};

class ReferenceParameter {
public:
    std::string name;
    DataType type;
    std::string source;
    
    ReferenceParameter(const std::string& n, DataType t, const std::string& s)
        : name(n), type(t), source(s) {}
};

// Stage definition
class Stage : public ASTNode {
public:
    std::string name;
    std::vector<std::string> decorators;
    std::vector<std::unique_ptr<Statement>> statements;
    std::map<std::string, std::string> decoratorArgs;
    
    explicit Stage(const std::string& n) : name(n) {}
    void accept(ASTVisitor* visitor) override;
};

// Pipeline definition
class Pipeline : public ASTNode {
public:
    std::string name;
    std::vector<InputParameter> inputs;
    std::vector<OutputParameter> outputs;
    std::vector<ReferenceParameter> references;
    std::vector<std::unique_ptr<Stage>> stages;
    
    explicit Pipeline(const std::string& n) : name(n) {}
    void accept(ASTVisitor* visitor) override;
};

// Import statement
class ImportStatement : public ASTNode {
public:
    std::vector<std::string> modulePath;
    std::string alias;  // Optional
    
    void accept(ASTVisitor* visitor) override;
};

// Program root
class Program : public ASTNode {
public:
    std::vector<std::unique_ptr<ImportStatement>> imports;
    std::vector<std::unique_ptr<Pipeline>> pipelines;
    std::vector<std::unique_ptr<ASTNode>> declarations;  // Other top-level declarations
    
    void accept(ASTVisitor* visitor) override;
};

// Visitor pattern for AST traversal
class ASTVisitor {
public:
    virtual ~ASTVisitor() = default;
    
    // Visit methods for each node type
    virtual void visitProgram(Program* node) = 0;
    virtual void visitPipeline(Pipeline* node) = 0;
    virtual void visitStage(Stage* node) = 0;
    virtual void visitImportStatement(ImportStatement* node) = 0;
    
    // Literals
    virtual void visitNumberLiteral(NumberLiteral* node) = 0;
    virtual void visitStringLiteral(StringLiteral* node) = 0;
    virtual void visitBooleanLiteral(BooleanLiteral* node) = 0;
    virtual void visitDNASequenceLiteral(DNASequenceLiteral* node) = 0;
    
    // Expressions
    virtual void visitIdentifier(Identifier* node) = 0;
    virtual void visitArrayExpression(ArrayExpression* node) = 0;
    virtual void visitMemberAccess(MemberAccess* node) = 0;
    virtual void visitIndexAccess(IndexAccess* node) = 0;
    virtual void visitFunctionCall(FunctionCall* node) = 0;
    virtual void visitParallelMap(ParallelMap* node) = 0;
    
    // Statements
    virtual void visitAssignment(Assignment* node) = 0;
    virtual void visitEmitStatement(EmitStatement* node) = 0;
    virtual void visitPrintStatement(PrintStatement* node) = 0;
    virtual void visitAlertStatement(AlertStatement* node) = 0;
    virtual void visitRecommendationStatement(RecommendationStatement* node) = 0;
    virtual void visitIfStatement(IfStatement* node) = 0;
    virtual void visitForStatement(ForStatement* node) = 0;
    virtual void visitReportBlock(ReportBlock* node) = 0;
};

// Utility functions
DataType stringToDataType(const std::string& type);
std::string dataTypeToString(DataType type);

// AST printer for debugging
class ASTPrinter : public ASTVisitor {
public:
    void print(ASTNode* node);
    
    // Implement all visit methods
    void visitProgram(Program* node) override;
    void visitPipeline(Pipeline* node) override;
    void visitStage(Stage* node) override;
    void visitImportStatement(ImportStatement* node) override;
    void visitNumberLiteral(NumberLiteral* node) override;
    void visitStringLiteral(StringLiteral* node) override;
    void visitBooleanLiteral(BooleanLiteral* node) override;
    void visitDNASequenceLiteral(DNASequenceLiteral* node) override;
    void visitIdentifier(Identifier* node) override;
    void visitArrayExpression(ArrayExpression* node) override;
    void visitMemberAccess(MemberAccess* node) override;
    void visitIndexAccess(IndexAccess* node) override;
    void visitFunctionCall(FunctionCall* node) override;
    void visitParallelMap(ParallelMap* node) override;
    void visitAssignment(Assignment* node) override;
    void visitEmitStatement(EmitStatement* node) override;
    void visitPrintStatement(PrintStatement* node) override;
    void visitAlertStatement(AlertStatement* node) override;
    void visitRecommendationStatement(RecommendationStatement* node) override;
    void visitIfStatement(IfStatement* node) override;
    void visitForStatement(ForStatement* node) override;
    void visitReportBlock(ReportBlock* node) override;
    
private:
    int indentLevel = 0;
    void indent();
};

} // namespace biogpu

#endif // BIOGPU_AST_H