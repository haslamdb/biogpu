// src/parser/ast.cpp
#include "biogpu/parser/ast.h"
#include <iostream>
#include <iomanip>

namespace biogpu {

// Accept methods for visitor pattern
void NumberLiteral::accept(ASTVisitor* visitor) { visitor->visitNumberLiteral(this); }
void StringLiteral::accept(ASTVisitor* visitor) { visitor->visitStringLiteral(this); }
void BooleanLiteral::accept(ASTVisitor* visitor) { visitor->visitBooleanLiteral(this); }
void DNASequenceLiteral::accept(ASTVisitor* visitor) { visitor->visitDNASequenceLiteral(this); }
void Identifier::accept(ASTVisitor* visitor) { visitor->visitIdentifier(this); }
void ArrayExpression::accept(ASTVisitor* visitor) { visitor->visitArrayExpression(this); }
void MemberAccess::accept(ASTVisitor* visitor) { visitor->visitMemberAccess(this); }
void IndexAccess::accept(ASTVisitor* visitor) { visitor->visitIndexAccess(this); }
void FunctionCall::accept(ASTVisitor* visitor) { visitor->visitFunctionCall(this); }
void ParallelMap::accept(ASTVisitor* visitor) { visitor->visitParallelMap(this); }
void Assignment::accept(ASTVisitor* visitor) { visitor->visitAssignment(this); }
void EmitStatement::accept(ASTVisitor* visitor) { visitor->visitEmitStatement(this); }
void PrintStatement::accept(ASTVisitor* visitor) { visitor->visitPrintStatement(this); }
void AlertStatement::accept(ASTVisitor* visitor) { visitor->visitAlertStatement(this); }
void RecommendationStatement::accept(ASTVisitor* visitor) { visitor->visitRecommendationStatement(this); }
void IfStatement::accept(ASTVisitor* visitor) { visitor->visitIfStatement(this); }
void ForStatement::accept(ASTVisitor* visitor) { visitor->visitForStatement(this); }
void ReportBlock::accept(ASTVisitor* visitor) { visitor->visitReportBlock(this); }
void Stage::accept(ASTVisitor* visitor) { visitor->visitStage(this); }
void Pipeline::accept(ASTVisitor* visitor) { visitor->visitPipeline(this); }
void ImportStatement::accept(ASTVisitor* visitor) { visitor->visitImportStatement(this); }
void Program::accept(ASTVisitor* visitor) { visitor->visitProgram(this); }

// Utility functions
DataType stringToDataType(const std::string& type) {
    if (type == "fastq_file") return DataType::FASTQ_FILE;
    if (type == "fasta_file") return DataType::FASTA_FILE;
    if (type == "csv_file") return DataType::CSV_FILE;
    if (type == "json") return DataType::JSON;
    if (type == "pdf") return DataType::PDF;
    if (type == "string") return DataType::STRING;
    if (type == "int") return DataType::INT;
    if (type == "float") return DataType::FLOAT;
    if (type == "boolean") return DataType::BOOLEAN;
    if (type == "array") return DataType::ARRAY;
    if (type == "map") return DataType::MAP;
    if (type == "genome_database") return DataType::GENOME_DATABASE;
    if (type == "mutation_database") return DataType::MUTATION_DATABASE;
    if (type == "gene_database") return DataType::GENE_DATABASE;
    return DataType::UNKNOWN;
}

std::string dataTypeToString(DataType type) {
    switch (type) {
        case DataType::FASTQ_FILE: return "fastq_file";
        case DataType::FASTA_FILE: return "fasta_file";
        case DataType::CSV_FILE: return "csv_file";
        case DataType::JSON: return "json";
        case DataType::PDF: return "pdf";
        case DataType::STRING: return "string";
        case DataType::INT: return "int";
        case DataType::FLOAT: return "float";
        case DataType::BOOLEAN: return "boolean";
        case DataType::ARRAY: return "array";
        case DataType::MAP: return "map";
        case DataType::GENOME_DATABASE: return "genome_database";
        case DataType::MUTATION_DATABASE: return "mutation_database";
        case DataType::GENE_DATABASE: return "gene_database";
        default: return "unknown";
    }
}

// AST Printer implementation
void ASTPrinter::print(ASTNode* node) {
    indentLevel = 0;
    node->accept(this);
}

void ASTPrinter::indent() {
    for (int i = 0; i < indentLevel; ++i) {
        std::cout << "  ";
    }
}

void ASTPrinter::visitProgram(Program* node) {
    std::cout << "Program {" << std::endl;
    indentLevel++;
    
    if (!node->imports.empty()) {
        indent();
        std::cout << "Imports:" << std::endl;
        indentLevel++;
        for (auto& import : node->imports) {
            import->accept(this);
        }
        indentLevel--;
    }
    
    if (!node->pipelines.empty()) {
        indent();
        std::cout << "Pipelines:" << std::endl;
        indentLevel++;
        for (auto& pipeline : node->pipelines) {
            pipeline->accept(this);
        }
        indentLevel--;
    }
    
    indentLevel--;
    std::cout << "}" << std::endl;
}

void ASTPrinter::visitPipeline(Pipeline* node) {
    indent();
    std::cout << "Pipeline " << node->name << " {" << std::endl;
    indentLevel++;
    
    // Print inputs
    if (!node->inputs.empty()) {
        indent();
        std::cout << "Inputs:" << std::endl;
        indentLevel++;
        for (const auto& input : node->inputs) {
            indent();
            std::cout << input.name << ": " << dataTypeToString(input.type);
            if (input.optional) std::cout << " (optional)";
            std::cout << std::endl;
        }
        indentLevel--;
    }
    
    // Print outputs
    if (!node->outputs.empty()) {
        indent();
        std::cout << "Outputs:" << std::endl;
        indentLevel++;
        for (const auto& output : node->outputs) {
            indent();
            std::cout << output.name << ": " << dataTypeToString(output.type) << std::endl;
        }
        indentLevel--;
    }
    
    // Print stages
    if (!node->stages.empty()) {
        indent();
        std::cout << "Stages:" << std::endl;
        indentLevel++;
        for (auto& stage : node->stages) {
            stage->accept(this);
        }
        indentLevel--;
    }
    
    indentLevel--;
    indent();
    std::cout << "}" << std::endl;
}

void ASTPrinter::visitStage(Stage* node) {
    indent();
    
    // Print decorators
    for (const auto& decorator : node->decorators) {
        std::cout << decorator << " ";
    }
    
    std::cout << "Stage " << node->name << " {" << std::endl;
    indentLevel++;
    
    for (auto& stmt : node->statements) {
        stmt->accept(this);
    }
    
    indentLevel--;
    indent();
    std::cout << "}" << std::endl;
}

void ASTPrinter::visitImportStatement(ImportStatement* node) {
    indent();
    std::cout << "import ";
    for (size_t i = 0; i < node->modulePath.size(); ++i) {
        if (i > 0) std::cout << ".";
        std::cout << node->modulePath[i];
    }
    if (!node->alias.empty()) {
        std::cout << " as " << node->alias;
    }
    std::cout << std::endl;
}

void ASTPrinter::visitNumberLiteral(NumberLiteral* node) {
    std::cout << node->value;
}

void ASTPrinter::visitStringLiteral(StringLiteral* node) {
    std::cout << "\"" << node->value << "\"";
}

void ASTPrinter::visitBooleanLiteral(BooleanLiteral* node) {
    std::cout << (node->value ? "true" : "false");
}

void ASTPrinter::visitDNASequenceLiteral(DNASequenceLiteral* node) {
    std::cout << "DNA:" << node->sequence;
}

void ASTPrinter::visitIdentifier(Identifier* node) {
    std::cout << node->name;
}

void ASTPrinter::visitArrayExpression(ArrayExpression* node) {
    std::cout << "[";
    for (size_t i = 0; i < node->elements.size(); ++i) {
        if (i > 0) std::cout << ", ";
        node->elements[i]->accept(this);
    }
    std::cout << "]";
}

void ASTPrinter::visitMemberAccess(MemberAccess* node) {
    node->object->accept(this);
    std::cout << "." << node->member;
}

void ASTPrinter::visitIndexAccess(IndexAccess* node) {
    node->array->accept(this);
    std::cout << "[";
    node->index->accept(this);
    std::cout << "]";
}

void ASTPrinter::visitFunctionCall(FunctionCall* node) {
    std::cout << node->functionName << "(";
    for (size_t i = 0; i < node->arguments.size(); ++i) {
        if (i > 0) std::cout << ", ";
        node->arguments[i]->accept(this);
    }
    std::cout << ")";
    
    if (!node->config.empty()) {
        std::cout << " {" << std::endl;
        indentLevel++;
        for (const auto& [key, value] : node->config) {
            indent();
            std::cout << key << ": ";
            value->accept(this);
            std::cout << std::endl;
        }
        indentLevel--;
        indent();
        std::cout << "}";
    }
}

void ASTPrinter::visitParallelMap(ParallelMap* node) {
    std::cout << "parallel_map(";
    node->input->accept(this);
    std::cout << ", ";
    node->database->accept(this);
    std::cout << ")";
}

void ASTPrinter::visitAssignment(Assignment* node) {
    indent();
    std::cout << node->variable << " = ";
    node->value->accept(this);
    std::cout << std::endl;
}

void ASTPrinter::visitEmitStatement(EmitStatement* node) {
    indent();
    std::cout << "emit: ";
    for (size_t i = 0; i < node->variables.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << node->variables[i];
    }
    std::cout << std::endl;
}

void ASTPrinter::visitPrintStatement(PrintStatement* node) {
    indent();
    std::cout << "print \"" << node->message << "\"" << std::endl;
}

void ASTPrinter::visitAlertStatement(AlertStatement* node) {
    indent();
    std::cout << "alert: \"" << node->message << "\"" << std::endl;
}

void ASTPrinter::visitRecommendationStatement(RecommendationStatement* node) {
    indent();
    std::cout << "recommendation: \"" << node->message << "\"" << std::endl;
}

void ASTPrinter::visitIfStatement(IfStatement* node) {
    indent();
    std::cout << "if ";
    node->condition->accept(this);
    std::cout << " {" << std::endl;
    
    indentLevel++;
    for (auto& stmt : node->thenBody) {
        stmt->accept(this);
    }
    indentLevel--;
    
    if (!node->elseBody.empty()) {
        indent();
        std::cout << "} else {" << std::endl;
        indentLevel++;
        for (auto& stmt : node->elseBody) {
            stmt->accept(this);
        }
        indentLevel--;
    }
    
    indent();
    std::cout << "}" << std::endl;
}

void ASTPrinter::visitForStatement(ForStatement* node) {
    indent();
    std::cout << "for " << node->variable << " in ";
    node->iterable->accept(this);
    std::cout << " {" << std::endl;
    
    indentLevel++;
    for (auto& stmt : node->body) {
        stmt->accept(this);
    }
    indentLevel--;
    
    indent();
    std::cout << "}" << std::endl;
}

void ASTPrinter::visitReportBlock(ReportBlock* node) {
    indent();
    std::cout << "report {" << std::endl;
    
    indentLevel++;
    for (auto& stmt : node->statements) {
        stmt->accept(this);
    }
    indentLevel--;
    
    indent();
    std::cout << "}" << std::endl;
}

} // namespace biogpu