
// Generated from BioGPU.g4 by ANTLR 4.10.1

#pragma once


#include "antlr4-runtime.h"
#include "BioGPUParser.h"



/**
 * This class defines an abstract visitor for a parse tree
 * produced by BioGPUParser.
 */
class  BioGPUVisitor : public antlr4::tree::AbstractParseTreeVisitor {
public:

  /**
   * Visit parse trees produced by BioGPUParser.
   */
    virtual std::any visitProgram(BioGPUParser::ProgramContext *context) = 0;

    virtual std::any visitImportStatement(BioGPUParser::ImportStatementContext *context) = 0;

    virtual std::any visitModulePath(BioGPUParser::ModulePathContext *context) = 0;

    virtual std::any visitPipeline(BioGPUParser::PipelineContext *context) = 0;

    virtual std::any visitPipelineBody(BioGPUParser::PipelineBodyContext *context) = 0;

    virtual std::any visitInputDecl(BioGPUParser::InputDeclContext *context) = 0;

    virtual std::any visitInputParam(BioGPUParser::InputParamContext *context) = 0;

    virtual std::any visitOutputDecl(BioGPUParser::OutputDeclContext *context) = 0;

    virtual std::any visitOutputParam(BioGPUParser::OutputParamContext *context) = 0;

    virtual std::any visitReferencesDecl(BioGPUParser::ReferencesDeclContext *context) = 0;

    virtual std::any visitReferenceParam(BioGPUParser::ReferenceParamContext *context) = 0;

    virtual std::any visitStage(BioGPUParser::StageContext *context) = 0;

    virtual std::any visitDecorator(BioGPUParser::DecoratorContext *context) = 0;

    virtual std::any visitDecoratorArgs(BioGPUParser::DecoratorArgsContext *context) = 0;

    virtual std::any visitDecoratorArg(BioGPUParser::DecoratorArgContext *context) = 0;

    virtual std::any visitStageBody(BioGPUParser::StageBodyContext *context) = 0;

    virtual std::any visitStatement(BioGPUParser::StatementContext *context) = 0;

    virtual std::any visitAssignment(BioGPUParser::AssignmentContext *context) = 0;

    virtual std::any visitFunctionCall(BioGPUParser::FunctionCallContext *context) = 0;

    virtual std::any visitArgumentList(BioGPUParser::ArgumentListContext *context) = 0;

    virtual std::any visitConfigBlock(BioGPUParser::ConfigBlockContext *context) = 0;

    virtual std::any visitConfigParam(BioGPUParser::ConfigParamContext *context) = 0;

    virtual std::any visitParallelMap(BioGPUParser::ParallelMapContext *context) = 0;

    virtual std::any visitParallelConfig(BioGPUParser::ParallelConfigContext *context) = 0;

    virtual std::any visitIfStatement(BioGPUParser::IfStatementContext *context) = 0;

    virtual std::any visitForStatement(BioGPUParser::ForStatementContext *context) = 0;

    virtual std::any visitEmitStatement(BioGPUParser::EmitStatementContext *context) = 0;

    virtual std::any visitReportBlock(BioGPUParser::ReportBlockContext *context) = 0;

    virtual std::any visitReportStatement(BioGPUParser::ReportStatementContext *context) = 0;

    virtual std::any visitExpression(BioGPUParser::ExpressionContext *context) = 0;

    virtual std::any visitArray(BioGPUParser::ArrayContext *context) = 0;

    virtual std::any visitLiteral(BioGPUParser::LiteralContext *context) = 0;

    virtual std::any visitDataType(BioGPUParser::DataTypeContext *context) = 0;

    virtual std::any visitReferenceType(BioGPUParser::ReferenceTypeContext *context) = 0;


};

