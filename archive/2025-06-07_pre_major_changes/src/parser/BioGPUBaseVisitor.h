
// Generated from BioGPU.g4 by ANTLR 4.10.1

#pragma once


#include "antlr4-runtime.h"
#include "BioGPUVisitor.h"


/**
 * This class provides an empty implementation of BioGPUVisitor, which can be
 * extended to create a visitor which only needs to handle a subset of the available methods.
 */
class  BioGPUBaseVisitor : public BioGPUVisitor {
public:

  virtual std::any visitProgram(BioGPUParser::ProgramContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitImportStatement(BioGPUParser::ImportStatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitModulePath(BioGPUParser::ModulePathContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitPipeline(BioGPUParser::PipelineContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitPipelineBody(BioGPUParser::PipelineBodyContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitInputDecl(BioGPUParser::InputDeclContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitInputParam(BioGPUParser::InputParamContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitOutputDecl(BioGPUParser::OutputDeclContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitOutputParam(BioGPUParser::OutputParamContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitReferencesDecl(BioGPUParser::ReferencesDeclContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitReferenceParam(BioGPUParser::ReferenceParamContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitStage(BioGPUParser::StageContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitDecorator(BioGPUParser::DecoratorContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitDecoratorArgs(BioGPUParser::DecoratorArgsContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitDecoratorArg(BioGPUParser::DecoratorArgContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitStageBody(BioGPUParser::StageBodyContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitStatement(BioGPUParser::StatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitAssignment(BioGPUParser::AssignmentContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitFunctionCall(BioGPUParser::FunctionCallContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitArgumentList(BioGPUParser::ArgumentListContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConfigBlock(BioGPUParser::ConfigBlockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitConfigParam(BioGPUParser::ConfigParamContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitParallelMap(BioGPUParser::ParallelMapContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitParallelConfig(BioGPUParser::ParallelConfigContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitIfStatement(BioGPUParser::IfStatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitForStatement(BioGPUParser::ForStatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitEmitStatement(BioGPUParser::EmitStatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitReportBlock(BioGPUParser::ReportBlockContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitReportStatement(BioGPUParser::ReportStatementContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitExpression(BioGPUParser::ExpressionContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitArray(BioGPUParser::ArrayContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitLiteral(BioGPUParser::LiteralContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitDataType(BioGPUParser::DataTypeContext *ctx) override {
    return visitChildren(ctx);
  }

  virtual std::any visitReferenceType(BioGPUParser::ReferenceTypeContext *ctx) override {
    return visitChildren(ctx);
  }


};

