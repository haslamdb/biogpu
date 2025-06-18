
// Generated from BioGPU.g4 by ANTLR 4.10.1

#pragma once


#include "antlr4-runtime.h"




class  BioGPUParser : public antlr4::Parser {
public:
  enum {
    T__0 = 1, PIPELINE = 2, STAGE = 3, INPUT = 4, OUTPUT = 5, REFERENCES = 6, 
    IMPORT = 7, AS = 8, IF = 9, ELSE = 10, FOR = 11, FOREACH = 12, IN = 13, 
    EMIT = 14, REPORT = 15, PRINT = 16, ALERT = 17, RECOMMENDATION = 18, 
    PARALLEL_MAP = 19, OPTIONAL = 20, FASTQ_FILE = 21, FASTA_FILE = 22, 
    CSV_FILE = 23, JSON = 24, PDF = 25, STRING_TYPE = 26, INT_TYPE = 27, 
    FLOAT_TYPE = 28, BOOLEAN_TYPE = 29, ARRAY = 30, MAP = 31, GENOME_DATABASE = 32, 
    MUTATION_DATABASE = 33, GENE_DATABASE = 34, GPU_KERNEL = 35, PARALLEL = 36, 
    ASSIGN = 37, DOT = 38, COMMA = 39, SEMI = 40, COLON = 41, LPAREN = 42, 
    RPAREN = 43, LBRACE = 44, RBRACE = 45, LBRACK = 46, RBRACK = 47, LT = 48, 
    GT = 49, BOOLEAN = 50, NUMBER = 51, STRING = 52, DNA_SEQUENCE = 53, 
    ID = 54, COMMENT = 55, BLOCK_COMMENT = 56, WS = 57
  };

  enum {
    RuleProgram = 0, RuleImportStatement = 1, RuleModulePath = 2, RulePipeline = 3, 
    RulePipelineBody = 4, RuleInputDecl = 5, RuleInputParam = 6, RuleOutputDecl = 7, 
    RuleOutputParam = 8, RuleReferencesDecl = 9, RuleReferenceParam = 10, 
    RuleStage = 11, RuleDecorator = 12, RuleDecoratorArgs = 13, RuleDecoratorArg = 14, 
    RuleStageBody = 15, RuleStatement = 16, RuleAssignment = 17, RuleFunctionCall = 18, 
    RuleArgumentList = 19, RuleConfigBlock = 20, RuleConfigParam = 21, RuleParallelMap = 22, 
    RuleParallelConfig = 23, RuleIfStatement = 24, RuleForStatement = 25, 
    RuleEmitStatement = 26, RuleReportBlock = 27, RuleReportStatement = 28, 
    RuleExpression = 29, RuleArray = 30, RuleLiteral = 31, RuleDataType = 32, 
    RuleReferenceType = 33
  };

  explicit BioGPUParser(antlr4::TokenStream *input);

  BioGPUParser(antlr4::TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options);

  ~BioGPUParser() override;

  std::string getGrammarFileName() const override;

  const antlr4::atn::ATN& getATN() const override;

  const std::vector<std::string>& getRuleNames() const override;

  const antlr4::dfa::Vocabulary& getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;


  class ProgramContext;
  class ImportStatementContext;
  class ModulePathContext;
  class PipelineContext;
  class PipelineBodyContext;
  class InputDeclContext;
  class InputParamContext;
  class OutputDeclContext;
  class OutputParamContext;
  class ReferencesDeclContext;
  class ReferenceParamContext;
  class StageContext;
  class DecoratorContext;
  class DecoratorArgsContext;
  class DecoratorArgContext;
  class StageBodyContext;
  class StatementContext;
  class AssignmentContext;
  class FunctionCallContext;
  class ArgumentListContext;
  class ConfigBlockContext;
  class ConfigParamContext;
  class ParallelMapContext;
  class ParallelConfigContext;
  class IfStatementContext;
  class ForStatementContext;
  class EmitStatementContext;
  class ReportBlockContext;
  class ReportStatementContext;
  class ExpressionContext;
  class ArrayContext;
  class LiteralContext;
  class DataTypeContext;
  class ReferenceTypeContext; 

  class  ProgramContext : public antlr4::ParserRuleContext {
  public:
    ProgramContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *EOF();
    std::vector<ImportStatementContext *> importStatement();
    ImportStatementContext* importStatement(size_t i);
    std::vector<PipelineContext *> pipeline();
    PipelineContext* pipeline(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ProgramContext* program();

  class  ImportStatementContext : public antlr4::ParserRuleContext {
  public:
    ImportStatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *IMPORT();
    ModulePathContext *modulePath();
    antlr4::tree::TerminalNode *SEMI();
    antlr4::tree::TerminalNode *AS();
    antlr4::tree::TerminalNode *ID();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ImportStatementContext* importStatement();

  class  ModulePathContext : public antlr4::ParserRuleContext {
  public:
    ModulePathContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<antlr4::tree::TerminalNode *> ID();
    antlr4::tree::TerminalNode* ID(size_t i);
    std::vector<antlr4::tree::TerminalNode *> DOT();
    antlr4::tree::TerminalNode* DOT(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ModulePathContext* modulePath();

  class  PipelineContext : public antlr4::ParserRuleContext {
  public:
    PipelineContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *PIPELINE();
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *LBRACE();
    PipelineBodyContext *pipelineBody();
    antlr4::tree::TerminalNode *RBRACE();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  PipelineContext* pipeline();

  class  PipelineBodyContext : public antlr4::ParserRuleContext {
  public:
    PipelineBodyContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    InputDeclContext *inputDecl();
    OutputDeclContext *outputDecl();
    ReferencesDeclContext *referencesDecl();
    std::vector<StageContext *> stage();
    StageContext* stage(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  PipelineBodyContext* pipelineBody();

  class  InputDeclContext : public antlr4::ParserRuleContext {
  public:
    InputDeclContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *INPUT();
    antlr4::tree::TerminalNode *COLON();
    std::vector<InputParamContext *> inputParam();
    InputParamContext* inputParam(size_t i);
    antlr4::tree::TerminalNode *LBRACE();
    antlr4::tree::TerminalNode *RBRACE();
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  InputDeclContext* inputDecl();

  class  InputParamContext : public antlr4::ParserRuleContext {
  public:
    InputParamContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *COLON();
    DataTypeContext *dataType();
    antlr4::tree::TerminalNode *OPTIONAL();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  InputParamContext* inputParam();

  class  OutputDeclContext : public antlr4::ParserRuleContext {
  public:
    OutputDeclContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *OUTPUT();
    antlr4::tree::TerminalNode *COLON();
    std::vector<OutputParamContext *> outputParam();
    OutputParamContext* outputParam(size_t i);
    antlr4::tree::TerminalNode *LBRACE();
    antlr4::tree::TerminalNode *RBRACE();
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  OutputDeclContext* outputDecl();

  class  OutputParamContext : public antlr4::ParserRuleContext {
  public:
    OutputParamContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *COLON();
    DataTypeContext *dataType();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  OutputParamContext* outputParam();

  class  ReferencesDeclContext : public antlr4::ParserRuleContext {
  public:
    ReferencesDeclContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *REFERENCES();
    antlr4::tree::TerminalNode *COLON();
    antlr4::tree::TerminalNode *LBRACE();
    std::vector<ReferenceParamContext *> referenceParam();
    ReferenceParamContext* referenceParam(size_t i);
    antlr4::tree::TerminalNode *RBRACE();
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ReferencesDeclContext* referencesDecl();

  class  ReferenceParamContext : public antlr4::ParserRuleContext {
  public:
    ReferenceParamContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *COLON();
    ReferenceTypeContext *referenceType();
    antlr4::tree::TerminalNode *LPAREN();
    antlr4::tree::TerminalNode *STRING();
    antlr4::tree::TerminalNode *RPAREN();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ReferenceParamContext* referenceParam();

  class  StageContext : public antlr4::ParserRuleContext {
  public:
    StageContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *STAGE();
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *LBRACE();
    StageBodyContext *stageBody();
    antlr4::tree::TerminalNode *RBRACE();
    std::vector<DecoratorContext *> decorator();
    DecoratorContext* decorator(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  StageContext* stage();

  class  DecoratorContext : public antlr4::ParserRuleContext {
  public:
    DecoratorContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *LPAREN();
    DecoratorArgsContext *decoratorArgs();
    antlr4::tree::TerminalNode *RPAREN();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DecoratorContext* decorator();

  class  DecoratorArgsContext : public antlr4::ParserRuleContext {
  public:
    DecoratorArgsContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<DecoratorArgContext *> decoratorArg();
    DecoratorArgContext* decoratorArg(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DecoratorArgsContext* decoratorArgs();

  class  DecoratorArgContext : public antlr4::ParserRuleContext {
  public:
    DecoratorArgContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<antlr4::tree::TerminalNode *> ID();
    antlr4::tree::TerminalNode* ID(size_t i);
    antlr4::tree::TerminalNode *COLON();
    LiteralContext *literal();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DecoratorArgContext* decoratorArg();

  class  StageBodyContext : public antlr4::ParserRuleContext {
  public:
    StageBodyContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<StatementContext *> statement();
    StatementContext* statement(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  StageBodyContext* stageBody();

  class  StatementContext : public antlr4::ParserRuleContext {
  public:
    StatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    AssignmentContext *assignment();
    FunctionCallContext *functionCall();
    ParallelMapContext *parallelMap();
    IfStatementContext *ifStatement();
    ForStatementContext *forStatement();
    EmitStatementContext *emitStatement();
    ReportBlockContext *reportBlock();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  StatementContext* statement();

  class  AssignmentContext : public antlr4::ParserRuleContext {
  public:
    AssignmentContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *ASSIGN();
    ExpressionContext *expression();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  AssignmentContext* assignment();

  class  FunctionCallContext : public antlr4::ParserRuleContext {
  public:
    FunctionCallContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *LPAREN();
    antlr4::tree::TerminalNode *RPAREN();
    ArgumentListContext *argumentList();
    antlr4::tree::TerminalNode *LBRACE();
    ConfigBlockContext *configBlock();
    antlr4::tree::TerminalNode *RBRACE();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  FunctionCallContext* functionCall();

  class  ArgumentListContext : public antlr4::ParserRuleContext {
  public:
    ArgumentListContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<ExpressionContext *> expression();
    ExpressionContext* expression(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ArgumentListContext* argumentList();

  class  ConfigBlockContext : public antlr4::ParserRuleContext {
  public:
    ConfigBlockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<ConfigParamContext *> configParam();
    ConfigParamContext* configParam(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ConfigBlockContext* configBlock();

  class  ConfigParamContext : public antlr4::ParserRuleContext {
  public:
    ConfigParamContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *COLON();
    LiteralContext *literal();
    ExpressionContext *expression();
    ArrayContext *array();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ConfigParamContext* configParam();

  class  ParallelMapContext : public antlr4::ParserRuleContext {
  public:
    ParallelMapContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *PARALLEL_MAP();
    antlr4::tree::TerminalNode *LPAREN();
    std::vector<ExpressionContext *> expression();
    ExpressionContext* expression(size_t i);
    antlr4::tree::TerminalNode *COMMA();
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *LBRACE();
    ParallelConfigContext *parallelConfig();
    antlr4::tree::TerminalNode *RBRACE();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ParallelMapContext* parallelMap();

  class  ParallelConfigContext : public antlr4::ParserRuleContext {
  public:
    ParallelConfigContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    std::vector<ConfigParamContext *> configParam();
    ConfigParamContext* configParam(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ParallelConfigContext* parallelConfig();

  class  IfStatementContext : public antlr4::ParserRuleContext {
  public:
    IfStatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *IF();
    ExpressionContext *expression();
    std::vector<antlr4::tree::TerminalNode *> LBRACE();
    antlr4::tree::TerminalNode* LBRACE(size_t i);
    std::vector<StageBodyContext *> stageBody();
    StageBodyContext* stageBody(size_t i);
    std::vector<antlr4::tree::TerminalNode *> RBRACE();
    antlr4::tree::TerminalNode* RBRACE(size_t i);
    antlr4::tree::TerminalNode *ELSE();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  IfStatementContext* ifStatement();

  class  ForStatementContext : public antlr4::ParserRuleContext {
  public:
    ForStatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *FOR();
    antlr4::tree::TerminalNode *ID();
    antlr4::tree::TerminalNode *IN();
    ExpressionContext *expression();
    antlr4::tree::TerminalNode *LBRACE();
    StageBodyContext *stageBody();
    antlr4::tree::TerminalNode *RBRACE();
    antlr4::tree::TerminalNode *FOREACH();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ForStatementContext* forStatement();

  class  EmitStatementContext : public antlr4::ParserRuleContext {
  public:
    EmitStatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *EMIT();
    antlr4::tree::TerminalNode *COLON();
    std::vector<antlr4::tree::TerminalNode *> ID();
    antlr4::tree::TerminalNode* ID(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  EmitStatementContext* emitStatement();

  class  ReportBlockContext : public antlr4::ParserRuleContext {
  public:
    ReportBlockContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *REPORT();
    antlr4::tree::TerminalNode *LBRACE();
    antlr4::tree::TerminalNode *RBRACE();
    std::vector<ReportStatementContext *> reportStatement();
    ReportStatementContext* reportStatement(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ReportBlockContext* reportBlock();

  class  ReportStatementContext : public antlr4::ParserRuleContext {
  public:
    ReportStatementContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *PRINT();
    antlr4::tree::TerminalNode *STRING();
    antlr4::tree::TerminalNode *ALERT();
    antlr4::tree::TerminalNode *COLON();
    antlr4::tree::TerminalNode *RECOMMENDATION();
    IfStatementContext *ifStatement();
    ForStatementContext *forStatement();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ReportStatementContext* reportStatement();

  class  ExpressionContext : public antlr4::ParserRuleContext {
  public:
    ExpressionContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *ID();
    LiteralContext *literal();
    FunctionCallContext *functionCall();
    antlr4::tree::TerminalNode *LPAREN();
    std::vector<ExpressionContext *> expression();
    ExpressionContext* expression(size_t i);
    antlr4::tree::TerminalNode *RPAREN();
    antlr4::tree::TerminalNode *DOT();
    antlr4::tree::TerminalNode *LBRACK();
    antlr4::tree::TerminalNode *RBRACK();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ExpressionContext* expression();
  ExpressionContext* expression(int precedence);
  class  ArrayContext : public antlr4::ParserRuleContext {
  public:
    ArrayContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *LBRACK();
    antlr4::tree::TerminalNode *RBRACK();
    std::vector<LiteralContext *> literal();
    LiteralContext* literal(size_t i);
    std::vector<antlr4::tree::TerminalNode *> ID();
    antlr4::tree::TerminalNode* ID(size_t i);
    std::vector<antlr4::tree::TerminalNode *> COMMA();
    antlr4::tree::TerminalNode* COMMA(size_t i);


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ArrayContext* array();

  class  LiteralContext : public antlr4::ParserRuleContext {
  public:
    LiteralContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *NUMBER();
    antlr4::tree::TerminalNode *STRING();
    antlr4::tree::TerminalNode *BOOLEAN();
    antlr4::tree::TerminalNode *DNA_SEQUENCE();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  LiteralContext* literal();

  class  DataTypeContext : public antlr4::ParserRuleContext {
  public:
    DataTypeContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *FASTQ_FILE();
    antlr4::tree::TerminalNode *FASTA_FILE();
    antlr4::tree::TerminalNode *CSV_FILE();
    antlr4::tree::TerminalNode *JSON();
    antlr4::tree::TerminalNode *PDF();
    antlr4::tree::TerminalNode *STRING_TYPE();
    antlr4::tree::TerminalNode *INT_TYPE();
    antlr4::tree::TerminalNode *FLOAT_TYPE();
    antlr4::tree::TerminalNode *BOOLEAN_TYPE();
    antlr4::tree::TerminalNode *ARRAY();
    antlr4::tree::TerminalNode *LT();
    std::vector<DataTypeContext *> dataType();
    DataTypeContext* dataType(size_t i);
    antlr4::tree::TerminalNode *GT();
    antlr4::tree::TerminalNode *MAP();
    antlr4::tree::TerminalNode *COMMA();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  DataTypeContext* dataType();

  class  ReferenceTypeContext : public antlr4::ParserRuleContext {
  public:
    ReferenceTypeContext(antlr4::ParserRuleContext *parent, size_t invokingState);
    virtual size_t getRuleIndex() const override;
    antlr4::tree::TerminalNode *GENOME_DATABASE();
    antlr4::tree::TerminalNode *MUTATION_DATABASE();
    antlr4::tree::TerminalNode *GENE_DATABASE();


    virtual std::any accept(antlr4::tree::ParseTreeVisitor *visitor) override;
   
  };

  ReferenceTypeContext* referenceType();


  bool sempred(antlr4::RuleContext *_localctx, size_t ruleIndex, size_t predicateIndex) override;

  bool expressionSempred(ExpressionContext *_localctx, size_t predicateIndex);

  // By default the static state used to implement the parser is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

private:
};

