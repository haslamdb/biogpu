
// Generated from BioGPU.g4 by ANTLR 4.10.1

#pragma once


#include "antlr4-runtime.h"




class  BioGPULexer : public antlr4::Lexer {
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

  explicit BioGPULexer(antlr4::CharStream *input);

  ~BioGPULexer() override;


  std::string getGrammarFileName() const override;

  const std::vector<std::string>& getRuleNames() const override;

  const std::vector<std::string>& getChannelNames() const override;

  const std::vector<std::string>& getModeNames() const override;

  const antlr4::dfa::Vocabulary& getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;

  const antlr4::atn::ATN& getATN() const override;

  // By default the static state used to implement the lexer is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

private:

  // Individual action functions triggered by action() above.

  // Individual semantic predicate functions triggered by sempred() above.

};

