
// Generated from BioGPU.g4 by ANTLR 4.10.1


#include "BioGPUVisitor.h"

#include "BioGPUParser.h"


using namespace antlrcpp;

using namespace antlr4;

namespace {

struct BioGPUParserStaticData final {
  BioGPUParserStaticData(std::vector<std::string> ruleNames,
                        std::vector<std::string> literalNames,
                        std::vector<std::string> symbolicNames)
      : ruleNames(std::move(ruleNames)), literalNames(std::move(literalNames)),
        symbolicNames(std::move(symbolicNames)),
        vocabulary(this->literalNames, this->symbolicNames) {}

  BioGPUParserStaticData(const BioGPUParserStaticData&) = delete;
  BioGPUParserStaticData(BioGPUParserStaticData&&) = delete;
  BioGPUParserStaticData& operator=(const BioGPUParserStaticData&) = delete;
  BioGPUParserStaticData& operator=(BioGPUParserStaticData&&) = delete;

  std::vector<antlr4::dfa::DFA> decisionToDFA;
  antlr4::atn::PredictionContextCache sharedContextCache;
  const std::vector<std::string> ruleNames;
  const std::vector<std::string> literalNames;
  const std::vector<std::string> symbolicNames;
  const antlr4::dfa::Vocabulary vocabulary;
  antlr4::atn::SerializedATNView serializedATN;
  std::unique_ptr<antlr4::atn::ATN> atn;
};

std::once_flag biogpuParserOnceFlag;
BioGPUParserStaticData *biogpuParserStaticData = nullptr;

void biogpuParserInitialize() {
  assert(biogpuParserStaticData == nullptr);
  auto staticData = std::make_unique<BioGPUParserStaticData>(
    std::vector<std::string>{
      "program", "importStatement", "modulePath", "pipeline", "pipelineBody", 
      "inputDecl", "inputParam", "outputDecl", "outputParam", "referencesDecl", 
      "referenceParam", "stage", "decorator", "decoratorArgs", "decoratorArg", 
      "stageBody", "statement", "assignment", "functionCall", "argumentList", 
      "configBlock", "configParam", "parallelMap", "parallelConfig", "ifStatement", 
      "forStatement", "emitStatement", "reportBlock", "reportStatement", 
      "expression", "array", "literal", "dataType", "referenceType"
    },
    std::vector<std::string>{
      "", "'@'", "'pipeline'", "'stage'", "'input'", "'output'", "'references'", 
      "'import'", "'as'", "'if'", "'else'", "'for'", "'foreach'", "'in'", 
      "'emit'", "'report'", "'print'", "'alert'", "'recommendation'", "'parallel_map'", 
      "'optional'", "'fastq_file'", "'fasta_file'", "'csv_file'", "'json'", 
      "'pdf'", "'string'", "'int'", "'float'", "'boolean'", "'array'", "'map'", 
      "'genome_database'", "'mutation_database'", "'gene_database'", "'@gpu_kernel'", 
      "'@parallel'", "'='", "'.'", "','", "';'", "':'", "'('", "')'", "'{'", 
      "'}'", "'['", "']'", "'<'", "'>'"
    },
    std::vector<std::string>{
      "", "", "PIPELINE", "STAGE", "INPUT", "OUTPUT", "REFERENCES", "IMPORT", 
      "AS", "IF", "ELSE", "FOR", "FOREACH", "IN", "EMIT", "REPORT", "PRINT", 
      "ALERT", "RECOMMENDATION", "PARALLEL_MAP", "OPTIONAL", "FASTQ_FILE", 
      "FASTA_FILE", "CSV_FILE", "JSON", "PDF", "STRING_TYPE", "INT_TYPE", 
      "FLOAT_TYPE", "BOOLEAN_TYPE", "ARRAY", "MAP", "GENOME_DATABASE", "MUTATION_DATABASE", 
      "GENE_DATABASE", "GPU_KERNEL", "PARALLEL", "ASSIGN", "DOT", "COMMA", 
      "SEMI", "COLON", "LPAREN", "RPAREN", "LBRACE", "RBRACE", "LBRACK", 
      "RBRACK", "LT", "GT", "BOOLEAN", "NUMBER", "STRING", "DNA_SEQUENCE", 
      "ID", "COMMENT", "BLOCK_COMMENT", "WS"
    }
  );
  static const int32_t serializedATNSegment[] = {
  	4,1,57,417,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,6,2,
  	7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,13,7,13,2,14,7,
  	14,2,15,7,15,2,16,7,16,2,17,7,17,2,18,7,18,2,19,7,19,2,20,7,20,2,21,7,
  	21,2,22,7,22,2,23,7,23,2,24,7,24,2,25,7,25,2,26,7,26,2,27,7,27,2,28,7,
  	28,2,29,7,29,2,30,7,30,2,31,7,31,2,32,7,32,2,33,7,33,1,0,5,0,70,8,0,10,
  	0,12,0,73,9,0,1,0,5,0,76,8,0,10,0,12,0,79,9,0,1,0,1,0,1,1,1,1,1,1,1,1,
  	3,1,87,8,1,1,1,1,1,1,2,1,2,1,2,5,2,94,8,2,10,2,12,2,97,9,2,1,3,1,3,1,
  	3,1,3,1,3,1,3,1,4,3,4,106,8,4,1,4,3,4,109,8,4,1,4,3,4,112,8,4,1,4,4,4,
  	115,8,4,11,4,12,4,116,1,5,1,5,1,5,1,5,1,5,1,5,1,5,5,5,126,8,5,10,5,12,
  	5,129,9,5,1,5,1,5,3,5,133,8,5,1,6,1,6,1,6,1,6,3,6,139,8,6,1,7,1,7,1,7,
  	1,7,1,7,1,7,1,7,5,7,148,8,7,10,7,12,7,151,9,7,1,7,1,7,3,7,155,8,7,1,8,
  	1,8,1,8,1,8,1,9,1,9,1,9,1,9,1,9,1,9,5,9,167,8,9,10,9,12,9,170,9,9,1,9,
  	1,9,1,10,1,10,1,10,1,10,1,10,1,10,1,10,1,11,5,11,182,8,11,10,11,12,11,
  	185,9,11,1,11,1,11,1,11,1,11,1,11,1,11,1,12,1,12,1,12,1,12,1,12,1,12,
  	3,12,199,8,12,1,13,1,13,1,13,5,13,204,8,13,10,13,12,13,207,9,13,1,14,
  	1,14,1,14,1,14,3,14,213,8,14,1,14,3,14,216,8,14,1,15,4,15,219,8,15,11,
  	15,12,15,220,1,16,1,16,1,16,1,16,1,16,1,16,1,16,3,16,230,8,16,1,17,1,
  	17,1,17,1,17,1,18,1,18,1,18,3,18,239,8,18,1,18,1,18,1,18,1,18,1,18,3,
  	18,246,8,18,1,19,1,19,1,19,5,19,251,8,19,10,19,12,19,254,9,19,1,20,1,
  	20,1,20,5,20,259,8,20,10,20,12,20,262,9,20,1,21,1,21,1,21,1,21,1,21,3,
  	21,269,8,21,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,22,1,23,1,
  	23,1,23,5,23,284,8,23,10,23,12,23,287,9,23,1,24,1,24,1,24,1,24,1,24,1,
  	24,1,24,1,24,1,24,1,24,3,24,299,8,24,1,25,1,25,1,25,1,25,1,25,1,25,1,
  	25,1,25,1,25,1,25,1,25,1,25,1,25,1,25,1,25,1,25,3,25,317,8,25,1,26,1,
  	26,1,26,1,26,1,26,5,26,324,8,26,10,26,12,26,327,9,26,1,27,1,27,1,27,4,
  	27,332,8,27,11,27,12,27,333,1,27,1,27,1,28,1,28,1,28,1,28,1,28,1,28,1,
  	28,1,28,1,28,1,28,3,28,348,8,28,1,29,1,29,1,29,1,29,1,29,1,29,1,29,1,
  	29,3,29,358,8,29,1,29,1,29,1,29,1,29,1,29,1,29,1,29,1,29,5,29,368,8,29,
  	10,29,12,29,371,9,29,1,30,1,30,1,30,3,30,376,8,30,1,30,1,30,1,30,3,30,
  	381,8,30,5,30,383,8,30,10,30,12,30,386,9,30,1,30,1,30,1,31,1,31,1,32,
  	1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,1,32,
  	1,32,1,32,1,32,1,32,1,32,1,32,3,32,413,8,32,1,33,1,33,1,33,0,1,58,34,
  	0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,
  	50,52,54,56,58,60,62,64,66,0,2,1,0,50,53,1,0,32,34,441,0,71,1,0,0,0,2,
  	82,1,0,0,0,4,90,1,0,0,0,6,98,1,0,0,0,8,105,1,0,0,0,10,118,1,0,0,0,12,
  	134,1,0,0,0,14,140,1,0,0,0,16,156,1,0,0,0,18,160,1,0,0,0,20,173,1,0,0,
  	0,22,183,1,0,0,0,24,192,1,0,0,0,26,200,1,0,0,0,28,215,1,0,0,0,30,218,
  	1,0,0,0,32,229,1,0,0,0,34,231,1,0,0,0,36,235,1,0,0,0,38,247,1,0,0,0,40,
  	255,1,0,0,0,42,263,1,0,0,0,44,270,1,0,0,0,46,280,1,0,0,0,48,288,1,0,0,
  	0,50,316,1,0,0,0,52,318,1,0,0,0,54,328,1,0,0,0,56,347,1,0,0,0,58,357,
  	1,0,0,0,60,372,1,0,0,0,62,389,1,0,0,0,64,412,1,0,0,0,66,414,1,0,0,0,68,
  	70,3,2,1,0,69,68,1,0,0,0,70,73,1,0,0,0,71,69,1,0,0,0,71,72,1,0,0,0,72,
  	77,1,0,0,0,73,71,1,0,0,0,74,76,3,6,3,0,75,74,1,0,0,0,76,79,1,0,0,0,77,
  	75,1,0,0,0,77,78,1,0,0,0,78,80,1,0,0,0,79,77,1,0,0,0,80,81,5,0,0,1,81,
  	1,1,0,0,0,82,83,5,7,0,0,83,86,3,4,2,0,84,85,5,8,0,0,85,87,5,54,0,0,86,
  	84,1,0,0,0,86,87,1,0,0,0,87,88,1,0,0,0,88,89,5,40,0,0,89,3,1,0,0,0,90,
  	95,5,54,0,0,91,92,5,38,0,0,92,94,5,54,0,0,93,91,1,0,0,0,94,97,1,0,0,0,
  	95,93,1,0,0,0,95,96,1,0,0,0,96,5,1,0,0,0,97,95,1,0,0,0,98,99,5,2,0,0,
  	99,100,5,54,0,0,100,101,5,44,0,0,101,102,3,8,4,0,102,103,5,45,0,0,103,
  	7,1,0,0,0,104,106,3,10,5,0,105,104,1,0,0,0,105,106,1,0,0,0,106,108,1,
  	0,0,0,107,109,3,14,7,0,108,107,1,0,0,0,108,109,1,0,0,0,109,111,1,0,0,
  	0,110,112,3,18,9,0,111,110,1,0,0,0,111,112,1,0,0,0,112,114,1,0,0,0,113,
  	115,3,22,11,0,114,113,1,0,0,0,115,116,1,0,0,0,116,114,1,0,0,0,116,117,
  	1,0,0,0,117,9,1,0,0,0,118,119,5,4,0,0,119,132,5,41,0,0,120,133,3,12,6,
  	0,121,122,5,44,0,0,122,127,3,12,6,0,123,124,5,39,0,0,124,126,3,12,6,0,
  	125,123,1,0,0,0,126,129,1,0,0,0,127,125,1,0,0,0,127,128,1,0,0,0,128,130,
  	1,0,0,0,129,127,1,0,0,0,130,131,5,45,0,0,131,133,1,0,0,0,132,120,1,0,
  	0,0,132,121,1,0,0,0,133,11,1,0,0,0,134,135,5,54,0,0,135,136,5,41,0,0,
  	136,138,3,64,32,0,137,139,5,20,0,0,138,137,1,0,0,0,138,139,1,0,0,0,139,
  	13,1,0,0,0,140,141,5,5,0,0,141,154,5,41,0,0,142,155,3,16,8,0,143,144,
  	5,44,0,0,144,149,3,16,8,0,145,146,5,39,0,0,146,148,3,16,8,0,147,145,1,
  	0,0,0,148,151,1,0,0,0,149,147,1,0,0,0,149,150,1,0,0,0,150,152,1,0,0,0,
  	151,149,1,0,0,0,152,153,5,45,0,0,153,155,1,0,0,0,154,142,1,0,0,0,154,
  	143,1,0,0,0,155,15,1,0,0,0,156,157,5,54,0,0,157,158,5,41,0,0,158,159,
  	3,64,32,0,159,17,1,0,0,0,160,161,5,6,0,0,161,162,5,41,0,0,162,163,5,44,
  	0,0,163,168,3,20,10,0,164,165,5,39,0,0,165,167,3,20,10,0,166,164,1,0,
  	0,0,167,170,1,0,0,0,168,166,1,0,0,0,168,169,1,0,0,0,169,171,1,0,0,0,170,
  	168,1,0,0,0,171,172,5,45,0,0,172,19,1,0,0,0,173,174,5,54,0,0,174,175,
  	5,41,0,0,175,176,3,66,33,0,176,177,5,42,0,0,177,178,5,52,0,0,178,179,
  	5,43,0,0,179,21,1,0,0,0,180,182,3,24,12,0,181,180,1,0,0,0,182,185,1,0,
  	0,0,183,181,1,0,0,0,183,184,1,0,0,0,184,186,1,0,0,0,185,183,1,0,0,0,186,
  	187,5,3,0,0,187,188,5,54,0,0,188,189,5,44,0,0,189,190,3,30,15,0,190,191,
  	5,45,0,0,191,23,1,0,0,0,192,193,5,1,0,0,193,198,5,54,0,0,194,195,5,42,
  	0,0,195,196,3,26,13,0,196,197,5,43,0,0,197,199,1,0,0,0,198,194,1,0,0,
  	0,198,199,1,0,0,0,199,25,1,0,0,0,200,205,3,28,14,0,201,202,5,39,0,0,202,
  	204,3,28,14,0,203,201,1,0,0,0,204,207,1,0,0,0,205,203,1,0,0,0,205,206,
  	1,0,0,0,206,27,1,0,0,0,207,205,1,0,0,0,208,209,5,54,0,0,209,212,5,41,
  	0,0,210,213,3,62,31,0,211,213,5,54,0,0,212,210,1,0,0,0,212,211,1,0,0,
  	0,213,216,1,0,0,0,214,216,3,62,31,0,215,208,1,0,0,0,215,214,1,0,0,0,216,
  	29,1,0,0,0,217,219,3,32,16,0,218,217,1,0,0,0,219,220,1,0,0,0,220,218,
  	1,0,0,0,220,221,1,0,0,0,221,31,1,0,0,0,222,230,3,34,17,0,223,230,3,36,
  	18,0,224,230,3,44,22,0,225,230,3,48,24,0,226,230,3,50,25,0,227,230,3,
  	52,26,0,228,230,3,54,27,0,229,222,1,0,0,0,229,223,1,0,0,0,229,224,1,0,
  	0,0,229,225,1,0,0,0,229,226,1,0,0,0,229,227,1,0,0,0,229,228,1,0,0,0,230,
  	33,1,0,0,0,231,232,5,54,0,0,232,233,5,37,0,0,233,234,3,58,29,0,234,35,
  	1,0,0,0,235,236,5,54,0,0,236,238,5,42,0,0,237,239,3,38,19,0,238,237,1,
  	0,0,0,238,239,1,0,0,0,239,240,1,0,0,0,240,245,5,43,0,0,241,242,5,44,0,
  	0,242,243,3,40,20,0,243,244,5,45,0,0,244,246,1,0,0,0,245,241,1,0,0,0,
  	245,246,1,0,0,0,246,37,1,0,0,0,247,252,3,58,29,0,248,249,5,39,0,0,249,
  	251,3,58,29,0,250,248,1,0,0,0,251,254,1,0,0,0,252,250,1,0,0,0,252,253,
  	1,0,0,0,253,39,1,0,0,0,254,252,1,0,0,0,255,260,3,42,21,0,256,257,5,39,
  	0,0,257,259,3,42,21,0,258,256,1,0,0,0,259,262,1,0,0,0,260,258,1,0,0,0,
  	260,261,1,0,0,0,261,41,1,0,0,0,262,260,1,0,0,0,263,264,5,54,0,0,264,268,
  	5,41,0,0,265,269,3,62,31,0,266,269,3,58,29,0,267,269,3,60,30,0,268,265,
  	1,0,0,0,268,266,1,0,0,0,268,267,1,0,0,0,269,43,1,0,0,0,270,271,5,19,0,
  	0,271,272,5,42,0,0,272,273,3,58,29,0,273,274,5,39,0,0,274,275,3,58,29,
  	0,275,276,5,43,0,0,276,277,5,44,0,0,277,278,3,46,23,0,278,279,5,45,0,
  	0,279,45,1,0,0,0,280,285,3,42,21,0,281,282,5,39,0,0,282,284,3,42,21,0,
  	283,281,1,0,0,0,284,287,1,0,0,0,285,283,1,0,0,0,285,286,1,0,0,0,286,47,
  	1,0,0,0,287,285,1,0,0,0,288,289,5,9,0,0,289,290,3,58,29,0,290,291,5,44,
  	0,0,291,292,3,30,15,0,292,298,5,45,0,0,293,294,5,10,0,0,294,295,5,44,
  	0,0,295,296,3,30,15,0,296,297,5,45,0,0,297,299,1,0,0,0,298,293,1,0,0,
  	0,298,299,1,0,0,0,299,49,1,0,0,0,300,301,5,11,0,0,301,302,5,54,0,0,302,
  	303,5,13,0,0,303,304,3,58,29,0,304,305,5,44,0,0,305,306,3,30,15,0,306,
  	307,5,45,0,0,307,317,1,0,0,0,308,309,5,12,0,0,309,310,5,54,0,0,310,311,
  	5,13,0,0,311,312,3,58,29,0,312,313,5,44,0,0,313,314,3,30,15,0,314,315,
  	5,45,0,0,315,317,1,0,0,0,316,300,1,0,0,0,316,308,1,0,0,0,317,51,1,0,0,
  	0,318,319,5,14,0,0,319,320,5,41,0,0,320,325,5,54,0,0,321,322,5,39,0,0,
  	322,324,5,54,0,0,323,321,1,0,0,0,324,327,1,0,0,0,325,323,1,0,0,0,325,
  	326,1,0,0,0,326,53,1,0,0,0,327,325,1,0,0,0,328,329,5,15,0,0,329,331,5,
  	44,0,0,330,332,3,56,28,0,331,330,1,0,0,0,332,333,1,0,0,0,333,331,1,0,
  	0,0,333,334,1,0,0,0,334,335,1,0,0,0,335,336,5,45,0,0,336,55,1,0,0,0,337,
  	338,5,16,0,0,338,348,5,52,0,0,339,340,5,17,0,0,340,341,5,41,0,0,341,348,
  	5,52,0,0,342,343,5,18,0,0,343,344,5,41,0,0,344,348,5,52,0,0,345,348,3,
  	48,24,0,346,348,3,50,25,0,347,337,1,0,0,0,347,339,1,0,0,0,347,342,1,0,
  	0,0,347,345,1,0,0,0,347,346,1,0,0,0,348,57,1,0,0,0,349,350,6,29,-1,0,
  	350,358,5,54,0,0,351,358,3,62,31,0,352,358,3,36,18,0,353,354,5,42,0,0,
  	354,355,3,58,29,0,355,356,5,43,0,0,356,358,1,0,0,0,357,349,1,0,0,0,357,
  	351,1,0,0,0,357,352,1,0,0,0,357,353,1,0,0,0,358,369,1,0,0,0,359,360,10,
  	3,0,0,360,361,5,38,0,0,361,368,5,54,0,0,362,363,10,2,0,0,363,364,5,46,
  	0,0,364,365,3,58,29,0,365,366,5,47,0,0,366,368,1,0,0,0,367,359,1,0,0,
  	0,367,362,1,0,0,0,368,371,1,0,0,0,369,367,1,0,0,0,369,370,1,0,0,0,370,
  	59,1,0,0,0,371,369,1,0,0,0,372,375,5,46,0,0,373,376,3,62,31,0,374,376,
  	5,54,0,0,375,373,1,0,0,0,375,374,1,0,0,0,376,384,1,0,0,0,377,380,5,39,
  	0,0,378,381,3,62,31,0,379,381,5,54,0,0,380,378,1,0,0,0,380,379,1,0,0,
  	0,381,383,1,0,0,0,382,377,1,0,0,0,383,386,1,0,0,0,384,382,1,0,0,0,384,
  	385,1,0,0,0,385,387,1,0,0,0,386,384,1,0,0,0,387,388,5,47,0,0,388,61,1,
  	0,0,0,389,390,7,0,0,0,390,63,1,0,0,0,391,413,5,21,0,0,392,413,5,22,0,
  	0,393,413,5,23,0,0,394,413,5,24,0,0,395,413,5,25,0,0,396,413,5,26,0,0,
  	397,413,5,27,0,0,398,413,5,28,0,0,399,413,5,29,0,0,400,401,5,30,0,0,401,
  	402,5,48,0,0,402,403,3,64,32,0,403,404,5,49,0,0,404,413,1,0,0,0,405,406,
  	5,31,0,0,406,407,5,48,0,0,407,408,3,64,32,0,408,409,5,39,0,0,409,410,
  	3,64,32,0,410,411,5,49,0,0,411,413,1,0,0,0,412,391,1,0,0,0,412,392,1,
  	0,0,0,412,393,1,0,0,0,412,394,1,0,0,0,412,395,1,0,0,0,412,396,1,0,0,0,
  	412,397,1,0,0,0,412,398,1,0,0,0,412,399,1,0,0,0,412,400,1,0,0,0,412,405,
  	1,0,0,0,413,65,1,0,0,0,414,415,7,1,0,0,415,67,1,0,0,0,39,71,77,86,95,
  	105,108,111,116,127,132,138,149,154,168,183,198,205,212,215,220,229,238,
  	245,252,260,268,285,298,316,325,333,347,357,367,369,375,380,384,412
  };
  staticData->serializedATN = antlr4::atn::SerializedATNView(serializedATNSegment, sizeof(serializedATNSegment) / sizeof(serializedATNSegment[0]));

  antlr4::atn::ATNDeserializer deserializer;
  staticData->atn = deserializer.deserialize(staticData->serializedATN);

  const size_t count = staticData->atn->getNumberOfDecisions();
  staticData->decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) { 
    staticData->decisionToDFA.emplace_back(staticData->atn->getDecisionState(i), i);
  }
  biogpuParserStaticData = staticData.release();
}

}

BioGPUParser::BioGPUParser(TokenStream *input) : BioGPUParser(input, antlr4::atn::ParserATNSimulatorOptions()) {}

BioGPUParser::BioGPUParser(TokenStream *input, const antlr4::atn::ParserATNSimulatorOptions &options) : Parser(input) {
  BioGPUParser::initialize();
  _interpreter = new atn::ParserATNSimulator(this, *biogpuParserStaticData->atn, biogpuParserStaticData->decisionToDFA, biogpuParserStaticData->sharedContextCache, options);
}

BioGPUParser::~BioGPUParser() {
  delete _interpreter;
}

const atn::ATN& BioGPUParser::getATN() const {
  return *biogpuParserStaticData->atn;
}

std::string BioGPUParser::getGrammarFileName() const {
  return "BioGPU.g4";
}

const std::vector<std::string>& BioGPUParser::getRuleNames() const {
  return biogpuParserStaticData->ruleNames;
}

const dfa::Vocabulary& BioGPUParser::getVocabulary() const {
  return biogpuParserStaticData->vocabulary;
}

antlr4::atn::SerializedATNView BioGPUParser::getSerializedATN() const {
  return biogpuParserStaticData->serializedATN;
}


//----------------- ProgramContext ------------------------------------------------------------------

BioGPUParser::ProgramContext::ProgramContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ProgramContext::EOF() {
  return getToken(BioGPUParser::EOF, 0);
}

std::vector<BioGPUParser::ImportStatementContext *> BioGPUParser::ProgramContext::importStatement() {
  return getRuleContexts<BioGPUParser::ImportStatementContext>();
}

BioGPUParser::ImportStatementContext* BioGPUParser::ProgramContext::importStatement(size_t i) {
  return getRuleContext<BioGPUParser::ImportStatementContext>(i);
}

std::vector<BioGPUParser::PipelineContext *> BioGPUParser::ProgramContext::pipeline() {
  return getRuleContexts<BioGPUParser::PipelineContext>();
}

BioGPUParser::PipelineContext* BioGPUParser::ProgramContext::pipeline(size_t i) {
  return getRuleContext<BioGPUParser::PipelineContext>(i);
}


size_t BioGPUParser::ProgramContext::getRuleIndex() const {
  return BioGPUParser::RuleProgram;
}


std::any BioGPUParser::ProgramContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitProgram(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ProgramContext* BioGPUParser::program() {
  ProgramContext *_localctx = _tracker.createInstance<ProgramContext>(_ctx, getState());
  enterRule(_localctx, 0, BioGPUParser::RuleProgram);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(71);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::IMPORT) {
      setState(68);
      importStatement();
      setState(73);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(77);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::PIPELINE) {
      setState(74);
      pipeline();
      setState(79);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(80);
    match(BioGPUParser::EOF);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ImportStatementContext ------------------------------------------------------------------

BioGPUParser::ImportStatementContext::ImportStatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ImportStatementContext::IMPORT() {
  return getToken(BioGPUParser::IMPORT, 0);
}

BioGPUParser::ModulePathContext* BioGPUParser::ImportStatementContext::modulePath() {
  return getRuleContext<BioGPUParser::ModulePathContext>(0);
}

tree::TerminalNode* BioGPUParser::ImportStatementContext::SEMI() {
  return getToken(BioGPUParser::SEMI, 0);
}

tree::TerminalNode* BioGPUParser::ImportStatementContext::AS() {
  return getToken(BioGPUParser::AS, 0);
}

tree::TerminalNode* BioGPUParser::ImportStatementContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}


size_t BioGPUParser::ImportStatementContext::getRuleIndex() const {
  return BioGPUParser::RuleImportStatement;
}


std::any BioGPUParser::ImportStatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitImportStatement(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ImportStatementContext* BioGPUParser::importStatement() {
  ImportStatementContext *_localctx = _tracker.createInstance<ImportStatementContext>(_ctx, getState());
  enterRule(_localctx, 2, BioGPUParser::RuleImportStatement);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(82);
    match(BioGPUParser::IMPORT);
    setState(83);
    modulePath();
    setState(86);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::AS) {
      setState(84);
      match(BioGPUParser::AS);
      setState(85);
      match(BioGPUParser::ID);
    }
    setState(88);
    match(BioGPUParser::SEMI);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ModulePathContext ------------------------------------------------------------------

BioGPUParser::ModulePathContext::ModulePathContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<tree::TerminalNode *> BioGPUParser::ModulePathContext::ID() {
  return getTokens(BioGPUParser::ID);
}

tree::TerminalNode* BioGPUParser::ModulePathContext::ID(size_t i) {
  return getToken(BioGPUParser::ID, i);
}

std::vector<tree::TerminalNode *> BioGPUParser::ModulePathContext::DOT() {
  return getTokens(BioGPUParser::DOT);
}

tree::TerminalNode* BioGPUParser::ModulePathContext::DOT(size_t i) {
  return getToken(BioGPUParser::DOT, i);
}


size_t BioGPUParser::ModulePathContext::getRuleIndex() const {
  return BioGPUParser::RuleModulePath;
}


std::any BioGPUParser::ModulePathContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitModulePath(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ModulePathContext* BioGPUParser::modulePath() {
  ModulePathContext *_localctx = _tracker.createInstance<ModulePathContext>(_ctx, getState());
  enterRule(_localctx, 4, BioGPUParser::RuleModulePath);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(90);
    match(BioGPUParser::ID);
    setState(95);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::DOT) {
      setState(91);
      match(BioGPUParser::DOT);
      setState(92);
      match(BioGPUParser::ID);
      setState(97);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- PipelineContext ------------------------------------------------------------------

BioGPUParser::PipelineContext::PipelineContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::PipelineContext::PIPELINE() {
  return getToken(BioGPUParser::PIPELINE, 0);
}

tree::TerminalNode* BioGPUParser::PipelineContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::PipelineContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

BioGPUParser::PipelineBodyContext* BioGPUParser::PipelineContext::pipelineBody() {
  return getRuleContext<BioGPUParser::PipelineBodyContext>(0);
}

tree::TerminalNode* BioGPUParser::PipelineContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}


size_t BioGPUParser::PipelineContext::getRuleIndex() const {
  return BioGPUParser::RulePipeline;
}


std::any BioGPUParser::PipelineContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitPipeline(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::PipelineContext* BioGPUParser::pipeline() {
  PipelineContext *_localctx = _tracker.createInstance<PipelineContext>(_ctx, getState());
  enterRule(_localctx, 6, BioGPUParser::RulePipeline);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(98);
    match(BioGPUParser::PIPELINE);
    setState(99);
    match(BioGPUParser::ID);
    setState(100);
    match(BioGPUParser::LBRACE);
    setState(101);
    pipelineBody();
    setState(102);
    match(BioGPUParser::RBRACE);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- PipelineBodyContext ------------------------------------------------------------------

BioGPUParser::PipelineBodyContext::PipelineBodyContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

BioGPUParser::InputDeclContext* BioGPUParser::PipelineBodyContext::inputDecl() {
  return getRuleContext<BioGPUParser::InputDeclContext>(0);
}

BioGPUParser::OutputDeclContext* BioGPUParser::PipelineBodyContext::outputDecl() {
  return getRuleContext<BioGPUParser::OutputDeclContext>(0);
}

BioGPUParser::ReferencesDeclContext* BioGPUParser::PipelineBodyContext::referencesDecl() {
  return getRuleContext<BioGPUParser::ReferencesDeclContext>(0);
}

std::vector<BioGPUParser::StageContext *> BioGPUParser::PipelineBodyContext::stage() {
  return getRuleContexts<BioGPUParser::StageContext>();
}

BioGPUParser::StageContext* BioGPUParser::PipelineBodyContext::stage(size_t i) {
  return getRuleContext<BioGPUParser::StageContext>(i);
}


size_t BioGPUParser::PipelineBodyContext::getRuleIndex() const {
  return BioGPUParser::RulePipelineBody;
}


std::any BioGPUParser::PipelineBodyContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitPipelineBody(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::PipelineBodyContext* BioGPUParser::pipelineBody() {
  PipelineBodyContext *_localctx = _tracker.createInstance<PipelineBodyContext>(_ctx, getState());
  enterRule(_localctx, 8, BioGPUParser::RulePipelineBody);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(105);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::INPUT) {
      setState(104);
      inputDecl();
    }
    setState(108);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::OUTPUT) {
      setState(107);
      outputDecl();
    }
    setState(111);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::REFERENCES) {
      setState(110);
      referencesDecl();
    }
    setState(114); 
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(113);
      stage();
      setState(116); 
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while (_la == BioGPUParser::T__0

    || _la == BioGPUParser::STAGE);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- InputDeclContext ------------------------------------------------------------------

BioGPUParser::InputDeclContext::InputDeclContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::InputDeclContext::INPUT() {
  return getToken(BioGPUParser::INPUT, 0);
}

tree::TerminalNode* BioGPUParser::InputDeclContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

std::vector<BioGPUParser::InputParamContext *> BioGPUParser::InputDeclContext::inputParam() {
  return getRuleContexts<BioGPUParser::InputParamContext>();
}

BioGPUParser::InputParamContext* BioGPUParser::InputDeclContext::inputParam(size_t i) {
  return getRuleContext<BioGPUParser::InputParamContext>(i);
}

tree::TerminalNode* BioGPUParser::InputDeclContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

tree::TerminalNode* BioGPUParser::InputDeclContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}

std::vector<tree::TerminalNode *> BioGPUParser::InputDeclContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::InputDeclContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::InputDeclContext::getRuleIndex() const {
  return BioGPUParser::RuleInputDecl;
}


std::any BioGPUParser::InputDeclContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitInputDecl(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::InputDeclContext* BioGPUParser::inputDecl() {
  InputDeclContext *_localctx = _tracker.createInstance<InputDeclContext>(_ctx, getState());
  enterRule(_localctx, 10, BioGPUParser::RuleInputDecl);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(118);
    match(BioGPUParser::INPUT);
    setState(119);
    match(BioGPUParser::COLON);
    setState(132);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::ID: {
        setState(120);
        inputParam();
        break;
      }

      case BioGPUParser::LBRACE: {
        setState(121);
        match(BioGPUParser::LBRACE);
        setState(122);
        inputParam();
        setState(127);
        _errHandler->sync(this);
        _la = _input->LA(1);
        while (_la == BioGPUParser::COMMA) {
          setState(123);
          match(BioGPUParser::COMMA);
          setState(124);
          inputParam();
          setState(129);
          _errHandler->sync(this);
          _la = _input->LA(1);
        }
        setState(130);
        match(BioGPUParser::RBRACE);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- InputParamContext ------------------------------------------------------------------

BioGPUParser::InputParamContext::InputParamContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::InputParamContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::InputParamContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

BioGPUParser::DataTypeContext* BioGPUParser::InputParamContext::dataType() {
  return getRuleContext<BioGPUParser::DataTypeContext>(0);
}

tree::TerminalNode* BioGPUParser::InputParamContext::OPTIONAL() {
  return getToken(BioGPUParser::OPTIONAL, 0);
}


size_t BioGPUParser::InputParamContext::getRuleIndex() const {
  return BioGPUParser::RuleInputParam;
}


std::any BioGPUParser::InputParamContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitInputParam(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::InputParamContext* BioGPUParser::inputParam() {
  InputParamContext *_localctx = _tracker.createInstance<InputParamContext>(_ctx, getState());
  enterRule(_localctx, 12, BioGPUParser::RuleInputParam);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(134);
    match(BioGPUParser::ID);
    setState(135);
    match(BioGPUParser::COLON);
    setState(136);
    dataType();
    setState(138);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::OPTIONAL) {
      setState(137);
      match(BioGPUParser::OPTIONAL);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- OutputDeclContext ------------------------------------------------------------------

BioGPUParser::OutputDeclContext::OutputDeclContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::OutputDeclContext::OUTPUT() {
  return getToken(BioGPUParser::OUTPUT, 0);
}

tree::TerminalNode* BioGPUParser::OutputDeclContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

std::vector<BioGPUParser::OutputParamContext *> BioGPUParser::OutputDeclContext::outputParam() {
  return getRuleContexts<BioGPUParser::OutputParamContext>();
}

BioGPUParser::OutputParamContext* BioGPUParser::OutputDeclContext::outputParam(size_t i) {
  return getRuleContext<BioGPUParser::OutputParamContext>(i);
}

tree::TerminalNode* BioGPUParser::OutputDeclContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

tree::TerminalNode* BioGPUParser::OutputDeclContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}

std::vector<tree::TerminalNode *> BioGPUParser::OutputDeclContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::OutputDeclContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::OutputDeclContext::getRuleIndex() const {
  return BioGPUParser::RuleOutputDecl;
}


std::any BioGPUParser::OutputDeclContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitOutputDecl(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::OutputDeclContext* BioGPUParser::outputDecl() {
  OutputDeclContext *_localctx = _tracker.createInstance<OutputDeclContext>(_ctx, getState());
  enterRule(_localctx, 14, BioGPUParser::RuleOutputDecl);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(140);
    match(BioGPUParser::OUTPUT);
    setState(141);
    match(BioGPUParser::COLON);
    setState(154);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::ID: {
        setState(142);
        outputParam();
        break;
      }

      case BioGPUParser::LBRACE: {
        setState(143);
        match(BioGPUParser::LBRACE);
        setState(144);
        outputParam();
        setState(149);
        _errHandler->sync(this);
        _la = _input->LA(1);
        while (_la == BioGPUParser::COMMA) {
          setState(145);
          match(BioGPUParser::COMMA);
          setState(146);
          outputParam();
          setState(151);
          _errHandler->sync(this);
          _la = _input->LA(1);
        }
        setState(152);
        match(BioGPUParser::RBRACE);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- OutputParamContext ------------------------------------------------------------------

BioGPUParser::OutputParamContext::OutputParamContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::OutputParamContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::OutputParamContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

BioGPUParser::DataTypeContext* BioGPUParser::OutputParamContext::dataType() {
  return getRuleContext<BioGPUParser::DataTypeContext>(0);
}


size_t BioGPUParser::OutputParamContext::getRuleIndex() const {
  return BioGPUParser::RuleOutputParam;
}


std::any BioGPUParser::OutputParamContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitOutputParam(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::OutputParamContext* BioGPUParser::outputParam() {
  OutputParamContext *_localctx = _tracker.createInstance<OutputParamContext>(_ctx, getState());
  enterRule(_localctx, 16, BioGPUParser::RuleOutputParam);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(156);
    match(BioGPUParser::ID);
    setState(157);
    match(BioGPUParser::COLON);
    setState(158);
    dataType();
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ReferencesDeclContext ------------------------------------------------------------------

BioGPUParser::ReferencesDeclContext::ReferencesDeclContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ReferencesDeclContext::REFERENCES() {
  return getToken(BioGPUParser::REFERENCES, 0);
}

tree::TerminalNode* BioGPUParser::ReferencesDeclContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

tree::TerminalNode* BioGPUParser::ReferencesDeclContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

std::vector<BioGPUParser::ReferenceParamContext *> BioGPUParser::ReferencesDeclContext::referenceParam() {
  return getRuleContexts<BioGPUParser::ReferenceParamContext>();
}

BioGPUParser::ReferenceParamContext* BioGPUParser::ReferencesDeclContext::referenceParam(size_t i) {
  return getRuleContext<BioGPUParser::ReferenceParamContext>(i);
}

tree::TerminalNode* BioGPUParser::ReferencesDeclContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}

std::vector<tree::TerminalNode *> BioGPUParser::ReferencesDeclContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::ReferencesDeclContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::ReferencesDeclContext::getRuleIndex() const {
  return BioGPUParser::RuleReferencesDecl;
}


std::any BioGPUParser::ReferencesDeclContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitReferencesDecl(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ReferencesDeclContext* BioGPUParser::referencesDecl() {
  ReferencesDeclContext *_localctx = _tracker.createInstance<ReferencesDeclContext>(_ctx, getState());
  enterRule(_localctx, 18, BioGPUParser::RuleReferencesDecl);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(160);
    match(BioGPUParser::REFERENCES);
    setState(161);
    match(BioGPUParser::COLON);
    setState(162);
    match(BioGPUParser::LBRACE);
    setState(163);
    referenceParam();
    setState(168);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(164);
      match(BioGPUParser::COMMA);
      setState(165);
      referenceParam();
      setState(170);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(171);
    match(BioGPUParser::RBRACE);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ReferenceParamContext ------------------------------------------------------------------

BioGPUParser::ReferenceParamContext::ReferenceParamContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ReferenceParamContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::ReferenceParamContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

BioGPUParser::ReferenceTypeContext* BioGPUParser::ReferenceParamContext::referenceType() {
  return getRuleContext<BioGPUParser::ReferenceTypeContext>(0);
}

tree::TerminalNode* BioGPUParser::ReferenceParamContext::LPAREN() {
  return getToken(BioGPUParser::LPAREN, 0);
}

tree::TerminalNode* BioGPUParser::ReferenceParamContext::STRING() {
  return getToken(BioGPUParser::STRING, 0);
}

tree::TerminalNode* BioGPUParser::ReferenceParamContext::RPAREN() {
  return getToken(BioGPUParser::RPAREN, 0);
}


size_t BioGPUParser::ReferenceParamContext::getRuleIndex() const {
  return BioGPUParser::RuleReferenceParam;
}


std::any BioGPUParser::ReferenceParamContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitReferenceParam(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ReferenceParamContext* BioGPUParser::referenceParam() {
  ReferenceParamContext *_localctx = _tracker.createInstance<ReferenceParamContext>(_ctx, getState());
  enterRule(_localctx, 20, BioGPUParser::RuleReferenceParam);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(173);
    match(BioGPUParser::ID);
    setState(174);
    match(BioGPUParser::COLON);
    setState(175);
    referenceType();
    setState(176);
    match(BioGPUParser::LPAREN);
    setState(177);
    match(BioGPUParser::STRING);
    setState(178);
    match(BioGPUParser::RPAREN);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- StageContext ------------------------------------------------------------------

BioGPUParser::StageContext::StageContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::StageContext::STAGE() {
  return getToken(BioGPUParser::STAGE, 0);
}

tree::TerminalNode* BioGPUParser::StageContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::StageContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

BioGPUParser::StageBodyContext* BioGPUParser::StageContext::stageBody() {
  return getRuleContext<BioGPUParser::StageBodyContext>(0);
}

tree::TerminalNode* BioGPUParser::StageContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}

std::vector<BioGPUParser::DecoratorContext *> BioGPUParser::StageContext::decorator() {
  return getRuleContexts<BioGPUParser::DecoratorContext>();
}

BioGPUParser::DecoratorContext* BioGPUParser::StageContext::decorator(size_t i) {
  return getRuleContext<BioGPUParser::DecoratorContext>(i);
}


size_t BioGPUParser::StageContext::getRuleIndex() const {
  return BioGPUParser::RuleStage;
}


std::any BioGPUParser::StageContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitStage(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::StageContext* BioGPUParser::stage() {
  StageContext *_localctx = _tracker.createInstance<StageContext>(_ctx, getState());
  enterRule(_localctx, 22, BioGPUParser::RuleStage);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(183);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::T__0) {
      setState(180);
      decorator();
      setState(185);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(186);
    match(BioGPUParser::STAGE);
    setState(187);
    match(BioGPUParser::ID);
    setState(188);
    match(BioGPUParser::LBRACE);
    setState(189);
    stageBody();
    setState(190);
    match(BioGPUParser::RBRACE);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DecoratorContext ------------------------------------------------------------------

BioGPUParser::DecoratorContext::DecoratorContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::DecoratorContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::DecoratorContext::LPAREN() {
  return getToken(BioGPUParser::LPAREN, 0);
}

BioGPUParser::DecoratorArgsContext* BioGPUParser::DecoratorContext::decoratorArgs() {
  return getRuleContext<BioGPUParser::DecoratorArgsContext>(0);
}

tree::TerminalNode* BioGPUParser::DecoratorContext::RPAREN() {
  return getToken(BioGPUParser::RPAREN, 0);
}


size_t BioGPUParser::DecoratorContext::getRuleIndex() const {
  return BioGPUParser::RuleDecorator;
}


std::any BioGPUParser::DecoratorContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitDecorator(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::DecoratorContext* BioGPUParser::decorator() {
  DecoratorContext *_localctx = _tracker.createInstance<DecoratorContext>(_ctx, getState());
  enterRule(_localctx, 24, BioGPUParser::RuleDecorator);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(192);
    match(BioGPUParser::T__0);
    setState(193);
    match(BioGPUParser::ID);
    setState(198);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::LPAREN) {
      setState(194);
      match(BioGPUParser::LPAREN);
      setState(195);
      decoratorArgs();
      setState(196);
      match(BioGPUParser::RPAREN);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DecoratorArgsContext ------------------------------------------------------------------

BioGPUParser::DecoratorArgsContext::DecoratorArgsContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<BioGPUParser::DecoratorArgContext *> BioGPUParser::DecoratorArgsContext::decoratorArg() {
  return getRuleContexts<BioGPUParser::DecoratorArgContext>();
}

BioGPUParser::DecoratorArgContext* BioGPUParser::DecoratorArgsContext::decoratorArg(size_t i) {
  return getRuleContext<BioGPUParser::DecoratorArgContext>(i);
}

std::vector<tree::TerminalNode *> BioGPUParser::DecoratorArgsContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::DecoratorArgsContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::DecoratorArgsContext::getRuleIndex() const {
  return BioGPUParser::RuleDecoratorArgs;
}


std::any BioGPUParser::DecoratorArgsContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitDecoratorArgs(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::DecoratorArgsContext* BioGPUParser::decoratorArgs() {
  DecoratorArgsContext *_localctx = _tracker.createInstance<DecoratorArgsContext>(_ctx, getState());
  enterRule(_localctx, 26, BioGPUParser::RuleDecoratorArgs);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(200);
    decoratorArg();
    setState(205);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(201);
      match(BioGPUParser::COMMA);
      setState(202);
      decoratorArg();
      setState(207);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DecoratorArgContext ------------------------------------------------------------------

BioGPUParser::DecoratorArgContext::DecoratorArgContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<tree::TerminalNode *> BioGPUParser::DecoratorArgContext::ID() {
  return getTokens(BioGPUParser::ID);
}

tree::TerminalNode* BioGPUParser::DecoratorArgContext::ID(size_t i) {
  return getToken(BioGPUParser::ID, i);
}

tree::TerminalNode* BioGPUParser::DecoratorArgContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

BioGPUParser::LiteralContext* BioGPUParser::DecoratorArgContext::literal() {
  return getRuleContext<BioGPUParser::LiteralContext>(0);
}


size_t BioGPUParser::DecoratorArgContext::getRuleIndex() const {
  return BioGPUParser::RuleDecoratorArg;
}


std::any BioGPUParser::DecoratorArgContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitDecoratorArg(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::DecoratorArgContext* BioGPUParser::decoratorArg() {
  DecoratorArgContext *_localctx = _tracker.createInstance<DecoratorArgContext>(_ctx, getState());
  enterRule(_localctx, 28, BioGPUParser::RuleDecoratorArg);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(215);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::ID: {
        enterOuterAlt(_localctx, 1);
        setState(208);
        match(BioGPUParser::ID);
        setState(209);
        match(BioGPUParser::COLON);
        setState(212);
        _errHandler->sync(this);
        switch (_input->LA(1)) {
          case BioGPUParser::BOOLEAN:
          case BioGPUParser::NUMBER:
          case BioGPUParser::STRING:
          case BioGPUParser::DNA_SEQUENCE: {
            setState(210);
            literal();
            break;
          }

          case BioGPUParser::ID: {
            setState(211);
            match(BioGPUParser::ID);
            break;
          }

        default:
          throw NoViableAltException(this);
        }
        break;
      }

      case BioGPUParser::BOOLEAN:
      case BioGPUParser::NUMBER:
      case BioGPUParser::STRING:
      case BioGPUParser::DNA_SEQUENCE: {
        enterOuterAlt(_localctx, 2);
        setState(214);
        literal();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- StageBodyContext ------------------------------------------------------------------

BioGPUParser::StageBodyContext::StageBodyContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<BioGPUParser::StatementContext *> BioGPUParser::StageBodyContext::statement() {
  return getRuleContexts<BioGPUParser::StatementContext>();
}

BioGPUParser::StatementContext* BioGPUParser::StageBodyContext::statement(size_t i) {
  return getRuleContext<BioGPUParser::StatementContext>(i);
}


size_t BioGPUParser::StageBodyContext::getRuleIndex() const {
  return BioGPUParser::RuleStageBody;
}


std::any BioGPUParser::StageBodyContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitStageBody(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::StageBodyContext* BioGPUParser::stageBody() {
  StageBodyContext *_localctx = _tracker.createInstance<StageBodyContext>(_ctx, getState());
  enterRule(_localctx, 30, BioGPUParser::RuleStageBody);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(218); 
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(217);
      statement();
      setState(220); 
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << BioGPUParser::IF)
      | (1ULL << BioGPUParser::FOR)
      | (1ULL << BioGPUParser::FOREACH)
      | (1ULL << BioGPUParser::EMIT)
      | (1ULL << BioGPUParser::REPORT)
      | (1ULL << BioGPUParser::PARALLEL_MAP)
      | (1ULL << BioGPUParser::ID))) != 0));
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- StatementContext ------------------------------------------------------------------

BioGPUParser::StatementContext::StatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

BioGPUParser::AssignmentContext* BioGPUParser::StatementContext::assignment() {
  return getRuleContext<BioGPUParser::AssignmentContext>(0);
}

BioGPUParser::FunctionCallContext* BioGPUParser::StatementContext::functionCall() {
  return getRuleContext<BioGPUParser::FunctionCallContext>(0);
}

BioGPUParser::ParallelMapContext* BioGPUParser::StatementContext::parallelMap() {
  return getRuleContext<BioGPUParser::ParallelMapContext>(0);
}

BioGPUParser::IfStatementContext* BioGPUParser::StatementContext::ifStatement() {
  return getRuleContext<BioGPUParser::IfStatementContext>(0);
}

BioGPUParser::ForStatementContext* BioGPUParser::StatementContext::forStatement() {
  return getRuleContext<BioGPUParser::ForStatementContext>(0);
}

BioGPUParser::EmitStatementContext* BioGPUParser::StatementContext::emitStatement() {
  return getRuleContext<BioGPUParser::EmitStatementContext>(0);
}

BioGPUParser::ReportBlockContext* BioGPUParser::StatementContext::reportBlock() {
  return getRuleContext<BioGPUParser::ReportBlockContext>(0);
}


size_t BioGPUParser::StatementContext::getRuleIndex() const {
  return BioGPUParser::RuleStatement;
}


std::any BioGPUParser::StatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitStatement(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::StatementContext* BioGPUParser::statement() {
  StatementContext *_localctx = _tracker.createInstance<StatementContext>(_ctx, getState());
  enterRule(_localctx, 32, BioGPUParser::RuleStatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(229);
    _errHandler->sync(this);
    switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 20, _ctx)) {
    case 1: {
      enterOuterAlt(_localctx, 1);
      setState(222);
      assignment();
      break;
    }

    case 2: {
      enterOuterAlt(_localctx, 2);
      setState(223);
      functionCall();
      break;
    }

    case 3: {
      enterOuterAlt(_localctx, 3);
      setState(224);
      parallelMap();
      break;
    }

    case 4: {
      enterOuterAlt(_localctx, 4);
      setState(225);
      ifStatement();
      break;
    }

    case 5: {
      enterOuterAlt(_localctx, 5);
      setState(226);
      forStatement();
      break;
    }

    case 6: {
      enterOuterAlt(_localctx, 6);
      setState(227);
      emitStatement();
      break;
    }

    case 7: {
      enterOuterAlt(_localctx, 7);
      setState(228);
      reportBlock();
      break;
    }

    default:
      break;
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- AssignmentContext ------------------------------------------------------------------

BioGPUParser::AssignmentContext::AssignmentContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::AssignmentContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::AssignmentContext::ASSIGN() {
  return getToken(BioGPUParser::ASSIGN, 0);
}

BioGPUParser::ExpressionContext* BioGPUParser::AssignmentContext::expression() {
  return getRuleContext<BioGPUParser::ExpressionContext>(0);
}


size_t BioGPUParser::AssignmentContext::getRuleIndex() const {
  return BioGPUParser::RuleAssignment;
}


std::any BioGPUParser::AssignmentContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitAssignment(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::AssignmentContext* BioGPUParser::assignment() {
  AssignmentContext *_localctx = _tracker.createInstance<AssignmentContext>(_ctx, getState());
  enterRule(_localctx, 34, BioGPUParser::RuleAssignment);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(231);
    match(BioGPUParser::ID);
    setState(232);
    match(BioGPUParser::ASSIGN);
    setState(233);
    expression(0);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- FunctionCallContext ------------------------------------------------------------------

BioGPUParser::FunctionCallContext::FunctionCallContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::FunctionCallContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::FunctionCallContext::LPAREN() {
  return getToken(BioGPUParser::LPAREN, 0);
}

tree::TerminalNode* BioGPUParser::FunctionCallContext::RPAREN() {
  return getToken(BioGPUParser::RPAREN, 0);
}

BioGPUParser::ArgumentListContext* BioGPUParser::FunctionCallContext::argumentList() {
  return getRuleContext<BioGPUParser::ArgumentListContext>(0);
}

tree::TerminalNode* BioGPUParser::FunctionCallContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

BioGPUParser::ConfigBlockContext* BioGPUParser::FunctionCallContext::configBlock() {
  return getRuleContext<BioGPUParser::ConfigBlockContext>(0);
}

tree::TerminalNode* BioGPUParser::FunctionCallContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}


size_t BioGPUParser::FunctionCallContext::getRuleIndex() const {
  return BioGPUParser::RuleFunctionCall;
}


std::any BioGPUParser::FunctionCallContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitFunctionCall(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::FunctionCallContext* BioGPUParser::functionCall() {
  FunctionCallContext *_localctx = _tracker.createInstance<FunctionCallContext>(_ctx, getState());
  enterRule(_localctx, 36, BioGPUParser::RuleFunctionCall);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(235);
    match(BioGPUParser::ID);
    setState(236);
    match(BioGPUParser::LPAREN);
    setState(238);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << BioGPUParser::LPAREN)
      | (1ULL << BioGPUParser::BOOLEAN)
      | (1ULL << BioGPUParser::NUMBER)
      | (1ULL << BioGPUParser::STRING)
      | (1ULL << BioGPUParser::DNA_SEQUENCE)
      | (1ULL << BioGPUParser::ID))) != 0)) {
      setState(237);
      argumentList();
    }
    setState(240);
    match(BioGPUParser::RPAREN);
    setState(245);
    _errHandler->sync(this);

    switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 22, _ctx)) {
    case 1: {
      setState(241);
      match(BioGPUParser::LBRACE);
      setState(242);
      configBlock();
      setState(243);
      match(BioGPUParser::RBRACE);
      break;
    }

    default:
      break;
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ArgumentListContext ------------------------------------------------------------------

BioGPUParser::ArgumentListContext::ArgumentListContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<BioGPUParser::ExpressionContext *> BioGPUParser::ArgumentListContext::expression() {
  return getRuleContexts<BioGPUParser::ExpressionContext>();
}

BioGPUParser::ExpressionContext* BioGPUParser::ArgumentListContext::expression(size_t i) {
  return getRuleContext<BioGPUParser::ExpressionContext>(i);
}

std::vector<tree::TerminalNode *> BioGPUParser::ArgumentListContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::ArgumentListContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::ArgumentListContext::getRuleIndex() const {
  return BioGPUParser::RuleArgumentList;
}


std::any BioGPUParser::ArgumentListContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitArgumentList(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ArgumentListContext* BioGPUParser::argumentList() {
  ArgumentListContext *_localctx = _tracker.createInstance<ArgumentListContext>(_ctx, getState());
  enterRule(_localctx, 38, BioGPUParser::RuleArgumentList);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(247);
    expression(0);
    setState(252);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(248);
      match(BioGPUParser::COMMA);
      setState(249);
      expression(0);
      setState(254);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ConfigBlockContext ------------------------------------------------------------------

BioGPUParser::ConfigBlockContext::ConfigBlockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<BioGPUParser::ConfigParamContext *> BioGPUParser::ConfigBlockContext::configParam() {
  return getRuleContexts<BioGPUParser::ConfigParamContext>();
}

BioGPUParser::ConfigParamContext* BioGPUParser::ConfigBlockContext::configParam(size_t i) {
  return getRuleContext<BioGPUParser::ConfigParamContext>(i);
}

std::vector<tree::TerminalNode *> BioGPUParser::ConfigBlockContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::ConfigBlockContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::ConfigBlockContext::getRuleIndex() const {
  return BioGPUParser::RuleConfigBlock;
}


std::any BioGPUParser::ConfigBlockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitConfigBlock(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ConfigBlockContext* BioGPUParser::configBlock() {
  ConfigBlockContext *_localctx = _tracker.createInstance<ConfigBlockContext>(_ctx, getState());
  enterRule(_localctx, 40, BioGPUParser::RuleConfigBlock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(255);
    configParam();
    setState(260);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(256);
      match(BioGPUParser::COMMA);
      setState(257);
      configParam();
      setState(262);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ConfigParamContext ------------------------------------------------------------------

BioGPUParser::ConfigParamContext::ConfigParamContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ConfigParamContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::ConfigParamContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

BioGPUParser::LiteralContext* BioGPUParser::ConfigParamContext::literal() {
  return getRuleContext<BioGPUParser::LiteralContext>(0);
}

BioGPUParser::ExpressionContext* BioGPUParser::ConfigParamContext::expression() {
  return getRuleContext<BioGPUParser::ExpressionContext>(0);
}

BioGPUParser::ArrayContext* BioGPUParser::ConfigParamContext::array() {
  return getRuleContext<BioGPUParser::ArrayContext>(0);
}


size_t BioGPUParser::ConfigParamContext::getRuleIndex() const {
  return BioGPUParser::RuleConfigParam;
}


std::any BioGPUParser::ConfigParamContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitConfigParam(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ConfigParamContext* BioGPUParser::configParam() {
  ConfigParamContext *_localctx = _tracker.createInstance<ConfigParamContext>(_ctx, getState());
  enterRule(_localctx, 42, BioGPUParser::RuleConfigParam);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(263);
    match(BioGPUParser::ID);
    setState(264);
    match(BioGPUParser::COLON);
    setState(268);
    _errHandler->sync(this);
    switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 25, _ctx)) {
    case 1: {
      setState(265);
      literal();
      break;
    }

    case 2: {
      setState(266);
      expression(0);
      break;
    }

    case 3: {
      setState(267);
      array();
      break;
    }

    default:
      break;
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ParallelMapContext ------------------------------------------------------------------

BioGPUParser::ParallelMapContext::ParallelMapContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ParallelMapContext::PARALLEL_MAP() {
  return getToken(BioGPUParser::PARALLEL_MAP, 0);
}

tree::TerminalNode* BioGPUParser::ParallelMapContext::LPAREN() {
  return getToken(BioGPUParser::LPAREN, 0);
}

std::vector<BioGPUParser::ExpressionContext *> BioGPUParser::ParallelMapContext::expression() {
  return getRuleContexts<BioGPUParser::ExpressionContext>();
}

BioGPUParser::ExpressionContext* BioGPUParser::ParallelMapContext::expression(size_t i) {
  return getRuleContext<BioGPUParser::ExpressionContext>(i);
}

tree::TerminalNode* BioGPUParser::ParallelMapContext::COMMA() {
  return getToken(BioGPUParser::COMMA, 0);
}

tree::TerminalNode* BioGPUParser::ParallelMapContext::RPAREN() {
  return getToken(BioGPUParser::RPAREN, 0);
}

tree::TerminalNode* BioGPUParser::ParallelMapContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

BioGPUParser::ParallelConfigContext* BioGPUParser::ParallelMapContext::parallelConfig() {
  return getRuleContext<BioGPUParser::ParallelConfigContext>(0);
}

tree::TerminalNode* BioGPUParser::ParallelMapContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}


size_t BioGPUParser::ParallelMapContext::getRuleIndex() const {
  return BioGPUParser::RuleParallelMap;
}


std::any BioGPUParser::ParallelMapContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitParallelMap(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ParallelMapContext* BioGPUParser::parallelMap() {
  ParallelMapContext *_localctx = _tracker.createInstance<ParallelMapContext>(_ctx, getState());
  enterRule(_localctx, 44, BioGPUParser::RuleParallelMap);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(270);
    match(BioGPUParser::PARALLEL_MAP);
    setState(271);
    match(BioGPUParser::LPAREN);
    setState(272);
    expression(0);
    setState(273);
    match(BioGPUParser::COMMA);
    setState(274);
    expression(0);
    setState(275);
    match(BioGPUParser::RPAREN);
    setState(276);
    match(BioGPUParser::LBRACE);
    setState(277);
    parallelConfig();
    setState(278);
    match(BioGPUParser::RBRACE);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ParallelConfigContext ------------------------------------------------------------------

BioGPUParser::ParallelConfigContext::ParallelConfigContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

std::vector<BioGPUParser::ConfigParamContext *> BioGPUParser::ParallelConfigContext::configParam() {
  return getRuleContexts<BioGPUParser::ConfigParamContext>();
}

BioGPUParser::ConfigParamContext* BioGPUParser::ParallelConfigContext::configParam(size_t i) {
  return getRuleContext<BioGPUParser::ConfigParamContext>(i);
}

std::vector<tree::TerminalNode *> BioGPUParser::ParallelConfigContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::ParallelConfigContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::ParallelConfigContext::getRuleIndex() const {
  return BioGPUParser::RuleParallelConfig;
}


std::any BioGPUParser::ParallelConfigContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitParallelConfig(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ParallelConfigContext* BioGPUParser::parallelConfig() {
  ParallelConfigContext *_localctx = _tracker.createInstance<ParallelConfigContext>(_ctx, getState());
  enterRule(_localctx, 46, BioGPUParser::RuleParallelConfig);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(280);
    configParam();
    setState(285);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(281);
      match(BioGPUParser::COMMA);
      setState(282);
      configParam();
      setState(287);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- IfStatementContext ------------------------------------------------------------------

BioGPUParser::IfStatementContext::IfStatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::IfStatementContext::IF() {
  return getToken(BioGPUParser::IF, 0);
}

BioGPUParser::ExpressionContext* BioGPUParser::IfStatementContext::expression() {
  return getRuleContext<BioGPUParser::ExpressionContext>(0);
}

std::vector<tree::TerminalNode *> BioGPUParser::IfStatementContext::LBRACE() {
  return getTokens(BioGPUParser::LBRACE);
}

tree::TerminalNode* BioGPUParser::IfStatementContext::LBRACE(size_t i) {
  return getToken(BioGPUParser::LBRACE, i);
}

std::vector<BioGPUParser::StageBodyContext *> BioGPUParser::IfStatementContext::stageBody() {
  return getRuleContexts<BioGPUParser::StageBodyContext>();
}

BioGPUParser::StageBodyContext* BioGPUParser::IfStatementContext::stageBody(size_t i) {
  return getRuleContext<BioGPUParser::StageBodyContext>(i);
}

std::vector<tree::TerminalNode *> BioGPUParser::IfStatementContext::RBRACE() {
  return getTokens(BioGPUParser::RBRACE);
}

tree::TerminalNode* BioGPUParser::IfStatementContext::RBRACE(size_t i) {
  return getToken(BioGPUParser::RBRACE, i);
}

tree::TerminalNode* BioGPUParser::IfStatementContext::ELSE() {
  return getToken(BioGPUParser::ELSE, 0);
}


size_t BioGPUParser::IfStatementContext::getRuleIndex() const {
  return BioGPUParser::RuleIfStatement;
}


std::any BioGPUParser::IfStatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitIfStatement(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::IfStatementContext* BioGPUParser::ifStatement() {
  IfStatementContext *_localctx = _tracker.createInstance<IfStatementContext>(_ctx, getState());
  enterRule(_localctx, 48, BioGPUParser::RuleIfStatement);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(288);
    match(BioGPUParser::IF);
    setState(289);
    expression(0);
    setState(290);
    match(BioGPUParser::LBRACE);
    setState(291);
    stageBody();
    setState(292);
    match(BioGPUParser::RBRACE);
    setState(298);
    _errHandler->sync(this);

    _la = _input->LA(1);
    if (_la == BioGPUParser::ELSE) {
      setState(293);
      match(BioGPUParser::ELSE);
      setState(294);
      match(BioGPUParser::LBRACE);
      setState(295);
      stageBody();
      setState(296);
      match(BioGPUParser::RBRACE);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ForStatementContext ------------------------------------------------------------------

BioGPUParser::ForStatementContext::ForStatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ForStatementContext::FOR() {
  return getToken(BioGPUParser::FOR, 0);
}

tree::TerminalNode* BioGPUParser::ForStatementContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

tree::TerminalNode* BioGPUParser::ForStatementContext::IN() {
  return getToken(BioGPUParser::IN, 0);
}

BioGPUParser::ExpressionContext* BioGPUParser::ForStatementContext::expression() {
  return getRuleContext<BioGPUParser::ExpressionContext>(0);
}

tree::TerminalNode* BioGPUParser::ForStatementContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

BioGPUParser::StageBodyContext* BioGPUParser::ForStatementContext::stageBody() {
  return getRuleContext<BioGPUParser::StageBodyContext>(0);
}

tree::TerminalNode* BioGPUParser::ForStatementContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}

tree::TerminalNode* BioGPUParser::ForStatementContext::FOREACH() {
  return getToken(BioGPUParser::FOREACH, 0);
}


size_t BioGPUParser::ForStatementContext::getRuleIndex() const {
  return BioGPUParser::RuleForStatement;
}


std::any BioGPUParser::ForStatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitForStatement(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ForStatementContext* BioGPUParser::forStatement() {
  ForStatementContext *_localctx = _tracker.createInstance<ForStatementContext>(_ctx, getState());
  enterRule(_localctx, 50, BioGPUParser::RuleForStatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(316);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::FOR: {
        enterOuterAlt(_localctx, 1);
        setState(300);
        match(BioGPUParser::FOR);
        setState(301);
        match(BioGPUParser::ID);
        setState(302);
        match(BioGPUParser::IN);
        setState(303);
        expression(0);
        setState(304);
        match(BioGPUParser::LBRACE);
        setState(305);
        stageBody();
        setState(306);
        match(BioGPUParser::RBRACE);
        break;
      }

      case BioGPUParser::FOREACH: {
        enterOuterAlt(_localctx, 2);
        setState(308);
        match(BioGPUParser::FOREACH);
        setState(309);
        match(BioGPUParser::ID);
        setState(310);
        match(BioGPUParser::IN);
        setState(311);
        expression(0);
        setState(312);
        match(BioGPUParser::LBRACE);
        setState(313);
        stageBody();
        setState(314);
        match(BioGPUParser::RBRACE);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- EmitStatementContext ------------------------------------------------------------------

BioGPUParser::EmitStatementContext::EmitStatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::EmitStatementContext::EMIT() {
  return getToken(BioGPUParser::EMIT, 0);
}

tree::TerminalNode* BioGPUParser::EmitStatementContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

std::vector<tree::TerminalNode *> BioGPUParser::EmitStatementContext::ID() {
  return getTokens(BioGPUParser::ID);
}

tree::TerminalNode* BioGPUParser::EmitStatementContext::ID(size_t i) {
  return getToken(BioGPUParser::ID, i);
}

std::vector<tree::TerminalNode *> BioGPUParser::EmitStatementContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::EmitStatementContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::EmitStatementContext::getRuleIndex() const {
  return BioGPUParser::RuleEmitStatement;
}


std::any BioGPUParser::EmitStatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitEmitStatement(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::EmitStatementContext* BioGPUParser::emitStatement() {
  EmitStatementContext *_localctx = _tracker.createInstance<EmitStatementContext>(_ctx, getState());
  enterRule(_localctx, 52, BioGPUParser::RuleEmitStatement);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(318);
    match(BioGPUParser::EMIT);
    setState(319);
    match(BioGPUParser::COLON);
    setState(320);
    match(BioGPUParser::ID);
    setState(325);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(321);
      match(BioGPUParser::COMMA);
      setState(322);
      match(BioGPUParser::ID);
      setState(327);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ReportBlockContext ------------------------------------------------------------------

BioGPUParser::ReportBlockContext::ReportBlockContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ReportBlockContext::REPORT() {
  return getToken(BioGPUParser::REPORT, 0);
}

tree::TerminalNode* BioGPUParser::ReportBlockContext::LBRACE() {
  return getToken(BioGPUParser::LBRACE, 0);
}

tree::TerminalNode* BioGPUParser::ReportBlockContext::RBRACE() {
  return getToken(BioGPUParser::RBRACE, 0);
}

std::vector<BioGPUParser::ReportStatementContext *> BioGPUParser::ReportBlockContext::reportStatement() {
  return getRuleContexts<BioGPUParser::ReportStatementContext>();
}

BioGPUParser::ReportStatementContext* BioGPUParser::ReportBlockContext::reportStatement(size_t i) {
  return getRuleContext<BioGPUParser::ReportStatementContext>(i);
}


size_t BioGPUParser::ReportBlockContext::getRuleIndex() const {
  return BioGPUParser::RuleReportBlock;
}


std::any BioGPUParser::ReportBlockContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitReportBlock(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ReportBlockContext* BioGPUParser::reportBlock() {
  ReportBlockContext *_localctx = _tracker.createInstance<ReportBlockContext>(_ctx, getState());
  enterRule(_localctx, 54, BioGPUParser::RuleReportBlock);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(328);
    match(BioGPUParser::REPORT);
    setState(329);
    match(BioGPUParser::LBRACE);
    setState(331); 
    _errHandler->sync(this);
    _la = _input->LA(1);
    do {
      setState(330);
      reportStatement();
      setState(333); 
      _errHandler->sync(this);
      _la = _input->LA(1);
    } while ((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << BioGPUParser::IF)
      | (1ULL << BioGPUParser::FOR)
      | (1ULL << BioGPUParser::FOREACH)
      | (1ULL << BioGPUParser::PRINT)
      | (1ULL << BioGPUParser::ALERT)
      | (1ULL << BioGPUParser::RECOMMENDATION))) != 0));
    setState(335);
    match(BioGPUParser::RBRACE);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ReportStatementContext ------------------------------------------------------------------

BioGPUParser::ReportStatementContext::ReportStatementContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ReportStatementContext::PRINT() {
  return getToken(BioGPUParser::PRINT, 0);
}

tree::TerminalNode* BioGPUParser::ReportStatementContext::STRING() {
  return getToken(BioGPUParser::STRING, 0);
}

tree::TerminalNode* BioGPUParser::ReportStatementContext::ALERT() {
  return getToken(BioGPUParser::ALERT, 0);
}

tree::TerminalNode* BioGPUParser::ReportStatementContext::COLON() {
  return getToken(BioGPUParser::COLON, 0);
}

tree::TerminalNode* BioGPUParser::ReportStatementContext::RECOMMENDATION() {
  return getToken(BioGPUParser::RECOMMENDATION, 0);
}

BioGPUParser::IfStatementContext* BioGPUParser::ReportStatementContext::ifStatement() {
  return getRuleContext<BioGPUParser::IfStatementContext>(0);
}

BioGPUParser::ForStatementContext* BioGPUParser::ReportStatementContext::forStatement() {
  return getRuleContext<BioGPUParser::ForStatementContext>(0);
}


size_t BioGPUParser::ReportStatementContext::getRuleIndex() const {
  return BioGPUParser::RuleReportStatement;
}


std::any BioGPUParser::ReportStatementContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitReportStatement(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ReportStatementContext* BioGPUParser::reportStatement() {
  ReportStatementContext *_localctx = _tracker.createInstance<ReportStatementContext>(_ctx, getState());
  enterRule(_localctx, 56, BioGPUParser::RuleReportStatement);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(347);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::PRINT: {
        enterOuterAlt(_localctx, 1);
        setState(337);
        match(BioGPUParser::PRINT);
        setState(338);
        match(BioGPUParser::STRING);
        break;
      }

      case BioGPUParser::ALERT: {
        enterOuterAlt(_localctx, 2);
        setState(339);
        match(BioGPUParser::ALERT);
        setState(340);
        match(BioGPUParser::COLON);
        setState(341);
        match(BioGPUParser::STRING);
        break;
      }

      case BioGPUParser::RECOMMENDATION: {
        enterOuterAlt(_localctx, 3);
        setState(342);
        match(BioGPUParser::RECOMMENDATION);
        setState(343);
        match(BioGPUParser::COLON);
        setState(344);
        match(BioGPUParser::STRING);
        break;
      }

      case BioGPUParser::IF: {
        enterOuterAlt(_localctx, 4);
        setState(345);
        ifStatement();
        break;
      }

      case BioGPUParser::FOR:
      case BioGPUParser::FOREACH: {
        enterOuterAlt(_localctx, 5);
        setState(346);
        forStatement();
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ExpressionContext ------------------------------------------------------------------

BioGPUParser::ExpressionContext::ExpressionContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ExpressionContext::ID() {
  return getToken(BioGPUParser::ID, 0);
}

BioGPUParser::LiteralContext* BioGPUParser::ExpressionContext::literal() {
  return getRuleContext<BioGPUParser::LiteralContext>(0);
}

BioGPUParser::FunctionCallContext* BioGPUParser::ExpressionContext::functionCall() {
  return getRuleContext<BioGPUParser::FunctionCallContext>(0);
}

tree::TerminalNode* BioGPUParser::ExpressionContext::LPAREN() {
  return getToken(BioGPUParser::LPAREN, 0);
}

std::vector<BioGPUParser::ExpressionContext *> BioGPUParser::ExpressionContext::expression() {
  return getRuleContexts<BioGPUParser::ExpressionContext>();
}

BioGPUParser::ExpressionContext* BioGPUParser::ExpressionContext::expression(size_t i) {
  return getRuleContext<BioGPUParser::ExpressionContext>(i);
}

tree::TerminalNode* BioGPUParser::ExpressionContext::RPAREN() {
  return getToken(BioGPUParser::RPAREN, 0);
}

tree::TerminalNode* BioGPUParser::ExpressionContext::DOT() {
  return getToken(BioGPUParser::DOT, 0);
}

tree::TerminalNode* BioGPUParser::ExpressionContext::LBRACK() {
  return getToken(BioGPUParser::LBRACK, 0);
}

tree::TerminalNode* BioGPUParser::ExpressionContext::RBRACK() {
  return getToken(BioGPUParser::RBRACK, 0);
}


size_t BioGPUParser::ExpressionContext::getRuleIndex() const {
  return BioGPUParser::RuleExpression;
}


std::any BioGPUParser::ExpressionContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitExpression(this);
  else
    return visitor->visitChildren(this);
}


BioGPUParser::ExpressionContext* BioGPUParser::expression() {
   return expression(0);
}

BioGPUParser::ExpressionContext* BioGPUParser::expression(int precedence) {
  ParserRuleContext *parentContext = _ctx;
  size_t parentState = getState();
  BioGPUParser::ExpressionContext *_localctx = _tracker.createInstance<ExpressionContext>(_ctx, parentState);
  BioGPUParser::ExpressionContext *previousContext = _localctx;
  (void)previousContext; // Silence compiler, in case the context is not used by generated code.
  size_t startState = 58;
  enterRecursionRule(_localctx, 58, BioGPUParser::RuleExpression, precedence);

    

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    unrollRecursionContexts(parentContext);
  });
  try {
    size_t alt;
    enterOuterAlt(_localctx, 1);
    setState(357);
    _errHandler->sync(this);
    switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 32, _ctx)) {
    case 1: {
      setState(350);
      match(BioGPUParser::ID);
      break;
    }

    case 2: {
      setState(351);
      literal();
      break;
    }

    case 3: {
      setState(352);
      functionCall();
      break;
    }

    case 4: {
      setState(353);
      match(BioGPUParser::LPAREN);
      setState(354);
      expression(0);
      setState(355);
      match(BioGPUParser::RPAREN);
      break;
    }

    default:
      break;
    }
    _ctx->stop = _input->LT(-1);
    setState(369);
    _errHandler->sync(this);
    alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 34, _ctx);
    while (alt != 2 && alt != atn::ATN::INVALID_ALT_NUMBER) {
      if (alt == 1) {
        if (!_parseListeners.empty())
          triggerExitRuleEvent();
        previousContext = _localctx;
        setState(367);
        _errHandler->sync(this);
        switch (getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 33, _ctx)) {
        case 1: {
          _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
          pushNewRecursionContext(_localctx, startState, RuleExpression);
          setState(359);

          if (!(precpred(_ctx, 3))) throw FailedPredicateException(this, "precpred(_ctx, 3)");
          setState(360);
          match(BioGPUParser::DOT);
          setState(361);
          match(BioGPUParser::ID);
          break;
        }

        case 2: {
          _localctx = _tracker.createInstance<ExpressionContext>(parentContext, parentState);
          pushNewRecursionContext(_localctx, startState, RuleExpression);
          setState(362);

          if (!(precpred(_ctx, 2))) throw FailedPredicateException(this, "precpred(_ctx, 2)");
          setState(363);
          match(BioGPUParser::LBRACK);
          setState(364);
          expression(0);
          setState(365);
          match(BioGPUParser::RBRACK);
          break;
        }

        default:
          break;
        } 
      }
      setState(371);
      _errHandler->sync(this);
      alt = getInterpreter<atn::ParserATNSimulator>()->adaptivePredict(_input, 34, _ctx);
    }
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }
  return _localctx;
}

//----------------- ArrayContext ------------------------------------------------------------------

BioGPUParser::ArrayContext::ArrayContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ArrayContext::LBRACK() {
  return getToken(BioGPUParser::LBRACK, 0);
}

tree::TerminalNode* BioGPUParser::ArrayContext::RBRACK() {
  return getToken(BioGPUParser::RBRACK, 0);
}

std::vector<BioGPUParser::LiteralContext *> BioGPUParser::ArrayContext::literal() {
  return getRuleContexts<BioGPUParser::LiteralContext>();
}

BioGPUParser::LiteralContext* BioGPUParser::ArrayContext::literal(size_t i) {
  return getRuleContext<BioGPUParser::LiteralContext>(i);
}

std::vector<tree::TerminalNode *> BioGPUParser::ArrayContext::ID() {
  return getTokens(BioGPUParser::ID);
}

tree::TerminalNode* BioGPUParser::ArrayContext::ID(size_t i) {
  return getToken(BioGPUParser::ID, i);
}

std::vector<tree::TerminalNode *> BioGPUParser::ArrayContext::COMMA() {
  return getTokens(BioGPUParser::COMMA);
}

tree::TerminalNode* BioGPUParser::ArrayContext::COMMA(size_t i) {
  return getToken(BioGPUParser::COMMA, i);
}


size_t BioGPUParser::ArrayContext::getRuleIndex() const {
  return BioGPUParser::RuleArray;
}


std::any BioGPUParser::ArrayContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitArray(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ArrayContext* BioGPUParser::array() {
  ArrayContext *_localctx = _tracker.createInstance<ArrayContext>(_ctx, getState());
  enterRule(_localctx, 60, BioGPUParser::RuleArray);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(372);
    match(BioGPUParser::LBRACK);
    setState(375);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::BOOLEAN:
      case BioGPUParser::NUMBER:
      case BioGPUParser::STRING:
      case BioGPUParser::DNA_SEQUENCE: {
        setState(373);
        literal();
        break;
      }

      case BioGPUParser::ID: {
        setState(374);
        match(BioGPUParser::ID);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
    setState(384);
    _errHandler->sync(this);
    _la = _input->LA(1);
    while (_la == BioGPUParser::COMMA) {
      setState(377);
      match(BioGPUParser::COMMA);
      setState(380);
      _errHandler->sync(this);
      switch (_input->LA(1)) {
        case BioGPUParser::BOOLEAN:
        case BioGPUParser::NUMBER:
        case BioGPUParser::STRING:
        case BioGPUParser::DNA_SEQUENCE: {
          setState(378);
          literal();
          break;
        }

        case BioGPUParser::ID: {
          setState(379);
          match(BioGPUParser::ID);
          break;
        }

      default:
        throw NoViableAltException(this);
      }
      setState(386);
      _errHandler->sync(this);
      _la = _input->LA(1);
    }
    setState(387);
    match(BioGPUParser::RBRACK);
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- LiteralContext ------------------------------------------------------------------

BioGPUParser::LiteralContext::LiteralContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::LiteralContext::NUMBER() {
  return getToken(BioGPUParser::NUMBER, 0);
}

tree::TerminalNode* BioGPUParser::LiteralContext::STRING() {
  return getToken(BioGPUParser::STRING, 0);
}

tree::TerminalNode* BioGPUParser::LiteralContext::BOOLEAN() {
  return getToken(BioGPUParser::BOOLEAN, 0);
}

tree::TerminalNode* BioGPUParser::LiteralContext::DNA_SEQUENCE() {
  return getToken(BioGPUParser::DNA_SEQUENCE, 0);
}


size_t BioGPUParser::LiteralContext::getRuleIndex() const {
  return BioGPUParser::RuleLiteral;
}


std::any BioGPUParser::LiteralContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitLiteral(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::LiteralContext* BioGPUParser::literal() {
  LiteralContext *_localctx = _tracker.createInstance<LiteralContext>(_ctx, getState());
  enterRule(_localctx, 62, BioGPUParser::RuleLiteral);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(389);
    _la = _input->LA(1);
    if (!((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << BioGPUParser::BOOLEAN)
      | (1ULL << BioGPUParser::NUMBER)
      | (1ULL << BioGPUParser::STRING)
      | (1ULL << BioGPUParser::DNA_SEQUENCE))) != 0))) {
    _errHandler->recoverInline(this);
    }
    else {
      _errHandler->reportMatch(this);
      consume();
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- DataTypeContext ------------------------------------------------------------------

BioGPUParser::DataTypeContext::DataTypeContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::DataTypeContext::FASTQ_FILE() {
  return getToken(BioGPUParser::FASTQ_FILE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::FASTA_FILE() {
  return getToken(BioGPUParser::FASTA_FILE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::CSV_FILE() {
  return getToken(BioGPUParser::CSV_FILE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::JSON() {
  return getToken(BioGPUParser::JSON, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::PDF() {
  return getToken(BioGPUParser::PDF, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::STRING_TYPE() {
  return getToken(BioGPUParser::STRING_TYPE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::INT_TYPE() {
  return getToken(BioGPUParser::INT_TYPE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::FLOAT_TYPE() {
  return getToken(BioGPUParser::FLOAT_TYPE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::BOOLEAN_TYPE() {
  return getToken(BioGPUParser::BOOLEAN_TYPE, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::ARRAY() {
  return getToken(BioGPUParser::ARRAY, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::LT() {
  return getToken(BioGPUParser::LT, 0);
}

std::vector<BioGPUParser::DataTypeContext *> BioGPUParser::DataTypeContext::dataType() {
  return getRuleContexts<BioGPUParser::DataTypeContext>();
}

BioGPUParser::DataTypeContext* BioGPUParser::DataTypeContext::dataType(size_t i) {
  return getRuleContext<BioGPUParser::DataTypeContext>(i);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::GT() {
  return getToken(BioGPUParser::GT, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::MAP() {
  return getToken(BioGPUParser::MAP, 0);
}

tree::TerminalNode* BioGPUParser::DataTypeContext::COMMA() {
  return getToken(BioGPUParser::COMMA, 0);
}


size_t BioGPUParser::DataTypeContext::getRuleIndex() const {
  return BioGPUParser::RuleDataType;
}


std::any BioGPUParser::DataTypeContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitDataType(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::DataTypeContext* BioGPUParser::dataType() {
  DataTypeContext *_localctx = _tracker.createInstance<DataTypeContext>(_ctx, getState());
  enterRule(_localctx, 64, BioGPUParser::RuleDataType);

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    setState(412);
    _errHandler->sync(this);
    switch (_input->LA(1)) {
      case BioGPUParser::FASTQ_FILE: {
        enterOuterAlt(_localctx, 1);
        setState(391);
        match(BioGPUParser::FASTQ_FILE);
        break;
      }

      case BioGPUParser::FASTA_FILE: {
        enterOuterAlt(_localctx, 2);
        setState(392);
        match(BioGPUParser::FASTA_FILE);
        break;
      }

      case BioGPUParser::CSV_FILE: {
        enterOuterAlt(_localctx, 3);
        setState(393);
        match(BioGPUParser::CSV_FILE);
        break;
      }

      case BioGPUParser::JSON: {
        enterOuterAlt(_localctx, 4);
        setState(394);
        match(BioGPUParser::JSON);
        break;
      }

      case BioGPUParser::PDF: {
        enterOuterAlt(_localctx, 5);
        setState(395);
        match(BioGPUParser::PDF);
        break;
      }

      case BioGPUParser::STRING_TYPE: {
        enterOuterAlt(_localctx, 6);
        setState(396);
        match(BioGPUParser::STRING_TYPE);
        break;
      }

      case BioGPUParser::INT_TYPE: {
        enterOuterAlt(_localctx, 7);
        setState(397);
        match(BioGPUParser::INT_TYPE);
        break;
      }

      case BioGPUParser::FLOAT_TYPE: {
        enterOuterAlt(_localctx, 8);
        setState(398);
        match(BioGPUParser::FLOAT_TYPE);
        break;
      }

      case BioGPUParser::BOOLEAN_TYPE: {
        enterOuterAlt(_localctx, 9);
        setState(399);
        match(BioGPUParser::BOOLEAN_TYPE);
        break;
      }

      case BioGPUParser::ARRAY: {
        enterOuterAlt(_localctx, 10);
        setState(400);
        match(BioGPUParser::ARRAY);
        setState(401);
        match(BioGPUParser::LT);
        setState(402);
        dataType();
        setState(403);
        match(BioGPUParser::GT);
        break;
      }

      case BioGPUParser::MAP: {
        enterOuterAlt(_localctx, 11);
        setState(405);
        match(BioGPUParser::MAP);
        setState(406);
        match(BioGPUParser::LT);
        setState(407);
        dataType();
        setState(408);
        match(BioGPUParser::COMMA);
        setState(409);
        dataType();
        setState(410);
        match(BioGPUParser::GT);
        break;
      }

    default:
      throw NoViableAltException(this);
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

//----------------- ReferenceTypeContext ------------------------------------------------------------------

BioGPUParser::ReferenceTypeContext::ReferenceTypeContext(ParserRuleContext *parent, size_t invokingState)
  : ParserRuleContext(parent, invokingState) {
}

tree::TerminalNode* BioGPUParser::ReferenceTypeContext::GENOME_DATABASE() {
  return getToken(BioGPUParser::GENOME_DATABASE, 0);
}

tree::TerminalNode* BioGPUParser::ReferenceTypeContext::MUTATION_DATABASE() {
  return getToken(BioGPUParser::MUTATION_DATABASE, 0);
}

tree::TerminalNode* BioGPUParser::ReferenceTypeContext::GENE_DATABASE() {
  return getToken(BioGPUParser::GENE_DATABASE, 0);
}


size_t BioGPUParser::ReferenceTypeContext::getRuleIndex() const {
  return BioGPUParser::RuleReferenceType;
}


std::any BioGPUParser::ReferenceTypeContext::accept(tree::ParseTreeVisitor *visitor) {
  if (auto parserVisitor = dynamic_cast<BioGPUVisitor*>(visitor))
    return parserVisitor->visitReferenceType(this);
  else
    return visitor->visitChildren(this);
}

BioGPUParser::ReferenceTypeContext* BioGPUParser::referenceType() {
  ReferenceTypeContext *_localctx = _tracker.createInstance<ReferenceTypeContext>(_ctx, getState());
  enterRule(_localctx, 66, BioGPUParser::RuleReferenceType);
  size_t _la = 0;

#if __cplusplus > 201703L
  auto onExit = finally([=, this] {
#else
  auto onExit = finally([=] {
#endif
    exitRule();
  });
  try {
    enterOuterAlt(_localctx, 1);
    setState(414);
    _la = _input->LA(1);
    if (!((((_la & ~ 0x3fULL) == 0) &&
      ((1ULL << _la) & ((1ULL << BioGPUParser::GENOME_DATABASE)
      | (1ULL << BioGPUParser::MUTATION_DATABASE)
      | (1ULL << BioGPUParser::GENE_DATABASE))) != 0))) {
    _errHandler->recoverInline(this);
    }
    else {
      _errHandler->reportMatch(this);
      consume();
    }
   
  }
  catch (RecognitionException &e) {
    _errHandler->reportError(this, e);
    _localctx->exception = std::current_exception();
    _errHandler->recover(this, _localctx->exception);
  }

  return _localctx;
}

bool BioGPUParser::sempred(RuleContext *context, size_t ruleIndex, size_t predicateIndex) {
  switch (ruleIndex) {
    case 29: return expressionSempred(antlrcpp::downCast<ExpressionContext *>(context), predicateIndex);

  default:
    break;
  }
  return true;
}

bool BioGPUParser::expressionSempred(ExpressionContext *_localctx, size_t predicateIndex) {
  switch (predicateIndex) {
    case 0: return precpred(_ctx, 3);
    case 1: return precpred(_ctx, 2);

  default:
    break;
  }
  return true;
}

void BioGPUParser::initialize() {
  std::call_once(biogpuParserOnceFlag, biogpuParserInitialize);
}
