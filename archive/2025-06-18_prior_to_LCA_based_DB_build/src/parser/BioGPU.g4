grammar BioGPU;

// Parser rules
program
    : (importStatement)* pipeline* EOF
    ;

importStatement
    : 'import' modulePath ('as' ID)? ';'
    ;

modulePath
    : ID ('.' ID)*
    ;

pipeline
    : 'pipeline' ID '{' pipelineBody '}'
    ;

pipelineBody
    : inputDecl? outputDecl? referencesDecl? stage+
    ;

inputDecl
    : 'input' ':' (inputParam | '{' inputParam (',' inputParam)* '}')
    ;

inputParam
    : ID ':' dataType ('optional')?
    ;

outputDecl
    : 'output' ':' (outputParam | '{' outputParam (',' outputParam)* '}')
    ;

outputParam
    : ID ':' dataType
    ;

referencesDecl
    : 'references' ':' '{' referenceParam (',' referenceParam)* '}'
    ;

referenceParam
    : ID ':' referenceType '(' STRING ')'
    ;

stage
    : decorator* 'stage' ID '{' stageBody '}'
    ;

decorator
    : '@' ID ('(' decoratorArgs ')')?
    ;

decoratorArgs
    : decoratorArg (',' decoratorArg)*
    ;

decoratorArg
    : ID ':' (literal | ID)
    | literal
    ;

stageBody
    : statement+
    ;

statement
    : assignment
    | functionCall
    | parallelMap
    | ifStatement
    | forStatement
    | emitStatement
    | reportBlock
    ;

assignment
    : ID '=' expression
    ;

functionCall
    : ID '(' argumentList? ')' ('{' configBlock '}')?
    ;

argumentList
    : expression (',' expression)*
    ;

configBlock
    : configParam (',' configParam)*
    ;

configParam
    : ID ':' (literal | expression | array)
    ;

parallelMap
    : 'parallel_map' '(' expression ',' expression ')' '{' parallelConfig '}'
    ;

parallelConfig
    : configParam (',' configParam)*
    ;

ifStatement
    : 'if' expression '{' stageBody '}' ('else' '{' stageBody '}')?
    ;

forStatement
    : 'for' ID 'in' expression '{' stageBody '}'
    | 'foreach' ID 'in' expression '{' stageBody '}'
    ;

emitStatement
    : 'emit' ':' ID (',' ID)*
    ;

reportBlock
    : 'report' '{' reportStatement+ '}'
    ;

reportStatement
    : 'print' STRING
    | 'alert' ':' STRING
    | 'recommendation' ':' STRING
    | ifStatement
    | forStatement
    ;

expression
    : ID
    | literal
    | functionCall
    | expression '.' ID
    | expression '[' expression ']'
    | '(' expression ')'
    ;

array
    : '[' (literal | ID) (',' (literal | ID))* ']'
    ;

literal
    : NUMBER
    | STRING
    | BOOLEAN
    | DNA_SEQUENCE
    ;

dataType
    : 'fastq_file'
    | 'fasta_file'
    | 'csv_file'
    | 'json'
    | 'pdf'
    | 'string'
    | 'int'
    | 'float'
    | 'boolean'
    | 'array' '<' dataType '>'
    | 'map' '<' dataType ',' dataType '>'
    ;

referenceType
    : 'genome_database'
    | 'mutation_database'
    | 'gene_database'
    ;

// Lexer rules
// Keywords
PIPELINE : 'pipeline';
STAGE : 'stage';
INPUT : 'input';
OUTPUT : 'output';
REFERENCES : 'references';
IMPORT : 'import';
AS : 'as';
IF : 'if';
ELSE : 'else';
FOR : 'for';
FOREACH : 'foreach';
IN : 'in';
EMIT : 'emit';
REPORT : 'report';
PRINT : 'print';
ALERT : 'alert';
RECOMMENDATION : 'recommendation';
PARALLEL_MAP : 'parallel_map';
OPTIONAL : 'optional';

// Data types
FASTQ_FILE : 'fastq_file';
FASTA_FILE : 'fasta_file';
CSV_FILE : 'csv_file';
JSON : 'json';
PDF : 'pdf';
STRING_TYPE : 'string';
INT_TYPE : 'int';
FLOAT_TYPE : 'float';
BOOLEAN_TYPE : 'boolean';
ARRAY : 'array';
MAP : 'map';

// Reference types
GENOME_DATABASE : 'genome_database';
MUTATION_DATABASE : 'mutation_database';
GENE_DATABASE : 'gene_database';

// Decorators
GPU_KERNEL : '@gpu_kernel';
PARALLEL : '@parallel';

// Operators
ASSIGN : '=';
DOT : '.';
COMMA : ',';
SEMI : ';';
COLON : ':';
LPAREN : '(';
RPAREN : ')';
LBRACE : '{';
RBRACE : '}';
LBRACK : '[';
RBRACK : ']';
LT : '<';
GT : '>';

// Literals
BOOLEAN
    : 'true'
    | 'false'
    ;

NUMBER
    : '-'? INT ('.' INT)? EXP?
    ;

fragment INT
    : '0' | [1-9] [0-9]*
    ;

fragment EXP
    : [Ee] [+\-]? INT
    ;

STRING
    : '"' ( ESC | ~[\\"] )* '"'
    ;

fragment ESC
    : '\\' [abfnrtv"'\\]
    | UNICODE
    ;

fragment UNICODE
    : '\\' 'u' HEX HEX HEX HEX
    ;

fragment HEX
    : [0-9a-fA-F]
    ;

DNA_SEQUENCE
    : 'DNA:' [ACGTN]+
    ;

ID
    : [a-zA-Z_] [a-zA-Z_0-9]*
    ;

// Comments and whitespace
COMMENT
    : '#' ~[\r\n]* -> skip
    ;

BLOCK_COMMENT
    : '/*' .*? '*/' -> skip
    ;

WS
    : [ \t\r\n]+ -> skip
    ;