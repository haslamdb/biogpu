#include <iostream>
#include <fstream>
#include "antlr4-runtime.h"
#include "BioGPULexer.h"
#include "BioGPUParser.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.biogpu>" << std::endl;
        return 1;
    }
    
    std::ifstream stream(argv[1]);
    antlr4::ANTLRInputStream input(stream);
    BioGPULexer lexer(&input);
    antlr4::CommonTokenStream tokens(&lexer);
    BioGPUParser parser(&tokens);
    
    auto tree = parser.program();
    std::cout << "Parse tree: " << tree->toStringTree(&parser) << std::endl;
    
    return 0;
}
