#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

int main(int argc, char* argv[]) {
    std::cout << "Simple heap corruption test" << std::endl;
    
    // Test basic string operations that might cause heap corruption
    std::vector<std::string> test_strings;
    
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            std::string arg(argv[i]);
            test_strings.push_back(arg);
            std::cout << "arg[" << i << "] = " << arg << std::endl;
        }
    }
    
    // Test some allocations
    for (int i = 0; i < 1000; i++) {
        std::string test = "test_string_" + std::to_string(i);
        test_strings.push_back(test);
    }
    
    std::cout << "Created " << test_strings.size() << " strings without crash" << std::endl;
    
    if (argc >= 3 && std::string(argv[1]) == "build") {
        std::cout << "Build command simulation - would process: " << argv[2] << std::endl;
    }
    
    return 0;
}
