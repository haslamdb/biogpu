// tests/test_translated_search.cu
// Test program for 6-frame translation and protein search

#include <iostream>
#include <cuda_runtime.h>

// Declare the translated search interface
extern "C" {
    void* create_translated_search_engine(int batch_size);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
}

int main(int argc, char** argv) {
    std::cout << "=== Testing Translated Search Engine ===" << std::endl;
    
    // Create engine
    void* engine = create_translated_search_engine(1000);
    if (!engine) {
        std::cerr << "Failed to create translated search engine" << std::endl;
        return 1;
    }
    
    std::cout << "Translated search engine created successfully" << std::endl;
    
    // Test loading protein database (if path provided)
    if (argc > 1) {
        std::cout << "Loading protein database from: " << argv[1] << std::endl;
        int result = load_protein_database(engine, argv[1]);
        if (result == 0) {
            std::cout << "Protein database loaded successfully" << std::endl;
        } else {
            std::cout << "Failed to load protein database" << std::endl;
        }
    }
    
    // Cleanup
    destroy_translated_search_engine(engine);
    std::cout << "Translated search engine destroyed" << std::endl;
    
    std::cout << "=== Test Complete ===" << std::endl;
    return 0;
}