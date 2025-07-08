#include <iostream>
#include <string>
#include <vector>
#include <sstream>

std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

std::string extractGeneNameFromHeader(const std::string& header) {
    // First check for AMRProt.fa format: 0|AAA16360.1|1|1|gene_name|...
    if (header.find('|') != std::string::npos) {
        std::vector<std::string> parts = splitString(header, '|');
        std::cout << "Header parts (" << parts.size() << "):" << std::endl;
        for (size_t i = 0; i < parts.size(); i++) {
            std::cout << "  [" << i << "] = '" << parts[i] << "'" << std::endl;
        }
        
        if (parts.size() >= 5 && parts[0].length() <= 2) {
            // Return the gene name from parts[4]
            std::cout << "Returning gene name from position 4: '" << parts[4] << "'" << std::endl;
            return parts[4];
        }
    }
    
    return "unknown";
}

int main() {
    std::string test_headers[] = {
        ">0|AAA16360.1|1|1|stxA2b|stxA2b||1|stxA2b|STX2|Shiga_toxin_Stx2b_subunit_A",
        ">0|AAA16361.1|1|1|stxB2b|stxB2b||1|stxB2b|STX2|Shiga_toxin_Stx2b_subunit_B",
        "0|AAA16360.1|1|1|stxA2b|stxA2b||1|stxA2b|STX2|Shiga_toxin_Stx2b_subunit_A"
    };
    
    for (const auto& header : test_headers) {
        std::cout << "\nTesting header: " << header << std::endl;
        std::string gene_name = extractGeneNameFromHeader(header);
        std::cout << "Extracted gene name: '" << gene_name << "'" << std::endl;
    }
    
    return 0;
}