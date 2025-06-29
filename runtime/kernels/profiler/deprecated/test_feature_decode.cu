// test_feature_decode.cu
// Simple program to decode and display feature flags from minimizers

#include "gpu_kraken_types.h"
#include <iostream>
#include <iomanip>
#include <bitset>

void decode_feature_flags(uint32_t feature_flags) {
    std::cout << "Feature flags: 0x" << std::hex << feature_flags 
              << " (binary: " << std::bitset<32>(feature_flags) << ")" << std::dec << std::endl;
    
    // Decode bits 0-2: GC content
    uint32_t gc_category = MinimizerFlags::get_gc_content_category(feature_flags);
    std::cout << "  GC Content Category: " << gc_category << " (";
    switch(gc_category) {
        case 0: std::cout << "0-12.5%"; break;
        case 1: std::cout << "12.5-25%"; break;
        case 2: std::cout << "25-37.5%"; break;
        case 3: std::cout << "37.5-50%"; break;
        case 4: std::cout << "50-62.5%"; break;
        case 5: std::cout << "62.5-75%"; break;
        case 6: std::cout << "75-87.5%"; break;
        case 7: std::cout << "87.5-100%"; break;
    }
    std::cout << ")" << std::endl;
    
    // Decode bits 3-5: Sequence complexity
    uint32_t complexity = MinimizerFlags::get_complexity_score(feature_flags);
    std::cout << "  Sequence Complexity: " << complexity << " (";
    if (complexity < 3) std::cout << "low";
    else if (complexity < 6) std::cout << "medium";
    else std::cout << "high";
    std::cout << ")" << std::endl;
    
    // Decode bit 6: Position bias
    bool position_bias = MinimizerFlags::has_position_bias(feature_flags);
    std::cout << "  Position Bias: " << (position_bias ? "clustered" : "uniform") << std::endl;
    
    // Decode bit 7: Contamination risk
    bool contam_risk = MinimizerFlags::has_contamination_risk(feature_flags);
    std::cout << "  Contamination Risk: " << (contam_risk ? "yes" : "no") << std::endl;
    
    // Decode bits 13-15: Co-occurrence score (new feature)
    uint8_t cooccur = MinimizerFlags::get_cooccurrence_score(feature_flags);
    std::cout << "  Co-occurrence Score: " << (int)cooccur << std::endl;
    
    // Decode bits 16-18: Taxonomic ambiguity (new feature)
    uint8_t tax_ambig = MinimizerFlags::get_taxonomic_ambiguity(feature_flags);
    std::cout << "  Taxonomic Ambiguity: " << (int)tax_ambig << std::endl;
    
    // Decode bits 19-21: Context complexity (new feature)
    uint8_t ctx_complex = MinimizerFlags::get_context_complexity(feature_flags);
    std::cout << "  Context Complexity: " << (int)ctx_complex << std::endl;
    
    // Decode bits 22-23: Read coherence (new feature)
    uint8_t read_coher = MinimizerFlags::get_read_coherence(feature_flags);
    std::cout << "  Read Coherence: " << (int)read_coher << std::endl;
    
    // Decode bits 28-30: Contamination type (moved from 13-15)
    uint8_t contam_type = (feature_flags >> 28) & 0x7;
    if (contam_type > 0) {
        std::cout << "  Contamination Type: " << (int)contam_type << " (";
        switch(contam_type) {
            case 1: std::cout << "human"; break;
            case 2: std::cout << "adapter"; break;
            case 3: std::cout << "barcode"; break;
            case 4: std::cout << "environmental"; break;
            case 5: std::cout << "phix"; break;
            case 6: std::cout << "synthetic"; break;
            case 7: std::cout << "other"; break;
        }
        std::cout << ")" << std::endl;
    }
}

int main() {
    std::cout << "=== MINIMIZER FEATURE FLAG DECODER ===\n" << std::endl;
    
    // Test with the feature flags from the debug output
    uint32_t test_flags[] = {61, 60, 54, 53, 52};
    
    std::cout << "Decoding feature flags from test output:\n" << std::endl;
    
    for (int i = 0; i < 5; i++) {
        std::cout << "Minimizer " << i << " feature flags:" << std::endl;
        decode_feature_flags(test_flags[i]);
        std::cout << std::endl;
    }
    
    // Test with some synthetic examples
    std::cout << "\nTesting with synthetic examples:\n" << std::endl;
    
    // Example 1: High GC, high complexity, clustered
    uint32_t example1 = 0;
    example1 = MinimizerFlags::set_gc_content_category(example1, 6);  // 75-87.5% GC
    example1 = MinimizerFlags::set_complexity_score(example1, 7);     // High complexity
    example1 = MinimizerFlags::set_position_bias(example1, true);     // Clustered
    example1 = MinimizerFlags::set_cooccurrence_score(example1, 5);   // High co-occurrence
    std::cout << "Example 1 - High GC, high complexity, clustered:" << std::endl;
    decode_feature_flags(example1);
    
    std::cout << "\n";
    
    // Example 2: Low GC, low complexity, contamination risk
    uint32_t example2 = 0;
    example2 = MinimizerFlags::set_gc_content_category(example2, 1);  // 12.5-25% GC
    example2 = MinimizerFlags::set_complexity_score(example2, 2);     // Low complexity
    example2 = MinimizerFlags::set_contamination_risk(example2, true); // Has contamination risk
    example2 |= (2 << 28);    // Adapter contamination (bits 28-30)
    std::cout << "Example 2 - Low GC, low complexity, contamination:" << std::endl;
    decode_feature_flags(example2);
    
    return 0;
}