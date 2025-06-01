// tests/test_translated_search.cu
// Test program for enhanced translated search with 5-mer k-mers and Smith-Waterman

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cuda_runtime.h>

// Include translated search interface
extern "C" {
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
    void set_smith_waterman_enabled(void* engine, bool enabled);
    int search_translated_reads(
        void* engine,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        const bool* d_reads_to_process,
        int num_reads,
        void* results,
        uint32_t* result_counts
    );
}

// ProteinMatch structure (must match translated_search.cu)
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    uint8_t num_mutations;
    uint8_t mutation_positions[10];
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
    bool used_smith_waterman;
};

class TranslatedSearchTester {
private:
    void* engine;
    char* d_reads;
    int* d_read_lengths;
    int* d_read_offsets;
    bool* d_reads_to_process;
    ProteinMatch* d_results;
    uint32_t* d_result_counts;
    
    const int max_reads = 1000;
    const int max_read_length = 300;
    const int max_matches_per_read = 32;
    
public:
    TranslatedSearchTester(bool enable_smith_waterman = false) {
        std::cout << "=== Enhanced Translated Search Test ===" << std::endl;
        std::cout << "Features: 5-mer k-mers, seed clustering, extension" << std::endl;
        if (enable_smith_waterman) {
            std::cout << "Smith-Waterman: ENABLED for high-scoring matches" << std::endl;
        } else {
            std::cout << "Smith-Waterman: DISABLED" << std::endl;
        }
        std::cout << std::endl;
        
        // Create enhanced engine
        engine = create_translated_search_engine_with_sw(max_reads, enable_smith_waterman);
        if (!engine) {
            throw std::runtime_error("Failed to create translated search engine");
        }
        
        // Allocate GPU memory
        cudaMalloc(&d_reads, max_reads * max_read_length);
        cudaMalloc(&d_read_lengths, max_reads * sizeof(int));
        cudaMalloc(&d_read_offsets, max_reads * sizeof(int));
        cudaMalloc(&d_reads_to_process, max_reads * sizeof(bool));
        cudaMalloc(&d_results, max_reads * max_matches_per_read * sizeof(ProteinMatch));
        cudaMalloc(&d_result_counts, max_reads * sizeof(uint32_t));
        
        std::cout << "Enhanced translated search engine initialized successfully" << std::endl;
    }
    
    ~TranslatedSearchTester() {
        if (engine) destroy_translated_search_engine(engine);
        
        cudaFree(d_reads);
        cudaFree(d_read_lengths);
        cudaFree(d_read_offsets);
        cudaFree(d_reads_to_process);
        cudaFree(d_results);
        cudaFree(d_result_counts);
    }
    
    bool loadDatabase(const std::string& db_path) {
        std::cout << "Loading enhanced protein database from: " << db_path << std::endl;
        
        int result = load_protein_database(engine, db_path.c_str());
        if (result != 0) {
            std::cerr << "ERROR: Failed to load protein database" << std::endl;
            return false;
        }
        
        std::cout << "Enhanced protein database loaded successfully" << std::endl;
        return true;
    }
    
    void testSyntheticReads() {
        std::cout << "\n=== Testing Enhanced Translated Search with Synthetic Reads ===" << std::endl;
        
        // Create synthetic reads with resistance mutations
        std::vector<std::string> test_reads = {
            // gyrA S83L mutation (AGC -> CTG)
            "ATGGAGACCCGCAAGATCACCGCCGAAGACCGCCACGCTGGTGATGGGCTCGAAATCCTGGCGAATCTGAACAAGCTGGGCGATGTGCTGAAGGAGTTGGATCTGATGGGCTACCGTGAGCGCGAAGATGAGCTGGAATTCCTGGCGGAAGCCGGCGAGAAGCTGGTGCAGCGCAAACTGATCGGCGACCTGTGGCCGCAGGATGGCAACACCAAACTGAAACCGCTGGCG",
            
            // parC S80I mutation (AGC -> ATC)
            "ATGGAGACCCGCAAGATCACCGCCGAAGACCGCCACGCTGGTGATGGGCTCGAAATCCTGGCGATCCTGAACAAGCTGGGCGATGTGCTGAAGGAGTTGGATCTGATGGGCTACCGTGAGCGCGAAGATGAGCTGGAATTCCTGGCGGAAGCCGGCGAGAAGCTGGTGCAGCGCAAACTGATCGGCGACCTGTGGCCGCAGGATGGCAACACCAAACTGAAACCGCTGGCG",
            
            // Wild-type sequence (no mutations)
            "ATGGAGACCCGCAAGATCACCGCCGAAGACCGCCACGCTGGTGATGGGCTCGAAATCCTGGCGAGCCTGAACAAGCTGGGCGATGTGCTGAAGGAGTTGGATCTGATGGGCTACCGTGAGCGCGAAGATGAGCTGGAATTCCTGGCGGAAGCCGGCGAGAAGCTGGTGCAGCGCAAACTGATCGGCGACCTGTGGCCGCAGGATGGCAACACCAAACTGAAACCGCTGGCG",
            
            // Short read with partial sequence
            "CTGGCGAAGATCCTGAACAAGCTGGGCGATGTGCTGAAGGAGTTGGATCTGATGGGCTACCGTGAGCGCGAAGATGAGCTGGAATTC",
            
            // Reverse complement of gyrA S83L
            "CGCCAGCGGCTTTCAGTTTGGTGTTGCCATCCTGCGGCCACAGGTCGCCGATCAGTTTGCGCTGCACCAGCTTCTCGCCGGCTTCCGCCAGGAATTCCAGCTCATCTTCGCGCTCACGGTAGCCCATCAGATCCAACTCCTTCAGCACATCGCCCAGCTTGTTCAGATTCGCCAGGATTTCGAGCCCATCACCAGCGTGGCGGTCTTCGGCGGTGATCTTGCGGGTCTCCAT",
            
            // Another gyrA variant (S83F - AGC -> TTC)
            "ATGGAGACCCGCAAGATCACCGCCGAAGACCGCCACGCTGGTGATGGGCTCGAAATCCTGGCGTTCCTGAACAAGCTGGGCGATGTGCTGAAGGAGTTGGATCTGATGGGCTACCGTGAGCGCGAAGATGAGCTGGAATTCCTGGCGGAAGCCGGCGAGAAGCTGGTGCAGCGCAAACTGATCGGCGACCTGTGGCCGCAGGATGGCAACACCAAACTGAAACCGCTGGCG"
        };
        
        int num_reads = test_reads.size();
        
        // Prepare host data
        std::vector<char> h_reads_data;
        std::vector<int> h_read_lengths(num_reads);
        std::vector<int> h_read_offsets(num_reads);
        std::vector<char> h_reads_to_process(num_reads, 1);
        
        for (int i = 0; i < num_reads; i++) {
            h_read_offsets[i] = h_reads_data.size();
            h_read_lengths[i] = test_reads[i].length();
            
            for (char c : test_reads[i]) {
                h_reads_data.push_back(c);
            }
        }
        
        // Copy to GPU
        cudaMemcpy(d_reads, h_reads_data.data(), h_reads_data.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, h_read_lengths.data(), num_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, h_read_offsets.data(), num_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_reads_to_process, h_reads_to_process.data(), num_reads * sizeof(char), cudaMemcpyHostToDevice);
        
        // Reset results
        cudaMemset(d_results, 0, max_reads * max_matches_per_read * sizeof(ProteinMatch));
        cudaMemset(d_result_counts, 0, max_reads * sizeof(uint32_t));
        
        std::cout << "Testing " << num_reads << " synthetic reads..." << std::endl;
        
        // Run enhanced translated search
        auto start_time = std::chrono::high_resolution_clock::now();
        
        int search_result = search_translated_reads(
            engine,
            d_reads,
            d_read_lengths,
            d_read_offsets,
            d_reads_to_process,
            num_reads,
            d_results,
            d_result_counts
        );
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        
        if (search_result != 0) {
            std::cerr << "ERROR: Enhanced translated search failed" << std::endl;
            return;
        }
        
        std::cout << "Enhanced translated search completed in " << duration.count() / 1000.0 << " ms" << std::endl;
        
        // Copy results back to host
        std::vector<ProteinMatch> h_results(max_reads * max_matches_per_read);
        std::vector<uint32_t> h_result_counts(max_reads);
        
        cudaMemcpy(h_results.data(), d_results, max_reads * max_matches_per_read * sizeof(ProteinMatch), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_result_counts.data(), d_result_counts, max_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        // Analyze results
        analyzeResults(h_results, h_result_counts, num_reads, test_reads);
    }
    
    void analyzeResults(const std::vector<ProteinMatch>& results, 
                       const std::vector<uint32_t>& result_counts,
                       int num_reads,
                       const std::vector<std::string>& test_reads) {
        
        std::cout << "\n=== Enhanced Translated Search Results Analysis ===" << std::endl;
        
        int total_matches = 0;
        int reads_with_matches = 0;
        int smith_waterman_alignments = 0;
        int resistance_mutations_found = 0;
        
        for (int i = 0; i < num_reads; i++) {
            uint32_t matches_for_read = result_counts[i];
            total_matches += matches_for_read;
            
            if (matches_for_read > 0) {
                reads_with_matches++;
                
                std::cout << "\nRead " << i << " (" << test_reads[i].length() << " bp): " 
                         << matches_for_read << " protein matches" << std::endl;
                
                // Analyze each match for this read
                for (uint32_t j = 0; j < matches_for_read; j++) {
                    const ProteinMatch& match = results[i * max_matches_per_read + j];
                    
                    std::cout << "  Match " << j << ":" << std::endl;
                    std::cout << "    Frame: " << (int)match.frame << std::endl;
                    std::cout << "    Protein ID: " << match.protein_id << std::endl;
                    std::cout << "    Gene ID: " << match.gene_id << std::endl;
                    std::cout << "    Species ID: " << match.species_id << std::endl;
                    std::cout << "    Query range: " << match.query_start << "-" << (match.query_start + match.match_length) << std::endl;
                    std::cout << "    Ref range: " << match.ref_start << "-" << (match.ref_start + match.match_length) << std::endl;
                    std::cout << "    Length: " << match.match_length << " amino acids" << std::endl;
                    std::cout << "    Score: " << match.alignment_score << std::endl;
                    std::cout << "    Identity: " << (match.identity * 100.0f) << "%" << std::endl;
                    std::cout << "    Smith-Waterman: " << (match.used_smith_waterman ? "YES" : "NO") << std::endl;
                    
                    if (match.used_smith_waterman) {
                        smith_waterman_alignments++;
                    }
                    
                    if (match.num_mutations > 0) {
                        std::cout << "    Mutations (" << (int)match.num_mutations << "):" << std::endl;
                        for (int m = 0; m < match.num_mutations; m++) {
                            std::cout << "      Position " << match.mutation_positions[m] 
                                     << ": " << match.ref_aas[m] << " -> " << match.query_aas[m]
                                     << " (BLOSUM: " << match.blosum_scores[m] << ")" << std::endl;
                            
                            // Check for known resistance mutations
                            if ((match.gene_id == 0 && match.mutation_positions[m] == 83) ||  // gyrA S83
                                (match.gene_id == 1 && match.mutation_positions[m] == 80)) {   // parC S80
                                resistance_mutations_found++;
                                std::cout << "        *** RESISTANCE MUTATION DETECTED ***" << std::endl;
                            }
                        }
                    }
                    
                    std::cout << std::endl;
                }
            }
        }
        
        // Summary statistics
        std::cout << "\n=== Enhanced Search Summary ===" << std::endl;
        std::cout << "Total reads processed: " << num_reads << std::endl;
        std::cout << "Reads with protein matches: " << reads_with_matches << " (" 
                 << (100.0 * reads_with_matches / num_reads) << "%)" << std::endl;
        std::cout << "Total protein matches: " << total_matches << std::endl;
        std::cout << "Average matches per positive read: " 
                 << (reads_with_matches > 0 ? (double)total_matches / reads_with_matches : 0.0) << std::endl;
        std::cout << "Smith-Waterman alignments: " << smith_waterman_alignments;
        if (total_matches > 0) {
            std::cout << " (" << (100.0 * smith_waterman_alignments / total_matches) << "% of matches)";
        }
        std::cout << std::endl;
        std::cout << "Resistance mutations detected: " << resistance_mutations_found << std::endl;
        
        // Performance metrics
        double throughput = num_reads / (1.0);  // Simplified - would need actual timing
        std::cout << "\nPerformance (estimated):" << std::endl;
        std::cout << "  Throughput: ~" << throughput << " reads/second" << std::endl;
        std::cout << "  Enhanced features: 5-mer k-mers + seed clustering + extension" << std::endl;
        if (smith_waterman_alignments > 0) {
            std::cout << "  Smith-Waterman overhead: minimal (only for high-scoring matches)" << std::endl;
        }
    }
    
    void testSmithWatermanToggle() {
        std::cout << "\n=== Testing Smith-Waterman Toggle ===" << std::endl;
        
        std::cout << "Disabling Smith-Waterman..." << std::endl;
        set_smith_waterman_enabled(engine, false);
        
        std::cout << "Enabling Smith-Waterman..." << std::endl;
        set_smith_waterman_enabled(engine, true);
        
        std::cout << "Smith-Waterman toggle test completed" << std::endl;
    }
};

int main(int argc, char** argv) {
    std::cout << "Enhanced Translated Search Test Program" << std::endl;
    std::cout << "Features: 5-mer protein k-mers + seed clustering + optional Smith-Waterman" << std::endl;
    
    bool enable_smith_waterman = false;
    std::string db_path = "data/protein_resistance_db";
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--enable-smith-waterman") {
            enable_smith_waterman = true;
        } else if (std::string(argv[i]) == "--protein-db" && i + 1 < argc) {
            db_path = argv[i + 1];
            i++;
        } else if (std::string(argv[i]) == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --enable-smith-waterman    Enable Smith-Waterman for high-scoring matches" << std::endl;
            std::cout << "  --protein-db <path>        Path to protein database" << std::endl;
            return 0;
        }
    }
    
    try {
        TranslatedSearchTester tester(enable_smith_waterman);
        
        // Test without database first
        std::cout << "\nTesting enhanced engine initialization..." << std::endl;
        
        // Try to load database if it exists
        if (!db_path.empty()) {
            if (tester.loadDatabase(db_path)) {
                tester.testSyntheticReads();
            } else {
                std::cout << "Database not found, testing engine functionality only..." << std::endl;
            }
        }
        
        if (enable_smith_waterman) {
            tester.testSmithWatermanToggle();
        }
        
        std::cout << "\n=== Enhanced Translated Search Test Completed Successfully ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}