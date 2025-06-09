// Key changes to add to clean_resistance_pipeline_main.cpp
// Add these configuration flags to the CleanResistancePipeline class

// In the private section of CleanResistancePipeline class, add:
    // Configuration flags
    bool use_bloom_filter = true;        // Enable bloom filtering by default
    bool use_smith_waterman = true;      // Enable Smith-Waterman by default
    
// Update the constructor to accept these flags:
public:
    CleanResistancePipeline(bool enable_bloom = true, bool enable_sw = true) 
        : use_bloom_filter(enable_bloom), use_smith_waterman(enable_sw) {
        
        std::cout << "Initializing Clean Resistance Detection Pipeline\n";
        std::cout << "  Bloom filter: " << (use_bloom_filter ? "ENABLED" : "DISABLED") << "\n";
        std::cout << "  Smith-Waterman: " << (use_smith_waterman ? "ENABLED" : "DISABLED") << "\n";
        
        // Only create bloom filter if enabled
        if (use_bloom_filter) {
            bloom_filter = create_bloom_filter(kmer_length);
            if (!bloom_filter) {
                throw std::runtime_error("Failed to create bloom filter");
            }
        } else {
            bloom_filter = nullptr;
        }
        
        // Create translated search engine with configurable Smith-Waterman
        translated_search_engine = create_translated_search_engine_with_sw(batch_size, use_smith_waterman);
        if (!translated_search_engine) {
            throw std::runtime_error("Failed to create translated search engine");
        }
        
        // ... rest of constructor remains the same ...
    }

// Update the processBatch method to conditionally use bloom filtering:
    void processBatch(const std::vector<FastqRecord>& batch_r1,
                     const std::vector<FastqRecord>& batch_r2,
                     std::ofstream& json_output,
                     bool& first_result) {
        
        // ... existing code up to bloom filtering section ...
        
        // Stage 1: Bloom filter screening (if enabled)
        if (use_bloom_filter && bloom_filter) {
            int bloom_result_r1 = bloom_filter_screen_reads_with_rc(
                bloom_filter, d_reads_r1, d_lengths_r1, d_offsets_r1,
                num_reads, d_bloom_passes_r1, d_bloom_kmers_r1,
                bloom_min_kmers, false, nullptr
            );
            
            int bloom_result_r2 = bloom_filter_screen_reads_with_rc(
                bloom_filter, d_reads_r2, d_lengths_r2, d_offsets_r2,
                num_reads, d_bloom_passes_r2, d_bloom_kmers_r2,
                bloom_min_kmers, true, nullptr
            );
            
            // Get bloom results
            std::vector<uint8_t> h_bloom_r1(num_reads), h_bloom_r2(num_reads);
            CUDA_CHECK(cudaMemcpy(h_bloom_r1.data(), d_bloom_passes_r1, num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_bloom_r2.data(), d_bloom_passes_r2, num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
            
            // Count bloom passes
            for (int i = 0; i < num_reads; i++) {
                if (h_bloom_r1[i] || h_bloom_r2[i]) {
                    stats.bloom_passed++;
                }
            }
        } else {
            // If bloom filtering is disabled, mark all reads as passed
            CUDA_CHECK(cudaMemset(d_bloom_passes_r1, 1, num_reads * sizeof(bool)));
            CUDA_CHECK(cudaMemset(d_bloom_passes_r2, 1, num_reads * sizeof(bool)));
            stats.bloom_passed += num_reads;
        }
        
        // ... rest of the method remains the same ...
    }

// Add setter methods for runtime configuration:
    void setBloomFilterEnabled(bool enabled) {
        if (enabled && !bloom_filter) {
            bloom_filter = create_bloom_filter(kmer_length);
            if (!bloom_filter) {
                std::cerr << "Warning: Failed to create bloom filter\n";
                use_bloom_filter = false;
                return;
            }
        } else if (!enabled && bloom_filter) {
            destroy_bloom_filter(bloom_filter);
            bloom_filter = nullptr;
        }
        use_bloom_filter = enabled;
        std::cout << "Bloom filter " << (enabled ? "ENABLED" : "DISABLED") << "\n";
    }
    
    void setSmithWatermanEnabled(bool enabled) {
        use_smith_waterman = enabled;
        if (translated_search_engine) {
            set_smith_waterman_enabled(translated_search_engine, enabled);
        }
        std::cout << "Smith-Waterman " << (enabled ? "ENABLED" : "DISABLED") << "\n";
    }

// Update the main function to accept command line arguments:
int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <r1.fq.gz> <r2.fq.gz> <nucleotide_index> <protein_db> <fq_csv> <output_prefix> [--no-bloom] [--no-sw]\n";
        return 1;
    }
    
    // Parse optional flags
    bool use_bloom = true;
    bool use_sw = true;
    
    for (int i = 7; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--no-bloom") {
            use_bloom = false;
        } else if (arg == "--no-sw") {
            use_sw = false;
        }
    }
    
    try {
        CleanResistancePipeline pipeline(use_bloom, use_sw);
        
        std::string r1_path(argv[1]);
        std::string r2_path(argv[2]);
        std::string nucleotide_index(argv[3]);
        std::string protein_db(argv[4]);
        std::string fq_csv(argv[5]);
        std::string output_prefix(argv[6]);
        
        pipeline.loadDatabases(nucleotide_index, protein_db, fq_csv);
        pipeline.processReads(r1_path, r2_path, output_prefix);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}