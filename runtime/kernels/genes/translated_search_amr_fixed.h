#ifndef TRANSLATED_SEARCH_AMR_FIXED_H
#define TRANSLATED_SEARCH_AMR_FIXED_H

#include "../shared/gpu_memory_pool.h"
#include <memory>
#include <string>

namespace BioGPU {

// Forward declarations
class TranslatedSearchEngineImpl;

// Improved translated search engine with proper memory management
class TranslatedSearchEngine {
public:
    // Configuration structure
    struct Config {
        int batch_size = 50000;
        bool enable_smith_waterman = false;
        int gpu_device = 0;
        bool use_memory_pool = true;
        size_t max_memory_gb = 4;  // Maximum GPU memory to use
        
        Config() = default;
    };
    
    // Constructor/Destructor
    explicit TranslatedSearchEngine(const Config& config = Config());
    ~TranslatedSearchEngine();
    
    // Disable copy, enable move
    TranslatedSearchEngine(const TranslatedSearchEngine&) = delete;
    TranslatedSearchEngine& operator=(const TranslatedSearchEngine&) = delete;
    TranslatedSearchEngine(TranslatedSearchEngine&&) noexcept;
    TranslatedSearchEngine& operator=(TranslatedSearchEngine&&) noexcept;
    
    // Initialize engine
    bool initialize();
    
    // Load protein database (only loads once)
    bool loadProteinDatabase(const std::string& db_path);
    
    // Check if database is loaded
    bool isDatabaseLoaded() const;
    
    // Search methods
    struct SearchResult {
        uint32_t read_id;
        int8_t frame;
        uint32_t protein_id;
        uint32_t gene_id;
        uint16_t query_start;
        uint16_t ref_start;
        uint16_t match_length;
        float alignment_score;
        float identity;
        uint32_t gene_length;
        uint16_t coverage_start;
        uint16_t coverage_end;
        bool used_smith_waterman;
        bool concordant;
        std::string query_peptide;
    };
    
    // Process a batch of reads
    bool processBatch(
        const char* d_reads,
        const uint32_t* d_read_offsets,
        const uint32_t* d_read_lengths,
        uint32_t num_reads,
        uint32_t max_read_length,
        std::vector<SearchResult>& results
    );
    
    // Get statistics
    struct Statistics {
        size_t total_reads_processed = 0;
        size_t total_matches_found = 0;
        size_t database_load_time_ms = 0;
        size_t total_search_time_ms = 0;
        size_t gpu_memory_used_mb = 0;
        size_t gpu_memory_peak_mb = 0;
    };
    
    Statistics getStatistics() const;
    
    // Reset engine state (keeps database loaded)
    void reset();
    
    // Clear GPU errors
    static void clearGPUErrors();
    
private:
    std::unique_ptr<TranslatedSearchEngineImpl> impl_;
};

// Singleton manager for shared protein database
class ProteinDatabaseManager {
public:
    static ProteinDatabaseManager& getInstance();
    
    // Load database (thread-safe, loads only once)
    bool loadDatabase(const std::string& db_path, int gpu_device);
    
    // Get database handle for a specific GPU
    void* getDatabaseHandle(int gpu_device);
    
    // Check if database is loaded for a GPU
    bool isDatabaseLoaded(int gpu_device) const;
    
    // Unload database from a specific GPU
    void unloadDatabase(int gpu_device);
    
    // Get database info
    struct DatabaseInfo {
        uint32_t num_proteins;
        uint32_t num_kmers;
        size_t memory_usage_bytes;
        std::string load_path;
    };
    
    DatabaseInfo getDatabaseInfo(int gpu_device) const;
    
private:
    ProteinDatabaseManager() = default;
    ~ProteinDatabaseManager();
    
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace BioGPU

// C-style wrapper functions for compatibility
extern "C" {
    void* create_translated_search_engine_fixed(int batch_size, bool enable_sw, int gpu_device);
    void destroy_translated_search_engine_fixed(void* engine);
    int load_protein_database_fixed(void* engine, const char* db_path);
    int search_translated_batch_fixed(
        void* engine,
        const char* d_reads,
        const uint32_t* d_read_offsets,
        const uint32_t* d_read_lengths,
        uint32_t num_reads,
        uint32_t max_read_length,
        void* d_matches,
        uint32_t* d_match_counts,
        uint32_t max_matches_per_read
    );
    void clear_gpu_errors_fixed();
}

#endif // TRANSLATED_SEARCH_AMR_FIXED_H