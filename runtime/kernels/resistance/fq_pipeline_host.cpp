// fq_pipeline_host.cpp
// Main host code for FQ resistance detection pipeline with enhanced translated search support

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <zlib.h>
#include <chrono>
#include <cuda_runtime.h>
#include "fq_mutation_detector.cuh"
#include "hdf5_alignment_writer.h"

// Include Bloom filter declarations
extern "C" {
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int bloom_filter_screen_reads(void* filter, const char* d_reads, const int* d_read_lengths,
                                  const int* d_read_offsets, int num_reads, bool* d_read_passes,
                                  int* d_kmers_found, int min_kmers_threshold);
    int save_bloom_filter(void* filter, const char* filename);
    int load_bloom_filter(void* filter, const char* filename);
    
    // RC-aware functions
    int bloom_filter_screen_reads_with_rc(
        void* filter,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        bool* d_read_passes,
        int* d_kmers_found,
        int min_kmers_threshold,
        bool check_rc,
        int* d_debug_stats
    );
    
    bool bloom_filter_has_rc(void* filter);
}

// Enhanced translated search declarations
extern "C" {
    void* create_translated_search_engine(int batch_size);
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

// Enhanced ProteinMatch structure from translated_search.cu
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

// Debug macros
#define DEBUG 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { fprintf(stderr, "[DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); fflush(stderr); }

// CUDA error checking
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "[CUDA ERROR] %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

// FASTQ record structure
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Read batch for GPU processing
struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

class FastqReader {
private:
    gzFile file;
    bool is_open;
    
public:
    FastqReader(const std::string& filename) : is_open(false) {
        DEBUG_PRINT("Opening FASTQ file: %s", filename.c_str());
        file = gzopen(filename.c_str(), "r");
        is_open = (file != NULL);
        if (!is_open) {
            DEBUG_PRINT("ERROR: Failed to open FASTQ file: %s", filename.c_str());
        } else {
            DEBUG_PRINT("Successfully opened FASTQ file");
        }
    }
    
    ~FastqReader() {
        if (is_open) {
            gzclose(file);
        }
    }
    
    bool readRecord(FastqRecord& record) {
        if (!is_open) return false;
        
        const int buffer_size = 1024;
        char buffer[buffer_size];
        
        // Read header
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.header = std::string(buffer);
        if (record.header.empty() || record.header[0] != '@') return false;
        
        // Read sequence
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.sequence = std::string(buffer);
        record.sequence.pop_back(); // Remove newline
        
        // Read plus line
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        
        // Read quality
        if (gzgets(file, buffer, buffer_size) == NULL) return false;
        record.quality = std::string(buffer);
        record.quality.pop_back(); // Remove newline
        
        return true;
    }
    
    bool isOpen() const { return is_open; }
};

class FQResistancePipeline {
private:
    // GPU memory for reads
    char* d_reads_r1;
    char* d_reads_r2;
    int* d_lengths_r1;
    int* d_lengths_r2;
    int* d_offsets_r1;
    int* d_offsets_r2;
    
    // Bloom filter memory
    void* bloom_filter;
    bool* d_bloom_passes_r1;
    bool* d_bloom_passes_r2;
    int* d_bloom_kmers_r1;
    int* d_bloom_kmers_r2;
    int* d_bloom_debug_stats;
    
    // Results for nucleotide search
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    AlignmentResult* d_results;
    uint32_t* d_result_count;
    
    // Enhanced translated search engine and results
    void* translated_search_engine;
    ProteinMatch* d_protein_matches;
    uint32_t* d_protein_match_counts;
    bool enable_translated_search;
    bool enable_smith_waterman;
    std::string protein_db_path;
    
    // HDF5 writer
    HDF5AlignmentWriter* hdf5_writer;
    
    // Batch parameters
    const int batch_size = 10000;
    const int max_read_length = 300;
    const int bloom_min_kmers = 3;
    const int kmer_length = 15;
    const int max_protein_matches_per_read = 32;
    
    FQMutationDetectorCUDA detector;
    std::string current_index_path;
    
    // Statistics
    int total_reads_processed = 0;
    int total_reads_passed_bloom = 0;
    int total_candidates_found = 0;
    int total_mutations_found = 0;
    int total_protein_matches_found = 0;
    int total_protein_resistance_found = 0;
    int total_smith_waterman_alignments = 0;
    
public:
    FQResistancePipeline(bool use_translated_search = false, bool use_smith_waterman = false) 
        : enable_translated_search(use_translated_search), enable_smith_waterman(use_smith_waterman) {
        DEBUG_PRINT("Initializing Enhanced FQ Resistance Pipeline");
        DEBUG_PRINT("Batch size: %d, Max read length: %d", batch_size, max_read_length);
        DEBUG_PRINT("Translated search: %s", enable_translated_search ? "ENABLED" : "DISABLED");
        DEBUG_PRINT("Smith-Waterman: %s", enable_smith_waterman ? "ENABLED" : "DISABLED");
        
        // Initialize HDF5 writer to null
        hdf5_writer = nullptr;
        
        // Create Bloom filter
        bloom_filter = create_bloom_filter(kmer_length);
        if (!bloom_filter) {
            std::cerr << "ERROR: Failed to create Bloom filter" << std::endl;
            exit(1);
        }
        DEBUG_PRINT("Bloom filter created successfully");
        
        // Initialize enhanced translated search engine
        translated_search_engine = nullptr;
        d_protein_matches = nullptr;
        d_protein_match_counts = nullptr;
        
        if (enable_translated_search) {
            translated_search_engine = create_translated_search_engine_with_sw(batch_size, enable_smith_waterman);
            if (!translated_search_engine) {
                std::cerr << "ERROR: Failed to create enhanced translated search engine" << std::endl;
                exit(1);
            }
            DEBUG_PRINT("Enhanced translated search engine created (SW: %s)", 
                       enable_smith_waterman ? "enabled" : "disabled");
            
            // Allocate memory for protein search results
            CUDA_CHECK(cudaMalloc(&d_protein_matches, 
                                  batch_size * max_protein_matches_per_read * sizeof(ProteinMatch)));
            CUDA_CHECK(cudaMalloc(&d_protein_match_counts, batch_size * sizeof(uint32_t)));
        }
        
        // Allocate GPU memory for batch processing
        size_t read_buffer_size = batch_size * max_read_length;
        DEBUG_PRINT("Allocating GPU memory: read buffer size = %zu bytes", read_buffer_size);
        
        CUDA_CHECK(cudaMalloc(&d_reads_r1, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_reads_r2, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_lengths_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_lengths_r2, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r2, batch_size * sizeof(int)));
        
        // Allocate Bloom filter results
        CUDA_CHECK(cudaMalloc(&d_bloom_passes_r1, batch_size * sizeof(bool)));
        CUDA_CHECK(cudaMalloc(&d_bloom_passes_r2, batch_size * sizeof(bool)));
        CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r2, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_bloom_debug_stats, 4 * sizeof(int)));
        
        // Allocate memory for nucleotide search results
        size_t candidates_size = batch_size * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch);
        size_t results_size = batch_size * MAX_CANDIDATES_PER_READ * sizeof(AlignmentResult);
        DEBUG_PRINT("Allocating GPU memory: candidates = %zu bytes, results = %zu bytes", 
                    candidates_size, results_size);
        
        CUDA_CHECK(cudaMalloc(&d_candidates, candidates_size));
        CUDA_CHECK(cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_results, results_size));
        CUDA_CHECK(cudaMalloc(&d_result_count, sizeof(uint32_t)));
        
        DEBUG_PRINT("GPU memory allocation completed successfully");
    }
    
    ~FQResistancePipeline() {
        // Destroy HDF5 writer
        if (hdf5_writer) {
            delete hdf5_writer;
        }
        
        // Destroy enhanced translated search engine
        if (translated_search_engine) {
            destroy_translated_search_engine(translated_search_engine);
        }
        
        // Destroy Bloom filter
        if (bloom_filter) {
            destroy_bloom_filter(bloom_filter);
        }
        
        // Free GPU memory
        cudaFree(d_reads_r1);
        cudaFree(d_reads_r2);
        cudaFree(d_lengths_r1);
        cudaFree(d_lengths_r2);
        cudaFree(d_offsets_r1);
        cudaFree(d_offsets_r2);
        cudaFree(d_bloom_passes_r1);
        cudaFree(d_bloom_passes_r2);
        cudaFree(d_bloom_kmers_r1);
        cudaFree(d_bloom_kmers_r2);
        cudaFree(d_bloom_debug_stats);
        cudaFree(d_candidates);
        cudaFree(d_candidate_counts);
        cudaFree(d_results);
        cudaFree(d_result_count);
        
        if (d_protein_matches) cudaFree(d_protein_matches);
        if (d_protein_match_counts) cudaFree(d_protein_match_counts);
    }
    
    void setProteinDatabase(const std::string& db_path) {
        protein_db_path = db_path;
        if (translated_search_engine && !protein_db_path.empty()) {
            DEBUG_PRINT("Loading enhanced protein database from: %s", protein_db_path.c_str());
            int result = load_protein_database(translated_search_engine, protein_db_path.c_str());
            if (result != 0) {
                std::cerr << "WARNING: Failed to load enhanced protein database" << std::endl;
                enable_translated_search = false;
            } else {
                std::cout << "Enhanced protein database loaded successfully (5-mer k-mers)" << std::endl;
                if (enable_smith_waterman) {
                    std::cout << "Smith-Waterman alignment enabled for high-scoring matches" << std::endl;
                }
            }
        }
    }
    
    void setSmithWatermanEnabled(bool enabled) {
        enable_smith_waterman = enabled;
        if (translated_search_engine) {
            set_smith_waterman_enabled(translated_search_engine, enabled);
            DEBUG_PRINT("Smith-Waterman alignment %s", enabled ? "ENABLED" : "DISABLED");
        }
    }
    
    void loadIndex(const std::string& index_path) {
        DEBUG_PRINT("Loading index from: %s", index_path.c_str());
        
        // Store index path for HDF5 initialization
        current_index_path = index_path;
        
        // Load k-mer index
        detector.loadIndex(index_path.c_str());
        DEBUG_PRINT("K-mer index loaded");
        
        // Build or load Bloom filter
        std::string bloom_path = index_path + "/bloom_filter.bin";
        
        // Try to load existing Bloom filter
        if (load_bloom_filter(bloom_filter, bloom_path.c_str()) == 0) {
            std::cout << "Loaded pre-built Bloom filter from: " << bloom_path << std::endl;
            
            // Check if it has RC k-mers
            bool has_rc = bloom_filter_has_rc(bloom_filter);
            std::cout << "Bloom filter contains RC k-mers: " << (has_rc ? "YES" : "NO") << std::endl;
            
            if (!has_rc) {
                std::cout << "WARNING: Bloom filter does not contain RC k-mers. Rebuilding..." << std::endl;
                // Force rebuild with RC
                if (detector.d_kmer_sorted && detector.num_kmers > 0) {
                    if (build_bloom_filter_from_index(bloom_filter, detector.d_kmer_sorted, detector.num_kmers) == 0) {
                        std::cout << "Bloom filter rebuilt with RC k-mers" << std::endl;
                        save_bloom_filter(bloom_filter, bloom_path.c_str());
                    }
                }
            }
        } else {
            // Build Bloom filter from k-mer index
            std::cout << "Building Bloom filter from k-mer index (with RC k-mers)..." << std::endl;
            
            if (detector.d_kmer_sorted && detector.num_kmers > 0) {
                if (build_bloom_filter_from_index(bloom_filter, detector.d_kmer_sorted, detector.num_kmers) == 0) {
                    std::cout << "Bloom filter built successfully with " << detector.num_kmers << " k-mers (including RC)" << std::endl;
                    
                    // Save for future use
                    save_bloom_filter(bloom_filter, bloom_path.c_str());
                    std::cout << "Saved Bloom filter to: " << bloom_path << std::endl;
                } else {
                    std::cerr << "WARNING: Failed to build Bloom filter from index" << std::endl;
                }
            } else {
                std::cerr << "WARNING: No k-mers available to build Bloom filter" << std::endl;
            }
        }
        
        DEBUG_PRINT("Index loading completed");
    }
    
    ReadBatch prepareBatch(const std::vector<FastqRecord>& records) {
        ReadBatch batch;
        batch.num_reads = records.size();
        
        // Calculate total size needed
        batch.total_bases = 0;
        for (const auto& record : records) {
            batch.total_bases += record.sequence.length();
        }
        
        // Allocate host memory
        batch.sequences = new char[batch.total_bases];
        batch.lengths = new int[batch.num_reads];
        batch.offsets = new int[batch.num_reads];
        
        // Copy sequences
        int offset = 0;
        for (int i = 0; i < batch.num_reads; i++) {
            const std::string& seq = records[i].sequence;
            memcpy(batch.sequences + offset, seq.c_str(), seq.length());
            batch.lengths[i] = seq.length();
            batch.offsets[i] = offset;
            offset += seq.length();
        }
        
        return batch;
    }
    
    void processPairedReads(const std::string& r1_path, const std::string& r2_path, 
                           const std::string& output_path) {
        DEBUG_PRINT("Starting enhanced paired read processing");
        DEBUG_PRINT("R1: %s", r1_path.c_str());
        DEBUG_PRINT("R2: %s", r2_path.c_str());
        DEBUG_PRINT("Output: %s", output_path.c_str());
        DEBUG_PRINT("5-mer protein search: %s", enable_translated_search ? "ENABLED" : "DISABLED");
        DEBUG_PRINT("Smith-Waterman: %s", enable_smith_waterman ? "ENABLED" : "DISABLED");
        
        // Create HDF5 output file
        std::string hdf5_path = output_path.substr(0, output_path.find_last_of('.')) + ".h5";
        hdf5_writer = new HDF5AlignmentWriter(hdf5_path);
        hdf5_writer->initialize(current_index_path, r1_path, r2_path);
        DEBUG_PRINT("HDF5 output initialized: %s", hdf5_path.c_str());
        
        FastqReader reader1(r1_path);
        FastqReader reader2(r2_path);
        
        if (!reader1.isOpen() || !reader2.isOpen()) {
            std::cerr << "ERROR: Failed to open input files" << std::endl;
            DEBUG_PRINT("Reader1 open: %s, Reader2 open: %s", 
                       reader1.isOpen() ? "YES" : "NO",
                       reader2.isOpen() ? "YES" : "NO");
            return;
        }
        
        DEBUG_PRINT("Both FASTQ files opened successfully");
        
        std::ofstream output(output_path);
        output << "{\n";
        output << "  \"sample\": \"" << r1_path << "\",\n";
        output << "  \"hdf5_output\": \"" << hdf5_path << "\",\n";
        output << "  \"pipeline_version\": \"0.4.0-enhanced\",\n";
        output << "  \"translated_search_enabled\": " << (enable_translated_search ? "true" : "false") << ",\n";
        output << "  \"smith_waterman_enabled\": " << (enable_smith_waterman ? "true" : "false") << ",\n";
        output << "  \"protein_kmer_size\": 5,\n";
        if (enable_translated_search && !protein_db_path.empty()) {
            output << "  \"protein_database\": \"" << protein_db_path << "\",\n";
        }
        output << "  \"mutations\": [\n";
        
        bool first_mutation = true;
        auto pipeline_start = std::chrono::high_resolution_clock::now();
        
        while (true) {
            // Read batch of paired reads
            std::vector<FastqRecord> batch_r1, batch_r2;
            
            DEBUG_PRINT("Reading batch %d (up to %d reads)", total_reads_processed/batch_size + 1, batch_size);
            
            for (int i = 0; i < batch_size; i++) {
                FastqRecord rec1, rec2;
                bool got1 = reader1.readRecord(rec1);
                bool got2 = reader2.readRecord(rec2);
                
                if (!got1 || !got2) {
                    if (i == 0) {
                        DEBUG_PRINT("No more reads available");
                    } else {
                        DEBUG_PRINT("Batch partially filled: %d reads", i);
                    }
                    break;
                }
                
                batch_r1.push_back(rec1);
                batch_r2.push_back(rec2);
            }
            
            if (batch_r1.empty()) {
                DEBUG_PRINT("Empty batch, ending processing");
                break;
            }
            
            DEBUG_PRINT("Batch contains %zu read pairs", batch_r1.size());
            
            // Prepare batches
            ReadBatch batch1 = prepareBatch(batch_r1);
            ReadBatch batch2 = prepareBatch(batch_r2);
            
            // Transfer to GPU
            DEBUG_PRINT("Transferring batch to GPU: %d reads, %d total bases (R1)", 
                       batch1.num_reads, batch1.total_bases);
            
            CUDA_CHECK(cudaMemcpy(d_reads_r1, batch1.sequences, batch1.total_bases, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lengths_r1, batch1.lengths, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_offsets_r1, batch1.offsets, batch1.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            
            CUDA_CHECK(cudaMemcpy(d_reads_r2, batch2.sequences, batch2.total_bases, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_lengths_r2, batch2.lengths, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_offsets_r2, batch2.offsets, batch2.num_reads * sizeof(int), cudaMemcpyHostToDevice));
            
            // Reset result counters
            CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
            CUDA_CHECK(cudaMemset(d_candidate_counts, 0, batch1.num_reads * sizeof(uint32_t)));
            if (enable_translated_search) {
                CUDA_CHECK(cudaMemset(d_protein_match_counts, 0, batch1.num_reads * sizeof(uint32_t)));
            }
            
            DEBUG_PRINT("GPU transfer completed");
            
            // =========================
            // Stage 0: Bloom Filter Pre-screening with RC support
            // =========================
            auto bloom_start = std::chrono::high_resolution_clock::now();
            DEBUG_PRINT("Stage 0: Bloom filter pre-screening for %d reads", batch1.num_reads);
            
            // Reset debug statistics
            CUDA_CHECK(cudaMemset(d_bloom_debug_stats, 0, 4 * sizeof(int)));
            
            // Screen R1 reads (forward orientation expected)
            int bloom_result = bloom_filter_screen_reads_with_rc(
                bloom_filter,
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                batch1.num_reads,
                d_bloom_passes_r1, d_bloom_kmers_r1,
                bloom_min_kmers,
                false,  // R1: don't check RC
                d_bloom_debug_stats
            );
            
            if (bloom_result != 0) {
                DEBUG_PRINT("WARNING: Bloom filter screening failed with code %d", bloom_result);
            }
            
            // Get R1 Bloom filter results
            std::vector<char> h_bloom_passes_raw(batch1.num_reads);
            std::vector<int> h_bloom_kmers(batch1.num_reads);
            CUDA_CHECK(cudaMemcpy(h_bloom_passes_raw.data(), d_bloom_passes_r1, 
                                 batch1.num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
            std::vector<bool> h_bloom_passes(batch1.num_reads);
            for (int i = 0; i < batch1.num_reads; i++) {
                h_bloom_passes[i] = h_bloom_passes_raw[i];
            }
            CUDA_CHECK(cudaMemcpy(h_bloom_kmers.data(), d_bloom_kmers_r1, 
                                 batch1.num_reads * sizeof(int), cudaMemcpyDeviceToHost));
            
            // Get R1 statistics
            int r1_stats[4];
            CUDA_CHECK(cudaMemcpy(r1_stats, d_bloom_debug_stats, 4 * sizeof(int), cudaMemcpyDeviceToHost));
            
            // Count R1 reads that passed
            int reads_passed_bloom_r1 = 0;
            for (int i = 0; i < batch1.num_reads; i++) {
                if (h_bloom_passes[i]) reads_passed_bloom_r1++;
            }
            
            DEBUG_PRINT("R1 Bloom filter results: %d/%d reads passed (%.1f%%)", 
                       reads_passed_bloom_r1, batch1.num_reads, 
                       100.0 * reads_passed_bloom_r1 / batch1.num_reads);
            DEBUG_PRINT("  R1 stats - Forward hits: %d, RC hits: %d", r1_stats[1], r1_stats[2]);
            
            // Screen R2 reads with RC support
            DEBUG_PRINT("Stage 0b: Bloom filter pre-screening for R2 reads (with RC)");
            
            // Reset debug statistics
            CUDA_CHECK(cudaMemset(d_bloom_debug_stats, 0, 4 * sizeof(int)));
            
            bloom_result = bloom_filter_screen_reads_with_rc(
                bloom_filter,
                d_reads_r2, d_lengths_r2, d_offsets_r2,
                batch2.num_reads,
                d_bloom_passes_r2, d_bloom_kmers_r2,
                bloom_min_kmers,
                true,   // R2: CHECK RC!
                d_bloom_debug_stats
            );
            
            // Get R2 Bloom results
            std::vector<char> h_bloom_passes_r2_raw(batch2.num_reads);
            CUDA_CHECK(cudaMemcpy(h_bloom_passes_r2_raw.data(), d_bloom_passes_r2, 
                                 batch2.num_reads * sizeof(bool), cudaMemcpyDeviceToHost));
            std::vector<bool> h_bloom_passes_r2(batch2.num_reads);
            for (int i = 0; i < batch2.num_reads; i++) {
                h_bloom_passes_r2[i] = h_bloom_passes_r2_raw[i];
            }
            
            // Get R2 statistics
            int r2_stats[4];
            CUDA_CHECK(cudaMemcpy(r2_stats, d_bloom_debug_stats, 4 * sizeof(int), cudaMemcpyDeviceToHost));
            
            int reads_passed_bloom_r2 = 0;
            for (int i = 0; i < batch2.num_reads; i++) {
                if (h_bloom_passes_r2[i]) reads_passed_bloom_r2++;
            }
            
            DEBUG_PRINT("R2 Bloom filter results: %d/%d reads passed (%.1f%%)", 
                       reads_passed_bloom_r2, batch2.num_reads, 
                       100.0 * reads_passed_bloom_r2 / batch2.num_reads);
            DEBUG_PRINT("  R2 stats - Forward hits: %d, RC hits: %d", r2_stats[1], r2_stats[2]);
            
            auto bloom_end = std::chrono::high_resolution_clock::now();
            auto bloom_time = std::chrono::duration_cast<std::chrono::microseconds>(bloom_end - bloom_start).count();
            
            // Combine R1 and R2 Bloom results - a read pair passes if EITHER read passes
            int reads_passed_bloom_combined = 0;
            std::vector<bool> h_pair_passes(batch1.num_reads);
            for (int i = 0; i < batch1.num_reads; i++) {
                h_pair_passes[i] = h_bloom_passes[i] || h_bloom_passes_r2[i];
                if (h_pair_passes[i]) reads_passed_bloom_combined++;
            }
            
            DEBUG_PRINT("Combined Bloom results: %d/%d read pairs passed (%.1f ms)", 
                       reads_passed_bloom_combined, batch1.num_reads, bloom_time / 1000.0);
            
            total_reads_passed_bloom += reads_passed_bloom_combined;
            
            // Only proceed if some read pairs passed Bloom filter
            if (reads_passed_bloom_combined > 0) {
                // We'll process both R1 and R2, then combine results
                std::map<int, AlignmentResult> best_results_per_read_pair;
                
                // =========================
                // Process R1 reads with nucleotide search
                // =========================
                DEBUG_PRINT("Processing R1 reads...");
                
                // Stage 1: K-mer Filtering for R1
                auto kmer_start = std::chrono::high_resolution_clock::now();
                DEBUG_PRINT("Stage 1: K-mer filtering for R1 reads");
                
                launch_kmer_filter(
                    d_reads_r1, d_lengths_r1, d_offsets_r1,
                    detector, batch1.num_reads,
                    d_candidates, d_candidate_counts,
                    false  // R1: don't check RC in k-mer screening
                );
                CUDA_CHECK(cudaDeviceSynchronize());
                
                // Check candidate counts for R1
                std::vector<uint32_t> h_candidate_counts_r1(batch1.num_reads);
                CUDA_CHECK(cudaMemcpy(h_candidate_counts_r1.data(), d_candidate_counts, 
                                     batch1.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                
                int total_candidates_r1 = 0;
                for (int i = 0; i < batch1.num_reads; i++) {
                    if (h_pair_passes[i] && h_candidate_counts_r1[i] > 0) {
                        total_candidates_r1 += h_candidate_counts_r1[i];
                    }
                }
                
                DEBUG_PRINT("R1 k-mer filtering: %d total candidates", total_candidates_r1);
                
                // Stage 2: Alignment for R1
                if (total_candidates_r1 > 0) {
                    DEBUG_PRINT("Stage 2: Position-weighted alignment for R1");
                    
                    CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
                    
                    launch_position_weighted_alignment(
                        d_reads_r1, d_lengths_r1, d_offsets_r1,
                        d_candidates, d_candidate_counts,
                        detector, batch1.num_reads,
                        d_results, d_result_count
                    );
                    CUDA_CHECK(cudaDeviceSynchronize());
                    
                    // Get R1 results
                    uint32_t num_results_r1;
                    CUDA_CHECK(cudaMemcpy(&num_results_r1, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                    
                    if (num_results_r1 > 0) {
                        std::vector<AlignmentResult> results_r1(num_results_r1);
                        cudaMemcpy(results_r1.data(), d_results, num_results_r1 * sizeof(AlignmentResult), 
                                  cudaMemcpyDeviceToHost);
                        
                        // Add to HDF5
                        hdf5_writer->addAlignmentBatch(results_r1.data(), num_results_r1, 
                                                     total_reads_processed);
                        
                        // Store R1 results for JSON output
                        for (const auto& result : results_r1) {
                            if (result.num_mutations_detected > 0 && h_pair_passes[result.read_id]) {
                                int read_pair_id = result.read_id;
                                
                                if (best_results_per_read_pair.find(read_pair_id) == best_results_per_read_pair.end() ||
                                    result.alignment_score > best_results_per_read_pair[read_pair_id].alignment_score) {
                                    best_results_per_read_pair[read_pair_id] = result;
                                }
                            }
                        }
                    }
                }
                
                // =========================
                // Process R2 reads with nucleotide search (with RC support)
                // =========================
                DEBUG_PRINT("Processing R2 reads (with RC support)...");
                
                // Reset candidate counts for R2
                CUDA_CHECK(cudaMemset(d_candidate_counts, 0, batch2.num_reads * sizeof(uint32_t)));
                
                // Stage 1: K-mer Filtering for R2 with RC
                DEBUG_PRINT("Stage 1: K-mer filtering for R2 reads (with RC)");
                
                launch_kmer_filter(
                    d_reads_r2, d_lengths_r2, d_offsets_r2,
                    detector, batch2.num_reads,
                    d_candidates, d_candidate_counts,
                    true   // R2: CHECK RC in k-mer screening!
                );
                CUDA_CHECK(cudaDeviceSynchronize());
                
                // Check candidate counts for R2
                std::vector<uint32_t> h_candidate_counts_r2(batch2.num_reads);
                CUDA_CHECK(cudaMemcpy(h_candidate_counts_r2.data(), d_candidate_counts, 
                                     batch2.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                
                int total_candidates_r2 = 0;
                for (int i = 0; i < batch2.num_reads; i++) {
                    if (h_pair_passes[i] && h_candidate_counts_r2[i] > 0) {
                        total_candidates_r2 += h_candidate_counts_r2[i];
                    }
                }
                
                DEBUG_PRINT("R2 k-mer filtering: %d total candidates", total_candidates_r2);
                
                // Stage 2: Alignment for R2
                if (total_candidates_r2 > 0) {
                    DEBUG_PRINT("Stage 2: Position-weighted alignment for R2");
                    
                    CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
                    
                    launch_position_weighted_alignment(
                        d_reads_r2, d_lengths_r2, d_offsets_r2,
                        d_candidates, d_candidate_counts,
                        detector, batch2.num_reads,
                        d_results, d_result_count
                    );
                    CUDA_CHECK(cudaDeviceSynchronize());
                    
                    // Get R2 results
                    uint32_t num_results_r2;
                    CUDA_CHECK(cudaMemcpy(&num_results_r2, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                    
                    if (num_results_r2 > 0) {
                        std::vector<AlignmentResult> results_r2(num_results_r2);
                        cudaMemcpy(results_r2.data(), d_results, num_results_r2 * sizeof(AlignmentResult), 
                                  cudaMemcpyDeviceToHost);
                        
                        // Add to HDF5
                        hdf5_writer->addAlignmentBatch(results_r2.data(), num_results_r2, 
                                                     total_reads_processed);
                        
                        // Combine with R1 results for JSON output, preferring higher scores
                        for (const auto& result : results_r2) {
                            if (result.num_mutations_detected > 0 && h_pair_passes[result.read_id]) {
                                int read_pair_id = result.read_id;
                                
                                // If R2 has better score than R1, use R2
                                if (best_results_per_read_pair.find(read_pair_id) == best_results_per_read_pair.end() ||
                                    result.alignment_score > best_results_per_read_pair[read_pair_id].alignment_score) {
                                    best_results_per_read_pair[read_pair_id] = result;
                                }
                            }
                        }
                    }
                }
                
                auto kmer_end = std::chrono::high_resolution_clock::now();
                auto kmer_time = std::chrono::duration_cast<std::chrono::microseconds>(kmer_end - kmer_start).count();
                
                DEBUG_PRINT("Combined nucleotide results: %zu read pairs with mutations (%.1f ms)", 
                           best_results_per_read_pair.size(), kmer_time / 1000.0);
                
                total_candidates_found += total_candidates_r1 + total_candidates_r2;
                
                // =========================
                // Stage 3: Enhanced 5-mer Translated Search (if enabled)
                // =========================
                if (enable_translated_search && translated_search_engine) {
                    auto trans_start = std::chrono::high_resolution_clock::now();
                    DEBUG_PRINT("Stage 3: Enhanced 6-frame translated search (5-mer k-mers, SW: %s)", 
                               enable_smith_waterman ? "enabled" : "disabled");
                    
                    // Search both R1 and R2 for protein matches
                    int trans_result = search_translated_reads(
                        translated_search_engine,
                        d_reads_r1, d_lengths_r1, d_offsets_r1,
                        d_bloom_passes_r1,  // Only search reads that passed bloom
                        batch1.num_reads,
                        d_protein_matches,
                        d_protein_match_counts
                    );
                    
                    if (trans_result == 0) {
                        // Get protein match counts
                        std::vector<uint32_t> h_protein_match_counts(batch1.num_reads);
                        CUDA_CHECK(cudaMemcpy(h_protein_match_counts.data(), d_protein_match_counts,
                                             batch1.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        
                        // Count total matches and Smith-Waterman alignments
                        int batch_protein_matches = 0;
                        int batch_sw_alignments = 0;
                        for (int i = 0; i < batch1.num_reads; i++) {
                            if (h_pair_passes[i]) {
                                batch_protein_matches += h_protein_match_counts[i];
                            }
                        }
                        
                        if (batch_protein_matches > 0) {
                            // Copy protein matches to host
                            std::vector<ProteinMatch> h_protein_matches(batch1.num_reads * max_protein_matches_per_read);
                            CUDA_CHECK(cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                                                 batch1.num_reads * max_protein_matches_per_read * sizeof(ProteinMatch),
                                                 cudaMemcpyDeviceToHost));
                            
                            // Add to HDF5
                            hdf5_writer->addTranslatedResults(h_protein_matches.data(), h_protein_match_counts.data(),
                                                            batch1.num_reads, total_reads_processed);
                            
                            // Count resistance mutations and Smith-Waterman usage
                            for (int i = 0; i < batch1.num_reads; i++) {
                                if (h_pair_passes[i] && h_protein_match_counts[i] > 0) {
                                    for (uint32_t j = 0; j < h_protein_match_counts[i]; j++) {
                                        const ProteinMatch& pm = h_protein_matches[i * max_protein_matches_per_read + j];
                                        total_protein_matches_found++;
                                        
                                        if (pm.used_smith_waterman) {
                                            batch_sw_alignments++;
                                        }
                                        
                                        // Check if any mutations are at resistance positions
                                        bool has_resistance = false;
                                        for (int m = 0; m < pm.num_mutations; m++) {
                                            if (pm.gene_id == 0 && (pm.mutation_positions[m] == 83 || pm.mutation_positions[m] == 87)) {
                                                has_resistance = true;
                                                break;
                                            }
                                            if (pm.gene_id == 1 && (pm.mutation_positions[m] == 80 || pm.mutation_positions[m] == 84)) {
                                                has_resistance = true;
                                                break;
                                            }
                                        }
                                        if (has_resistance) {
                                            total_protein_resistance_found++;
                                        }
                                    }
                                }
                            }
                        }
                        
                        total_smith_waterman_alignments += batch_sw_alignments;
                        DEBUG_PRINT("R1 enhanced translated search: %d protein matches found (%d used SW)", 
                                   batch_protein_matches, batch_sw_alignments);
                    }
                    
                    // Also search R2 (optional, could skip for speed)
                    CUDA_CHECK(cudaMemset(d_protein_match_counts, 0, batch2.num_reads * sizeof(uint32_t)));
                    
                    trans_result = search_translated_reads(
                        translated_search_engine,
                        d_reads_r2, d_lengths_r2, d_offsets_r2,
                        d_bloom_passes_r2,
                        batch2.num_reads,
                        d_protein_matches,
                        d_protein_match_counts
                    );
                    
                    if (trans_result == 0) {
                        std::vector<uint32_t> h_protein_match_counts(batch2.num_reads);
                        CUDA_CHECK(cudaMemcpy(h_protein_match_counts.data(), d_protein_match_counts,
                                             batch2.num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
                        
                        int batch_protein_matches_r2 = 0;
                        int batch_sw_alignments_r2 = 0;
                        for (int i = 0; i < batch2.num_reads; i++) {
                            if (h_pair_passes[i]) {
                                batch_protein_matches_r2 += h_protein_match_counts[i];
                            }
                        }
                        
                        if (batch_protein_matches_r2 > 0) {
                            std::vector<ProteinMatch> h_protein_matches(batch2.num_reads * max_protein_matches_per_read);
                            CUDA_CHECK(cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                                                 batch2.num_reads * max_protein_matches_per_read * sizeof(ProteinMatch),
                                                 cudaMemcpyDeviceToHost));
                            
                            hdf5_writer->addTranslatedResults(h_protein_matches.data(), h_protein_match_counts.data(),
                                                            batch2.num_reads, total_reads_processed);
                            
                            // Count SW usage for R2
                            for (int i = 0; i < batch2.num_reads; i++) {
                                if (h_pair_passes[i] && h_protein_match_counts[i] > 0) {
                                    for (uint32_t j = 0; j < h_protein_match_counts[i]; j++) {
                                        const ProteinMatch& pm = h_protein_matches[i * max_protein_matches_per_read + j];
                                        if (pm.used_smith_waterman) {
                                            batch_sw_alignments_r2++;
                                        }
                                    }
                                }
                            }
                            
                            total_protein_matches_found += batch_protein_matches_r2;
                            total_smith_waterman_alignments += batch_sw_alignments_r2;
                        }
                        
                        DEBUG_PRINT("R2 enhanced translated search: %d protein matches found (%d used SW)", 
                                   batch_protein_matches_r2, batch_sw_alignments_r2);
                    }
                    
                    auto trans_end = std::chrono::high_resolution_clock::now();
                    auto trans_time = std::chrono::duration_cast<std::chrono::microseconds>(trans_end - trans_start).count();
                    DEBUG_PRINT("Enhanced translated search completed (%.1f ms)", trans_time / 1000.0);
                }
                
                // Output combined results to JSON
                for (const auto& pair : best_results_per_read_pair) {
                    const auto& result = pair.second;
                    if (!first_mutation) output << ",\n";
                    output << "    {\n";
                    output << "      \"read_pair\": " << (result.read_id + total_reads_processed) << ",\n";
                    output << "      \"gene_id\": " << result.gene_id << ",\n";
                    output << "      \"species_id\": " << result.species_id << ",\n";
                    output << "      \"alignment_score\": " << result.alignment_score << ",\n";
                    output << "      \"identity\": " << result.identity << ",\n";
                    output << "      \"mutations_detected\": " << (int)result.num_mutations_detected << ",\n";
                    output << "      \"source\": \"" << (result.start_pos < 1000 ? "R1" : "R2") << "\"\n";
                    output << "    }";
                    first_mutation = false;
                    total_mutations_found++;
                }
            }
            
            total_reads_processed += batch1.num_reads;
            
            // Cleanup batch memory
            delete[] batch1.sequences;
            delete[] batch1.lengths;
            delete[] batch1.offsets;
            delete[] batch2.sequences;
            delete[] batch2.lengths;
            delete[] batch2.offsets;
            
            // Progress update
            if (total_reads_processed % 100000 == 0) {
                std::cout << "Processed " << total_reads_processed << " read pairs..." << std::endl;
                std::cout << "  Bloom filter pass rate: " << (100.0 * total_reads_passed_bloom / total_reads_processed) << "%" << std::endl;
                if (enable_translated_search) {
                    std::cout << "  Protein matches found: " << total_protein_matches_found << std::endl;
                    if (enable_smith_waterman) {
                        std::cout << "  Smith-Waterman alignments: " << total_smith_waterman_alignments << std::endl;
                    }
                }
            }
        }
        
        auto pipeline_end = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(pipeline_end - pipeline_start).count();
        
        // Finalize HDF5 output
        hdf5_writer->finalize(output_path);
        
        output << "\n  ],\n";
        output << "  \"summary\": {\n";
        output << "    \"total_reads\": " << total_reads_processed << ",\n";
        output << "    \"reads_passed_bloom\": " << total_reads_passed_bloom << ",\n";
        output << "    \"bloom_pass_rate\": " << (100.0 * total_reads_passed_bloom / total_reads_processed) << ",\n";
        output << "    \"bloom_reduction\": " << (100.0 * (total_reads_processed - total_reads_passed_bloom) / total_reads_processed) << ",\n";
        output << "    \"total_candidates\": " << total_candidates_found << ",\n";
        output << "    \"mutations_found\": " << total_mutations_found << ",\n";
        if (enable_translated_search) {
            output << "    \"protein_matches_found\": " << total_protein_matches_found << ",\n";
            output << "    \"protein_resistance_found\": " << total_protein_resistance_found << ",\n";
            output << "    \"protein_kmer_size\": 5,\n";
            if (enable_smith_waterman) {
                output << "    \"smith_waterman_alignments\": " << total_smith_waterman_alignments << ",\n";
            }
        }
        output << "    \"processing_time_seconds\": " << total_time << "\n";
        output << "  }\n";
        output << "}\n";
        output.close();
        
        std::cout << "\n=== Enhanced Processing Complete ===" << std::endl;
        std::cout << "Total reads: " << total_reads_processed << std::endl;
        std::cout << "Bloom filter:" << std::endl;
        std::cout << "  Passed: " << total_reads_passed_bloom 
                  << " (" << (100.0 * total_reads_passed_bloom / total_reads_processed) << "%)" << std::endl;
        std::cout << "  Filtered out: " << (total_reads_processed - total_reads_passed_bloom)
                  << " (" << (100.0 * (total_reads_processed - total_reads_passed_bloom) / total_reads_processed) << "%)" << std::endl;
        std::cout << "Candidates found: " << total_candidates_found << std::endl;
        std::cout << "Mutations found: " << total_mutations_found << std::endl;
        if (enable_translated_search) {
            std::cout << "Enhanced protein search (5-mer k-mers):" << std::endl;
            std::cout << "  Protein matches: " << total_protein_matches_found << std::endl;
            std::cout << "  Protein resistance mutations: " << total_protein_resistance_found << std::endl;
            if (enable_smith_waterman) {
                std::cout << "  Smith-Waterman alignments: " << total_smith_waterman_alignments << std::endl;
                double sw_rate = total_smith_waterman_alignments > 0 ? 
                    (100.0 * total_smith_waterman_alignments / total_protein_matches_found) : 0.0;
                std::cout << "  SW usage rate: " << sw_rate << "%" << std::endl;
            }
        }
        std::cout << "Total time: " << total_time << " seconds" << std::endl;
        std::cout << "Throughput: " << (total_reads_processed / (double)total_time) << " reads/second" << std::endl;
        std::cout << "\nHDF5 output: " << hdf5_path << std::endl;
    }
};

// Main entry point with enhanced command line options
int main(int argc, char** argv) {
    DEBUG_PRINT("=== Enhanced FQ Pipeline GPU Starting (v0.4.0) ===");
    DEBUG_PRINT("Features: Bloom filter + K-mer enrichment + Enhanced 5-mer protein search + Optional Smith-Waterman");
    
    if (argc < 5 || argc > 10) {
        std::cerr << "Usage: " << argv[0] << " <index_path> <reads_R1.fastq.gz> <reads_R2.fastq.gz> <output.json> [options]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --enable-translated-search     Enable 6-frame translated search with 5-mer k-mers\n";
        std::cerr << "  --enable-smith-waterman        Enable Smith-Waterman for high-scoring protein matches\n";
        std::cerr << "  --protein-db <path>            Path to protein resistance database\n";
        return 1;
    }
    
    std::string index_path = argv[1];
    std::string r1_path = argv[2];
    std::string r2_path = argv[3];
    std::string output_path = argv[4];
    
    // Parse enhanced command line options
    bool enable_translated_search = false;
    bool enable_smith_waterman = false;
    std::string protein_db_path;
    
    for (int i = 5; i < argc; i++) {
        if (std::string(argv[i]) == "--enable-translated-search") {
            enable_translated_search = true;
        } else if (std::string(argv[i]) == "--enable-smith-waterman") {
            enable_smith_waterman = true;
        } else if (std::string(argv[i]) == "--protein-db" && i + 1 < argc) {
            protein_db_path = argv[i + 1];
            i++; // Skip next argument
        }
    }
    
    // Display configuration
    std::cout << "=== Enhanced Fluoroquinolone Resistance Detection Pipeline (v0.4.0) ===" << std::endl;
    std::cout << "Features: Bloom filter + K-mer enrichment + Alignment + HDF5 output";
    if (enable_translated_search) {
        std::cout << " + Enhanced 5-mer protein search";
        if (enable_smith_waterman) {
            std::cout << " + Smith-Waterman";
        }
    }
    std::cout << std::endl;
    std::cout << "Index: " << index_path << std::endl;
    std::cout << "R1 reads: " << r1_path << std::endl;
    std::cout << "R2 reads: " << r2_path << std::endl;
    std::cout << "Output: " << output_path << std::endl;
    if (enable_translated_search) {
        std::cout << "Enhanced translated search: ENABLED (5-mer k-mers)" << std::endl;
        if (enable_smith_waterman) {
            std::cout << "Smith-Waterman alignment: ENABLED (for high-scoring matches)" << std::endl;
        }
        if (!protein_db_path.empty()) {
            std::cout << "Protein database: " << protein_db_path << std::endl;
        }
    }
    
    // Check input files exist
    DEBUG_PRINT("Checking input files...");
    std::ifstream idx_check(index_path);
    if (!idx_check.good()) {
        std::cerr << "ERROR: Index file not found: " << index_path << std::endl;
        return 1;
    }
    idx_check.close();
    
    std::ifstream r1_check(r1_path);
    if (!r1_check.good()) {
        std::cerr << "ERROR: R1 file not found: " << r1_path << std::endl;
        return 1;
    }
    r1_check.close();
    
    std::ifstream r2_check(r2_path);
    if (!r2_check.good()) {
        std::cerr << "ERROR: R2 file not found: " << r2_path << std::endl;
        return 1;
    }
    r2_check.close();
    
    // Check CUDA device
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << std::endl;
    std::cout << "Compute capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "Memory: " << prop.totalGlobalMem / (1024*1024*1024) << " GB" << std::endl;
    std::cout << std::endl;
    
    // Run enhanced pipeline
    try {
        DEBUG_PRINT("Creating enhanced pipeline instance");
        FQResistancePipeline pipeline(enable_translated_search, enable_smith_waterman);
        
        DEBUG_PRINT("Loading index");
        pipeline.loadIndex(index_path);
        
        // Load protein database if provided
        if (enable_translated_search && !protein_db_path.empty()) {
            pipeline.setProteinDatabase(protein_db_path);
        }
        
        DEBUG_PRINT("Starting enhanced paired read processing");
        pipeline.processPairedReads(r1_path, r2_path, output_path);
        
        DEBUG_PRINT("=== Enhanced FQ Pipeline GPU Completed Successfully ===");
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Exception caught: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "ERROR: Unknown exception caught" << std::endl;
        return 1;
    }
    
    return 0;
}