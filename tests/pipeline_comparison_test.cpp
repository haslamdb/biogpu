// pipeline_comparison_test.cpp
// Comprehensive test to compare different filtering combinations in the FQ resistance pipeline
// Tests: 1) Full pipeline, 2) Skip bloom filter, 3) Skip k-mer matching

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <cuda_runtime.h>
#include <zlib.h>

// Include necessary headers from your pipeline
#include "fq_mutation_detector.cuh"

// External function declarations
extern "C" {
    // Bloom filter functions
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int bloom_filter_screen_reads_with_rc(void* filter, const char* d_reads, const int* d_read_lengths,
                                         const int* d_read_offsets, int num_reads, bool* d_read_passes,
                                         int* d_kmers_found, int min_kmers_threshold, bool check_rc,
                                         int* d_debug_stats);
    
    // K-mer matching functions
    void launch_kmer_filter(const char* d_reads, const int* d_read_lengths, const int* d_read_offsets,
                           FQMutationDetectorCUDA& detector, int num_reads, CandidateMatch* d_candidates,
                           uint32_t* d_candidate_counts, bool check_reverse_complement);
    
    // Translated search functions
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
    int search_translated_reads(void* engine, const char* d_reads, const int* d_read_lengths,
                               const int* d_read_offsets, const bool* d_reads_to_process,
                               int num_reads, void* results, uint32_t* result_counts);
}

// Structure to hold pipeline statistics
struct PipelineStats {
    int total_reads = 0;
    int reads_after_bloom = 0;
    int reads_after_kmer = 0;
    int total_protein_matches = 0;
    int unique_proteins_found = 0;
    int resistance_mutations_found = 0;
    double bloom_time_ms = 0.0;
    double kmer_time_ms = 0.0;
    double translation_time_ms = 0.0;
    double total_time_ms = 0.0;
    
    // Additional metrics
    double reads_per_second = 0.0;
    double bloom_reduction_rate = 0.0;
    double kmer_reduction_rate = 0.0;
    double protein_hit_rate = 0.0;
};

// Protein match structure (matching your translated_search.cu)
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
    char query_peptide[51];
    bool is_qrdr_alignment;
};

// FASTQ record structure
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

// Function to read FASTQ records
std::vector<FastqRecord> readFastqBatch(gzFile file, int batch_size) {
    std::vector<FastqRecord> records;
    const int buffer_size = 1024;
    char buffer[buffer_size];
    
    for (int i = 0; i < batch_size; i++) {
        FastqRecord record;
        
        // Read header
        if (gzgets(file, buffer, buffer_size) == NULL) break;
        record.header = std::string(buffer);
        if (record.header.empty() || record.header[0] != '@') break;
        
        // Read sequence
        if (gzgets(file, buffer, buffer_size) == NULL) break;
        record.sequence = std::string(buffer);
        record.sequence.pop_back(); // Remove newline
        
        // Read plus line
        if (gzgets(file, buffer, buffer_size) == NULL) break;
        
        // Read quality
        if (gzgets(file, buffer, buffer_size) == NULL) break;
        record.quality = std::string(buffer);
        record.quality.pop_back(); // Remove newline
        
        records.push_back(record);
    }
    
    return records;
}

// Pipeline execution modes
enum PipelineMode {
    FULL_PIPELINE,      // Bloom -> K-mer -> Translation
    SKIP_BLOOM,         // K-mer -> Translation
    SKIP_KMER,          // Bloom -> Translation
    TRANSLATION_ONLY    // Translation only (baseline)
};

// Function to run pipeline with specific configuration
PipelineStats runPipelineTest(
    const std::string& index_path,
    const std::string& r1_path,
    const std::string& r2_path,
    const std::string& protein_db_path,
    PipelineMode mode,
    int batch_size = 10000,
    bool enable_smith_waterman = false
) {
    PipelineStats stats;
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Initialize components based on mode
    void* bloom_filter = nullptr;
    FQMutationDetectorCUDA* detector = nullptr;
    void* translated_search_engine = nullptr;
    
    // Create bloom filter if needed
    if (mode == FULL_PIPELINE || mode == SKIP_KMER) {
        bloom_filter = create_bloom_filter(15); // 15-mer
        if (!bloom_filter) {
            std::cerr << "Failed to create bloom filter" << std::endl;
            return stats;
        }
    }
    
    // Create k-mer detector if needed
    if (mode == FULL_PIPELINE || mode == SKIP_BLOOM) {
        detector = new FQMutationDetectorCUDA();
        detector->loadIndex(index_path.c_str());
        
        // Build bloom filter from k-mer index if using full pipeline
        if (mode == FULL_PIPELINE && bloom_filter && detector->d_kmer_sorted && detector->num_kmers > 0) {
            build_bloom_filter_from_index(bloom_filter, detector->d_kmer_sorted, detector->num_kmers);
        }
    }
    
    // Always create translated search engine with Smith-Waterman enabled
    translated_search_engine = create_translated_search_engine_with_sw(batch_size, true);  // Always enable SW
    if (!translated_search_engine) {
        std::cerr << "Failed to create translated search engine" << std::endl;
        return stats;
    }
    
    // Load protein database
    if (load_protein_database(translated_search_engine, protein_db_path.c_str()) != 0) {
        std::cerr << "Failed to load protein database" << std::endl;
        return stats;
    }
    
    // Open FASTQ files
    gzFile r1_file = gzopen(r1_path.c_str(), "r");
    gzFile r2_file = gzopen(r2_path.c_str(), "r");
    
    if (!r1_file || !r2_file) {
        std::cerr << "Failed to open FASTQ files" << std::endl;
        return stats;
    }
    
    // Allocate GPU memory
    char* d_reads_r1;
    char* d_reads_r2;
    int* d_lengths_r1;
    int* d_lengths_r2;
    int* d_offsets_r1;
    int* d_offsets_r2;
    bool* d_bloom_passes_r1;
    bool* d_bloom_passes_r2;
    int* d_bloom_kmers_r1;
    int* d_bloom_kmers_r2;
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    ProteinMatch* d_protein_matches;
    uint32_t* d_protein_match_counts;
    
    size_t read_buffer_size = batch_size * 300; // Max 300bp per read
    cudaMalloc(&d_reads_r1, read_buffer_size);
    cudaMalloc(&d_reads_r2, read_buffer_size);
    cudaMalloc(&d_lengths_r1, batch_size * sizeof(int));
    cudaMalloc(&d_lengths_r2, batch_size * sizeof(int));
    cudaMalloc(&d_offsets_r1, batch_size * sizeof(int));
    cudaMalloc(&d_offsets_r2, batch_size * sizeof(int));
    cudaMalloc(&d_bloom_passes_r1, batch_size * sizeof(bool));
    cudaMalloc(&d_bloom_passes_r2, batch_size * sizeof(bool));
    cudaMalloc(&d_bloom_kmers_r1, batch_size * sizeof(int));
    cudaMalloc(&d_bloom_kmers_r2, batch_size * sizeof(int));
    cudaMalloc(&d_candidates, batch_size * 64 * sizeof(CandidateMatch));
    cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t));
    cudaMalloc(&d_protein_matches, batch_size * 32 * sizeof(ProteinMatch));
    cudaMalloc(&d_protein_match_counts, batch_size * sizeof(uint32_t));
    
    // Process batches
    std::set<uint32_t> all_proteins_found;
    
    while (true) {
        // Read batch
        auto batch_r1 = readFastqBatch(r1_file, batch_size);
        auto batch_r2 = readFastqBatch(r2_file, batch_size);
        
        if (batch_r1.empty()) break;
        
        int num_reads = batch_r1.size();
        stats.total_reads += num_reads;
        
        // Prepare batch data
        std::vector<char> h_reads_r1;
        std::vector<int> h_lengths_r1(num_reads);
        std::vector<int> h_offsets_r1(num_reads);
        
        int offset = 0;
        for (int i = 0; i < num_reads; i++) {
            h_offsets_r1[i] = offset;
            h_lengths_r1[i] = batch_r1[i].sequence.length();
            h_reads_r1.insert(h_reads_r1.end(), batch_r1[i].sequence.begin(), batch_r1[i].sequence.end());
            offset += batch_r1[i].sequence.length();
        }
        
        // Copy to GPU
        cudaMemcpy(d_reads_r1, h_reads_r1.data(), h_reads_r1.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(d_lengths_r1, h_lengths_r1.data(), num_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_offsets_r1, h_offsets_r1.data(), num_reads * sizeof(int), cudaMemcpyHostToDevice);
        
        // Initialize passes array (all true by default)
        std::vector<uint8_t> h_passes(num_reads, 1);
        cudaMemcpy(d_bloom_passes_r1, h_passes.data(), num_reads * sizeof(bool), cudaMemcpyHostToDevice);
        
        // Stage 1: Bloom filter (if enabled)
        if (mode == FULL_PIPELINE || mode == SKIP_KMER) {
            auto bloom_start = std::chrono::high_resolution_clock::now();
            
            bloom_filter_screen_reads_with_rc(
                bloom_filter,
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                num_reads,
                d_bloom_passes_r1, d_bloom_kmers_r1,
                3,     // min_kmers_threshold
                false, // check_rc
                nullptr
            );
            
            cudaDeviceSynchronize();
            auto bloom_end = std::chrono::high_resolution_clock::now();
            stats.bloom_time_ms += std::chrono::duration<double, std::milli>(bloom_end - bloom_start).count();
            
            // Get bloom results
            cudaMemcpy(h_passes.data(), d_bloom_passes_r1, num_reads * sizeof(bool), cudaMemcpyDeviceToHost);
            int passes_bloom = std::count(h_passes.begin(), h_passes.end(), (uint8_t)1);
            stats.reads_after_bloom += passes_bloom;
        } else {
            // If skipping bloom, all reads pass
            stats.reads_after_bloom += num_reads;
        }
        
        // Stage 2: K-mer matching (if enabled)
        if (mode == FULL_PIPELINE || mode == SKIP_BLOOM) {
            auto kmer_start = std::chrono::high_resolution_clock::now();
            
            cudaMemset(d_candidate_counts, 0, num_reads * sizeof(uint32_t));
            
            launch_kmer_filter(
                d_reads_r1, d_lengths_r1, d_offsets_r1,
                *detector, num_reads,
                d_candidates, d_candidate_counts,
                false // check_reverse_complement
            );
            
            cudaDeviceSynchronize();
            auto kmer_end = std::chrono::high_resolution_clock::now();
            stats.kmer_time_ms += std::chrono::duration<double, std::milli>(kmer_end - kmer_start).count();
            
            // Get k-mer results
            std::vector<uint32_t> h_candidate_counts(num_reads);
            cudaMemcpy(h_candidate_counts.data(), d_candidate_counts, num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            
            // Update passes based on k-mer hits
            for (int i = 0; i < num_reads; i++) {
                if (mode == FULL_PIPELINE) {
                    // Both bloom and k-mer must pass
                    h_passes[i] = h_passes[i] && (h_candidate_counts[i] > 0) ? 1 : 0;
                } else {
                    // Only k-mer check
                    h_passes[i] = (h_candidate_counts[i] > 0) ? 1 : 0;
                }
            }
            
            int passes_kmer = std::count(h_passes.begin(), h_passes.end(), (uint8_t)1);
            stats.reads_after_kmer += passes_kmer;
        } else {
            // If skipping k-mer, use bloom results (or all if translation only)
            stats.reads_after_kmer += (mode == SKIP_KMER) ? stats.reads_after_bloom : num_reads;
        }
        
        // Update d_bloom_passes_r1 with final passes for translated search
        cudaMemcpy(d_bloom_passes_r1, h_passes.data(), num_reads * sizeof(bool), cudaMemcpyHostToDevice);
        
        // Stage 3: Translated search
        auto trans_start = std::chrono::high_resolution_clock::now();
        
        cudaMemset(d_protein_match_counts, 0, num_reads * sizeof(uint32_t));
        
        search_translated_reads(
            translated_search_engine,
            d_reads_r1, d_lengths_r1, d_offsets_r1,
            d_bloom_passes_r1, // Only process reads that passed previous filters
            num_reads,
            d_protein_matches,
            d_protein_match_counts
        );
        
        cudaDeviceSynchronize();
        auto trans_end = std::chrono::high_resolution_clock::now();
        stats.translation_time_ms += std::chrono::duration<double, std::milli>(trans_end - trans_start).count();
        
        // Get protein results
        std::vector<uint32_t> h_protein_match_counts(num_reads);
        std::vector<ProteinMatch> h_protein_matches(num_reads * 32);
        
        cudaMemcpy(h_protein_match_counts.data(), d_protein_match_counts, 
                   num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                   num_reads * 32 * sizeof(ProteinMatch), cudaMemcpyDeviceToHost);
        
        // Analyze protein matches
        for (int i = 0; i < num_reads; i++) {
            if (h_protein_match_counts[i] > 0) {
                for (uint32_t j = 0; j < h_protein_match_counts[i]; j++) {
                    const ProteinMatch& match = h_protein_matches[i * 32 + j];
                    all_proteins_found.insert(match.protein_id);
                    stats.total_protein_matches++;
                    
                    // Count resistance mutations
                    if (match.num_mutations > 0) {
                        for (int m = 0; m < match.num_mutations; m++) {
                            // Check for known resistance positions
                            if ((match.gene_id == 0 && (match.mutation_positions[m] == 83 || match.mutation_positions[m] == 87)) ||
                                (match.gene_id == 1 && (match.mutation_positions[m] == 80 || match.mutation_positions[m] == 84))) {
                                stats.resistance_mutations_found++;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Calculate final statistics
    auto total_end = std::chrono::high_resolution_clock::now();
    stats.total_time_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    stats.unique_proteins_found = all_proteins_found.size();
    
    if (stats.total_reads > 0) {
        stats.reads_per_second = (stats.total_reads * 1000.0) / stats.total_time_ms;
        stats.bloom_reduction_rate = 1.0 - (double)stats.reads_after_bloom / stats.total_reads;
        stats.kmer_reduction_rate = 1.0 - (double)stats.reads_after_kmer / stats.total_reads;
        stats.protein_hit_rate = (double)stats.total_protein_matches / stats.total_reads;
    }
    
    // Cleanup
    gzclose(r1_file);
    gzclose(r2_file);
    
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
    cudaFree(d_candidates);
    cudaFree(d_candidate_counts);
    cudaFree(d_protein_matches);
    cudaFree(d_protein_match_counts);
    
    if (bloom_filter) destroy_bloom_filter(bloom_filter);
    if (detector) delete detector;
    if (translated_search_engine) destroy_translated_search_engine(translated_search_engine);
    
    return stats;
}

// Function to print comparison results
void printComparisonResults(const std::map<std::string, PipelineStats>& results) {
    std::cout << "\n=== PIPELINE COMPARISON RESULTS ===" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    
    // Header
    std::cout << "\n" << std::left << std::setw(25) << "Pipeline Mode"
              << std::right << std::setw(12) << "Total Time"
              << std::setw(12) << "Reads/sec"
              << std::setw(12) << "Proteins"
              << std::setw(12) << "Mutations"
              << std::setw(12) << "Hit Rate"
              << std::endl;
    std::cout << std::string(85, '-') << std::endl;
    
    // Results
    for (const auto& [mode, stats] : results) {
        std::cout << std::left << std::setw(25) << mode
                  << std::right << std::setw(12) << stats.total_time_ms << " ms"
                  << std::setw(12) << (int)stats.reads_per_second
                  << std::setw(12) << stats.unique_proteins_found
                  << std::setw(12) << stats.resistance_mutations_found
                  << std::setw(11) << stats.protein_hit_rate * 100 << "%"
                  << std::endl;
    }
    
    // Detailed timing breakdown
    std::cout << "\n=== TIMING BREAKDOWN ===" << std::endl;
    std::cout << std::left << std::setw(25) << "Pipeline Mode"
              << std::right << std::setw(15) << "Bloom (ms)"
              << std::setw(15) << "K-mer (ms)"
              << std::setw(15) << "Translation (ms)"
              << std::setw(15) << "% Translation"
              << std::endl;
    std::cout << std::string(85, '-') << std::endl;
    
    for (const auto& [mode, stats] : results) {
        double trans_percent = (stats.translation_time_ms / stats.total_time_ms) * 100;
        std::cout << std::left << std::setw(25) << mode
                  << std::right << std::setw(15) << stats.bloom_time_ms
                  << std::setw(15) << stats.kmer_time_ms
                  << std::setw(15) << stats.translation_time_ms
                  << std::setw(14) << trans_percent << "%"
                  << std::endl;
    }
    
    // Filtering efficiency
    std::cout << "\n=== FILTERING EFFICIENCY ===" << std::endl;
    std::cout << std::left << std::setw(25) << "Pipeline Mode"
              << std::right << std::setw(15) << "Total Reads"
              << std::setw(15) << "After Bloom"
              << std::setw(15) << "After K-mer"
              << std::setw(15) << "Final %"
              << std::endl;
    std::cout << std::string(85, '-') << std::endl;
    
    for (const auto& [mode, stats] : results) {
        double final_percent = (double)stats.reads_after_kmer / stats.total_reads * 100;
        std::cout << std::left << std::setw(25) << mode
                  << std::right << std::setw(15) << stats.total_reads
                  << std::setw(15) << stats.reads_after_bloom
                  << std::setw(15) << stats.reads_after_kmer
                  << std::setw(14) << final_percent << "%"
                  << std::endl;
    }
    
    // Speedup analysis
    if (results.count("Full Pipeline") && results.count("Translation Only")) {
        std::cout << "\n=== SPEEDUP ANALYSIS ===" << std::endl;
        const auto& full = results.at("Full Pipeline");
        const auto& trans_only = results.at("Translation Only");
        
        for (const auto& [mode, stats] : results) {
            if (mode != "Full Pipeline") {
                double speedup = full.total_time_ms / stats.total_time_ms;
                double protein_ratio = (double)stats.unique_proteins_found / full.unique_proteins_found;
                
                std::cout << mode << " vs Full Pipeline:" << std::endl;
                std::cout << "  Speedup: " << std::setprecision(2) << speedup << "x" << std::endl;
                std::cout << "  Protein detection ratio: " << std::setprecision(1) << protein_ratio * 100 << "%" << std::endl;
                std::cout << "  Time saved: " << (full.total_time_ms - stats.total_time_ms) << " ms" << std::endl;
                std::cout << std::endl;
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <index_path> <r1.fastq.gz> <r2.fastq.gz> <protein_db> <batch_size> [--smith-waterman]" << std::endl;
        return 1;
    }
    
    std::string index_path = argv[1];
    std::string r1_path = argv[2];
    std::string r2_path = argv[3];
    std::string protein_db_path = argv[4];
    int batch_size = std::stoi(argv[5]);
    bool enable_sw = (argc > 6 && std::string(argv[6]) == "--smith-waterman");
    
    std::cout << "=== FQ Resistance Pipeline Comparison Test ===" << std::endl;
    std::cout << "Index: " << index_path << std::endl;
    std::cout << "R1: " << r1_path << std::endl;
    std::cout << "R2: " << r2_path << std::endl;
    std::cout << "Protein DB: " << protein_db_path << std::endl;
    std::cout << "Batch size: " << batch_size << std::endl;
    std::cout << "Smith-Waterman: enabled (for all tests)" << std::endl;
    std::cout << std::endl;
    
    // Run tests
    std::map<std::string, PipelineStats> results;
    
    std::cout << "Running Full Pipeline (Bloom -> K-mer -> Translation)..." << std::endl;
    results["Full Pipeline"] = runPipelineTest(index_path, r1_path, r2_path, protein_db_path, 
                                              FULL_PIPELINE, batch_size, enable_sw);
    
    std::cout << "Running Skip Bloom (K-mer -> Translation)..." << std::endl;
    results["Skip Bloom"] = runPipelineTest(index_path, r1_path, r2_path, protein_db_path, 
                                           SKIP_BLOOM, batch_size, enable_sw);
    
    std::cout << "Running Skip K-mer (Bloom -> Translation)..." << std::endl;
    results["Skip K-mer"] = runPipelineTest(index_path, r1_path, r2_path, protein_db_path, 
                                           SKIP_KMER, batch_size, enable_sw);
    
    std::cout << "Running Translation Only..." << std::endl;
    results["Translation Only"] = runPipelineTest(index_path, r1_path, r2_path, protein_db_path, 
                                                 TRANSLATION_ONLY, batch_size, enable_sw);
    
    // Print comparison results
    printComparisonResults(results);
    
    // Recommendations
    std::cout << "\n=== RECOMMENDATIONS ===" << std::endl;
    
    // Find fastest mode that maintains good protein detection
    const auto& full_stats = results["Full Pipeline"];
    double min_time = full_stats.total_time_ms;
    std::string best_mode = "Full Pipeline";
    
    for (const auto& [mode, stats] : results) {
        double protein_ratio = (double)stats.unique_proteins_found / full_stats.unique_proteins_found;
        if (protein_ratio >= 0.95 && stats.total_time_ms < min_time) {
            min_time = stats.total_time_ms;
            best_mode = mode;
        }
    }
    
    std::cout << "Best configuration (95% protein retention): " << best_mode << std::endl;
    
    // Check if k-mer is redundant
    if (results.count("Full Pipeline") && results.count("Skip K-mer")) {
        const auto& full = results.at("Full Pipeline");
        const auto& skip_kmer = results.at("Skip K-mer");
        
        double read_diff = std::abs((double)full.reads_after_kmer - skip_kmer.reads_after_kmer) / full.total_reads;
        if (read_diff < 0.01) { // Less than 1% difference
            std::cout << "\nK-mer filtering appears redundant (<1% additional filtering)" << std::endl;
            std::cout << "Consider using Bloom -> Translation pipeline for better performance" << std::endl;
        }
    }
    
    return 0;
}