// clean_resistance_pipeline_main.cpp
// Clean integrated pipeline for GPU-accelerated fluoroquinolone resistance detection
// Uses existing CUDA kernels for nucleotide k-mer matching and translated protein search

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cuda_runtime.h>
#include <zlib.h>

#include "fq_mutation_detector.cuh"
#include "hdf5_alignment_writer.h"

// External CUDA kernel declarations
extern "C" {
    // Bloom filter functions
    void* create_bloom_filter(int kmer_length);
    void destroy_bloom_filter(void* filter);
    int build_bloom_filter_from_index(void* filter, const uint64_t* d_kmers, uint32_t num_kmers);
    int bloom_filter_screen_reads_with_rc(void* filter, const char* d_reads, const int* d_read_lengths,
                                          const int* d_read_offsets, int num_reads, bool* d_read_passes,
                                          int* d_kmers_found, int min_kmers_threshold, bool check_rc,
                                          int* d_debug_stats);
    
    // Translated search functions
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
    void set_smith_waterman_enabled(void* engine, bool enabled);
    int search_translated_reads(void* engine, const char* d_reads, const int* d_read_lengths,
                                const int* d_read_offsets, const bool* d_reads_to_process,
                                int num_reads, void* results, uint32_t* result_counts);
    
    // K-mer filtering functions
    void launch_kmer_filter(const char* d_reads, const int* d_read_lengths, const int* d_read_offsets,
                           FQMutationDetectorCUDA& detector, int num_reads, CandidateMatch* d_candidates,
                           uint32_t* d_candidate_counts, bool check_reverse_complement);
    
    void launch_position_weighted_alignment(const char* d_reads, const int* d_read_lengths,
                                          const int* d_read_offsets, const CandidateMatch* d_candidates,
                                          const uint32_t* d_candidate_counts, FQMutationDetectorCUDA& detector,
                                          int num_reads, AlignmentResult* d_results, uint32_t* d_result_count);
}

// Structure definitions (matching your existing code)
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

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

// CUDA error checking
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "[CUDA ERROR] %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        exit(1); \
    } \
} while(0)

class CleanResistancePipeline {
private:
    // GPU memory for reads
    char* d_reads_r1;
    char* d_reads_r2;
    int* d_lengths_r1;
    int* d_lengths_r2;
    int* d_offsets_r1;
    int* d_offsets_r2;
    
    // Bloom filter
    void* bloom_filter;
    bool* d_bloom_passes_r1;
    bool* d_bloom_passes_r2;
    int* d_bloom_kmers_r1;
    int* d_bloom_kmers_r2;
    
    // Nucleotide search results
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    AlignmentResult* d_results;
    uint32_t* d_result_count;
    
    // Protein search
    void* translated_search_engine;
    ProteinMatch* d_protein_matches;
    uint32_t* d_protein_match_counts;
    
    // Components
    FQMutationDetectorCUDA detector;
    HDF5AlignmentWriter* hdf5_writer;
    
    // Configuration
    const int batch_size = 10000;
    const int max_read_length = 300;
    const int bloom_min_kmers = 3;
    const int kmer_length = 15;
    const int max_candidates_per_read = 64;
    const int max_protein_matches_per_read = 32;
    
    // Statistics
    struct Stats {
        int total_reads = 0;
        int bloom_passed = 0;
        int nucleotide_matches = 0;
        int protein_matches = 0;
        int resistance_mutations = 0;
    } stats;

public:
    CleanResistancePipeline() {
        std::cout << "Initializing Clean Resistance Detection Pipeline\n";
        
        // Create bloom filter
        bloom_filter = create_bloom_filter(kmer_length);
        if (!bloom_filter) {
            throw std::runtime_error("Failed to create bloom filter");
        }
        
        // Create translated search engine with Smith-Waterman
        translated_search_engine = create_translated_search_engine_with_sw(batch_size, true);
        if (!translated_search_engine) {
            throw std::runtime_error("Failed to create translated search engine");
        }
        
        // Allocate GPU memory
        size_t read_buffer_size = batch_size * max_read_length;
        CUDA_CHECK(cudaMalloc(&d_reads_r1, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_reads_r2, read_buffer_size));
        CUDA_CHECK(cudaMalloc(&d_lengths_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_lengths_r2, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_offsets_r2, batch_size * sizeof(int)));
        
        // Bloom filter results
        CUDA_CHECK(cudaMalloc(&d_bloom_passes_r1, batch_size * sizeof(bool)));
        CUDA_CHECK(cudaMalloc(&d_bloom_passes_r2, batch_size * sizeof(bool)));
        CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r1, batch_size * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_bloom_kmers_r2, batch_size * sizeof(int)));
        
        // Nucleotide search results
        CUDA_CHECK(cudaMalloc(&d_candidates, batch_size * max_candidates_per_read * sizeof(CandidateMatch)));
        CUDA_CHECK(cudaMalloc(&d_candidate_counts, batch_size * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_results, batch_size * max_candidates_per_read * sizeof(AlignmentResult)));
        CUDA_CHECK(cudaMalloc(&d_result_count, sizeof(uint32_t)));
        
        // Protein search results
        CUDA_CHECK(cudaMalloc(&d_protein_matches, batch_size * max_protein_matches_per_read * sizeof(ProteinMatch)));
        CUDA_CHECK(cudaMalloc(&d_protein_match_counts, batch_size * sizeof(uint32_t)));
        
        hdf5_writer = nullptr;
    }
    
    ~CleanResistancePipeline() {
        // Clean up
        if (bloom_filter) destroy_bloom_filter(bloom_filter);
        if (translated_search_engine) destroy_translated_search_engine(translated_search_engine);
        if (hdf5_writer) delete hdf5_writer;
        
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
        cudaFree(d_candidates);
        cudaFree(d_candidate_counts);
        cudaFree(d_results);
        cudaFree(d_result_count);
        cudaFree(d_protein_matches);
        cudaFree(d_protein_match_counts);
    }
    
    void loadDatabases(const std::string& nucleotide_index_path, 
                      const std::string& protein_db_path) {
        std::cout << "Loading databases...\n";
        
        // Load nucleotide k-mer index
        detector.loadIndex(nucleotide_index_path.c_str());
        
        // Build bloom filter from k-mer index
        if (detector.d_kmer_sorted && detector.num_kmers > 0) {
            if (build_bloom_filter_from_index(bloom_filter, detector.d_kmer_sorted, 
                                            detector.num_kmers) == 0) {
                std::cout << "Bloom filter built with " << detector.num_kmers << " k-mers\n";
            }
        }
        
        // Load protein database
        if (load_protein_database(translated_search_engine, protein_db_path.c_str()) == 0) {
            std::cout << "Protein database loaded from " << protein_db_path << "\n";
        } else {
            std::cerr << "Warning: Failed to load protein database\n";
        }
    }
    
    void processReads(const std::string& r1_path, const std::string& r2_path,
                     const std::string& output_path) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Initialize HDF5 output
        std::string hdf5_path = output_path + ".h5";
        hdf5_writer = new HDF5AlignmentWriter(hdf5_path);
        hdf5_writer->initialize("", r1_path, r2_path);
        
        // Open FASTQ files
        gzFile gz_r1 = gzopen(r1_path.c_str(), "r");
        gzFile gz_r2 = gzopen(r2_path.c_str(), "r");
        
        if (!gz_r1 || !gz_r2) {
            throw std::runtime_error("Failed to open input files");
        }
        
        // Open output JSON
        std::ofstream json_output(output_path + ".json");
        json_output << "{\n";
        json_output << "  \"pipeline\": \"clean_resistance_detection\",\n";
        json_output << "  \"version\": \"1.0\",\n";
        json_output << "  \"results\": [\n";
        
        bool first_result = true;
        
        // Process in batches
        while (true) {
            std::vector<FastqRecord> batch_r1, batch_r2;
            
            if (!readBatch(gz_r1, gz_r2, batch_r1, batch_r2, batch_size)) {
                break;
            }
            
            if (batch_r1.empty()) break;
            
            // Process batch
            processBatch(batch_r1, batch_r2, json_output, first_result);
            
            // Progress update
            if (stats.total_reads % 100000 == 0) {
                std::cout << "Processed " << stats.total_reads << " reads...\n";
            }
        }
        
        // Close files
        gzclose(gz_r1);
        gzclose(gz_r2);
        
        // Finalize output
        json_output << "\n  ],\n";
        json_output << "  \"summary\": {\n";
        json_output << "    \"total_reads\": " << stats.total_reads << ",\n";
        json_output << "    \"bloom_passed\": " << stats.bloom_passed << ",\n";
        json_output << "    \"nucleotide_matches\": " << stats.nucleotide_matches << ",\n";
        json_output << "    \"protein_matches\": " << stats.protein_matches << ",\n";
        json_output << "    \"resistance_mutations\": " << stats.resistance_mutations << "\n";
        json_output << "  }\n";
        json_output << "}\n";
        json_output.close();
        
        hdf5_writer->finalize(output_path + "_summary.json");
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Print summary
        std::cout << "\n=== PROCESSING COMPLETE ===\n";
        std::cout << "Total reads: " << stats.total_reads << "\n";
        std::cout << "Bloom filter passed: " << stats.bloom_passed << " ("
                  << (100.0 * stats.bloom_passed / stats.total_reads) << "%)\n";
        std::cout << "Nucleotide matches: " << stats.nucleotide_matches << "\n";
        std::cout << "Protein matches: " << stats.protein_matches << "\n";
        std::cout << "Resistance mutations: " << stats.resistance_mutations << "\n";
        std::cout << "Processing time: " << duration.count() << " seconds\n";
        std::cout << "Output: " << output_path << ".json and " << hdf5_path << "\n";
    }

private:
    bool readBatch(gzFile gz_r1, gzFile gz_r2, std::vector<FastqRecord>& batch_r1,
                  std::vector<FastqRecord>& batch_r2, int max_size) {
        char buffer[1024];
        
        for (int i = 0; i < max_size; i++) {
            FastqRecord rec1, rec2;
            
            // Read R1
            if (gzgets(gz_r1, buffer, 1024) == NULL) return i > 0;
            rec1.header = std::string(buffer);
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
            rec1.sequence = std::string(buffer);
            rec1.sequence.pop_back(); // Remove newline
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false; // +
            if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
            rec1.quality = std::string(buffer);
            rec1.quality.pop_back();
            
            // Read R2
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.header = std::string(buffer);
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.sequence = std::string(buffer);
            rec2.sequence.pop_back();
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false; // +
            if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
            rec2.quality = std::string(buffer);
            rec2.quality.pop_back();
            
            batch_r1.push_back(rec1);
            batch_r2.push_back(rec2);
        }
        
        return true;
    }
    
    ReadBatch prepareBatch(const std::vector<FastqRecord>& records) {
        ReadBatch batch;
        batch.num_reads = records.size();
        batch.total_bases = 0;
        
        for (const auto& rec : records) {
            batch.total_bases += rec.sequence.length();
        }
        
        batch.sequences = new char[batch.total_bases];
        batch.lengths = new int[batch.num_reads];
        batch.offsets = new int[batch.num_reads];
        
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
    
    void processBatch(const std::vector<FastqRecord>& batch_r1,
                     const std::vector<FastqRecord>& batch_r2,
                     std::ofstream& json_output,
                     bool& first_result) {
        int num_reads = batch_r1.size();
        stats.total_reads += num_reads;
        
        // Prepare batches
        ReadBatch r1_batch = prepareBatch(batch_r1);
        ReadBatch r2_batch = prepareBatch(batch_r2);
        
        // Transfer to GPU
        CUDA_CHECK(cudaMemcpy(d_reads_r1, r1_batch.sequences, r1_batch.total_bases, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_lengths_r1, r1_batch.lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_offsets_r1, r1_batch.offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        
        CUDA_CHECK(cudaMemcpy(d_reads_r2, r2_batch.sequences, r2_batch.total_bases, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_lengths_r2, r2_batch.lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_offsets_r2, r2_batch.offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice));
        
        // Reset counters
        CUDA_CHECK(cudaMemset(d_result_count, 0, sizeof(uint32_t)));
        CUDA_CHECK(cudaMemset(d_candidate_counts, 0, num_reads * sizeof(uint32_t)));
        CUDA_CHECK(cudaMemset(d_protein_match_counts, 0, num_reads * sizeof(uint32_t)));
        
        // Stage 1: Bloom filter screening
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
        
        // Stage 2: Nucleotide k-mer matching
        launch_kmer_filter(d_reads_r1, d_lengths_r1, d_offsets_r1,
                          detector, num_reads, d_candidates, d_candidate_counts, false);
        
        launch_position_weighted_alignment(d_reads_r1, d_lengths_r1, d_offsets_r1,
                                         d_candidates, d_candidate_counts,
                                         detector, num_reads, d_results, d_result_count);
        
        // Get nucleotide results
        uint32_t num_results;
        CUDA_CHECK(cudaMemcpy(&num_results, d_result_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
        
        if (num_results > 0) {
            std::vector<AlignmentResult> h_results(num_results);
            CUDA_CHECK(cudaMemcpy(h_results.data(), d_results, 
                                 num_results * sizeof(AlignmentResult), cudaMemcpyDeviceToHost));
            
            stats.nucleotide_matches += num_results;
            
            // Add to HDF5
            hdf5_writer->addAlignmentBatch(h_results.data(), num_results, 
                                         stats.total_reads - num_reads);
        }
        
        // Stage 3: Protein search
        search_translated_reads(translated_search_engine,
                               d_reads_r1, d_lengths_r1, d_offsets_r1,
                               d_bloom_passes_r1, num_reads,
                               d_protein_matches, d_protein_match_counts);
        
        // Get protein results
        std::vector<uint32_t> h_protein_counts(num_reads);
        CUDA_CHECK(cudaMemcpy(h_protein_counts.data(), d_protein_match_counts,
                             num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost));
        
        int total_protein_matches = 0;
        for (int i = 0; i < num_reads; i++) {
            total_protein_matches += h_protein_counts[i];
        }
        
        if (total_protein_matches > 0) {
            std::vector<ProteinMatch> h_protein_matches(num_reads * max_protein_matches_per_read);
            CUDA_CHECK(cudaMemcpy(h_protein_matches.data(), d_protein_matches,
                                 num_reads * max_protein_matches_per_read * sizeof(ProteinMatch),
                                 cudaMemcpyDeviceToHost));
            
            // Process protein matches
            for (int i = 0; i < num_reads; i++) {
                for (uint32_t j = 0; j < h_protein_counts[i]; j++) {
                    ProteinMatch& match = h_protein_matches[i * max_protein_matches_per_read + j];
                    
                    // Count ALL protein matches, not just Smith-Waterman
                    stats.protein_matches++;
                    
                    // Debug: print first few matches
                    if (stats.protein_matches <= 10) {
                        std::cout << "DEBUG: Protein match " << stats.protein_matches << ": "
                                  << "gene_id=" << match.gene_id 
                                  << ", identity=" << match.identity
                                  << ", num_mutations=" << (int)match.num_mutations
                                  << ", match_length=" << match.match_length
                                  << ", used_sw=" << match.used_smith_waterman << "\n";
                    }
                    
                    // Check for resistance mutations in ALL matches
                    if (match.num_mutations > 0) {
                        stats.resistance_mutations++;
                            
                            // Output to JSON
                            if (!first_result) json_output << ",\n";
                            json_output << "    {\n";
                            json_output << "      \"read_id\": " << (stats.total_reads - num_reads + i) << ",\n";
                            json_output << "      \"gene_id\": " << match.gene_id << ",\n";
                            json_output << "      \"species_id\": " << match.species_id << ",\n";
                            json_output << "      \"frame\": " << (int)match.frame << ",\n";
                            json_output << "      \"alignment_score\": " << match.alignment_score << ",\n";
                            json_output << "      \"identity\": " << match.identity << ",\n";
                            json_output << "      \"mutations\": [\n";
                            
                            for (int k = 0; k < match.num_mutations; k++) {
                                if (k > 0) json_output << ",\n";
                                json_output << "        {\n";
                                json_output << "          \"position\": " << (match.ref_start + match.mutation_positions[k]) << ",\n";
                                json_output << "          \"ref_aa\": \"" << match.ref_aas[k] << "\",\n";
                                json_output << "          \"query_aa\": \"" << match.query_aas[k] << "\"\n";
                                json_output << "        }";
                            }
                            
                            json_output << "\n      ],\n";
                            json_output << "      \"peptide\": \"" << match.query_peptide << "\"\n";
                            json_output << "    }";
                            first_result = false;
                        }
                }
            }
            
            // Add to HDF5
            hdf5_writer->addTranslatedResults(h_protein_matches.data(), h_protein_counts.data(),
                                            num_reads, stats.total_reads - num_reads);
        }
        
        // Clean up batch memory
        delete[] r1_batch.sequences;
        delete[] r1_batch.lengths;
        delete[] r1_batch.offsets;
        delete[] r2_batch.sequences;
        delete[] r2_batch.lengths;
        delete[] r2_batch.offsets;
    }
};

// Main function
int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <nucleotide_index> <protein_db> <reads_r1.fastq.gz> <reads_r2.fastq.gz> [output_prefix]\n";
        std::cerr << "\nRequired:\n";
        std::cerr << "  nucleotide_index: Directory containing k-mer index\n";
        std::cerr << "  protein_db: Directory containing protein database\n";
        std::cerr << "  reads_r1.fastq.gz: Forward reads\n";
        std::cerr << "  reads_r2.fastq.gz: Reverse reads\n";
        std::cerr << "\nOptional:\n";
        std::cerr << "  output_prefix: Output file prefix (default: resistance_results)\n";
        return 1;
    }
    
    std::string nucleotide_index = argv[1];
    std::string protein_db = argv[2];
    std::string r1_path = argv[3];
    std::string r2_path = argv[4];
    std::string output_prefix = (argc > 5) ? argv[5] : "resistance_results";
    
    std::cout << "=== Clean Fluoroquinolone Resistance Detection Pipeline ===\n";
    std::cout << "Nucleotide index: " << nucleotide_index << "\n";
    std::cout << "Protein database: " << protein_db << "\n";
    std::cout << "R1 reads: " << r1_path << "\n";
    std::cout << "R2 reads: " << r2_path << "\n";
    std::cout << "Output prefix: " << output_prefix << "\n\n";
    
    // Check CUDA device
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) {
        std::cerr << "No CUDA devices found!\n";
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << "\n";
    std::cout << "Memory: " << (prop.totalGlobalMem / (1024*1024*1024)) << " GB\n\n";
    
    try {
        // Create and run pipeline
        CleanResistancePipeline pipeline;
        pipeline.loadDatabases(nucleotide_index, protein_db);
        pipeline.processReads(r1_path, r2_path, output_prefix);
        
        std::cout << "\nPipeline completed successfully!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}