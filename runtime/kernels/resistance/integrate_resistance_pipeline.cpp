// integrate_resistance_pipeline.cpp
// Complete integration module for GPU-accelerated fluoroquinolone resistance detection
// Combines nucleotide k-mer search, translated protein search, and pileup-based variant calling

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cuda_runtime.h>
#include <zlib.h>
#include "fq_mutation_detector.cuh"
#include "hdf5_alignment_writer.h"

// Forward declarations for structures from other modules
struct ResistanceCall {
    uint32_t gene_id;
    uint16_t position;
    char wildtype_aa;
    char observed_aa;
    float allele_frequency;
    uint32_t supporting_reads;
    uint32_t total_depth;
    float confidence_score;
    char drug_affected[32];
    bool is_resistance_mutation;
};

struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;    // Position in translated frame
    uint16_t ref_start;      // Position in reference protein
    uint16_t match_length;
    float alignment_score;
    float identity;
    uint8_t num_mutations;
    uint8_t mutation_positions[10];  // Up to 10 mutations
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
    bool used_smith_waterman;  // Flag indicating if SW was used
    char query_peptide[51];  // Store aligned peptide sequence (up to 50 AA + null terminator)
    bool is_qrdr_alignment;  // Flag for QRDR region alignment
};

// External declarations for translated search
extern "C" {
    void* create_translated_search_engine_with_sw(int batch_size, bool enable_sw);
    void destroy_translated_search_engine(void* engine);
    int load_protein_database(void* engine, const char* db_path);
    void set_smith_waterman_enabled(void* engine, bool enabled);
    int search_translated_reads(void* engine, const char* d_reads, const int* d_read_lengths,
                                const int* d_read_offsets, const bool* d_reads_to_process,
                                int num_reads, void* results, uint32_t* result_counts);
}

// External declarations for resistance detection GPU
extern "C" {
    void* create_resistance_detector_gpu(int max_reads, int num_genes, int num_mutations);
    void destroy_resistance_detector_gpu(void* detector);
    int detect_resistance_gpu(void* detector, const char** translated_reads, int* read_lengths,
                             int num_reads, void* results, float min_score, float min_depth,
                             float min_allele_freq);
}

// Structure to hold integrated results
struct IntegratedResistanceResult {
    // From nucleotide alignment
    uint32_t read_id;
    std::string read_header;
    
    // From protein search
    int frame;
    uint32_t gene_id;
    std::string gene_name;
    uint32_t species_id;
    std::string species_name;
    
    // From resistance detection
    std::vector<ResistanceCall> resistance_calls;
    float max_confidence;
    bool is_resistant;
    std::string resistance_summary;
    
    // Additional metadata
    float alignment_score;
    float identity;
    int num_mutations;
    std::string aligned_peptide;
};

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

class IntegratedResistancePipeline {
private:
    // Pipeline components
    FQMutationDetectorCUDA* nucleotide_detector;
    void* translated_search_engine;
    void* resistance_detector_gpu;
    HDF5AlignmentWriter* hdf5_writer;
    
    // GPU memory for reads
    char* d_reads;
    int* d_read_lengths;
    int* d_read_offsets;
    bool* d_reads_to_process;
    
    // Results buffers
    CandidateMatch* d_candidates;
    uint32_t* d_candidate_counts;
    ProteinMatch* d_protein_matches;
    uint32_t* d_protein_match_counts;
    
    // Mappings
    std::map<uint32_t, std::string> gene_id_to_name;
    std::map<uint32_t, std::string> species_id_to_name;
    std::map<std::string, uint32_t> gene_name_to_id;
    
    // Configuration
    struct Config {
        bool enable_nucleotide_search = true;
        bool enable_protein_search = true;
        bool enable_pileup_calling = true;
        bool enable_smith_waterman = true;
        float min_alignment_score = 50.0f;
        float min_allele_frequency = 0.1f;
        float min_depth = 5.0f;
        float min_confidence = 0.8f;
        int batch_size = 10000;
        int max_candidates_per_read = 64;
        int max_protein_matches_per_read = 32;
    } config;
    
    // Statistics
    struct Stats {
        int total_reads = 0;
        int reads_with_hits = 0;
        int resistance_calls = 0;
        int high_confidence_calls = 0;
        int nucleotide_matches = 0;
        int protein_matches = 0;
        int pileup_variants = 0;
        std::chrono::milliseconds total_time;
        std::chrono::milliseconds nucleotide_time;
        std::chrono::milliseconds protein_time;
        std::chrono::milliseconds pileup_time;
    } stats;
    
public:
    IntegratedResistancePipeline() {
        nucleotide_detector = nullptr;
        translated_search_engine = nullptr;
        resistance_detector_gpu = nullptr;
        hdf5_writer = nullptr;
        
        // Allocate GPU memory
        size_t read_buffer_size = config.batch_size * 300; // Max 300 bp per read
        cudaMalloc(&d_reads, read_buffer_size);
        cudaMalloc(&d_read_lengths, config.batch_size * sizeof(int));
        cudaMalloc(&d_read_offsets, config.batch_size * sizeof(int));
        cudaMalloc(&d_reads_to_process, config.batch_size * sizeof(bool));
        
        // Allocate results buffers
        cudaMalloc(&d_candidates, config.batch_size * config.max_candidates_per_read * sizeof(CandidateMatch));
        cudaMalloc(&d_candidate_counts, config.batch_size * sizeof(uint32_t));
        cudaMalloc(&d_protein_matches, config.batch_size * config.max_protein_matches_per_read * sizeof(ProteinMatch));
        cudaMalloc(&d_protein_match_counts, config.batch_size * sizeof(uint32_t));
    }
    
    ~IntegratedResistancePipeline() {
        cleanup();
        
        // Free GPU memory
        cudaFree(d_reads);
        cudaFree(d_read_lengths);
        cudaFree(d_read_offsets);
        cudaFree(d_reads_to_process);
        cudaFree(d_candidates);
        cudaFree(d_candidate_counts);
        cudaFree(d_protein_matches);
        cudaFree(d_protein_match_counts);
    }
    
    void initialize(const std::string& index_path, 
                   const std::string& protein_db_path,
                   const std::string& resistance_db_path,
                   int max_batch_size = 10000) {
        
        std::cout << "Initializing integrated resistance detection pipeline..." << std::endl;
        
        config.batch_size = max_batch_size;
        
        // Initialize nucleotide search (existing k-mer index)
        if (config.enable_nucleotide_search) {
            nucleotide_detector = new FQMutationDetectorCUDA();
            nucleotide_detector->loadIndex(index_path.c_str());
            std::cout << "Loaded nucleotide k-mer index from: " << index_path << std::endl;
        }
        
        // Initialize protein search with Smith-Waterman
        if (config.enable_protein_search) {
            translated_search_engine = create_translated_search_engine_with_sw(
                max_batch_size, config.enable_smith_waterman);
            
            if (load_protein_database(translated_search_engine, protein_db_path.c_str()) != 0) {
                std::cerr << "ERROR: Failed to load protein database" << std::endl;
                config.enable_protein_search = false;
            } else {
                std::cout << "Loaded protein database from: " << protein_db_path << std::endl;
            }
        }
        
        // Initialize GPU resistance detector for pileup-based calling
        if (config.enable_pileup_calling) {
            // Load resistance database metadata
            loadResistanceMetadata(resistance_db_path);
            
            // Create GPU detector
            int num_genes = gene_id_to_name.size();
            int num_mutations = 100;  // Estimate
            resistance_detector_gpu = create_resistance_detector_gpu(
                max_batch_size, num_genes, num_mutations);
            
            // Load reference sequences and known mutations
            loadResistanceDatabase(resistance_db_path);
            
            std::cout << "Loaded resistance database with " << num_genes 
                     << " genes from: " << resistance_db_path << std::endl;
        }
        
        std::cout << "Pipeline initialized successfully" << std::endl;
        std::cout << "Configuration:" << std::endl;
        std::cout << "  Nucleotide search: " << (config.enable_nucleotide_search ? "ENABLED" : "DISABLED") << std::endl;
        std::cout << "  Protein search: " << (config.enable_protein_search ? "ENABLED" : "DISABLED") << std::endl;
        std::cout << "  Smith-Waterman: " << (config.enable_smith_waterman ? "ENABLED" : "DISABLED") << std::endl;
        std::cout << "  Pileup calling: " << (config.enable_pileup_calling ? "ENABLED" : "DISABLED") << std::endl;
    }
    
    std::vector<IntegratedResistanceResult> processReads(
        const std::string& r1_path,
        const std::string& r2_path,
        const std::string& output_path
    ) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::vector<IntegratedResistanceResult> all_results;
        
        // Initialize HDF5 output
        std::string hdf5_path = output_path + ".h5";
        hdf5_writer = new HDF5AlignmentWriter(hdf5_path);
        hdf5_writer->initialize(output_path, r1_path, r2_path);
        
        // Open FASTQ files
        gzFile gz_r1 = gzopen(r1_path.c_str(), "r");
        gzFile gz_r2 = gzopen(r2_path.c_str(), "r");
        
        if (!gz_r1 || !gz_r2) {
            std::cerr << "ERROR: Failed to open input files" << std::endl;
            return all_results;
        }
        
        // Process in batches
        int batch_num = 0;
        while (true) {
            // Read batch
            std::vector<FastqRecord> batch_r1, batch_r2;
            if (!readBatch(gz_r1, gz_r2, batch_r1, batch_r2, config.batch_size)) {
                break;
            }
            
            if (batch_r1.empty()) break;
            
            stats.total_reads += batch_r1.size();
            
            // Process batch through complete pipeline
            auto batch_results = processBatch(batch_r1, batch_r2, batch_num);
            
            // Accumulate results
            all_results.insert(all_results.end(), batch_results.begin(), batch_results.end());
            
            // Write intermediate results to HDF5
            writeBatchResults(batch_results, batch_num);
            
            // Progress update
            if (stats.total_reads % 100000 == 0) {
                std::cout << "Processed " << stats.total_reads << " reads..." << std::endl;
                printIntermediateStats();
            }
            
            batch_num++;
        }
        
        // Close files
        gzclose(gz_r1);
        gzclose(gz_r2);
        
        // Finalize
        hdf5_writer->finalize(output_path + "_summary.json");
        
        auto end_time = std::chrono::high_resolution_clock::now();
        stats.total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        // Print summary
        printSummary(all_results);
        
        return all_results;
    }
    
    void setConfig(const std::string& key, float value) {
        if (key == "min_alignment_score") config.min_alignment_score = value;
        else if (key == "min_allele_frequency") config.min_allele_frequency = value;
        else if (key == "min_depth") config.min_depth = value;
        else if (key == "min_confidence") config.min_confidence = value;
    }
    
private:
    bool readBatch(gzFile gz_r1, gzFile gz_r2, 
                  std::vector<FastqRecord>& batch_r1,
                  std::vector<FastqRecord>& batch_r2,
                  int batch_size) {
        char buffer[1024];
        
        for (int i = 0; i < batch_size; i++) {
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
    
    std::vector<IntegratedResistanceResult> processBatch(
        const std::vector<FastqRecord>& batch_r1,
        const std::vector<FastqRecord>& batch_r2,
        int batch_num
    ) {
        std::vector<IntegratedResistanceResult> results;
        int num_reads = batch_r1.size();
        
        // Prepare batch for GPU
        ReadBatch batch = prepareBatch(batch_r1, batch_r2);
        
        // Transfer to GPU
        cudaMemcpy(d_reads, batch.sequences, batch.total_bases, cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, batch.lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, batch.offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice);
        
        // Initialize all reads as "to process"
        std::vector<uint8_t> h_reads_to_process(num_reads, 1);
        cudaMemcpy(d_reads_to_process, h_reads_to_process.data(), 
                  num_reads * sizeof(bool), cudaMemcpyHostToDevice);
        
        // Stage 1: Nucleotide k-mer filtering (optional)
        std::vector<CandidateMatch> nucleotide_candidates;
        if (config.enable_nucleotide_search) {
            auto nuc_start = std::chrono::high_resolution_clock::now();
            nucleotide_candidates = runNucleotideSearch(batch, num_reads);
            auto nuc_end = std::chrono::high_resolution_clock::now();
            stats.nucleotide_time += std::chrono::duration_cast<std::chrono::milliseconds>(nuc_end - nuc_start);
        }
        
        // Stage 2: 6-frame translation and protein search
        std::vector<ProteinMatch> protein_matches;
        if (config.enable_protein_search) {
            auto prot_start = std::chrono::high_resolution_clock::now();
            protein_matches = runProteinSearch(num_reads);
            auto prot_end = std::chrono::high_resolution_clock::now();
            stats.protein_time += std::chrono::duration_cast<std::chrono::milliseconds>(prot_end - prot_start);
        }
        
        // Stage 3: Extract translated sequences for positive hits
        std::vector<std::string> translated_sequences;
        std::vector<int> read_mapping;  // Maps translated seq to original read
        
        for (size_t i = 0; i < protein_matches.size(); i++) {
            const auto& match = protein_matches[i];
            if (match.alignment_score >= config.min_alignment_score) {
                // Extract the aligned peptide sequence
                std::string peptide(match.query_peptide);
                translated_sequences.push_back(peptide);
                read_mapping.push_back(match.read_id);
                stats.protein_matches++;
            }
        }
        
        // Stage 4: Pileup-based variant calling (if we have translated sequences)
        std::vector<ResistanceCall> resistance_calls;
        if (config.enable_pileup_calling && !translated_sequences.empty()) {
            auto pileup_start = std::chrono::high_resolution_clock::now();
            resistance_calls = runPileupCalling(translated_sequences);
            auto pileup_end = std::chrono::high_resolution_clock::now();
            stats.pileup_time += std::chrono::duration_cast<std::chrono::milliseconds>(pileup_end - pileup_start);
        }
        
        // Stage 5: Integrate all results
        results = integrateResults(batch_r1, batch_r2, nucleotide_candidates, 
                                 protein_matches, resistance_calls, read_mapping, batch_num);
        
        // Cleanup
        delete[] batch.sequences;
        delete[] batch.lengths;
        delete[] batch.offsets;
        
        return results;
    }
    
    ReadBatch prepareBatch(const std::vector<FastqRecord>& batch_r1,
                          const std::vector<FastqRecord>& batch_r2) {
        ReadBatch batch;
        batch.num_reads = batch_r1.size() * 2; // R1 and R2
        
        // Calculate total size needed
        batch.total_bases = 0;
        for (const auto& record : batch_r1) {
            batch.total_bases += record.sequence.length();
        }
        for (const auto& record : batch_r2) {
            batch.total_bases += record.sequence.length();
        }
        
        // Allocate host memory
        batch.sequences = new char[batch.total_bases];
        batch.lengths = new int[batch.num_reads];
        batch.offsets = new int[batch.num_reads];
        
        // Copy sequences (R1 first, then R2)
        int offset = 0;
        for (int i = 0; i < batch_r1.size(); i++) {
            const std::string& seq = batch_r1[i].sequence;
            memcpy(batch.sequences + offset, seq.c_str(), seq.length());
            batch.lengths[i] = seq.length();
            batch.offsets[i] = offset;
            offset += seq.length();
        }
        
        for (int i = 0; i < batch_r2.size(); i++) {
            const std::string& seq = batch_r2[i].sequence;
            memcpy(batch.sequences + offset, seq.c_str(), seq.length());
            batch.lengths[batch_r1.size() + i] = seq.length();
            batch.offsets[batch_r1.size() + i] = offset;
            offset += seq.length();
        }
        
        return batch;
    }
    
    std::vector<CandidateMatch> runNucleotideSearch(const ReadBatch& batch, int num_reads) {
        // Reset counts
        cudaMemset(d_candidate_counts, 0, num_reads * sizeof(uint32_t));
        
        // Launch k-mer filter kernel
        extern void launch_kmer_filter(const char*, const int*, const int*, 
                                     FQMutationDetectorCUDA&, int, CandidateMatch*, 
                                     uint32_t*, bool);
        
        launch_kmer_filter(d_reads, d_read_lengths, d_read_offsets,
                         *nucleotide_detector, num_reads,
                         d_candidates, d_candidate_counts, false);
        
        cudaDeviceSynchronize();
        
        // Copy results back
        std::vector<uint32_t> h_counts(num_reads);
        cudaMemcpy(h_counts.data(), d_candidate_counts, 
                  num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::vector<CandidateMatch> all_candidates;
        for (int i = 0; i < num_reads; i++) {
            if (h_counts[i] > 0) {
                std::vector<CandidateMatch> read_candidates(h_counts[i]);
                cudaMemcpy(read_candidates.data(), 
                          d_candidates + i * config.max_candidates_per_read,
                          h_counts[i] * sizeof(CandidateMatch), 
                          cudaMemcpyDeviceToHost);
                
                all_candidates.insert(all_candidates.end(), 
                                    read_candidates.begin(), 
                                    read_candidates.end());
                
                stats.nucleotide_matches += h_counts[i];
            }
        }
        
        return all_candidates;
    }
    
    std::vector<ProteinMatch> runProteinSearch(int num_reads) {
        // Reset counts
        cudaMemset(d_protein_match_counts, 0, num_reads * sizeof(uint32_t));
        
        // Run translated search
        search_translated_reads(
            translated_search_engine,
            d_reads, d_read_lengths, d_read_offsets, d_reads_to_process,
            num_reads, d_protein_matches, d_protein_match_counts
        );
        
        cudaDeviceSynchronize();
        
        // Copy results back
        std::vector<uint32_t> h_counts(num_reads);
        cudaMemcpy(h_counts.data(), d_protein_match_counts,
                  num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        std::vector<ProteinMatch> all_matches;
        for (int i = 0; i < num_reads; i++) {
            if (h_counts[i] > 0) {
                std::vector<ProteinMatch> read_matches(h_counts[i]);
                cudaMemcpy(read_matches.data(),
                          d_protein_matches + i * config.max_protein_matches_per_read,
                          h_counts[i] * sizeof(ProteinMatch),
                          cudaMemcpyDeviceToHost);
                
                all_matches.insert(all_matches.end(),
                                 read_matches.begin(),
                                 read_matches.end());
            }
        }
        
        return all_matches;
    }
    
    std::vector<ResistanceCall> runPileupCalling(
        const std::vector<std::string>& translated_sequences
    ) {
        if (!resistance_detector_gpu || translated_sequences.empty()) {
            return std::vector<ResistanceCall>();
        }
        
        // Convert to format expected by GPU detector
        std::vector<const char*> seq_ptrs;
        std::vector<int> seq_lengths;
        
        for (const auto& seq : translated_sequences) {
            seq_ptrs.push_back(seq.c_str());
            seq_lengths.push_back(seq.length());
        }
        
        // Allocate space for results
        std::vector<ResistanceCall> calls(1000); // Max 1000 calls
        
        int num_calls = detect_resistance_gpu(
            resistance_detector_gpu,
            seq_ptrs.data(),
            seq_lengths.data(),
            translated_sequences.size(),
            calls.data(),
            config.min_alignment_score,
            config.min_depth,
            config.min_allele_frequency
        );
        
        calls.resize(num_calls);
        stats.pileup_variants += num_calls;
        
        return calls;
    }
    
    std::vector<IntegratedResistanceResult> integrateResults(
        const std::vector<FastqRecord>& batch_r1,
        const std::vector<FastqRecord>& batch_r2,
        const std::vector<CandidateMatch>& nucleotide_candidates,
        const std::vector<ProteinMatch>& protein_matches,
        const std::vector<ResistanceCall>& resistance_calls,
        const std::vector<int>& read_mapping,
        int batch_num
    ) {
        std::vector<IntegratedResistanceResult> results;
        
        // Group results by read
        std::map<int, IntegratedResistanceResult> read_results;
        
        // Process protein matches
        for (const auto& match : protein_matches) {
            if (match.alignment_score < config.min_alignment_score) continue;
            
            int read_id = match.read_id;
            auto& result = read_results[read_id];
            
            result.read_id = batch_num * config.batch_size + read_id;
            result.read_header = (read_id < batch_r1.size()) ? 
                                batch_r1[read_id].header : 
                                batch_r2[read_id - batch_r1.size()].header;
            
            result.frame = match.frame;
            result.gene_id = match.gene_id;
            result.gene_name = getGeneName(match.gene_id);
            result.species_id = match.species_id;
            result.species_name = getSpeciesName(match.species_id);
            result.alignment_score = match.alignment_score;
            result.identity = match.identity;
            result.num_mutations = match.num_mutations;
            result.aligned_peptide = std::string(match.query_peptide);
            
            // Check if this alignment used Smith-Waterman
            if (match.used_smith_waterman) {
                result.resistance_summary += "Smith-Waterman alignment detected mutations: ";
                for (int i = 0; i < match.num_mutations; i++) {
                    result.resistance_summary += match.ref_aas[i];
                    result.resistance_summary += std::to_string(match.mutation_positions[i] + 1);
                    result.resistance_summary += match.query_aas[i];
                    if (i < match.num_mutations - 1) result.resistance_summary += ", ";
                }
                result.resistance_summary += "; ";
            }
        }
        
        // Add resistance calls from pileup
        for (const auto& call : resistance_calls) {
            // Find which read this call came from
            for (size_t i = 0; i < read_mapping.size(); i++) {
                int read_id = read_mapping[i];
                auto& result = read_results[read_id];
                
                if (call.confidence_score >= config.min_confidence) {
                    result.resistance_calls.push_back(call);
                    result.is_resistant = true;
                    result.max_confidence = std::max(result.max_confidence, call.confidence_score);
                    
                    result.resistance_summary += getGeneName(call.gene_id) + " ";
                    result.resistance_summary += call.wildtype_aa;
                    result.resistance_summary += std::to_string(call.position);
                    result.resistance_summary += call.observed_aa;
                    result.resistance_summary += " (";
                    result.resistance_summary += std::to_string(int(call.allele_frequency * 100));
                    result.resistance_summary += "%, ";
                    result.resistance_summary += call.drug_affected;
                    result.resistance_summary += "); ";
                    
                    stats.resistance_calls++;
                    if (call.confidence_score >= 0.9) {
                        stats.high_confidence_calls++;
                    }
                }
            }
        }
        
        // Convert map to vector
        for (auto& pair : read_results) {
            if (pair.second.alignment_score >= config.min_alignment_score ||
                !pair.second.resistance_calls.empty()) {
                results.push_back(pair.second);
                stats.reads_with_hits++;
            }
        }
        
        return results;
    }
    
    void writeBatchResults(const std::vector<IntegratedResistanceResult>& results,
                          int batch_num) {
        // Convert to alignment results for HDF5
        std::vector<AlignmentResult> alignments;
        std::vector<ProteinMatch> protein_matches;
        
        for (const auto& result : results) {
            // Create simplified alignment result
            AlignmentResult align;
            align.read_id = result.read_id;
            align.gene_id = result.gene_id;
            align.species_id = result.species_id;
            align.alignment_score = result.alignment_score;
            align.identity = result.identity;
            align.num_mutations_detected = result.num_mutations;
            alignments.push_back(align);
            
            // Create protein match
            ProteinMatch pmatch;
            pmatch.read_id = result.read_id;
            pmatch.frame = result.frame;
            pmatch.gene_id = result.gene_id;
            pmatch.species_id = result.species_id;
            pmatch.alignment_score = result.alignment_score;
            pmatch.identity = result.identity;
            pmatch.num_mutations = result.num_mutations;
            strncpy(pmatch.query_peptide, result.aligned_peptide.c_str(), 50);
            protein_matches.push_back(pmatch);
        }
        
        if (!alignments.empty()) {
            hdf5_writer->addAlignmentBatch(alignments.data(), alignments.size(), 
                                         batch_num * config.batch_size);
        }
        
        if (!protein_matches.empty()) {
            std::vector<uint32_t> counts(protein_matches.size(), 1);
            hdf5_writer->addTranslatedResults(protein_matches.data(), counts.data(),
                                            protein_matches.size(), 
                                            batch_num * config.batch_size);
        }
    }
    
    void loadResistanceMetadata(const std::string& db_path) {
        // Load gene and species mappings from resistance_db.json
        std::ifstream metadata_file(db_path + "/resistance_db.json");
        if (!metadata_file.good()) {
            std::cerr << "WARNING: Cannot load resistance metadata from " << db_path + "/resistance_db.json" << std::endl;
            return;
        }
        
        // Simple JSON parsing (in production, use a proper JSON library)
        std::string line;
        while (std::getline(metadata_file, line)) {
            if (line.find("\"gene_map\"") != std::string::npos) {
                // Parse gene mappings
                while (std::getline(metadata_file, line) && line.find("}") == std::string::npos) {
                    size_t id_pos = line.find("\"");
                    size_t id_end = line.find("\"", id_pos + 1);
                    size_t name_pos = line.find("\"", id_end + 1);
                    size_t name_end = line.find("\"", name_pos + 1);
                    
                    if (id_pos != std::string::npos && name_pos != std::string::npos) {
                        uint32_t id = std::stoi(line.substr(id_pos + 1, id_end - id_pos - 1));
                        std::string name = line.substr(name_pos + 1, name_end - name_pos - 1);
                        gene_id_to_name[id] = name;
                        gene_name_to_id[name] = id;
                    }
                }
            }
            if (line.find("\"species_map\"") != std::string::npos) {
                // Parse species mappings
                while (std::getline(metadata_file, line) && line.find("}") == std::string::npos) {
                    size_t id_pos = line.find("\"");
                    size_t id_end = line.find("\"", id_pos + 1);
                    size_t name_pos = line.find("\"", id_end + 1);
                    size_t name_end = line.find("\"", name_pos + 1);
                    
                    if (id_pos != std::string::npos && name_pos != std::string::npos) {
                        uint32_t id = std::stoi(line.substr(id_pos + 1, id_end - id_pos - 1));
                        std::string name = line.substr(name_pos + 1, name_end - name_pos - 1);
                        species_id_to_name[id] = name;
                    }
                }
            }
        }
    }
    
    void loadResistanceDatabase(const std::string& db_path) {
        // Load resistance mutations from resistance_db.json
        std::ifstream db_file(db_path + "/resistance_db.json");
        if (!db_file.good()) {
            std::cerr << "ERROR: Cannot load resistance database from " << db_path + "/resistance_db.json" << std::endl;
            return;
        }
        
        // Load binary resistance catalog if available
        std::string catalog_path = db_path + "/resistance_catalog.bin";
        std::ifstream catalog_file(catalog_path, std::ios::binary);
        if (catalog_file.good()) {
            // Read number of mutations
            uint32_t num_mutations;
            catalog_file.read(reinterpret_cast<char*>(&num_mutations), sizeof(uint32_t));
            std::cout << "Loading " << num_mutations << " mutations from binary catalog" << std::endl;
            
            // TODO: Load actual mutation data and pass to GPU detector
            catalog_file.close();
        }
        
        // Load protein sequences from protein database
        std::string protein_metadata_path = db_path + "/protein/metadata.json";
        std::ifstream protein_meta(protein_metadata_path);
        if (protein_meta.good()) {
            std::cout << "Found protein metadata at " << protein_metadata_path << std::endl;
            protein_meta.close();
        }
        
        std::cout << "Loaded " << gene_id_to_name.size() << " genes and "
                 << species_id_to_name.size() << " species" << std::endl;
    }
    
    void cleanup() {
        if (nucleotide_detector) {
            delete nucleotide_detector;
            nucleotide_detector = nullptr;
        }
        
        if (translated_search_engine) {
            destroy_translated_search_engine(translated_search_engine);
            translated_search_engine = nullptr;
        }
        
        if (resistance_detector_gpu) {
            destroy_resistance_detector_gpu(resistance_detector_gpu);
            resistance_detector_gpu = nullptr;
        }
        
        if (hdf5_writer) {
            delete hdf5_writer;
            hdf5_writer = nullptr;
        }
    }
    
    std::string getGeneName(uint32_t gene_id) {
        auto it = gene_id_to_name.find(gene_id);
        if (it != gene_id_to_name.end()) {
            return it->second;
        }
        return "Gene" + std::to_string(gene_id);
    }
    
    std::string getSpeciesName(uint32_t species_id) {
        auto it = species_id_to_name.find(species_id);
        if (it != species_id_to_name.end()) {
            return it->second;
        }
        return "Species" + std::to_string(species_id);
    }
    
    void printIntermediateStats() {
        std::cout << "  Nucleotide matches: " << stats.nucleotide_matches << std::endl;
        std::cout << "  Protein matches: " << stats.protein_matches << std::endl;
        std::cout << "  Resistance calls: " << stats.resistance_calls << std::endl;
    }
    
    void printSummary(const std::vector<IntegratedResistanceResult>& all_results) {
        std::cout << "\n=== INTEGRATED RESISTANCE DETECTION SUMMARY ===" << std::endl;
        std::cout << "Total reads processed: " << stats.total_reads << std::endl;
        std::cout << "Reads with hits: " << stats.reads_with_hits << std::endl;
        std::cout << "Total resistance calls: " << stats.resistance_calls << std::endl;
        std::cout << "High confidence calls: " << stats.high_confidence_calls << std::endl;
        
        // Summarize by gene
        std::map<std::string, int> gene_counts;
        std::map<std::string, std::set<std::string>> gene_mutations;
        
        for (const auto& result : all_results) {
            if (result.is_resistant) {
                gene_counts[result.gene_name]++;
                for (const auto& call : result.resistance_calls) {
                    std::string mutation = std::string(1, call.wildtype_aa) + 
                                         std::to_string(call.position) + 
                                         call.observed_aa;
                    gene_mutations[result.gene_name].insert(mutation);
                }
            }
        }
        
        std::cout << "\nResistance by gene:" << std::endl;
        for (const auto& pair : gene_counts) {
            std::cout << "  " << pair.first << ": " << pair.second << " reads";
            std::cout << " (mutations: ";
            bool first = true;
            for (const auto& mut : gene_mutations[pair.first]) {
                if (!first) std::cout << ", ";
                std::cout << mut;
                first = false;
            }
            std::cout << ")" << std::endl;
        }
        
        // Timing breakdown
        std::cout << "\nProcessing time breakdown:" << std::endl;
        std::cout << "  Total time: " << stats.total_time.count() << " ms" << std::endl;
        if (config.enable_nucleotide_search) {
            std::cout << "  Nucleotide search: " << stats.nucleotide_time.count() << " ms" << std::endl;
        }
        if (config.enable_protein_search) {
            std::cout << "  Protein search: " << stats.protein_time.count() << " ms" << std::endl;
        }
        if (config.enable_pileup_calling) {
            std::cout << "  Pileup calling: " << stats.pileup_time.count() << " ms" << std::endl;
        }
        
        // Throughput
        double reads_per_second = stats.total_reads / (stats.total_time.count() / 1000.0);
        std::cout << "\nThroughput: " << reads_per_second << " reads/second" << std::endl;
    }
};

// Main function for testing
int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <kmer_index> <protein_db> <resistance_db> "
                 << "<reads_r1.fastq.gz> <reads_r2.fastq.gz> [output_prefix]" << std::endl;
        std::cerr << "\nRequired inputs:" << std::endl;
        std::cerr << "  kmer_index: Directory containing nucleotide k-mer index" << std::endl;
        std::cerr << "  protein_db: Directory containing protein database and 5-mer index" << std::endl;
        std::cerr << "  resistance_db: Directory containing resistance mutations database" << std::endl;
        std::cerr << "  reads_r1.fastq.gz: Forward reads (gzipped FASTQ)" << std::endl;
        std::cerr << "  reads_r2.fastq.gz: Reverse reads (gzipped FASTQ)" << std::endl;
        std::cerr << "\nOptional:" << std::endl;
        std::cerr << "  output_prefix: Prefix for output files (default: resistance_report)" << std::endl;
        return 1;
    }
    
    std::string kmer_index = argv[1];
    std::string protein_db = argv[2];
    std::string resistance_db = argv[3];
    std::string r1_path = argv[4];
    std::string r2_path = argv[5];
    std::string output_prefix = (argc > 6) ? argv[6] : "resistance_report";
    
    std::cout << "=== Integrated Fluoroquinolone Resistance Detection Pipeline ===" << std::endl;
    std::cout << "Version: 1.0.0" << std::endl;
    std::cout << "K-mer index: " << kmer_index << std::endl;
    std::cout << "Protein database: " << protein_db << std::endl;
    std::cout << "Resistance database: " << resistance_db << std::endl;
    std::cout << "Input R1: " << r1_path << std::endl;
    std::cout << "Input R2: " << r2_path << std::endl;
    std::cout << "Output prefix: " << output_prefix << std::endl;
    
    try {
        // Create pipeline
        IntegratedResistancePipeline pipeline;
        
        // Initialize with databases
        pipeline.initialize(kmer_index, protein_db, resistance_db);
        
        // Process reads
        auto results = pipeline.processReads(r1_path, r2_path, output_prefix);
        
        // Write detailed results to JSON
        std::string json_path = output_prefix + "_detailed.json";
        std::ofstream json_file(json_path);
        json_file << "{\n";
        json_file << "  \"total_reads_with_resistance\": " << results.size() << ",\n";
        json_file << "  \"resistance_calls\": [\n";
        
        bool first = true;
        for (const auto& result : results) {
            if (result.is_resistant) {
                if (!first) json_file << ",\n";
                json_file << "    {\n";
                json_file << "      \"read_id\": " << result.read_id << ",\n";
                json_file << "      \"gene\": \"" << result.gene_name << "\",\n";
                json_file << "      \"species\": \"" << result.species_name << "\",\n";
                json_file << "      \"confidence\": " << result.max_confidence << ",\n";
                json_file << "      \"summary\": \"" << result.resistance_summary << "\"\n";
                json_file << "    }";
                first = false;
            }
        }
        
        json_file << "\n  ]\n";
        json_file << "}\n";
        json_file.close();
        
        std::cout << "\nDetailed results written to: " << json_path << std::endl;
        std::cout << "HDF5 output: " << output_prefix << ".h5" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}