// amr_detection_pipeline.cpp
#include "amr_detection_pipeline.h"
#include "amr_detection_kernels.h"
#include "translated_search_amr.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <cstdio>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cctype>
#include <limits>
#include <cuda_runtime.h>

// CUDA error checking macro
#define CHECK_CUDA_ERROR(msg) do { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                  << " (" << msg << ") - " << cudaGetErrorString(err) << std::endl; \
    } \
} while(0)

// Genetic code initialization is now in amr_detection_kernels_wrapper.cu


// Structure to match ProteinMatch from translated_search_amr.cu
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint16_t query_start;    // Position in translated frame
    uint16_t ref_start;      // Position in reference protein
    uint16_t match_length;
    float alignment_score;
    float identity;
    // Coverage tracking fields
    uint32_t gene_length;     // Total gene length
    uint16_t coverage_start;  // Start of covered region
    uint16_t coverage_end;    // End of covered region
    // Remove mutation-specific fields for AMR gene detection
    // (mutations are not needed for gene presence/absence)
    bool used_smith_waterman;  // Flag indicating if SW was used
    bool concordant;           // Flag for paired-end concordance
    char query_peptide[51];  // Store aligned peptide sequence (up to 50 AA + null terminator)
};

AMRDetectionPipeline::AMRDetectionPipeline(const AMRDetectionConfig& cfg) 
    : config(cfg), current_batch_size(0), total_reads_processed(0),
      translated_search_engine(nullptr), search_engine_initialized(false),
      engine_capacity((cfg.reads_per_batch * 2) + (cfg.reads_per_batch / 5)),
      reads_processed_checkpoint(0) {
    
    // Initialize performance tracking
    processing_start_time = std::chrono::steady_clock::now();
    last_performance_report = processing_start_time;
    
    // Initialize CUDA
    initializeGeneticCode();
    
    // Initialize pointers
    d_reads = nullptr;
    d_read_offsets = nullptr;
    d_read_lengths = nullptr;
    d_read_ids = nullptr;
    d_minimizers = nullptr;
    d_minimizer_counts = nullptr;
    d_minimizer_offsets = nullptr;
    d_bloom_filter = nullptr;
    d_read_passes_filter = nullptr;
    d_amr_hits = nullptr;
    d_hit_counts = nullptr;
    d_coverage_stats = nullptr;
}

AMRDetectionPipeline::~AMRDetectionPipeline() {
    // Clean up search engine if initialized
    if (translated_search_engine) {
        destroy_translated_search_engine(translated_search_engine);
        translated_search_engine = nullptr;
    }
    freeGPUMemory();
}


bool AMRDetectionPipeline::initialize(const std::string& amr_db_path) {
    std::cout << "Initializing AMR detection pipeline..." << std::endl;
    
    // Load AMR database
    amr_db = std::make_unique<NCBIAMRDatabaseLoader>();
    
    // Parse database path - expecting "dna.fasta,protein.fasta" format
    size_t comma_pos = amr_db_path.find(',');
    if (comma_pos == std::string::npos) {
        std::cerr << "Database path should be 'dna.fasta,protein.fasta'" << std::endl;
        return false;
    }
    
    std::string dna_path = amr_db_path.substr(0, comma_pos);
    std::string protein_path = amr_db_path.substr(comma_pos + 1);
    
    if (!amr_db->loadFromFastaFiles(dna_path, protein_path)) {
        std::cerr << "Failed to load AMR database" << std::endl;
        return false;
    }
    
    amr_db->printDatabaseStats();
    
    
    // Allocate GPU memory
    allocateGPUMemory();
    
    // Build bloom filter only if enabled
    if (config.use_bloom_filter) {
        buildBloomFilter();
    }
    
    // Initialize coverage statistics
    initializeCoverageStats();
    
    // Create translated search engine with Smith-Waterman
    // In AMRDetectionPipeline constructor initialization, ensure consistency:
    
    // Use the SW version as in the resistance pipeline
    translated_search_engine = create_translated_search_engine_with_sw(engine_capacity, true);
    if (!translated_search_engine) {
        std::cerr << "Failed to create translated search engine" << std::endl;
        return false;
    }
    
    // NOTE: The resistance version doesn't have initialize_genetic_code_gpu
    // The genetic code is initialized as a __constant__ in the CUDA file
    
    // Load protein database once
    std::string protein_db_path = config.protein_db_path;
    if (protein_db_path.empty()) {
        protein_db_path = "amr_protein_db";
    }
    
    if (load_protein_database(translated_search_engine, protein_db_path.c_str()) == 0) {
        search_engine_initialized = true;
        std::cout << "Protein database loaded successfully from " << protein_db_path << std::endl;
    } else {
        std::cerr << "Failed to load protein database from " << protein_db_path << std::endl;
        std::cerr << "Please ensure you have run: ./build_amr_protein_db AMRProt.fa " << protein_db_path << std::endl;
        destroy_translated_search_engine(translated_search_engine);
        translated_search_engine = nullptr;
        return false;
    }
    
    // Initialize coverage statistics after database is loaded
    initializeCoverageStats();
    std::cout << "Coverage statistics initialized" << std::endl;
    
    return true;
}

void AMRDetectionPipeline::allocateGPUMemory() {
    // Allocate for maximum batch size
    // For paired-end reads, we process R1 and R2 separately, so need 2x batch size
    // Add 10% safety margin to avoid boundary issues
    size_t max_batch = (config.reads_per_batch * 2) + (config.reads_per_batch / 5);  // 220,000 instead of 200,000
    size_t max_read_len = config.max_read_length;
    
    
    // Read data - no longer merging reads, just process them separately
    size_t max_single_read_len = max_read_len;
    
    cudaError_t err;
    err = cudaMalloc(&d_reads, max_batch * max_single_read_len);
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for reads: " << cudaGetErrorString(err) << std::endl;
        std::cerr << "Requested size: " << (max_batch * max_single_read_len) << " bytes" << std::endl;
        throw std::runtime_error("GPU allocation failed");
    }
    
    err = cudaMalloc(&d_read_offsets, max_batch * sizeof(int));
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for read offsets: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("GPU allocation failed");
    }
    
    err = cudaMalloc(&d_read_lengths, max_batch * sizeof(int));
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for read lengths: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("GPU allocation failed");
    }
    
    err = cudaMalloc(&d_read_ids, max_batch * sizeof(uint32_t));
    if (err != cudaSuccess) {
        std::cerr << "Failed to allocate GPU memory for read IDs: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("GPU allocation failed");
    }
    
    // Minimizers (estimate ~100 minimizers per read)
    size_t max_minimizers = max_batch * 1000;
    cudaMalloc(&d_minimizers, max_minimizers * sizeof(Minimizer));
    cudaMalloc(&d_minimizer_counts, max_batch * sizeof(uint32_t));
    cudaMalloc(&d_minimizer_offsets, (max_batch + 1) * sizeof(uint32_t));  // +1 for offset array
    
    // Bloom filter (only if enabled)
    if (config.use_bloom_filter) {
        size_t bloom_words = (config.bloom_filter_size + 63) / 64;
        cudaMalloc(&d_bloom_filter, bloom_words * sizeof(uint64_t));
        cudaMemset(d_bloom_filter, 0, bloom_words * sizeof(uint64_t));
    }
    
    // Bloom filter results (always allocate for consistency)
    cudaMalloc(&d_read_passes_filter, max_batch * sizeof(bool));
    
    // Hits (max MAX_MATCHES_PER_READ per read)
    cudaMalloc(&d_amr_hits, max_batch * MAX_MATCHES_PER_READ * sizeof(AMRHit));
    cudaMalloc(&d_hit_counts, max_batch * sizeof(uint32_t));
    
    // Coverage statistics
    size_t num_genes = amr_db->getNumGenes();
    cudaMalloc(&d_coverage_stats, num_genes * sizeof(AMRCoverageStats));
    cudaMemset(d_coverage_stats, 0, num_genes * sizeof(AMRCoverageStats));
}

void AMRDetectionPipeline::freeGPUMemory() {
    if (d_reads) cudaFree(d_reads);
    if (d_read_offsets) cudaFree(d_read_offsets);
    if (d_read_lengths) cudaFree(d_read_lengths);
    if (d_read_ids) cudaFree(d_read_ids);
    if (d_minimizers) cudaFree(d_minimizers);
    if (d_minimizer_counts) cudaFree(d_minimizer_counts);
    if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
    if (d_bloom_filter) cudaFree(d_bloom_filter);
    if (d_read_passes_filter) cudaFree(d_read_passes_filter);
    if (d_amr_hits) cudaFree(d_amr_hits);
    if (d_hit_counts) cudaFree(d_hit_counts);
    if (d_coverage_stats) {
        // Free position counts arrays
        freeCoverageStats();
        cudaFree(d_coverage_stats);
    }
}

void AMRDetectionPipeline::initializeCoverageStats() {
    uint32_t num_genes = amr_db->getNumGenes();
    h_coverage_stats.resize(num_genes);
    
    // Get gene information to know lengths
    std::vector<AMRGeneEntry> gene_entries(num_genes);
    cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(),
               num_genes * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    // Calculate total memory needed for position counts
    size_t total_position_memory = 0;
    for (uint32_t i = 0; i < num_genes; i++) {
        h_coverage_stats[i].gene_length = gene_entries[i].protein_length;
        total_position_memory += gene_entries[i].protein_length * sizeof(uint32_t);
    }
    
    // Allocate one big block for all position counts
    uint32_t* d_all_position_counts;
    cudaMalloc(&d_all_position_counts, total_position_memory);
    cudaMemset(d_all_position_counts, 0, total_position_memory);
    
    // Set up pointers for each gene
    size_t offset = 0;
    for (uint32_t i = 0; i < num_genes; i++) {
        h_coverage_stats[i].position_counts = d_all_position_counts + offset;
        offset += gene_entries[i].protein_length;
    }
    
    // Copy coverage stats to GPU
    cudaMemcpy(d_coverage_stats, h_coverage_stats.data(),
               num_genes * sizeof(AMRCoverageStats), cudaMemcpyHostToDevice);
}

void AMRDetectionPipeline::freeCoverageStats() {
    if (!h_coverage_stats.empty() && h_coverage_stats[0].position_counts) {
        // Free the big block of position counts
        cudaFree(h_coverage_stats[0].position_counts);
    }
}

void AMRDetectionPipeline::buildBloomFilter() {
    std::cout << "Building bloom filter from AMR sequences..." << std::endl;
    
    // Get sequence data from database
    char* d_amr_sequences = amr_db->getGPUDNASequences();
    AMRGeneEntry* d_gene_entries = amr_db->getGPUGeneEntries();
    uint32_t num_sequences = amr_db->getNumGenes();
    
    // Create offset and length arrays
    std::vector<uint32_t> offsets(num_sequences);
    std::vector<uint32_t> lengths(num_sequences);
    
    // Copy gene entries to get offsets and lengths
    std::vector<AMRGeneEntry> gene_entries(num_sequences);
    cudaMemcpy(gene_entries.data(), d_gene_entries, 
               num_sequences * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    for (uint32_t i = 0; i < num_sequences; i++) {
        offsets[i] = gene_entries[i].gpu_offset_dna;
        lengths[i] = gene_entries[i].dna_length;
    }
    
    // Copy to GPU
    uint32_t* d_offsets;
    uint32_t* d_lengths;
    cudaMalloc(&d_offsets, num_sequences * sizeof(uint32_t));
    cudaMalloc(&d_lengths, num_sequences * sizeof(uint32_t));
    cudaMemcpy(d_offsets, offsets.data(), num_sequences * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lengths, lengths.data(), num_sequences * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch kernel
    launch_build_bloom_filter_kernel(
        d_amr_sequences,
        d_offsets,
        d_lengths,
        d_bloom_filter,
        num_sequences,
        config.bloom_filter_size,
        config.kmer_length
    );
    
    cudaDeviceSynchronize();
    
    // Cleanup
    cudaFree(d_offsets);
    cudaFree(d_lengths);
    
    std::cout << "Bloom filter built successfully" << std::endl;
}

void AMRDetectionPipeline::copyReadsToGPU(const std::vector<std::string>& reads) {
    
    // Validate input
    if (reads.empty()) {
        std::cerr << "ERROR: reads vector is empty!" << std::endl;
        return;
    }
    
    // Prepare reads in GPU-friendly format
    std::vector<char> concatenated_reads;
    std::vector<int> offsets;
    std::vector<int> lengths;
    
    for (const auto& read : reads) {
        // Skip empty reads
        if (read.empty()) {
            std::cerr << "WARNING: Skipping empty read" << std::endl;
            continue;
        }
        offsets.push_back(concatenated_reads.size());
        lengths.push_back(read.length());
        concatenated_reads.insert(concatenated_reads.end(), read.begin(), read.end());
    }
    
    
    // Validate sizes
    if (offsets.size() != lengths.size()) {
        std::cerr << "ERROR: Size mismatch - offsets: " << offsets.size() 
                  << ", lengths: " << lengths.size() << std::endl;
    }
    
    // Copy to GPU
    cudaMemcpy(d_reads, concatenated_reads.data(), 
               concatenated_reads.size(), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERROR("copying reads to GPU");
    
    cudaMemcpy(d_read_offsets, offsets.data(), 
               offsets.size() * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERROR("copying read offsets to GPU");
    
    cudaMemcpy(d_read_lengths, lengths.data(), 
               lengths.size() * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_CUDA_ERROR("copying read lengths to GPU");
}

void AMRDetectionPipeline::processBatch(const std::vector<std::string>& reads,
                                       const std::vector<std::string>& read_ids) {
    
    // At the start of processBatch(), add:
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "WARNING: Clearing existing CUDA error: " << cudaGetErrorString(err) << std::endl;
        // Clear the error but continue
    }
    
    if (reads.empty()) return;
    
    // Reset current_batch_size
    current_batch_size = reads.size();
    
    total_reads_processed += current_batch_size;
    std::cout << "Processing batch of " << current_batch_size << " reads" << std::endl;
    
    // Performance reporting
    auto current_time = std::chrono::steady_clock::now();
    auto time_since_last_report = std::chrono::duration_cast<std::chrono::seconds>(
        current_time - last_performance_report).count();
    
    if (time_since_last_report >= 60) {  // Report every minute
        auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            current_time - processing_start_time).count();
        uint64_t reads_since_checkpoint = total_reads_processed - reads_processed_checkpoint;
        
        if (time_since_last_report > 0) {
            double reads_per_minute = reads_since_checkpoint * 60.0 / time_since_last_report;
            double reads_per_second = reads_since_checkpoint / static_cast<double>(time_since_last_report);
            
            std::cout << "\n=== Performance Report ===" << std::endl;
            std::cout << "Total reads processed: " << total_reads_processed << std::endl;
            std::cout << "Reads in last period: " << reads_since_checkpoint << std::endl;
            std::cout << "Throughput: " << std::fixed << std::setprecision(1) 
                      << reads_per_minute << " reads/minute (" 
                      << reads_per_second << " reads/second)" << std::endl;
            std::cout << "Total elapsed time: " << total_elapsed << " seconds" << std::endl;
            std::cout << "========================\n" << std::endl;
        }
        
        reads_processed_checkpoint = total_reads_processed;
        last_performance_report = current_time;
    }
    
    // Copy reads to GPU with timing
    auto copy_start = std::chrono::high_resolution_clock::now();
    copyReadsToGPU(reads);
    auto copy_end = std::chrono::high_resolution_clock::now();
    auto copy_duration = std::chrono::duration_cast<std::chrono::milliseconds>(copy_end - copy_start);
    if (copy_duration.count() > 100) {  // Report if takes more than 100ms
        std::cout << "GPU copy time: " << copy_duration.count() << "ms" << std::endl;
    }
    
    // After copyReadsToGPU(), add:
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "ERROR: CUDA error after copyReadsToGPU: " << cudaGetErrorString(err) << std::endl;
        return; // Don't continue with corrupted state
    }
    
    
    // Skip bloom filter completely if disabled (default behavior)
    if (!config.use_bloom_filter) {
        // Don't allocate or use any bloom filter memory
        // Pass nullptr to indicate all reads should be processed
        performTranslatedAlignment();
        
        // Get hits from this batch IMMEDIATELY after alignment
        auto batch_hits = getAMRHits();
        
        // Accumulate hits for EM if enabled
        if (config.use_em && !batch_hits.empty()) {
            accumulated_hits.insert(accumulated_hits.end(), batch_hits.begin(), batch_hits.end());
            std::cout << "Accumulated " << batch_hits.size() << " hits (total: " 
                      << accumulated_hits.size() << ")" << std::endl;
        }
        
        // Update coverage statistics with hits from this batch
        updateCoverageStatsFromHits();
        
        std::cout << "Batch completed successfully" << std::endl;
        return;
    }
    
    // Only do bloom filter operations if explicitly enabled
    generateMinimizers();
    screenWithBloomFilter();
    performTranslatedAlignment();
    
    // Get hits from this batch IMMEDIATELY after alignment
    auto batch_hits = getAMRHits();
    
    // Accumulate hits for EM if enabled
    if (config.use_em && !batch_hits.empty()) {
        accumulated_hits.insert(accumulated_hits.end(), batch_hits.begin(), batch_hits.end());
        std::cout << "Accumulated " << batch_hits.size() << " hits (total: " 
                  << accumulated_hits.size() << ")" << std::endl;
    }
    
    // Update coverage statistics with hits from this batch
    updateCoverageStatsFromHits();
    
    std::cout << "Batch completed successfully" << std::endl;
}

void AMRDetectionPipeline::generateMinimizers() {
    // Calculate minimizer offsets
    std::vector<uint32_t> offsets(current_batch_size + 1);
    offsets[0] = 0;
    for (int i = 1; i <= current_batch_size; i++) {
        offsets[i] = offsets[i-1] + 1000;  // Max minimizers per read
    }
    
    cudaMemcpy(d_minimizer_offsets, offsets.data(), 
               offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch kernel
    launch_generate_minimizers_kernel(
        d_reads,
        d_read_offsets,
        d_read_lengths,
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        current_batch_size,
        config.minimizer_k,
        config.minimizer_w
    );
    
    cudaDeviceSynchronize();
}

void AMRDetectionPipeline::screenWithBloomFilter() {
    // Launch kernel using the class member boolean array
    launch_screen_minimizers_kernel(
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        d_bloom_filter,
        config.bloom_filter_size,
        d_read_passes_filter,
        current_batch_size
    );
    
    cudaDeviceSynchronize();
    
    // Count how many passed
    std::vector<uint8_t> passes_filter(current_batch_size);
    cudaMemcpy(passes_filter.data(), d_read_passes_filter, 
               current_batch_size * sizeof(bool), cudaMemcpyDeviceToHost);
    
    int passed = std::count(passes_filter.begin(), passes_filter.end(), (uint8_t)1);
    std::cout << "Bloom filter: " << passed << "/" << current_batch_size 
              << " reads passed (" << (100.0 * passed / current_batch_size) << "%)" << std::endl;
}

void AMRDetectionPipeline::performTranslatedAlignment() {
    auto align_start = std::chrono::high_resolution_clock::now();
    
    // Check if search engine is initialized
    if (!translated_search_engine || !search_engine_initialized) {
        std::cerr << "Translated search engine not initialized" << std::endl;
        return;
    }
    
    // Add validation
    if (current_batch_size == 0) {
        std::cerr << "ERROR: current_batch_size is 0" << std::endl;
        return;
    }
    
    // Get engine capacity
    if (current_batch_size > engine_capacity) {
        std::cerr << "ERROR: current_batch_size (" << current_batch_size 
                  << ") exceeds engine capacity (" << engine_capacity << ")" << std::endl;
        return;
    }
    
    // Validate concatenated read data
    std::vector<int> h_read_lengths(current_batch_size);
    cudaMemcpy(h_read_lengths.data(), d_read_lengths, 
               current_batch_size * sizeof(int), cudaMemcpyDeviceToHost);
    
    size_t total_bases = 0;
    for (int len : h_read_lengths) {
        if (len <= 0 || len > config.max_read_length) {
            std::cerr << "WARNING: Invalid read length: " << len << std::endl;
        }
        total_bases += len;
    }
    
    
    // Load gene entries from database if not already loaded
    if (gene_entries.empty()) {
        uint32_t num_genes = amr_db->getNumGenes();
        if (num_genes == 0) {
            std::cerr << "No genes loaded in AMR database" << std::endl;
            return;
        }
        
        AMRGeneEntry* gpu_entries = amr_db->getGPUGeneEntries();
        if (!gpu_entries) {
            std::cerr << "GPU gene entries pointer is NULL" << std::endl;
            return;
        }
        
        gene_entries.resize(num_genes);
        
        // Ensure all previous GPU operations are complete
        cudaDeviceSynchronize();
        
        // Check for any existing errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error before gene entries copy: " << cudaGetErrorString(err) << std::endl;
            // Clear the error
            cudaGetLastError();
        }
        
        cudaMemcpy(gene_entries.data(), gpu_entries,
                   num_genes * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
                   
        // Check for errors after copy
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error after gene entries copy: " << cudaGetErrorString(err) << std::endl;
        }
    }
    
    // Allocate space for results
    ProteinMatch* d_protein_matches = nullptr;
    const size_t MAX_MATCHES_PER_READ = 32;  // Must match what's in the kernel
    size_t protein_match_size = current_batch_size * MAX_MATCHES_PER_READ * sizeof(ProteinMatch);
    
    
    cudaError_t alloc_err = cudaMalloc(&d_protein_matches, protein_match_size);
    if (alloc_err != cudaSuccess) {
        std::cerr << "ERROR: Failed to allocate protein matches: " << cudaGetErrorString(alloc_err) << std::endl;
        return;
    }
    
    // Initialize hit counts to zero
    cudaMemset(d_hit_counts, 0, current_batch_size * sizeof(uint32_t));
    
    // Debug: Check all pointers before calling search
    if (!d_reads || !d_read_lengths || !d_read_offsets || !d_read_ids || 
        !d_protein_matches || !d_hit_counts) {
        std::cerr << "ERROR: One or more GPU pointers are NULL:" << std::endl;
        std::cerr << "  d_reads: " << (void*)d_reads << std::endl;
        std::cerr << "  d_read_lengths: " << (void*)d_read_lengths << std::endl;
        std::cerr << "  d_read_offsets: " << (void*)d_read_offsets << std::endl;
        std::cerr << "  d_read_ids: " << (void*)d_read_ids << std::endl;
        std::cerr << "  d_protein_matches: " << (void*)d_protein_matches << std::endl;
        std::cerr << "  d_hit_counts: " << (void*)d_hit_counts << std::endl;
        cudaFree(d_protein_matches);
        return;
    }
    
    // Clear any existing CUDA errors before the call
    cudaGetLastError();
    
    // Bounds checking
    size_t allocated_max = (config.reads_per_batch * 2) + (config.reads_per_batch / 5);
    if (current_batch_size > allocated_max) {
        std::cerr << "ERROR: Batch size " << current_batch_size 
                  << " exceeds allocated maximum " << allocated_max << std::endl;
        return;
    }
    
    
    
    // Use the class member search engine
    // Pass bloom filter results if enabled, otherwise nullptr to process all reads
    bool* reads_filter = config.use_bloom_filter ? d_read_passes_filter : (bool*)nullptr;
    
    
    
    
    // Clear any existing CUDA errors before the call
    cudaError_t pre_err = cudaGetLastError();
    if (pre_err != cudaSuccess) {
        std::cerr << "WARNING: Existing CUDA error before search: " << cudaGetErrorString(pre_err) << std::endl;
    }
    
    
    // Clear any pre-existing CUDA errors before the search
    cudaError_t clear_err = cudaGetLastError();
    if (clear_err != cudaSuccess) {
        std::cerr << "WARNING: Clearing pre-existing CUDA error before search: " 
                  << cudaGetErrorString(clear_err) << std::endl;
    }
    
    // Call search
    int result = search_translated_reads(translated_search_engine, d_reads, d_read_lengths, 
                                       d_read_offsets, reads_filter, 
                                       current_batch_size, d_protein_matches, d_hit_counts);
    
    // Now check for errors
    cudaError_t post_err = cudaGetLastError();
    if (post_err != cudaSuccess) {
        std::cerr << "CUDA error after search_translated_reads: " 
                  << cudaGetErrorString(post_err) << std::endl;
    }
    
    if (result != 0) {
        std::cerr << "search_translated_reads returned error code: " << result << std::endl;
    }
    
    cudaDeviceSynchronize();
    
    // Convert ProteinMatch results to AMRHit format
    std::vector<uint32_t> hit_counts(current_batch_size);
    cudaMemcpy(hit_counts.data(), d_hit_counts, 
               current_batch_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    std::vector<ProteinMatch> protein_matches(current_batch_size * MAX_MATCHES_PER_READ);
    cudaMemcpy(protein_matches.data(), d_protein_matches,
               protein_matches.size() * sizeof(ProteinMatch), cudaMemcpyDeviceToHost);
    
    // Apply paired-end concordance scoring if we have paired read info
    if (!paired_read_info.empty()) {
        applyPairedConcordanceScoring(protein_matches.data(), hit_counts.data(), current_batch_size);
    }
    
    // Convert to AMRHit format and process results with concordance information
    std::vector<AMRHit> amr_hits;
    for (int i = 0; i < current_batch_size; i++) {
        if (hit_counts[i] > 0) {
            // Get paired read info if available
            const PairedReadInfo* pairInfo = nullptr;
            if (!paired_read_info.empty() && i < paired_read_info.size()) {
                pairInfo = &paired_read_info[i];
            }
            
            
            for (uint32_t j = 0; j < hit_counts[i]; j++) {
                ProteinMatch& pm = protein_matches[i * MAX_MATCHES_PER_READ + j];
                AMRHit hit = {};
                
                hit.read_id = pm.read_id;
                hit.gene_id = pm.gene_id;  // Use gene_id, not protein_id
                hit.ref_start = pm.ref_start;
                hit.ref_end = pm.ref_start + pm.match_length;
                hit.read_start = pm.query_start;
                hit.read_end = pm.query_start + pm.match_length;
                hit.identity = pm.identity;
                hit.coverage = (float)pm.match_length / gene_entries[pm.gene_id].protein_length;
                hit.frame = pm.frame;
                hit.is_complete_gene = (pm.ref_start == 0 && 
                                       hit.ref_end >= gene_entries[pm.gene_id].protein_length * 0.95);
                hit.concordant = pm.concordant;  // Copy concordance information
                
                strncpy(hit.gene_name, gene_entries[pm.gene_id].gene_name, 63);
                strncpy(hit.drug_class, gene_entries[pm.gene_id].class_, 31);
                strncpy(hit.gene_family, gene_entries[pm.gene_id].gene_family, 31);
                hit.gene_family[31] = '\0';
                
                // Include concordance information in debug output
                if (pairInfo) {
                    printf("Read %d %s: Gene %d - Score: %.2f %s\n",
                        pairInfo->read_idx,
                        pairInfo->is_read2 ? "R2" : "R1",
                        pm.gene_id,
                        pm.alignment_score,
                        pm.concordant ? "[CONCORDANT]" : "[DISCORDANT]");
                }
                
                amr_hits.push_back(hit);
            }
        }
    }
    
    // Copy converted hits back to GPU in the correct layout
    if (!amr_hits.empty()) {
        // Copy hit counts to GPU
        std::vector<uint32_t> gpu_hit_counts(current_batch_size, 0);
        
        // Count hits per read for GPU storage
        for (const auto& hit : amr_hits) {
            if (hit.read_id < current_batch_size) {
                gpu_hit_counts[hit.read_id]++;
            }
        }
        
        // Copy hits to GPU in the correct layout (grouped by read)
        std::vector<AMRHit> gpu_hits_layout(current_batch_size * MAX_MATCHES_PER_READ);
        std::vector<uint32_t> hit_index(current_batch_size, 0);
        
        for (const auto& hit : amr_hits) {
            uint32_t read_id = hit.read_id;
            uint32_t idx = hit_index[read_id]++;
            if (idx < MAX_MATCHES_PER_READ) {
                gpu_hits_layout[read_id * MAX_MATCHES_PER_READ + idx] = hit;
            }
        }
        
        // Copy to GPU
        cudaMemcpy(d_amr_hits, gpu_hits_layout.data(), 
                   gpu_hits_layout.size() * sizeof(AMRHit), cudaMemcpyHostToDevice);
        cudaMemcpy(d_hit_counts, gpu_hit_counts.data(), 
                   gpu_hit_counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        std::cout << "Copied " << amr_hits.size() << " hits to GPU" << std::endl;
    }
    
    // Cleanup (but don't destroy the search engine - it's reused)
    cudaFree(d_protein_matches);
    
    // Ensure all GPU operations complete and check for errors
    cudaDeviceSynchronize();
    cudaError_t final_err = cudaGetLastError();
    if (final_err != cudaSuccess) {
        std::cerr << "CUDA error after translated search: " << cudaGetErrorString(final_err) << std::endl;
        
        // Try to recover without destroying the engine
        // First, clear the error
        cudaGetLastError();
        
        // If we get repeated errors, consider resetting the engine
        static int consecutive_errors = 0;
        consecutive_errors++;
        
        if (consecutive_errors > 5) {
            std::cerr << "Too many consecutive CUDA errors, resetting translated search engine" << std::endl;
            resetTranslatedSearchEngine();
            consecutive_errors = 0;
        }
    } else {
        // Reset error counter on success
        static int consecutive_errors = 0;
        consecutive_errors = 0;
    }
    
    // Report hits
    int total_hits = 0;
    for (auto count : hit_counts) {
        total_hits += count;
    }
    
    std::cout << "Found " << total_hits << " AMR hits using translated search" << std::endl;
    
    // Report timing
    auto align_end = std::chrono::high_resolution_clock::now();
    auto align_duration = std::chrono::duration_cast<std::chrono::milliseconds>(align_end - align_start);
    if (align_duration.count() > 500) {  // Report if takes more than 500ms
        std::cout << "Translated alignment time: " << align_duration.count() << "ms" << std::endl;
    }
}

void AMRDetectionPipeline::extendAlignments() {
    // Get protein info
    char* d_amr_proteins = amr_db->getGPUProteinSequences();
    uint32_t num_proteins = amr_db->getNumGenes();
    
    // Create protein arrays (reuse from previous step if possible)
    std::vector<AMRGeneEntry> gene_entries(num_proteins);
    AMRGeneEntry* d_gene_entries = amr_db->getGPUGeneEntries();
    cudaMemcpy(gene_entries.data(), d_gene_entries, 
               num_proteins * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    std::vector<uint32_t> protein_offsets(num_proteins);
    std::vector<uint32_t> protein_lengths(num_proteins);
    
    for (uint32_t i = 0; i < num_proteins; i++) {
        protein_offsets[i] = gene_entries[i].gpu_offset_protein;
        protein_lengths[i] = gene_entries[i].protein_length;
    }
    
    uint32_t* d_protein_offsets;
    uint32_t* d_protein_lengths;
    cudaMalloc(&d_protein_offsets, num_proteins * sizeof(uint32_t));
    cudaMalloc(&d_protein_lengths, num_proteins * sizeof(uint32_t));
    cudaMemcpy(d_protein_offsets, protein_offsets.data(), 
               num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_protein_lengths, protein_lengths.data(), 
               num_proteins * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch extension kernel
    launch_extend_alignments_kernel(
        d_reads,
        d_read_offsets,
        d_read_lengths,
        d_minimizers,
        d_minimizer_offsets,
        d_minimizer_counts,
        d_amr_hits,
        d_hit_counts,
        d_amr_proteins,
        d_protein_offsets,
        d_protein_lengths,
        current_batch_size
    );
    
    cudaDeviceSynchronize();
    
    // Cleanup
    cudaFree(d_protein_offsets);
    cudaFree(d_protein_lengths);
}

void AMRDetectionPipeline::updateCoverageStatsFromHits() {
    // Get hits from current batch
    std::vector<uint32_t> hit_counts(current_batch_size);
    cudaMemcpy(hit_counts.data(), d_hit_counts, 
               current_batch_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    // Update stats for each hit
    std::vector<AMRHit> batch_hits(current_batch_size * MAX_MATCHES_PER_READ);
    cudaMemcpy(batch_hits.data(), d_amr_hits, 
               batch_hits.size() * sizeof(AMRHit), cudaMemcpyDeviceToHost);
    
    // Track actual positions covered for each gene using sets
    // Static to persist across batches within a sample
    static std::map<uint32_t, std::set<uint32_t>> gene_covered_positions;
    static bool first_batch = true;
    
    // Clear on first batch of a new sample
    if (first_batch && total_reads_processed == 0) {
        gene_covered_positions.clear();
        first_batch = false;
    }
    
    int hits_processed = 0;
    for (int i = 0; i < current_batch_size; i++) {
        for (uint32_t j = 0; j < hit_counts[i]; j++) {
            AMRHit& hit = batch_hits[i * MAX_MATCHES_PER_READ + j];
            
            if (hit.gene_id < h_coverage_stats.size()) {
                // Update read count and bases mapped
                h_coverage_stats[hit.gene_id].total_reads++;
                h_coverage_stats[hit.gene_id].total_bases_mapped += 
                    (hit.ref_end - hit.ref_start) * 3; // Convert AA to nucleotides
                
                // Track actual positions covered (in amino acid coordinates)
                for (uint16_t pos = hit.ref_start; pos < hit.ref_end; pos++) {
                    gene_covered_positions[hit.gene_id].insert(pos);
                }
                
                // Update covered_positions with actual count of unique positions
                h_coverage_stats[hit.gene_id].covered_positions = 
                    gene_covered_positions[hit.gene_id].size();
                
                hits_processed++;
            }
        }
    }
    
    // Note: The static map will persist across batches within a sample
    // It should be cleared when clearResults() is called between samples
    
    // Reset first_batch flag after processing hits
    if (hits_processed > 0 && total_reads_processed > 0) {
        first_batch = true;  // Reset for next sample
    }
    
    std::cout << "Updated coverage stats with " << hits_processed << " hits from current batch" << std::endl;
}

void AMRDetectionPipeline::finalizeCoverageStats() {
    uint32_t num_genes = amr_db->getNumGenes();
    
    std::cout << "\n=== Finalizing Coverage Stats ===" << std::endl;
    std::cout << "Total reads processed: " << total_reads_processed << std::endl;
    
    // Get gene information
    std::vector<AMRGeneEntry> gene_entries(num_genes);
    cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(),
               num_genes * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    // Calculate final metrics
    int genes_with_coverage = 0;
    double total_rpk = 0.0; // For TPM calculation
    
    // First pass: calculate RPK (reads per kilobase) for each gene
    for (uint32_t i = 0; i < num_genes; i++) {
        auto& stats = h_coverage_stats[i];
        
        if (stats.total_reads > 0) {
            genes_with_coverage++;
            stats.gene_length = gene_entries[i].protein_length;
            
            // Calculate percent coverage (already set by updateCoverageStatsFromHits)
            // covered_positions now contains actual unique positions covered
            if (stats.gene_length > 0) {
                stats.percent_coverage = (float)stats.covered_positions / stats.gene_length * 100.0f;
            }
            
            // Calculate mean depth over covered regions (not entire gene)
            if (stats.covered_positions > 0) {
                // total_bases_mapped is in nucleotides, covered_positions is in amino acids
                // So we need to convert: divide by 3 to get amino acid depth
                stats.mean_depth = (float)stats.total_bases_mapped / (stats.covered_positions * 3);
            } else {
                stats.mean_depth = 0.0f;
            }
            
            // Calculate RPKM and RPK for TPM
            if (total_reads_processed > 0) {
                // Gene length in kilobases (protein length * 3 for nucleotides / 1000)
                double gene_length_kb = (stats.gene_length * 3.0) / 1000.0;
                
                if (gene_length_kb > 0) {
                    // RPKM = (reads * 10^6) / (gene_length_kb * total_reads)
                    stats.rpkm = (stats.total_reads * 1e6) / (gene_length_kb * total_reads_processed);
                    
                    // RPK for TPM calculation
                    double rpk = stats.total_reads / gene_length_kb;
                    total_rpk += rpk;
                }
            }
            
            // Debug first few
            if (genes_with_coverage <= 5) {
                std::cout << "Gene " << i << " (" << gene_entries[i].gene_name << "): "
                          << stats.total_reads << " reads, "
                          << stats.percent_coverage << "% coverage, "
                          << stats.mean_depth << " depth, "
                          << "RPKM=" << stats.rpkm << std::endl;
            }
        }
    }
    
    std::cout << "Genes with coverage: " << genes_with_coverage << std::endl;
    std::cout << "Total RPK sum: " << total_rpk << std::endl;
    
    // REMOVED: resolveAmbiguousAssignmentsEM() call - this should NOT be here!
    // EM should only be called explicitly from main after all batches
    
    // Second pass: calculate TPM
    if (total_rpk > 0) {
        for (uint32_t i = 0; i < num_genes; i++) {
            auto& stats = h_coverage_stats[i];
            
            if (stats.total_reads > 0 && stats.gene_length > 0) {
                double gene_length_kb = (stats.gene_length * 3.0) / 1000.0;
                if (gene_length_kb > 0) {
                    double rpk = stats.total_reads / gene_length_kb;
                    stats.tpm = (rpk / total_rpk) * 1e6;
                    
                    // Debug TPM calculation for first few genes
                    if (i < 5 && stats.tpm > 0) {
                        std::cout << "Gene " << i << " TPM calculation: "
                                  << "reads=" << stats.total_reads 
                                  << ", length_kb=" << gene_length_kb
                                  << ", rpk=" << rpk
                                  << ", tpm=" << stats.tpm << std::endl;
                    }
                }
            }
        }
    }
    
    // Copy final stats to GPU if needed
    cudaMemcpy(d_coverage_stats, h_coverage_stats.data(),
               num_genes * sizeof(AMRCoverageStats), cudaMemcpyHostToDevice);
}

void AMRDetectionPipeline::calculateCoverageStats() {
    // This is now just a wrapper that calls finalizeCoverageStats
    finalizeCoverageStats();
}

std::vector<AMRHit> AMRDetectionPipeline::getAMRHits() {
    std::vector<AMRHit> all_hits;
    std::vector<uint32_t> hit_counts(current_batch_size);
    
    cudaMemcpy(hit_counts.data(), d_hit_counts, 
               current_batch_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    // Count total hits
    int total_hits = 0;
    for (auto count : hit_counts) {
        total_hits += count;
    }
    
    std::cout << "\n=== getAMRHits Debug ===" << std::endl;
    std::cout << "Current batch size: " << current_batch_size << std::endl;
    std::cout << "Total hits in GPU memory: " << total_hits << std::endl;
    
    if (total_hits > 0) {
        // Copy all hits - use MAX_MATCHES_PER_READ instead of hardcoded 10
        std::vector<AMRHit> batch_hits(current_batch_size * MAX_MATCHES_PER_READ);
        cudaMemcpy(batch_hits.data(), d_amr_hits, 
                   batch_hits.size() * sizeof(AMRHit), cudaMemcpyDeviceToHost);
        
        // Extract actual hits and debug first few
        int hits_shown = 0;
        for (int i = 0; i < current_batch_size; i++) {
            for (uint32_t j = 0; j < hit_counts[i]; j++) {
                AMRHit& hit = batch_hits[i * MAX_MATCHES_PER_READ + j];
                all_hits.push_back(hit);
                
                // Debug first few hits
                if (hits_shown < 5) {
                    std::cout << "Hit " << hits_shown << ": Gene=" << hit.gene_name 
                              << " Identity=" << hit.identity 
                              << " Coverage=" << hit.coverage 
                              << " Drug=" << hit.drug_class << std::endl;
                    hits_shown++;
                }
            }
        }
    }
    
    std::cout << "Total hits returned: " << all_hits.size() << std::endl;
    std::cout << "========================\n" << std::endl;
    
    return all_hits;
}

std::vector<AMRCoverageStats> AMRDetectionPipeline::getCoverageStats() {
    uint32_t num_genes = amr_db->getNumGenes();
    std::vector<AMRCoverageStats> stats(num_genes);
    
    cudaMemcpy(stats.data(), d_coverage_stats, 
               num_genes * sizeof(AMRCoverageStats), cudaMemcpyDeviceToHost);
    
    // Debug coverage stats
    int genes_with_coverage = 0;
    for (uint32_t i = 0; i < num_genes; i++) {
        if (stats[i].total_reads > 0) {
            genes_with_coverage++;
            if (genes_with_coverage <= 5) {
                std::cout << "Coverage Stats - Gene " << i << ": reads=" << stats[i].total_reads 
                          << " coverage=" << stats[i].percent_coverage << "%" << std::endl;
            }
        }
    }
    std::cout << "Total genes with coverage > 0: " << genes_with_coverage << " out of " << num_genes << std::endl;
    
    return stats;
}

void AMRDetectionPipeline::writeResults(const std::string& output_prefix) {
    // Write hits to TSV file
    std::string hits_file = output_prefix + "_hits.tsv";
    std::ofstream out(hits_file);
    
    out << "read_id\tgene_name\tdrug_class\tidentity\tcoverage\t"
        << "ref_start\tref_end\tframe\tcomplete_gene\tconcordant\n";
    
    auto hits = getAMRHits();
    for (const auto& hit : hits) {
        out << hit.read_id << "\t"
            << hit.gene_name << "\t"
            << hit.drug_class << "\t"
            << hit.identity << "\t"
            << hit.coverage << "\t"
            << hit.ref_start << "\t"
            << hit.ref_end << "\t"
            << (int)hit.frame << "\t"
            << (hit.is_complete_gene ? "Y" : "N") << "\t"
            << (hit.concordant ? "Y" : "N") << "\n";
    }
    
    out.close();
    
    // Write coverage statistics
    std::string coverage_file = output_prefix + "_coverage.tsv";
    out.open(coverage_file);
    
    out << "gene_id\tgene_name\ttotal_reads\ttotal_bases\tcovered_positions\t"
        << "gene_length\tpercent_coverage\tmean_depth\tRPKM\tTPM\n";
    
    auto coverage_stats = getCoverageStats();
    std::vector<AMRGeneEntry> gene_entries(amr_db->getNumGenes());
    cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(), 
               gene_entries.size() * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    for (size_t i = 0; i < coverage_stats.size(); i++) {
        const auto& stats = coverage_stats[i];
        if (stats.total_reads > 0) {
            out << i << "\t"
                << gene_entries[i].gene_name << "\t"
                << stats.total_reads << "\t"
                << stats.total_bases_mapped << "\t"
                << stats.covered_positions << "\t"
                << stats.gene_length << "\t"
                << std::fixed << std::setprecision(2)
                << stats.percent_coverage << "\t"
                << stats.mean_depth << "\t"
                << stats.rpkm << "\t"
                << stats.tpm << "\n";
        }
    }
    
    out.close();
    
    std::cout << "Results written to " << hits_file << " and " << coverage_file << std::endl;
}

void AMRDetectionPipeline::generateClinicalReport(const std::string& output_file) {
    std::ofstream report(output_file);
    
    report << "=== Clinical Diagnostic Report: Antibiotic Resistance Gene Detection ===\n\n";
    report << "Purpose: Detection of antibiotic resistance genes in patient microbiome sample\n";
    report << "Clinical Application: Guide antibiotic selection based on detected resistance genes\n\n";
    
    // Get coverage statistics and gene information
    auto coverage_stats = getCoverageStats();
    std::vector<AMRGeneEntry> gene_entries(amr_db->getNumGenes());
    cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(),
               gene_entries.size() * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    // Group detected genes by drug class
    std::map<std::string, std::vector<GeneAbundance>> genes_by_class;
    
    for (size_t i = 0; i < coverage_stats.size(); i++) {
        const auto& stats = coverage_stats[i];
        const auto& gene = gene_entries[i];
        
        // Report genes with >90% coverage as "detected"
        if (stats.percent_coverage > 90.0f && stats.total_reads > 0) {
            GeneAbundance abundance;
            abundance.gene_id = i;
            strncpy(abundance.gene_name, gene.gene_name, 63);
            strncpy(abundance.drug_class, gene.class_, 31);
            abundance.read_count = stats.total_reads;
            abundance.rpkm = stats.rpkm;
            abundance.tpm = stats.tpm;
            abundance.coverage_depth = stats.mean_depth;
            abundance.coverage_breadth = stats.percent_coverage;
            
            genes_by_class[gene.class_].push_back(abundance);
        }
    }
    
    // Summary of detected resistance
    report << "=== SUMMARY OF DETECTED RESISTANCE GENES ===\n\n";
    report << "Total antibiotic classes with resistance detected: " << genes_by_class.size() << "\n\n";
    
    // Report by drug class
    for (const auto& [drug_class, genes] : genes_by_class) {
        report << "Antibiotic Class: " << drug_class << "\n";
        report << "Clinical Implication: Resistance detected - consider alternative antibiotics\n";
        report << "Number of resistance genes detected: " << genes.size() << "\n\n";
        
        report << "Detected Genes:\n";
        report << std::left << std::setw(20) << "Gene" 
               << std::setw(15) << "Read Count" 
               << std::setw(15) << "Coverage %" 
               << std::setw(15) << "Depth"
               << std::setw(15) << "TPM" << "\n";
        report << std::string(80, '-') << "\n";
        
        for (const auto& gene : genes) {
            report << std::left << std::setw(20) << gene.gene_name
                   << std::setw(15) << gene.read_count
                   << std::setw(15) << std::fixed << std::setprecision(1) << gene.coverage_breadth
                   << std::setw(15) << std::fixed << std::setprecision(1) << gene.coverage_depth
                   << std::setw(15) << std::fixed << std::setprecision(2) << gene.tpm << "\n";
        }
        report << "\n";
    }
    
    // Clinical recommendations
    report << "=== CLINICAL RECOMMENDATIONS ===\n\n";
    if (genes_by_class.empty()) {
        report << "No resistance genes detected with high confidence (>90% coverage).\n";
        report << "Standard antibiotic therapy may be appropriate.\n";
    } else {
        report << "Resistance genes detected. Consider the following:\n\n";
        
        for (const auto& [drug_class, genes] : genes_by_class) {
            report << "- " << drug_class << " antibiotics: RESISTANCE DETECTED\n";
            report << "  Detected genes: ";
            for (size_t i = 0; i < genes.size(); i++) {
                if (i > 0) report << ", ";
                report << genes[i].gene_name;
            }
            report << "\n  Clinical Action: Avoid " << drug_class << " antibiotics\n\n";
        }
    }
    
    // Technical details
    report << "\n=== TECHNICAL DETAILS ===\n";
    report << "Total reads processed: " << total_reads_processed << "\n";
    report << "Detection threshold: >90% gene coverage\n";
    report << "Note: This analysis detects known resistance genes present in the patient's microbiome\n";
    
    report.close();
    std::cout << "Clinical report written to " << output_file << std::endl;
}

void AMRDetectionPipeline::calculateAbundanceMetrics() {
    if (total_reads_processed == 0) return;
    
    uint32_t num_genes = amr_db->getNumGenes();
    
    // Get current coverage stats from GPU
    cudaMemcpy(h_coverage_stats.data(), d_coverage_stats,
               num_genes * sizeof(AMRCoverageStats), cudaMemcpyDeviceToHost);
    
    // Get gene information
    std::vector<AMRGeneEntry> gene_entries(num_genes);
    cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(),
               num_genes * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    // Calculate total mapped reads for TPM normalization
    uint64_t total_mapped_reads = 0;
    for (uint32_t i = 0; i < num_genes; i++) {
        total_mapped_reads += h_coverage_stats[i].total_reads;
    }
    
    // First pass: Calculate RPKM and the sum for TPM normalization
    double tpm_sum = 0.0;
    for (uint32_t i = 0; i < num_genes; i++) {
        auto& stats = h_coverage_stats[i];
        
        // Set gene length from database
        stats.gene_length = gene_entries[i].protein_length;
        
        if (stats.gene_length > 0) {
            // Calculate percent coverage
            stats.percent_coverage = (float)stats.covered_positions / stats.gene_length * 100.0f;
            
            // Calculate mean depth over covered regions
            if (stats.covered_positions > 0) {
                // total_bases_mapped is in nucleotides, covered_positions is in amino acids
                // So we need to convert: divide by 3 to get amino acid depth
                stats.mean_depth = (float)stats.total_bases_mapped / (stats.covered_positions * 3);
            } else {
                stats.mean_depth = 0.0f;
            }
            
            // Calculate RPKM: (reads * 10^9) / (gene_length * total_reads)
            // Note: gene_length is in amino acids, so multiply by 3 for nucleotides
            double gene_length_kb = (stats.gene_length * 3) / 1000.0;
            double total_reads_millions = total_reads_processed / 1000000.0;
            
            if (gene_length_kb > 0 && total_reads_millions > 0) {
                stats.rpkm = (float)(stats.total_reads / (gene_length_kb * total_reads_millions));
            } else {
                stats.rpkm = 0.0f;
            }
            
            // Calculate reads per kilobase for TPM
            double rpk = stats.total_reads / gene_length_kb;
            tpm_sum += rpk;
        }
    }
    
    // Second pass: Calculate TPM
    if (tpm_sum > 0) {
        for (uint32_t i = 0; i < num_genes; i++) {
            auto& stats = h_coverage_stats[i];
            if (stats.gene_length > 0) {
                double gene_length_kb = (stats.gene_length * 3) / 1000.0;
                double rpk = stats.total_reads / gene_length_kb;
                stats.tpm = (float)(rpk / tpm_sum * 1000000.0);
            } else {
                stats.tpm = 0.0f;
            }
        }
    }
    
    // Copy updated stats back to GPU
    cudaMemcpy(d_coverage_stats, h_coverage_stats.data(),
               num_genes * sizeof(AMRCoverageStats), cudaMemcpyHostToDevice);
}

void AMRDetectionPipeline::exportAbundanceTable(const std::string& output_file) {
    std::ofstream out(output_file);
    
    // Write header
    out << "gene_name\tdrug_class\tread_count\trpkm\ttpm\tcoverage_depth\tcoverage_breadth\n";
    
    // Get coverage statistics
    auto coverage_stats = getCoverageStats();
    
    // Get gene information
    std::vector<AMRGeneEntry> gene_entries(amr_db->getNumGenes());
    cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(),
               gene_entries.size() * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
    // Export abundance data for each gene
    for (size_t i = 0; i < coverage_stats.size(); i++) {
        const auto& stats = coverage_stats[i];
        const auto& gene = gene_entries[i];
        
        // Only export genes with detected reads
        if (stats.total_reads > 0) {
            GeneAbundance abundance;
            abundance.gene_id = i;
            strncpy(abundance.gene_name, gene.gene_name, 63);
            strncpy(abundance.drug_class, gene.class_, 31);
            abundance.read_count = stats.total_reads;
            abundance.rpkm = stats.rpkm;
            abundance.tpm = stats.tpm;
            abundance.coverage_depth = stats.mean_depth;
            abundance.coverage_breadth = stats.percent_coverage;
            
            // Write to file
            out << abundance.gene_name << "\t"
                << abundance.drug_class << "\t"
                << abundance.read_count << "\t"
                << std::fixed << std::setprecision(2)
                << abundance.rpkm << "\t"
                << abundance.tpm << "\t"
                << abundance.coverage_depth << "\t"
                << abundance.coverage_breadth << "\n";
        }
    }
    
    out.close();
    std::cout << "Abundance table written to " << output_file << std::endl;
}

void AMRDetectionPipeline::resetTranslatedSearchEngine() {
    std::cout << "Resetting translated search engine..." << std::endl;
    
    // Add more thorough cleanup
    if (translated_search_engine) {
        // Ensure all operations complete
        cudaDeviceSynchronize();
        
        // Clear any errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cout << "Clearing CUDA error before engine reset: " << cudaGetErrorString(err) << std::endl;
        }
        
        destroy_translated_search_engine(translated_search_engine);
        translated_search_engine = nullptr;
        search_engine_initialized = false;
    }
    
    // Wait a moment for GPU to settle
    cudaDeviceSynchronize();
    
    // Create new engine with proper size
    int engine_batch_size = (config.reads_per_batch * 2) + (config.reads_per_batch / 5);
    std::cout << "Creating new translated search engine with batch size: " 
              << engine_batch_size << std::endl;
    
    translated_search_engine = create_translated_search_engine_with_sw(engine_batch_size, true);
    if (!translated_search_engine) {
        std::cerr << "Failed to create new translated search engine during reset" << std::endl;
        return;
    }
    
    // Reload protein database
    std::string protein_db_path = config.protein_db_path;
    if (protein_db_path.empty()) {
        protein_db_path = "amr_protein_db";
    }
    
    if (load_protein_database(translated_search_engine, protein_db_path.c_str()) == 0) {
        search_engine_initialized = true;
        std::cout << "Translated search engine reset successfully" << std::endl;
    } else {
        std::cerr << "Failed to reload protein database after reset" << std::endl;
        destroy_translated_search_engine(translated_search_engine);
        translated_search_engine = nullptr;
    }
}

void AMRDetectionPipeline::processBatchPaired(const std::vector<std::string>& reads1,
                                              const std::vector<std::string>& reads2,
                                              const std::vector<std::string>& read_ids) {
    
    // Enhanced paired-end processing strategy:
    // 1. Process R1 and R2 reads maintaining their pair relationship
    // 2. During alignment phase:
    //    - Search both R1 and R2 against the database
    //    - Track which genes each read maps to
    //    - Identify concordant pairs (both reads map to same gene)
    // 3. Scoring adjustments:
    //    - Concordant pairs: Boost alignment scores significantly (2.0x)
    //    - Discordant pairs: Penalize scores (0.5x)
    //    - Use fragment length distribution for additional scoring
    // 4. For EM algorithm:
    //    - Treat read pairs as single units for assignment
    //    - Prefer assigning pairs to genes where both reads map
    // 5. Benefits:
    //    - Increased specificity for closely related genes
    //    - Reduced false positives from spurious single-read alignments
    //    - Better resolution of multi-mapping reads
    
    // Clear paired read info from previous batch
    paired_read_info.clear();
    
    // Reset batch size 
    current_batch_size = 0;
    
    
    if (reads1.empty() && reads2.empty()) return;
    
    // Update total reads processed for paired-end
    total_reads_processed += std::min(reads1.size(), reads2.size());
    
    // Performance reporting
    auto current_time = std::chrono::steady_clock::now();
    auto time_since_last_report = std::chrono::duration_cast<std::chrono::seconds>(
        current_time - last_performance_report).count();
    
    if (time_since_last_report >= 60) {  // Report every minute
        auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            current_time - processing_start_time).count();
        uint64_t reads_since_checkpoint = total_reads_processed - reads_processed_checkpoint;
        
        if (time_since_last_report > 0) {
            double reads_per_minute = reads_since_checkpoint * 60.0 / time_since_last_report;
            double reads_per_second = reads_since_checkpoint / static_cast<double>(time_since_last_report);
            
            std::cout << "\n=== Performance Report (Paired-End) ===" << std::endl;
            std::cout << "Total read pairs processed: " << total_reads_processed << std::endl;
            std::cout << "Read pairs in last period: " << reads_since_checkpoint << std::endl;
            std::cout << "Throughput: " << std::fixed << std::setprecision(1) 
                      << reads_per_minute << " read pairs/minute (" 
                      << reads_per_second << " read pairs/second)" << std::endl;
            std::cout << "Total elapsed time: " << total_elapsed << " seconds" << std::endl;
            std::cout << "=====================================\n" << std::endl;
        }
        
        reads_processed_checkpoint = total_reads_processed;
        last_performance_report = current_time;
    }
    
    // Pre-allocate for 2x reads (R1 and R2 separate)
    int num_pairs = std::min(reads1.size(), reads2.size());
    int max_reads = num_pairs * 2;
    paired_read_info.clear();
    paired_read_info.reserve(max_reads);
    
    // Prepare reads for GPU processing
    std::vector<std::string> all_reads;
    std::vector<std::string> all_read_ids;
    all_reads.reserve(max_reads);
    all_read_ids.reserve(max_reads);
    
    // Process each read pair
    for (int i = 0; i < num_pairs; i++) {
        // Process R1
        if (!reads1[i].empty()) {
            all_reads.push_back(reads1[i]);
            all_read_ids.push_back(read_ids[i] + "_R1");
            
            uint32_t r1_idx = all_reads.size() - 1;
            uint32_t r2_idx = r1_idx + 1; // Will be the R2 index
            
            paired_read_info.push_back({
                static_cast<uint32_t>(i),     // Original pair index
                false,                        // This is R1
                r2_idx                       // Points to R2
            });
        }
        
        // Process R2
        if (!reads2[i].empty()) {
            all_reads.push_back(reads2[i]);
            all_read_ids.push_back(read_ids[i] + "_R2");
            
            uint32_t r2_idx = all_reads.size() - 1;
            uint32_t r1_idx = r2_idx - 1; // Points back to R1
            
            paired_read_info.push_back({
                static_cast<uint32_t>(i),     // Original pair index
                true,                         // This is R2
                r1_idx                       // Points to R1
            });
        }
    }
    
    
    // Process all reads together (R1 and R2 separately)
    // Ensure current_batch_size doesn't exceed engine capacity
    int engine_batch_size = (config.reads_per_batch * 2) + (config.reads_per_batch / 5);  // Match allocation
    
    if (all_reads.size() > engine_batch_size) {
        // Process in smaller chunks
        std::cout << "Batch size " << all_reads.size() << " exceeds engine capacity " 
                  << engine_batch_size << ", processing in chunks" << std::endl;
        
        size_t chunk_size = engine_batch_size;
        for (size_t i = 0; i < all_reads.size(); i += chunk_size) {
            size_t end = std::min(i + chunk_size, all_reads.size());
            std::vector<std::string> chunk_reads(all_reads.begin() + i, all_reads.begin() + end);
            std::vector<std::string> chunk_ids(all_read_ids.begin() + i, all_read_ids.begin() + end);
            
            std::cout << "Processing chunk " << (i/chunk_size + 1) << ": reads " << i << "-" << end << std::endl;
            current_batch_size = chunk_reads.size();
            processBatch(chunk_reads, chunk_ids);
        }
        return;
    }
    
    // Process normally if within limits
    current_batch_size = all_reads.size();
    if (current_batch_size > 0) {
        processBatch(all_reads, all_read_ids);
        
        // After processing, build pair alignment information
        current_pair_alignments.clear();
        current_pair_alignments.reserve(num_pairs);

        // Get the hits from this batch
        auto batch_hits = getAMRHits();

        // Build alignment information for each pair
        for (int pair_idx = 0; pair_idx < num_pairs; pair_idx++) {
            ReadPairAlignment pair_align;
            pair_align.pair_id = pair_idx;
            pair_align.r1_read_idx = pair_idx * 2;
            pair_align.r2_read_idx = pair_idx * 2 + 1;
            
            // Collect R1 and R2 hits for this pair
            for (const auto& hit : batch_hits) {
                if (hit.read_id == pair_align.r1_read_idx) {
                    pair_align.r1_gene_ids.push_back(hit.gene_id);
                    pair_align.r1_scores.push_back(hit.identity * (hit.ref_end - hit.ref_start));
                } else if (hit.read_id == pair_align.r2_read_idx) {
                    pair_align.r2_gene_ids.push_back(hit.gene_id);
                    pair_align.r2_scores.push_back(hit.identity * (hit.ref_end - hit.ref_start));
                }
            }
            
            // Find concordant genes (present in both R1 and R2)
            for (size_t i = 0; i < pair_align.r1_gene_ids.size(); i++) {
                for (size_t j = 0; j < pair_align.r2_gene_ids.size(); j++) {
                    if (pair_align.r1_gene_ids[i] == pair_align.r2_gene_ids[j]) {
                        pair_align.concordant_genes.push_back(pair_align.r1_gene_ids[i]);
                        pair_align.concordant_scores.push_back(
                            pair_align.r1_scores[i] + pair_align.r2_scores[j]
                        );
                    }
                }
            }
            
            if (pair_align.hasAlignment()) {
                current_pair_alignments.push_back(pair_align);
            }
        }

        std::cout << "Pair alignment summary:" << std::endl;
        std::cout << "  Total pairs: " << num_pairs << std::endl;
        std::cout << "  Pairs with alignments: " << current_pair_alignments.size() << std::endl;
        int concordant_count = 0;
        for (const auto& pa : current_pair_alignments) {
            if (pa.isConcordant()) concordant_count++;
        }
        std::cout << "  Concordant pairs: " << concordant_count << std::endl;
        
        // Apply integrated paired-end scoring
        integratePairedEndScoring();
    }
}

void AMRDetectionPipeline::processPairedBatch(const std::vector<ReadPairData>& read_pairs) {
    if (read_pairs.empty()) {
        return;
    }
    
    // Track read processing performance
    total_reads_processed += read_pairs.size();
    
    // Performance reporting
    auto current_time = std::chrono::steady_clock::now();
    double time_since_last_report = std::chrono::duration<double>(
        current_time - last_performance_report).count();
    
    if (time_since_last_report >= 60.0) {  // Report every minute
        double total_elapsed = std::chrono::duration<double>(
            current_time - processing_start_time).count();
        uint64_t reads_since_checkpoint = total_reads_processed - reads_processed_checkpoint;
        
        if (time_since_last_report > 0) {
            double reads_per_minute = reads_since_checkpoint * 60.0 / time_since_last_report;
            double reads_per_second = reads_since_checkpoint / static_cast<double>(time_since_last_report);
            
            std::cout << "\n=== Performance Report (Paired-End Improved) ===" << std::endl;
            std::cout << "Total read pairs processed: " << total_reads_processed << std::endl;
            std::cout << "Read pairs in last period: " << reads_since_checkpoint << std::endl;
            std::cout << "Throughput: " << std::fixed << std::setprecision(1) 
                      << reads_per_minute << " read pairs/minute (" 
                      << reads_per_second << " read pairs/second)" << std::endl;
            std::cout << "Total elapsed time: " << total_elapsed << " seconds" << std::endl;
            std::cout << "=====================================\n" << std::endl;
        }
        
        reads_processed_checkpoint = total_reads_processed;
        last_performance_report = current_time;
    }
    
    // Clear previous batch data
    paired_read_info.clear();
    current_pair_alignments.clear();
    
    // Prepare separate reads for processing
    std::vector<std::string> all_reads;
    std::vector<std::string> all_read_ids;
    all_reads.reserve(read_pairs.size() * 2);
    all_read_ids.reserve(read_pairs.size() * 2);
    paired_read_info.reserve(read_pairs.size() * 2);
    
    // Process each read pair, keeping R1 and R2 separate
    for (size_t i = 0; i < read_pairs.size(); i++) {
        const auto& pair = read_pairs[i];
        
        // Add R1
        if (!pair.read1_seq.empty()) {
            all_reads.push_back(pair.read1_seq);
            all_read_ids.push_back(pair.read1_id);
            
            uint32_t r1_idx = all_reads.size() - 1;
            paired_read_info.push_back({
                pair.pair_index,              // Original pair index
                false,                        // This is R1
                r1_idx + 1                   // Points to R2 (will be next)
            });
        }
        
        // Add R2
        if (!pair.read2_seq.empty()) {
            all_reads.push_back(pair.read2_seq);
            all_read_ids.push_back(pair.read2_id);
            
            uint32_t r2_idx = all_reads.size() - 1;
            paired_read_info.push_back({
                pair.pair_index,              // Original pair index
                true,                         // This is R2
                r2_idx - 1                   // Points to R1 (previous)
            });
        }
    }
    
    // Check batch size limits
    int engine_batch_size = (config.reads_per_batch * 2) + (config.reads_per_batch / 5);
    
    if (all_reads.size() > engine_batch_size) {
        // Process in smaller chunks
        std::cout << "Batch size " << all_reads.size() << " exceeds engine capacity " 
                  << engine_batch_size << ", processing in chunks" << std::endl;
        
        size_t chunk_size = engine_batch_size;
        for (size_t i = 0; i < all_reads.size(); i += chunk_size) {
            size_t end = std::min(i + chunk_size, all_reads.size());
            std::vector<std::string> chunk_reads(all_reads.begin() + i, all_reads.begin() + end);
            std::vector<std::string> chunk_ids(all_read_ids.begin() + i, all_read_ids.begin() + end);
            
            std::cout << "Processing chunk " << (i/chunk_size + 1) << ": reads " << i << "-" << end << std::endl;
            current_batch_size = chunk_reads.size();
            processBatch(chunk_reads, chunk_ids);
        }
        return;
    }
    
    // Process all reads in single batch
    current_batch_size = all_reads.size();
    if (current_batch_size > 0) {
        processBatch(all_reads, all_read_ids);
        
        // After processing, build enhanced pair alignment information
        auto batch_hits = getAMRHits();
        
        // Update hits with paired-end information
        updateHitsWithPairedInfo(batch_hits, read_pairs);
        
        buildEnhancedPairAlignments(batch_hits, read_pairs.size());
        
        // Estimate fragment length distribution from concordant pairs
        if (!fragment_dist_estimated && current_pair_alignments.size() >= 100) {
            estimateFragmentLengthDistribution(current_pair_alignments);
        }
        
        // Update alignment scores based on fragment length probabilities
        if (fragment_dist_estimated) {
            updatePairAlignmentScores(current_pair_alignments);
        }
        
        // Apply integrated paired-end scoring with fragment length info
        integratePairedEndScoring();
        
        // Store the updated hits back to GPU for further processing
        updateGPUHitsWithPairedInfo(batch_hits);
    }
}

void AMRDetectionPipeline::buildEnhancedPairAlignments(
    const std::vector<AMRHit>& batch_hits, 
    size_t num_pairs
) {
    current_pair_alignments.clear();
    current_pair_alignments.reserve(num_pairs);
    
    // Build alignment information for each pair
    for (size_t pair_idx = 0; pair_idx < num_pairs; pair_idx++) {
        ReadPairAlignment pair_align;
        pair_align.pair_id = pair_idx;
        pair_align.r1_read_idx = pair_idx * 2;
        pair_align.r2_read_idx = pair_idx * 2 + 1;
        
        // Collect R1 and R2 hits with position information
        for (const auto& hit : batch_hits) {
            if (hit.read_id == pair_align.r1_read_idx) {
                pair_align.r1_gene_ids.push_back(hit.gene_id);
                pair_align.r1_scores.push_back(hit.identity * (hit.ref_end - hit.ref_start));
                pair_align.r1_positions.push_back(hit.ref_start);
                pair_align.r1_frames.push_back(hit.frame);
            } else if (hit.read_id == pair_align.r2_read_idx) {
                pair_align.r2_gene_ids.push_back(hit.gene_id);
                pair_align.r2_scores.push_back(hit.identity * (hit.ref_end - hit.ref_start));
                pair_align.r2_positions.push_back(hit.ref_start);
                pair_align.r2_frames.push_back(hit.frame);
            }
        }
        
        // Find concordant genes and calculate fragment lengths
        for (size_t i = 0; i < pair_align.r1_gene_ids.size(); i++) {
            for (size_t j = 0; j < pair_align.r2_gene_ids.size(); j++) {
                if (pair_align.r1_gene_ids[i] == pair_align.r2_gene_ids[j]) {
                    uint32_t gene_id = pair_align.r1_gene_ids[i];
                    
                    // Get gene length for fragment calculation
                    AMRGeneEntry gene_entry = getGeneEntry(gene_id);
                    int32_t gene_length = gene_entry.protein_length;
                    
                    // Calculate fragment length
                    int32_t frag_len = pair_align.calculateFragmentLength(
                        pair_align.r1_positions[i],
                        pair_align.r2_positions[j],
                        gene_length,
                        150,  // Assume 150bp reads (can be made configurable)
                        150
                    );
                    
                    if (frag_len > 0) {
                        pair_align.concordant_genes.push_back(gene_id);
                        pair_align.concordant_scores.push_back(
                            pair_align.r1_scores[i] + pair_align.r2_scores[j]
                        );
                        pair_align.fragment_lengths.push_back(frag_len);
                        
                        // Calculate fragment probability if distribution is known
                        float frag_prob = fragment_dist_estimated ? 
                            calculateFragmentLengthProbability(frag_len) : 1.0f;
                        pair_align.fragment_probabilities.push_back(frag_prob);
                    }
                }
            }
        }
        
        // Find best concordant match
        if (pair_align.isConcordant()) {
            float best_score = 0.0f;
            for (size_t i = 0; i < pair_align.concordant_genes.size(); i++) {
                float combined_score = pair_align.concordant_scores[i] * 
                                     pair_align.fragment_probabilities[i];
                if (combined_score > best_score) {
                    best_score = combined_score;
                    pair_align.best_gene_id = pair_align.concordant_genes[i];
                    pair_align.best_fragment_length = pair_align.fragment_lengths[i];
                    pair_align.best_concordant_score = combined_score;
                }
            }
        }
        
        if (pair_align.hasAlignment()) {
            current_pair_alignments.push_back(pair_align);
        }
    }
    
    // Report alignment statistics
    std::cout << "Enhanced pair alignment summary:" << std::endl;
    std::cout << "  Total pairs: " << num_pairs << std::endl;
    std::cout << "  Pairs with alignments: " << current_pair_alignments.size() << std::endl;
    
    int concordant_count = 0;
    float avg_fragment_length = 0.0f;
    int valid_fragments = 0;
    
    for (const auto& pa : current_pair_alignments) {
        if (pa.isConcordant()) {
            concordant_count++;
            if (pa.best_fragment_length > 0) {
                avg_fragment_length += pa.best_fragment_length;
                valid_fragments++;
            }
        }
    }
    
    std::cout << "  Concordant pairs: " << concordant_count << std::endl;
    if (valid_fragments > 0) {
        avg_fragment_length /= valid_fragments;
        std::cout << "  Average fragment length: " << avg_fragment_length << " bp" << std::endl;
    }
}

void AMRDetectionPipeline::estimateFragmentLengthDistribution(
    const std::vector<ReadPairAlignment>& alignments
) {
    std::vector<float> fragment_lengths;
    
    // Collect fragment lengths from high-confidence concordant pairs
    for (const auto& align : alignments) {
        if (align.isConcordant() && align.best_fragment_length > 0) {
            // Only use high-scoring pairs for distribution estimation
            if (align.best_concordant_score > 100.0f) {
                fragment_lengths.push_back(static_cast<float>(align.best_fragment_length));
            }
        }
    }
    
    if (fragment_lengths.size() < 50) {
        std::cout << "Not enough concordant pairs for fragment length estimation" << std::endl;
        return;
    }
    
    // Calculate mean
    float sum = 0.0f;
    for (float len : fragment_lengths) {
        sum += len;
    }
    estimated_fragment_mean = sum / fragment_lengths.size();
    
    // Calculate standard deviation
    float sq_sum = 0.0f;
    for (float len : fragment_lengths) {
        float diff = len - estimated_fragment_mean;
        sq_sum += diff * diff;
    }
    estimated_fragment_std = std::sqrt(sq_sum / fragment_lengths.size());
    
    // Remove outliers and recalculate
    std::vector<float> filtered_lengths;
    for (float len : fragment_lengths) {
        float z_score = std::abs(len - estimated_fragment_mean) / estimated_fragment_std;
        if (z_score <= 2.0f) {  // Keep within 2 standard deviations
            filtered_lengths.push_back(len);
        }
    }
    
    if (filtered_lengths.size() >= 30) {
        // Recalculate with filtered data
        sum = 0.0f;
        for (float len : filtered_lengths) {
            sum += len;
        }
        estimated_fragment_mean = sum / filtered_lengths.size();
        
        sq_sum = 0.0f;
        for (float len : filtered_lengths) {
            float diff = len - estimated_fragment_mean;
            sq_sum += diff * diff;
        }
        estimated_fragment_std = std::sqrt(sq_sum / filtered_lengths.size());
    }
    
    fragment_dist_estimated = true;
    
    std::cout << "Fragment length distribution estimated:" << std::endl;
    std::cout << "  Mean: " << estimated_fragment_mean << " bp" << std::endl;
    std::cout << "  Std dev: " << estimated_fragment_std << " bp" << std::endl;
    std::cout << "  Based on " << filtered_lengths.size() << " concordant pairs" << std::endl;
}

float AMRDetectionPipeline::calculateFragmentLengthProbability(int32_t fragment_length) const {
    if (!fragment_dist_estimated || fragment_length <= 0) {
        return 1.0f;  // Neutral probability if no distribution
    }
    
    // Calculate z-score
    float z = (fragment_length - estimated_fragment_mean) / estimated_fragment_std;
    
    // Probability from normal distribution
    // Using simplified formula: exp(-0.5 * z^2)
    float prob = std::exp(-0.5f * z * z);
    
    // Clamp to reasonable range
    return std::max(0.01f, std::min(1.0f, prob));
}

void AMRDetectionPipeline::updatePairAlignmentScores(
    std::vector<ReadPairAlignment>& alignments
) {
    for (auto& align : alignments) {
        // Recalculate fragment probabilities with estimated distribution
        for (size_t i = 0; i < align.fragment_lengths.size(); i++) {
            align.fragment_probabilities[i] = 
                calculateFragmentLengthProbability(align.fragment_lengths[i]);
        }
        
        // Update best concordant match
        if (align.isConcordant()) {
            float best_score = 0.0f;
            for (size_t i = 0; i < align.concordant_genes.size(); i++) {
                float combined_score = align.concordant_scores[i] * 
                                     align.fragment_probabilities[i] *
                                     config.concordance_bonus;  // Apply concordance bonus
                                     
                if (combined_score > best_score) {
                    best_score = combined_score;
                    align.best_gene_id = align.concordant_genes[i];
                    align.best_fragment_length = align.fragment_lengths[i];
                    align.best_concordant_score = combined_score;
                }
            }
        }
    }
}

void AMRDetectionPipeline::applyPairedConcordanceScoring(
    void* matches_ptr,
    uint32_t* match_counts,
    int num_reads
) {
    // Cast void pointer to ProteinMatch*
    ProteinMatch* matches = static_cast<ProteinMatch*>(matches_ptr);
    
    // Only apply concordance scoring if we have paired read info
    if (paired_read_info.empty()) {
        return;
    }
    
    // Create a map of read pairs to their matches
    struct GeneAssignment {
        uint32_t gene_id;
        float score;
        int read_idx;
    };
    
    std::map<uint32_t, std::vector<GeneAssignment>> pairAssignments;
    
    // Collect all assignments by pair
    for (int i = 0; i < num_reads; i++) {
        if (match_counts[i] == 0) continue;
        
        // Make sure we have paired info for this read
        if (i >= paired_read_info.size()) continue;
        
        uint32_t pair_idx = paired_read_info[i].read_idx;
        
        for (uint32_t j = 0; j < match_counts[i]; j++) {
            ProteinMatch& match = matches[i * MAX_MATCHES_PER_READ + j];
            pairAssignments[pair_idx].push_back({
                match.gene_id,
                match.alignment_score,
                i
            });
        }
    }
    
    // Apply concordance bonus (increased stringency)
    const float CONCORDANCE_BONUS = 2.0f;  // 100% score boost for concordant pairs (increased from 1.5f)
    const float DISCORD_PENALTY = 0.5f;    // 50% penalty for discordant pairs (decreased from 0.8f)
    
    for (auto& pair : pairAssignments) {
        auto& assignments = pair.second;
        
        // Remove species_id from the pairing logic
        // Just use gene_id for concordance checking
        std::map<uint32_t, std::vector<int>> genePairs;
        
        for (size_t i = 0; i < assignments.size(); i++) {
            genePairs[assignments[i].gene_id].push_back(i);
        }
        
        // Apply bonuses/penalties
        for (auto& g_pair : genePairs) {
            if (g_pair.second.size() >= 2) {
                // Concordant: both reads map to same gene
                for (int idx : g_pair.second) {
                    int read_idx = assignments[idx].read_idx;
                    
                    // Find the match in the results
                    for (uint32_t j = 0; j < match_counts[read_idx]; j++) {
                        ProteinMatch& match = matches[read_idx * MAX_MATCHES_PER_READ + j];
                        if (match.gene_id == assignments[idx].gene_id) {
                            match.alignment_score *= CONCORDANCE_BONUS;
                            match.concordant = true;
                            break;
                        }
                    }
                }
            } else {
                // Discordant: only one read maps to this gene
                for (int idx : g_pair.second) {
                    int read_idx = assignments[idx].read_idx;
                    for (uint32_t j = 0; j < match_counts[read_idx]; j++) {
                        ProteinMatch& match = matches[read_idx * MAX_MATCHES_PER_READ + j];
                        if (match.gene_id == assignments[idx].gene_id) {
                            match.alignment_score *= DISCORD_PENALTY;
                            match.concordant = false;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Re-sort matches by score within each read after applying bonuses
    for (int i = 0; i < num_reads; i++) {
        if (match_counts[i] > 1) {
            std::sort(
                matches + i * MAX_MATCHES_PER_READ,
                matches + i * MAX_MATCHES_PER_READ + match_counts[i],
                [](const ProteinMatch& a, const ProteinMatch& b) {
                    return a.alignment_score > b.alignment_score;
                }
            );
        }
    }
}

// Add these implementations
void AMRDetectionPipeline::setAccumulatedHits(const std::vector<AMRHit>& hits) {
    accumulated_hits = hits;
    std::cout << "Set " << accumulated_hits.size() << " accumulated hits for EM processing" << std::endl;
}

std::vector<AMRHit> AMRDetectionPipeline::getAllAccumulatedHits() const {
    return accumulated_hits;
}

void AMRDetectionPipeline::addBatchHits(const std::vector<AMRHit>& batch_hits) {
    accumulated_hits.insert(accumulated_hits.end(), batch_hits.begin(), batch_hits.end());
    std::cout << "Added " << batch_hits.size() << " hits to accumulated total (now " 
              << accumulated_hits.size() << " total)" << std::endl;
}

void AMRDetectionPipeline::resolveAmbiguousAssignmentsEM() {
    std::cout << "=== EM ENTRY POINT REACHED ===" << std::endl;
    std::cout << "\n=== STARTING EM ALGORITHM ===" << std::endl;
    std::cout << "EM enabled: " << (config.use_em ? "YES" : "NO") << std::endl;
    std::cout << "Accumulated hits: " << accumulated_hits.size() << std::endl;
    
    if (!config.use_em) {
        std::cout << "EM algorithm is disabled - skipping" << std::endl;
        return;
    }
    
    std::cout << "\n=== Resolving Ambiguous Read Assignments with Kallisto-style EM ===" << std::endl;
    
    // This method should only be called when EM is enabled and we have accumulated hits
    if (accumulated_hits.empty()) {
        std::cout << "No accumulated hits to resolve" << std::endl;
        return;
    }
    
    std::cout << "Input: " << accumulated_hits.size() << " total accumulated hits" << std::endl;
    
    // Check if we have paired-end data
    bool has_paired_data = false;
    for (const auto& hit : accumulated_hits) {
        if (hit.pair_id != UINT32_MAX && hit.mate_read_id != UINT32_MAX) {
            has_paired_data = true;
            break;
        }
    }
    
    if (has_paired_data) {
        std::cout << "Detected paired-end data - using paired-end aware EM" << std::endl;
        
        // Build paired read assignments
        buildPairedReadAssignments(accumulated_hits);
        
        if (paired_read_assignments.empty()) {
            std::cout << "No paired assignments to resolve" << std::endl;
            return;
        }
        
        // Run paired-end aware EM
        runPairedEndEM();
    } else {
        // Fall back to single-end EM
        std::cout << "Using single-end EM algorithm" << std::endl;
        
        // Build read assignments with quality filtering using accumulated hits
        buildReadAssignments(accumulated_hits);
        
        if (read_assignments.empty()) {
            std::cout << "No high-quality assignments to resolve" << std::endl;
            return;
        }
        
        // Run the full EM algorithm
        runKallistoStyleEM();
    }
    
    // Update coverage statistics based on EM results
    updateCoverageStatsFromKallistoEM();
    
    // Generate reports
    reportEMResults();
    analyzeBetaLactamaseAssignments();
    
    std::cout << "Kallisto-style EM resolution complete" << std::endl;
}

std::string AMRDetectionPipeline::extractGeneFamily(const std::string& gene_name) {
    std::string name = gene_name;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    // Beta-lactamase families (highly important for clinical decisions)
    if (name.find("ctx-m") == 0) return "CTX-M";
    if (name.find("tem-") == 0 || name.find("tem") == 0) return "TEM";
    if (name.find("shv-") == 0 || name.find("shv") == 0) return "SHV";
    if (name.find("oxa-") == 0 || name.find("oxa") == 0) return "OXA";
    if (name.find("kpc-") == 0 || name.find("kpc") == 0) return "KPC";
    if (name.find("ndm-") == 0 || name.find("ndm") == 0) return "NDM";
    if (name.find("vim-") == 0 || name.find("vim") == 0) return "VIM";
    if (name.find("imp-") == 0 || name.find("imp") == 0) return "IMP";
    
    // Quinolone resistance families
    if (name.find("qnr") == 0) {
        // Extract qnrA, qnrB, qnrC, etc.
        size_t alpha_pos = name.find_first_of("abcdefghijklmnopqrstuvwxyz", 3);
        if (alpha_pos != std::string::npos && alpha_pos < 5) {
            return "qnr" + std::string(1, std::toupper(name[alpha_pos]));
        }
        return "qnr";
    }
    if (name.find("qep") == 0) return "qepA";
    if (name.find("aac(6')-ib-cr") != std::string::npos) return "aac6-Ib-cr";
    
    // Aminoglycoside resistance families
    if (name.find("aac(6')") == 0) return "aac6";
    if (name.find("aac(3)") == 0) return "aac3";
    if (name.find("ant(2\")") == 0) return "ant2";
    if (name.find("ant(3\")") == 0) return "ant3";
    if (name.find("aph(3')") == 0) return "aph3";
    if (name.find("aph(6)") == 0) return "aph6";
    if (name.find("arma") == 0) return "armA";
    if (name.find("rmta") == 0) return "rmtA";
    if (name.find("rmtb") == 0) return "rmtB";
    
    // Tetracycline resistance
    if (name.find("tet(") == 0) {
        size_t paren_pos = name.find('(');
        size_t close_paren = name.find(')', paren_pos);
        if (paren_pos != std::string::npos && close_paren != std::string::npos) {
            return "tet" + name.substr(paren_pos, close_paren - paren_pos + 1);
        }
        return "tet";
    }
    
    // Macrolide resistance
    if (name.find("erm") == 0) return "erm";
    if (name.find("mef") == 0) return "mef";
    if (name.find("mph") == 0) return "mph";
    
    // Vancomycin resistance
    if (name.find("van") == 0) {
        if (name.length() > 3) {
            return "van" + std::string(1, std::toupper(name[3]));
        }
        return "van";
    }
    
    // Colistin resistance (very important clinically)
    if (name.find("mcr-") == 0 || name.find("mcr") == 0) return "mcr";
    
    // Sulfonamide resistance
    if (name.find("sul") == 0) return "sul";
    
    // Generic pattern extraction for unknown families
    // Extract before dash or number
    size_t dash_pos = name.find('-');
    if (dash_pos != std::string::npos) {
        return name.substr(0, dash_pos);
    }
    
    // Remove trailing numbers
    size_t num_pos = name.find_first_of("0123456789");
    if (num_pos != std::string::npos && num_pos > 0) {
        return name.substr(0, num_pos);
    }
    
    return name;  // Return as-is if no pattern matched
}

void AMRDetectionPipeline::buildGeneFamiliesMap() {
    gene_families_map.clear();
    
    std::cout << "Building gene families map..." << std::endl;
    
    for (size_t i = 0; i < gene_entries.size(); i++) {
        std::string family = extractGeneFamily(gene_entries[i].gene_name);
        gene_families_map[family].push_back(i);
    }
    
    std::cout << "Built " << gene_families_map.size() << " gene families:" << std::endl;
    
    // Report family sizes (helpful for debugging)
    for (const auto& [family, gene_ids] : gene_families_map) {
        if (gene_ids.size() > 1) {  // Only report multi-gene families
            std::cout << "  " << family << ": " << gene_ids.size() << " variants" << std::endl;
            
            // Show first few variants for important families
            if (gene_ids.size() <= 5) {
                std::cout << "    Variants: ";
                for (size_t i = 0; i < gene_ids.size(); i++) {
                    if (i > 0) std::cout << ", ";
                    std::cout << gene_entries[gene_ids[i]].gene_name;
                }
                std::cout << std::endl;
            }
        }
    }
}

float AMRDetectionPipeline::calculateAssignmentScore(const AMRHit& hit) {
    // Multi-factor scoring for assignment quality
    float identity_score = hit.identity;
    float coverage_score = hit.coverage;
    
    // Length penalty: favor longer alignments
    float length_penalty = std::min(1.0f, (float)(hit.ref_end - hit.ref_start) / 50.0f);
    
    // Concordance bonus for paired-end reads
    float concordance_bonus = hit.concordant ? 1.2f : 1.0f;
    
    // Frame penalty: prefer certain frames (optional)
    float frame_penalty = (abs(hit.frame) <= 3) ? 1.0f : 0.9f;
    
    return identity_score * coverage_score * length_penalty * concordance_bonus * frame_penalty;
}

bool AMRDetectionPipeline::isHighQualityHit(const AMRHit& hit) {
    // Apply quality filters
    
    // Smith-Waterman score filter (primary filter)
    // Note: We need to add alignment_score to AMRHit structure or use a proxy
    // Use identity * alignment_length as proxy (ignore gene coverage for quality)
    float alignment_length = hit.ref_end - hit.ref_start;
    float proxy_score = hit.identity * alignment_length;
    
    // Debug output for first few hits
    // Note: This will show debug for every call, consider adding a member variable to track this
    static thread_local int debug_count = 0;
    if (debug_count < 5) {
        std::cout << "  Hit debug: Gene=" << hit.gene_name 
                  << " Identity=" << hit.identity 
                  << " Coverage=" << hit.coverage 
                  << " Length=" << alignment_length
                  << " ProxyScore=" << proxy_score
                  << " (MIN_SCORE=" << MIN_SMITH_WATERMAN_SCORE << ")" << std::endl;
        debug_count++;
    }
    
    if (proxy_score < MIN_SMITH_WATERMAN_SCORE) {
        return false;
    }
    
    // Minimum identity threshold
    if (hit.identity < 0.85f) {
        return false;
    }
    
    // Minimum coverage threshold (only apply if configured)
    if (config.min_hit_coverage >= 0.0f && hit.coverage < config.min_hit_coverage) {
        return false;
    }
    
    // Minimum alignment length
    if ((hit.ref_end - hit.ref_start) < 20) {
        return false;
    }
    
    return true;
}

void AMRDetectionPipeline::buildReadAssignments(const std::vector<AMRHit>& all_hits) {
    read_assignments.clear();
    std::map<uint32_t, ReadAssignment> read_map;
    
    std::cout << "Building read assignments from " << all_hits.size() << " hits..." << std::endl;
    
    // First, build a map of pair_id to hits for concordance checking
    std::map<uint32_t, std::vector<const AMRHit*>> pair_hits;
    for (const auto& hit : all_hits) {
        pair_hits[hit.pair_id].push_back(&hit);
    }
    
    int high_quality_hits = 0;
    int filtered_hits = 0;
    int concordant_bonuses = 0;
    
    // Group hits by read_id and filter for quality
    for (const auto& hit : all_hits) {
        if (isHighQualityHit(hit)) {
            high_quality_hits++;
            
            auto& assignment = read_map[hit.read_id];
            assignment.read_id = hit.read_id;
            assignment.candidate_genes.push_back(hit.gene_id);
            
            float score = calculateAssignmentScore(hit);
            
            // Apply paired-end concordance bonus if applicable
            if (hit.concordant && hit.pair_score > 0) {
                score *= config.concordance_bonus;
                concordant_bonuses++;
            }
            
            // Check for fragment length probability if we have it
            if (fragment_dist_estimated && hit.concordant) {
                // Find the corresponding pair alignment
                for (const auto& pair_align : read_pair_alignments) {
                    if ((hit.is_read2 ? pair_align.r2_read_idx : pair_align.r1_read_idx) == hit.read_id) {
                        // Find this gene in the concordant list
                        for (size_t i = 0; i < pair_align.concordant_genes.size(); i++) {
                            if (pair_align.concordant_genes[i] == hit.gene_id) {
                                float frag_prob = pair_align.fragment_probabilities[i];
                                score *= frag_prob;
                                break;
                            }
                        }
                        break;
                    }
                }
            }
            
            assignment.alignment_scores.push_back(score);
            assignment.high_quality = true;
        } else {
            filtered_hits++;
        }
    }
    
    std::cout << "Quality filtering results:" << std::endl;
    std::cout << "  High quality hits: " << high_quality_hits << std::endl;
    std::cout << "  Filtered out: " << filtered_hits << std::endl;
    std::cout << "  Concordant bonuses applied: " << concordant_bonuses << std::endl;
    if (fragment_dist_estimated) {
        std::cout << "  Fragment length distribution: mean=" << estimated_fragment_mean 
                  << " bp, std=" << estimated_fragment_std << " bp" << std::endl;
    }
    
    // Convert to vector and calculate total scores
    for (auto& [read_id, assignment] : read_map) {
        if (assignment.candidate_genes.size() >= 1) {
            // Calculate total score for normalization
            assignment.total_score = 0.0f;
            for (float score : assignment.alignment_scores) {
                assignment.total_score += score;
            }
            
            // Initialize probabilities uniformly for now
            assignment.assignment_probabilities.resize(assignment.candidate_genes.size(), 
                                                     1.0f / assignment.candidate_genes.size());
            
            read_assignments.push_back(assignment);
        }
    }
    
    // Report assignment statistics
    int multi_assignment_reads = 0;
    int total_assignments = 0;
    
    for (const auto& assignment : read_assignments) {
        total_assignments += assignment.candidate_genes.size();
        if (assignment.candidate_genes.size() > 1) {
            multi_assignment_reads++;
        }
    }
    
    std::cout << "Assignment statistics:" << std::endl;
    std::cout << "  Reads with assignments: " << read_assignments.size() << std::endl;
    std::cout << "  Reads with multiple candidates: " << multi_assignment_reads << std::endl;
    std::cout << "  Average candidates per read: " << 
                 (read_assignments.empty() ? 0.0f : (float)total_assignments / read_assignments.size()) << std::endl;
}

void AMRDetectionPipeline::buildPairedReadAssignments(const std::vector<AMRHit>& all_hits) {
    paired_read_assignments.clear();
    std::map<uint32_t, PairedReadAssignment> pair_map;
    
    std::cout << "Building paired read assignments from " << all_hits.size() << " hits..." << std::endl;
    
    // Group hits by pair_id
    std::map<uint32_t, std::vector<const AMRHit*>> hits_by_pair;
    for (const auto& hit : all_hits) {
        if (isHighQualityHit(hit)) {
            hits_by_pair[hit.pair_id].push_back(&hit);
        }
    }
    
    std::cout << "Found " << hits_by_pair.size() << " read pairs with hits" << std::endl;
    
    // Build assignments for each pair
    for (const auto& [pair_id, pair_hits] : hits_by_pair) {
        PairedReadAssignment pair_assign;
        pair_assign.pair_id = pair_id;
        
        // Separate R1 and R2 hits
        std::vector<const AMRHit*> r1_hits, r2_hits;
        for (const auto* hit : pair_hits) {
            if (hit->is_read2) {
                r2_hits.push_back(hit);
                pair_assign.r2_read_id = hit->read_id;
            } else {
                r1_hits.push_back(hit);
                pair_assign.r1_read_id = hit->read_id;
            }
        }
        
        // Process R1 hits
        std::set<uint32_t> r1_genes;
        for (const auto* hit : r1_hits) {
            float score = calculateAssignmentScore(*hit);
            pair_assign.r1_gene_scores[hit->gene_id] = score;
            r1_genes.insert(hit->gene_id);
        }
        
        // Process R2 hits
        std::set<uint32_t> r2_genes;
        for (const auto* hit : r2_hits) {
            float score = calculateAssignmentScore(*hit);
            pair_assign.r2_gene_scores[hit->gene_id] = score;
            r2_genes.insert(hit->gene_id);
        }
        
        // Find compatible genes (union of R1 and R2 genes)
        std::set<uint32_t> all_genes;
        all_genes.insert(r1_genes.begin(), r1_genes.end());
        all_genes.insert(r2_genes.begin(), r2_genes.end());
        pair_assign.compatible_genes.assign(all_genes.begin(), all_genes.end());
        
        // Calculate concordance scores and fragment information
        for (uint32_t gene_id : all_genes) {
            if (pair_assign.isConcordant(gene_id)) {
                // Both reads hit this gene - calculate combined score
                float r1_score = pair_assign.r1_gene_scores[gene_id];
                float r2_score = pair_assign.r2_gene_scores[gene_id];
                
                // Find the actual hits to get position information
                const AMRHit* r1_hit = nullptr;
                const AMRHit* r2_hit = nullptr;
                for (const auto* hit : r1_hits) {
                    if (hit->gene_id == gene_id) {
                        r1_hit = hit;
                        break;
                    }
                }
                for (const auto* hit : r2_hits) {
                    if (hit->gene_id == gene_id) {
                        r2_hit = hit;
                        break;
                    }
                }
                
                if (r1_hit && r2_hit) {
                    // Calculate fragment length
                    AMRGeneEntry gene_entry = getGeneEntry(gene_id);
                    int32_t gene_length = gene_entry.protein_length;
                    
                    // For protein alignments, positions are in amino acids
                    int32_t r1_nt_start = r1_hit->ref_start * 3;
                    int32_t r2_nt_end = r2_hit->ref_end * 3;
                    
                    if (r2_nt_end > r1_nt_start) {
                        int32_t fragment_length = r2_nt_end - r1_nt_start + 150; // Add read length
                        pair_assign.fragment_lengths[gene_id] = fragment_length;
                        
                        // Calculate fragment probability
                        float frag_prob = calculateFragmentLengthProbability(fragment_length);
                        pair_assign.fragment_probabilities[gene_id] = frag_prob;
                        
                        // Combined score with concordance bonus and fragment probability
                        float combined_score = (r1_score + r2_score) * config.concordance_bonus * frag_prob;
                        pair_assign.concordance_scores[gene_id] = combined_score;
                    }
                }
            }
        }
        
        // Initialize uniform probabilities
        if (!pair_assign.compatible_genes.empty()) {
            float uniform_prob = 1.0f / pair_assign.compatible_genes.size();
            for (uint32_t gene_id : pair_assign.compatible_genes) {
                pair_assign.assignment_probabilities[gene_id] = uniform_prob;
            }
            paired_read_assignments.push_back(pair_assign);
        }
    }
    
    // Report statistics
    int concordant_pairs = 0;
    int discordant_pairs = 0;
    int single_end_pairs = 0;
    
    for (const auto& pair : paired_read_assignments) {
        bool has_concordant = false;
        for (uint32_t gene_id : pair.compatible_genes) {
            if (pair.isConcordant(gene_id)) {
                has_concordant = true;
                break;
            }
        }
        
        if (has_concordant) {
            concordant_pairs++;
        } else if (pair.r1_gene_scores.empty() || pair.r2_gene_scores.empty()) {
            single_end_pairs++;
        } else {
            discordant_pairs++;
        }
    }
    
    std::cout << "Paired assignment statistics:" << std::endl;
    std::cout << "  Total pairs with assignments: " << paired_read_assignments.size() << std::endl;
    std::cout << "  Concordant pairs: " << concordant_pairs << std::endl;
    std::cout << "  Discordant pairs: " << discordant_pairs << std::endl;
    std::cout << "  Single-end pairs: " << single_end_pairs << std::endl;
}

float AMRDetectionPipeline::calculateNameSimilarity(const std::string& name1, const std::string& name2) {
    if (name1 == name2) return 1.0f;
    
    // Parse gene names into family and variant
    auto [family1, variant1] = parseGeneName(name1);
    auto [family2, variant2] = parseGeneName(name2);
    
    if (family1 != family2) {
        // Different families have low similarity
        return 0.05f;
    }
    
    // Same family, different variants
    if (variant1.empty() || variant2.empty()) {
        return 0.7f;  // One has no variant info
    }
    
    if (variant1 == variant2) {
        return 0.95f;  // Same variant (should be identical)
    }
    
    // Different variants of same family
    // Special cases for known relationships
    if (family1 == "CTX-M") {
        // CTX-M variants have complex relationships
        return calculateCTXMSimilarity(variant1, variant2);
    } else if (family1 == "TEM") {
        return calculateTEMSimilarity(variant1, variant2);
    } else if (family1 == "qnrA" || family1 == "qnrB" || family1 == "qnrC") {
        return 0.8f;  // qnr variants are quite similar
    }
    
    // Default similarity for same family, different variants
    return 0.7f;
}

std::pair<std::string, std::string> AMRDetectionPipeline::parseGeneName(const std::string& name) {
    std::string family = extractGeneFamily(name);
    
    // Extract variant part
    size_t family_len = family.length();
    if (name.length() > family_len) {
        std::string remainder = name.substr(family_len);
        
        // Remove common separators
        if (!remainder.empty() && (remainder[0] == '-' || remainder[0] == '_')) {
            remainder = remainder.substr(1);
        }
        
        return {family, remainder};
    }
    
    return {family, ""};
}

float AMRDetectionPipeline::calculateCTXMSimilarity(const std::string& var1, const std::string& var2) {
    // CTX-M variants grouped by phylogenetic relationships
    static const std::map<std::string, int> ctxm_groups = {
        {"1", 1}, {"3", 1}, {"10", 1}, {"11", 1}, {"12", 1}, {"15", 1}, {"23", 1}, {"28", 1},
        {"2", 2}, {"4", 2}, {"5", 2}, {"6", 2}, {"7", 2}, {"20", 2},
        {"8", 3}, {"40", 3}, {"41", 3},
        {"9", 4}, {"13", 4}, {"14", 4}, {"16", 4}, {"17", 4}, {"18", 4}, {"19", 4}, {"21", 4},
        {"25", 5}, {"26", 5}, {"27", 5}, {"39", 5}, {"42", 5}
    };
    
    auto it1 = ctxm_groups.find(var1);
    auto it2 = ctxm_groups.find(var2);
    
    if (it1 != ctxm_groups.end() && it2 != ctxm_groups.end()) {
        if (it1->second == it2->second) {
            return 0.9f;  // Same phylogenetic group
        } else {
            return 0.6f;  // Different groups
        }
    }
    
    return 0.7f;  // Unknown variants
}

float AMRDetectionPipeline::calculateTEMSimilarity(const std::string& var1, const std::string& var2) {
    // TEM variants: TEM-1 is ancestral, others are derivatives
    if (var1 == "1" || var2 == "1") {
        return 0.8f;  // High similarity to TEM-1
    }
    
    // Other TEM variants
    return 0.75f;
}

void AMRDetectionPipeline::buildSequenceSimilarityMatrix() {
    size_t n = gene_entries.size();
    gene_similarity_matrix.assign(n, std::vector<float>(n, 0.0f));
    
    std::cout << "Building sequence similarity matrix for " << n << " genes..." << std::endl;
    
    // Calculate pairwise similarities
    for (size_t i = 0; i < n; i++) {
        gene_similarity_matrix[i][i] = 1.0f;  // Self-similarity
        
        for (size_t j = i + 1; j < n; j++) {
            float similarity = calculateNameSimilarity(gene_entries[i].gene_name, 
                                                     gene_entries[j].gene_name);
            gene_similarity_matrix[i][j] = similarity;
            gene_similarity_matrix[j][i] = similarity;
        }
        
        // Progress reporting for large matrices
        if (i % 100 == 0 && i > 0) {
            std::cout << "  Processed " << i << "/" << n << " genes" << std::endl;
        }
    }
    
    // Report high-similarity pairs for verification
    std::cout << "High similarity gene pairs (>0.8):" << std::endl;
    int high_sim_count = 0;
    for (size_t i = 0; i < n && high_sim_count < 20; i++) {
        for (size_t j = i + 1; j < n && high_sim_count < 20; j++) {
            if (gene_similarity_matrix[i][j] > 0.8f) {
                std::cout << "  " << gene_entries[i].gene_name << " <-> " 
                          << gene_entries[j].gene_name << ": " 
                          << std::fixed << std::setprecision(2) 
                          << gene_similarity_matrix[i][j] << std::endl;
                high_sim_count++;
            }
        }
    }
}

float AMRDetectionPipeline::getGeneFamilyPrior(const std::string& family, const std::string& gene_name) {
    // Known common variants get higher priors based on clinical prevalence
    static const std::map<std::string, float> common_variants = {
        // CTX-M family (very common ESBLs)
        {"CTX-M-15", 3.0f}, {"CTX-M-1", 2.5f}, {"CTX-M-14", 2.0f}, 
        {"CTX-M-27", 1.8f}, {"CTX-M-3", 1.5f}, {"CTX-M-2", 1.3f},
        
        // TEM family
        {"TEM-1", 2.5f}, {"TEM-52", 1.8f}, {"TEM-116", 1.5f}, {"TEM-135", 1.3f},
        
        // SHV family
        {"SHV-1", 2.0f}, {"SHV-12", 1.5f}, {"SHV-5", 1.3f}, {"SHV-2", 1.2f},
        
        // Carbapenemases (high priority)
        {"KPC-2", 2.0f}, {"KPC-3", 1.8f}, {"NDM-1", 2.2f}, {"NDM-5", 1.5f},
        {"OXA-48", 1.8f}, {"OXA-181", 1.3f}, {"VIM-1", 1.5f}, {"IMP-1", 1.3f},
        
        // Quinolone resistance
        {"qnrA1", 1.5f}, {"qnrB1", 1.8f}, {"qnrB2", 1.3f}, {"qnrB4", 1.2f},
        {"qnrB6", 1.1f}, {"qnrC", 1.0f}, {"qnrS1", 1.4f},
        
        // Aminoglycoside resistance
        {"aac(6')-Ib", 1.8f}, {"aac(6')-Ib-cr", 1.5f}, {"aac(3)-Ia", 1.3f},
        {"ant(2\")-Ia", 1.2f}, {"aph(3')-Ia", 1.4f}, {"armA", 1.6f}, {"rmtB", 1.3f},
        
        // Colistin resistance (rare but critical)
        {"mcr-1", 0.8f}, {"mcr-2", 0.3f}, {"mcr-3", 0.2f},
        
        // Vancomycin resistance (rare but critical)
        {"vanA", 0.5f}, {"vanB", 0.3f}, {"vanC", 0.2f}
    };
    
    auto it = common_variants.find(gene_name);
    if (it != common_variants.end()) {
        return it->second;
    }
    
    // Family-level priors for variants not specifically listed
    static const std::map<std::string, float> family_priors = {
        {"CTX-M", 1.5f}, {"TEM", 1.3f}, {"SHV", 1.2f},
        {"KPC", 1.0f}, {"NDM", 1.0f}, {"OXA", 1.1f}, {"VIM", 0.8f}, {"IMP", 0.7f},
        {"qnrA", 1.0f}, {"qnrB", 1.2f}, {"qnrC", 0.8f}, {"qnrS", 1.0f},
        {"aac6", 1.2f}, {"aac3", 1.0f}, {"ant2", 0.9f}, {"ant3", 0.8f},
        {"aph3", 1.0f}, {"armA", 1.0f}, {"rmtA", 0.8f}, {"rmtB", 0.9f},
        {"tet", 1.1f}, {"sul", 1.0f}, {"erm", 0.9f}, {"mef", 0.8f},
        {"mcr", 0.3f}, {"vanA", 0.3f}, {"vanB", 0.2f}
    };
    
    auto family_it = family_priors.find(family);
    if (family_it != family_priors.end()) {
        return family_it->second;
    }
    
    return 1.0f;  // Default prior
}

void AMRDetectionPipeline::initializeGeneAbundances() {
    gene_abundance_info.clear();
    
    std::cout << "Initializing gene abundances with clinical priors..." << std::endl;
    
    for (size_t i = 0; i < gene_entries.size(); i++) {
        GeneAbundanceInfo info;
        info.gene_id = i;
        info.gene_name = gene_entries[i].gene_name;
        info.gene_family = extractGeneFamily(gene_entries[i].gene_name);
        
        // Effective length (protein length converted to nucleotides)
        info.effective_length = gene_entries[i].protein_length * 3.0f;
        
        // Clinical prevalence prior
        info.prior_weight = getGeneFamilyPrior(info.gene_family, info.gene_name);
        
        // Initialize abundance with prior (will be updated by EM)
        info.abundance = info.prior_weight;
        info.total_reads = 0.0f;
        
        gene_abundance_info[i] = info;
    }
    
    // Report high-priority genes
    std::cout << "High-priority genes (prior > 1.5):" << std::endl;
    int count = 0;
    for (const auto& [gene_id, info] : gene_abundance_info) {
        if (info.prior_weight > 1.5f && count < 15) {
            std::cout << "  " << info.gene_name << " (" << info.gene_family 
                      << "): prior = " << std::fixed << std::setprecision(1) 
                      << info.prior_weight << std::endl;
            count++;
        }
    }
    
    std::cout << "Initialized " << gene_abundance_info.size() << " genes" << std::endl;
}

void AMRDetectionPipeline::normalizeInitialAbundances() {
    // Normalize abundances within families to prevent dominance
    std::map<std::string, float> family_totals;
    std::map<std::string, int> family_counts;
    
    // Calculate family totals
    for (const auto& [gene_id, info] : gene_abundance_info) {
        family_totals[info.gene_family] += info.abundance;
        family_counts[info.gene_family]++;
    }
    
    // Apply gentle normalization within families with multiple members
    for (auto& [gene_id, info] : gene_abundance_info) {
        if (family_counts[info.gene_family] > 1) {
            float family_avg = family_totals[info.gene_family] / family_counts[info.gene_family];
            
            // Gentle normalization: move 30% toward family average
            info.abundance = 0.7f * info.abundance + 0.3f * family_avg;
        }
    }
    
    std::cout << "Applied family-level abundance normalization" << std::endl;
}

void AMRDetectionPipeline::updateAssignmentProbabilities() {
    // E-step: Update assignment probabilities based on current abundances
    
    for (auto& assignment : read_assignments) {
        std::vector<float> weights;
        weights.reserve(assignment.candidate_genes.size());
        
        for (size_t i = 0; i < assignment.candidate_genes.size(); i++) {
            uint32_t gene_id = assignment.candidate_genes[i];
            float alignment_score = assignment.alignment_scores[i];
            
            const auto& gene_info = gene_abundance_info.at(gene_id);
            
            // Kallisto-style weight: abundance / effective_length
            float abundance_weight = gene_info.abundance / (gene_info.effective_length + 1e-10f);
            
            // Sequence similarity boost
            float similarity_boost = 1.0f;
            if (SIMILARITY_WEIGHT > 0 && assignment.candidate_genes.size() > 1) {
                float max_similarity_contribution = 0.0f;
                
                for (size_t j = 0; j < assignment.candidate_genes.size(); j++) {
                    if (i != j) {
                        uint32_t other_gene = assignment.candidate_genes[j];
                        float similarity = gene_similarity_matrix[gene_id][other_gene];
                        float other_abundance = gene_abundance_info.at(other_gene).abundance;
                        
                        max_similarity_contribution = std::max(max_similarity_contribution, 
                                                             similarity * other_abundance);
                    }
                }
                
                similarity_boost = 1.0f + SIMILARITY_WEIGHT * max_similarity_contribution;
            }
            
            // Combined weight
            float weight = alignment_score * abundance_weight * similarity_boost;
            weights.push_back(weight);
        }
        
        // Normalize to probabilities
        float total_weight = 0.0f;
        for (float w : weights) {
            total_weight += w;
        }
        
        if (total_weight > 0) {
            for (size_t i = 0; i < weights.size(); i++) {
                assignment.assignment_probabilities[i] = weights[i] / total_weight;
            }
        } else {
            // Fallback to uniform distribution
            float uniform_prob = 1.0f / assignment.candidate_genes.size();
            std::fill(assignment.assignment_probabilities.begin(), 
                     assignment.assignment_probabilities.end(), uniform_prob);
        }
    }
}

void AMRDetectionPipeline::updateGeneAbundances() {
    // M-step: Update gene abundances based on assignment probabilities
    
    // Reset read counts
    for (auto& [gene_id, info] : gene_abundance_info) {
        info.total_reads = 0.0f;
    }
    
    // Accumulate fractional read assignments
    for (const auto& assignment : read_assignments) {
        for (size_t i = 0; i < assignment.candidate_genes.size(); i++) {
            uint32_t gene_id = assignment.candidate_genes[i];
            float probability = assignment.assignment_probabilities[i];
            
            gene_abundance_info[gene_id].total_reads += probability;
        }
    }
    
    // Update abundances (RPKM-like calculation)
    for (auto& [gene_id, info] : gene_abundance_info) {
        // Include prior in abundance calculation
        float prior_contribution = 0.1f * info.prior_weight;  // 10% weight to prior
        
        // RPKM = (reads * 1000) / (gene_length_kb)
        float rpkm = (info.total_reads + prior_contribution) / (info.effective_length / 1000.0f);
        
        info.abundance = rpkm;
    }
}

void AMRDetectionPipeline::runPairedEndEM() {
    std::cout << "\n=== Running Paired-End Aware EM Algorithm ===" << std::endl;
    
    if (paired_read_assignments.empty()) {
        std::cout << "No paired assignments to process" << std::endl;
        return;
    }
    
    // Initialize gene families, similarity matrix, and abundances
    buildGeneFamiliesMap();
    buildSequenceSimilarityMatrix();
    initializeGeneAbundances();
    normalizeInitialAbundances();
    
    // Store previous abundances for convergence checking
    std::map<uint32_t, float> prev_abundances;
    
    // EM iterations
    for (int iter = 0; iter < config.em_iterations; iter++) {
        std::cout << "\nEM Iteration " << (iter + 1) << "/" << config.em_iterations << std::endl;
        
        // Store current abundances for convergence check
        for (const auto& [gene_id, info] : gene_abundance_info) {
            prev_abundances[gene_id] = info.abundance;
        }
        
        // E-step: Update assignment probabilities
        updatePairedAssignmentProbabilities();
        
        // M-step: Update gene abundances
        updateGeneAbundancesFromPairs();
        
        // Apply family constraints if needed
        applyFamilyConstraints();
        
        // Check convergence
        float max_change = 0.0f;
        for (const auto& [gene_id, info] : gene_abundance_info) {
            float change = std::abs(info.abundance - prev_abundances[gene_id]);
            max_change = std::max(max_change, change);
        }
        
        std::cout << "Max abundance change: " << max_change << std::endl;
        
        if (max_change < config.em_convergence) {
            std::cout << "EM converged after " << (iter + 1) << " iterations" << std::endl;
            break;
        }
    }
}

void AMRDetectionPipeline::updatePairedAssignmentProbabilities() {
    // E-step for paired reads
    
    for (auto& pair_assign : paired_read_assignments) {
        std::map<uint32_t, float> gene_weights;
        float total_weight = 0.0f;
        
        for (uint32_t gene_id : pair_assign.compatible_genes) {
            float weight = calculatePairProbability(pair_assign, gene_id);
            gene_weights[gene_id] = weight;
            total_weight += weight;
        }
        
        // Normalize to probabilities
        if (total_weight > 0) {
            for (auto& [gene_id, weight] : gene_weights) {
                pair_assign.assignment_probabilities[gene_id] = weight / total_weight;
            }
        } else {
            // Fallback to uniform
            float uniform_prob = 1.0f / pair_assign.compatible_genes.size();
            for (uint32_t gene_id : pair_assign.compatible_genes) {
                pair_assign.assignment_probabilities[gene_id] = uniform_prob;
            }
        }
        
        pair_assign.total_probability = total_weight;
    }
}

float AMRDetectionPipeline::calculatePairProbability(
    const PairedReadAssignment& pair, 
    uint32_t gene_id
) {
    const auto& gene_info = gene_abundance_info.at(gene_id);
    
    // Base probability from abundance
    float abundance_weight = gene_info.abundance / (gene_info.effective_length + 1e-10f);
    
    // Alignment scores
    float alignment_weight = 0.0f;
    
    if (pair.isConcordant(gene_id)) {
        // Concordant pair - use combined score with fragment length probability
        if (pair.concordance_scores.count(gene_id)) {
            alignment_weight = pair.concordance_scores.at(gene_id);
        } else {
            // Fallback to sum of individual scores
            alignment_weight = pair.getCombinedScore(gene_id) * config.concordance_bonus;
        }
        
        // Apply fragment length probability if available
        if (pair.fragment_probabilities.count(gene_id)) {
            alignment_weight *= pair.fragment_probabilities.at(gene_id);
        }
    } else {
        // Discordant or single-end - use individual scores with penalty
        alignment_weight = pair.getCombinedScore(gene_id) * config.discord_penalty;
    }
    
    // Sequence similarity boost (if multiple candidates)
    float similarity_boost = 1.0f;
    if (SIMILARITY_WEIGHT > 0 && pair.compatible_genes.size() > 1) {
        float max_similarity_contribution = 0.0f;
        
        for (uint32_t other_gene : pair.compatible_genes) {
            if (other_gene != gene_id) {
                float similarity = gene_similarity_matrix[gene_id][other_gene];
                float other_abundance = gene_abundance_info.at(other_gene).abundance;
                
                max_similarity_contribution = std::max(max_similarity_contribution, 
                                                     similarity * other_abundance);
            }
        }
        
        similarity_boost = 1.0f + SIMILARITY_WEIGHT * max_similarity_contribution;
    }
    
    return alignment_weight * abundance_weight * similarity_boost;
}

void AMRDetectionPipeline::updateGeneAbundancesFromPairs() {
    // M-step for paired reads
    
    // Reset read counts
    for (auto& [gene_id, info] : gene_abundance_info) {
        info.total_reads = 0.0f;
    }
    
    // Accumulate fractional pair assignments
    for (const auto& pair_assign : paired_read_assignments) {
        for (const auto& [gene_id, prob] : pair_assign.assignment_probabilities) {
            // Each pair counts as 2 reads
            gene_abundance_info[gene_id].total_reads += 2.0f * prob;
        }
    }
    
    // Also include unpaired reads if we have them
    for (const auto& assignment : read_assignments) {
        // Check if this read is part of a pair that was already counted
        bool is_paired = false;
        for (const auto& pair : paired_read_assignments) {
            if (assignment.read_id == pair.r1_read_id || 
                assignment.read_id == pair.r2_read_id) {
                is_paired = true;
                break;
            }
        }
        
        if (!is_paired) {
            // This is an orphaned read, count it separately
            for (size_t i = 0; i < assignment.candidate_genes.size(); i++) {
                uint32_t gene_id = assignment.candidate_genes[i];
                float probability = assignment.assignment_probabilities[i];
                gene_abundance_info[gene_id].total_reads += probability;
            }
        }
    }
    
    // Update abundances (RPKM-like calculation)
    float total_assigned_reads = 0.0f;
    for (const auto& [gene_id, info] : gene_abundance_info) {
        total_assigned_reads += info.total_reads;
    }
    
    for (auto& [gene_id, info] : gene_abundance_info) {
        // Include prior in abundance calculation
        float prior_contribution = 0.1f * info.prior_weight;
        
        // RPKM = (reads * 1e6) / (gene_length_kb * total_reads)
        float rpkm = ((info.total_reads + prior_contribution) * 1e6f) / 
                     ((info.effective_length / 1000.0f) * (total_assigned_reads + 1.0f));
        
        info.abundance = rpkm;
    }
}

void AMRDetectionPipeline::applyFamilyConstraints() {
    // Optional: Apply constraints within gene families
    const float FAMILY_SMOOTHING = 0.05f;  // 5% smoothing within families
    
    for (const auto& [family_name, gene_ids] : gene_families_map) {
        if (gene_ids.size() <= 1) continue;
        
        // Calculate family statistics
        float family_total_abundance = 0.0f;
        float family_total_reads = 0.0f;
        
        for (uint32_t gene_id : gene_ids) {
            family_total_abundance += gene_abundance_info[gene_id].abundance;
            family_total_reads += gene_abundance_info[gene_id].total_reads;
        }
        
        if (family_total_reads > 0) {
            float family_mean_abundance = family_total_abundance / gene_ids.size();
            
            // Apply gentle smoothing toward family mean
            for (uint32_t gene_id : gene_ids) {
                auto& info = gene_abundance_info[gene_id];
                info.abundance = (1.0f - FAMILY_SMOOTHING) * info.abundance + 
                                FAMILY_SMOOTHING * family_mean_abundance;
            }
        }
    }
}

void AMRDetectionPipeline::runKallistoStyleEM() {
    std::cout << "\n=== Running Kallisto-style EM Algorithm ===" << std::endl;
    std::cout << "Processing " << read_assignments.size() << " read assignments" << std::endl;
    
    if (read_assignments.empty()) {
        std::cout << "No read assignments available for EM" << std::endl;
        return;
    }
    
    // Initialize
    buildGeneFamiliesMap();
    buildSequenceSimilarityMatrix();
    initializeGeneAbundances();
    normalizeInitialAbundances();
    
    std::cout << "Starting EM iterations..." << std::endl;
    
    float prev_likelihood = -std::numeric_limits<float>::infinity();
    
    for (int iter = 0; iter < MAX_EM_ITERATIONS; iter++) {
        // E-step
        updateAssignmentProbabilities();
        
        // M-step
        updateGeneAbundances();
        
        // Optional family constraints
        if (iter % 5 == 0) {  // Apply every 5 iterations
            applyFamilyConstraints();
        }
        
        // Check convergence
        float max_change = 0.0f;
        for (const auto& [gene_id, info] : gene_abundance_info) {
            // Track changes in genes with significant reads
            if (info.total_reads > 0.1f) {
                max_change = std::max(max_change, std::abs(info.abundance));
            }
        }
        
        if (iter % 10 == 0 || iter < 5) {
            std::cout << "EM iteration " << iter << ", max abundance change: " 
                      << std::fixed << std::setprecision(4) << max_change << std::endl;
        }
        
        // Convergence check
        if (iter > 5 && max_change < EM_CONVERGENCE_THRESHOLD) {
            std::cout << "EM converged after " << iter + 1 << " iterations" << std::endl;
            break;
        }
    }
    
    std::cout << "EM algorithm completed" << std::endl;
    
    // Update coverage statistics and report results
    updateCoverageStatsFromKallistoEM();
    reportEMResults();
}

void AMRDetectionPipeline::updateCoverageStatsFromKallistoEM() {
    // Update coverage statistics based on Kallisto-style EM results
    std::cout << "Updating coverage statistics from EM results..." << std::endl;
    
    for (const auto& [gene_id, info] : gene_abundance_info) {
        if (gene_id < h_coverage_stats.size() && info.total_reads > 0.01f) {
            // Update read count with EM result
            h_coverage_stats[gene_id].total_reads = std::round(info.total_reads);
            
            // Calculate TPM and RPKM based on EM abundances
            h_coverage_stats[gene_id].rpkm = info.abundance;  // This is already RPKM-like
            
            // Proportional update of other metrics
            if (h_coverage_stats[gene_id].total_reads > 0) {
                float original_reads = h_coverage_stats[gene_id].total_reads;
                float scale_factor = info.total_reads / std::max(1.0f, original_reads);
                
                h_coverage_stats[gene_id].total_bases_mapped = 
                    std::round(h_coverage_stats[gene_id].total_bases_mapped * scale_factor);
            }
        }
    }
}

void AMRDetectionPipeline::reportEMResults() {
    std::cout << "\n=== EM Assignment Results by Gene Family ===" << std::endl;
    
    // Sort families by total abundance
    std::vector<std::pair<std::string, float>> family_abundances;
    
    for (const auto& [family_name, gene_ids] : gene_families_map) {
        float family_total = 0.0f;
        for (uint32_t gene_id : gene_ids) {
            family_total += gene_abundance_info[gene_id].total_reads;
        }
        
        if (family_total > 0.1f) {  // Only report families with significant reads
            family_abundances.push_back({family_name, family_total});
        }
    }
    
    std::sort(family_abundances.begin(), family_abundances.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Report top families
    for (const auto& [family_name, family_total] : family_abundances) {
        const auto& gene_ids = gene_families_map[family_name];
        
        std::vector<std::pair<float, std::string>> family_genes;
        for (uint32_t gene_id : gene_ids) {
            const auto& info = gene_abundance_info[gene_id];
            if (info.total_reads > 0.01f) {
                family_genes.push_back({info.total_reads, info.gene_name});
            }
        }
        
        if (!family_genes.empty()) {
            std::sort(family_genes.begin(), family_genes.end(),
                     [](const auto& a, const auto& b) { return a.first > b.first; });
            
            std::cout << "\n" << family_name << " family (total: " 
                      << std::fixed << std::setprecision(1) << family_total << " reads):" << std::endl;
            
            for (const auto& [reads, name] : family_genes) {
                float percentage = (reads / family_total) * 100.0f;
                std::cout << "  " << std::setw(15) << std::left << name << ": " 
                          << std::setw(6) << std::fixed << std::setprecision(1) << reads 
                          << " reads (" << std::setw(4) << percentage << "%)" << std::endl;
            }
        }
    }
}

void AMRDetectionPipeline::analyzeBetaLactamaseAssignments() {
    std::cout << "\n=== Beta-lactamase Family Analysis ===" << std::endl;
    
    // Use accumulated hits for analysis (not just current batch)
    auto hits = accumulated_hits;
    std::map<std::string, std::vector<AMRHit>> beta_lactamase_hits;
    
    for (const auto& hit : hits) {
        if (hit.gene_id < gene_entries.size()) {
            std::string gene_name = gene_entries[hit.gene_id].gene_name;
            std::string family = extractGeneFamily(gene_name);
            
            // Check if it's a beta-lactamase family
            if (family == "CTX-M" || family == "TEM" || family == "SHV" || 
                family == "OXA" || family == "KPC" || family == "NDM" || 
                family == "VIM" || family == "IMP") {
                beta_lactamase_hits[family].push_back(hit);
            }
        }
    }
    
    for (const auto& [family, family_hits] : beta_lactamase_hits) {
        std::map<std::string, int> variant_counts;
        std::set<uint32_t> unique_reads;
        std::map<std::string, float> variant_em_reads;
        
        // Count raw hits
        for (const auto& hit : family_hits) {
            std::string gene_name = gene_entries[hit.gene_id].gene_name;
            variant_counts[gene_name]++;
            unique_reads.insert(hit.read_id);
        }
        
        // Get EM-assigned reads for this family
        for (uint32_t gene_id : gene_families_map[family]) {
            const auto& info = gene_abundance_info[gene_id];
            if (info.total_reads > 0.01f) {
                variant_em_reads[info.gene_name] = info.total_reads;
            }
        }
        
        std::cout << "\n" << family << " family analysis:" << std::endl;
        std::cout << "  Raw hits: " << family_hits.size() << std::endl;
        std::cout << "  Unique reads: " << unique_reads.size() << std::endl;
        std::cout << "  Avg hits per read: " << std::fixed << std::setprecision(2) 
                  << (unique_reads.empty() ? 0.0f : (float)family_hits.size() / unique_reads.size()) << std::endl;
        
        std::cout << "  Raw variant distribution:" << std::endl;
        for (const auto& [variant, count] : variant_counts) {
            std::cout << "    " << std::setw(15) << std::left << variant 
                      << ": " << std::setw(4) << count << " hits" << std::endl;
        }
        
        if (!variant_em_reads.empty()) {
            std::cout << "  EM-resolved distribution:" << std::endl;
            for (const auto& [variant, em_reads] : variant_em_reads) {
                std::cout << "    " << std::setw(15) << std::left << variant 
                          << ": " << std::setw(6) << std::fixed << std::setprecision(1) 
                          << em_reads << " reads" << std::endl;
            }
        }
    }
}

void AMRDetectionPipeline::updateHitsWithPairedInfo(
    std::vector<AMRHit>& hits, 
    const std::vector<ReadPairData>& read_pairs
) {
    // Create a map for quick lookup of pair information
    std::map<uint32_t, uint32_t> read_to_pair;
    std::map<uint32_t, uint32_t> read_to_mate;
    
    for (const auto& pair : read_pairs) {
        // R1 and R2 indices
        uint32_t r1_idx = pair.pair_index * 2;
        uint32_t r2_idx = pair.pair_index * 2 + 1;
        
        read_to_pair[r1_idx] = pair.pair_index;
        read_to_pair[r2_idx] = pair.pair_index;
        read_to_mate[r1_idx] = r2_idx;
        read_to_mate[r2_idx] = r1_idx;
    }
    
    // Update each hit with paired information
    for (auto& hit : hits) {
        if (read_to_pair.count(hit.read_id)) {
            hit.pair_id = read_to_pair[hit.read_id];
            hit.is_read2 = (hit.read_id % 2 == 1);  // R2 has odd indices
            hit.mate_read_id = read_to_mate[hit.read_id];
        } else {
            // Single-end or orphaned read
            hit.pair_id = UINT32_MAX;
            hit.is_read2 = false;
            hit.mate_read_id = UINT32_MAX;
        }
        
        // Initialize pair_score to 0 (will be updated later)
        hit.pair_score = 0.0f;
    }
    
    // Now check for concordant hits
    std::map<std::pair<uint32_t, uint32_t>, std::vector<AMRHit*>> pair_gene_hits;
    
    for (auto& hit : hits) {
        if (hit.pair_id != UINT32_MAX) {
            pair_gene_hits[{hit.pair_id, hit.gene_id}].push_back(&hit);
        }
    }
    
    // Mark concordant hits and calculate pair scores
    for (auto& [pair_gene, hit_ptrs] : pair_gene_hits) {
        if (hit_ptrs.size() >= 2) {
            // Check if we have both R1 and R2
            bool has_r1 = false, has_r2 = false;
            float r1_score = 0.0f, r2_score = 0.0f;
            
            for (auto* hit : hit_ptrs) {
                float score = hit->identity * (hit->ref_end - hit->ref_start);
                if (hit->is_read2) {
                    has_r2 = true;
                    r2_score = std::max(r2_score, score);
                } else {
                    has_r1 = true;
                    r1_score = std::max(r1_score, score);
                }
            }
            
            if (has_r1 && has_r2) {
                // This is a concordant pair
                float combined_score = (r1_score + r2_score) * config.concordance_bonus;
                
                for (auto* hit : hit_ptrs) {
                    hit->concordant = true;
                    hit->pair_score = combined_score;
                }
            }
        }
    }
    
    std::cout << "Updated " << hits.size() << " hits with paired-end information" << std::endl;
}

void AMRDetectionPipeline::updateGPUHitsWithPairedInfo(const std::vector<AMRHit>& hits) {
    // Copy updated hits back to GPU
    if (!hits.empty() && d_amr_hits) {
        std::vector<AMRHit> gpu_hits(current_batch_size * MAX_MATCHES_PER_READ);
        
        // First copy existing hits from GPU
        cudaMemcpy(gpu_hits.data(), d_amr_hits, 
                   gpu_hits.size() * sizeof(AMRHit), cudaMemcpyDeviceToHost);
        
        // Update with our modified hits
        size_t hit_idx = 0;
        for (int i = 0; i < current_batch_size && hit_idx < hits.size(); i++) {
            for (uint32_t j = 0; j < MAX_MATCHES_PER_READ && hit_idx < hits.size(); j++) {
                if (hits[hit_idx].read_id == i) {
                    gpu_hits[i * MAX_MATCHES_PER_READ + j] = hits[hit_idx];
                    hit_idx++;
                } else {
                    break;
                }
            }
        }
        
        // Copy back to GPU
        cudaMemcpy(d_amr_hits, gpu_hits.data(), 
                   gpu_hits.size() * sizeof(AMRHit), cudaMemcpyHostToDevice);
    }
}

void AMRDetectionPipeline::integratePairedEndScoring() {
    if (current_pair_alignments.empty()) return;
    
    std::cout << "\nIntegrating paired-end scoring..." << std::endl;
    
    // Create a map of read_id to best gene assignment considering pairs
    std::map<uint32_t, uint32_t> read_to_best_gene;
    std::map<uint32_t, float> read_to_best_score;
    
    for (const auto& pair_align : current_pair_alignments) {
        if (pair_align.isConcordant()) {
            // For concordant pairs, strongly prefer the concordant gene
            uint32_t best_gene = pair_align.concordant_genes[0];
            float best_score = pair_align.concordant_scores[0];
            
            for (size_t i = 1; i < pair_align.concordant_genes.size(); i++) {
                if (pair_align.concordant_scores[i] > best_score) {
                    best_gene = pair_align.concordant_genes[i];
                    best_score = pair_align.concordant_scores[i];
                }
            }
            
            // Apply concordance bonus
            best_score *= config.concordance_bonus;
            
            read_to_best_gene[pair_align.r1_read_idx] = best_gene;
            read_to_best_gene[pair_align.r2_read_idx] = best_gene;
            read_to_best_score[pair_align.r1_read_idx] = best_score;
            read_to_best_score[pair_align.r2_read_idx] = best_score;
            
        } else {
            // For discordant pairs, use individual best hits with penalty
            if (!pair_align.r1_gene_ids.empty()) {
                size_t best_idx = 0;
                for (size_t i = 1; i < pair_align.r1_scores.size(); i++) {
                    if (pair_align.r1_scores[i] > pair_align.r1_scores[best_idx]) {
                        best_idx = i;
                    }
                }
                read_to_best_gene[pair_align.r1_read_idx] = pair_align.r1_gene_ids[best_idx];
                read_to_best_score[pair_align.r1_read_idx] = 
                    pair_align.r1_scores[best_idx] * config.discord_penalty;
            }
            
            if (!pair_align.r2_gene_ids.empty()) {
                size_t best_idx = 0;
                for (size_t i = 1; i < pair_align.r2_scores.size(); i++) {
                    if (pair_align.r2_scores[i] > pair_align.r2_scores[best_idx]) {
                        best_idx = i;
                    }
                }
                read_to_best_gene[pair_align.r2_read_idx] = pair_align.r2_gene_ids[best_idx];
                read_to_best_score[pair_align.r2_read_idx] = 
                    pair_align.r2_scores[best_idx] * config.discord_penalty;
            }
        }
    }
    
    std::cout << "Paired-end scoring applied to " << read_to_best_gene.size() 
              << " reads" << std::endl;
}
