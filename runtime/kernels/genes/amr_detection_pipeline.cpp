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
    uint32_t species_id;
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
      engine_capacity((cfg.reads_per_batch * 2) + (cfg.reads_per_batch / 5)) {
    
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

void AMRDetectionPipeline::validateGeneMappings() {
    std::cout << "\n=== Gene Mapping Validation ===" << std::endl;
    
    // Get first 10 gene entries
    std::vector<AMRGeneEntry> entries = getGeneEntries();
    std::cout << "First 10 gene entries from database:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(10), entries.size()); i++) {
        std::cout << "  [" << i << "] gene_name: " << entries[i].gene_name
                  << ", gene_family: " << entries[i].gene_family
                  << ", class: " << entries[i].class_ << std::endl;
    }
    
    // Check for common gene families
    std::map<std::string, int> family_counts;
    for (const auto& entry : entries) {
        family_counts[entry.gene_family]++;
    }
    
    std::cout << "\nGene family distribution:" << std::endl;
    for (const auto& [family, count] : family_counts) {
        if (count > 1) {  // Only show families with multiple genes
            std::cout << "  " << family << ": " << count << " genes" << std::endl;
        }
    }
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
    
    // Validate gene mappings
    validateGeneMappings();
    
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
    std::cout << "Creating search engine for batch size: " << engine_capacity << std::endl;
    std::cout << "This should match max GPU allocation: " << ((config.reads_per_batch * 2) + (config.reads_per_batch / 5)) << std::endl;
    
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
    
    return true;
}

void AMRDetectionPipeline::allocateGPUMemory() {
    // Allocate for maximum batch size
    // For paired-end reads, we process R1 and R2 separately, so need 2x batch size
    // Add 10% safety margin to avoid boundary issues
    size_t max_batch = (config.reads_per_batch * 2) + (config.reads_per_batch / 5);  // 220,000 instead of 200,000
    size_t max_read_len = config.max_read_length;
    
    std::cout << "=== MEMORY ALLOCATION DEBUG ===" << std::endl;
    std::cout << "config.reads_per_batch: " << config.reads_per_batch << std::endl;
    std::cout << "max_batch calculated: " << max_batch << std::endl;
    std::cout << "Allocating for reads: " << (max_batch * max_read_len) << " bytes" << std::endl;
    std::cout << "Allocating read_offsets: " << (max_batch * sizeof(int)) << " bytes" << std::endl;
    std::cout << "Allocating read_lengths: " << (max_batch * sizeof(int)) << " bytes" << std::endl;
    
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
    // Debug input reads
    std::cout << "DEBUG copyReadsToGPU: Processing " << reads.size() << " reads" << std::endl;
    for (int i = 0; i < std::min(5, (int)reads.size()); i++) {
        std::cout << "  Read " << i << ": length=" << reads[i].length() << std::endl;
    }
    
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
    
    // Debug built arrays
    std::cout << "DEBUG: Built " << offsets.size() << " offsets, " << lengths.size() << " lengths" << std::endl;
    for (int i = 0; i < std::min(5, (int)lengths.size()); i++) {
        std::cout << "  Length[" << i << "]=" << lengths[i] << std::endl;
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
    // At the start of processBatch():
    std::cout << "\n### processBatch ENTRY ###" << std::endl;
    std::cout << "current_batch_size before processing: " << current_batch_size << std::endl;
    std::cout << "total_reads_processed before: " << total_reads_processed << std::endl;
    std::cout << "Input reads.size(): " << reads.size() << std::endl;
    
    // At the start of processBatch(), add:
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "WARNING: Clearing existing CUDA error: " << cudaGetErrorString(err) << std::endl;
        // Clear the error but continue
    }
    
    if (reads.empty()) return;
    
    // Reset current_batch_size
    current_batch_size = reads.size();
    std::cout << "Set current_batch_size to: " << current_batch_size << std::endl;
    
    total_reads_processed += current_batch_size;
    std::cout << "Processing batch of " << current_batch_size << " reads" << std::endl;
    
    // Copy reads to GPU
    copyReadsToGPU(reads);
    
    // After copyReadsToGPU(), add:
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "ERROR: CUDA error after copyReadsToGPU: " << cudaGetErrorString(err) << std::endl;
        return; // Don't continue with corrupted state
    }
    
    // Debug: Verify the offsets are correct
    std::vector<int> h_offsets(current_batch_size);
    std::vector<int> h_lengths(current_batch_size);
    cudaMemcpy(h_offsets.data(), d_read_offsets, current_batch_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_lengths.data(), d_read_lengths, current_batch_size * sizeof(int), cudaMemcpyDeviceToHost);
    
    std::cout << "First 5 read offsets and lengths:" << std::endl;
    for (int i = 0; i < std::min(5, current_batch_size); i++) {
        std::cout << "  Read " << i << ": offset=" << h_offsets[i] << ", length=" << h_lengths[i] << std::endl;
    }
    
    // Check if offsets are monotonically increasing
    for (int i = 1; i < std::min(100, current_batch_size); i++) {
        if (h_offsets[i] <= h_offsets[i-1]) {
            std::cerr << "ERROR: Non-monotonic offsets at index " << i 
                      << ": " << h_offsets[i-1] << " -> " << h_offsets[i] << std::endl;
            break;
        }
    }
    
    // Skip bloom filter completely if disabled (default behavior)
    if (!config.use_bloom_filter) {
        // Don't allocate or use any bloom filter memory
        // Pass nullptr to indicate all reads should be processed
        performTranslatedAlignment();
        
        // Skip per-batch statistics to avoid memory corruption
        // Statistics will be collected at the end for all batches
        std::cout << "Batch completed successfully" << std::endl;
        
        // Skip per-batch calculations to avoid memory corruption
        // These will be done once at the end for all batches
        return;
    }
    
    // Only do bloom filter operations if explicitly enabled
    generateMinimizers();
    screenWithBloomFilter();
    performTranslatedAlignment();
    
    // Skip per-batch statistics to avoid memory corruption
    // Statistics will be collected at the end for all batches
    std::cout << "Batch completed successfully" << std::endl;
    
    // Skip per-batch calculations to avoid memory corruption
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
    
    std::cout << "DEBUG: Processing " << current_batch_size << " reads with " 
              << total_bases << " total bases" << std::endl;
    
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
    
    std::cout << "Allocating protein matches: " << protein_match_size << " bytes" << std::endl;
    std::cout << "  ProteinMatch size: " << sizeof(ProteinMatch) << " bytes" << std::endl;
    std::cout << "  Max matches per read: " << MAX_MATCHES_PER_READ << std::endl;
    
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
    
    // Check GPU memory before kernel call
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "GPU memory before translated search: " << free_mem << "/" << total_mem << " bytes free" << std::endl;
    
    std::cout << "=== Translated Search Debug ===" << std::endl;
    std::cout << "Current batch size: " << current_batch_size << std::endl;
    std::cout << "Engine capacity: " << engine_capacity << std::endl;
    std::cout << "Using bloom filter: " << (config.use_bloom_filter ? "YES" : "NO") << std::endl;
    
    // Use the class member search engine
    // Pass bloom filter results if enabled, otherwise nullptr to process all reads
    bool* reads_filter = config.use_bloom_filter ? d_read_passes_filter : (bool*)nullptr;
    
    // Debug: Print all parameters being passed
    std::cout << "\n=== DEBUG: search_translated_reads parameters ===" << std::endl;
    std::cout << "translated_search_engine: " << translated_search_engine << std::endl;
    std::cout << "d_reads: " << (void*)d_reads << std::endl;
    std::cout << "d_read_lengths: " << (void*)d_read_lengths << std::endl;
    std::cout << "d_read_offsets: " << (void*)d_read_offsets << std::endl;
    std::cout << "bloom filter: " << (reads_filter ? (void*)reads_filter : nullptr) << std::endl;
    std::cout << "current_batch_size: " << current_batch_size << std::endl;
    std::cout << "d_protein_matches: " << (void*)d_protein_matches << std::endl;
    std::cout << "d_hit_counts: " << (void*)d_hit_counts << std::endl;
    std::cout << "=== END DEBUG ===" << std::endl;
    
    // Debug: Verify GPU allocations are valid
    std::cout << "\n=== DEBUG: Verifying GPU memory before search ===" << std::endl;
    
    // Test if we can read from the GPU pointers
    int test_length;
    cudaError_t err = cudaMemcpy(&test_length, d_read_lengths, sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cerr << "ERROR: Cannot read from d_read_lengths: " << cudaGetErrorString(err) << std::endl;
    } else {
        std::cout << "First read length: " << test_length << std::endl;
    }
    
    
    // Clear any existing CUDA errors before the call
    cudaError_t pre_err = cudaGetLastError();
    if (pre_err != cudaSuccess) {
        std::cerr << "WARNING: Existing CUDA error before search: " << cudaGetErrorString(pre_err) << std::endl;
    }
    
    std::cout << "Calling search_translated_reads..." << std::endl;
    
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
            
            // Debug output for first few reads with hits
            if (i < 5 && hit_counts[i] > 0) {
                for (uint32_t j = 0; j < std::min(hit_counts[i], 3u); j++) {
                    ProteinMatch& pm = protein_matches[i * MAX_MATCHES_PER_READ + j];
                    std::cout << "DEBUG: Read " << i << " Match " << j << ":" << std::endl;
                    std::cout << "  protein_id=" << pm.protein_id << ", gene_id=" << pm.gene_id 
                              << " (expecting gene_id from FASTA header)" << std::endl;
                    std::cout << "  identity=" << pm.identity << ", match_length=" << pm.match_length << std::endl;
                    std::cout << "  ref_start=" << pm.ref_start << ", query_start=" << pm.query_start << std::endl;
                    std::cout << "  alignment_score=" << pm.alignment_score << std::endl;
                    
                    // Check what gene this is
                    if (pm.gene_id < gene_entries.size()) {
                        std::cout << "  gene_name=" << gene_entries[pm.gene_id].gene_name 
                                  << ", gene_family=" << gene_entries[pm.gene_id].gene_family
                                  << ", drug_class=" << gene_entries[pm.gene_id].class_ << std::endl;
                    }
                }
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
                    printf("Read %d %s: Gene %d (Species %d) - Score: %.2f %s\n",
                        pairInfo->read_idx,
                        pairInfo->is_read2 ? "R2" : "R1",
                        pm.gene_id,
                        pm.species_id,
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

void AMRDetectionPipeline::calculateCoverageStats() {
    uint32_t num_genes = amr_db->getNumGenes();
    
    // Update coverage with current batch
    launch_update_coverage_stats_kernel(
        d_amr_hits,
        d_hit_counts,
        d_coverage_stats,
        current_batch_size,
        num_genes
    );
    
    cudaDeviceSynchronize();
    
    // Finalize statistics
    AMRGeneEntry* d_gene_entries = amr_db->getGPUGeneEntries();
    
    launch_finalize_coverage_stats_kernel(
        d_coverage_stats,
        d_gene_entries,
        num_genes
    );
    
    cudaDeviceSynchronize();
}

std::vector<AMRHit> AMRDetectionPipeline::getAMRHits() {
    std::vector<AMRHit> all_hits;
    std::vector<uint32_t> hit_counts(current_batch_size);
    
    cudaMemcpy(hit_counts.data(), d_hit_counts, 
               current_batch_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    // Debug: Print hit counts
    int total_hits = 0;
    for (auto count : hit_counts) {
        total_hits += count;
    }
    std::cout << "DEBUG getAMRHits: current_batch_size=" << current_batch_size 
              << ", total_hits=" << total_hits << std::endl;
    
    if (total_hits > 0) {
        // Copy all hits - use MAX_MATCHES_PER_READ instead of hardcoded 10
        std::vector<AMRHit> batch_hits(current_batch_size * MAX_MATCHES_PER_READ);
        cudaMemcpy(batch_hits.data(), d_amr_hits, 
                   batch_hits.size() * sizeof(AMRHit), cudaMemcpyDeviceToHost);
        
        // Extract actual hits and debug first few
        for (int i = 0; i < current_batch_size; i++) {
            for (uint32_t j = 0; j < hit_counts[i]; j++) {
                AMRHit& hit = batch_hits[i * MAX_MATCHES_PER_READ + j];
                if (all_hits.size() < 5) {  // Debug first 5 hits
                    std::cout << "  Hit " << all_hits.size() << ": gene_id=" << hit.gene_id 
                              << ", identity=" << hit.identity << std::endl;
                }
                all_hits.push_back(hit);
            }
        }
    }
    
    return all_hits;
}

std::vector<AMRCoverageStats> AMRDetectionPipeline::getCoverageStats() {
    uint32_t num_genes = amr_db->getNumGenes();
    std::vector<AMRCoverageStats> stats(num_genes);
    
    cudaMemcpy(stats.data(), d_coverage_stats, 
               num_genes * sizeof(AMRCoverageStats), cudaMemcpyDeviceToHost);
    
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
            
            // Calculate mean depth
            if (stats.covered_positions > 0) {
                stats.mean_depth = (float)stats.total_bases_mapped / stats.covered_positions;
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
    // At the very start of processBatchPaired():
    std::cout << "### RESETTING PIPELINE STATE ###" << std::endl;
    
    // Clear paired read info from previous batch
    paired_read_info.clear();
    
    // Reset batch size 
    current_batch_size = 0;
    
    std::cout << "Pipeline state reset completed" << std::endl;
    
    std::cout << "\n### processBatchPaired ENTRY ###" << std::endl;
    std::cout << "Input validation:" << std::endl;
    std::cout << "  reads1.size(): " << reads1.size() << std::endl;
    std::cout << "  reads2.size(): " << reads2.size() << std::endl;
    std::cout << "  read_ids.size(): " << read_ids.size() << std::endl;

    // Check for empty reads in input
    int empty_r1 = 0, empty_r2 = 0;
    for (const auto& read : reads1) {
        if (read.empty()) empty_r1++;
    }
    for (const auto& read : reads2) {
        if (read.empty()) empty_r2++;
    }

    if (empty_r1 > 0 || empty_r2 > 0) {
        std::cout << "WARNING: Found " << empty_r1 << " empty R1 and " 
                  << empty_r2 << " empty R2 reads in input!" << std::endl;
    }
    
    if (reads1.empty() && reads2.empty()) return;
    
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
    
    // Debug output after building all_reads vector
    std::cout << "DEBUG processBatchPaired: Built " << all_reads.size() << " total reads from " 
              << num_pairs << " pairs" << std::endl;
              
    // Check for empty reads
    int empty_count = 0;
    for (const auto& read : all_reads) {
        if (read.empty()) empty_count++;
    }
    if (empty_count > 0) {
        std::cout << "WARNING: Found " << empty_count << " empty reads in all_reads vector!" << std::endl;
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
    
    // Apply concordance bonus
    const float CONCORDANCE_BONUS = 1.5f;  // 50% score boost for concordant pairs
    const float DISCORD_PENALTY = 0.8f;    // 20% penalty for discordant pairs
    
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
