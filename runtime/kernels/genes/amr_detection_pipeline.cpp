// amr_detection_pipeline.cpp
#include "amr_detection_pipeline.h"
#include "amr_detection_kernels.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <cuda_runtime.h>

// Genetic code initialization is now in amr_detection_kernels_wrapper.cu

AMRDetectionPipeline::AMRDetectionPipeline(const AMRDetectionConfig& cfg) 
    : config(cfg), current_batch_size(0) {
    
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
    d_amr_hits = nullptr;
    d_hit_counts = nullptr;
    d_coverage_stats = nullptr;
}

AMRDetectionPipeline::~AMRDetectionPipeline() {
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
    
    // Build bloom filter
    buildBloomFilter();
    
    // Initialize coverage statistics
    h_coverage_stats.resize(amr_db->getNumGenes());
    cudaMemset(d_coverage_stats, 0, 
               amr_db->getNumGenes() * sizeof(AMRCoverageStats));
    
    return true;
}

void AMRDetectionPipeline::allocateGPUMemory() {
    // Allocate for maximum batch size
    size_t max_batch = config.reads_per_batch;
    size_t max_read_len = config.max_read_length;
    
    // Read data
    cudaMalloc(&d_reads, max_batch * max_read_len);
    cudaMalloc(&d_read_offsets, max_batch * sizeof(int));
    cudaMalloc(&d_read_lengths, max_batch * sizeof(int));
    cudaMalloc(&d_read_ids, max_batch * sizeof(uint32_t));
    
    // Minimizers (estimate ~100 minimizers per read)
    size_t max_minimizers = max_batch * 1000;
    cudaMalloc(&d_minimizers, max_minimizers * sizeof(Minimizer));
    cudaMalloc(&d_minimizer_counts, max_batch * sizeof(uint32_t));
    cudaMalloc(&d_minimizer_offsets, max_batch * sizeof(uint32_t));
    
    // Bloom filter
    size_t bloom_words = (config.bloom_filter_size + 63) / 64;
    cudaMalloc(&d_bloom_filter, bloom_words * sizeof(uint64_t));
    cudaMemset(d_bloom_filter, 0, bloom_words * sizeof(uint64_t));
    
    // Hits (max 10 per read)
    cudaMalloc(&d_amr_hits, max_batch * 10 * sizeof(AMRHit));
    cudaMalloc(&d_hit_counts, max_batch * sizeof(uint32_t));
    
    // Coverage statistics
    size_t num_genes = amr_db->getNumGenes();
    cudaMalloc(&d_coverage_stats, num_genes * sizeof(AMRCoverageStats));
    
    // Allocate position counts for each gene
    for (size_t i = 0; i < num_genes; i++) {
        uint16_t gene_length = 1000;  // Max protein length, adjust as needed
        uint32_t* position_counts;
        cudaMalloc(&position_counts, gene_length * sizeof(uint32_t));
        cudaMemset(position_counts, 0, gene_length * sizeof(uint32_t));
        
        // Set the pointer in coverage stats
        cudaMemcpy(&d_coverage_stats[i].position_counts, &position_counts,
                   sizeof(uint32_t*), cudaMemcpyHostToDevice);
    }
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
    if (d_amr_hits) cudaFree(d_amr_hits);
    if (d_hit_counts) cudaFree(d_hit_counts);
    
    // Free position counts in coverage stats
    if (d_coverage_stats) {
        // Note: In practice, you'd need to free individual position_counts arrays
        cudaFree(d_coverage_stats);
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
    // Prepare reads in GPU-friendly format
    std::vector<char> concatenated_reads;
    std::vector<int> offsets;
    std::vector<int> lengths;
    
    for (const auto& read : reads) {
        offsets.push_back(concatenated_reads.size());
        lengths.push_back(read.length());
        concatenated_reads.insert(concatenated_reads.end(), read.begin(), read.end());
    }
    
    // Copy to GPU
    cudaMemcpy(d_reads, concatenated_reads.data(), 
               concatenated_reads.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_read_offsets, offsets.data(), 
               offsets.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_read_lengths, lengths.data(), 
               lengths.size() * sizeof(int), cudaMemcpyHostToDevice);
}

void AMRDetectionPipeline::processBatch(const std::vector<std::string>& reads,
                                       const std::vector<std::string>& read_ids) {
    if (reads.empty()) return;
    
    current_batch_size = reads.size();
    std::cout << "Processing batch of " << current_batch_size << " reads" << std::endl;
    
    // Copy reads to GPU
    copyReadsToGPU(reads);
    
    // Generate minimizers
    generateMinimizers();
    
    // Screen with bloom filter
    screenWithBloomFilter();
    
    // Perform translated alignment
    performTranslatedAlignment();
    
    // Extend alignments using minimizer information
    extendAlignments();
    
    // Update coverage statistics
    calculateCoverageStats();
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
    // Allocate filter results
    bool* d_passes_filter;
    cudaMalloc(&d_passes_filter, current_batch_size * sizeof(bool));
    
    // Launch kernel
    launch_screen_minimizers_kernel(
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        d_bloom_filter,
        config.bloom_filter_size,
        d_passes_filter,
        current_batch_size
    );
    
    cudaDeviceSynchronize();
    
    // Count how many passed
    std::vector<uint8_t> passes_filter(current_batch_size);
    cudaMemcpy(passes_filter.data(), d_passes_filter, 
               current_batch_size * sizeof(bool), cudaMemcpyDeviceToHost);
    
    int passed = std::count(passes_filter.begin(), passes_filter.end(), (uint8_t)1);
    std::cout << "Bloom filter: " << passed << "/" << current_batch_size 
              << " reads passed (" << (100.0 * passed / current_batch_size) << "%)" << std::endl;
    
    // Store for alignment stage
    cudaMemcpy(d_read_ids, d_passes_filter, current_batch_size * sizeof(bool), 
               cudaMemcpyDeviceToDevice);
    
    cudaFree(d_passes_filter);
}

void AMRDetectionPipeline::performTranslatedAlignment() {
    // Get protein database info
    char* d_amr_proteins = amr_db->getGPUProteinSequences();
    AMRGeneEntry* d_gene_entries = amr_db->getGPUGeneEntries();
    uint32_t num_proteins = amr_db->getNumGenes();
    
    // Create protein offset and length arrays
    std::vector<uint32_t> protein_offsets(num_proteins);
    std::vector<uint32_t> protein_lengths(num_proteins);
    std::vector<AMRGeneEntry> gene_entries(num_proteins);
    
    cudaMemcpy(gene_entries.data(), d_gene_entries, 
               num_proteins * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    
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
    
    // Launch alignment kernel
    launch_translated_alignment_kernel(
        d_reads,
        d_read_offsets,
        d_read_lengths,
        (bool*)d_read_ids,  // Reusing as filter results
        d_amr_proteins,
        d_protein_offsets,
        d_protein_lengths,
        d_gene_entries,
        d_amr_hits,
        d_hit_counts,
        current_batch_size,
        num_proteins,
        config
    );
    
    cudaDeviceSynchronize();
    
    // Cleanup
    cudaFree(d_protein_offsets);
    cudaFree(d_protein_lengths);
    
    // Report hits
    std::vector<uint32_t> hit_counts(current_batch_size);
    cudaMemcpy(hit_counts.data(), d_hit_counts, 
               current_batch_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    int total_hits = 0;
    for (auto count : hit_counts) {
        total_hits += count;
    }
    
    std::cout << "Found " << total_hits << " AMR hits" << std::endl;
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
    
    // Calculate total hits
    int total_hits = 0;
    for (auto count : hit_counts) {
        total_hits += count;
    }
    
    if (total_hits > 0) {
        // Copy all hits
        std::vector<AMRHit> batch_hits(current_batch_size * 10);
        cudaMemcpy(batch_hits.data(), d_amr_hits, 
                   batch_hits.size() * sizeof(AMRHit), cudaMemcpyDeviceToHost);
        
        // Extract actual hits
        for (int i = 0; i < current_batch_size; i++) {
            for (uint32_t j = 0; j < hit_counts[i]; j++) {
                all_hits.push_back(batch_hits[i * 10 + j]);
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
        << "ref_start\tref_end\tframe\tcomplete_gene\n";
    
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
            << (hit.is_complete_gene ? "Y" : "N") << "\n";
    }
    
    out.close();
    
    // Write coverage statistics
    std::string coverage_file = output_prefix + "_coverage.tsv";
    out.open(coverage_file);
    
    out << "gene_id\tgene_name\ttotal_reads\tmean_coverage\t"
        << "covered_positions\tgene_length\tcoverage_uniformity\n";
    
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
                << stats.mean_coverage << "\t"
                << stats.covered_positions << "\t"
                << stats.gene_length << "\t"
                << stats.coverage_uniformity << "\n";
        }
    }
    
    out.close();
    
    std::cout << "Results written to " << hits_file << " and " << coverage_file << std::endl;
}

void AMRDetectionPipeline::generateClinicalReport(const std::string& output_file) {
    std::ofstream report(output_file);
    
    report << "=== AMR Detection Clinical Report ===\n\n";
    
    // Get hits and coverage
    auto hits = getAMRHits();
    auto coverage_stats = getCoverageStats();
    
    // Group by drug class
    std::map<std::string, std::vector<AMRHit>> hits_by_class;
    for (const auto& hit : hits) {
        hits_by_class[hit.drug_class].push_back(hit);
    }
    
    // Report by drug class
    for (const auto& [drug_class, class_hits] : hits_by_class) {
        report << "Drug Class: " << drug_class << "\n";
        report << "Number of resistance genes detected: " << class_hits.size() << "\n";
        
        // Find unique genes
        std::set<std::string> unique_genes;
        for (const auto& hit : class_hits) {
            unique_genes.insert(hit.gene_name);
        }
        
        report << "Unique genes: ";
        for (const auto& gene : unique_genes) {
            report << gene << " ";
        }
        report << "\n\n";
    }
    
    // High-confidence complete gene detections
    report << "=== High-Confidence AMR Genes (>95% coverage, >95% identity) ===\n";
    for (const auto& hit : hits) {
        if (hit.coverage >= 0.95f && hit.identity >= 0.95f) {
            report << hit.gene_name << " (" << hit.drug_class << "): "
                   << "Coverage=" << hit.coverage * 100 << "%, "
                   << "Identity=" << hit.identity * 100 << "%\n";
        }
    }
    
    report.close();
    std::cout << "Clinical report written to " << output_file << std::endl;
}
