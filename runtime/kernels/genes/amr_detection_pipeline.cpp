// amr_detection_pipeline.cpp
#include "amr_detection_pipeline.h"
#include "amr_detection_kernels.h"
#include "translated_search_amr.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <cstdio>
#include <iomanip>
#include <cuda_runtime.h>

// Genetic code initialization is now in amr_detection_kernels_wrapper.cu

AMRDetectionPipeline::AMRDetectionPipeline(const AMRDetectionConfig& cfg) 
    : config(cfg), current_batch_size(0), total_reads_processed(0) {
    
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
    initializeCoverageStats();
    
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
    total_reads_processed += current_batch_size;
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
    
    // Calculate abundance metrics (RPKM/TPM)
    calculateAbundanceMetrics();
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
    // Create translated search engine with Smith-Waterman enabled
    void* search_engine = create_translated_search_engine_with_sw(current_batch_size, true);
    
    // Use pre-built protein database directory
    // This should be created by running: ./build_amr_protein_db AMRProt.fa amr_protein_db/
    std::string protein_db_path = config.protein_db_path;
    if (protein_db_path.empty()) {
        // Default location
        protein_db_path = "amr_protein_db";
    }
    
    // Load protein database into search engine
    if (load_protein_database(search_engine, protein_db_path.c_str()) != 0) {
        std::cerr << "Failed to load protein database from " << protein_db_path << std::endl;
        std::cerr << "Please ensure you have run: ./build_amr_protein_db AMRProt.fa " << protein_db_path << std::endl;
        destroy_translated_search_engine(search_engine);
        return;
    }
    
    // Load gene entries from database if not already loaded
    if (gene_entries.empty()) {
        uint32_t num_genes = amr_db->getNumGenes();
        gene_entries.resize(num_genes);
        cudaMemcpy(gene_entries.data(), amr_db->getGPUGeneEntries(),
                   num_genes * sizeof(AMRGeneEntry), cudaMemcpyDeviceToHost);
    }
    
    // Structure to match ProteinMatch from translated_search_amr.cu
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
    
    // Allocate space for results
    ProteinMatch* d_protein_matches;
    cudaMalloc(&d_protein_matches, current_batch_size * 10 * sizeof(ProteinMatch));
    
    // Use the search engine
    search_translated_reads(search_engine, d_reads, d_read_lengths, 
                          d_read_offsets, (bool*)d_read_ids, 
                          current_batch_size, d_protein_matches, d_hit_counts);
    
    cudaDeviceSynchronize();
    
    // Convert ProteinMatch results to AMRHit format
    std::vector<uint32_t> hit_counts(current_batch_size);
    cudaMemcpy(hit_counts.data(), d_hit_counts, 
               current_batch_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    
    std::vector<ProteinMatch> protein_matches(current_batch_size * 10);
    cudaMemcpy(protein_matches.data(), d_protein_matches,
               protein_matches.size() * sizeof(ProteinMatch), cudaMemcpyDeviceToHost);
    
    // Convert to AMRHit format
    std::vector<AMRHit> amr_hits;
    for (int i = 0; i < current_batch_size; i++) {
        for (uint32_t j = 0; j < hit_counts[i]; j++) {
            ProteinMatch& pm = protein_matches[i * 10 + j];
            AMRHit hit = {};
            
            hit.read_id = pm.read_id;
            hit.gene_id = pm.protein_id;  // Using protein_id as gene_id
            hit.ref_start = pm.ref_start;
            hit.ref_end = pm.ref_start + pm.match_length;
            hit.read_start = pm.query_start;
            hit.read_end = pm.query_start + pm.match_length;
            hit.identity = pm.identity;
            hit.coverage = (float)pm.match_length / gene_entries[pm.protein_id].protein_length;
            hit.frame = pm.frame;
            hit.num_mutations = pm.num_mutations;
            hit.is_complete_gene = (pm.ref_start == 0 && 
                                   hit.ref_end >= gene_entries[pm.protein_id].protein_length * 0.95);
            
            strncpy(hit.gene_name, gene_entries[pm.protein_id].gene_name, 63);
            strncpy(hit.drug_class, gene_entries[pm.protein_id].class_, 31);
            
            amr_hits.push_back(hit);
        }
    }
    
    // Copy converted hits back to GPU
    if (!amr_hits.empty()) {
        cudaMemcpy(d_amr_hits, amr_hits.data(), 
                   amr_hits.size() * sizeof(AMRHit), cudaMemcpyHostToDevice);
    }
    
    // Cleanup
    cudaFree(d_protein_matches);
    destroy_translated_search_engine(search_engine);
    
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
