// runtime/kernels/genes/amr_detection_pipeline_v2.cpp
#include "amr_detection_pipeline_v2.h"
#include "amr_detection_kernels_wrapper.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cuda_runtime.h>

namespace BioGPU {

// AMRResults implementation
void AMRResults::writeReport(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open output file " << filename << std::endl;
        return;
    }
    
    out << "# AMR Detection Report\n";
    out << "# Total sequences processed: " << sequences_processed << "\n";
    out << "# Total bases processed: " << bases_processed << "\n";
    out << "# Processing time: " << processing_time_seconds << " seconds\n";
    out << "# Total AMR hits: " << hits.size() << "\n";
    out << "# High-confidence hits: " << high_confidence_hits.size() << "\n\n";
    
    // Header
    out << "read_id\tgene_name\tdrug_class\tidentity\tcoverage\tframe\t"
        << "ref_start\tref_end\tread_start\tread_end\tmutations\tcomplete_gene\n";
    
    // Write all hits
    for (const auto& hit : hits) {
        out << hit.read_id << "\t"
            << hit.gene_name << "\t"
            << hit.drug_class << "\t"
            << std::fixed << std::setprecision(2) << hit.identity * 100 << "%\t"
            << std::fixed << std::setprecision(2) << hit.coverage * 100 << "%\t"
            << (int)hit.frame << "\t"
            << hit.ref_start << "\t"
            << hit.ref_end << "\t"
            << hit.read_start << "\t"
            << hit.read_end << "\t"
            << (int)hit.num_mutations << "\t"
            << (hit.is_complete_gene ? "Yes" : "No") << "\n";
    }
}

void AMRResults::writeSummary(std::ostream& out) const {
    out << "AMR Detection Summary:\n";
    out << "  Total sequences: " << sequences_processed << "\n";
    out << "  Total AMR hits: " << hits.size() << "\n";
    out << "  High-confidence hits: " << high_confidence_hits.size() << "\n";
    
    // Count unique genes hit
    std::set<std::string> unique_genes;
    std::map<std::string, int> drug_class_counts;
    
    for (const auto& hit : high_confidence_hits) {
        unique_genes.insert(hit.gene_name);
        drug_class_counts[hit.drug_class]++;
    }
    
    out << "  Unique AMR genes detected: " << unique_genes.size() << "\n";
    out << "  Drug classes with resistance:\n";
    for (const auto& [drug_class, count] : drug_class_counts) {
        out << "    " << drug_class << ": " << count << " hits\n";
    }
}

void AMRResults::writeClinicalReport(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open clinical report file " << filename << std::endl;
        return;
    }
    
    out << "=== AMR Detection Clinical Report ===\n\n";
    
    // High-confidence section
    out << "=== High-Confidence AMR Genes (>95% coverage, >95% identity) ===\n\n";
    
    std::map<std::string, std::vector<const AMRHit*>> genes_by_drug_class;
    for (const auto& hit : high_confidence_hits) {
        genes_by_drug_class[hit.drug_class].push_back(&hit);
    }
    
    for (const auto& [drug_class, gene_hits] : genes_by_drug_class) {
        out << "Drug Class: " << drug_class << "\n";
        
        // Get unique genes for this drug class
        std::set<std::string> unique_genes;
        for (const auto* hit : gene_hits) {
            unique_genes.insert(hit->gene_name);
        }
        
        for (const auto& gene : unique_genes) {
            // Count reads supporting this gene
            int supporting_reads = 0;
            float max_identity = 0;
            float max_coverage = 0;
            
            for (const auto* hit : gene_hits) {
                if (hit->gene_name == gene) {
                    supporting_reads++;
                    max_identity = std::max(max_identity, hit->identity);
                    max_coverage = std::max(max_coverage, hit->coverage);
                }
            }
            
            out << "  - " << gene << " (";
            out << supporting_reads << " reads, ";
            out << "max identity: " << std::fixed << std::setprecision(1) 
                << max_identity * 100 << "%, ";
            out << "max coverage: " << std::fixed << std::setprecision(1) 
                << max_coverage * 100 << "%)\n";
        }
        out << "\n";
    }
    
    // Summary statistics
    out << "=== Summary ===\n";
    out << "Total AMR genes detected: " << unique_genes.size() << "\n";
    out << "Total supporting reads: " << high_confidence_hits.size() << "\n";
    out << "Drug classes affected: " << genes_by_drug_class.size() << "\n";
}

// AMRDetectionPipeline implementation
AMRDetectionPipeline::AMRDetectionPipeline(const GenesConfig& config)
    : PipelineBase("AMR Detection", config), genes_config(config) {
    accumulated_results = std::make_unique<AMRResults>();
}

void AMRDetectionPipeline::setDatabasePath(const std::string& nucl_path, const std::string& prot_path) {
    amr_db = std::make_unique<NCBIAMRDatabaseLoader>();
    // Database will be loaded during initialization
}

bool AMRDetectionPipeline::initializePipeline() {
    log("Initializing AMR detection pipeline");
    
    // Load AMR database
    if (!loadAMRDatabase()) {
        logError("Failed to load AMR database");
        return false;
    }
    
    // Estimate memory requirements
    size_t estimated_sequences = config.max_sequences_per_batch;
    size_t estimated_bases = config.max_bases_per_batch;
    size_t required_memory = estimateMemoryRequirements(estimated_sequences, estimated_bases);
    
    // Allocate workspace
    if (!allocateWorkspace(required_memory)) {
        logError("Failed to allocate GPU workspace");
        return false;
    }
    
    // Allocate AMR-specific GPU memory
    size_t offset = 0;
    
    // Minimizers (estimate 1 minimizer per 5 bases)
    max_minimizers_per_batch = estimated_bases / 5;
    d_minimizers = (Minimizer*)((char*)d_workspace + offset);
    offset += max_minimizers_per_batch * sizeof(Minimizer);
    
    d_minimizer_counts = (uint32_t*)((char*)d_workspace + offset);
    offset += estimated_sequences * sizeof(uint32_t);
    
    d_minimizer_offsets = (uint32_t*)((char*)d_workspace + offset);
    offset += (estimated_sequences + 1) * sizeof(uint32_t);
    
    // Hits (estimate up to 10 hits per sequence)
    max_hits_per_batch = estimated_sequences * 10;
    d_amr_hits = (void*)((char*)d_workspace + offset);
    offset += max_hits_per_batch * sizeof(AMRResults::AMRHit);
    
    d_hit_counts = (uint32_t*)((char*)d_workspace + offset);
    offset += estimated_sequences * sizeof(uint32_t);
    
    // Coverage stats (one per AMR gene)
    d_coverage_stats = (uint32_t*)((char*)d_workspace + offset);
    offset += amr_db->getNumGenes() * sizeof(uint32_t) * 4;  // Simplified coverage
    
    // Build bloom filter
    if (!buildBloomFilter()) {
        logError("Failed to build bloom filter");
        return false;
    }
    
    log("AMR detection pipeline initialized successfully");
    return true;
}

void AMRDetectionPipeline::cleanupPipeline() {
    // Free bloom filter
    if (d_bloom_filter) {
        cudaFree(d_bloom_filter);
        d_bloom_filter = nullptr;
    }
    
    // Clear results
    accumulated_results.reset();
    
    // Base class cleanup handles workspace and streams
}

bool AMRDetectionPipeline::processBatchImpl(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) {
    if (!gpu_buffer || gpu_buffer->getCurrentSequences() == 0) {
        return true;  // Nothing to process
    }
    
    log("Processing batch with " + std::to_string(gpu_buffer->getCurrentSequences()) + " sequences");
    
    // Clear previous batch results
    clearBatchResults();
    
    // 1. Generate minimizers
    launchGenerateMinimizersKernel(
        gpu_buffer->getSequences(),
        gpu_buffer->getOffsets(),
        gpu_buffer->getLengths(),
        gpu_buffer->getCurrentSequences(),
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        stream
    );
    
    // 2. Screen with bloom filter
    uint8_t* d_passed_filter = (uint8_t*)((char*)d_workspace + workspace_size - 
                                          gpu_buffer->getCurrentSequences());
    
    launchScreenWithBloomFilterKernel(
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        d_bloom_filter,
        gpu_buffer->getCurrentSequences(),
        d_passed_filter,
        stream
    );
    
    // 3. Perform translated search on sequences that passed bloom filter
    launchTranslatedSearchKernel(
        gpu_buffer->getSequences(),
        gpu_buffer->getOffsets(),
        gpu_buffer->getLengths(),
        d_passed_filter,
        gpu_buffer->getCurrentSequences(),
        amr_db->getProteinSequencesGPU(),
        amr_db->getProteinOffsetsGPU(),
        amr_db->getNumProteins(),
        d_amr_hits,
        d_hit_counts,
        stream
    );
    
    // 4. Process hits and update results
    processHits(stream);
    
    // 5. Update coverage statistics
    updateCoverageStats(stream);
    
    // Synchronize to ensure all operations complete
    cudaStreamSynchronize(stream);
    
    return true;
}

bool AMRDetectionPipeline::loadAMRDatabase() {
    if (!amr_db) {
        amr_db = std::make_unique<NCBIAMRDatabaseLoader>();
    }
    
    log("Loading AMR database from " + genes_config.amr_cds_path + " and " + 
        genes_config.amr_prot_path);
    
    if (!amr_db->loadDatabase(genes_config.amr_cds_path, genes_config.amr_prot_path)) {
        logError("Failed to load AMR database");
        return false;
    }
    
    log("Loaded " + std::to_string(amr_db->getNumGenes()) + " AMR genes with " + 
        std::to_string(amr_db->getNumProteins()) + " protein sequences");
    
    // Transfer database to GPU
    if (!amr_db->transferToGPU()) {
        logError("Failed to transfer AMR database to GPU");
        return false;
    }
    
    return true;
}

bool AMRDetectionPipeline::buildBloomFilter() {
    log("Building bloom filter for k-mer screening");
    
    // Allocate bloom filter
    if (cudaMalloc(&d_bloom_filter, bloom_filter_size) != cudaSuccess) {
        logError("Failed to allocate bloom filter memory");
        return false;
    }
    
    // Initialize bloom filter to zeros
    cudaMemset(d_bloom_filter, 0, bloom_filter_size);
    
    // Build bloom filter from AMR database k-mers
    launch_build_bloom_filter_kernel(
        amr_db->getNucleotideSequencesGPU(),
        amr_db->getNucleotideOffsetsGPU(),
        amr_db->getNucleotideLengthsGPU(),
        amr_db->getNumGenes(),
        d_bloom_filter,
        bloom_filter_size,
        genes_config.kmer_size,
        bloom_k
    );
    
    cudaDeviceSynchronize();
    log("Bloom filter built successfully");
    return true;
}

void AMRDetectionPipeline::launchGenerateMinimizersKernel(
    const char* d_sequences,
    const int* d_offsets,
    const int* d_lengths,
    int num_sequences,
    Minimizer* d_minimizers,
    uint32_t* d_minimizer_counts,
    uint32_t* d_minimizer_offsets,
    cudaStream_t stream) {
    
    // Call the wrapper function from amr_detection_kernels_wrapper.cu
    launch_generate_minimizers_kernel(
        d_sequences,
        d_offsets,
        d_lengths,
        num_sequences,
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        genes_config.minimizer_k,
        genes_config.minimizer_w,
        max_minimizers_per_batch,
        stream
    );
}

void AMRDetectionPipeline::launchScreenWithBloomFilterKernel(
    const Minimizer* d_minimizers,
    const uint32_t* d_minimizer_counts,
    const uint32_t* d_minimizer_offsets,
    const uint64_t* d_bloom_filter,
    int num_sequences,
    uint8_t* d_passed_filter,
    cudaStream_t stream) {
    
    launch_screen_reads_kernel(
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        d_bloom_filter,
        bloom_filter_size,
        bloom_k,
        num_sequences,
        d_passed_filter,
        genes_config.min_exact_matches,
        stream
    );
}

void AMRDetectionPipeline::launchTranslatedSearchKernel(
    const char* d_sequences,
    const int* d_offsets,
    const int* d_lengths,
    const uint8_t* d_passed_filter,
    int num_sequences,
    const char* d_protein_db,
    const int* d_protein_offsets,
    int num_proteins,
    void* d_hits,
    uint32_t* d_hit_counts,
    cudaStream_t stream) {
    
    launch_extend_alignments_kernel(
        d_sequences,
        d_offsets,
        d_lengths,
        d_passed_filter,
        num_sequences,
        d_minimizers,
        d_minimizer_counts,
        d_minimizer_offsets,
        d_protein_db,
        d_protein_offsets,
        num_proteins,
        d_hits,
        d_hit_counts,
        genes_config.max_mismatches,
        genes_config.translate_all_frames,
        stream
    );
}

void AMRDetectionPipeline::processHits(cudaStream_t stream) {
    // Copy hit counts from GPU
    std::vector<uint32_t> hit_counts(config.max_sequences_per_batch);
    cudaMemcpyAsync(hit_counts.data(), d_hit_counts, 
                    config.max_sequences_per_batch * sizeof(uint32_t),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    // Count total hits
    size_t total_hits = 0;
    for (size_t i = 0; i < config.max_sequences_per_batch; i++) {
        total_hits += hit_counts[i];
    }
    
    if (total_hits == 0) {
        return;  // No hits to process
    }
    
    // Copy hits from GPU
    std::vector<AMRResults::AMRHit> batch_hits(total_hits);
    cudaMemcpyAsync(batch_hits.data(), d_amr_hits,
                    total_hits * sizeof(AMRResults::AMRHit),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    // Process hits and add to accumulated results
    for (const auto& hit : batch_hits) {
        // Get gene metadata from database
        auto gene_info = amr_db->getGeneInfo(hit.gene_id);
        
        AMRResults::AMRHit processed_hit = hit;
        processed_hit.gene_name = gene_info.gene_name;
        processed_hit.drug_class = gene_info.drug_class;
        
        accumulated_results->hits.push_back(processed_hit);
        
        // Check if high-confidence hit
        if (hit.identity >= genes_config.min_gene_identity &&
            hit.coverage >= genes_config.min_gene_coverage) {
            accumulated_results->high_confidence_hits.push_back(processed_hit);
        }
    }
    
    log("Processed " + std::to_string(total_hits) + " AMR hits");
}

void AMRDetectionPipeline::updateCoverageStats(cudaStream_t stream) {
    // This would update coverage statistics for each AMR gene
    // For now, this is a placeholder
    // In a full implementation, this would track coverage depth across each gene
}

size_t AMRDetectionPipeline::estimateMemoryRequirements(size_t num_sequences, size_t total_bases) {
    size_t memory_needed = 0;
    
    // Minimizers
    size_t max_minimizers = total_bases / 5;  // Rough estimate
    memory_needed += max_minimizers * sizeof(Minimizer);
    memory_needed += num_sequences * sizeof(uint32_t);  // counts
    memory_needed += (num_sequences + 1) * sizeof(uint32_t);  // offsets
    
    // Hits
    size_t max_hits = num_sequences * 10;  // Up to 10 hits per sequence
    memory_needed += max_hits * sizeof(AMRResults::AMRHit);
    memory_needed += num_sequences * sizeof(uint32_t);  // hit counts
    
    // Coverage stats
    memory_needed += amr_db->getNumGenes() * sizeof(uint32_t) * 4;
    
    // Temporary arrays
    memory_needed += num_sequences * sizeof(uint8_t);  // passed filter
    
    return memory_needed;
}

void AMRDetectionPipeline::clearBatchResults() {
    // Clear GPU arrays
    cudaMemset(d_hit_counts, 0, config.max_sequences_per_batch * sizeof(uint32_t));
    cudaMemset(d_minimizer_counts, 0, config.max_sequences_per_batch * sizeof(uint32_t));
}

} // namespace BioGPU