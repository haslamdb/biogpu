// hierarchical_profiler_pipeline.cu - Updated profiler with hierarchical database
#include <iostream>
#include <chrono>
#include <iomanip>
#include <map>
#include <algorithm>
#include <fstream>
#include <cuda_runtime.h>
#include "fastq_processing.h"
#include "minimizer_extractor.h"
#include "hierarchical_gpu_database.h"

using namespace biogpu;

// Enhanced abundance profile with hierarchical database support
struct HierarchicalAbundanceProfile {
    std::map<uint32_t, uint64_t> taxon_read_counts;
    std::map<uint32_t, double> taxon_relative_abundance;
    uint64_t total_reads = 0;
    uint64_t classified_reads = 0;
    uint64_t unclassified_reads = 0;
    
    // Database performance metrics
    HierarchicalDBStats db_stats;
    double avg_batch_time_ms = 0.0;
    size_t batches_processed = 0;
};

class HierarchicalProfilerPipeline {
private:
    std::unique_ptr<HierarchicalGPUDatabase> database;
    std::unique_ptr<MinimizerExtractor> extractor;
    
    // GPU memory for batch processing
    uint64_t* d_minimizer_hashes = nullptr;
    uint32_t* d_minimizer_taxons = nullptr;
    uint32_t* d_minimizer_offsets = nullptr;
    uint32_t* d_minimizer_counts = nullptr;
    uint32_t* d_read_classifications = nullptr;
    float* d_confidence_scores = nullptr;
    
    // For abundance counting
    uint32_t* d_taxon_counts = nullptr;
    size_t max_taxon_id = 10000;  // Adjust based on your database
    
    size_t allocated_reads = 0;
    size_t allocated_minimizers = 0;
    
    // Abundance tracking
    HierarchicalAbundanceProfile abundance;
    
    // Taxon name mapping (would be loaded from taxonomy file)
    std::map<uint32_t, std::string> taxon_names = {
        {100, "Escherichia_coli"},
        {101, "Klebsiella_pneumoniae"},
        {102, "Enterococcus_faecalis"},
        {103, "Staphylococcus_aureus"},
        {104, "Streptococcus_pneumoniae"},
        {105, "Haemophilus_influenzae"},
        {106, "Proteus_mirabilis"},
        {107, "Pseudomonas_aeruginosa"},
        {108, "Clostridioides_difficile"},
        // Add more as needed
    };
    
    void allocate_gpu_memory(size_t num_reads, size_t num_minimizers) {
        if (num_reads > allocated_reads) {
            if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
            if (d_minimizer_counts) cudaFree(d_minimizer_counts);
            if (d_read_classifications) cudaFree(d_read_classifications);
            if (d_confidence_scores) cudaFree(d_confidence_scores);
            
            cudaMalloc(&d_minimizer_offsets, (num_reads + 1) * sizeof(uint32_t));
            cudaMalloc(&d_minimizer_counts, num_reads * sizeof(uint32_t));
            cudaMalloc(&d_read_classifications, num_reads * sizeof(uint32_t));
            cudaMalloc(&d_confidence_scores, num_reads * sizeof(float));
            
            allocated_reads = num_reads;
        }
        
        if (num_minimizers > allocated_minimizers) {
            if (d_minimizer_hashes) cudaFree(d_minimizer_hashes);
            if (d_minimizer_taxons) cudaFree(d_minimizer_taxons);
            
            cudaMalloc(&d_minimizer_hashes, num_minimizers * sizeof(uint64_t));
            cudaMalloc(&d_minimizer_taxons, num_minimizers * sizeof(uint32_t));
            
            allocated_minimizers = num_minimizers;
        }
    }
    
public:
    HierarchicalProfilerPipeline(const std::string& db_path, 
                                 size_t max_memory_gb = 10,
                                 int k = 31, int m = 15) {
        
        // Initialize hierarchical database with memory limit
        HierarchicalDBConfig db_config;
        db_config.max_gpu_memory_gb = max_memory_gb;
        db_config.tier_size_mb = 512;  // 512MB tiers
        db_config.cache_tiers = std::max(1UL, max_memory_gb * 2);  // About 2 tiers per GB
        
        database = std::make_unique<HierarchicalGPUDatabase>(db_config);
        database->load_database(db_path);
        
        extractor = std::make_unique<MinimizerExtractor>(k, m);
        
        // Allocate taxon counting array on GPU
        cudaMalloc(&d_taxon_counts, max_taxon_id * sizeof(uint32_t));
        cudaMemset(d_taxon_counts, 0, max_taxon_id * sizeof(uint32_t));
        
        std::cout << "Hierarchical profiler initialized:\n";
        std::cout << "  Database path: " << db_path << "\n";
        std::cout << "  Max GPU memory: " << max_memory_gb << " GB\n";
        std::cout << "  Database tiers: " << database->get_tier_count() << "\n";
    }
    
    ~HierarchicalProfilerPipeline() {
        if (d_minimizer_hashes) cudaFree(d_minimizer_hashes);
        if (d_minimizer_taxons) cudaFree(d_minimizer_taxons);
        if (d_minimizer_offsets) cudaFree(d_minimizer_offsets);
        if (d_minimizer_counts) cudaFree(d_minimizer_counts);
        if (d_read_classifications) cudaFree(d_read_classifications);
        if (d_confidence_scores) cudaFree(d_confidence_scores);
        if (d_taxon_counts) cudaFree(d_taxon_counts);
    }
    
    void process_batch_with_abundance(const std::vector<std::string>& sequences,
                                     std::vector<uint32_t>& classifications,
                                     std::vector<float>& confidences,
                                     bool update_abundance = true) {
        size_t num_reads = sequences.size();
        if (num_reads == 0) return;
        
        auto batch_start = std::chrono::high_resolution_clock::now();
        
        // Step 1: Extract minimizers on GPU
        auto minimizer_results = extractor->extract_minimizers(sequences);
        
        // Prepare data for GPU processing
        std::vector<uint64_t> all_hashes;
        std::vector<uint32_t> offsets(num_reads + 1);
        std::vector<uint32_t> counts(num_reads);
        
        uint32_t current_offset = 0;
        for (size_t i = 0; i < num_reads; i++) {
            offsets[i] = current_offset;
            counts[i] = minimizer_results[i].size();
            
            for (const auto& m : minimizer_results[i]) {
                all_hashes.push_back(m.hash);
            }
            
            current_offset += counts[i];
        }
        offsets[num_reads] = current_offset;
        
        auto extract_time = std::chrono::high_resolution_clock::now();
        
        // Allocate GPU memory
        allocate_gpu_memory(num_reads, all_hashes.size());
        
        // Transfer to GPU
        cudaMemcpy(d_minimizer_hashes, all_hashes.data(), 
                   all_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_offsets, offsets.data(), 
                   offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_counts, counts.data(), 
                   counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        auto transfer_time = std::chrono::high_resolution_clock::now();
        
        // Step 2: Hierarchical database lookup
        database->lookup_batch_gpu(d_minimizer_hashes, d_minimizer_taxons, all_hashes.size());
        
        auto lookup_time = std::chrono::high_resolution_clock::now();
        
        // Step 3: GPU classification
        int block_size = 256;
        int num_blocks = (num_reads + block_size - 1) / block_size;
        
        classify_reads_kernel<<<num_blocks, block_size>>>(
            d_minimizer_taxons, d_minimizer_offsets, d_minimizer_counts,
            nullptr, d_read_classifications, d_confidence_scores, num_reads
        );
        
        // Update abundance counts on GPU
        if (update_abundance) {
            update_abundance_kernel<<<num_blocks, block_size>>>(
                d_read_classifications, d_taxon_counts, num_reads, max_taxon_id
            );
        }
        
        cudaDeviceSynchronize();
        auto classify_time = std::chrono::high_resolution_clock::now();
        
        // Transfer results back
        classifications.resize(num_reads);
        confidences.resize(num_reads);
        
        cudaMemcpy(classifications.data(), d_read_classifications, 
                   num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(confidences.data(), d_confidence_scores, 
                   num_reads * sizeof(float), cudaMemcpyDeviceToHost);
        
        // Update host-side abundance tracking
        if (update_abundance) {
            abundance.total_reads += num_reads;
            for (auto taxon : classifications) {
                if (taxon > 0) {
                    abundance.taxon_read_counts[taxon]++;
                    abundance.classified_reads++;
                } else {
                    abundance.unclassified_reads++;
                }
            }
        }
        
        auto batch_end = std::chrono::high_resolution_clock::now();
        
        // Update batch timing statistics
        auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(batch_end - batch_start);
        abundance.avg_batch_time_ms = (abundance.avg_batch_time_ms * abundance.batches_processed + 
                                      batch_duration.count()) / (abundance.batches_processed + 1);
        abundance.batches_processed++;
        
        // Update database statistics
        abundance.db_stats = database->get_statistics();
        
        // Print detailed timing breakdown
        auto extract_us = std::chrono::duration_cast<std::chrono::microseconds>(extract_time - batch_start).count();
        auto transfer_us = std::chrono::duration_cast<std::chrono::microseconds>(transfer_time - extract_time).count();
        auto lookup_us = std::chrono::duration_cast<std::chrono::microseconds>(lookup_time - transfer_time).count();
        auto classify_us = std::chrono::duration_cast<std::chrono::microseconds>(classify_time - lookup_time).count();
        auto return_us = std::chrono::duration_cast<std::chrono::microseconds>(batch_end - classify_time).count();
        
        std::cout << "Batch " << abundance.batches_processed << " timing (μs): "
                  << "extract=" << extract_us 
                  << " transfer=" << transfer_us
                  << " lookup=" << lookup_us
                  << " classify=" << classify_us
                  << " return=" << return_us
                  << " total=" << batch_duration.count() * 1000
                  << " minimizers=" << all_hashes.size()
                  << " cache_hit=" << std::fixed << std::setprecision(1) 
                  << (abundance.db_stats.get_cache_hit_rate() * 100.0) << "%\n";
    }
    
    // Wrapper methods for backward compatibility
    void process_batch(const std::vector<std::string>& sequences,
                      std::vector<uint32_t>& classifications,
                      std::vector<float>& confidences) {
        process_batch_with_abundance(sequences, classifications, confidences, true);
    }
    
    void process_paired_batch(const std::vector<std::string>& sequences1,
                             const std::vector<std::string>& sequences2,
                             std::vector<uint32_t>& classifications,
                             std::vector<float>& confidences) {
        // Combine paired reads
        std::vector<std::string> combined_sequences;
        combined_sequences.reserve(sequences1.size() + sequences2.size());
        
        for (size_t i = 0; i < sequences1.size(); i++) {
            combined_sequences.push_back(sequences1[i]);
            if (i < sequences2.size()) {
                combined_sequences.push_back(sequences2[i]);
            }
        }
        
        // Process combined sequences
        std::vector<uint32_t> combined_classifications;
        std::vector<float> combined_confidences;
        
        process_batch_with_abundance(combined_sequences, combined_classifications, combined_confidences, true);
        
        // Merge classifications for paired reads (take best classification)
        classifications.resize(sequences1.size());
        confidences.resize(sequences1.size());
        
        for (size_t i = 0; i < sequences1.size(); i++) {
            size_t idx1 = i * 2;
            size_t idx2 = i * 2 + 1;
            
            if (idx2 < combined_classifications.size()) {
                // Take classification with higher confidence
                if (combined_confidences[idx1] >= combined_confidences[idx2]) {
                    classifications[i] = combined_classifications[idx1];
                    confidences[i] = combined_confidences[idx1];
                } else {
                    classifications[i] = combined_classifications[idx2];
                    confidences[i] = combined_confidences[idx2];
                }
            } else {
                classifications[i] = combined_classifications[idx1];
                confidences[i] = combined_confidences[idx1];
            }
        }
    }
    
    // Get current abundance profile from GPU
    HierarchicalAbundanceProfile get_abundance_profile() {
        // Copy GPU counts to host
        std::vector<uint32_t> host_counts(max_taxon_id);
        cudaMemcpy(host_counts.data(), d_taxon_counts, 
                   max_taxon_id * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        // Update abundance profile with GPU counts
        for (size_t i = 0; i < max_taxon_id; i++) {
            if (host_counts[i] > 0) {
                abundance.taxon_read_counts[i] = host_counts[i];
            }
        }
        
        // Calculate relative abundances
        if (abundance.classified_reads > 0) {
            for (const auto& [taxon, count] : abundance.taxon_read_counts) {
                abundance.taxon_relative_abundance[taxon] = 
                    100.0 * count / abundance.classified_reads;
            }
        }
        
        // Update database statistics
        abundance.db_stats = database->get_statistics();
        
        return abundance;
    }
    
    // Enhanced abundance report with hierarchical database stats
    void print_abundance_report(std::ostream& out = std::cout) {
        auto profile = get_abundance_profile();
        
        out << "\n=== Hierarchical Database Profiling Results ===\n";
        out << "Total reads: " << profile.total_reads << "\n";
        out << "Classified: " << profile.classified_reads 
            << " (" << std::fixed << std::setprecision(2) 
            << (100.0 * profile.classified_reads / profile.total_reads) << "%)\n";
        out << "Unclassified: " << profile.unclassified_reads 
            << " (" << (100.0 * profile.unclassified_reads / profile.total_reads) << "%)\n";
        
        out << "\n=== Database Performance ===\n";
        out << "Cache hit rate: " << std::fixed << std::setprecision(1) 
            << (profile.db_stats.get_cache_hit_rate() * 100.0) << "%\n";
        out << "Total lookups: " << profile.db_stats.total_lookups << "\n";
        out << "Tier loads: " << profile.db_stats.tier_loads << "\n";
        out << "Tier evictions: " << profile.db_stats.tier_evictions << "\n";
        out << "Average lookup time: " << profile.db_stats.avg_lookup_time_us << " μs\n";
        out << "Average batch time: " << profile.avg_batch_time_ms << " ms\n";
        out << "Batches processed: " << profile.batches_processed << "\n";
        
        // Sort by abundance
        std::vector<std::pair<uint32_t, uint64_t>> sorted_taxa;
        for (const auto& [taxon, count] : profile.taxon_read_counts) {
            sorted_taxa.push_back({taxon, count});
        }
        std::sort(sorted_taxa.begin(), sorted_taxa.end(), 
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        out << "\n=== Top Taxa by Read Count ===\n";
        out << std::left << std::setw(30) << "Species" 
            << std::right << std::setw(12) << "Reads" 
            << std::setw(10) << "Percent" << "\n";
        out << std::string(52, '-') << "\n";
        
        int shown = 0;
        for (const auto& [taxon, count] : sorted_taxa) {
            if (shown++ >= 20) break;  // Show top 20
            
            std::string name = taxon_names.count(taxon) ? 
                taxon_names[taxon] : "Taxon_" + std::to_string(taxon);
            
            out << std::left << std::setw(30) << name
                << std::right << std::setw(12) << count
                << std::setw(9) << std::fixed << std::setprecision(2) 
                << profile.taxon_relative_abundance[taxon] << "%\n";
        }
        
        if (sorted_taxa.size() > 20) {
            out << "... and " << (sorted_taxa.size() - 20) << " more taxa\n";
        }
        
        // Print database detailed statistics
        out << "\n";
        database->print_statistics();
    }
    
    // Optimize database cache for better performance
    void optimize_cache() {
        std::cout << "Optimizing hierarchical database cache...\n";
        database->preload_frequent_tiers(3);  // Preload top 3 most accessed tiers
    }
    
    // Reset abundance counts (useful for processing multiple samples)
    void reset_abundance() {
        abundance = HierarchicalAbundanceProfile();
        cudaMemset(d_taxon_counts, 0, max_taxon_id * sizeof(uint32_t));
        database->clear_cache();  // Clear database cache for fresh start
    }
    
    // Public getters for accessing private members
    const std::map<uint32_t, std::string>& get_taxon_names() const {
        return taxon_names;
    }
    
    HierarchicalGPUDatabase* get_database() {
        return database.get();
    }
};

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <hierarchical_db_path> <reads.fastq> [reads2.fastq] [options]\n";
        std::cerr << "\nOptions:\n";
        std::cerr << "  --memory <GB>         Maximum GPU memory to use (default: 10)\n";
        std::cerr << "  --output-prefix PREFIX Output prefix for results (default: hierarchical_profile)\n";
        std::cerr << "  --optimize-cache      Preload frequent database tiers\n";
        std::cerr << "\nExample:\n";
        std::cerr << "  " << argv[0] << " pathogen_db data/sample_R1.fastq data/sample_R2.fastq --memory 8\n";
        return 1;
    }
    
    std::string db_path = argv[1];
    std::string fastq_file1 = argv[2];
    std::string fastq_file2;
    std::string output_prefix = "hierarchical_profile";
    size_t max_memory_gb = 10;
    bool optimize_cache = false;
    
    // Parse arguments
    bool paired_end = false;
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--memory" && i + 1 < argc) {
            max_memory_gb = std::stoul(argv[++i]);
        } else if (arg == "--output-prefix" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (arg == "--optimize-cache") {
            optimize_cache = true;
        } else if (arg.find(".fastq") != std::string::npos || 
                   arg.find(".fq") != std::string::npos) {
            fastq_file2 = arg;
            paired_end = true;
        }
    }
    
    try {
        // Initialize pipeline
        HierarchicalProfilerPipeline pipeline(db_path, max_memory_gb);
        
        if (optimize_cache) {
            pipeline.optimize_cache();
        }
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (paired_end) {
            std::cout << "Running paired-end hierarchical profiling...\n";
            std::cout << "R1: " << fastq_file1 << "\n";
            std::cout << "R2: " << fastq_file2 << "\n";
            std::cout << "Max GPU memory: " << max_memory_gb << " GB\n";
            
            FastqReader reader1(fastq_file1);
            FastqReader reader2(fastq_file2);
            ReadBatch batch1, batch2;
            
            size_t batch_num = 0;
            while (reader1.read_batch(batch1, 10000) && reader2.read_batch(batch2, 10000)) {
                if (batch1.size() != batch2.size()) {
                    std::cerr << "Warning: Batch size mismatch between paired files\n";
                    break;
                }
                
                std::vector<uint32_t> classifications;
                std::vector<float> confidences;
                
                pipeline.process_paired_batch(
                    batch1.sequences, batch2.sequences, 
                    classifications, confidences);
                
                if (++batch_num % 5 == 0) {  // More frequent updates for hierarchical profiling
                    auto profile = pipeline.get_abundance_profile();
                    std::cout << "\nProgress: Processed " << profile.total_reads << " read pairs, "
                              << profile.classified_reads << " classified ("
                              << std::fixed << std::setprecision(1)
                              << (100.0 * profile.classified_reads / profile.total_reads) 
                              << "%), cache hit rate: " 
                              << (profile.db_stats.get_cache_hit_rate() * 100.0) << "%\n";
                }
            }
        } else {
            std::cout << "Running single-end hierarchical profiling...\n";
            std::cout << "Input: " << fastq_file1 << "\n";
            std::cout << "Max GPU memory: " << max_memory_gb << " GB\n";
            
            FastqReader reader(fastq_file1);
            ReadBatch batch;
            
            size_t batch_num = 0;
            while (reader.read_batch(batch, 10000)) {
                std::vector<uint32_t> classifications;
                std::vector<float> confidences;
                
                pipeline.process_batch_with_abundance(
                    batch.sequences, classifications, confidences, true);
                
                if (++batch_num % 5 == 0) {  // More frequent updates
                    auto profile = pipeline.get_abundance_profile();
                    std::cout << "\nProgress: Processed " << profile.total_reads << " reads, "
                              << profile.classified_reads << " classified ("
                              << std::fixed << std::setprecision(1)
                              << (100.0 * profile.classified_reads / profile.total_reads) 
                              << "%), cache hit rate: " 
                              << (profile.db_stats.get_cache_hit_rate() * 100.0) << "%\n";
                }
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Print final results
        pipeline.print_abundance_report();
        
        std::cout << "\n=== Performance Summary ===\n";
        std::cout << "Total processing time: " << duration.count() << " seconds\n";
        
        auto profile = pipeline.get_abundance_profile();
        if (paired_end) {
            std::cout << "Read pairs per second: " 
                      << (profile.total_reads * 1.0 / duration.count()) << "\n";
        } else {
            std::cout << "Reads per second: " 
                      << (profile.total_reads * 1.0 / duration.count()) << "\n";
        }
        
        std::cout << "Average batch processing time: " << profile.avg_batch_time_ms << " ms\n";
        std::cout << "Database cache efficiency: " 
                  << std::fixed << std::setprecision(1) 
                  << (profile.db_stats.get_cache_hit_rate() * 100.0) << "%\n";
        
        // Calculate effective throughput considering cache misses
        double effective_throughput = profile.total_reads / duration.count();
        double cache_penalty = (1.0 - profile.db_stats.get_cache_hit_rate()) * 0.1; // 10% penalty per miss
        double optimized_throughput = effective_throughput * (1.0 - cache_penalty);
        
        std::cout << "Effective throughput (accounting for cache): " 
                  << optimized_throughput << " reads/sec\n";
        
        // Memory usage recommendations
        auto db_stats = profile.db_stats;
        if (db_stats.get_cache_hit_rate() < 0.8) {
            std::cout << "\n⚠️  Low cache hit rate detected (" 
                      << (db_stats.get_cache_hit_rate() * 100.0) << "%)\n";
            std::cout << "Consider increasing --memory to improve performance\n";
            std::cout << "Recommended: --memory " << (max_memory_gb * 2) << "\n";
        } else if (db_stats.get_cache_hit_rate() > 0.95) {
            std::cout << "\n✅ Excellent cache performance (" 
                      << (db_stats.get_cache_hit_rate() * 100.0) << "%)\n";
            std::cout << "Current memory allocation is optimal\n";
        }
        
        // Save results
        std::string tsv_file = output_prefix + "_abundance.tsv";
        std::ofstream tsv_out(tsv_file);
        tsv_out << "taxon_id\ttaxon_name\tread_count\trelative_abundance\n";
        
        for (const auto& [taxon, count] : profile.taxon_read_counts) {
            const auto& taxon_names = pipeline.get_taxon_names();
            std::string name = taxon_names.count(taxon) ? 
                taxon_names.at(taxon) : "Unknown";
            
            tsv_out << taxon << "\t" << name << "\t" << count << "\t"
                    << std::fixed << std::setprecision(4) 
                    << profile.taxon_relative_abundance[taxon] << "\n";
        }
        tsv_out.close();
        std::cout << "\nSaved abundance profile to: " << tsv_file << "\n";
        
        // Save detailed report with performance metrics
        std::string report_file = output_prefix + "_detailed_report.txt";
        std::ofstream report(report_file);
        pipeline.print_abundance_report(report);
        
        // Add performance section to report
        report << "\n=== Performance Analysis ===\n";
        report << "Hardware: " << max_memory_gb << " GB GPU memory limit\n";
        report << "Database tiers: " << pipeline.get_database()->get_tier_count() << "\n";
        report << "Processing mode: " << (paired_end ? "Paired-end" : "Single-end") << "\n";
        report << "Batches processed: " << profile.batches_processed << "\n";
        report << "Tier loads during processing: " << db_stats.tier_loads << "\n";
        report << "Tier evictions during processing: " << db_stats.tier_evictions << "\n";
        
        if (db_stats.tier_evictions > 0) {
            report << "\nMemory pressure detected - " << db_stats.tier_evictions 
                   << " tier evictions occurred.\n";
            report << "Consider increasing memory allocation for better performance.\n";
        }
        
        report.close();
        std::cout << "Saved detailed report to: " << report_file << "\n";
        
        // Save performance metrics for analysis
        std::string perf_file = output_prefix + "_performance.json";
        std::ofstream perf(perf_file);
        perf << "{\n";
        perf << "  \"total_reads\": " << profile.total_reads << ",\n";
        perf << "  \"processing_time_seconds\": " << duration.count() << ",\n";
        perf << "  \"reads_per_second\": " << (profile.total_reads * 1.0 / duration.count()) << ",\n";
        perf << "  \"cache_hit_rate\": " << profile.db_stats.get_cache_hit_rate() << ",\n";
        perf << "  \"average_batch_time_ms\": " << profile.avg_batch_time_ms << ",\n";
        perf << "  \"tier_loads\": " << db_stats.tier_loads << ",\n";
        perf << "  \"tier_evictions\": " << db_stats.tier_evictions << ",\n";
        perf << "  \"memory_limit_gb\": " << max_memory_gb << ",\n";
        perf << "  \"database_tiers\": " << pipeline.get_database()->get_tier_count() << ",\n";
        perf << "  \"classification_rate\": " << (100.0 * profile.classified_reads / profile.total_reads) << "\n";
        perf << "}\n";
        perf.close();
        std::cout << "Saved performance metrics to: " << perf_file << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}