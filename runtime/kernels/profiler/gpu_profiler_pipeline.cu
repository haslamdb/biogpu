// gpu_profiler_pipeline.cu - GPU profiler with abundance estimation
#include <iostream>
#include <chrono>
#include <iomanip>
#include <map>
#include <algorithm>
#include <fstream>
#include <cuda_runtime.h>
#include "fastq_processing.h"
#include "minimizer_extractor.h"
#include "gpu_kmer_database.h"

using namespace biogpu;

// Structure to hold abundance results
struct AbundanceProfile {
    std::map<uint32_t, uint64_t> taxon_read_counts;
    std::map<uint32_t, double> taxon_relative_abundance;
    uint64_t total_reads = 0;
    uint64_t classified_reads = 0;
    uint64_t unclassified_reads = 0;
};

class GPUProfilerPipeline {
private:
    std::unique_ptr<GPUKmerDatabase> database;
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
    AbundanceProfile abundance;
    
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
    GPUProfilerPipeline(const std::string& db_file, int k = 31, int m = 15) {
        database = std::make_unique<GPUKmerDatabase>();
        database->load_binary(db_file);
        
        extractor = std::make_unique<MinimizerExtractor>(k, m);
        
        // Allocate taxon counting array on GPU
        cudaMalloc(&d_taxon_counts, max_taxon_id * sizeof(uint32_t));
        cudaMemset(d_taxon_counts, 0, max_taxon_id * sizeof(uint32_t));
    }
    
    ~GPUProfilerPipeline() {
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
        
        auto start = std::chrono::high_resolution_clock::now();
        
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
        
        // Step 2: GPU database lookup
        database->lookup_batch_gpu(d_minimizer_hashes, d_minimizer_taxons, all_hashes.size());
        
        auto lookup_time = std::chrono::high_resolution_clock::now();
        
        // Step 3: GPU classification
        int block_size = 256;
        int num_blocks = (num_reads + block_size - 1) / block_size;
        
        classify_reads_kernel<<<num_blocks, block_size>>>(
            d_minimizer_taxons, d_minimizer_offsets, d_minimizer_counts,
            nullptr, d_read_classifications, d_confidence_scores, num_reads
        );
        
        // NEW: Update abundance counts on GPU
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
        
        auto end = std::chrono::high_resolution_clock::now();
        
        // Print timing breakdown
        auto extract_ms = std::chrono::duration_cast<std::chrono::microseconds>(extract_time - start).count();
        auto transfer_ms = std::chrono::duration_cast<std::chrono::microseconds>(transfer_time - extract_time).count();
        auto lookup_ms = std::chrono::duration_cast<std::chrono::microseconds>(lookup_time - transfer_time).count();
        auto classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(classify_time - lookup_time).count();
        auto return_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - classify_time).count();
        
        std::cout << "Batch timing (μs): extract=" << extract_ms 
                  << " transfer=" << transfer_ms
                  << " lookup=" << lookup_ms
                  << " classify=" << classify_ms
                  << " return=" << return_ms
                  << " total=" << extract_ms + transfer_ms + lookup_ms + classify_ms + return_ms
                  << " minimizers=" << all_hashes.size() << "\n";
    }
    
    void process_paired_batch_with_abundance(const std::vector<std::string>& sequences1,
                                            const std::vector<std::string>& sequences2,
                                            std::vector<uint32_t>& classifications,
                                            std::vector<float>& confidences,
                                            bool update_abundance = true) {
        // Process both forward and reverse reads together
        size_t num_pairs = sequences1.size();
        if (num_pairs == 0 || sequences1.size() != sequences2.size()) return;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Extract minimizers from both reads
        auto minimizer_results1 = extractor->extract_minimizers(sequences1);
        auto minimizer_results2 = extractor->extract_minimizers(sequences2);
        
        // Combine minimizers from paired reads
        std::vector<uint64_t> all_hashes;
        std::vector<uint32_t> offsets(num_pairs + 1);
        std::vector<uint32_t> counts(num_pairs);
        
        uint32_t current_offset = 0;
        for (size_t i = 0; i < num_pairs; i++) {
            offsets[i] = current_offset;
            
            // Add minimizers from both reads
            for (const auto& m : minimizer_results1[i]) {
                all_hashes.push_back(m.hash);
            }
            for (const auto& m : minimizer_results2[i]) {
                all_hashes.push_back(m.hash);
            }
            
            counts[i] = minimizer_results1[i].size() + minimizer_results2[i].size();
            current_offset += counts[i];
        }
        offsets[num_pairs] = current_offset;
        
        auto extract_time = std::chrono::high_resolution_clock::now();
        
        // Allocate GPU memory
        allocate_gpu_memory(num_pairs, all_hashes.size());
        
        // Transfer to GPU
        cudaMemcpy(d_minimizer_hashes, all_hashes.data(), 
                   all_hashes.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_offsets, offsets.data(), 
                   offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_minimizer_counts, counts.data(), 
                   counts.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        auto transfer_time = std::chrono::high_resolution_clock::now();
        
        // GPU database lookup
        database->lookup_batch_gpu(d_minimizer_hashes, d_minimizer_taxons, all_hashes.size());
        
        auto lookup_time = std::chrono::high_resolution_clock::now();
        
        // GPU classification
        int block_size = 256;
        int num_blocks = (num_pairs + block_size - 1) / block_size;
        
        classify_reads_kernel<<<num_blocks, block_size>>>(
            d_minimizer_taxons, d_minimizer_offsets, d_minimizer_counts,
            nullptr, d_read_classifications, d_confidence_scores, num_pairs
        );
        
        // Update abundance counts on GPU
        if (update_abundance) {
            update_abundance_kernel<<<num_blocks, block_size>>>(
                d_read_classifications, d_taxon_counts, num_pairs, max_taxon_id
            );
        }
        
        cudaDeviceSynchronize();
        auto classify_time = std::chrono::high_resolution_clock::now();
        
        // Transfer results back
        classifications.resize(num_pairs);
        confidences.resize(num_pairs);
        
        cudaMemcpy(classifications.data(), d_read_classifications, 
                   num_pairs * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(confidences.data(), d_confidence_scores, 
                   num_pairs * sizeof(float), cudaMemcpyDeviceToHost);
        
        // Update host-side abundance tracking
        if (update_abundance) {
            abundance.total_reads += num_pairs;
            for (auto taxon : classifications) {
                if (taxon > 0) {
                    abundance.taxon_read_counts[taxon]++;
                    abundance.classified_reads++;
                } else {
                    abundance.unclassified_reads++;
                }
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        
        auto extract_ms = std::chrono::duration_cast<std::chrono::microseconds>(extract_time - start).count();
        auto transfer_ms = std::chrono::duration_cast<std::chrono::microseconds>(transfer_time - extract_time).count();
        auto lookup_ms = std::chrono::duration_cast<std::chrono::microseconds>(lookup_time - transfer_time).count();
        auto classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(classify_time - lookup_time).count();
        auto return_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - classify_time).count();
        
        std::cout << "Paired batch timing (μs): extract=" << extract_ms 
                  << " transfer=" << transfer_ms
                  << " lookup=" << lookup_ms
                  << " classify=" << classify_ms
                  << " return=" << return_ms
                  << " total=" << extract_ms + transfer_ms + lookup_ms + classify_ms + return_ms
                  << " minimizers=" << all_hashes.size() << "\n";
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
        process_paired_batch_with_abundance(sequences1, sequences2, 
                                          classifications, confidences, true);
    }
    
    // Get current abundance profile from GPU
    AbundanceProfile get_abundance_profile() {
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
        
        return abundance;
    }
    
    // Print abundance report
    void print_abundance_report(std::ostream& out = std::cout) {
        auto profile = get_abundance_profile();
        
        out << "\n=== Abundance Profile ===\n";
        out << "Total reads: " << profile.total_reads << "\n";
        out << "Classified: " << profile.classified_reads 
            << " (" << std::fixed << std::setprecision(2) 
            << (100.0 * profile.classified_reads / profile.total_reads) << "%)\n";
        out << "Unclassified: " << profile.unclassified_reads 
            << " (" << (100.0 * profile.unclassified_reads / profile.total_reads) << "%)\n\n";
        
        // Sort by abundance
        std::vector<std::pair<uint32_t, uint64_t>> sorted_taxa;
        for (const auto& [taxon, count] : profile.taxon_read_counts) {
            sorted_taxa.push_back({taxon, count});
        }
        std::sort(sorted_taxa.begin(), sorted_taxa.end(), 
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        out << "Top Taxa by Read Count:\n";
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
    }
    
    // Save abundance profile to file
    void save_abundance_tsv(const std::string& filename) {
        auto profile = get_abundance_profile();
        std::ofstream out(filename);
        
        out << "taxon_id\ttaxon_name\tread_count\trelative_abundance\n";
        
        for (const auto& [taxon, count] : profile.taxon_read_counts) {
            std::string name = taxon_names.count(taxon) ? 
                taxon_names[taxon] : "Unknown";
            
            out << taxon << "\t" << name << "\t" << count << "\t"
                << std::fixed << std::setprecision(4) 
                << profile.taxon_relative_abundance[taxon] << "\n";
        }
    }
    
    // Reset abundance counts (useful for processing multiple samples)
    void reset_abundance() {
        abundance = AbundanceProfile();
        cudaMemset(d_taxon_counts, 0, max_taxon_id * sizeof(uint32_t));
    }
};

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <database.bin> <reads.fastq> [reads2.fastq] [--output-prefix PREFIX]\n";
        std::cerr << "  Single-end mode: provide one FASTQ file\n";
        std::cerr << "  Paired-end mode: provide two FASTQ files (R1 and R2)\n";
        std::cerr << "  --output-prefix: Prefix for output files (default: microbiome_profile)\n";
        return 1;
    }
    
    std::string db_file = argv[1];
    std::string fastq_file1 = argv[2];
    std::string fastq_file2;
    std::string output_prefix = "microbiome_profile";
    
    // Parse arguments
    bool paired_end = false;
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--output-prefix" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (arg.find(".fastq") != std::string::npos || 
                   arg.find(".fq") != std::string::npos) {
            fastq_file2 = arg;
            paired_end = true;
        }
    }
    
    try {
        // Initialize pipeline
        GPUProfilerPipeline pipeline(db_file);
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (paired_end) {
            std::cout << "Running paired-end profiling...\n";
            std::cout << "R1: " << fastq_file1 << "\n";
            std::cout << "R2: " << fastq_file2 << "\n";
            
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
                
                pipeline.process_paired_batch_with_abundance(
                    batch1.sequences, batch2.sequences, 
                    classifications, confidences, true);
                
                if (++batch_num % 10 == 0) {
                    auto profile = pipeline.get_abundance_profile();
                    std::cout << "Processed " << profile.total_reads << " read pairs, "
                              << profile.classified_reads << " classified ("
                              << std::fixed << std::setprecision(1)
                              << (100.0 * profile.classified_reads / profile.total_reads) 
                              << "%)\n";
                }
            }
        } else {
            std::cout << "Running single-end profiling...\n";
            std::cout << "Input: " << fastq_file1 << "\n";
            
            FastqReader reader(fastq_file1);
            ReadBatch batch;
            
            size_t batch_num = 0;
            while (reader.read_batch(batch, 10000)) {
                std::vector<uint32_t> classifications;
                std::vector<float> confidences;
                
                pipeline.process_batch_with_abundance(
                    batch.sequences, classifications, confidences, true);
                
                if (++batch_num % 10 == 0) {
                    auto profile = pipeline.get_abundance_profile();
                    std::cout << "Processed " << profile.total_reads << " reads, "
                              << profile.classified_reads << " classified ("
                              << std::fixed << std::setprecision(1)
                              << (100.0 * profile.classified_reads / profile.total_reads) 
                              << "%)\n";
                }
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        // Print final results
        pipeline.print_abundance_report();
        
        std::cout << "\nProcessing time: " << duration.count() << " seconds\n";
        
        auto profile = pipeline.get_abundance_profile();
        if (paired_end) {
            std::cout << "Pairs per second: " 
                      << (profile.total_reads * 1.0 / duration.count()) << "\n";
        } else {
            std::cout << "Reads per second: " 
                      << (profile.total_reads * 1.0 / duration.count()) << "\n";
        }
        
        // Save results
        std::string tsv_file = output_prefix + "_abundance.tsv";
        pipeline.save_abundance_tsv(tsv_file);
        std::cout << "\nSaved abundance profile to: " << tsv_file << "\n";
        
        // Save detailed report
        std::string report_file = output_prefix + "_report.txt";
        std::ofstream report(report_file);
        pipeline.print_abundance_report(report);
        std::cout << "Saved detailed report to: " << report_file << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}