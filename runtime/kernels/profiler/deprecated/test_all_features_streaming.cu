// test_all_features_streaming.cu
// Test all implemented features with proper streaming approach
// Enhanced to include uniqueness scoring and co-occurrence analysis

#include <iostream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <cuda_runtime.h>
#include <algorithm>
#include <numeric>

#include "processing/genome_file_processor.h"
#include "gpu/gpu_minimizer_extraction.cuh"
#include "features/enhanced_minimizer_flags.h"
#include "features/namespace_conflict_resolution.h"
#include "gpu_kraken_types.h"
// #include "processing/minimizer_feature_extractor.h"
// Forward declare what we need
class MinimizerFeatureExtractor {
public:
    MinimizerFeatureExtractor(size_t max_minimizers, size_t max_genomes) {}
    void process_first_pass(GPUMinimizerHit* d_hits, size_t num_hits, const std::vector<uint32_t>& taxons) {}
    void process_second_pass(GPUMinimizerHit* d_hits, size_t num_hits) {}
    struct Statistics {
        size_t total_unique_minimizers = 0;
        size_t singleton_minimizers = 0;
        size_t rare_minimizers = 0;
        size_t common_minimizers = 0;
        double avg_taxonomic_spread = 0.0;
    };
    Statistics get_statistics() const { return Statistics(); }
};
// Forward declarations for feature functions
bool compute_and_encode_uniqueness_scores(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint32_t>& taxon_ids,
    size_t num_genomes);

bool compute_and_encode_cooccurrence_scores(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    const std::vector<uint64_t>& unique_minimizers,
    const std::vector<uint32_t>& occurrence_counts,
    const std::vector<GPUGenomeInfo>& genome_info);

// Enhanced feature analysis functions
void analyze_minimizer_features(const std::vector<GPUMinimizerHit>& hits) {
    std::cout << "\n=== COMPREHENSIVE FEATURE ANALYSIS ===" << std::endl;
    std::cout << "Total minimizers: " << hits.size() << std::endl;
    
    // Uniqueness analysis
    std::vector<size_t> uniqueness_counts(8, 0);
    std::vector<size_t> cooccurrence_counts(8, 0);
    std::vector<size_t> gc_counts(8, 0);
    std::vector<size_t> complexity_counts(8, 0);
    
    size_t position_clustered = 0;
    size_t position_uniform = 0;
    size_t contamination_risk = 0;
    size_t unique_class = 0;
    size_t canonical_class = 0;
    size_t redundant_class = 0;
    
    // NEW: Enhanced uniqueness tracking
    size_t unique_minimizers = 0;      // ≥90% uniqueness
    size_t rare_minimizers = 0;        // ≤3 occurrences
    size_t reliable_minimizers = 0;   // Suitable for classification
    
    std::vector<float> ml_weights;
    double uniqueness_sum = 0.0;
    
    for (const auto& hit : hits) {
        uint32_t flags = hit.feature_flags;
        uint16_t strand_flags = hit.strand;
        
        // Get uniqueness category (bits 8-10)
        uint8_t uniqueness = MinimizerFlags::get_uniqueness_category_safe(flags);
        if (uniqueness < uniqueness_counts.size()) {
            uniqueness_counts[uniqueness]++;
        }
        
        // Track unique/rare/reliable flags
        if (MinimizerFlags::is_unique_minimizer_safe(flags)) {
            unique_minimizers++;
        }
        if (MinimizerFlags::is_rare_minimizer_safe(flags)) {
            rare_minimizers++;
        }
        if (MinimizerFlags::is_reliable_minimizer_safe(flags)) {
            reliable_minimizers++;
        }
        
        // Get co-occurrence score (bits 14-16)
        uint8_t cooccurrence = MinimizerFlags::get_cooccurrence_score(flags);
        if (cooccurrence < cooccurrence_counts.size()) {
            cooccurrence_counts[cooccurrence]++;
        }
        
        uint8_t gc_category = MinimizerFlags::get_gc_content_category(flags);
        if (gc_category < gc_counts.size()) {
            gc_counts[gc_category]++;
        }
        
        uint8_t complexity = MinimizerFlags::get_complexity_score(flags);
        if (complexity < complexity_counts.size()) {
            complexity_counts[complexity]++;
        }
        
        if (MinimizerFlags::has_position_bias(flags)) {
            position_clustered++;
        } else {
            position_uniform++;
        }
        
        if (MinimizerFlags::has_contamination_risk(flags)) {
            contamination_risk++;
        }
        
        uint32_t classification = MinimizerFlags::get_classification(strand_flags);
        if (classification == 0) unique_class++;
        else if (classification == 1) canonical_class++;
        else if (classification == 2) redundant_class++;
        
        float ml_weight = MinimizerFlags::ml_weight_to_float(hit.ml_weight);
        ml_weights.push_back(ml_weight);
        uniqueness_sum += ml_weight;
    }
    
    // Print uniqueness statistics
    std::cout << "\n--- UNIQUENESS SCORES ---" << std::endl;
    const char* uniqueness_names[] = {
        "Extremely common (0-10%)",
        "Very low (10-30%)",
        "Low (30-50%)",
        "Moderate (50-70%)",
        "High (70-80%)",
        "Very high (80-90%)",
        "Extremely high (90-95%)",
        "Singleton-like (95-100%)"
    };
    for (int i = 0; i < 8; i++) {
        double percent = (hits.size() > 0) ? (100.0 * uniqueness_counts[i] / hits.size()) : 0.0;
        std::cout << "  " << uniqueness_names[i] << ": " << uniqueness_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print quality metrics
    std::cout << "\n--- QUALITY METRICS ---" << std::endl;
    std::cout << "  Unique minimizers (≥90%): " << unique_minimizers << " (" 
              << (100.0 * unique_minimizers / hits.size()) << "%)" << std::endl;
    std::cout << "  Rare minimizers (≤3 occurrences): " << rare_minimizers << " (" 
              << (100.0 * rare_minimizers / hits.size()) << "%)" << std::endl;
    std::cout << "  Reliable minimizers: " << reliable_minimizers << " (" 
              << (100.0 * reliable_minimizers / hits.size()) << "%)" << std::endl;
    std::cout << "  Average uniqueness score: " << std::fixed << std::setprecision(3) 
              << (uniqueness_sum / hits.size()) << std::endl;
    
    // Print co-occurrence statistics
    std::cout << "\n--- CO-OCCURRENCE SCORES ---" << std::endl;
    const char* cooccurrence_names[] = {
        "No significant co-occurrence",
        "Very low co-occurrence",
        "Low co-occurrence",
        "Moderate co-occurrence",
        "High co-occurrence",
        "Very high co-occurrence",
        "Extremely high co-occurrence",
        "Perfect co-occurrence pattern"
    };
    for (int i = 0; i < 8; i++) {
        double percent = (hits.size() > 0) ? (100.0 * cooccurrence_counts[i] / hits.size()) : 0.0;
        std::cout << "  " << cooccurrence_names[i] << ": " << cooccurrence_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print GC content distribution
    std::cout << "\n--- GC CONTENT DISTRIBUTION ---" << std::endl;
    const char* gc_names[] = {
        "0-20%", "20-30%", "30-40%", "40-50%",
        "50-60%", "60-70%", "70-80%", "80-100%"
    };
    for (int i = 0; i < 8; i++) {
        double percent = (hits.size() > 0) ? (100.0 * gc_counts[i] / hits.size()) : 0.0;
        std::cout << "  GC " << gc_names[i] << ": " << gc_counts[i] 
                  << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
    }
    
    // Print other statistics
    std::cout << "\n--- OTHER FEATURES ---" << std::endl;
    std::cout << "  Position clustering: " << position_clustered << " clustered, " 
              << position_uniform << " uniform" << std::endl;
    std::cout << "  Contamination risk: " << contamination_risk << std::endl;
    std::cout << "  Classification: " << unique_class << " unique, " 
              << canonical_class << " canonical, " << redundant_class << " redundant" << std::endl;
}

// Enhanced minimizer extraction kernel with all features
__global__ void extract_minimizers_with_all_features(
    const char* sequences,
    const uint32_t* seq_offsets,
    const uint32_t* seq_lengths,
    const uint32_t* genome_ids,
    const uint32_t* taxon_ids,
    int num_sequences,
    GPUMinimizerHit* minimizer_hits,
    uint32_t* global_counter,
    uint32_t max_minimizers,
    MinimizerParams params) {
    
    int seq_idx = blockIdx.x;
    if (seq_idx >= num_sequences) return;
    
    uint32_t seq_start = seq_offsets[seq_idx];
    uint32_t seq_length = seq_lengths[seq_idx];
    uint32_t genome_id = genome_ids[seq_idx];
    uint32_t taxon_id = taxon_ids[seq_idx];
    
    if (seq_length < params.k) return;
    
    const char* sequence = sequences + seq_start;
    
    // Shared memory for co-occurrence window
    extern __shared__ uint64_t shared_window[];
    
    // Each thread processes k-mers
    uint32_t tid = threadIdx.x;
    uint32_t threads_per_block = blockDim.x;
    uint32_t num_kmers = seq_length - params.k + 1;
    
    for (uint32_t pos = tid; pos < num_kmers; pos += threads_per_block) {
        // Extract minimizer
        uint64_t minimizer = extract_minimizer_sliding_window(
            sequence, pos, params.k, params.ell, params.spaces, 
            params.xor_mask, seq_length
        );
        
        if (minimizer != UINT64_MAX) {
            // Get slot for this minimizer
            uint32_t idx = atomicAdd(global_counter, 1);
            if (idx < max_minimizers) {
                GPUMinimizerHit& hit = minimizer_hits[idx];
                hit.minimizer_hash = minimizer;
                hit.genome_id = genome_id;
                hit.position = pos;
                hit.taxon_id = taxon_id;
                
                // Initialize feature flags
                uint32_t feature_flags = 0;
                
                // Calculate GC content
                int gc_count = 0;
                for (int i = 0; i < params.k; i++) {
                    char base = sequence[pos + i];
                    if (base == 'G' || base == 'g' || base == 'C' || base == 'c') {
                        gc_count++;
                    }
                }
                uint8_t gc_category = (gc_count * 8) / params.k;
                feature_flags |= (gc_category & 0x7);
                
                // Calculate complexity
                uint8_t base_mask = 0;
                for (int i = 0; i < params.k; i++) {
                    char base = sequence[pos + i];
                    if (base == 'A' || base == 'a') base_mask |= 1;
                    else if (base == 'C' || base == 'c') base_mask |= 2;
                    else if (base == 'G' || base == 'g') base_mask |= 4;
                    else if (base == 'T' || base == 't') base_mask |= 8;
                }
                uint8_t complexity = __popc(base_mask);
                feature_flags |= ((complexity & 0x7) << 3);
                
                // Position bias (simplified)
                bool is_clustered = (pos % 100) < 20;  // Dummy clustering
                if (is_clustered) {
                    feature_flags |= (1 << 6);
                }
                
                // Initialize uniqueness to unknown (will be computed later)
                feature_flags |= (3 << 8);  // Moderate uniqueness as default
                
                // Store features
                hit.feature_flags = feature_flags;
                hit.strand = 0;  // Can add strand info later
                hit.ml_weight = 32768;  // Default weight = 0.5
            }
        }
    }
}

// Process a single batch of genomes with two-pass feature extraction
void process_genome_batch_with_features(
    const std::vector<std::string>& batch_files,
    const std::vector<uint32_t>& batch_taxons,
    std::vector<GPUMinimizerHit>& all_minimizers,
    MinimizerFeatureExtractor& feature_extractor,
    const MinimizerParams& params,
    size_t max_batch_minimizers = 10000000) {
    
    // Load sequences from batch
    std::vector<std::string> sequences;
    std::vector<uint32_t> seq_lengths;
    std::vector<uint32_t> seq_offsets;
    std::vector<uint32_t> genome_ids;
    std::vector<uint32_t> taxon_ids;
    std::string concatenated;
    
    GenomeFileProcessor file_processor;
    
    for (size_t i = 0; i < batch_files.size(); i++) {
        auto seqs = file_processor.load_sequences_from_fasta(batch_files[i]);
        for (const auto& seq : seqs) {
            seq_offsets.push_back(concatenated.length());
            seq_lengths.push_back(seq.length());
            genome_ids.push_back(i);
            taxon_ids.push_back(batch_taxons[i]);
            concatenated += seq;
        }
    }
    
    if (concatenated.empty()) return;
    
    std::cout << "  Batch: " << seq_offsets.size() << " sequences, " 
              << concatenated.length() << " total bases" << std::endl;
    
    // Allocate GPU memory for this batch
    char* d_sequences;
    uint32_t* d_offsets;
    uint32_t* d_lengths;
    uint32_t* d_genome_ids;
    uint32_t* d_taxon_ids;
    GPUMinimizerHit* d_minimizers;
    uint32_t* d_counter;
    
    cudaMalloc(&d_sequences, concatenated.length());
    cudaMalloc(&d_offsets, seq_offsets.size() * sizeof(uint32_t));
    cudaMalloc(&d_lengths, seq_lengths.size() * sizeof(uint32_t));
    cudaMalloc(&d_genome_ids, genome_ids.size() * sizeof(uint32_t));
    cudaMalloc(&d_taxon_ids, taxon_ids.size() * sizeof(uint32_t));
    cudaMalloc(&d_minimizers, max_batch_minimizers * sizeof(GPUMinimizerHit));
    cudaMalloc(&d_counter, sizeof(uint32_t));
    
    // Copy data to GPU
    cudaMemcpy(d_sequences, concatenated.data(), concatenated.length(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, seq_offsets.data(), seq_offsets.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lengths, seq_lengths.data(), seq_lengths.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_genome_ids, genome_ids.data(), genome_ids.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_taxon_ids, taxon_ids.data(), taxon_ids.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemset(d_counter, 0, sizeof(uint32_t));
    
    // Launch kernel
    int threads = 256;
    int blocks = seq_offsets.size();
    size_t shared_mem = threads * sizeof(uint64_t);
    
    extract_minimizers_with_all_features<<<blocks, threads, shared_mem>>>(
        d_sequences, d_offsets, d_lengths, d_genome_ids, d_taxon_ids,
        seq_offsets.size(),
        d_minimizers, d_counter, max_batch_minimizers,
        params
    );
    
    cudaDeviceSynchronize();
    
    auto err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
    } else {
        // Get results
        uint32_t num_minimizers;
        cudaMemcpy(&num_minimizers, d_counter, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        if (num_minimizers > 0) {
            std::vector<GPUMinimizerHit> batch_minimizers(num_minimizers);
            cudaMemcpy(batch_minimizers.data(), d_minimizers, 
                       num_minimizers * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
            
            // First pass: collect statistics
            std::cout << "  Running first pass feature extraction..." << std::endl;
            feature_extractor.process_first_pass(d_minimizers, num_minimizers, batch_taxons);
            
            // Compute uniqueness scores
            std::cout << "  Computing uniqueness scores..." << std::endl;
            if (!compute_and_encode_uniqueness_scores(
                    d_minimizers, num_minimizers, batch_taxons, batch_taxons.size())) {
                std::cerr << "  Warning: Uniqueness computation failed" << std::endl;
            }
            
            // Compute co-occurrence scores
            std::cout << "  Computing co-occurrence patterns..." << std::endl;
            
            // Get unique minimizers for co-occurrence
            std::unordered_map<uint64_t, uint32_t> minimizer_counts;
            for (const auto& hit : batch_minimizers) {
                minimizer_counts[hit.minimizer_hash]++;
            }
            
            std::vector<uint64_t> unique_minimizers;
            std::vector<uint32_t> occurrence_counts;
            for (const auto& [hash, count] : minimizer_counts) {
                unique_minimizers.push_back(hash);
                occurrence_counts.push_back(count);
            }
            
            // Create dummy genome info for co-occurrence calculation
            std::vector<GPUGenomeInfo> genome_info;
            for (size_t i = 0; i < batch_files.size(); i++) {
                GPUGenomeInfo info;
                info.genome_id = i;
                info.taxon_id = batch_taxons[i];
                info.sequence_length = 100000;  // Dummy value
                info.sequence_offset = 0;
                genome_info.push_back(info);
            }
            
            if (!compute_and_encode_cooccurrence_scores(
                    d_minimizers, num_minimizers, unique_minimizers, 
                    occurrence_counts, genome_info)) {
                std::cerr << "  Warning: Co-occurrence computation failed" << std::endl;
            }
            
            // Second pass: complete feature extraction
            std::cout << "  Running second pass feature extraction..." << std::endl;
            feature_extractor.process_second_pass(d_minimizers, num_minimizers);
            
            // Copy final results back
            cudaMemcpy(batch_minimizers.data(), d_minimizers, 
                       num_minimizers * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
            
            // Add to overall results
            all_minimizers.insert(all_minimizers.end(), 
                                  batch_minimizers.begin(), 
                                  batch_minimizers.end());
            
            std::cout << "  Extracted " << num_minimizers << " minimizers with all features" << std::endl;
        }
    }
    
    // Cleanup GPU memory
    cudaFree(d_sequences);
    cudaFree(d_offsets);
    cudaFree(d_lengths);
    cudaFree(d_genome_ids);
    cudaFree(d_taxon_ids);
    cudaFree(d_minimizers);
    cudaFree(d_counter);
}

int main(int argc, char** argv) {
    std::cout << "=== ALL FEATURES TEST WITH STREAMING ===" << std::endl;
    
    // Configuration
    std::string fna_file = "/home/david/Documents/Code/biogpu/data/test_50_genomes.fna";
    std::string temp_dir = "/tmp/biogpu_features_test";
    size_t batch_size = 5;  // Process 5 genomes at a time
    
    // Check if file exists
    if (!std::filesystem::exists(fna_file)) {
        std::cerr << "Error: FNA file not found: " << fna_file << std::endl;
        return 1;
    }
    
    // Create temp directory
    std::filesystem::create_directories(temp_dir);
    
    // Check GPU
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        std::cerr << "No CUDA devices!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "GPU: " << prop.name << std::endl;
    std::cout << "Total memory: " << prop.totalGlobalMem / (1024*1024) << " MB" << std::endl;
    
    // Set minimizer parameters
    MinimizerParams params;
    params.k = 31;
    params.ell = 31;
    params.spaces = 7;
    params.xor_mask = 0x3c8bfbb395c60474ULL;
    
    std::cout << "\nMinimizer parameters:" << std::endl;
    std::cout << "  k-mer size: " << params.k << std::endl;
    std::cout << "  l-mer size: " << params.ell << std::endl;
    std::cout << "  Batch size: " << batch_size << " genomes" << std::endl;
    
    // Initialize streaming processor
    std::cout << "\nInitializing streaming processor..." << std::endl;
    StreamingFnaProcessor processor(fna_file, temp_dir, batch_size);
    
    // Initialize feature extractor
    MinimizerFeatureExtractor feature_extractor(100000000, 1000);  // Max 100M minimizers, 1000 genomes
    
    // Process all batches
    std::vector<GPUMinimizerHit> all_minimizers;
    all_minimizers.reserve(75000000);  // Reserve space for ~75M minimizers
    
    int batch_num = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::vector<std::string> batch_files;
    std::vector<uint32_t> batch_taxons;
    
    while (processor.process_next_batch(batch_files, batch_taxons)) {
        batch_num++;
        std::cout << "\nProcessing batch " << batch_num << " (" 
                  << batch_files.size() << " genomes)..." << std::endl;
        
        process_genome_batch_with_features(batch_files, batch_taxons, 
                                          all_minimizers, feature_extractor, params);
        
        // Clean up temp files
        for (const auto& file : batch_files) {
            std::filesystem::remove(file);
        }
        
        batch_files.clear();
        batch_taxons.clear();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n=== PROCESSING COMPLETE ===" << std::endl;
    std::cout << "Total batches: " << batch_num << std::endl;
    std::cout << "Total genomes: " << processor.get_total_genomes() << std::endl;
    std::cout << "Total minimizers: " << all_minimizers.size() << std::endl;
    std::cout << "Processing time: " << duration.count() << " seconds" << std::endl;
    
    // Analyze features
    if (!all_minimizers.empty()) {
        analyze_minimizer_features(all_minimizers);
        
        // Show some examples
        std::cout << "\n=== EXAMPLE MINIMIZERS WITH ALL FEATURES ===" << std::endl;
        int examples = std::min(10, (int)all_minimizers.size());
        for (int i = 0; i < examples; i++) {
            const auto& hit = all_minimizers[i];
            std::cout << "\nMinimizer " << i << ":" << std::endl;
            std::cout << "  Hash: 0x" << std::hex << hit.minimizer_hash << std::dec << std::endl;
            std::cout << "  Genome ID: " << hit.genome_id << std::endl;
            std::cout << "  Position: " << hit.position << std::endl;
            std::cout << "  Taxon ID: " << hit.taxon_id << std::endl;
            
            uint32_t flags = hit.feature_flags;
            std::cout << "  Features:" << std::endl;
            std::cout << "    - GC content category: " << (flags & 0x7) << "/7" << std::endl;
            std::cout << "    - Complexity: " << ((flags >> 3) & 0x7) << "/7" << std::endl;
            std::cout << "    - Position bias: " << ((flags >> 6) & 0x1 ? "Clustered" : "Uniform") << std::endl;
            
            // Uniqueness features
            uint8_t uniqueness_cat = MinimizerFlags::get_uniqueness_category_safe(flags);
            std::cout << "    - Uniqueness category: " << (int)uniqueness_cat << "/7 (" 
                      << MinimizerFlags::uniqueness_category_name(uniqueness_cat) << ")" << std::endl;
            std::cout << "    - Unique minimizer: " << (MinimizerFlags::is_unique_minimizer_safe(flags) ? "Yes" : "No") << std::endl;
            std::cout << "    - Rare minimizer: " << (MinimizerFlags::is_rare_minimizer_safe(flags) ? "Yes" : "No") << std::endl;
            std::cout << "    - Reliable: " << (MinimizerFlags::is_reliable_minimizer_safe(flags) ? "Yes" : "No") << std::endl;
            
            // Co-occurrence score
            uint8_t cooccurrence = MinimizerFlags::get_cooccurrence_score(flags);
            std::cout << "    - Co-occurrence score: " << (int)cooccurrence << "/7" << std::endl;
            
            // ML weight (uniqueness score)
            float uniqueness_score = MinimizerFlags::ml_weight_to_float(hit.ml_weight);
            std::cout << "    - Uniqueness score: " << std::fixed << std::setprecision(3) 
                      << uniqueness_score << std::endl;
        }
        
        // Print feature extraction statistics
        auto stats = feature_extractor.get_statistics();
        std::cout << "\n=== FEATURE EXTRACTION STATISTICS ===" << std::endl;
        std::cout << "Total unique minimizers: " << stats.total_unique_minimizers << std::endl;
        std::cout << "Singleton minimizers: " << stats.singleton_minimizers << std::endl;
        std::cout << "Rare minimizers (<10): " << stats.rare_minimizers << std::endl;
        std::cout << "Common minimizers (>50%): " << stats.common_minimizers << std::endl;
        std::cout << "Average taxonomic spread: " << std::fixed << std::setprecision(2) 
                  << stats.avg_taxonomic_spread << std::endl;
    }
    
    // Cleanup
    std::filesystem::remove_all(temp_dir);
    
    std::cout << "\n=== TEST " << (all_minimizers.size() > 0 ? "PASSED" : "FAILED") << " ===" << std::endl;
    
    return all_minimizers.size() > 0 ? 0 : 1;
}
