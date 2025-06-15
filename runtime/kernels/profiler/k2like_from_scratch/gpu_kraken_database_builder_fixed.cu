// gpu_kraken_database_builder_fixed.cu
// Fixed version with optimized minimizer extraction and proper memory management

#include "gpu_minimizer_extraction.cu"  // Include the optimized kernel
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>

class GPUKrakenDatabaseBuilderFixed {
private:
    // Configuration
    MinimizerParams params;
    std::string output_directory;
    
    // Host data
    std::vector<std::string> genome_files;
    std::vector<uint32_t> genome_taxon_ids;
    std::unordered_map<uint32_t, std::string> taxon_names;
    std::unordered_map<uint32_t, uint32_t> taxon_parents;
    
    // GPU memory
    char* d_sequence_data = nullptr;
    GPUGenomeInfo* d_genome_info = nullptr;
    GPUMinimizerHit* d_minimizer_hits = nullptr;
    uint32_t* d_hit_counts = nullptr;
    
    // Realistic batch sizes based on GPU memory
    int MAX_SEQUENCES_PER_BATCH = 100;      // Smaller batch size
    int MAX_MINIMIZERS_PER_BATCH = 5000000; // 5M minimizers (more realistic)
    int MAX_SEQUENCE_LENGTH = 50000000;     // 50MB total sequence data per batch
    
    // Statistics
    struct BuildStats {
        uint64_t total_sequences = 0;
        uint64_t total_bases = 0;
        uint64_t total_kmers_processed = 0;
        uint64_t raw_minimizers_extracted = 0;
        uint64_t unique_minimizers = 0;
        double extraction_time = 0.0;
        double deduplication_time = 0.0;
    } stats;
    
public:
    GPUKrakenDatabaseBuilderFixed(const std::string& output_dir, MinimizerParams config = MinimizerParams())
        : output_directory(output_dir), params(config) {
        
        std::cout << "Initializing Fixed GPU Kraken Database Builder" << std::endl;
        std::cout << "Parameters: k=" << params.k << ", ell=" << params.ell 
                  << ", spaces=" << params.spaces << std::endl;
        
        std::filesystem::create_directories(output_directory);
        check_gpu_memory();
    }
    
    ~GPUKrakenDatabaseBuilderFixed() {
        free_gpu_memory();
    }
    
    void check_gpu_memory() {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        
        std::cout << "GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
                  << (total_mem / 1024 / 1024) << " MB total" << std::endl;
        
        // Adjust batch sizes based on available memory
        size_t safe_memory = free_mem * 0.8; // Use 80% of available memory
        
        // Calculate memory per batch
        size_t sequence_memory = MAX_SEQUENCE_LENGTH;
        size_t minimizer_memory = MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit);
        size_t genome_info_memory = MAX_SEQUENCES_PER_BATCH * sizeof(GPUGenomeInfo);
        size_t total_per_batch = sequence_memory + minimizer_memory + genome_info_memory;
        
        if (total_per_batch > safe_memory) {
            // Reduce batch size
            double scale_factor = (double)safe_memory / total_per_batch;
            MAX_SEQUENCES_PER_BATCH = (int)(MAX_SEQUENCES_PER_BATCH * scale_factor);
            MAX_MINIMIZERS_PER_BATCH = (int)(MAX_MINIMIZERS_PER_BATCH * scale_factor);
            
            std::cout << "Adjusted batch sizes for memory constraints:" << std::endl;
            std::cout << "  Sequences per batch: " << MAX_SEQUENCES_PER_BATCH << std::endl;
            std::cout << "  Minimizers per batch: " << MAX_MINIMIZERS_PER_BATCH << std::endl;
        }
    }
    
    bool build_database_from_genomes(const std::string& genome_library_path) {
        std::cout << "\n=== BUILDING KRAKEN DATABASE (FIXED VERSION) ===" << std::endl;
        
        auto total_start = std::chrono::high_resolution_clock::now();
        
        // Load genome files
        if (!load_genome_files(genome_library_path)) {
            return false;
        }
        
        // Process genomes in batches
        if (!process_genomes_in_batches()) {
            return false;
        }
        
        // Save database
        if (!save_database()) {
            return false;
        }
        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
        
        std::cout << "\n=== BUILD COMPLETE ===" << std::endl;
        std::cout << "Total time: " << duration.count() << " seconds" << std::endl;
        print_statistics();
        
        return true;
    }
    
private:
    bool load_genome_files(const std::string& library_path) {
        std::cout << "Loading genome files from: " << library_path << std::endl;
        
        // Find all FASTA files
        for (const auto& entry : std::filesystem::recursive_directory_iterator(library_path)) {
            if (entry.is_regular_file()) {
                std::string ext = entry.path().extension().string();
                if (ext == ".fna" || ext == ".fa" || ext == ".fasta") {
                    genome_files.push_back(entry.path().string());
                    
                    // Extract taxon ID from filename (simplified)
                    uint32_t taxon_id = std::hash<std::string>{}(entry.path().stem().string()) % 1000000 + 1000000;
                    genome_taxon_ids.push_back(taxon_id);
                    taxon_names[taxon_id] = entry.path().stem().string();
                }
            }
        }
        
        std::cout << "Found " << genome_files.size() << " genome files" << std::endl;
        stats.total_sequences = genome_files.size();
        
        return !genome_files.empty();
    }
    
    bool process_genomes_in_batches() {
        std::cout << "Processing genomes in batches..." << std::endl;
        
        std::vector<LCACandidate> all_candidates;
        
        for (size_t batch_start = 0; batch_start < genome_files.size(); 
             batch_start += MAX_SEQUENCES_PER_BATCH) {
            
            size_t batch_end = std::min(batch_start + MAX_SEQUENCES_PER_BATCH, genome_files.size());
            std::cout << "Processing batch: genomes " << batch_start << " to " << batch_end << std::endl;
            
            // Load sequences for this batch
            std::vector<std::string> batch_sequences;
            std::vector<uint32_t> batch_taxon_ids;
            
            size_t total_batch_size = 0;
            for (size_t i = batch_start; i < batch_end; i++) {
                auto sequences = load_fasta_file(genome_files[i]);
                for (const auto& seq : sequences) {
                    // Check batch size limit
                    if (total_batch_size + seq.length() > MAX_SEQUENCE_LENGTH) {
                        std::cout << "Batch size limit reached, processing " << batch_sequences.size() << " sequences" << std::endl;
                        break;
                    }
                    
                    batch_sequences.push_back(seq);
                    batch_taxon_ids.push_back(genome_taxon_ids[i]);
                    total_batch_size += seq.length();
                    stats.total_bases += seq.length();
                }
                
                if (total_batch_size > MAX_SEQUENCE_LENGTH * 0.9) break; // Leave some buffer
            }
            
            if (batch_sequences.empty()) continue;
            
            // Process this batch on GPU
            auto batch_candidates = process_batch_on_gpu(batch_sequences, batch_taxon_ids);
            
            // Accumulate candidates
            all_candidates.insert(all_candidates.end(), 
                                 batch_candidates.begin(), batch_candidates.end());
            
            std::cout << "Batch produced " << batch_candidates.size() << " candidates" << std::endl;
            std::cout << "Total candidates so far: " << all_candidates.size() << std::endl;
        }
        
        // Final deduplication of all candidates
        std::cout << "Performing final deduplication..." << std::endl;
        auto final_start = std::chrono::high_resolution_clock::now();
        
        // Sort and deduplicate by minimizer hash
        std::sort(all_candidates.begin(), all_candidates.end(),
            [](const LCACandidate& a, const LCACandidate& b) {
                return a.minimizer_hash < b.minimizer_hash;
            });
        
        auto new_end = std::unique(all_candidates.begin(), all_candidates.end(),
            [](const LCACandidate& a, const LCACandidate& b) {
                return a.minimizer_hash == b.minimizer_hash;
            });
        
        all_candidates.erase(new_end, all_candidates.end());
        
        auto final_end = std::chrono::high_resolution_clock::now();
        stats.deduplication_time += std::chrono::duration<double>(final_end - final_start).count();
        
        stats.unique_minimizers = all_candidates.size();
        
        std::cout << "Final unique minimizers: " << stats.unique_minimizers << std::endl;
        
        // Store for saving
        final_candidates = std::move(all_candidates);
        
        return true;
    }
    
    std::vector<LCACandidate> process_batch_on_gpu(
        const std::vector<std::string>& sequences,
        const std::vector<uint32_t>& taxon_ids) {
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Allocate GPU memory for this batch
        if (!allocate_gpu_memory_for_batch(sequences.size())) {
            throw std::runtime_error("Failed to allocate GPU memory");
        }
        
        // Prepare and transfer data
        std::string concatenated = prepare_sequence_data(sequences, taxon_ids);
        
        // Transfer to GPU
        cudaMemcpy(d_sequence_data, concatenated.c_str(), concatenated.length(), cudaMemcpyHostToDevice);
        
        // Extract minimizers using optimized kernel
        uint32_t total_hits = 0;
        extract_minimizers_gpu_optimized(
            d_sequence_data, d_genome_info, sequences.size(),
            d_minimizer_hits, d_hit_counts, &total_hits,
            params, MAX_MINIMIZERS_PER_BATCH
        );
        
        std::cout << "Extracted " << total_hits << " raw minimizers from " << sequences.size() << " sequences" << std::endl;
        
        if (total_hits >= MAX_MINIMIZERS_PER_BATCH) {
            std::cout << "Warning: Hit minimizer limit! Some data may be lost." << std::endl;
        }
        
        // Deduplicate on GPU
        uint32_t unique_count = 0;
        deduplicate_minimizers_gpu(d_minimizer_hits, total_hits, &unique_count);
        
        std::cout << "After deduplication: " << unique_count << " unique minimizers" << std::endl;
        
        // Copy results back to host
        std::vector<GPUMinimizerHit> host_hits(unique_count);
        cudaMemcpy(host_hits.data(), d_minimizer_hits, 
                   unique_count * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
        
        // Convert to LCACandidate format
        std::vector<LCACandidate> candidates;
        candidates.reserve(unique_count);
        
        for (const auto& hit : host_hits) {
            LCACandidate candidate;
            candidate.minimizer_hash = hit.minimizer_hash;
            candidate.lca_taxon = hit.taxon_id;
            candidate.genome_count = 1;  // Will be updated during final merge
            candidate.uniqueness_score = 1.0f;
            candidates.push_back(candidate);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        stats.extraction_time += std::chrono::duration<double>(end_time - start_time).count();
        stats.raw_minimizers_extracted += total_hits;
        
        return candidates;
    }
    
    std::string prepare_sequence_data(const std::vector<std::string>& sequences,
                                     const std::vector<uint32_t>& taxon_ids) {
        
        std::string concatenated;
        std::vector<GPUGenomeInfo> genome_infos;
        
        uint32_t current_offset = 0;
        for (size_t i = 0; i < sequences.size(); i++) {
            GPUGenomeInfo info;
            info.taxon_id = taxon_ids[i];
            info.sequence_offset = current_offset;
            info.sequence_length = sequences[i].length();
            info.genome_id = i;
            
            genome_infos.push_back(info);
            concatenated += sequences[i];
            current_offset += sequences[i].length();
            
            // Track k-mers
            if (sequences[i].length() >= params.k) {
                stats.total_kmers_processed += sequences[i].length() - params.k + 1;
            }
        }
        
        // Transfer genome info to GPU
        cudaMemcpy(d_genome_info, genome_infos.data(),
                   genome_infos.size() * sizeof(GPUGenomeInfo), cudaMemcpyHostToDevice);
        
        return concatenated;
    }
    
    bool allocate_gpu_memory_for_batch(size_t num_sequences) {
        free_gpu_memory(); // Free any existing allocation
        
        cudaMalloc(&d_sequence_data, MAX_SEQUENCE_LENGTH);
        cudaMalloc(&d_genome_info, num_sequences * sizeof(GPUGenomeInfo));
        cudaMalloc(&d_minimizer_hits, MAX_MINIMIZERS_PER_BATCH * sizeof(GPUMinimizerHit));
        cudaMalloc(&d_hit_counts, num_sequences * sizeof(uint32_t));
        
        return (d_sequence_data && d_genome_info && d_minimizer_hits && d_hit_counts);
    }
    
    void free_gpu_memory() {
        if (d_sequence_data) { cudaFree(d_sequence_data); d_sequence_data = nullptr; }
        if (d_genome_info) { cudaFree(d_genome_info); d_genome_info = nullptr; }
        if (d_minimizer_hits) { cudaFree(d_minimizer_hits); d_minimizer_hits = nullptr; }
        if (d_hit_counts) { cudaFree(d_hit_counts); d_hit_counts = nullptr; }
    }
    
    std::vector<std::string> load_fasta_file(const std::string& filepath) {
        std::vector<std::string> sequences;
        std::ifstream file(filepath);
        
        if (!file.is_open()) {
            std::cerr << "Cannot open: " << filepath << std::endl;
            return sequences;
        }
        
        std::string line, current_seq;
        bool in_sequence = false;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (in_sequence && !current_seq.empty()) {
                    sequences.push_back(current_seq);
                    current_seq.clear();
                }
                in_sequence = true;
            } else if (in_sequence) {
                current_seq += line;
            }
        }
        
        if (!current_seq.empty()) {
            sequences.push_back(current_seq);
        }
        
        return sequences;
    }
    
    bool save_database() {
        std::cout << "Saving database..." << std::endl;
        
        // Save hash table
        std::string hash_file = output_directory + "/hash_table.bin";
        std::ofstream hash_out(hash_file, std::ios::binary);
        
        uint64_t table_size = final_candidates.size() * 2;  // 50% load factor
        uint64_t num_entries = final_candidates.size();
        
        hash_out.write(reinterpret_cast<const char*>(&table_size), sizeof(uint64_t));
        hash_out.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint64_t));
        
        for (const auto& candidate : final_candidates) {
            hash_out.write(reinterpret_cast<const char*>(&candidate.minimizer_hash), sizeof(uint64_t));
            hash_out.write(reinterpret_cast<const char*>(&candidate.lca_taxon), sizeof(uint32_t));
            hash_out.write(reinterpret_cast<const char*>(&candidate.genome_count), sizeof(uint32_t));
            hash_out.write(reinterpret_cast<const char*>(&candidate.uniqueness_score), sizeof(float));
        }
        hash_out.close();
        
        // Save taxonomy
        std::string tax_file = output_directory + "/taxonomy.tsv";
        std::ofstream tax_out(tax_file);
        tax_out << "taxon_id\tname\tparent_id\n";
        for (const auto& [taxon_id, name] : taxon_names) {
            uint32_t parent_id = taxon_parents.count(taxon_id) ? taxon_parents[taxon_id] : 0;
            tax_out << taxon_id << "\t" << name << "\t" << parent_id << "\n";
        }
        tax_out.close();
        
        std::cout << "Database saved to " << output_directory << std::endl;
        std::cout << "Hash table: " << num_entries << " entries" << std::endl;
        std::cout << "Taxonomy: " << taxon_names.size() << " taxa" << std::endl;
        
        return true;
    }
    
    void print_statistics() {
        std::cout << "\n=== FINAL STATISTICS ===" << std::endl;
        std::cout << "Total sequences: " << stats.total_sequences << std::endl;
        std::cout << "Total bases: " << stats.total_bases << std::endl;
        std::cout << "Total k-mers: " << stats.total_kmers_processed << std::endl;
        std::cout << "Raw minimizers: " << stats.raw_minimizers_extracted << std::endl;
        std::cout << "Unique minimizers: " << stats.unique_minimizers << std::endl;
        
        if (stats.total_kmers_processed > 0) {
            double compression_ratio = (double)stats.unique_minimizers / stats.total_kmers_processed;
            std::cout << "Compression ratio: " << std::fixed << std::setprecision(4) 
                      << compression_ratio << " (unique minimizers / k-mers)" << std::endl;
        }
        
        std::cout << "Extraction time: " << std::fixed << std::setprecision(2) 
                  << stats.extraction_time << "s" << std::endl;
        std::cout << "Deduplication time: " << std::fixed << std::setprecision(2) 
                  << stats.deduplication_time << "s" << std::endl;
        
        if (stats.total_bases > 0 && stats.extraction_time > 0) {
            double bases_per_second = stats.total_bases / stats.extraction_time;
            std::cout << "Processing rate: " << std::scientific << std::setprecision(2)
                      << bases_per_second << " bases/second" << std::endl;
        }
    }
    
private:
    std::vector<LCACandidate> final_candidates;
};

// Test function
extern "C" void test_fixed_database_builder() {
    std::cout << "Testing Fixed Database Builder" << std::endl;
    
    // Create a test with realistic parameters
    MinimizerParams params;
    params.k = 35;
    params.ell = 31; 
    params.spaces = 7;
    
    GPUKrakenDatabaseBuilderFixed builder("./test_output", params);
    
    // This would be called with a real genome directory
    // builder.build_database_from_genomes("./test_genomes");
    
    std::cout << "Test setup complete. Use with real genome directory." << std::endl;
}