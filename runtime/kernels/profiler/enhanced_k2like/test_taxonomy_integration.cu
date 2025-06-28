// enhanced_k2like/test_taxonomy_integration.cu
// Test script to validate taxonomy flow through minimizer extraction

#include <iostream>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <cuda_runtime.h>

#include "core/gpu_database_builder_core.h"
#include "processing/genome_file_processor.h"
#include "taxonomy/taxonomy_processor.h"
#include "gpu/gpu_database_kernels.h"
#include "memory/gpu_memory_manager.h"

void print_minimizer_taxonomy_stats(const std::vector<GPUMinimizerHit>& hits,
                                   const std::unordered_map<uint32_t, std::string>& taxon_names) {
    // Count minimizers per taxon
    std::unordered_map<uint32_t, uint32_t> taxon_counts;
    std::unordered_map<uint32_t, std::vector<uint64_t>> taxon_minimizers;
    
    for (const auto& hit : hits) {
        taxon_counts[hit.taxon_id]++;
        taxon_minimizers[hit.taxon_id].push_back(hit.minimizer_hash);
    }
    
    std::cout << "\n=== MINIMIZER TAXONOMY STATISTICS ===" << std::endl;
    std::cout << "Total unique taxa: " << taxon_counts.size() << std::endl;
    std::cout << "\nMinimizers per taxon:" << std::endl;
    
    for (const auto& [taxon_id, count] : taxon_counts) {
        std::string name = "Unknown";
        auto it = taxon_names.find(taxon_id);
        if (it != taxon_names.end()) {
            name = it->second;
        }
        
        std::cout << "  Taxon " << taxon_id << " (" << name << "): " 
                  << count << " minimizers" << std::endl;
        
        // Show first few minimizer hashes
        if (count > 0) {
            std::cout << "    First minimizers: ";
            for (size_t i = 0; i < std::min(size_t(3), taxon_minimizers[taxon_id].size()); i++) {
                std::cout << std::hex << taxon_minimizers[taxon_id][i] << std::dec << " ";
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    std::cout << "=== TAXONOMY INTEGRATION TEST ===" << std::endl;
    
    // Check for required files
    std::string data_dir = "../../../../data/";
    std::string nodes_file = data_dir + "nodes.dmp";
    std::string names_file = data_dir + "names.dmp";
    std::string fna_file = data_dir + "test_50_genomes.fna";
    std::string temp_dir = "temp_taxonomy_test";
    
    if (!std::filesystem::exists(nodes_file) || 
        !std::filesystem::exists(names_file) || 
        !std::filesystem::exists(fna_file)) {
        std::cerr << "Required files not found in data/ directory" << std::endl;
        std::cerr << "Expected: nodes.dmp, names.dmp, test_50_genomes.fna" << std::endl;
        return 1;
    }
    
    // Step 1: Load taxonomy
    std::cout << "\n1. Loading NCBI taxonomy..." << std::endl;
    EnhancedNCBITaxonomyProcessor taxonomy_processor;
    
    if (!taxonomy_processor.load_ncbi_taxonomy(nodes_file, names_file)) {
        std::cerr << "Failed to load taxonomy" << std::endl;
        return 1;
    }
    
    std::cout << "✓ Taxonomy loaded: " << taxonomy_processor.get_taxonomy_size() << " taxa" << std::endl;
    
    // Step 2: Process FNA file with streaming processor
    std::cout << "\n2. Processing FNA file..." << std::endl;
    std::filesystem::create_directories(temp_dir);
    
    StreamingFnaProcessor fna_processor(fna_file, temp_dir, 10); // Process 10 genomes at a time
    
    std::vector<std::string> batch_files;
    std::vector<uint32_t> batch_taxons;
    std::vector<std::string> all_genome_files;
    std::vector<uint32_t> all_taxon_ids;
    
    int batch_num = 0;
    while (fna_processor.process_next_batch(batch_files, batch_taxons)) {
        batch_num++;
        std::cout << "  Batch " << batch_num << ": " << batch_files.size() << " genomes" << std::endl;
        
        // Show taxon IDs for first batch
        if (batch_num == 1) {
            std::cout << "  First batch taxon IDs: ";
            for (size_t i = 0; i < std::min(size_t(5), batch_taxons.size()); i++) {
                std::cout << batch_taxons[i] << " ";
            }
            std::cout << std::endl;
        }
        
        all_genome_files.insert(all_genome_files.end(), batch_files.begin(), batch_files.end());
        all_taxon_ids.insert(all_taxon_ids.end(), batch_taxons.begin(), batch_taxons.end());
        
        batch_files.clear();
        batch_taxons.clear();
    }
    
    std::cout << "✓ Processed " << all_genome_files.size() << " genomes" << std::endl;
    std::cout << "  Unique taxon IDs: " << std::set<uint32_t>(all_taxon_ids.begin(), all_taxon_ids.end()).size() << std::endl;
    
    // Step 3: Initialize GPU memory
    std::cout << "\n3. Initializing GPU memory..." << std::endl;
    MemoryConfig mem_config;
    mem_config.minimizer_capacity = 1000000; // 1M minimizers
    mem_config.sequence_batch_size = 50;
    
    GPUMemoryManager memory_manager(mem_config);
    if (!memory_manager.initialize()) {
        std::cerr << "Failed to initialize GPU memory" << std::endl;
        return 1;
    }
    
    // Step 4: Load sequences and prepare for GPU
    std::cout << "\n4. Loading sequences for GPU processing..." << std::endl;
    std::vector<GPUGenomeInfo> genome_info;
    std::string concatenated_sequences;
    size_t sequence_offset = 0;
    
    for (size_t i = 0; i < all_genome_files.size(); i++) {
        std::ifstream file(all_genome_files[i]);
        if (!file.is_open()) continue;
        
        std::string line, sequence;
        while (std::getline(file, line)) {
            if (line[0] != '>') {
                sequence += line;
            }
        }
        file.close();
        
        if (!sequence.empty()) {
            GPUGenomeInfo info;
            info.genome_id = i;
            info.sequence_offset = sequence_offset;
            info.sequence_length = sequence.length();
            info.taxon_id = all_taxon_ids[i];
            info.minimizer_count = 0;
            
            genome_info.push_back(info);
            concatenated_sequences += sequence + '\0';
            sequence_offset += sequence.length() + 1;
        }
    }
    
    std::cout << "✓ Loaded " << genome_info.size() << " sequences" << std::endl;
    std::cout << "  Total sequence length: " << concatenated_sequences.length() << " bp" << std::endl;
    
    // Debug: Show first few genome taxon IDs
    std::cout << "  First 10 genome info (genome_id -> taxon_id): " << std::endl;
    for (size_t i = 0; i < std::min(size_t(10), genome_info.size()); i++) {
        std::cout << "    Genome " << genome_info[i].genome_id 
                  << " -> Taxon " << genome_info[i].taxon_id << std::endl;
    }
    std::cout << std::endl;
    
    // Step 5: Allocate GPU memory and copy data
    std::cout << "\n5. Copying data to GPU..." << std::endl;
    
    if (!memory_manager.allocate_sequence_memory(genome_info.size(), concatenated_sequences.length())) {
        std::cerr << "Failed to allocate sequence memory" << std::endl;
        return 1;
    }
    
    if (!memory_manager.allocate_minimizer_memory(mem_config.minimizer_capacity)) {
        std::cerr << "Failed to allocate minimizer memory" << std::endl;
        return 1;
    }
    
    // Copy to GPU
    char* d_sequences = memory_manager.get_sequence_buffer();
    GPUGenomeInfo* d_genome_info = memory_manager.get_genome_info_buffer();
    GPUMinimizerHit* d_minimizers = memory_manager.get_minimizer_buffer();
    uint32_t* d_hit_counter = memory_manager.get_global_counter();
    
    cudaMemcpy(d_sequences, concatenated_sequences.c_str(), 
               concatenated_sequences.length(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_genome_info, genome_info.data(), 
               genome_info.size() * sizeof(GPUGenomeInfo), cudaMemcpyHostToDevice);
    cudaMemset(d_hit_counter, 0, sizeof(uint32_t));
    
    // Debug: Verify GPU data by copying back first few genomes
    std::vector<GPUGenomeInfo> verify_genomes(5);
    cudaMemcpy(verify_genomes.data(), d_genome_info, 
               std::min(size_t(5), genome_info.size()) * sizeof(GPUGenomeInfo), 
               cudaMemcpyDeviceToHost);
    std::cout << "  GPU verification - first 5 genomes:" << std::endl;
    for (size_t i = 0; i < verify_genomes.size() && i < genome_info.size(); i++) {
        std::cout << "    GPU Genome " << verify_genomes[i].genome_id 
                  << " -> Taxon " << verify_genomes[i].taxon_id << std::endl;
    }
    
    // Step 6: Extract minimizers with taxonomy
    std::cout << "\n6. Extracting minimizers on GPU..." << std::endl;
    
    GPUBatchData batch_data;
    batch_data.d_sequence_data = d_sequences;
    batch_data.d_genome_info = d_genome_info;
    batch_data.d_minimizer_hits = d_minimizers;
    batch_data.d_global_counter = d_hit_counter;
    batch_data.max_genomes = genome_info.size();
    batch_data.max_minimizers = mem_config.minimizer_capacity;
    batch_data.sequence_buffer_size = concatenated_sequences.length();
    
    MinimizerParams params;
    params.k = 31;
    params.ell = 31;
    params.spaces = 7;
    params.xor_mask = 0x3c8bfbb395c60474ULL;
    
    uint32_t total_hits = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (!launch_minimizer_extraction_kernel(batch_data, params, &total_hits)) {
        std::cerr << "Failed to launch minimizer extraction kernel" << std::endl;
        return 1;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "✓ Extracted " << total_hits << " minimizers in " << duration << " seconds" << std::endl;
    
    // Step 7: Copy results back and analyze
    std::cout << "\n7. Analyzing results..." << std::endl;
    
    std::vector<GPUMinimizerHit> minimizer_hits(total_hits);
    cudaMemcpy(minimizer_hits.data(), d_minimizers, 
               total_hits * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
    
    // Get taxon names for display
    std::unordered_map<uint32_t, std::string> taxon_names;
    for (uint32_t taxon_id : all_taxon_ids) {
        taxon_names[taxon_id] = taxonomy_processor.get_scientific_name(taxon_id);
    }
    
    // Display taxonomy statistics
    print_minimizer_taxonomy_stats(minimizer_hits, taxon_names);
    
    // Debug: Show detailed info for first few minimizers
    std::cout << "\nDEBUG: First 10 minimizer details:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(10), minimizer_hits.size()); i++) {
        const auto& hit = minimizer_hits[i];
        std::cout << "  Minimizer " << i << ": "
                  << "hash=" << std::hex << hit.minimizer_hash << std::dec
                  << ", genome_id=" << hit.genome_id
                  << ", taxon_id=" << hit.taxon_id
                  << ", position=" << hit.position << std::endl;
    }
    
    // Step 8: Test phylogenetic calculations
    std::cout << "\n8. Testing phylogenetic calculations..." << std::endl;
    
    // Get unique taxa from results
    std::set<uint32_t> unique_taxa;
    for (const auto& hit : minimizer_hits) {
        unique_taxa.insert(hit.taxon_id);
    }
    
    if (unique_taxa.size() >= 2) {
        std::vector<uint32_t> taxa_list(unique_taxa.begin(), unique_taxa.end());
        
        // Compute LCA
        uint32_t lca = taxonomy_processor.compute_lca_of_species(taxa_list);
        std::cout << "  LCA of all taxa: " << lca << " (" 
                  << taxonomy_processor.get_scientific_name(lca) << ")" << std::endl;
        
        // Calculate phylogenetic spread
        uint8_t spread = taxonomy_processor.calculate_phylogenetic_spread(taxa_list, lca);
        std::cout << "  Phylogenetic spread: " << (int)spread << std::endl;
        
        // Show distances from some taxa to LCA
        std::cout << "  Distances to LCA:" << std::endl;
        for (size_t i = 0; i < std::min(size_t(5), taxa_list.size()); i++) {
            uint8_t dist = taxonomy_processor.calculate_distance_to_lca(taxa_list[i], lca);
            std::cout << "    " << taxon_names[taxa_list[i]] << ": " << (int)dist << " steps" << std::endl;
        }
    }
    
    // Cleanup
    std::filesystem::remove_all(temp_dir);
    
    std::cout << "\n✓ Taxonomy integration test completed successfully!" << std::endl;
    
    return 0;
}
