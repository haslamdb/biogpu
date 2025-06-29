// test_direct_minimizer.cu
// Direct test of minimizer extraction from FNA files

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include "gpu_kraken_types.h"
#include "gpu/gpu_database_kernels.h"
#include "memory/gpu_memory_manager.h"
#include "processing/genome_file_processor.h"

// Simple function to read a FASTA file
std::vector<std::string> read_fasta_sequences(const std::string& filename) {
    std::vector<std::string> sequences;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return sequences;
    }
    
    std::string line, current_seq;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                sequences.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }
    
    if (!current_seq.empty()) {
        sequences.push_back(current_seq);
    }
    
    return sequences;
}

int main(int argc, char** argv) {
    std::cout << "=== Direct Minimizer Extraction Test ===" << std::endl;
    
    // Get FNA file
    std::string fna_file = "/home/david/Documents/Code/biogpu/data/type_strain_reference_genomes/GCF_000006765.1_Pseudomonas_aeruginosa_PAO1_g-proteobacteria_genomic.fna";
    if (argc > 1) {
        fna_file = argv[1];
    }
    
    std::cout << "\n1. Reading FNA file: " << fna_file << std::endl;
    
    // Read sequences
    std::vector<std::string> sequences = read_fasta_sequences(fna_file);
    if (sequences.empty()) {
        std::cerr << "No sequences found in file!" << std::endl;
        return 1;
    }
    
    std::cout << "   Found " << sequences.size() << " sequences" << std::endl;
    
    // Concatenate sequences and prepare genome info
    std::string concatenated;
    std::vector<GPUGenomeInfo> genome_info;
    
    size_t offset = 0;
    for (size_t i = 0; i < sequences.size(); i++) {
        GPUGenomeInfo info;
        info.genome_id = i;
        info.taxon_id = 1000 + i;  // Dummy taxon IDs
        info.sequence_offset = offset;
        info.sequence_length = sequences[i].length();
        genome_info.push_back(info);
        
        concatenated += sequences[i];
        offset += sequences[i].length();
        
        std::cout << "   Sequence " << i << ": " << sequences[i].length() << " bp" << std::endl;
    }
    
    std::cout << "   Total length: " << concatenated.length() << " bp" << std::endl;
    
    // Initialize GPU memory manager
    std::cout << "\n2. Initializing GPU..." << std::endl;
    
    MemoryConfig mem_config;
    mem_config.minimizer_capacity = 1000000;  // 1M minimizers
    mem_config.sequence_batch_size = sequences.size();
    
    GPUMemoryManager memory_manager(mem_config);
    
    // Allocate GPU memory
    if (!memory_manager.allocate_sequence_memory(sequences.size(), concatenated.length())) {
        std::cerr << "Failed to allocate sequence memory!" << std::endl;
        return 1;
    }
    
    if (!memory_manager.allocate_minimizer_memory(mem_config.minimizer_capacity)) {
        std::cerr << "Failed to allocate minimizer memory!" << std::endl;
        return 1;
    }
    
    std::cout << "✓ GPU memory allocated" << std::endl;
    
    // Get GPU buffers
    char* d_sequence_data = memory_manager.get_sequence_buffer();
    GPUGenomeInfo* d_genome_info = memory_manager.get_genome_info_buffer();
    GPUMinimizerHit* d_minimizer_hits = memory_manager.get_minimizer_buffer();
    uint32_t* d_hit_counter = memory_manager.get_global_counter();
    
    // Copy data to GPU
    std::cout << "\n3. Copying data to GPU..." << std::endl;
    cudaMemcpy(d_sequence_data, concatenated.c_str(), concatenated.length(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_genome_info, genome_info.data(), genome_info.size() * sizeof(GPUGenomeInfo), cudaMemcpyHostToDevice);
    cudaMemset(d_hit_counter, 0, sizeof(uint32_t));
    
    // Prepare batch data
    GPUBatchData batch_data;
    batch_data.d_sequence_data = d_sequence_data;
    batch_data.d_genome_info = d_genome_info;
    batch_data.d_minimizer_hits = d_minimizer_hits;
    batch_data.d_global_counter = d_hit_counter;
    batch_data.sequence_buffer_size = concatenated.length();
    batch_data.max_genomes = genome_info.size();
    batch_data.max_minimizers = mem_config.minimizer_capacity;
    
    // Set minimizer parameters
    MinimizerParams params;
    params.k = 31;
    params.ell = 31;
    params.spaces = 7;
    params.xor_mask = 0x3c8bfbb395c60474ULL;
    
    std::cout << "✓ Data copied to GPU" << std::endl;
    std::cout << "   K-mer size: " << params.k << std::endl;
    std::cout << "   L-mer size: " << params.ell << std::endl;
    std::cout << "   Spaces: " << params.spaces << std::endl;
    
    // Launch minimizer extraction
    std::cout << "\n4. Extracting minimizers..." << std::endl;
    
    uint32_t total_hits = 0;
    bool success = launch_improved_minimizer_kernel(
        batch_data,
        params,
        0,  // min_clear_hash_value
        0xe37e28c4271b5a2dULL,  // toggle_mask
        &total_hits
    );
    
    if (!success) {
        std::cerr << "Minimizer extraction kernel failed!" << std::endl;
        return 1;
    }
    
    std::cout << "✓ Kernel executed successfully" << std::endl;
    std::cout << "   Total minimizers extracted: " << total_hits << std::endl;
    
    // Copy results back
    if (total_hits > 0) {
        std::cout << "\n5. Retrieving results..." << std::endl;
        
        std::vector<GPUMinimizerHit> minimizer_hits(total_hits);
        cudaMemcpy(minimizer_hits.data(), d_minimizer_hits, 
                   total_hits * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
        
        // Analyze results
        std::cout << "\n6. Minimizer Statistics:" << std::endl;
        
        // Count minimizers per genome
        std::vector<int> minimizers_per_genome(genome_info.size(), 0);
        for (const auto& hit : minimizer_hits) {
            if (hit.genome_id < genome_info.size()) {
                minimizers_per_genome[hit.genome_id]++;
            }
        }
        
        for (size_t i = 0; i < genome_info.size(); i++) {
            std::cout << "   Genome " << i << ": " << minimizers_per_genome[i] << " minimizers";
            if (genome_info[i].sequence_length > 0) {
                double density = (double)minimizers_per_genome[i] / genome_info[i].sequence_length * 1000;
                std::cout << " (density: " << density << " per kb)";
            }
            std::cout << std::endl;
        }
        
        // Show first few minimizers
        std::cout << "\n7. First 10 minimizers:" << std::endl;
        for (size_t i = 0; i < std::min(size_t(10), minimizer_hits.size()); i++) {
            const auto& hit = minimizer_hits[i];
            std::cout << "   [" << i << "] Hash: 0x" << std::hex << hit.minimizer_hash << std::dec
                      << ", Genome: " << hit.genome_id 
                      << ", Position: " << hit.position
                      << ", Taxon: " << hit.taxon_id << std::endl;
        }
        
        std::cout << "\n✓ MINIMIZER EXTRACTION SUCCESSFUL!" << std::endl;
        std::cout << "   Extracted " << total_hits << " minimizers from " 
                  << concatenated.length() << " bases" << std::endl;
        std::cout << "   Average density: " << (double)total_hits / concatenated.length() * 1000 
                  << " minimizers per kb" << std::endl;
    } else {
        std::cout << "\n✗ No minimizers extracted!" << std::endl;
    }
    
    // Cleanup is handled by destructor
    
    return (total_hits > 0) ? 0 : 1;
}