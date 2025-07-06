#pragma once

#include <cuda_runtime.h>
#include <cstdint>
#include <memory>

namespace BioGPU {
namespace Unified {

class GPUSequenceBuffer;

struct KmerGeneratorConfig {
    int kmer_size_nucleotide = 31;
    int kmer_size_protein = 15;
    int minimizer_k = 15;
    int minimizer_w = 10;
    bool generate_reverse_complement = true;
    bool generate_minimizers = true;
    size_t max_kmers_per_batch = 100000000;  // 100M
};

class SharedKmerGenerator {
public:
    explicit SharedKmerGenerator(const KmerGeneratorConfig& config = KmerGeneratorConfig());
    ~SharedKmerGenerator();
    
    // Generate k-mers for nucleotide sequences
    void generateNucleotideKmers(
        GPUSequenceBuffer* buffer, 
        cudaStream_t stream);
    
    // Generate k-mers for protein sequences (translated)
    void generateProteinKmers(
        GPUSequenceBuffer* buffer,
        cudaStream_t stream);
    
    // Generate minimizers for taxonomy
    void generateMinimizers(
        GPUSequenceBuffer* buffer,
        cudaStream_t stream);
    
    // Accessors for generated data
    uint64_t* getNucleotideKmers() const { return d_nucleotide_kmers_; }
    uint64_t* getProteinKmers() const { return d_protein_kmers_; }
    uint64_t* getMinimizers() const { return d_minimizers_; }
    
    uint32_t* getNucleotideKmerCounts() const { return d_nucleotide_kmer_counts_; }
    uint32_t* getProteinKmerCounts() const { return d_protein_kmer_counts_; }
    uint32_t* getMinimizerCounts() const { return d_minimizer_counts_; }
    
    // Get bloom filter for initial screening
    void* getBloomFilter() const { return bloom_filter_; }
    
    // Configuration
    const KmerGeneratorConfig& getConfig() const { return config_; }
    
    // Reset counts for new batch
    void resetCounts(cudaStream_t stream);

private:
    KmerGeneratorConfig config_;
    
    // GPU memory for different k-mer types
    uint64_t* d_nucleotide_kmers_;
    uint64_t* d_protein_kmers_;
    uint64_t* d_minimizers_;
    
    uint32_t* d_nucleotide_kmer_counts_;
    uint32_t* d_protein_kmer_counts_;
    uint32_t* d_minimizer_counts_;
    
    // Shared bloom filter
    void* bloom_filter_;
    
    // Temporary buffers
    uint8_t* d_translated_sequences_;
    size_t translated_buffer_size_;
    
    // Helper methods
    void allocateMemory();
    void freeMemory();
};

} // namespace Unified
} // namespace BioGPU