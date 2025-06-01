// bloom_filter.cu
// GPU-accelerated Bloom filter for pre-screening reads before k-mer analysis
// Optimized for fluoroquinolone resistance detection pipeline

#ifndef BLOOM_FILTER_CU
#define BLOOM_FILTER_CU

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace cg = cooperative_groups;

// Debug macros
#define DEBUG_BLOOM 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG_BLOOM) { printf("[BLOOM DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

// Bloom filter constants
#define BLOOM_FILTER_SIZE_BITS 26  // 2^26 bits = 64MB bloom filter
#define BLOOM_FILTER_SIZE (1ULL << BLOOM_FILTER_SIZE_BITS)
#define BLOOM_FILTER_MASK (BLOOM_FILTER_SIZE - 1)
#define NUM_HASH_FUNCTIONS 3
#define BITS_PER_WORD 64
#define BLOOM_WORDS (BLOOM_FILTER_SIZE / BITS_PER_WORD)

// MurmurHash3 constants
#define MURMUR_C1 0xcc9e2d51
#define MURMUR_C2 0x1b873593
#define MURMUR_R1 15
#define MURMUR_R2 13
#define MURMUR_M 5
#define MURMUR_N 0xe6546b64

// Device functions for hashing
__device__ inline uint32_t rotl32(uint32_t x, int8_t r) {
    return (x << r) | (x >> (32 - r));
}

__device__ inline uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    uint32_t h = seed;
    const uint32_t* blocks = (const uint32_t*)key;
    int nblocks = len / 4;
    
    for (int i = 0; i < nblocks; i++) {
        uint32_t k = blocks[i];
        k *= MURMUR_C1;
        k = rotl32(k, MURMUR_R1);
        k *= MURMUR_C2;
        
        h ^= k;
        h = rotl32(h, MURMUR_R2);
        h = h * MURMUR_M + MURMUR_N;
    }
    
    const uint8_t* tail = (const uint8_t*)(key + nblocks * 4);
    uint32_t k1 = 0;
    
    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= MURMUR_C1;
                k1 = rotl32(k1, MURMUR_R1);
                k1 *= MURMUR_C2;
                h ^= k1;
    }
    
    h ^= len;
    h ^= (h >> 16);
    h *= 0x85ebca6b;
    h ^= (h >> 13);
    h *= 0xc2b2ae35;
    h ^= (h >> 16);
    
    return h;
}

// Generate multiple hash values from a k-mer
__device__ inline void generate_bloom_hashes(uint64_t kmer, uint32_t* hashes, int k) {
    // Convert k-mer to bytes for hashing
    uint8_t kmer_bytes[8];
    for (int i = 0; i < 8; i++) {
        kmer_bytes[i] = (kmer >> (i * 8)) & 0xFF;
    }
    
    // Generate NUM_HASH_FUNCTIONS different hash values
    for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
        hashes[i] = murmur3_32(kmer_bytes, sizeof(uint64_t), i * 0x9747b28c) & BLOOM_FILTER_MASK;
    }
}

// Encode base to 2-bit representation
__device__ inline uint8_t encode_base_bloom(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // Invalid base
    }
}

// Encode k-mer from sequence
__device__ inline uint64_t encode_kmer_bloom(const char* seq, int pos, int k) {
    uint64_t kmer = 0;
    for (int i = 0; i < k; i++) {
        uint8_t base = encode_base_bloom(seq[pos + i]);
        if (base == 4) return UINT64_MAX; // Invalid k-mer
        kmer = (kmer << 2) | base;
    }
    return kmer;
}

// Kernel to build Bloom filter from k-mer index
__global__ void build_bloom_filter_kernel(
    const uint64_t* kmers,
    const uint32_t num_kmers,
    uint64_t* bloom_filter,
    const int kmer_length
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    if (tid == 0) {
        DEBUG_PRINT("[BUILD] Building Bloom filter from %d k-mers", num_kmers);
    }
    
    for (uint32_t idx = tid; idx < num_kmers; idx += stride) {
        uint64_t kmer = kmers[idx];
        if (kmer == UINT64_MAX) continue; // Skip invalid k-mers
        
        uint32_t hashes[NUM_HASH_FUNCTIONS];
        generate_bloom_hashes(kmer, hashes, kmer_length);
        
        // Set bits in Bloom filter
        for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
            uint32_t bit_pos = hashes[i];
            uint32_t word_idx = bit_pos / BITS_PER_WORD;
            uint32_t bit_offset = bit_pos % BITS_PER_WORD;
            uint64_t mask = 1ULL << bit_offset;
            
            atomicOr((unsigned long long*)&bloom_filter[word_idx], (unsigned long long)mask);
        }
        
        // Debug first few k-mers
        if (idx < 5) {
            DEBUG_PRINT("[BUILD] K-mer %d: %llu, hashes: %u %u %u", 
                       idx, kmer, hashes[0], hashes[1], hashes[2]);
        }
    }
    
    if (tid == 0) {
        DEBUG_PRINT("[BUILD] Bloom filter construction completed");
    }
}

// Kernel to pre-screen reads using Bloom filter
__global__ void bloom_filter_screen_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const int num_reads,
    const uint64_t* bloom_filter,
    const int kmer_length,
    const int min_kmers_threshold,
    bool* read_passes_filter,
    int* kmers_found_per_read
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    // Shared memory for Bloom filter chunk (improves memory access pattern)
    __shared__ uint64_t shared_bloom[512]; // Cache 512 words = 32KB
    
    if (tid == 0) {
        DEBUG_PRINT("[SCREEN] Screening %d reads with Bloom filter", num_reads);
    }
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        
        if (read_len < kmer_length) {
            read_passes_filter[read_idx] = false;
            kmers_found_per_read[read_idx] = 0;
            continue;
        }
        
        int kmers_in_bloom = 0;
        int total_kmers_checked = 0;
        
        // Check all k-mers in the read
        for (int pos = 0; pos <= read_len - kmer_length; pos++) {
            uint64_t kmer = encode_kmer_bloom(read, pos, kmer_length);
            if (kmer == UINT64_MAX) continue; // Skip invalid k-mers
            
            total_kmers_checked++;
            
            // Generate hash values
            uint32_t hashes[NUM_HASH_FUNCTIONS];
            generate_bloom_hashes(kmer, hashes, kmer_length);
            
            // Check if all hash positions are set in Bloom filter
            bool possibly_in_set = true;
            for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
                uint32_t bit_pos = hashes[i];
                uint32_t word_idx = bit_pos / BITS_PER_WORD;
                uint32_t bit_offset = bit_pos % BITS_PER_WORD;
                uint64_t mask = 1ULL << bit_offset;
                
                // Load from shared memory if cached, otherwise from global
                uint64_t word;
                if (word_idx < 512 && threadIdx.x < 512) {
                    word = shared_bloom[word_idx];
                } else {
                    word = bloom_filter[word_idx];
                }
                
                if ((word & mask) == 0) {
                    possibly_in_set = false;
                    break;
                }
            }
            
            if (possibly_in_set) {
                kmers_in_bloom++;
            }
        }
        
        // Determine if read passes filter
        bool passes = (kmers_in_bloom >= min_kmers_threshold);
        read_passes_filter[read_idx] = passes;
        kmers_found_per_read[read_idx] = kmers_in_bloom;
        
        // Debug info for first few reads
        if (read_idx < 10) {
            DEBUG_PRINT("[SCREEN] Read %d: %d/%d k-mers in Bloom, passes: %s", 
                       read_idx, kmers_in_bloom, total_kmers_checked, 
                       passes ? "YES" : "NO");
        }
    }
    
    if (tid == 0) {
        DEBUG_PRINT("[SCREEN] Bloom filter screening completed");
    }
}

// Host class for Bloom filter management
class BloomFilterGPU {
private:
    uint64_t* d_bloom_filter;
    size_t bloom_size_bytes;
    int kmer_length;
    bool initialized;
    
public:
    BloomFilterGPU(int k = 15) : kmer_length(k), initialized(false) {
        bloom_size_bytes = BLOOM_WORDS * sizeof(uint64_t);
        d_bloom_filter = nullptr;
        
        DEBUG_PRINT("Initializing Bloom filter: size=%zu MB, k=%d", 
                   bloom_size_bytes / (1024*1024), kmer_length);
    }
    
    ~BloomFilterGPU() {
        if (d_bloom_filter) {
            cudaFree(d_bloom_filter);
        }
    }
    
    bool initialize() {
        // Allocate GPU memory for Bloom filter
        cudaError_t err = cudaMalloc(&d_bloom_filter, bloom_size_bytes);
        if (err != cudaSuccess) {
            std::cerr << "Failed to allocate Bloom filter: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        // Initialize to zero
        cudaMemset(d_bloom_filter, 0, bloom_size_bytes);
        
        initialized = true;
        return true;
    }
    
    bool buildFromKmers(const uint64_t* d_kmers, uint32_t num_kmers) {
        if (!initialized) {
            std::cerr << "Bloom filter not initialized" << std::endl;
            return false;
        }
        
        DEBUG_PRINT("Building Bloom filter from %u k-mers", num_kmers);
        
        // Launch kernel to build Bloom filter
        int block_size = 256;
        int grid_size = (num_kmers + block_size - 1) / block_size;
        
        build_bloom_filter_kernel<<<grid_size, block_size>>>(
            d_kmers, num_kmers, d_bloom_filter, kmer_length
        );
        
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Bloom filter build kernel error: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        cudaDeviceSynchronize();
        
        // Calculate and report filter statistics
        calculateStatistics();
        
        return true;
    }
    
    bool screenReads(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        bool* d_read_passes,
        int* d_kmers_found,
        int min_kmers_threshold = 3
    ) {
        if (!initialized) {
            std::cerr << "Bloom filter not initialized" << std::endl;
            return false;
        }
        
        DEBUG_PRINT("Screening %d reads (min k-mers threshold: %d)", num_reads, min_kmers_threshold);
        
        // Launch screening kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        bloom_filter_screen_kernel<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets, num_reads,
            d_bloom_filter, kmer_length, min_kmers_threshold,
            d_read_passes, d_kmers_found
        );
        
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Bloom filter screen kernel error: " << cudaGetErrorString(err) << std::endl;
            return false;
        }
        
        cudaDeviceSynchronize();
        return true;
    }
    
    bool saveToFile(const std::string& filename) {
        if (!initialized) return false;
        
        std::vector<uint64_t> h_bloom(BLOOM_WORDS);
        cudaMemcpy(h_bloom.data(), d_bloom_filter, bloom_size_bytes, cudaMemcpyDeviceToHost);
        
        std::ofstream file(filename, std::ios::binary);
        if (!file.good()) return false;
        
        // Write header
        file.write(reinterpret_cast<const char*>(&kmer_length), sizeof(int));
        int size_bits = BLOOM_FILTER_SIZE_BITS;
        file.write(reinterpret_cast<const char*>(&size_bits), sizeof(int));
        
        // Write Bloom filter data
        file.write(reinterpret_cast<const char*>(h_bloom.data()), bloom_size_bytes);
        
        file.close();
        return true;
    }
    
    bool loadFromFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.good()) return false;
        
        // Read header
        int file_kmer_length, file_size_bits;
        file.read(reinterpret_cast<char*>(&file_kmer_length), sizeof(int));
        file.read(reinterpret_cast<char*>(&file_size_bits), sizeof(int));
        
        if (file_kmer_length != kmer_length || file_size_bits != BLOOM_FILTER_SIZE_BITS) {
            std::cerr << "Bloom filter file parameters don't match" << std::endl;
            return false;
        }
        
        // Read Bloom filter data
        std::vector<uint64_t> h_bloom(BLOOM_WORDS);
        file.read(reinterpret_cast<char*>(h_bloom.data()), bloom_size_bytes);
        file.close();
        
        if (!initialized) initialize();
        
        cudaMemcpy(d_bloom_filter, h_bloom.data(), bloom_size_bytes, cudaMemcpyHostToDevice);
        
        calculateStatistics();
        return true;
    }
    
    void calculateStatistics() {
        // Copy filter to host for analysis
        std::vector<uint64_t> h_bloom(BLOOM_WORDS);
        cudaMemcpy(h_bloom.data(), d_bloom_filter, bloom_size_bytes, cudaMemcpyDeviceToHost);
        
        // Count set bits
        uint64_t bits_set = 0;
        for (uint64_t word : h_bloom) {
            bits_set += __builtin_popcountll(word);
        }
        
        double fill_ratio = (double)bits_set / BLOOM_FILTER_SIZE;
        double false_positive_rate = pow(fill_ratio, NUM_HASH_FUNCTIONS);
        
        std::cout << "Bloom Filter Statistics:" << std::endl;
        std::cout << "  Size: " << bloom_size_bytes / (1024*1024) << " MB" << std::endl;
        std::cout << "  Bits set: " << bits_set << " / " << BLOOM_FILTER_SIZE << std::endl;
        std::cout << "  Fill ratio: " << fill_ratio * 100 << "%" << std::endl;
        std::cout << "  Estimated FP rate: " << false_positive_rate * 100 << "%" << std::endl;
    }
    
    uint64_t* getDevicePointer() { return d_bloom_filter; }
    bool isInitialized() const { return initialized; }
};

// Integration function for existing pipeline
extern "C" {
    void* create_bloom_filter(int kmer_length) {
        BloomFilterGPU* filter = new BloomFilterGPU(kmer_length);
        if (filter->initialize()) {
            return filter;
        }
        delete filter;
        return nullptr;
    }
    
    void destroy_bloom_filter(void* filter) {
        if (filter) {
            delete static_cast<BloomFilterGPU*>(filter);
        }
    }
    
    int build_bloom_filter_from_index(
        void* filter,
        const uint64_t* d_kmers,
        uint32_t num_kmers
    ) {
        if (!filter) return -1;
        BloomFilterGPU* bf = static_cast<BloomFilterGPU*>(filter);
        return bf->buildFromKmers(d_kmers, num_kmers) ? 0 : -1;
    }
    
    int bloom_filter_screen_reads(
        void* filter,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        bool* d_read_passes,
        int* d_kmers_found,
        int min_kmers_threshold
    ) {
        if (!filter) return -1;
        BloomFilterGPU* bf = static_cast<BloomFilterGPU*>(filter);
        return bf->screenReads(d_reads, d_read_lengths, d_read_offsets, 
                              num_reads, d_read_passes, d_kmers_found, 
                              min_kmers_threshold) ? 0 : -1;
    }
    
    int save_bloom_filter(void* filter, const char* filename) {
        if (!filter) return -1;
        BloomFilterGPU* bf = static_cast<BloomFilterGPU*>(filter);
        return bf->saveToFile(filename) ? 0 : -1;
    }
    
    int load_bloom_filter(void* filter, const char* filename) {
        if (!filter) return -1;
        BloomFilterGPU* bf = static_cast<BloomFilterGPU*>(filter);
        return bf->loadFromFile(filename) ? 0 : -1;
    }
}

#endif // BLOOM_FILTER_CU