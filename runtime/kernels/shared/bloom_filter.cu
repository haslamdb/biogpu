// bloom_filter.cu
// GPU-accelerated Bloom filter for pre-screening reads before k-mer analysis
// Optimized for fluoroquinolone resistance detection pipeline
// Enhanced with reverse complement support for R2 reads

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
#define DEBUG_BLOOM 0
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

// Base complement for reverse complement
__device__ inline uint8_t complement_base_bloom(uint8_t base) {
    switch(base) {
        case 0: return 3; // A -> T
        case 1: return 2; // C -> G
        case 2: return 1; // G -> C
        case 3: return 0; // T -> A
        default: return 4; // N -> N
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

// Compute reverse complement of encoded k-mer
__device__ inline uint64_t reverse_complement_kmer(uint64_t kmer, int k) {
    uint64_t rc_kmer = 0;
    for (int i = 0; i < k; i++) {
        uint8_t base = kmer & 3;
        rc_kmer = (rc_kmer << 2) | complement_base_bloom(base);
        kmer >>= 2;
    }
    return rc_kmer;
}

// Encode k-mer and its reverse complement in one pass
__device__ inline void encode_kmer_with_rc_bloom(
    const char* seq, 
    int pos, 
    int k,
    uint64_t* forward_kmer,
    uint64_t* rc_kmer
) {
    *forward_kmer = 0;
    *rc_kmer = 0;
    
    // Build forward k-mer
    for (int i = 0; i < k; i++) {
        uint8_t base = encode_base_bloom(seq[pos + i]);
        if (base == 4) {
            *forward_kmer = UINT64_MAX;
            *rc_kmer = UINT64_MAX;
            return;
        }
        *forward_kmer = (*forward_kmer << 2) | base;
    }
    
    // Build reverse complement
    uint64_t temp = *forward_kmer;
    for (int i = 0; i < k; i++) {
        uint8_t base = temp & 3;
        *rc_kmer = (*rc_kmer << 2) | complement_base_bloom(base);
        temp >>= 2;
    }
}

// Enhanced kernel to build Bloom filter with both forward and RC k-mers
__global__ void build_bloom_filter_kernel_with_rc(
    const uint64_t* kmers,
    const uint32_t num_kmers,
    uint64_t* bloom_filter,
    const int kmer_length,
    const bool include_rc  // New parameter to control RC inclusion
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    if (tid == 0) {
        DEBUG_PRINT("[BUILD] Building Bloom filter from %d k-mers (include_rc: %s)", 
                   num_kmers, include_rc ? "YES" : "NO");
    }
    
    for (uint32_t idx = tid; idx < num_kmers; idx += stride) {
        uint64_t kmer = kmers[idx];
        if (kmer == UINT64_MAX) continue; // Skip invalid k-mers
        
        // Add forward k-mer
        uint32_t hashes[NUM_HASH_FUNCTIONS];
        generate_bloom_hashes(kmer, hashes, kmer_length);
        
        for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
            uint32_t bit_pos = hashes[i];
            uint32_t word_idx = bit_pos / BITS_PER_WORD;
            uint32_t bit_offset = bit_pos % BITS_PER_WORD;
            uint64_t mask = 1ULL << bit_offset;
            
            atomicOr((unsigned long long*)&bloom_filter[word_idx], (unsigned long long)mask);
        }
        
        // Also add reverse complement if requested
        if (include_rc) {
            uint64_t rc_kmer = reverse_complement_kmer(kmer, kmer_length);
            generate_bloom_hashes(rc_kmer, hashes, kmer_length);
            
            for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
                uint32_t bit_pos = hashes[i];
                uint32_t word_idx = bit_pos / BITS_PER_WORD;
                uint32_t bit_offset = bit_pos % BITS_PER_WORD;
                uint64_t mask = 1ULL << bit_offset;
                
                atomicOr((unsigned long long*)&bloom_filter[word_idx], (unsigned long long)mask);
            }
        }
        
        // Debug first few k-mers
        if (idx < 5) {
            DEBUG_PRINT("[BUILD] K-mer %d: %llu (forward), %llu (RC)", 
                       idx, kmer, include_rc ? reverse_complement_kmer(kmer, kmer_length) : 0);
        }
    }
    
    if (tid == 0) {
        DEBUG_PRINT("[BUILD] Bloom filter construction completed");
    }
}

// Enhanced kernel to pre-screen reads with RC support
__global__ void bloom_filter_screen_kernel_with_rc(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const int num_reads,
    const uint64_t* bloom_filter,
    const int kmer_length,
    const int min_kmers_threshold,
    bool* read_passes_filter,
    int* kmers_found_per_read,
    const bool check_rc,        // Check reverse complement for R2 reads
    int* debug_stats = nullptr  // Optional statistics [total, forward_hits, rc_hits, both_hits]
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    // Shared memory for statistics
    __shared__ int shared_stats[4];
    if (threadIdx.x < 4 && debug_stats != nullptr) {
        shared_stats[threadIdx.x] = 0;
    }
    __syncthreads();
    
    if (tid == 0) {
        DEBUG_PRINT("[SCREEN] Screening %d reads with Bloom filter (check_rc: %s)", 
                   num_reads, check_rc ? "YES" : "NO");
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
        int forward_hits = 0;
        int rc_hits = 0;
        int both_hits = 0;
        
        // Check all k-mers in the read
        for (int pos = 0; pos <= read_len - kmer_length; pos++) {
            uint64_t forward_kmer, rc_kmer;
            
            if (check_rc) {
                // For R2 reads, compute both forward and RC
                encode_kmer_with_rc_bloom(read, pos, kmer_length, &forward_kmer, &rc_kmer);
                if (forward_kmer == UINT64_MAX) continue;
            } else {
                // For R1 reads, just compute forward
                forward_kmer = encode_kmer_bloom(read, pos, kmer_length);
                if (forward_kmer == UINT64_MAX) continue;
                rc_kmer = UINT64_MAX;
            }
            
            total_kmers_checked++;
            
            // Check forward k-mer
            uint32_t hashes[NUM_HASH_FUNCTIONS];
            generate_bloom_hashes(forward_kmer, hashes, kmer_length);
            
            bool forward_in_set = true;
            for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
                uint32_t bit_pos = hashes[i];
                uint32_t word_idx = bit_pos / BITS_PER_WORD;
                uint32_t bit_offset = bit_pos % BITS_PER_WORD;
                uint64_t mask = 1ULL << bit_offset;
                
                if ((bloom_filter[word_idx] & mask) == 0) {
                    forward_in_set = false;
                    break;
                }
            }
            
            // Check RC k-mer if needed
            bool rc_in_set = false;
            if (check_rc && !forward_in_set) {
                generate_bloom_hashes(rc_kmer, hashes, kmer_length);
                rc_in_set = true;
                
                for (int i = 0; i < NUM_HASH_FUNCTIONS; i++) {
                    uint32_t bit_pos = hashes[i];
                    uint32_t word_idx = bit_pos / BITS_PER_WORD;
                    uint32_t bit_offset = bit_pos % BITS_PER_WORD;
                    uint64_t mask = 1ULL << bit_offset;
                    
                    if ((bloom_filter[word_idx] & mask) == 0) {
                        rc_in_set = false;
                        break;
                    }
                }
            }
            
            if (forward_in_set || rc_in_set) {
                kmers_in_bloom++;
                if (forward_in_set && rc_in_set) both_hits++;
                else if (forward_in_set) forward_hits++;
                else rc_hits++;
            }
        }
        
        // Determine if read passes filter
        bool passes = (kmers_in_bloom >= min_kmers_threshold);
        read_passes_filter[read_idx] = passes;
        kmers_found_per_read[read_idx] = kmers_in_bloom;
        
        // Update statistics
        if (debug_stats != nullptr && passes) {
            atomicAdd(&shared_stats[0], 1); // total passing reads
            if (forward_hits > 0) atomicAdd(&shared_stats[1], 1);
            if (rc_hits > 0) atomicAdd(&shared_stats[2], 1);
            if (both_hits > 0) atomicAdd(&shared_stats[3], 1);
        }
        
        // Debug info for first few reads
        if (read_idx < 10) {
            DEBUG_PRINT("[SCREEN] Read %d: %d/%d k-mers in Bloom (F:%d, RC:%d, Both:%d), passes: %s", 
                       read_idx, kmers_in_bloom, total_kmers_checked, 
                       forward_hits, rc_hits, both_hits,
                       passes ? "YES" : "NO");
        }
    }
    
    // Write statistics
    __syncthreads();
    if (threadIdx.x < 4 && debug_stats != nullptr) {
        atomicAdd(&debug_stats[threadIdx.x], shared_stats[threadIdx.x]);
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
    bool has_rc_kmers;  // Track if filter includes RC k-mers
    
public:
    BloomFilterGPU(int k = 15) : kmer_length(k), initialized(false), has_rc_kmers(false) {
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
    
    bool buildFromKmers(const uint64_t* d_kmers, uint32_t num_kmers, bool include_rc = true) {
        if (!initialized) {
            std::cerr << "Bloom filter not initialized" << std::endl;
            return false;
        }
        
        has_rc_kmers = include_rc;
        DEBUG_PRINT("Building Bloom filter from %u k-mers (include RC: %s)", 
                   num_kmers, include_rc ? "YES" : "NO");
        
        // Launch kernel to build Bloom filter
        int block_size = 256;
        int grid_size = (num_kmers + block_size - 1) / block_size;
        
        build_bloom_filter_kernel_with_rc<<<grid_size, block_size>>>(
            d_kmers, num_kmers, d_bloom_filter, kmer_length, include_rc
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
        int min_kmers_threshold = 3,
        bool check_rc = false,
        int* d_debug_stats = nullptr
    ) {
        if (!initialized) {
            std::cerr << "Bloom filter not initialized" << std::endl;
            return false;
        }
        
        DEBUG_PRINT("Screening %d reads (min k-mers: %d, check RC: %s, filter has RC: %s)", 
                   num_reads, min_kmers_threshold, 
                   check_rc ? "YES" : "NO",
                   has_rc_kmers ? "YES" : "NO");
        
        // Launch screening kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        bloom_filter_screen_kernel_with_rc<<<grid_size, block_size>>>(
            d_reads, d_read_lengths, d_read_offsets, num_reads,
            d_bloom_filter, kmer_length, min_kmers_threshold,
            d_read_passes, d_kmers_found, check_rc, d_debug_stats
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
        file.write(reinterpret_cast<const char*>(&has_rc_kmers), sizeof(bool));
        
        // Write Bloom filter data
        file.write(reinterpret_cast<const char*>(h_bloom.data()), bloom_size_bytes);
        
        file.close();
        DEBUG_PRINT("Saved Bloom filter (has RC: %s)", has_rc_kmers ? "YES" : "NO");
        return true;
    }
    
    bool loadFromFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.good()) return false;
        
        // Read header
        int file_kmer_length, file_size_bits;
        file.read(reinterpret_cast<char*>(&file_kmer_length), sizeof(int));
        file.read(reinterpret_cast<char*>(&file_size_bits), sizeof(int));
        
        // Try to read RC flag (may not exist in old files)
        bool file_has_rc = false;
        file.read(reinterpret_cast<char*>(&file_has_rc), sizeof(bool));
        if (!file.good()) {
            // Old file format, assume no RC
            file_has_rc = false;
            file.clear();
            file.seekg(sizeof(int) * 2); // Reset to after old header
        }
        
        if (file_kmer_length != kmer_length || file_size_bits != BLOOM_FILTER_SIZE_BITS) {
            std::cerr << "Bloom filter file parameters don't match" << std::endl;
            return false;
        }
        
        has_rc_kmers = file_has_rc;
        
        // Read Bloom filter data
        std::vector<uint64_t> h_bloom(BLOOM_WORDS);
        file.read(reinterpret_cast<char*>(h_bloom.data()), bloom_size_bytes);
        file.close();
        
        if (!initialized) initialize();
        
        cudaMemcpy(d_bloom_filter, h_bloom.data(), bloom_size_bytes, cudaMemcpyHostToDevice);
        
        DEBUG_PRINT("Loaded Bloom filter (has RC: %s)", has_rc_kmers ? "YES" : "NO");
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
        std::cout << "  Contains RC k-mers: " << (has_rc_kmers ? "YES" : "NO") << std::endl;
    }
    
    uint64_t* getDevicePointer() { return d_bloom_filter; }
    bool isInitialized() const { return initialized; }
    bool hasReverseComplement() const { return has_rc_kmers; }
};

// Integration functions for existing pipeline
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
        // Build with RC k-mers by default for better R2 support
        return bf->buildFromKmers(d_kmers, num_kmers, true) ? 0 : -1;
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
        // Default behavior: don't check RC (for backward compatibility)
        return bf->screenReads(d_reads, d_read_lengths, d_read_offsets, 
                              num_reads, d_read_passes, d_kmers_found, 
                              min_kmers_threshold, false, nullptr) ? 0 : -1;
    }
    
    // New function with RC support
    int bloom_filter_screen_reads_with_rc(
        void* filter,
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        bool* d_read_passes,
        int* d_kmers_found,
        int min_kmers_threshold,
        bool check_rc,
        int* d_debug_stats
    ) {
        if (!filter) return -1;
        BloomFilterGPU* bf = static_cast<BloomFilterGPU*>(filter);
        return bf->screenReads(d_reads, d_read_lengths, d_read_offsets, 
                              num_reads, d_read_passes, d_kmers_found, 
                              min_kmers_threshold, check_rc, d_debug_stats) ? 0 : -1;
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
    
    bool bloom_filter_has_rc(void* filter) {
        if (!filter) return false;
        BloomFilterGPU* bf = static_cast<BloomFilterGPU*>(filter);
        return bf->hasReverseComplement();
    }
}

#endif // BLOOM_FILTER_CU