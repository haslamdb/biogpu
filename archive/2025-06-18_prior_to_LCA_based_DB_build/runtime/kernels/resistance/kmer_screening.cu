// fixed_kmer_screening.cu
// Fixed implementation of k-mer screening for FQ resistance detection
// Enhanced with reverse complement support for R2 reads

#include "fq_mutation_detector.cuh"
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <cooperative_groups.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

// Debug macros
#define DEBUG 0
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { printf("[GPU DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

namespace cg = cooperative_groups;

// Enhanced device functions
__device__ inline uint8_t encode_base(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // Invalid base
    }
}

__device__ inline char decode_base(uint8_t encoded) {
    switch(encoded) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

// Device function for base complement
__device__ inline uint8_t complement_base(uint8_t base) {
    switch(base) {
        case 0: return 3; // A -> T
        case 1: return 2; // C -> G
        case 2: return 1; // G -> C
        case 3: return 0; // T -> A
        default: return 4; // N -> N
    }
}

__device__ inline uint64_t encode_kmer(const char* seq, int pos, int k) {
    uint64_t kmer = 0;
    for (int i = 0; i < k; i++) {
        uint8_t base = encode_base(seq[pos + i]);
        if (base == 4) return UINT64_MAX; // Invalid k-mer
        kmer = (kmer << 2) | base;
    }
    return kmer;
}

// New function: encode k-mer and its reverse complement
__device__ inline void encode_kmer_with_rc(
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
        uint8_t base = encode_base(seq[pos + i]);
        if (base == 4) {
            *forward_kmer = UINT64_MAX;
            *rc_kmer = UINT64_MAX;
            return;
        }
        *forward_kmer = (*forward_kmer << 2) | base;
    }
    
    // Build reverse complement directly
    for (int i = k - 1; i >= 0; i--) {
        uint8_t base = encode_base(seq[pos + i]);
        uint8_t comp = complement_base(base);
        *rc_kmer = (*rc_kmer << 2) | comp;
    }
}

__device__ inline void decode_kmer(uint64_t encoded_kmer, char* output, int k) {
    for (int i = k - 1; i >= 0; i--) {
        output[i] = decode_base(encoded_kmer & 3);
        encoded_kmer >>= 2;
    }
    output[k] = '\0';
}

// Binary search for k-mer in sorted index
__device__ int binary_search_kmer(
    const uint64_t* sorted_kmers,
    uint64_t target_kmer,
    int num_kmers
) {
    int left = 0;
    int right = num_kmers - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        uint64_t mid_kmer = sorted_kmers[mid];
        
        if (mid_kmer == target_kmer) {
            // Find first occurrence
            while (mid > 0 && sorted_kmers[mid - 1] == target_kmer) {
                mid--;
            }
            return mid;
        } else if (mid_kmer < target_kmer) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    return -1; // Not found
}

// Helper device function to process k-mer matches
__device__ void process_kmer_match(
    int kmer_pos,
    const KmerEntry* kmer_index,
    const uint64_t* sorted_kmers,
    const uint32_t* kmer_start_positions,
    const uint32_t num_unique_kmers,
    const uint32_t total_kmer_entries,
    CandidateMatch* local_candidates,
    uint32_t* local_count,
    uint32_t* seen_candidates,
    uint32_t* num_seen,
    const uint32_t max_candidates_per_read
) {
    // Get range of entries for this k-mer
    uint32_t start_entry = kmer_start_positions[kmer_pos];
    uint32_t end_entry = (kmer_pos + 1 < num_unique_kmers) ? 
                        kmer_start_positions[kmer_pos + 1] : total_kmer_entries;
    
    // Process all entries for this k-mer
    for (uint32_t entry_idx = start_entry; entry_idx < end_entry; entry_idx++) {
        if (entry_idx >= total_kmer_entries) break;
        
        const KmerEntry& entry = kmer_index[entry_idx];
        
        // Create unique identifier for this candidate
        uint32_t candidate_id = (entry.species_id << 16) | entry.seq_id;
        
        // Check if we've already seen this candidate
        bool already_seen = false;
        for (uint32_t i = 0; i < *num_seen; i++) {
            if (seen_candidates[i] == candidate_id) {
                // Update existing candidate
                for (uint32_t j = 0; j < *local_count; j++) {
                    if (local_candidates[j].species_id == entry.species_id &&
                        local_candidates[j].seq_id == entry.seq_id) {
                        local_candidates[j].kmer_hits++;
                        already_seen = true;
                        break;
                    }
                }
                break;
            }
        }
        
        // Add new candidate if not seen
        if (!already_seen && *local_count < max_candidates_per_read && 
            *num_seen < max_candidates_per_read) {
            
            local_candidates[*local_count].gene_id = entry.gene_id;
            local_candidates[*local_count].species_id = entry.species_id;
            local_candidates[*local_count].seq_id = entry.seq_id;
            local_candidates[*local_count].kmer_hits = 1;
            
            seen_candidates[*num_seen] = candidate_id;
            (*num_seen)++;
            (*local_count)++;
        }
    }
}

// Enhanced k-mer filtering kernel with reverse complement support
__global__ void enhanced_kmer_filter_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const KmerEntry* kmer_index,
    const uint64_t* sorted_kmers,
    const uint32_t* kmer_start_positions,
    const uint32_t num_unique_kmers,
    const uint32_t total_kmer_entries,
    const int num_reads,
    const int kmer_length,
    CandidateMatch* candidates,
    uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read,
    const bool check_reverse_complement = false,  // New parameter for R2 reads
    int* debug_stats = nullptr                    // Optional statistics collection
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    // Shared memory for statistics
    __shared__ int shared_stats[4];
    if (threadIdx.x < 4 && debug_stats != nullptr) {
        shared_stats[threadIdx.x] = 0;
    }
    __syncthreads();
    
    // Debug from first thread
    if (tid == 0) {
        DEBUG_PRINT("[KERNEL] Enhanced k-mer filter started");
        DEBUG_PRINT("  num_reads: %d", num_reads);
        DEBUG_PRINT("  num_unique_kmers: %d", num_unique_kmers);
        DEBUG_PRINT("  total_kmer_entries: %d", total_kmer_entries);
        DEBUG_PRINT("  kmer_length: %d", kmer_length);
        DEBUG_PRINT("  check_reverse_complement: %s", check_reverse_complement ? "YES" : "NO");
        
        // Debug first few k-mers
        for (int i = 0; i < min(5, num_unique_kmers); i++) {
            char kmer_str[32];
            decode_kmer(sorted_kmers[i], kmer_str, kmer_length);
            DEBUG_PRINT("  sorted_kmer[%d]: %llu (%s)", i, sorted_kmers[i], kmer_str);
        }
    }
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        
        if (read_len < kmer_length) {
            candidate_counts[read_idx] = 0;
            continue;
        }
        
        // Local storage for candidates (avoid too many global memory writes)
        CandidateMatch local_candidates[MAX_CANDIDATES_PER_READ];
        uint32_t local_count = 0;
        
        // Track unique candidates to avoid duplicates
        uint32_t seen_candidates[MAX_CANDIDATES_PER_READ];
        uint32_t num_seen = 0;
        
        int kmers_tested = 0;
        int kmers_found = 0;
        int forward_hits = 0;
        int rc_hits = 0;
        
        // Scan all k-mers in the read
        for (int pos = 0; pos <= read_len - kmer_length; pos++) {
            uint64_t forward_kmer, rc_kmer;
            
            if (check_reverse_complement) {
                // For R2 reads, compute both forward and RC
                encode_kmer_with_rc(read, pos, kmer_length, &forward_kmer, &rc_kmer);
                if (forward_kmer == UINT64_MAX) continue;
            } else {
                // For R1 reads, just compute forward
                forward_kmer = encode_kmer(read, pos, kmer_length);
                if (forward_kmer == UINT64_MAX) continue;
                rc_kmer = UINT64_MAX; // Not needed
            }
            
            kmers_tested++;
            
            // Try forward k-mer first
            int kmer_pos = binary_search_kmer(sorted_kmers, forward_kmer, num_unique_kmers);
            bool found_forward = (kmer_pos >= 0);
            
            // Try reverse complement if checking R2 or if forward not found
            int rc_pos = -1;
            if (check_reverse_complement && !found_forward) {
                rc_pos = binary_search_kmer(sorted_kmers, rc_kmer, num_unique_kmers);
            }
            
            if (found_forward || rc_pos >= 0) {
                kmers_found++;
                
                // Debug for first few reads
                if (read_idx < 3 && kmers_found <= 3) {
                    char kmer_str[32];
                    if (found_forward) {
                        decode_kmer(forward_kmer, kmer_str, kmer_length);
                        DEBUG_PRINT("[KERNEL] Read %d: Found FORWARD k-mer %s at pos %d (index pos %d)", 
                                   read_idx, kmer_str, pos, kmer_pos);
                        forward_hits++;
                    } else {
                        decode_kmer(rc_kmer, kmer_str, kmer_length);
                        DEBUG_PRINT("[KERNEL] Read %d: Found RC k-mer %s at pos %d (index pos %d)", 
                                   read_idx, kmer_str, pos, rc_pos);
                        rc_hits++;
                    }
                }
                
                // Process whichever was found (prefer forward for consistency)
                if (found_forward) {
                    process_kmer_match(kmer_pos, kmer_index, sorted_kmers, kmer_start_positions,
                                     num_unique_kmers, total_kmer_entries,
                                     local_candidates, &local_count, 
                                     seen_candidates, &num_seen, max_candidates_per_read);
                } else if (rc_pos >= 0) {
                    process_kmer_match(rc_pos, kmer_index, sorted_kmers, kmer_start_positions,
                                     num_unique_kmers, total_kmer_entries,
                                     local_candidates, &local_count,
                                     seen_candidates, &num_seen, max_candidates_per_read);
                }
            }
        }
        
        // Update statistics
        if (debug_stats != nullptr && kmers_tested > 0) {
            atomicAdd(&shared_stats[0], 1); // reads processed
            if (forward_hits > 0) atomicAdd(&shared_stats[1], 1); // reads with forward hits
            if (rc_hits > 0) atomicAdd(&shared_stats[2], 1); // reads with RC hits
            if (check_reverse_complement) atomicAdd(&shared_stats[3], 1); // R2 reads processed
        }
        
        // Debug summary for first few reads
        if (read_idx < 5) {
            DEBUG_PRINT("[KERNEL] Read %d summary: %d kmers tested, %d found (%d forward, %d RC), %d candidates", 
                       read_idx, kmers_tested, kmers_found, forward_hits, rc_hits, local_count);
        }
        
        // Write results to global memory
        uint32_t global_offset = read_idx * max_candidates_per_read;
        for (uint32_t i = 0; i < local_count; i++) {
            candidates[global_offset + i] = local_candidates[i];
        }
        candidate_counts[read_idx] = local_count;
    }
    
    // Write statistics
    __syncthreads();
    if (threadIdx.x < 4 && debug_stats != nullptr) {
        atomicAdd(&debug_stats[threadIdx.x], shared_stats[threadIdx.x]);
    }
    
    // Final summary from first thread
    if (tid == 0) {
        DEBUG_PRINT("[KERNEL] Enhanced k-mer filter completed");
    }
}

// Host-side index loader for binary format
class BinaryIndexLoader {
public:
    struct LoadedIndex {
        std::vector<KmerEntry> kmer_entries;
        std::vector<uint64_t> sorted_kmers;
        std::vector<uint32_t> kmer_start_positions;
        uint32_t kmer_length;
        uint32_t num_unique_kmers;
        uint32_t total_entries;
    };
    
    static LoadedIndex loadFromBinary(const std::string& index_dir) {
        LoadedIndex index;
        
        // Load k-mer index
        std::string kmer_path = index_dir + "/kmer_index.bin";
        std::ifstream kmer_file(kmer_path, std::ios::binary);
        
        if (!kmer_file.good()) {
            std::cerr << "ERROR: Cannot open k-mer index: " << kmer_path << std::endl;
            return index;
        }
        
        // Read header
        uint32_t num_entries;
        kmer_file.read(reinterpret_cast<char*>(&num_entries), sizeof(uint32_t));
        kmer_file.read(reinterpret_cast<char*>(&index.kmer_length), sizeof(uint32_t));
        
        std::cout << "Loading " << num_entries << " k-mer entries (k=" << index.kmer_length << ")" << std::endl;
        
        // Read all k-mer entries
        index.kmer_entries.resize(num_entries);
        for (uint32_t i = 0; i < num_entries; i++) {
            KmerEntry entry;
            kmer_file.read(reinterpret_cast<char*>(&entry.kmer), sizeof(uint64_t));
            kmer_file.read(reinterpret_cast<char*>(&entry.gene_id), sizeof(uint32_t));
            kmer_file.read(reinterpret_cast<char*>(&entry.species_id), sizeof(uint32_t));
            kmer_file.read(reinterpret_cast<char*>(&entry.seq_id), sizeof(uint32_t));
            kmer_file.read(reinterpret_cast<char*>(&entry.position), sizeof(uint16_t));
            
            index.kmer_entries[i] = entry;
        }
        kmer_file.close();
        
        std::cout << "Loaded " << index.kmer_entries.size() << " k-mer entries" << std::endl;
        
        // Build sorted k-mer list and position mapping
        std::map<uint64_t, std::vector<uint32_t>> kmer_to_positions;
        
        for (uint32_t i = 0; i < index.kmer_entries.size(); i++) {
            uint64_t kmer = index.kmer_entries[i].kmer;
            kmer_to_positions[kmer].push_back(i);
        }
        
        // Create sorted arrays
        for (const auto& pair : kmer_to_positions) {
            index.sorted_kmers.push_back(pair.first);
            index.kmer_start_positions.push_back(pair.second[0]); // Position of first occurrence
        }
        
        index.num_unique_kmers = index.sorted_kmers.size();
        index.total_entries = index.kmer_entries.size();
        
        std::cout << "Created sorted k-mer list: " << index.num_unique_kmers << " unique k-mers" << std::endl;
        
        // Debug: show first few k-mers
        for (int i = 0; i < std::min(5, (int)index.sorted_kmers.size()); i++) {
            std::cout << "  sorted_kmer[" << i << "] = " << index.sorted_kmers[i] 
                     << " (start_pos: " << index.kmer_start_positions[i] << ")" << std::endl;
        }
        
        return index;
    }
};

// Enhanced FQ mutation detector with proper binary index loading
class EnhancedFQMutationDetector {
private:
    BinaryIndexLoader::LoadedIndex index;
    
    // Device memory
    KmerEntry* d_kmer_entries;
    uint64_t* d_sorted_kmers;
    uint32_t* d_kmer_start_positions;
    
public:
    EnhancedFQMutationDetector() {
        d_kmer_entries = nullptr;
        d_sorted_kmers = nullptr;
        d_kmer_start_positions = nullptr;
    }
    
    ~EnhancedFQMutationDetector() {
        if (d_kmer_entries) cudaFree(d_kmer_entries);
        if (d_sorted_kmers) cudaFree(d_sorted_kmers);
        if (d_kmer_start_positions) cudaFree(d_kmer_start_positions);
    }
    
    bool loadIndex(const std::string& index_dir) {
        std::cout << "Loading index from: " << index_dir << std::endl;
        
        // Load from binary files
        index = BinaryIndexLoader::loadFromBinary(index_dir);
        
        if (index.kmer_entries.empty()) {
            std::cerr << "ERROR: Failed to load index" << std::endl;
            return false;
        }
        
        // Allocate GPU memory
        size_t entries_size = index.total_entries * sizeof(KmerEntry);
        size_t sorted_size = index.num_unique_kmers * sizeof(uint64_t);
        size_t positions_size = index.num_unique_kmers * sizeof(uint32_t);
        
        std::cout << "Allocating GPU memory:" << std::endl;
        std::cout << "  K-mer entries: " << entries_size << " bytes" << std::endl;
        std::cout << "  Sorted k-mers: " << sorted_size << " bytes" << std::endl;
        std::cout << "  Start positions: " << positions_size << " bytes" << std::endl;
        
        cudaMalloc(&d_kmer_entries, entries_size);
        cudaMalloc(&d_sorted_kmers, sorted_size);
        cudaMalloc(&d_kmer_start_positions, positions_size);
        
        // Copy to GPU
        cudaMemcpy(d_kmer_entries, index.kmer_entries.data(), entries_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_sorted_kmers, index.sorted_kmers.data(), sorted_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_kmer_start_positions, index.kmer_start_positions.data(), positions_size, cudaMemcpyHostToDevice);
        
        // Check for CUDA errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA error: " << cudaGetErrorString(error) << std::endl;
            return false;
        }
        
        std::cout << "Index loaded successfully to GPU" << std::endl;
        return true;
    }
    
    void screenReads(
        const char* d_reads,
        const int* d_read_lengths,
        const int* d_read_offsets,
        int num_reads,
        CandidateMatch* d_candidates,
        uint32_t* d_candidate_counts,
        bool check_reverse_complement = false,
        int* d_debug_stats = nullptr
    ) {
        std::cout << "Screening " << num_reads << " reads";
        if (check_reverse_complement) {
            std::cout << " (with reverse complement)";
        }
        std::cout << "..." << std::endl;
        
        // Launch enhanced k-mer filter kernel
        int block_size = 256;
        int grid_size = (num_reads + block_size - 1) / block_size;
        
        std::cout << "Launching kernel: grid_size=" << grid_size 
                 << ", block_size=" << block_size << std::endl;
        
        enhanced_kmer_filter_kernel<<<grid_size, block_size>>>(
            d_reads,
            d_read_lengths,
            d_read_offsets,
            d_kmer_entries,
            d_sorted_kmers,
            d_kmer_start_positions,
            index.num_unique_kmers,
            index.total_entries,
            num_reads,
            index.kmer_length,
            d_candidates,
            d_candidate_counts,
            MAX_CANDIDATES_PER_READ,
            check_reverse_complement,
            d_debug_stats
        );
        
        // Check for kernel launch errors
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "Kernel launch error: " << cudaGetErrorString(error) << std::endl;
            return;
        }
        
        // Wait for completion
        cudaDeviceSynchronize();
        
        error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "Kernel execution error: " << cudaGetErrorString(error) << std::endl;
            return;
        }
        
        std::cout << "K-mer screening completed" << std::endl;
    }
    
    void debugResults(
        uint32_t* d_candidate_counts,
        CandidateMatch* d_candidates,
        int num_reads,
        int* d_debug_stats = nullptr
    ) {
        std::cout << "Debugging screening results..." << std::endl;
        
        // Copy results back to host
        std::vector<uint32_t> h_candidate_counts(num_reads);
        cudaMemcpy(h_candidate_counts.data(), d_candidate_counts, 
                   num_reads * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        // Get debug statistics if available
        if (d_debug_stats != nullptr) {
            int h_stats[4];
            cudaMemcpy(h_stats, d_debug_stats, 4 * sizeof(int), cudaMemcpyDeviceToHost);
            
            std::cout << "\nK-mer screening statistics:" << std::endl;
            std::cout << "  Total reads processed: " << h_stats[0] << std::endl;
            std::cout << "  Reads with forward hits: " << h_stats[1] << std::endl;
            std::cout << "  Reads with RC hits: " << h_stats[2] << std::endl;
            std::cout << "  R2 reads processed: " << h_stats[3] << std::endl;
        }
        
        // Analyze results
        int total_candidates = 0;
        int reads_with_candidates = 0;
        int max_candidates = 0;
        
        for (int i = 0; i < num_reads; i++) {
            uint32_t count = h_candidate_counts[i];
            total_candidates += count;
            if (count > 0) reads_with_candidates++;
            if (count > max_candidates) max_candidates = count;
        }
        
        std::cout << "\nScreening results summary:" << std::endl;
        std::cout << "  Total reads: " << num_reads << std::endl;
        std::cout << "  Reads with candidates: " << reads_with_candidates << std::endl;
        std::cout << "  Total candidates: " << total_candidates << std::endl;
        std::cout << "  Max candidates per read: " << max_candidates << std::endl;
        std::cout << "  Hit rate: " << (double)reads_with_candidates / num_reads << std::endl;
        
        // Show detailed results for first few reads
        std::cout << "\nDetailed results for first 10 reads:" << std::endl;
        for (int i = 0; i < std::min(10, num_reads); i++) {
            std::cout << "  Read " << i << ": " << h_candidate_counts[i] << " candidates" << std::endl;
        }
        
        // Show some candidate details
        if (total_candidates > 0) {
            std::cout << "\nSample candidate details:" << std::endl;
            
            std::vector<CandidateMatch> h_candidates(num_reads * MAX_CANDIDATES_PER_READ);
            cudaMemcpy(h_candidates.data(), d_candidates, 
                       num_reads * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch), 
                       cudaMemcpyDeviceToHost);
            
            int shown = 0;
            for (int read_idx = 0; read_idx < num_reads && shown < 10; read_idx++) {
                if (h_candidate_counts[read_idx] > 0) {
                    uint32_t offset = read_idx * MAX_CANDIDATES_PER_READ;
                    for (uint32_t i = 0; i < h_candidate_counts[read_idx] && shown < 10; i++) {
                        const CandidateMatch& candidate = h_candidates[offset + i];
                        std::cout << "    Read " << read_idx << ", Candidate " << i 
                                 << ": gene=" << candidate.gene_id 
                                 << ", species=" << candidate.species_id
                                 << ", seq=" << candidate.seq_id
                                 << ", hits=" << candidate.kmer_hits << std::endl;
                        shown++;
                    }
                }
            }
        }
    }
};

// Test function for validating k-mer screening
extern "C" {
    int test_kmer_screening(
        const char* index_dir,
        const char* reads_data,
        const int* read_lengths,
        const int* read_offsets,
        int num_reads,
        bool test_reverse_complement = false
    ) {
        std::cout << "=== Testing K-mer Screening ===" << std::endl;
        
        // Create detector and load index
        EnhancedFQMutationDetector detector;
        if (!detector.loadIndex(index_dir)) {
            std::cerr << "Failed to load index" << std::endl;
            return -1;
        }
        
        // Allocate GPU memory for reads
        size_t total_bases = 0;
        for (int i = 0; i < num_reads; i++) {
            total_bases += read_lengths[i];
        }
        
        char* d_reads;
        int* d_read_lengths;
        int* d_read_offsets;
        int* d_debug_stats;
        
        cudaMalloc(&d_reads, total_bases);
        cudaMalloc(&d_read_lengths, num_reads * sizeof(int));
        cudaMalloc(&d_read_offsets, num_reads * sizeof(int));
        cudaMalloc(&d_debug_stats, 4 * sizeof(int));
        cudaMemset(d_debug_stats, 0, 4 * sizeof(int));
        
        cudaMemcpy(d_reads, reads_data, total_bases, cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_lengths, read_lengths, num_reads * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_read_offsets, read_offsets, num_reads * sizeof(int), cudaMemcpyHostToDevice);
        
        // Allocate memory for results
        CandidateMatch* d_candidates;
        uint32_t* d_candidate_counts;
        
        cudaMalloc(&d_candidates, num_reads * MAX_CANDIDATES_PER_READ * sizeof(CandidateMatch));
        cudaMalloc(&d_candidate_counts, num_reads * sizeof(uint32_t));
        cudaMemset(d_candidate_counts, 0, num_reads * sizeof(uint32_t));
        
        // Run screening
        detector.screenReads(d_reads, d_read_lengths, d_read_offsets, num_reads, 
                           d_candidates, d_candidate_counts, test_reverse_complement, d_debug_stats);
        
        // Debug results
        detector.debugResults(d_candidate_counts, d_candidates, num_reads, d_debug_stats);
        
        // Cleanup
        cudaFree(d_reads);
        cudaFree(d_read_lengths);
        cudaFree(d_read_offsets);
        cudaFree(d_candidates);
        cudaFree(d_candidate_counts);
        cudaFree(d_debug_stats);
        
        std::cout << "=== K-mer Screening Test Complete ===" << std::endl;
        return 0;
    }
}