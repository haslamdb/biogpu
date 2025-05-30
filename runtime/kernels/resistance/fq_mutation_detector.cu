#include "fq_mutation_detector.cuh"
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cooperative_groups.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
// #include <hdf5.h>  // Not needed for current implementation
// #include <hdf5_hl.h>
#include <vector>
#include <string>
#include <algorithm>

// Debug macros
#define DEBUG 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { printf("[GPU DEBUG] %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); }

namespace cg = cooperative_groups;

// Device functions
__device__ inline uint8_t encode_base(char base) {
    switch(base) {
        case 'A': case 'a': return BASE_A;
        case 'C': case 'c': return BASE_C;
        case 'G': case 'g': return BASE_G;
        case 'T': case 't': return BASE_T;
        default: return BASE_N;
    }
}

__device__ inline uint64_t encode_kmer(const char* seq, int pos) {
    uint64_t kmer = 0;
    for (int i = 0; i < KMER_LENGTH; i++) {
        uint8_t base = encode_base(seq[pos + i]);
        if (base == BASE_N) return UINT64_MAX; // Invalid k-mer
        kmer = (kmer << 2) | base;
    }
    return kmer;
}

// Stage 1: K-mer filtering kernel
__global__ void kmer_filter_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const KmerEntry* kmer_index,
    const uint64_t* kmer_sorted,
    const uint32_t* kmer_positions,
    const uint32_t num_kmers,
    const int num_reads,
    CandidateMatch* candidates,
    uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    // Debug print from first thread only
    if (tid == 0) {
        DEBUG_PRINT("[KERNEL] kmer_filter_kernel started: num_reads=%d, num_kmers=%d", 
                    num_reads, num_kmers);
        
        // Print first few k-mers from the index
        for (int i = 0; i < min(5, num_kmers); i++) {
            DEBUG_PRINT("[KERNEL] Index k-mer[%d] = %llu", i, kmer_sorted[i]);
        }
    }
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        
        if (read_len < KMER_LENGTH) continue;
        
        // Local candidate tracking
        CandidateMatch local_candidates[MAX_CANDIDATES_PER_READ];
        uint32_t local_count = 0;
        
        // Extract k-mers from read
        for (int pos = 0; pos <= read_len - KMER_LENGTH; pos++) {
            uint64_t kmer = encode_kmer(read, pos);
            if (kmer == UINT64_MAX) continue;
            
            // Binary search in sorted k-mer index
            uint32_t left = 0;
            uint32_t right = num_kmers;
            
            // Debug: print first few k-mers from first read
            if (read_idx == 0 && pos < 3) {
                DEBUG_PRINT("[KERNEL] Read 0, pos %d: k-mer = %llu", pos, kmer);
            }
            
            while (left < right) {
                uint32_t mid = (left + right) / 2;
                if (kmer_sorted[mid] < kmer) {
                    left = mid + 1;
                } else if (kmer_sorted[mid] > kmer) {
                    right = mid;
                } else {
                    // Found match - get all entries with this k-mer
                    uint32_t start = kmer_positions[mid];
                    uint32_t end = (mid + 1 < num_kmers) ? kmer_positions[mid + 1] : num_kmers;
                    
                    if (read_idx < 5) {
                        DEBUG_PRINT("[KERNEL] Read %d: Found k-mer match at position %d", read_idx, pos);
                    }
                    
                    for (uint32_t i = start; i < end && local_count < MAX_CANDIDATES_PER_READ; i++) {
                        const KmerEntry& entry = kmer_index[i];
                        
                        // Check if we already have this candidate
                        bool found = false;
                        for (uint32_t j = 0; j < local_count; j++) {
                            if (local_candidates[j].gene_id == entry.gene_id &&
                                local_candidates[j].species_id == entry.species_id &&
                                local_candidates[j].seq_id == entry.seq_id) {
                                local_candidates[j].kmer_hits++;
                                found = true;
                                break;
                            }
                        }
                        
                        if (!found && local_count < MAX_CANDIDATES_PER_READ) {
                            local_candidates[local_count].gene_id = entry.gene_id;
                            local_candidates[local_count].species_id = entry.species_id;
                            local_candidates[local_count].seq_id = entry.seq_id;
                            local_candidates[local_count].kmer_hits = 1;
                            local_count++;
                        }
                    }
                    break;
                }
            }
        }
        
        // Write candidates to global memory
        uint32_t global_offset = read_idx * max_candidates_per_read;
        for (uint32_t i = 0; i < local_count; i++) {
            candidates[global_offset + i] = local_candidates[i];
        }
        candidate_counts[read_idx] = local_count;
        
        // Debug: print first few reads with candidates
        if (read_idx < 5 && local_count > 0) {
            DEBUG_PRINT("[KERNEL] Read %d: found %d candidates", read_idx, local_count);
        }
    }
    
    // Final debug from first thread
    if (tid == 0) {
        DEBUG_PRINT("[KERNEL] kmer_filter_kernel completed");
    }
}

// Stage 2: Position-weighted alignment kernel
__global__ void position_weighted_alignment_kernel(
    const char* reads,
    const int* read_lengths,
    const int* read_offsets,
    const CandidateMatch* candidates,
    const uint32_t* candidate_counts,
    const uint32_t max_candidates_per_read,
    const char* reference_sequences,
    const uint32_t* ref_lengths,
    const uint32_t* ref_offsets,
    const float* position_weights,
    const uint8_t* mutation_masks,
    const MutationInfo* mutation_info,
    const uint32_t* mutation_counts,
    const AlignmentParams params,
    const int num_reads,
    AlignmentResult* results,
    uint32_t* result_count
) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    
    // Shared memory for dynamic programming matrix
    extern __shared__ float shared_mem[];
    
    for (int read_idx = tid; read_idx < num_reads; read_idx += stride) {
        const char* read = reads + read_offsets[read_idx];
        const int read_len = read_lengths[read_idx];
        uint32_t num_candidates = candidate_counts[read_idx];
        
        if (num_candidates == 0) continue;
        
        // Process each candidate for this read
        for (uint32_t cand_idx = 0; cand_idx < num_candidates; cand_idx++) {
            const CandidateMatch& candidate = candidates[read_idx * max_candidates_per_read + cand_idx];
            
            // Get reference sequence
            uint32_t ref_idx = candidate.seq_id;
            const char* ref_seq = reference_sequences + ref_offsets[ref_idx];
            const uint32_t ref_len = ref_lengths[ref_idx];
            
            // Get position weights and mutation mask for this sequence
            const float* seq_weights = position_weights + ref_offsets[ref_idx];
            const uint8_t* seq_mut_mask = mutation_masks + ref_offsets[ref_idx];
            
            // Sliding window alignment with position weights
            float best_score = -1000.0f;
            int best_pos = -1;
            int best_matches = 0;
            
            // Try different starting positions
            for (int start_pos = max(0, (int)ref_len - read_len - 50); 
                 start_pos < min((int)ref_len - read_len + 50, (int)ref_len); 
                 start_pos++) {
                
                float score = 0.0f;
                int matches = 0;
                
                // Simple scoring (full DP would be more complex)
                for (int i = 0; i < read_len && start_pos + i < ref_len; i++) {
                    uint8_t read_base = encode_base(read[i]);
                    uint8_t ref_base = encode_base(ref_seq[start_pos + i]);
                    
                    if (read_base == BASE_N || ref_base == BASE_N) continue;
                    
                    float pos_weight = seq_weights[start_pos + i];
                    
                    if (read_base == ref_base) {
                        score += params.match_score * pos_weight;
                        matches++;
                        
                        // Bonus for matching at mutation site
                        if (seq_mut_mask[start_pos + i]) {
                            score += params.mutation_match_bonus * pos_weight;
                        }
                    } else {
                        score += params.mismatch_score * pos_weight;
                        
                        // Penalty for mismatch at mutation site
                        if (seq_mut_mask[start_pos + i]) {
                            score -= params.mutation_mismatch_penalty * pos_weight;
                        }
                    }
                }
                
                if (score > best_score) {
                    best_score = score;
                    best_pos = start_pos;
                    best_matches = matches;
                }
            }
            
            // Check if alignment passes threshold
            float identity = (float)best_matches / read_len;
            if (best_score >= params.min_alignment_score && identity >= params.min_identity) {
                // Get result slot
                uint32_t result_idx = atomicAdd(result_count, 1);
                
                AlignmentResult& result = results[result_idx];
                result.read_id = read_idx;
                result.gene_id = candidate.gene_id;
                result.species_id = candidate.species_id;
                result.seq_id = candidate.seq_id;
                result.alignment_score = best_score;
                result.identity = identity;
                result.start_pos = best_pos;
                result.matches = best_matches;
                result.passes_threshold = true;
                
                // Check for mutations
                result.num_mutations_detected = 0;
                const MutationInfo* seq_mutations = mutation_info + ref_idx * MAX_MUTATIONS_PER_GENE;
                uint32_t num_mutations = mutation_counts[ref_idx];
                
                for (uint32_t m = 0; m < num_mutations && result.num_mutations_detected < MAX_MUTATIONS_PER_GENE; m++) {
                    const MutationInfo& mut = seq_mutations[m];
                    
                    // Check if mutation position is covered by alignment
                    if (best_pos <= mut.position && mut.position < best_pos + read_len) {
                        int read_pos = mut.position - best_pos;
                        uint8_t read_base = encode_base(read[read_pos]);
                        uint8_t ref_base = encode_base(ref_seq[mut.position]);
                        
                        // Store mutation detection result
                        result.mutations_detected[result.num_mutations_detected++] = 
                            (read_base == mut.mutant) ? 1 : 0;
                    }
                }
            }
        }
    }
}

// Implementation of FQMutationDetectorCUDA methods
FQMutationDetectorCUDA::FQMutationDetectorCUDA() {
        // Initialize default parameters
        align_params.match_score = 2.0f;
        align_params.mismatch_score = -1.0f;
        align_params.gap_open = -3.0f;
        align_params.gap_extend = -1.0f;
        align_params.mutation_match_bonus = 3.0f;
        align_params.mutation_mismatch_penalty = 2.0f;
        align_params.min_alignment_score = 50.0f;
        align_params.min_identity = 0.8f;
    }

FQMutationDetectorCUDA::~FQMutationDetectorCUDA() {
        // Free device memory
        if (d_kmer_index) cudaFree(d_kmer_index);
        if (d_kmer_sorted) cudaFree(d_kmer_sorted);
        if (d_kmer_positions) cudaFree(d_kmer_positions);
        if (d_reference_sequences) cudaFree(d_reference_sequences);
        if (d_ref_lengths) cudaFree(d_ref_lengths);
        if (d_ref_offsets) cudaFree(d_ref_offsets);
        if (d_position_weights) cudaFree(d_position_weights);
        if (d_mutation_masks) cudaFree(d_mutation_masks);
        if (d_mutation_info) cudaFree(d_mutation_info);
        if (d_mutation_counts) cudaFree(d_mutation_counts);
    }

// Forward declaration - implemented in fq_simple_loader.cu
extern void loadSimpleIndexImpl(const char* index_path, FQMutationDetectorCUDA& detector);

void FQMutationDetectorCUDA::loadIndex(const char* index_path) {
        DEBUG_PRINT("FQMutationDetectorCUDA::loadIndex called with: %s", index_path);
        
        // Check if file exists
        std::ifstream test_file(index_path);
        if (!test_file.good()) {
            std::cerr << "ERROR: Cannot open index file: " << index_path << std::endl;
            DEBUG_PRINT("File does not exist or cannot be opened");
            // Initialize to prevent crashes
            num_kmers = 0;
            num_sequences = 0;
            return;
        }
        test_file.close();
        
        // Check if this is a simple index (by file name or by trying to open)
        std::string path_str(index_path);
        if (path_str.find("simple") != std::string::npos) {
            DEBUG_PRINT("Using simplified index loader");
            // The actual implementation is in fq_simple_loader.cu
            // We'll keep the old implementation below for now
        }
        
        // Open HDF5 file
        hid_t file_id = H5Fopen(index_path, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            std::cerr << "ERROR: Failed to open HDF5 file: " << index_path << std::endl;
            return;
        }
        
        DEBUG_PRINT("HDF5 file opened successfully");
        
        // Read alignment parameters
        hid_t align_group = H5Gopen2(file_id, "alignment", H5P_DEFAULT);
        if (align_group >= 0) {
            // Read attributes
            H5LTget_attribute_float(align_group, ".", "gap_open", &align_params.gap_open);
            H5LTget_attribute_float(align_group, ".", "gap_extend", &align_params.gap_extend);
            H5LTget_attribute_float(align_group, ".", "mutation_match_bonus", &align_params.mutation_match_bonus);
            H5LTget_attribute_float(align_group, ".", "mutation_mismatch_penalty", &align_params.mutation_mismatch_penalty);
            H5LTget_attribute_float(align_group, ".", "min_alignment_score", &align_params.min_alignment_score);
            H5LTget_attribute_float(align_group, ".", "min_identity", &align_params.min_identity);
            
            // Also set match/mismatch scores
            align_params.match_score = 2.0f;
            align_params.mismatch_score = -3.0f;
            
            H5Gclose(align_group);
            DEBUG_PRINT("Loaded alignment parameters");
        }
        
        // Read k-mer index
        hid_t kmer_group = H5Gopen2(file_id, "kmer_index", H5P_DEFAULT);
        if (kmer_group >= 0) {
            // Get k-mer dataset
            hid_t kmer_dataset = H5Dopen2(kmer_group, "kmers", H5P_DEFAULT);
            if (kmer_dataset >= 0) {
                // Get dataset dimensions
                hid_t space = H5Dget_space(kmer_dataset);
                hsize_t dims[2];
                H5Sget_simple_extent_dims(space, dims, NULL);
                num_kmers = dims[0];
                int k = dims[1];  // k-mer length
                
                DEBUG_PRINT("Found %d k-mers of length %d", num_kmers, k);
                
                // Read k-mer strings - they are stored as individual characters
                // Get the datatype from the dataset itself
                hid_t dtype = H5Dget_type(kmer_dataset);
                hid_t native_type = H5Tget_native_type(dtype, H5T_DIR_ASCEND);
                
                char* kmer_chars = new char[num_kmers * k];
                herr_t status = H5Dread(kmer_dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, kmer_chars);
                if (status < 0) {
                    DEBUG_PRINT("Failed to read k-mer dataset: %d", status);
                } else {
                    DEBUG_PRINT("Successfully read k-mer dataset");
                }
                H5Tclose(native_type);
                H5Tclose(dtype);
                
                // Print first k-mer for debugging
                DEBUG_PRINT("First k-mer chars: ");
                for (int i = 0; i < k && i < 15; i++) {
                    DEBUG_PRINT("  [%d] = '%c' (ASCII %d)", i, kmer_chars[i], (int)kmer_chars[i]);
                }
                
                // Convert k-mer strings to encoded values
                std::vector<uint64_t> kmer_values;
                for (int i = 0; i < num_kmers; i++) {
                    uint64_t kmer_val = 0;
                    bool valid = true;
                    for (int j = 0; j < k; j++) {
                        char base = kmer_chars[i * k + j];
                        uint8_t encoded = 0;
                        switch(base) {
                            case 'A': case 'a': encoded = 0; break;
                            case 'C': case 'c': encoded = 1; break;
                            case 'G': case 'g': encoded = 2; break;
                            case 'T': case 't': encoded = 3; break;
                            default: 
                                if (i < 5) {
                                    DEBUG_PRINT("Invalid base '%c' (ASCII %d) at k-mer %d, pos %d", 
                                               base, (int)base, i, j);
                                }
                                valid = false; 
                                break;
                        }
                        if (!valid) break;
                        kmer_val = (kmer_val << 2) | encoded;
                    }
                    if (valid) {
                        kmer_values.push_back(kmer_val);
                    }
                }
                
                // Sort k-mers for binary search
                std::sort(kmer_values.begin(), kmer_values.end());
                num_kmers = kmer_values.size();
                
                DEBUG_PRINT("Encoded and sorted %d k-mers", num_kmers);
                
                // Print first few k-mers for debugging
                for (int i = 0; i < std::min(5, (int)num_kmers); i++) {
                    DEBUG_PRINT("Sorted k-mer[%d] = %llu", i, kmer_values[i]);
                }
                
                // Allocate GPU memory for k-mers
                cudaMalloc(&d_kmer_sorted, num_kmers * sizeof(uint64_t));
                cudaMemcpy(d_kmer_sorted, kmer_values.data(), num_kmers * sizeof(uint64_t), cudaMemcpyHostToDevice);
                
                // Create position array (for binary search)
                uint32_t* h_kmer_positions = new uint32_t[num_kmers];
                for (uint32_t i = 0; i < num_kmers; i++) {
                    h_kmer_positions[i] = i;
                }
                cudaMalloc(&d_kmer_positions, num_kmers * sizeof(uint32_t));
                cudaMemcpy(d_kmer_positions, h_kmer_positions, num_kmers * sizeof(uint32_t), cudaMemcpyHostToDevice);
                
                delete[] kmer_chars;
                delete[] h_kmer_positions;
                
                H5Sclose(space);
                H5Dclose(kmer_dataset);
            }
            H5Gclose(kmer_group);
        }
        
        // Read k-mer lookup from JSON file
        std::string json_path = std::string(index_path);
        size_t pos = json_path.rfind('/');
        if (pos != std::string::npos) {
            json_path = json_path.substr(0, pos + 1) + "kmer_lookup.json";
        }
        
        DEBUG_PRINT("Loading k-mer lookup from: %s", json_path.c_str());
        
        // For now, create a simple k-mer index
        // In a full implementation, we would parse the JSON and create proper KmerEntry structures
        KmerEntry* h_kmer_index = new KmerEntry[num_kmers];
        for (uint32_t i = 0; i < num_kmers; i++) {
            h_kmer_index[i].kmer = i;  // Just use index as placeholder
            h_kmer_index[i].gene_id = 0;
            h_kmer_index[i].species_id = 0;
            h_kmer_index[i].seq_id = 0;
            h_kmer_index[i].position = 0;
        }
        
        cudaMalloc(&d_kmer_index, num_kmers * sizeof(KmerEntry));
        cudaMemcpy(d_kmer_index, h_kmer_index, num_kmers * sizeof(KmerEntry), cudaMemcpyHostToDevice);
        delete[] h_kmer_index;
        
        // Load gene sequences (simplified for now)
        // Count total sequences and get total length
        num_sequences = 0;
        total_ref_length = 0;
        
        // List all gene groups
        hsize_t num_objs;
        H5Gget_num_objs(file_id, &num_objs);
        
        std::vector<std::string> gene_names;
        for (hsize_t i = 0; i < num_objs; i++) {
            char name[256];
            H5Gget_objname_by_idx(file_id, i, name, sizeof(name));
            std::string obj_name(name);
            
            // Skip non-gene groups
            if (obj_name != "alignment" && obj_name != "kmer_index") {
                gene_names.push_back(obj_name);
                
                // Open gene group and count sequences
                hid_t gene_group = H5Gopen2(file_id, name, H5P_DEFAULT);
                if (gene_group >= 0) {
                    hid_t lengths_dataset = H5Dopen2(gene_group, "lengths", H5P_DEFAULT);
                    if (lengths_dataset >= 0) {
                        hid_t space = H5Dget_space(lengths_dataset);
                        hsize_t dims[1];
                        H5Sget_simple_extent_dims(space, dims, NULL);
                        num_sequences += dims[0];
                        
                        // Read lengths to calculate total
                        int32_t* lengths = new int32_t[dims[0]];
                        H5Dread(lengths_dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, lengths);
                        for (hsize_t j = 0; j < dims[0]; j++) {
                            total_ref_length += lengths[j];
                        }
                        delete[] lengths;
                        
                        H5Sclose(space);
                        H5Dclose(lengths_dataset);
                    }
                    H5Gclose(gene_group);
                }
            }
        }
        
        DEBUG_PRINT("Found %d genes with %d total sequences, %d total bases", 
                    (int)gene_names.size(), num_sequences, total_ref_length);
        
        // Allocate GPU memory for reference sequences
        cudaMalloc(&d_reference_sequences, total_ref_length);
        cudaMalloc(&d_ref_lengths, num_sequences * sizeof(uint32_t));
        cudaMalloc(&d_ref_offsets, num_sequences * sizeof(uint32_t));
        cudaMalloc(&d_position_weights, total_ref_length * sizeof(float));
        cudaMalloc(&d_mutation_masks, total_ref_length);
        
        // Initialize with dummy data for now
        cudaMemset(d_reference_sequences, 0, total_ref_length);
        cudaMemset(d_position_weights, 0, total_ref_length * sizeof(float));
        cudaMemset(d_mutation_masks, 0, total_ref_length);
        
        // Allocate mutation info arrays
        cudaMalloc(&d_mutation_info, num_sequences * MAX_MUTATIONS_PER_GENE * sizeof(MutationInfo));
        cudaMalloc(&d_mutation_counts, num_sequences * sizeof(uint32_t));
        cudaMemset(d_mutation_info, 0, num_sequences * MAX_MUTATIONS_PER_GENE * sizeof(MutationInfo));
        cudaMemset(d_mutation_counts, 0, num_sequences * sizeof(uint32_t));
        
        // Close HDF5 file
        H5Fclose(file_id);
        
        std::cerr << "\n=== Index Loading Summary ===" << std::endl;
        std::cerr << "Loaded " << num_kmers << " k-mers" << std::endl;
        std::cerr << "Found " << gene_names.size() << " genes" << std::endl;
        std::cerr << "Total sequences: " << num_sequences << std::endl;
        std::cerr << "Total reference length: " << total_ref_length << " bases" << std::endl;
        std::cerr << "\nNOTE: K-mer to sequence mapping not fully implemented yet" << std::endl;
        std::cerr << "============================\n" << std::endl;
    }

void FQMutationDetectorCUDA::processReads(const char* r1_path, const char* r2_path, const char* output_path) {
        printf("Processing paired-end reads:\n");
        printf("  R1: %s\n", r1_path);
        printf("  R2: %s\n", r2_path);
        printf("  Output: %s\n", output_path);
        
        // TODO: Implement full pipeline
        // 1. Load reads from FASTQ files
        // 2. Transfer to GPU
        // 3. Run stage 1 k-mer filtering
        // 4. Run stage 2 alignment
        // 5. Merge results from R1 and R2
        // 6. Write results to output file
    }

// Main entry point for testing
extern "C" void run_fq_mutation_detection(
    const char* index_path,
    const char* r1_path,
    const char* r2_path,
    const char* output_path
) {
    FQMutationDetectorCUDA detector;
    detector.loadIndex(index_path);
    detector.processReads(r1_path, r2_path, output_path);
}