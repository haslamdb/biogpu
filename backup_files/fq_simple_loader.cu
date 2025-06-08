#include "fq_mutation_detector.cuh"
// HDF5 no longer needed - using binary format
#include <vector>
#include <algorithm>
#include <iostream>

// Debug macros
#define DEBUG 1
#define DEBUG_PRINT(fmt, ...) if(DEBUG) { printf("[LOADER DEBUG] " fmt "\n", ##__VA_ARGS__); }

// Simple structures for the unified index
struct ResistanceMarker {
    uint32_t marker_id;
    uint32_t position;
    char gene_name[32];
    char species[64];
    char mutation[16];
    char ref_aa;
    char mut_aa;
};

struct SimpleFQIndex {
    // Markers
    ResistanceMarker* markers;
    uint32_t num_markers;
    
    // Sequences
    char* sequence_data;
    uint32_t* marker_ids;
    uint32_t* seq_starts;
    uint32_t* seq_lengths;
    uint32_t num_sequences;
    
    // K-mers
    uint64_t* kmer_values;
    uint32_t* marker_id_data;
    uint32_t* marker_id_offsets;
    uint32_t num_kmers;
    uint32_t kmer_length;
};

// Load simplified index from HDF5
// COMMENTED OUT - This function loads HDF5 format, we now use binary format
/*
SimpleFQIndex loadSimpleIndex(const char* filename) {
    SimpleFQIndex index = {};
    
    DEBUG_PRINT("Loading simplified index from: %s", filename);
    
    // Open HDF5 file
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "ERROR: Failed to open HDF5 file: " << filename << std::endl;
        return index;
    }
    
    // Read metadata
    hid_t meta_group = H5Gopen2(file_id, "metadata", H5P_DEFAULT);
    if (meta_group >= 0) {
        H5LTget_attribute_uint(meta_group, ".", "num_markers", &index.num_markers);
        H5LTget_attribute_uint(meta_group, ".", "num_sequences", &index.num_sequences);
        H5LTget_attribute_uint(meta_group, ".", "num_kmers", &index.num_kmers);
        H5LTget_attribute_uint(meta_group, ".", "kmer_length", &index.kmer_length);
        H5Gclose(meta_group);
        
        DEBUG_PRINT("Metadata: %d markers, %d sequences, %d kmers (length %d)",
                    index.num_markers, index.num_sequences, index.num_kmers, index.kmer_length);
    }
    
    // Read markers
    if (index.num_markers > 0) {
        index.markers = new ResistanceMarker[index.num_markers];
        memset(index.markers, 0, index.num_markers * sizeof(ResistanceMarker));
        
        hid_t markers_group = H5Gopen2(file_id, "markers", H5P_DEFAULT);
        if (markers_group >= 0) {
            // Read marker IDs
            uint64_t* marker_ids = new uint64_t[index.num_markers];
            H5LTread_dataset(file_id, "markers/marker_id", H5T_NATIVE_INT64, marker_ids);
            
            // Read positions
            uint64_t* positions = new uint64_t[index.num_markers];
            H5LTread_dataset(file_id, "markers/position", H5T_NATIVE_INT64, positions);
            
            // Copy to structures
            for (uint32_t i = 0; i < index.num_markers; i++) {
                index.markers[i].marker_id = marker_ids[i];
                index.markers[i].position = positions[i];
            }
            
            delete[] marker_ids;
            delete[] positions;
            
            // Read string data (simplified - in real code, handle variable-length strings properly)
            // For now, just set test values
            for (uint32_t i = 0; i < index.num_markers; i++) {
                sprintf(index.markers[i].gene_name, "gene_%d", i);
                sprintf(index.markers[i].species, "species_%d", i);
                sprintf(index.markers[i].mutation, "mut_%d", i);
                index.markers[i].ref_aa = 'A';
                index.markers[i].mut_aa = 'T';
            }
            
            H5Gclose(markers_group);
        }
    }
    
    // Read k-mers
    if (index.num_kmers > 0) {
        // Read k-mer values
        index.kmer_values = new uint64_t[index.num_kmers];
        H5LTread_dataset(file_id, "kmers/kmer_encoded", H5T_NATIVE_UINT64, index.kmer_values);
        
        // Read marker ID offsets
        index.marker_id_offsets = new uint32_t[index.num_kmers + 1];
        uint64_t* temp_offsets = new uint64_t[index.num_kmers + 1];
        H5LTread_dataset(file_id, "kmers/marker_id_offsets", H5T_NATIVE_INT64, temp_offsets);
        for (uint32_t i = 0; i <= index.num_kmers; i++) {
            index.marker_id_offsets[i] = temp_offsets[i];
        }
        delete[] temp_offsets;
        
        // Read marker IDs
        uint32_t total_marker_ids = index.marker_id_offsets[index.num_kmers];
        index.marker_id_data = new uint32_t[total_marker_ids];
        uint64_t* temp_ids = new uint64_t[total_marker_ids];
        H5LTread_dataset(file_id, "kmers/marker_ids", H5T_NATIVE_INT64, temp_ids);
        for (uint32_t i = 0; i < total_marker_ids; i++) {
            index.marker_id_data[i] = temp_ids[i];
        }
        delete[] temp_ids;
        
        DEBUG_PRINT("Loaded %d k-mers with %d total marker associations",
                    index.num_kmers, total_marker_ids);
        
        // Debug: print first few k-mers
        for (uint32_t i = 0; i < std::min(5u, index.num_kmers); i++) {
            uint32_t start = index.marker_id_offsets[i];
            uint32_t end = index.marker_id_offsets[i + 1];
            DEBUG_PRINT("K-mer[%d] = %llu, maps to %d markers", 
                       i, index.kmer_values[i], end - start);
        }
    }
    
    // Read sequences
    if (index.num_sequences > 0) {
        // Read sequence metadata
        index.marker_ids = new uint32_t[index.num_sequences];
        index.seq_starts = new uint32_t[index.num_sequences];
        index.seq_lengths = new uint32_t[index.num_sequences];
        
        uint64_t* temp_data = new uint64_t[index.num_sequences];
        
        H5LTread_dataset(file_id, "sequences/marker_id", H5T_NATIVE_INT64, temp_data);
        for (uint32_t i = 0; i < index.num_sequences; i++) {
            index.marker_ids[i] = temp_data[i];
        }
        
        H5LTread_dataset(file_id, "sequences/seq_start", H5T_NATIVE_INT64, temp_data);
        for (uint32_t i = 0; i < index.num_sequences; i++) {
            index.seq_starts[i] = temp_data[i];
        }
        
        H5LTread_dataset(file_id, "sequences/seq_length", H5T_NATIVE_INT64, temp_data);
        for (uint32_t i = 0; i < index.num_sequences; i++) {
            index.seq_lengths[i] = temp_data[i];
        }
        
        delete[] temp_data;
        
        // Read concatenated sequence data
        // Get total length
        uint32_t total_seq_length = 0;
        for (uint32_t i = 0; i < index.num_sequences; i++) {
            total_seq_length = std::max(total_seq_length, 
                                       index.seq_starts[i] + index.seq_lengths[i]);
        }
        
        // For now, create dummy sequence data
        index.sequence_data = new char[total_seq_length + 1];
        memset(index.sequence_data, 'A', total_seq_length);
        index.sequence_data[total_seq_length] = '\0';
        
        DEBUG_PRINT("Loaded %d sequences, total length %d", 
                    index.num_sequences, total_seq_length);
    }
    
    H5Fclose(file_id);
    
    DEBUG_PRINT("Index loading complete!");
    return index;
}

// Free index memory
void freeSimpleIndex(SimpleFQIndex& index) {
    delete[] index.markers;
    delete[] index.sequence_data;
    delete[] index.marker_ids;
    delete[] index.seq_starts;
    delete[] index.seq_lengths;
    delete[] index.kmer_values;
    delete[] index.marker_id_data;
    delete[] index.marker_id_offsets;
    
    // Reset to empty
    index = {};
}

// Convert simple index to GPU-compatible format
void convertToGPUIndex(const SimpleFQIndex& simple_index, FQMutationDetectorCUDA& detector) {
    DEBUG_PRINT("Converting simple index to GPU format...");
    
    // Set basic parameters
    detector.num_kmers = simple_index.num_kmers;
    detector.num_sequences = simple_index.num_sequences;
    detector.total_ref_length = 0;
    
    // Calculate total reference length
    for (uint32_t i = 0; i < simple_index.num_sequences; i++) {
        detector.total_ref_length += simple_index.seq_lengths[i];
    }
    
    // Allocate and copy k-mers
    if (simple_index.num_kmers > 0) {
        cudaMalloc(&detector.d_kmer_sorted, simple_index.num_kmers * sizeof(uint64_t));
        cudaMemcpy(detector.d_kmer_sorted, simple_index.kmer_values, 
                   simple_index.num_kmers * sizeof(uint64_t), cudaMemcpyHostToDevice);
        
        // Create k-mer index entries
        KmerEntry* h_kmer_index = new KmerEntry[simple_index.num_kmers];
        for (uint32_t i = 0; i < simple_index.num_kmers; i++) {
            h_kmer_index[i].kmer = simple_index.kmer_values[i];
            
            // Get first marker for this k-mer (simplified)
            if (simple_index.marker_id_offsets[i] < simple_index.marker_id_offsets[i + 1]) {
                uint32_t marker_id = simple_index.marker_id_data[simple_index.marker_id_offsets[i]];
                h_kmer_index[i].gene_id = marker_id;  // Use marker_id as gene_id for now
                h_kmer_index[i].species_id = 0;
                h_kmer_index[i].seq_id = marker_id;
                h_kmer_index[i].position = 0;
            }
        }
        
        cudaMalloc(&detector.d_kmer_index, simple_index.num_kmers * sizeof(KmerEntry));
        cudaMemcpy(detector.d_kmer_index, h_kmer_index, 
                   simple_index.num_kmers * sizeof(KmerEntry), cudaMemcpyHostToDevice);
        
        // Create position array
        uint32_t* h_positions = new uint32_t[simple_index.num_kmers];
        for (uint32_t i = 0; i < simple_index.num_kmers; i++) {
            h_positions[i] = i;
        }
        cudaMalloc(&detector.d_kmer_positions, simple_index.num_kmers * sizeof(uint32_t));
        cudaMemcpy(detector.d_kmer_positions, h_positions, 
                   simple_index.num_kmers * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        delete[] h_kmer_index;
        delete[] h_positions;
    }
    
    // Allocate reference sequences
    if (detector.total_ref_length > 0) {
        cudaMalloc(&detector.d_reference_sequences, detector.total_ref_length);
        cudaMalloc(&detector.d_ref_lengths, simple_index.num_sequences * sizeof(uint32_t));
        cudaMalloc(&detector.d_ref_offsets, simple_index.num_sequences * sizeof(uint32_t));
        
        // Copy sequence data
        cudaMemcpy(detector.d_reference_sequences, simple_index.sequence_data,
                   detector.total_ref_length, cudaMemcpyHostToDevice);
        cudaMemcpy(detector.d_ref_lengths, simple_index.seq_lengths,
                   simple_index.num_sequences * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(detector.d_ref_offsets, simple_index.seq_starts,
                   simple_index.num_sequences * sizeof(uint32_t), cudaMemcpyHostToDevice);
        
        // Allocate position weights and mutation masks
        cudaMalloc(&detector.d_position_weights, detector.total_ref_length * sizeof(float));
        cudaMalloc(&detector.d_mutation_masks, detector.total_ref_length);
        
        // Initialize with defaults
        cudaMemset(detector.d_position_weights, 0, detector.total_ref_length * sizeof(float));
        cudaMemset(detector.d_mutation_masks, 0, detector.total_ref_length);
    }
    
    // Allocate mutation info
    cudaMalloc(&detector.d_mutation_info, simple_index.num_markers * sizeof(MutationInfo));
    cudaMalloc(&detector.d_mutation_counts, simple_index.num_sequences * sizeof(uint32_t));
    cudaMemset(detector.d_mutation_counts, 0, simple_index.num_sequences * sizeof(uint32_t));
    
    DEBUG_PRINT("GPU conversion complete: %d kmers, %d sequences", 
                detector.num_kmers, detector.num_sequences);
}

// Updated loadIndex function for FQMutationDetectorCUDA
// COMMENTED OUT - Using the version in fq_mutation_detector.cu that loads binary format
/*
void FQMutationDetectorCUDA::loadIndex(const char* index_path) {
    DEBUG_PRINT("Loading index using simplified loader: %s", index_path);
    
    // Load the simplified index
    SimpleFQIndex simple_index = loadSimpleIndex(index_path);
    
    if (simple_index.num_kmers == 0) {
        std::cerr << "ERROR: Failed to load index or index is empty" << std::endl;
        return;
    }
    
    // Convert to GPU format
    convertToGPUIndex(simple_index, *this);
    
    // Set alignment parameters
    align_params.match_score = 2.0f;
    align_params.mismatch_score = -3.0f;
    align_params.gap_open = -5.0f;
    align_params.gap_extend = -2.0f;
    align_params.mutation_match_bonus = 3.0f;
    align_params.mutation_mismatch_penalty = 2.0f;
    align_params.min_alignment_score = 50.0f;
    align_params.min_identity = 0.8f;
    
    // Clean up
    freeSimpleIndex(simple_index);
    
    std::cerr << "\n=== Simple Index Loaded Successfully ===" << std::endl;
    std::cerr << "K-mers: " << num_kmers << std::endl;
    std::cerr << "Sequences: " << num_sequences << std::endl;
    std::cerr << "Total reference length: " << total_ref_length << std::endl;
    std::cerr << "========================================\n" << std::endl;
}
*/