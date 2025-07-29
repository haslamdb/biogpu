// wavefront_seeding_kernel.cu
#include "wavefront_alignment_types.h"

/**
 * @brief Performs massively parallel seeding using a GPU-native FM-Index.
 *
 * This kernel is the first step in the Synchronous Wavefront Alignment. Each thread
 * takes a k-mer from a read, queries the FM-Index to find all exact matches in the
 * pangenome graph, and creates an initial "WavefrontState" for each match. This
 * collection of states represents "Tick 0" of the alignment process.
 *
 * @param reads A character array of all concatenated reads for the batch.
 * @param read_offsets Offsets for the start of each read in the `reads` array.
 * @param read_lengths Length of each read.
 * @param num_reads Total number of reads in the batch.
 * @param fm_index A pointer to the pre-built FM-Index on the GPU.
 * @param k The k-mer size to use for seeding.
 * @param initial_wavefront_states_out Output array to store the generated "Tick 0" states.
 * @param state_counts_out Output array to store the number of initial states generated for each read.
 */
__global__ void initial_seeding_kernel(
    const char* reads,
    const uint32_t* read_offsets,
    const uint32_t* read_lengths,
    const int num_reads,
    const FMIndexDevice* fm_index,
    const int k,
    WavefrontState* initial_wavefront_states_out,
    uint32_t* state_counts_out)
{
    // Each thread block can be responsible for one read.
    const int read_id = blockIdx.x;
    if (read_id >= num_reads) {
        return;
    }

    // Each thread within the block processes a starting position in the read.
    const int start_pos_in_read = threadIdx.x;
    const int read_len = read_lengths[read_id];
    const char* read_seq = &reads[read_offsets[read_id]];

    if (start_pos_in_read + k > read_len) {
        return;
    }

    // 1. Extract the k-mer from the read.
    const char* kmer = &read_seq[start_pos_in_read];

    // 2. Query the FM-Index to find all occurrences of this k-mer.
    //    (This is a complex operation in a real implementation).
    //    For this stub, we'll imagine it returns a list of locations.
    //
    //    uint32_t match_locations[10]; // Dummy array for matches
    //    int num_matches = fm_index_query(fm_index, kmer, k, match_locations);
    int num_matches = 1; // Dummy value

    // 3. For each match found, create an initial WavefrontState.
    for (int i = 0; i < num_matches; ++i) {
        // Decode the match location to get node_id and position_in_node
        // uint64_t pangenome_pos = match_locations[i];
        // uint32_t node_id = decode_node_id(pangenome_pos);
        // uint32_t pos_in_node = decode_pos_in_node(pangenome_pos);
        uint32_t node_id = 123;     // Dummy value
        uint32_t pos_in_node = 45;  // Dummy value

        // Atomically get an index in the output array.
        // This is a simplified approach; a real implementation would use a more
        // efficient parallel allocation strategy.
        uint32_t output_idx = atomicAdd(&state_counts_out[read_id], 1);

        // Populate the "Tick 0" state.
        WavefrontState& state = initial_wavefront_states_out[read_id * 1000 + output_idx]; // Assuming max 1000 seeds/read
        state.read_id = read_id;
        state.node_id = node_id;
        state.pos_in_node = pos_in_node + k -1; // Position at the end of the k-mer
        state.pos_in_read = start_pos_in_read + k - 1;
        state.score = 1.0f; // Initial score
        state.predecessor_idx = 0; // No predecessor at Tick 0
    }
}

