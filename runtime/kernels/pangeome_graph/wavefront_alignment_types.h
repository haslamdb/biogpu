// wavefront_alignment_types.h
#ifndef WAVEFRONT_ALIGNMENT_TYPES_H
#define WAVEFRONT_ALIGNMENT_TYPES_H

#include <cstdint>
#include "pangenome_graph_types.h" // Assuming this file exists from our previous plan

// --- FM-Index Structures (GPU-Native) ---

// A simplified representation of an FM-Index on the GPU.
// A real implementation would be more complex, likely involving bit-compressed arrays.
struct FMIndexDevice {
    const char* d_bwt; // Burrows-Wheeler Transformed string of the pangenome
    const uint32_t* d_suffix_array; // To locate positions in the original text
    const uint32_t* d_C_table;      // Checkpoint table for faster lookups
    const uint32_t* d_occ;          // Occurrence table
    uint64_t text_length;
};


// --- Wavefront Structures ---

// Represents the state of a single active alignment wavefront at a given tick.
// This is the core data structure that gets passed from one tick to the next.
struct WavefrontState {
    uint32_t read_id;       // Which read this wavefront belongs to
    uint32_t node_id;       // Current node in the pangenome graph
    uint32_t pos_in_node;   // Current position within that node's sequence
    uint16_t pos_in_read;   // How far along the read we've matched (this is the "tick")
    float score;            // The accumulated alignment score
    uint32_t predecessor_idx; // Index to the state in the previous tick's array that generated this one
};

#endif // WAVEFRONT_ALIGNMENT_TYPES_H
