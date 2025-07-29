// pangenome_graph_types.h
#ifndef PANGENOME_GRAPH_TYPES_H
#define PANGENOME_GRAPH_TYPES_H

#include <cstdint>
#include <vector>
#include <string>

// --- Core Graph Structures ---

// Represents a segment of sequence in the graph.
struct PangenomeNode {
    uint32_t node_id;
    uint64_t sequence_offset; // Offset in the global sequence buffer
    uint32_t sequence_length;
    uint32_t edge_start_index; // Index into the global edge list
    uint32_t edge_count;
    bool is_amr_gene;         // Flag for AMR gene nodes
    uint32_t amr_gene_id;     // If is_amr_gene, links to your existing gene ID
};

// Represents a labeled path through the graph (a genome or plasmid).
struct PangenomePath {
    uint32_t path_id;
    uint32_t species_id;      // Taxonomy ID (e.g., from NCBI)
    char species_name[128];
    char strain_name[128];
    bool is_plasmid;
    uint32_t node_list_start_index; // Index into the global path-node list
    uint32_t node_count;
};

// --- GPU-Optimized Database Representation ---

// This struct will be copied to the GPU and hold pointers to device memory.
struct PangenomeGraphDevice {
    // Node information
    PangenomeNode* d_nodes;
    uint32_t num_nodes;

    // Concatenated sequences of all nodes
    char* d_sequences;
    uint64_t total_sequence_length;

    // Graph topology (CSR format)
    uint32_t* d_edge_list; // Stores destination node_id for each edge
    uint32_t num_edges;

    // Path information
    PangenomePath* d_paths;
    uint32_t num_paths;
    uint32_t* d_path_node_lists; // Concatenated list of node_ids for all paths

    // Minimizer index for seeding
    uint64_t* d_minimizer_hashes;
    uint32_t* d_minimizer_offsets;
    uint32_t* d_minimizer_positions; // Encoded (node_id << 16 | position_in_node)
    uint32_t num_minimizers;
};

// --- Analysis Output ---

// Extends the existing AMRHit to include species attribution.
struct AttributedAMRHit {
    // Inherits fields from your existing AMRHit
    uint32_t read_id;
    uint32_t gene_id;
    float identity;
    float coverage;
    // ... other fields from AMRHit

    // New attribution fields
    uint32_t species_id;      // The attributed species taxonomy ID
    uint32_t path_id;         // The specific path (genome/plasmid) the hit was found on
    float attribution_confidence; // A score from 0.0 to 1.0
};


#endif // PANGENOME_GRAPH_TYPES_H
