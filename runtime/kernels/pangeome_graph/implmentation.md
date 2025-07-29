# BioGPU Pangenome Graph Module: Implementation Plan

This document outlines the design and implementation strategy for extending the BioGPU pipeline with a reference pangenome graph module. The goal is to enable highly accurate, GPU-accelerated species-level attribution for detected antibiotic resistance (AMR) genes.

## 1. Core Concept

The fundamental shift is from matching reads against a linear list of AMR genes to mapping reads onto a comprehensive, pre-built graph. This graph represents the known universe of bacterial genomes, plasmids, and their associated AMR genes. By leveraging paired-end read information, we can determine the most likely "path" a read pair takes through the graph, simultaneously identifying the AMR gene and the species carrying it.

## 2. System Architecture

The system is divided into two main components: an **Offline Graph Builder** and the **Online Analysis Pipeline**.

### 2.1. Offline: The Graph Compiler (`build_pangenome_db`)

This is a new, standalone C++ tool responsible for creating the graph database. It is run once to generate the necessary binary files for the analysis pipeline.

* **Inputs:**
    * A curated set of complete bacterial reference genomes (FASTA format).
    * A comprehensive database of known plasmid sequences (FASTA format).
    * The existing NCBI AMR gene sequences used by BioGPU (FASTA format).
* **Process:**
    1.  **Node Segmentation:** Deconstruct all input sequences into unique, overlapping segments (nodes).
    2.  **Edge Creation:** Create directed edges between nodes that are adjacent in the original sequences.
    3.  **Path Annotation:** Store the original genome/plasmid sequences as explicit paths (ordered lists of node IDs). Each path is annotated with metadata (e.g., `species_id`, `strain_name`, `is_plasmid`). AMR gene nodes are also specially flagged.
    4.  **Minimizer Indexing:** Generate a global minimizer index for all node sequences to enable rapid seeding.
* **Outputs:** A directory containing several GPU-optimized binary files:
    * `graph_nodes.bin`: Node metadata.
    * `graph_topology.csr`: Adjacency information in Compressed Sparse Row (CSR) format.
    * `graph_paths.bin`: Path metadata and node lists.
    * `graph_sequences.bin`: Concatenated DNA sequences of all nodes.
    * `graph_minimizer.idx`: The global minimizer index.