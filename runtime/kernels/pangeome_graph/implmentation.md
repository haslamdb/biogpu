BioGPU Pangenome Graph Module: Implementation Plan
This document outlines the design and implementation strategy for extending the BioGPU pipeline with a reference pangenome graph module. The goal is to enable highly accurate, GPU-accelerated species-level attribution for detected antibiotic resistance (AMR) genes.

1. Core Concept
The fundamental shift is from matching reads against a linear list of AMR genes to mapping reads onto a comprehensive, pre-built graph. This graph represents the known universe of bacterial genomes, plasmids, and their associated AMR genes. By leveraging paired-end read information, we can determine the most likely "path" a read pair takes through the graph, simultaneously identifying the AMR gene and the species carrying it.

2. System Architecture
The system is divided into two main components: an Offline Graph Builder and the Online Analysis Pipeline.

2.1. Offline: The Graph Compiler (build_pangenome_db)
This is a new, standalone C++ tool responsible for creating the graph database. It is run once to generate the necessary binary files for the analysis pipeline.

Inputs:

A curated set of complete bacterial reference genomes (FASTA format).

A comprehensive database of known plasmid sequences (FASTA format).

The existing NCBI AMR gene sequences used by BioGPU (FASTA format).

Process:

Node Segmentation: Deconstruct all input sequences into unique, overlapping segments (nodes).

Edge Creation: Create directed edges between nodes that are adjacent in the original sequences.

Path Annotation: Store the original genome/plasmid sequences as explicit paths (ordered lists of node IDs). Each path is annotated with metadata (e.g., species_id, strain_name, is_plasmid). AMR gene nodes are also specially flagged.

Minimizer Indexing: Generate a global minimizer index for all node sequences to enable rapid seeding.

Outputs: A directory containing several GPU-optimized binary files:

graph_nodes.bin: Node metadata.

graph_topology.csr: Adjacency information in Compressed Sparse Row (CSR) format.

graph_paths.bin: Path metadata and node lists.

graph_sequences.bin: Concatenated DNA sequences of all nodes.

graph_minimizer.idx: The global minimizer index.

2.2. Online: The Analysis Pipeline (amr_detection_pipeline)
This involves modifying the existing AMRDetectionPipeline to use the new graph database.

Input: Paired-end FASTQ files.

Process:

Initialization: The pipeline will load the pre-compiled pangenome graph files into GPU memory.

Minimizer Seeding: For each read, generate minimizers (using the existing generate_minimizers_kernel) and query the graph's minimizer index to find all potential node locations (seed hits).

Path Attribution (New CUDA Kernel): This is the core of the new functionality. For each read pair:

Collect all seed hits for both R1 and R2.

Identify "anchor pairs" where one read hits a node of interest (e.g., an AMR gene) and the mate hits a distinct, species-specific node.

Perform a rapid, path-constrained graph traversal between the two anchor nodes.

If a valid path is found within the expected fragment length, a high-confidence attribution is made.

A new AttributedAMRHit struct, containing gene_id and species_id, is generated.

EM Abundance Quantification: The existing EM algorithm is adapted to work on (gene_id, species_id) tuples. This resolves ambiguous mappings where a read pair is consistent with multiple species carrying the same gene, providing a final, quantified abundance for each species-specific resistance gene.

Output: The final report will now include species attribution for each detected resistance gene, along with its abundance.

3. Implementation Steps
Define Data Structures: Finalize the C++/CUDA structs for graph nodes, paths, and the new AttributedAMRHit.

Build the Offline Compiler: Create the build_pangenome_db tool. This will be CPU-based C++ code.

Modify the Pipeline:

Update AMRDetectionPipeline::initialize to load the new graph files.

Create a new GraphSeeding class/module to replace the bloom filter logic.

Implement the path_attribution_kernel.cu CUDA kernel.

Adapt the EM algorithm in amr_detection_pipeline.cpp.

Update Reporting: Modify the reporting tools to display the new species-level information.