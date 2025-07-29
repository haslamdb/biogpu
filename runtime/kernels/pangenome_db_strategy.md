# BioGPU Pangenome Database: Build Strategy & Design Notes

This document summarizes the key architectural decisions and resource estimates for constructing the pangenome graph database, which will serve as the reference for the BioGPU alignment engine.

## 1. Objective

The goal is to build a comprehensive, offline database that represents the known bacterial pangenome. This database will be a high-resolution sequence graph, enabling a novel GPU-based alignment algorithm to perform simultaneous gene detection and species attribution from metagenomic reads.

## 2. Database Construction Workflow

The database construction is a one-time, offline process performed by a new tool, `build_pangenome_db`.

### Inputs:
1.  **Bacterial Genomes:** All "complete" genomes from the NCBI RefSeq database.
2.  **Plasmids:** A comprehensive collection of known plasmid sequences (e.g., NCBI Plasmid RefSeq).
3.  **Target Genes:** The existing BioGPU AMR gene database to ensure these specific sequences are explicitly represented and annotated.

### Core Process:
1.  **Sequence Ingestion:** All input FASTA files are read and their sequences are concatenated into a single, massive text string. This text is the foundation of the index.
2.  **Suffix Array & FM-Index Construction:** The most critical step is building a GPU-native FM-Index from the concatenated text. This involves:
    * **Suffix Array (SA) Generation:** Creating an array of all suffixes of the text, sorted lexicographically. This is the most computationally and memory-intensive part of the entire offline build.
    * **Burrows-Wheeler Transform (BWT):** Derived from the SA, this permutation of the text has properties that make it highly compressible and searchable.
    * **Auxiliary Tables (C-table, Occ):** These tables are built from the BWT to enable extremely fast k-mer searching, which is the basis for the aligner's seeding stage.
3.  **Graph Generation & Annotation:** The sequence and its index are used to define the graph structure.
    * **Nodes & Edges:** The builder segments the sequences into nodes and records their adjacencies as edges.
    * **Path & Annotation Data:** The builder creates metadata files that define which nodes belong to which original genome/plasmid (Paths) and flags specific nodes that correspond to AMR genes.

## 3. The Homology Parameter: A Critical Design Choice

A key decision is how the builder handles similar, but not identical, sequences (e.g., `blaCTX-M-15` vs. `blaCTX-M-19`). This is controlled by a **sequence identity threshold** that dictates whether to merge sequences into a shared node or keep them separate.

### The Trade-Off:
* **High Aggregation (e.g., 98% identity):** Merges similar sequences into a shared path with small "bubbles" at variant sites. This creates a smaller, less complex graph but makes it harder to distinguish between close variants at the alignment stage.
* **Low Aggregation (e.g., 100% identity):** Only merges perfectly identical segments. Similar genes are represented as distinct, parallel paths. This creates a much larger, more complex graph but preserves all variant information with maximum fidelity.

### Recommendation for BioGPU:
For the proposed alignment method, we **strongly recommend the Low Aggregation / High-Resolution approach**. The philosophy is to build a graph that is a perfectly faithful representation of the input sequences. This provides the aligner with the maximum possible information. The task of resolving ambiguity (e.g., a read that could map to either `CTX-M-15` or `CTX-M-19`) is then handled by the downstream alignment algorithm and the EM statistical model, which can use global information from all reads to make the most probable assignment.

## 4. Resource & Time Estimates

These are order-of-magnitude estimates for building the database using all complete bacterial RefSeq genomes and plasmids on the specified workstation.

* **Input Data Size (Download):**
    * Complete Bacterial Genomes (RefSeq): **~150-200 GB**
    * Complete Plasmid Database (RefSeq): **~5-10 GB**
    * Total Raw Sequence: **~50-70 billion base pairs**.
* **Offline Build Time:**
    * This is a CPU- and RAM-bound process. With all 64 cores utilized, the Suffix Array construction will be the longest step.
    * **Estimated Build Time: 48 - 96 hours.** This is highly dependent on the efficiency of the Suffix Array construction algorithm (e.g., SA-IS). Using the 8TB NVMe drive for intermediate build files will accelerate I/O and likely place the build time in the lower end of this range.
* **Final Database Size (On-Disk Footprint):**
    * The final, compressed graph and FM-Index will be very large.
    * **Estimated Final Size: 1.0 - 1.5 TB.**
    * The largest components will be the Suffix Array and the checkpointed Occurrence table of the FM-Index. This database should be stored on the 8TB NVMe drive for fastest loading during the analysis phase.

## 5. Hardware Strategy

The specified workstation is exceptionally well-suited for this project. The recommended strategy for utilizing its components is as follows:

* **CPU (64-core Threadripper Pro):** Dedicate all cores to the offline `build_pangenome_db` tool. This is essential for processing the tens of billions of base pairs in a reasonable timeframe.
* **RAM (512 GB):** This is the key enabler for the offline build. It allows the large, intermediate data structures (like the full Suffix Array) to be held in memory, making the construction of a large-scale graph feasible on a single machine.
* **Storage (NVMe & SATA):**
    * **8TB NVMe SSD:** This is the primary "working drive." Use it to store the final 1.0-1.5 TB database. All intermediate files generated during the build process should also be written here to maximize I/O speed.
    * **4TB NVMe SSD (OS):** Keep this drive dedicated to the OS and applications for system responsiveness.
    * **100TB SATA Storage:** Use this for long-term archival. The initial multi-terabyte download of the full RefSeq database should reside here, as well as backups of the final pangenome database.
* **GPU (RTX A6000 & A5000):**
    * **RTX A6000 (48 GB VRAM):** This is the primary GPU for the online analysis. Its 48 GB of VRAM defines the upper limit on the size of the graph index that can be loaded for real-time alignment, making it the target for the full-scale pangenome.
    * **RTX A5000 (24 GB VRAM):** An excellent secondary GPU. Ideal for developing and testing the alignment kernels with a smaller, curated "mini" version of the pangenome graph without tying up the A6000.
