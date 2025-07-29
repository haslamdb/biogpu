# GPU Pangenome Alignment: The Synchronous Wavefront Method

Given unlimited computational tools, the fastest and most accurate way to align reads to a pangenome graph would be to move away from the traditional "seed-and-extend" paradigm. Instead, we can treat the pangenome graph as a circuit and the reads as a signal. The goal is to find the path through the circuit that **resonates** most strongly with the signal. This method, **Synchronous Wavefront Alignment**, is designed for the massively parallel architecture of a modern GPU.

## 1. The Core Idea: Signal Resonance

Imagine the entire pangenome graph is "etched" into the GPU's memory. Each node is a transistor and each edge is a wire. Instead of processing one read at a time, we apply the "signal" of **all reads simultaneously** to the circuit. We then watch to see which paths through the graph light up the brightest.

This approach consists of three main stages, all performed in parallel on the GPU.

## 2. Stage 1: Massively Parallel Seeding via GPU-Native FM-Index

The first step is to find all possible starting points for our signal. Instead of using minimizers, which are a heuristic to reduce the number of seeds, we find **all exact k-mer matches** for all reads at once.

* **How it Works:** We pre-compile the entire DNA sequence of the pangenome graph into a **GPU-native FM-index**. This is a compressed data structure that allows for extremely fast substring searches.
* **The Kernel:** A single, massive CUDA kernel takes all k-mers from all reads in the batch. Each thread processes a k-mer and queries the FM-index.
* **Output:** The result is a huge list of "activation points": `(read_id, position_in_read, node_id, position_in_node)`. This tells us every single place in the graph where a piece of a read has a perfect match.

## 3. Stage 2: Synchronous Wavefront Propagation

This is the core of the algorithm and where it diverges most from traditional methods. The alignment proceeds in synchronous "ticks," like a processor's clock cycle.

* **Initialization:** At Tick=0, every node in the graph that was identified in the seeding stage is "activated" with an initial score.
* **The Propagation Kernel:** In each subsequent tick, a single CUDA kernel is launched where every thread is responsible for a single **activated node**. The kernel does the following:
    1.  It looks at the next base in the corresponding sequencing read.
    2.  It checks all of its outgoing edges in the graph.
    3.  If the sequence of an edge and its destination node matches the next base(s) of the read, it "propagates" its score to that neighboring node.
    4.  The neighbor node becomes activated for the next tick.
* **The Wavefront:** This process creates a "wavefront" of activation scores that flows through the graph, cycle by cycle. The wavefronts naturally follow the paths where the read sequence matches the graph sequence. If a read doesn't match a particular branch, that part of the wavefront simply dies out.
* **Paired-End Pincer Movement:** For paired-end reads, we initiate two wavefronts simultaneously—one for Read 1 and one for Read 2. These wavefronts travel towards each other. A successful alignment is declared when two wavefronts originating from the same read pair meet on a valid path in the graph. This provides immense accuracy and immediately resolves most repeats.

## 4. Stage 3: Path Resonance and Traceback

After a set number of ticks (e.g., 150, the length of a read), some paths through the graph will have accumulated a very high score—they have "resonated" with the read signal.

* **Scoring:** The score isn't just a count. It's a quality-weighted score that can incorporate base quality scores from the FASTQ file and use a scoring matrix (like Smith-Waterman) for the local comparisons at each tick.
* **Traceback:** A final kernel identifies the nodes with the highest terminal scores. It then performs a parallel traceback, following the chain of predecessors that contributed the highest scores at each tick. This reconstructs the single best-scoring path for each read pair.

### Why This is the Fastest and Most Accurate Method

* **Speed:** It leverages the GPU's architecture perfectly. The wavefront propagation is a highly structured, data-parallel operation with predictable memory access patterns, minimizing the random memory lookups that plague traditional alignment. Every stage is a massive parallel operation across thousands of cores.
* **Accuracy:** It is inherently more accurate because it explores all possible alignments simultaneously. The best alignment isn't just the first one found; it's the one that emerges as the global optimum after considering all evidence. It naturally handles complex structural variations, insertions, and deletions, as these are simply different branches in the graph for the wavefront to explore.

This "Synchronous Wavefront Alignment" method represents a paradigm shift from the serial, CPU-centric logic of seed-and-extend to a massively parallel, signal-processing approach tailored for the future of computational hardware.
