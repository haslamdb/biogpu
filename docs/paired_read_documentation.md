# AMR Detection Pipeline Workflow Documentation

## Enhanced Paired-End Read Handling Implementation
*Date: July 9, 2025*

### Overview

For algorithms designed for accurate quantification, especially in the context of highly similar AMR genes, the way paired-end reads are handled is critical. The pipeline has been updated to treat R1 and R2 as a pair throughout the alignment/assignment process, rather than merging them into a single, longer read. This approach follows the methodology used by tools like Kallisto and Salmon, leveraging fragment length distribution and relative orientation information for improved accuracy.

### Why the Previous Merging Approach Was Suboptimal

The previous `mergePairedReads` function concatenated R1 and R2 with 'N's, which caused several issues:

1. **Loss of Insert Size Information**: The 'N's artificially inflated the fragment length and obscured the true insert size, which is a powerful statistical cue for assigning reads.

2. **Loss of Relative Orientation**: The merged read lost the explicit information about R1 and R2's individual orientations and their relative positions.

3. **Alignment Challenges**: Aligners might struggle with the 'N' gaps or treat them as long indels, potentially leading to less accurate initial alignments or hits.

4. **Reduced Specificity**: Multi-mapping reads that could be resolved by their pair were instead treated as single, ambiguous entities.

### Implementation Details

#### 1. Data Structure Enhancements

**New ReadPairData Structure**:
```cpp
struct ReadPairData {
    std::string read1_seq;
    std::string read2_seq;
    std::string read1_id;
    std::string read2_id;
    uint32_t pair_index;    // Original index in the batch
};
```

**Enhanced AMRHit Structure**:
```cpp
struct AMRHit {
    // ... existing fields ...
    // Enhanced paired-end tracking
    uint32_t pair_id;       // Original read pair ID
    bool is_read2;          // false for R1, true for R2
    uint32_t mate_read_id;  // ID of the mate read
    float pair_score;       // Combined score if concordant
};
```

**PairedReadAssignment Structure for EM**:
```cpp
struct PairedReadAssignment {
    uint32_t pair_id;
    uint32_t r1_read_id;
    uint32_t r2_read_id;
    
    // Compatible genes (where both reads can map consistently)
    std::vector<uint32_t> compatible_genes;
    
    // Individual hit information
    std::map<uint32_t, float> r1_gene_scores;  // gene_id -> score
    std::map<uint32_t, float> r2_gene_scores;  // gene_id -> score
    
    // Paired information
    std::map<uint32_t, float> concordance_scores;     // gene_id -> combined score
    std::map<uint32_t, int32_t> fragment_lengths;     // gene_id -> fragment length
    std::map<uint32_t, float> fragment_probabilities; // gene_id -> P(fragment|gene)
    
    // Assignment probabilities for each gene
    std::map<uint32_t, float> assignment_probabilities;
};
```

#### 2. Processing Workflow

**Main Processing Changes**:
1. Reads are no longer merged in `processSamplePaired`
2. New `processPairedBatch` method processes R1 and R2 separately
3. Paired information is maintained throughout the pipeline

**Key Methods**:
- `processPairedBatch()`: Main entry point for paired-end processing
- `updateHitsWithPairedInfo()`: Populates paired-end fields in AMRHit structures
- `buildEnhancedPairAlignments()`: Builds pair alignment information including fragment lengths

#### 3. Fragment Length Distribution

The pipeline now learns the fragment length distribution from the data:

```cpp
void estimateFragmentLengthDistribution(const std::vector<ReadPairAlignment>& alignments)
```

- Collects fragment lengths from high-confidence concordant pairs
- Calculates mean and standard deviation
- Removes outliers using z-score filtering (2 standard deviations)
- Used to calculate fragment length probabilities for scoring

#### 4. Paired-End Aware EM Algorithm

The EM algorithm has been enhanced to handle paired-end reads as fundamental units:

**E-Step (Expectation)**:
```cpp
void updatePairedAssignmentProbabilities()
```
Calculates assignment probabilities considering:
- Individual alignment scores for R1 and R2
- Concordance status (both reads mapping to same gene)
- Fragment length probability
- Gene abundance
- Sequence similarity between genes

**M-Step (Maximization)**:
```cpp
void updateGeneAbundancesFromPairs()
```
- Updates gene abundances based on paired assignments
- Each pair counts as 2 reads
- Handles orphaned reads separately

**Probability Calculation**:
```
P(pair|gene) = alignment_weight × abundance_weight × similarity_boost × fragment_probability
```

Where:
- `alignment_weight` includes concordance bonus/penalty
- `abundance_weight` is based on current gene abundance estimate
- `similarity_boost` accounts for sequence similarity between candidate genes
- `fragment_probability` is based on learned fragment length distribution

#### 5. Configuration Parameters

New configuration options for paired-end handling:

```cpp
// Paired-end concordance parameters
float concordance_bonus = 2.0f;         // Bonus for concordant pairs
float discord_penalty = 0.5f;           // Penalty for discordant pairs
int expected_fragment_length = 500;     // Expected fragment length
int fragment_length_stddev = 100;       // Standard deviation for fragment length
float max_fragment_length_zscore = 3.0f; // Max z-score for valid fragments
```

#### 6. Automatic Detection

The pipeline automatically detects paired-end data and selects the appropriate EM algorithm:

```cpp
// In resolveAmbiguousAssignmentsEM()
if (has_paired_data) {
    buildPairedReadAssignments(accumulated_hits);
    runPairedEndEM();
} else {
    buildReadAssignments(accumulated_hits);
    runKallistoStyleEM();
}
```

### Benefits of the New Approach

1. **Improved Specificity**: Multi-mapping reads can be resolved using paired-end evidence
2. **Better Quantification**: Fragment length information helps distinguish between similar genes
3. **Adaptive Learning**: Fragment length distribution is learned from the data itself
4. **Preservation of Information**: All paired-end information is maintained throughout the pipeline
5. **Concordance Scoring**: Concordant pairs receive appropriate bonuses, improving assignment accuracy

### Usage

The new paired-end handling is automatically enabled when processing paired-end reads. No changes to the command-line interface are required. The pipeline will:

1. Detect paired-end input
2. Process R1 and R2 separately
3. Learn fragment length distribution from concordant pairs
4. Use paired-end aware EM for quantification
5. Report concordance statistics in the output

### Performance Considerations

- Memory usage is slightly higher due to tracking both R1 and R2 separately
- Processing time may increase slightly due to additional paired-end calculations
- The benefits in accuracy, especially for highly similar AMR genes, outweigh the performance costs

### Future Enhancements

Potential improvements for future versions:
1. Support for different library preparation protocols (e.g., strand-specific)
2. More sophisticated fragment length models (e.g., mixture models)
3. Integration of insert size information from alignment BAM files
4. Adaptive concordance scoring based on data quality