# Phase 1 Enhanced K2-like Classifier

## Overview

The Phase 1 Enhanced Classifier builds upon the existing Kraken2-style paired-end GPU classifier with the following enhancements:

### Key Features

1. **Two-Stage Filtering**
   - Stage 1: Primary confidence threshold (default: 0.1)
   - Stage 2: Secondary confidence threshold with enhanced validation (default: 0.3)

2. **Phylogenetic Consistency Validation**
   - Weighted phylogenetic consistency scoring
   - Detection of taxonomic conflicts
   - Multiple weighting methods (exponential decay, linear, hybrid)

3. **Multi-Minimizer Validation**
   - Requires multiple distinct minimizers for robust classification
   - Configurable minimizer diversity requirements

4. **K-mer Position Tracking** (Forward Compatibility)
   - Tracks k-mer positions in reads for future post-hoc analysis
   - Stores minimizer hashes, positions, and confidence weights
   - Prepares data for future coverage validation

5. **Enhanced Confidence Calculation**
   - Composite confidence score using multiple factors:
     - Primary k-mer vote ratio (40%)
     - Phylogenetic consistency (30%)
     - Minimizer diversity (20%)
     - Taxonomic specificity (10%)

## Architecture

The enhanced classifier uses composition with the existing `PairedEndGPUKrakenClassifier`:

```
Phase1EnhancedClassifier
    ├── PairedEndGPUKrakenClassifier (base functionality)
    ├── Enhanced GPU kernels
    ├── K-mer tracking system
    └── Phylogenetic validation
```

## Building

Add the CMakeLists additions from `CMakeLists_addition_cc.txt` to the profiler's CMakeLists.txt:

```bash
cd /path/to/biogpu
# Add the enhanced classifier to CMakeLists.txt
# Then build:
mkdir build && cd build
cmake ..
make gpu_kraken_enhanced_pipeline test_enhanced_classifier
```

## Usage

### Basic Usage

```bash
# Test the enhanced classifier
./test_enhanced_classifier /path/to/database /path/to/reads.fastq

# With custom parameters
./test_enhanced_classifier /path/to/database /path/to/reads.fastq \
    --confidence 0.05 \
    --confidence2 0.25 \
    --min-kmers 5 \
    --min-minimizers 2

# Compare with base classifier
./test_enhanced_classifier /path/to/database /path/to/reads.fastq --compare
```

### Parameters

- `--confidence <float>`: Primary confidence threshold (Stage 1)
- `--confidence2 <float>`: Secondary confidence threshold (Stage 2)
- `--min-kmers <int>`: Minimum k-mers required for classification
- `--min-minimizers <int>`: Minimum distinct minimizers required
- `--no-phylo`: Disable phylogenetic validation
- `--no-tracking`: Disable k-mer position tracking
- `--max-reads <int>`: Maximum reads to process
- `--compare`: Compare results with base classifier

## Enhanced Classification Process

1. **Initial K-mer Processing**
   - Extract minimizers from reads
   - Look up LCA in hash table
   - Track distinct minimizers
   - Record k-mer positions (if tracking enabled)

2. **Stage 1 Validation**
   - Check primary confidence threshold
   - Verify minimum k-mer count
   - Verify minimizer diversity

3. **Stage 2 Enhanced Validation**
   - Calculate phylogenetic consistency score
   - Detect taxonomic conflicts
   - Compute composite confidence score

4. **Final Classification**
   - High confidence: Passed both stages
   - Low confidence: Passed Stage 1 only
   - Unclassified: Failed Stage 1

## Output Format

The enhanced classifier provides detailed results including:

```
Phase1ClassificationResult {
    taxon_id                        // Final classification
    primary_confidence              // Stage 1 confidence
    secondary_confidence            // Stage 2 enhanced confidence
    passed_stage1/stage2           // Stage validation results
    is_high_confidence_call        // Overall confidence level
    classified_kmers/total_kmers   // K-mer statistics
    distinct_minimizers_found      // Minimizer diversity
    phylogenetic_consistency_score // Phylo validation score
    has_taxonomic_conflicts        // Conflict detection
    kmer_tracking {                // Position tracking data
        minimizer_hashes
        read_positions
        contributing_taxa
        confidence_weights
    }
}
```

## Performance Considerations

- K-mer tracking adds memory overhead (~400 bytes per read with default settings)
- Phylogenetic validation adds computational overhead but improves accuracy
- Two-stage filtering reduces false positives at the cost of some sensitivity
- GPU memory usage scales with batch size and tracking parameters

## Integration with Existing Pipeline

The enhanced classifier can be used as a drop-in replacement or alongside the existing classifier:

```cpp
// Option 1: Standalone enhanced classifier
Phase1EnhancedParams params;
Phase1EnhancedClassifier classifier(params);
classifier.load_database(db_dir);
auto results = classifier.classify_enhanced(reads);

// Option 2: Use with existing classifier
PairedEndGPUKrakenClassifier* base = new PairedEndGPUKrakenClassifier();
Phase1EnhancedClassifier enhanced(base, params);
// ... use enhanced classifier
```

## Future Enhancements (Phase 2)

The current implementation is forward-compatible with planned Phase 2 features:
- Coverage validation using tracked k-mer positions
- Multi-genome mapping analysis
- Post-hoc positional coverage assessment
- Advanced contamination detection

## Troubleshooting

1. **"Enhanced classification requires database pointers to be set"**
   - The enhanced classifier needs access to GPU database structures
   - Use `set_database_pointers()` method if using custom integration

2. **High memory usage**
   - Reduce `max_kmers_to_track` parameter
   - Disable k-mer tracking with `--no-tracking`
   - Process smaller batches

3. **Different results from base classifier**
   - This is expected due to enhanced validation
   - Use `--compare` flag to analyze differences
   - Adjust confidence thresholds if needed