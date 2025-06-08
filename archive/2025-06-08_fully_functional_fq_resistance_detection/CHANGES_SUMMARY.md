# Changes Summary - 2025-06-08

## Critical Fixes Applied

### 1. Fixed QRDR Detection
**Problem**: The `is_qrdr_alignment` field was always `false` because it wasn't being set on the host side.

**Solution**: Added QRDR detection logic in `clean_resistance_pipeline_main.cpp` after protein matches are retrieved:
```cpp
// Check if this alignment covers QRDR region
bool is_qrdr_alignment = false;
if (g_fq_resistance_db) {
    // Check if any position in the alignment range is a QRDR position
    for (int pos = match.ref_start; pos < match.ref_start + match.match_length; pos++) {
        if (g_fq_resistance_db->isQRDRPosition(gene_name, pos)) {
            is_qrdr_alignment = true;
            break;
        }
    }
}
match.is_qrdr_alignment = is_qrdr_alignment;
```

### 2. Fixed Gene/Species ID Mapping
**Problem**: After disabling `fq_mutation_reporter`, the metadata mappings weren't being loaded because of a check for the reporter being non-null.

**Solution**: Modified `loadMappingsForReporter()` to always load mappings:
```cpp
void loadMappingsForReporter() {
    // Always load mappings since we need them for the main pipeline too
    if (protein_db_path.empty()) return;
    // ... rest of loading code
}
```

### 3. Created New Clinical Report Generator
**Problem**: The old `fq_mutation_reporter` was filtering out non-QRDR mutations, preventing comprehensive mutation analysis.

**Solution**: Created `clinical_fq_report_generator.cpp` that:
- Reports all mutations, not just QRDR ones
- Provides clinical interpretation with confidence scores
- Generates HTML, JSON, and text reports
- Clearly distinguishes between known FQ resistance and QRDR mutations

### 4. Added Performance Metrics
Added reads/second calculation to the pipeline summary:
```cpp
double reads_per_second = stats.total_reads / static_cast<double>(duration.count());
std::cout << "Performance: " << std::fixed << std::setprecision(0) 
          << reads_per_second << " reads/second\n";
```

## Results After Fixes

- **Species/Gene Detection**: Working correctly (e.g., "Escherichia_coli", "gyrA")
- **QRDR Detection**: Properly identifies QRDR alignments
- **FQ Resistance Detection**: Successfully detects known mutations (e.g., D87G)
- **Clinical Reporting**: Generates comprehensive reports with proper confidence levels
- **Performance**: ~16,667-21,739 reads/second on NVIDIA TITAN Xp

## Files Modified

1. `runtime/kernels/resistance/clean_resistance_pipeline_main.cpp`
2. `runtime/kernels/resistance/clinical_fq_report_generator.cpp` (new file)
3. `CMakeLists.txt` (added clinical report generator)

## Test Results

On 1M synthetic reads:
- Total mutations detected: 18,523
- FQ resistance mutations: 658
- High confidence resistance detected in E. coli
- D87G mutation in gyrA properly identified as high-level FQ resistance