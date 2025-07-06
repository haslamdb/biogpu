// runtime/kernels/resistance/gpu_diagnostic_adapters.cu
// Adapters for clinical diagnostic mutation detection in bacterial samples
// These adapters enable identification of antibiotic resistance mutations for treatment decisions

#include "gpu_diagnostic_adapters.h"
#include "fq_mutation_detector.cuh"
#include "../../common/gpu/gpu_sequence_buffer.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace BioGPU {

// Device constant memory for genetic code (bacterial code table 11)
__constant__ char d_genetic_code[64];

// Helper function to initialize genetic code on device
void initializeGeneticCode() {
    static bool initialized = false;
    if (!initialized) {
        // Standard bacterial genetic code (table 11)
        const char genetic_code[64] = {
            'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
            'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
            'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
            '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
        };
        cudaMemcpyToSymbol(d_genetic_code, genetic_code, 64 * sizeof(char));
        initialized = true;
    }
}

// Adapter kernel: Convert sequence data format for mutation detection
__global__ void convertSequenceFormat(
    const char* d_sequences,
    const int* d_offsets,
    const int* d_lengths,
    int num_sequences,
    char* d_legacy_sequences,
    int* d_legacy_lengths,
    int max_length) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_sequences) return;
    
    int start = d_offsets[tid];
    int length = d_lengths[tid];
    
    // Copy sequence to legacy format (padded fixed-length)
    int legacy_offset = tid * max_length;
    for (int i = 0; i < length && i < max_length; i++) {
        d_legacy_sequences[legacy_offset + i] = d_sequences[start + i];
    }
    
    // Pad remaining with 'N'
    for (int i = length; i < max_length; i++) {
        d_legacy_sequences[legacy_offset + i] = 'N';
    }
    
    d_legacy_lengths[tid] = min(length, max_length);
}

// Main mutation detection adapter for clinical diagnostics
void launchMutationDetectionKernel_v2(
    GPUSequenceBuffer& buffer,
    const FQMutationDatabase& known_mutations_db,
    DiagnosticResults& detected_mutations,
    cudaStream_t stream) {
    
    // Initialize genetic code if needed
    initializeGeneticCode();
    
    // Get buffer dimensions
    int num_sequences = buffer.getCurrentSequences();
    int max_length = 300;  // Standard read length for diagnostic sequencing
    
    // Allocate temporary buffer for legacy format
    char* d_legacy_sequences;
    int* d_legacy_lengths;
    size_t legacy_size = num_sequences * max_length * sizeof(char);
    
    cudaMalloc(&d_legacy_sequences, legacy_size);
    cudaMalloc(&d_legacy_lengths, num_sequences * sizeof(int));
    
    // Convert format
    int threads = 256;
    int blocks = (num_sequences + threads - 1) / threads;
    
    convertSequenceFormat<<<blocks, threads, 0, stream>>>(
        buffer.getSequences(),
        buffer.getOffsets(),
        buffer.getLengths(),
        num_sequences,
        d_legacy_sequences,
        d_legacy_lengths,
        max_length
    );
    
    // Call existing mutation detection kernel
    // This kernel identifies specific mutations that confer antibiotic resistance
    detect_resistance_mutations<<<blocks, threads, 0, stream>>>(
        d_legacy_sequences,
        d_legacy_lengths,
        num_sequences,
        known_mutations_db.d_mutations,
        known_mutations_db.num_mutations,
        detected_mutations.d_hits,
        detected_mutations.d_hit_counts,
        detected_mutations.d_confidence_scores
    );
    
    // Free temporary buffers
    cudaFree(d_legacy_sequences);
    cudaFree(d_legacy_lengths);
}

// K-mer based mutation screening adapter for rapid diagnostic screening
void launchKmerDiagnosticKernel_v2(
    GPUSequenceBuffer& buffer,
    const uint64_t* d_mutation_kmers,
    int num_mutation_kmers,
    int kmer_length,
    bool* d_positive_screens,
    cudaStream_t stream) {
    
    int num_sequences = buffer.getCurrentSequences();
    
    // Launch k-mer screening kernel
    // This performs rapid screening for known resistance-conferring k-mers
    int threads = 128;
    int blocks = (num_sequences + threads - 1) / threads;
    
    screen_mutation_kmers<<<blocks, threads, 0, stream>>>(
        buffer.getSequences(),
        buffer.getOffsets(),
        buffer.getLengths(),
        num_sequences,
        d_mutation_kmers,
        num_mutation_kmers,
        kmer_length,
        d_positive_screens
    );
}

// Translated protein mutation search adapter for accurate amino acid variant detection
void launchProteinMutationSearch_v2(
    GPUSequenceBuffer& buffer,
    const ProteinMutationDatabase& protein_db,
    ProteinDiagnosticResults& results,
    cudaStream_t stream) {
    
    int num_sequences = buffer.getCurrentSequences();
    
    // Shared memory size for codon translation
    size_t shmem_size = 64 * sizeof(char);  // Genetic code table
    
    // Launch protein translation and mutation search
    // This identifies amino acid changes that affect antibiotic binding sites
    int threads = 64;
    int blocks = (num_sequences + threads - 1) / threads;
    
    translate_and_detect_mutations<<<blocks, threads, shmem_size, stream>>>(
        buffer.getSequences(),
        buffer.getOffsets(),
        buffer.getLengths(),
        num_sequences,
        protein_db.d_reference_proteins,
        protein_db.d_protein_offsets,
        protein_db.d_mutation_positions,
        protein_db.num_proteins,
        results.d_amino_acid_variants,
        results.d_variant_counts,
        results.d_protein_coverages
    );
}

// QRDR (Quinolone Resistance Determining Region) specific adapter
void launchQRDRAnalysis_v2(
    GPUSequenceBuffer& buffer,
    const QRDRDatabase& qrdr_db,
    QRDRDiagnosticResults& results,
    cudaStream_t stream) {
    
    int num_sequences = buffer.getCurrentSequences();
    
    // QRDR regions are critical for fluoroquinolone resistance
    // This kernel specifically checks mutations in gyrA, gyrB, parC, parE
    int threads = 128;
    int blocks = (num_sequences + threads - 1) / threads;
    
    analyze_qrdr_regions<<<blocks, threads, 0, stream>>>(
        buffer.getSequences(),
        buffer.getOffsets(),
        buffer.getLengths(),
        num_sequences,
        qrdr_db.d_qrdr_regions,
        qrdr_db.d_critical_positions,
        qrdr_db.num_regions,
        results.d_qrdr_mutations,
        results.d_mutation_counts,
        results.d_resistance_predictions
    );
}

// Batch processing adapter for clinical laboratory workflow
void launchBatchDiagnosticAnalysis_v2(
    GPUSequenceBuffer& buffer,
    const ClinicalMutationPanel& panel,
    BatchDiagnosticResults& results,
    cudaStream_t stream) {
    
    // This adapter processes multiple diagnostic targets in a single pass
    // Useful for comprehensive antibiotic resistance panels
    
    // Stage 1: Rapid k-mer screening
    launchKmerDiagnosticKernel_v2(
        buffer,
        panel.d_screening_kmers,
        panel.num_screening_kmers,
        panel.kmer_length,
        results.d_initial_screens,
        stream
    );
    
    // Stage 2: Targeted mutation detection on positive screens
    launchMutationDetectionKernel_v2(
        buffer,
        panel.mutation_db,
        results.mutation_results,
        stream
    );
    
    // Stage 3: Protein-level analysis for confirmed mutations
    launchProteinMutationSearch_v2(
        buffer,
        panel.protein_db,
        results.protein_results,
        stream
    );
    
    // Stage 4: QRDR analysis for fluoroquinolone resistance
    launchQRDRAnalysis_v2(
        buffer,
        panel.qrdr_db,
        results.qrdr_results,
        stream
    );
}

// Allele frequency calculation adapter for population-level resistance monitoring
__global__ void calculateAlleleFrequencies(
    const ProteinDiagnosticResults results,
    int num_sequences,
    AlleleFrequencyData* d_allele_freqs,
    int num_positions) {
    
    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    if (pos >= num_positions) return;
    
    // Count amino acid variants at each position
    int aa_counts[20] = {0};  // 20 standard amino acids
    int total_coverage = 0;
    
    for (int seq = 0; seq < num_sequences; seq++) {
        int variant_idx = seq * num_positions + pos;
        char aa = results.d_amino_acid_variants[variant_idx];
        
        if (aa >= 'A' && aa <= 'Y') {
            // Map amino acid to index
            int aa_idx = aa_to_index(aa);
            if (aa_idx >= 0 && aa_idx < 20) {
                aa_counts[aa_idx]++;
                total_coverage++;
            }
        }
    }
    
    // Calculate frequencies
    if (total_coverage > 0) {
        for (int i = 0; i < 20; i++) {
            d_allele_freqs[pos].frequencies[i] = 
                (float)aa_counts[i] / total_coverage;
        }
        d_allele_freqs[pos].total_coverage = total_coverage;
    }
}

// Clinical report generation adapter
void generateClinicalDiagnosticReport_v2(
    const BatchDiagnosticResults& results,
    const PatientMetadata& metadata,
    ClinicalReport& report,
    cudaStream_t stream) {
    
    // Aggregate results for clinical interpretation
    // This helps physicians make informed antibiotic treatment decisions
    
    // Copy results to host for report generation
    int num_mutations = 0;
    cudaMemcpyAsync(&num_mutations, results.mutation_results.d_hit_counts, 
                    sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    if (num_mutations > 0) {
        // Allocate host memory for mutations
        DiagnosticHit* h_mutations = new DiagnosticHit[num_mutations];
        cudaMemcpyAsync(h_mutations, results.mutation_results.d_hits,
                        num_mutations * sizeof(DiagnosticHit),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
        
        // Generate clinical interpretation
        report.patient_id = metadata.patient_id;
        report.sample_date = metadata.sample_date;
        report.sample_type = metadata.sample_type;
        
        // Classify resistance profile
        for (int i = 0; i < num_mutations; i++) {
            const auto& mutation = h_mutations[i];
            
            // Determine clinical significance
            if (mutation.confidence > 0.95f) {
                ResistanceInterpretation interp;
                interp.gene = mutation.gene_name;
                interp.mutation = mutation.mutation_name;
                interp.drug_class = getDrugClass(mutation.gene_name);
                interp.resistance_level = getResistanceLevel(mutation);
                interp.treatment_recommendation = 
                    getTreatmentRecommendation(interp.drug_class, interp.resistance_level);
                
                report.interpretations.push_back(interp);
            }
        }
        
        delete[] h_mutations;
    }
    
    // Generate summary recommendation
    report.overall_recommendation = generateTreatmentSummary(report.interpretations);
    report.confidence_level = calculateOverallConfidence(results);
}

// Helper functions for clinical interpretation
const char* getDrugClass(const std::string& gene) {
    if (gene == "gyrA" || gene == "gyrB" || gene == "parC" || gene == "parE") {
        return "Fluoroquinolone";
    } else if (gene == "rpoB") {
        return "Rifamycin";
    } else if (gene == "katG" || gene == "inhA") {
        return "Isoniazid";
    }
    // Add more gene-drug mappings as needed
    return "Unknown";
}

ResistanceLevel getResistanceLevel(const DiagnosticHit& mutation) {
    // Determine resistance level based on mutation type and position
    if (mutation.is_primary_mutation && mutation.in_critical_region) {
        return ResistanceLevel::HIGH;
    } else if (mutation.is_primary_mutation) {
        return ResistanceLevel::MODERATE;
    }
    return ResistanceLevel::LOW;
}

std::string getTreatmentRecommendation(const std::string& drug_class, 
                                       ResistanceLevel level) {
    if (level == ResistanceLevel::HIGH) {
        return "Avoid " + drug_class + "; use alternative antibiotic class";
    } else if (level == ResistanceLevel::MODERATE) {
        return "Use " + drug_class + " with caution; consider higher dose or combination therapy";
    }
    return drug_class + " likely effective at standard dose";
}

} // namespace BioGPU