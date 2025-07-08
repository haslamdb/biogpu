// runtime/kernels/resistance/gpu_diagnostic_adapters.h
// Header for GPU diagnostic adapters - enabling clinical laboratories to identify
// bacterial mutations that affect antibiotic treatment decisions

#ifndef GPU_DIAGNOSTIC_ADAPTERS_H
#define GPU_DIAGNOSTIC_ADAPTERS_H

#include "../../common/gpu/gpu_sequence_buffer.h"
#include "fq_mutations_hardcoded.h"
#include <cuda_runtime.h>
#include <vector>
#include <string>

namespace BioGPU {

// Forward declarations
class GPUSequenceBuffer;

// Diagnostic result structures for clinical reporting
struct DiagnosticHit {
    uint32_t read_id;
    uint32_t mutation_id;
    char gene_name[32];
    char mutation_name[64];  // e.g., "GyrA_S83L"
    int position;
    char wildtype_aa;
    char mutant_aa;
    float confidence;        // Detection confidence (0-1)
    bool is_primary_mutation;
    bool in_critical_region;
    int supporting_reads;
};

struct DiagnosticResults {
    DiagnosticHit* d_hits;
    uint32_t* d_hit_counts;
    float* d_confidence_scores;
    size_t max_hits;
    
    DiagnosticResults(size_t max_hits_per_batch) : max_hits(max_hits_per_batch) {
        cudaMalloc(&d_hits, max_hits * sizeof(DiagnosticHit));
        cudaMalloc(&d_hit_counts, sizeof(uint32_t));
        cudaMalloc(&d_confidence_scores, max_hits * sizeof(float));
        cudaMemset(d_hit_counts, 0, sizeof(uint32_t));
    }
    
    ~DiagnosticResults() {
        cudaFree(d_hits);
        cudaFree(d_hit_counts);
        cudaFree(d_confidence_scores);
    }
};

// FQ Mutation Database for GPU
struct FQMutationDatabase {
    FQMutation* d_mutations;
    int num_mutations;
    
    FQMutationDatabase() {
        // Load hardcoded mutations
        const FQMutation* h_mutations = get_fq_mutations();
        num_mutations = get_num_fq_mutations();
        
        // Copy to GPU
        cudaMalloc(&d_mutations, num_mutations * sizeof(FQMutation));
        cudaMemcpy(d_mutations, h_mutations, num_mutations * sizeof(FQMutation), 
                   cudaMemcpyHostToDevice);
    }
    
    ~FQMutationDatabase() {
        cudaFree(d_mutations);
    }
};

// Protein mutation database for amino acid variant detection
struct ProteinMutationDatabase {
    char* d_reference_proteins;
    int* d_protein_offsets;
    int* d_mutation_positions;
    int num_proteins;
};

struct ProteinDiagnosticResults {
    char* d_amino_acid_variants;
    uint32_t* d_variant_counts;
    float* d_protein_coverages;
};

// QRDR specific structures
struct QRDRDatabase {
    struct QRDRRegion {
        char gene[32];
        char species[64];
        int start_pos;
        int end_pos;
        int critical_positions[10];
        int num_critical;
    };
    
    QRDRRegion* d_qrdr_regions;
    int* d_critical_positions;
    int num_regions;
};

struct QRDRDiagnosticResults {
    DiagnosticHit* d_qrdr_mutations;
    uint32_t* d_mutation_counts;
    float* d_resistance_predictions;  // Predicted resistance level (0-1)
};

// Clinical mutation panel for comprehensive testing
struct ClinicalMutationPanel {
    FQMutationDatabase mutation_db;
    ProteinMutationDatabase protein_db;
    QRDRDatabase qrdr_db;
    uint64_t* d_screening_kmers;
    int num_screening_kmers;
    int kmer_length;
};

// Batch diagnostic results
struct BatchDiagnosticResults {
    bool* d_initial_screens;
    DiagnosticResults mutation_results;
    ProteinDiagnosticResults protein_results;
    QRDRDiagnosticResults qrdr_results;
    
    BatchDiagnosticResults(size_t batch_size, size_t max_hits) 
        : mutation_results(max_hits) {
        cudaMalloc(&d_initial_screens, batch_size * sizeof(bool));
    }
    
    ~BatchDiagnosticResults() {
        cudaFree(d_initial_screens);
    }
};

// Allele frequency data for population-level monitoring
struct AlleleFrequencyData {
    float frequencies[20];  // Frequencies for 20 standard amino acids
    int total_coverage;
    int position;
    char reference_aa;
};

// Patient metadata for clinical context
struct PatientMetadata {
    std::string patient_id;
    std::string sample_date;
    std::string sample_type;  // e.g., "urine", "blood", "sputum"
    std::string suspected_organism;
    std::vector<std::string> current_antibiotics;
    std::vector<std::string> allergy_history;
};

// Resistance levels for clinical interpretation
enum class ResistanceLevel {
    NONE,
    LOW,
    MODERATE,
    HIGH,
    COMPLETE
};

// Clinical interpretation structure
struct ResistanceInterpretation {
    std::string gene;
    std::string mutation;
    std::string drug_class;
    ResistanceLevel resistance_level;
    std::string treatment_recommendation;
};

// Clinical report structure
struct ClinicalReport {
    std::string patient_id;
    std::string sample_date;
    std::string sample_type;
    std::vector<ResistanceInterpretation> interpretations;
    std::string overall_recommendation;
    float confidence_level;
    
    void generateHTML(const std::string& output_path) const;
    void generateJSON(const std::string& output_path) const;
    void generatePDF(const std::string& output_path) const;
};

// Main adapter functions for clinical diagnostic workflows

// Primary mutation detection for antibiotic resistance
void launchMutationDetectionKernel_v2(
    GPUSequenceBuffer& buffer,
    const FQMutationDatabase& known_mutations_db,
    DiagnosticResults& detected_mutations,
    cudaStream_t stream = 0
);

// Rapid k-mer based screening for known resistance markers
void launchKmerDiagnosticKernel_v2(
    GPUSequenceBuffer& buffer,
    const uint64_t* d_mutation_kmers,
    int num_mutation_kmers,
    int kmer_length,
    bool* d_positive_screens,
    cudaStream_t stream = 0
);

// Protein-level mutation analysis for amino acid changes
void launchProteinMutationSearch_v2(
    GPUSequenceBuffer& buffer,
    const ProteinMutationDatabase& protein_db,
    ProteinDiagnosticResults& results,
    cudaStream_t stream = 0
);

// QRDR-specific analysis for fluoroquinolone resistance
void launchQRDRAnalysis_v2(
    GPUSequenceBuffer& buffer,
    const QRDRDatabase& qrdr_db,
    QRDRDiagnosticResults& results,
    cudaStream_t stream = 0
);

// Comprehensive batch analysis for clinical laboratories
void launchBatchDiagnosticAnalysis_v2(
    GPUSequenceBuffer& buffer,
    const ClinicalMutationPanel& panel,
    BatchDiagnosticResults& results,
    cudaStream_t stream = 0
);

// Clinical report generation
void generateClinicalDiagnosticReport_v2(
    const BatchDiagnosticResults& results,
    const PatientMetadata& metadata,
    ClinicalReport& report,
    cudaStream_t stream = 0
);

// Helper functions for clinical interpretation
const char* getDrugClass(const std::string& gene);
ResistanceLevel getResistanceLevel(const DiagnosticHit& mutation);
std::string getTreatmentRecommendation(const std::string& drug_class, 
                                       ResistanceLevel level);
std::string generateTreatmentSummary(const std::vector<ResistanceInterpretation>& interps);
float calculateOverallConfidence(const BatchDiagnosticResults& results);

// Utility function to convert amino acid to array index
inline __device__ int aa_to_index(char aa) {
    // Convert single letter amino acid code to index (0-19)
    const char aa_order[] = "ACDEFGHIKLMNPQRSTVWY";
    for (int i = 0; i < 20; i++) {
        if (aa_order[i] == aa) return i;
    }
    return -1;
}

} // namespace BioGPU

#endif // GPU_DIAGNOSTIC_ADAPTERS_H