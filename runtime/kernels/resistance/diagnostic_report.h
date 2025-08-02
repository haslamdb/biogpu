#ifndef DIAGNOSTIC_REPORT_H
#define DIAGNOSTIC_REPORT_H

#include <string>
#include <vector>
#include <cstdint>

// Forward declaration of ProteinMatch if it's defined elsewhere
// Or include the header that defines it
// Assuming it's defined in a file like "protein_match.h"
// #include "protein_match.h" 

// If not, define a minimal version for the interface
#ifndef PROTEIN_MATCH_DEFINED
#define PROTEIN_MATCH_DEFINED
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    uint8_t num_mutations;
    uint8_t mutation_positions[10];
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
    bool used_smith_waterman;
    char query_peptide[51];
    bool is_qrdr_alignment;
};
#endif

class DiagnosticReportGenerator {
public:
    DiagnosticReportGenerator(const std::string& output_path);
    ~DiagnosticReportGenerator();

    void addPipelineStatistics(int total_reads, int reads_after_bloom, int reads_after_kmer,
                               int total_protein_alignments, int total_mutations,
                               int qrdr_mutations, int sw_alignments);

    void addAlignments(const std::vector<ProteinMatch>& alignments);
    void generateReport();

private:
    void writeHeader();
    void writePipelineSummary();
    void writeQRDRMutationAnalysis();
    void writeAlignmentDetails();
    void writeFooter();
    void loadMetadataMappings();

    std::ofstream report_file;
    std::string report_path;

    struct PipelineStats {
        int total_reads = 0;
        int reads_after_bloom = 0;
        int reads_after_kmer = 0;
        int total_protein_alignments = 0;
        int total_mutations_found = 0;
        int qrdr_mutations_found = 0;
        int smith_waterman_alignments = 0;
        double bloom_retention_rate = 0.0;
        double kmer_retention_rate = 0.0;
        double protein_hit_rate = 0.0;
    } stats;

    std::vector<ProteinMatch> top_alignments;
    std::vector<ProteinMatch> qrdr_alignments;
};

#endif // DIAGNOSTIC_REPORT_H
