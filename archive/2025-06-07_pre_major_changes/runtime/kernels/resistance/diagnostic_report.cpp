// diagnostic_report.cpp
// Enhanced diagnostic report generator for BioGPU pipeline troubleshooting
// Creates detailed alignment visualization and mutation analysis

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <chrono>

// Include the ProteinMatch structure (should match translated_search.cu)
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
    char query_peptide[51];  // Store aligned peptide sequence (up to 50 AA + null terminator)
    bool is_qrdr_alignment;  // Flag for QRDR region alignment
};

// Dynamic gene and species name mappings (loaded from protein database metadata)
std::map<uint32_t, std::string> gene_names;
std::map<uint32_t, std::string> species_names;

// Known QRDR positions for each gene
std::map<uint32_t, std::vector<int>> qrdr_positions = {
    {0, {83, 87}},        // gyrA
    {1, {80, 84}},        // parC  
    {2, {426, 447}},      // gyrB
    {3, {416, 420}}       // parE
};

// Wild-type amino acids at key positions (for mutation detection)
std::map<std::pair<uint32_t, int>, char> wildtype_residues = {
    {{0, 83}, 'S'},   // gyrA S83
    {{0, 87}, 'D'},   // gyrA D87
    {{1, 80}, 'S'},   // parC S80
    {{1, 84}, 'E'},   // parC E84
    {{2, 426}, 'G'},  // gyrB G426 (example)
    {{2, 447}, 'A'},  // gyrB A447 (example)
    {{3, 416}, 'G'},  // parE G416 (example)
    {{3, 420}, 'A'}   // parE A420 (example)
};

class DiagnosticReportGenerator {
private:
    std::ofstream report_file;
    std::string report_path;
    
    // Pipeline statistics
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
    
public:
    DiagnosticReportGenerator(const std::string& output_path) : report_path(output_path) {
        report_file.open(output_path);
        if (!report_file.good()) {
            std::cerr << "ERROR: Failed to create diagnostic report: " << output_path << std::endl;
            return;
        }
        
        // Load dynamic gene and species mappings from protein database metadata
        loadMetadataMappings();
        
        writeHeader();
    }
    
    ~DiagnosticReportGenerator() {
        if (report_file.is_open()) {
            writeFooter();
            report_file.close();
        }
    }
    
    void updatePipelineStats(int total_reads, int after_bloom, int after_kmer, 
                           int protein_alignments, int mutations_found) {
        stats.total_reads = total_reads;
        stats.reads_after_bloom = after_bloom;
        stats.reads_after_kmer = after_kmer;
        stats.total_protein_alignments = protein_alignments;
        stats.total_mutations_found = mutations_found;
        
        // Calculate rates
        if (total_reads > 0) {
            stats.bloom_retention_rate = (double)after_bloom / total_reads * 100.0;
            stats.kmer_retention_rate = (double)after_kmer / total_reads * 100.0;
            stats.protein_hit_rate = (double)protein_alignments / total_reads * 100.0;
        }
    }
    
    void addProteinMatch(const ProteinMatch& match) {
        top_alignments.push_back(match);
        
        // Count Smith-Waterman usage
        if (match.used_smith_waterman) {
            stats.smith_waterman_alignments++;
        }
        
        // Check if this is a QRDR alignment
        uint32_t gene_id = match.gene_id;
        bool is_qrdr = false;
        
        if (qrdr_positions.find(gene_id) != qrdr_positions.end()) {
            auto& qrdr_pos = qrdr_positions[gene_id];
            // Check if the alignment overlaps with any QRDR position
            for (int pos : qrdr_pos) {
                if (pos >= match.ref_start && pos < match.ref_start + match.match_length) {
                    is_qrdr = true;
                    break;
                }
            }
        }
        
        // Add to QRDR alignments if it overlaps QRDR regions
        if (is_qrdr) {
            qrdr_alignments.push_back(match);
        }
        
        // Count QRDR mutations
        if (qrdr_positions.find(gene_id) != qrdr_positions.end()) {
            for (int i = 0; i < match.num_mutations; i++) {
                int pos = match.mutation_positions[i];
                auto& qrdr_pos = qrdr_positions[gene_id];
                if (std::find(qrdr_pos.begin(), qrdr_pos.end(), pos) != qrdr_pos.end()) {
                    stats.qrdr_mutations_found++;
                }
            }
        }
    }
    
    void generateReport() {
        writePipelineStatistics();
        writeTopAlignments();
        writeMutationAnalysis();
        writeQRDRAnalysis();
        writeTroubleshootingAdvice();
    }
    
private:
    void loadMetadataMappings() {
        // Load gene and species mappings from protein database metadata
        std::ifstream metadata_file("data/wildtype_protein_db/metadata.json");
        if (!metadata_file.good()) {
            std::cerr << "WARNING: Could not load protein database metadata for diagnostic reporting" << std::endl;
            // Fallback mappings
            gene_names[0] = "gyrA";
            gene_names[1] = "parC"; 
            gene_names[2] = "parE";
            species_names[0] = "Species0";
            species_names[1] = "Species1";
            return;
        }
        
        std::string line;
        bool in_species_map = false;
        bool in_gene_map = false;
        
        while (std::getline(metadata_file, line)) {
            // Look for "species_map" section
            if (line.find("\"species_map\"") != std::string::npos) {
                in_species_map = true;
                in_gene_map = false;
                continue;
            }
            
            // Look for "gene_map" section
            if (line.find("\"gene_map\"") != std::string::npos) {
                in_gene_map = true;
                in_species_map = false;
                continue;
            }
            
            // End of current section
            if ((in_species_map || in_gene_map) && line.find("}") != std::string::npos) {
                in_species_map = false;
                in_gene_map = false;
                continue;
            }
            
            // Parse mappings: "0": "name",
            if ((in_species_map || in_gene_map) && line.find("\":") != std::string::npos) {
                size_t id_start = line.find("\"") + 1;
                size_t id_end = line.find("\"", id_start);
                size_t name_start = line.find("\"", id_end + 1) + 1;
                size_t name_end = line.find("\"", name_start);
                
                if (id_end != std::string::npos && name_end != std::string::npos) {
                    uint32_t id = std::stoi(line.substr(id_start, id_end - id_start));
                    std::string name = line.substr(name_start, name_end - name_start);
                    
                    if (in_species_map) {
                        species_names[id] = name;
                    } else if (in_gene_map) {
                        gene_names[id] = name;
                    }
                }
            }
        }
        
        metadata_file.close();
        
        // Ensure we have some fallback mappings
        if (gene_names.empty()) {
            gene_names[0] = "gyrA";
            gene_names[1] = "parC";
            gene_names[2] = "parE";
        }
        
        if (species_names.empty()) {
            species_names[0] = "Species0";
            species_names[1] = "Species1";
        }
    }
    
    void writeHeader() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        
        report_file << "================================\n";
        report_file << "BioGPU PIPELINE DIAGNOSTIC REPORT\n";
        report_file << "================================\n\n";
        report_file << "Generated: " << std::ctime(&time_t);
        report_file << "Pipeline Version: 0.4.0 (Enhanced Translated Search)\n";
        report_file << "Features: Bloom Filter + 15-mer nucleotide + 5-mer protein + Smith-Waterman\n\n";
    }
    
    void writePipelineStatistics() {
        report_file << "PIPELINE STATISTICS\n";
        report_file << "===================\n\n";
        
        report_file << "Read Processing:\n";
        report_file << "  Total input reads:           " << std::setw(10) << stats.total_reads << "\n";
        report_file << "  After Bloom filter:          " << std::setw(10) << stats.reads_after_bloom 
                   << " (" << std::fixed << std::setprecision(1) << stats.bloom_retention_rate << "%)\n";
        report_file << "  After k-mer enrichment:      " << std::setw(10) << stats.reads_after_kmer
                   << " (" << std::fixed << std::setprecision(1) << stats.kmer_retention_rate << "%)\n\n";
        
        report_file << "Protein Search Results:\n";
        report_file << "  Total protein alignments:    " << std::setw(10) << stats.total_protein_alignments << "\n";
        report_file << "  Protein hit rate:             " << std::setw(9) << std::fixed << std::setprecision(1) 
                   << stats.protein_hit_rate << "%\n";
        report_file << "  Smith-Waterman alignments:   " << std::setw(10) << stats.smith_waterman_alignments << "\n\n";
        
        report_file << "Mutation Detection:\n";
        report_file << "  Total mutations found:       " << std::setw(10) << stats.total_mutations_found << "\n";
        report_file << "  QRDR mutations found:         " << std::setw(10) << stats.qrdr_mutations_found << "\n\n";
        
        // Filter efficiency analysis
        int bloom_filtered = stats.total_reads - stats.reads_after_bloom;
        int kmer_filtered = stats.reads_after_bloom - stats.reads_after_kmer;
        
        report_file << "Filter Efficiency:\n";
        report_file << "  Bloom filter removed:         " << std::setw(10) << bloom_filtered 
                   << " (" << std::fixed << std::setprecision(1) 
                   << (double)bloom_filtered / stats.total_reads * 100.0 << "%)\n";
        report_file << "  K-mer filter removed:         " << std::setw(10) << kmer_filtered
                   << " (" << std::fixed << std::setprecision(1)
                   << (double)kmer_filtered / stats.reads_after_bloom * 100.0 << "%)\n\n";
    }
    
    void writeTopAlignments() {
        report_file << "TOP 20 QRDR ALIGNMENTS\n";
        report_file << "======================\n\n";
        
        if (qrdr_alignments.empty()) {
            report_file << "No QRDR alignments found.\n\n";
            return;
        }
        
        // Sort by alignment score (descending)
        std::sort(qrdr_alignments.begin(), qrdr_alignments.end(),
                 [](const ProteinMatch& a, const ProteinMatch& b) {
                     return a.alignment_score > b.alignment_score;
                 });
        
        report_file << "These represent alignments within Quinolone Resistance Determining Regions:\n\n";
        
        // Show top 20 QRDR alignments (same format as stdout)
        int num_to_show = std::min(20, (int)qrdr_alignments.size());
        
        for (int i = 0; i < num_to_show; i++) {
            const auto& match = qrdr_alignments[i];
            
            // Map gene IDs to names using dynamic mapping from metadata
            std::string gene_name;
            if (gene_names.count(match.gene_id)) {
                gene_name = gene_names[match.gene_id];
            } else {
                gene_name = "Gene" + std::to_string(match.gene_id);
            }
            
            // Map species IDs to names using dynamic mapping from metadata  
            std::string species_name;
            if (species_names.count(match.species_id)) {
                species_name = species_names[match.species_id];
            } else {
                species_name = "Species" + std::to_string(match.species_id);
            }
            
            report_file << "  " << (i + 1) << ". Read " << match.read_id 
                       << ", " << species_name << " " << gene_name 
                       << ", Frame " << (int)match.frame
                       << ", Score: " << std::fixed << std::setprecision(1) << match.alignment_score
                       << ", Identity: " << std::fixed << std::setprecision(1) << (match.identity * 100) << "%"
                       << ", Position: " << match.ref_start << "-" << (match.ref_start + match.match_length - 1)
                       << ", Length: " << match.match_length
                       << ", Peptide: " << match.query_peptide << "\n";
        }
        
        report_file << "\nTotal QRDR alignments found: " << qrdr_alignments.size() << "\n\n";
    }
    
    void writeAlignmentVisualization(const ProteinMatch& match, int rank) {
        std::string gene_name = gene_names.count(match.gene_id) ? 
                               gene_names[match.gene_id] : "Unknown";
        
        report_file << "Alignment #" << rank << " - " << gene_name 
                   << " (Score: " << std::fixed << std::setprecision(1) << match.alignment_score
                   << ", Identity: " << std::fixed << std::setprecision(1) << match.identity * 100 << "%)\n";
        
        report_file << "Read " << match.read_id << ", Frame " << (int)match.frame
                   << ", Protein " << match.protein_id << "\n";
        
        if (match.used_smith_waterman) {
            report_file << "*** SMITH-WATERMAN ALIGNMENT ***\n";
        }
        
        // Create alignment visualization
        std::string query_seq = "(Query sequence not available in current structure)";
        std::string ref_seq = "(Reference sequence not available in current structure)";
        
        report_file << "Position " << match.ref_start << "-" << (match.ref_start + match.match_length - 1) << ":\n";
        
        // Show mutations if any
        if (match.num_mutations > 0) {
            report_file << "MUTATIONS DETECTED (" << (int)match.num_mutations << "):\n";
            for (int i = 0; i < match.num_mutations; i++) {
                int pos = match.mutation_positions[i];
                char ref_aa = match.ref_aas[i];
                char query_aa = match.query_aas[i];
                float blosum = match.blosum_scores[i];
                
                report_file << "  Position " << (match.ref_start + pos) << ": " 
                           << ref_aa << " -> " << query_aa 
                           << " (BLOSUM: " << std::fixed << std::setprecision(1) << blosum << ")";
                
                // Check if this is a known QRDR position
                if (qrdr_positions.count(match.gene_id)) {
                    auto& qrdr_pos = qrdr_positions[match.gene_id];
                    int global_pos = match.ref_start + pos;
                    if (std::find(qrdr_pos.begin(), qrdr_pos.end(), global_pos) != qrdr_pos.end()) {
                        report_file << " *** QRDR MUTATION ***";
                    }
                }
                report_file << "\n";
            }
        } else {
            report_file << "No mutations detected in this alignment.\n";
        }
        
        // Simplified alignment view (since we don't have full sequences in the structure)
        report_file << "\nAlignment region:\n";
        report_file << "Query:  [" << match.query_start << "-" << (match.query_start + match.match_length - 1) << "] ";
        if (match.num_mutations > 0) {
            for (int i = 0; i < match.match_length && i < 20; i++) {
                bool is_mutation = false;
                char display_char = 'X';  // Placeholder since we don't have sequence
                
                for (int j = 0; j < match.num_mutations; j++) {
                    if (match.mutation_positions[j] == i) {
                        display_char = match.query_aas[j];
                        is_mutation = true;
                        break;
                    }
                }
                
                if (is_mutation) {
                    report_file << "\033[31m" << display_char << "\033[0m";  // Red for mutations
                } else {
                    report_file << display_char;
                }
            }
        } else {
            report_file << "(Perfect match - no sequence differences)";
        }
        report_file << "\n";
        
        report_file << "Ref:    [" << match.ref_start << "-" << (match.ref_start + match.match_length - 1) << "] ";
        if (match.num_mutations > 0) {
            for (int i = 0; i < match.match_length && i < 20; i++) {
                char display_char = 'X';  // Placeholder
                
                for (int j = 0; j < match.num_mutations; j++) {
                    if (match.mutation_positions[j] == i) {
                        display_char = match.ref_aas[j];
                        break;
                    }
                }
                report_file << display_char;
            }
        } else {
            report_file << "(Perfect match)";
        }
        report_file << "\n";
    }
    
    void writeMutationAnalysis() {
        report_file << "\nMUTATION ANALYSIS\n";
        report_file << "=================\n\n";
        
        // Count mutations by gene
        std::map<uint32_t, int> mutations_per_gene;
        std::map<std::pair<uint32_t, int>, int> position_mutations;
        
        for (const auto& match : top_alignments) {
            if (match.num_mutations > 0) {
                mutations_per_gene[match.gene_id] += match.num_mutations;
                
                for (int i = 0; i < match.num_mutations; i++) {
                    int global_pos = match.ref_start + match.mutation_positions[i];
                    position_mutations[{match.gene_id, global_pos}]++;
                }
            }
        }
        
        if (mutations_per_gene.empty()) {
            report_file << "❌ NO MUTATIONS DETECTED\n\n";
        } else {
            report_file << "Mutations by gene:\n";
            for (const auto& pair : mutations_per_gene) {
                std::string gene_name = gene_names.count(pair.first) ? 
                                      gene_names[pair.first] : "Unknown";
                report_file << "  " << gene_name << ": " << pair.second << " mutations\n";
            }
            report_file << "\n";
            
            report_file << "Mutations by position:\n";
            for (const auto& pair : position_mutations) {
                uint32_t gene_id = pair.first.first;
                int position = pair.first.second;
                int count = pair.second;
                
                std::string gene_name = gene_names.count(gene_id) ? 
                                      gene_names[gene_id] : "Unknown";
                
                report_file << "  " << gene_name << " position " << position << ": " 
                           << count << " occurrences";
                
                // Check if QRDR
                if (qrdr_positions.count(gene_id)) {
                    auto& qrdr_pos = qrdr_positions[gene_id];
                    if (std::find(qrdr_pos.begin(), qrdr_pos.end(), position) != qrdr_pos.end()) {
                        report_file << " (QRDR)";
                    }
                }
                report_file << "\n";
            }
        }
        report_file << "\n";
    }
    
    void writeQRDRAnalysis() {
        report_file << "QRDR RESISTANCE ANALYSIS\n";
        report_file << "========================\n\n";
        
        // Analyze coverage of QRDR regions
        std::map<uint32_t, bool> qrdr_covered;
        std::map<uint32_t, std::vector<int>> qrdr_mutations_found;
        
        for (const auto& match : top_alignments) {
            uint32_t gene_id = match.gene_id;
            if (qrdr_positions.count(gene_id)) {
                qrdr_covered[gene_id] = true;
                
                // Check for mutations in QRDR
                for (int i = 0; i < match.num_mutations; i++) {
                    int global_pos = match.ref_start + match.mutation_positions[i];
                    auto& qrdr_pos = qrdr_positions[gene_id];
                    if (std::find(qrdr_pos.begin(), qrdr_pos.end(), global_pos) != qrdr_pos.end()) {
                        qrdr_mutations_found[gene_id].push_back(global_pos);
                    }
                }
            }
        }
        
        report_file << "QRDR Coverage:\n";
        for (const auto& pair : qrdr_positions) {
            uint32_t gene_id = pair.first;
            std::string gene_name = gene_names.count(gene_id) ? 
                                  gene_names[gene_id] : "Unknown";
            
            bool covered = qrdr_covered.count(gene_id) && qrdr_covered[gene_id];
            report_file << "  " << gene_name << " QRDR: " 
                       << (covered ? "✓ COVERED" : "❌ NOT COVERED") << "\n";
            
            if (covered && qrdr_mutations_found.count(gene_id)) {
                auto& mutations = qrdr_mutations_found[gene_id];
                if (!mutations.empty()) {
                    report_file << "    Mutations found at positions: ";
                    for (size_t i = 0; i < mutations.size(); i++) {
                        if (i > 0) report_file << ", ";
                        report_file << mutations[i];
                    }
                    report_file << "\n";
                }
            }
        }
        
        report_file << "\nFluoroquinolone Resistance Assessment:\n";
        if (stats.qrdr_mutations_found > 0) {
            report_file << "✓ RESISTANCE MUTATIONS DETECTED (" << stats.qrdr_mutations_found << " total)\n";
            report_file << "  Recommendation: Avoid fluoroquinolone therapy\n";
        } else {
            report_file << "❌ NO RESISTANCE MUTATIONS DETECTED\n";
            report_file << "  Recommendation: Fluoroquinolones may be effective\n";
            report_file << "  Note: Check mutation detection logic if mutant reads were expected\n";
        }
        report_file << "\n";
    }
    
    void writeTroubleshootingAdvice() {
        report_file << "TROUBLESHOOTING ADVICE\n";
        report_file << "======================\n\n";
        
        if (stats.total_protein_alignments == 0) {
            report_file << "🔍 NO PROTEIN ALIGNMENTS FOUND\n";
            report_file << "Possible causes:\n";
            report_file << "- Protein database not loaded or corrupted\n";
            report_file << "- K-mer size mismatch between reads and database\n";
            report_file << "- Reads don't contain resistance genes\n";
            report_file << "- Identity threshold too stringent (current: 90%)\n\n";
        }
        
        
        if (stats.bloom_retention_rate < 0.1) {
            report_file << "⚠️  EXTREMELY LOW BLOOM FILTER RETENTION RATE (<0.1%)\n";
            report_file << "This may indicate reads don't contain any target sequences.\n";
            report_file << "This is expected for metagenomic data with low pathogen abundance.\n\n";
        }
        
        if (stats.protein_hit_rate < 0.1) {
            report_file << "⚠️  EXTREMELY LOW PROTEIN HIT RATE (<0.1%)\n";
            report_file << "This may indicate very low abundance of target organisms.\n";
            report_file << "This is expected for metagenomic data.\n\n";
        }
        
        if (stats.total_protein_alignments > 0) {
            report_file << "✅ Target sequences detected\n";
            report_file << "Review QRDR alignments above for resistance mutations.\n\n";
        }
    }
    
    void writeFooter() {
        report_file << "\nReport generated by BioGPU Enhanced Diagnostic System\n";
        report_file << "For technical support, check pipeline documentation\n";
        report_file << "================================\n";
    }
};

// C interface for integration with existing pipeline
extern "C" {
    void* create_diagnostic_reporter(const char* output_path) {
        return new DiagnosticReportGenerator(std::string(output_path));
    }
    
    void destroy_diagnostic_reporter(void* reporter) {
        if (reporter) {
            delete static_cast<DiagnosticReportGenerator*>(reporter);
        }
    }
    
    void update_pipeline_statistics(void* reporter, int total_reads, int after_bloom, 
                                  int after_kmer, int protein_alignments, int mutations_found) {
        if (reporter) {
            static_cast<DiagnosticReportGenerator*>(reporter)->updatePipelineStats(
                total_reads, after_bloom, after_kmer, protein_alignments, mutations_found);
        }
    }
    
    void add_protein_match_to_report(void* reporter, const ProteinMatch* match) {
        if (reporter && match) {
            static_cast<DiagnosticReportGenerator*>(reporter)->addProteinMatch(*match);
        }
    }
    
    void generate_diagnostic_report(void* reporter) {
        if (reporter) {
            static_cast<DiagnosticReportGenerator*>(reporter)->generateReport();
        }
    }
}