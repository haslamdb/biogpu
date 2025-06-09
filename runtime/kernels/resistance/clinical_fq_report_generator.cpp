// clinical_fq_report_generator.cpp
// Clinical-focused fluoroquinolone resistance diagnostic report generator
// Provides comprehensive analysis of FQ resistance for clinical decision-making

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <chrono>
#include <cmath>
#include "global_fq_resistance_mapper.h"

// Structure for clinical reporting
struct ClinicalFQReport {
    // Summary flags
    bool has_fq_resistance;
    bool has_qrdr_mutations;
    bool has_high_confidence_resistance;
    
    // Species-level summary
    struct SpeciesSummary {
        std::string species_name;
        int total_reads;
        int reads_with_mutations;
        int fq_resistance_mutations;
        int qrdr_mutations;
        std::map<std::string, std::vector<std::string>> gene_mutations; // gene -> list of mutations
        float max_identity_score;
        bool likely_resistant;
    };
    
    // Gene-level details
    struct GeneMutation {
        std::string species;
        std::string gene;
        std::string mutation_code;  // e.g., "S83L"
        int position;
        char wildtype_aa;
        char mutant_aa;
        int occurrence_count;
        float avg_identity;
        float avg_alignment_score;
        bool is_known_fq_resistance;
        bool is_in_qrdr;
        std::string clinical_significance;
    };
    
    // Collected data
    std::map<std::string, SpeciesSummary> species_summaries;
    std::vector<GeneMutation> all_mutations;
    std::map<std::string, int> mutation_counts;  // mutation_code -> count
    
    // Allele frequency data
    struct AlleleFrequencyEntry {
        std::string species;
        std::string gene;
        uint16_t position;
        uint32_t total_depth;
        char wildtype_aa;
        uint32_t wildtype_count;
        float wildtype_frequency;
        char dominant_mutant_aa;
        uint32_t dominant_mutant_count;
        float dominant_mutant_frequency;
        float total_resistant_frequency;
        bool has_resistance_mutation;
        std::string mutation_summary;
        std::map<char, uint32_t> aa_counts;
    };
    std::vector<AlleleFrequencyEntry> allele_frequencies;
    bool has_allele_frequency_data;
    
    // Statistics
    int total_reads_analyzed;
    int reads_with_protein_matches;
    int reads_with_any_mutations;
    int reads_with_fq_resistance;
    
    // Performance metrics
    double processing_time_seconds;
    double reads_per_second;
    
    // Clinical interpretation
    std::string overall_interpretation;
    std::vector<std::string> clinical_notes;
    float resistance_confidence;  // 0-1 scale
};

class ClinicalFQReportGenerator {
private:
    ClinicalFQReport report;
    std::string output_path;
    GlobalFQResistanceMapper& fq_mapper;
    
    // Known QRDR regions for major genes
    std::map<std::string, std::pair<int, int>> qrdr_regions = {
        {"gyrA", {67, 106}},   // E. coli numbering
        {"gyrB", {426, 464}},
        {"parC", {64, 102}},
        {"parE", {410, 450}},
        {"grlA", {64, 102}},   // S. aureus
        {"grlB", {410, 450}}
    };
    
    // Clinical significance mapping
    std::map<std::string, std::string> clinical_significance = {
        // High-level resistance mutations
        {"gyrA_S83L", "High-level FQ resistance"},
        {"gyrA_S83F", "High-level FQ resistance"},
        {"gyrA_D87N", "High-level FQ resistance"},
        {"gyrA_D87G", "High-level FQ resistance"},
        {"parC_S80I", "Moderate FQ resistance"},
        {"parC_E84K", "Moderate FQ resistance"},
        // Add more as needed
    };
    
public:
    ClinicalFQReportGenerator(const std::string& output) 
        : output_path(output), fq_mapper(GlobalFQResistanceMapper::getInstance()) {
        report.has_fq_resistance = false;
        report.has_qrdr_mutations = false;
        report.has_high_confidence_resistance = false;
        report.has_allele_frequency_data = false;
        report.total_reads_analyzed = 0;
        report.reads_with_protein_matches = 0;
        report.reads_with_any_mutations = 0;
        report.reads_with_fq_resistance = 0;
        report.resistance_confidence = 0.0f;
        report.processing_time_seconds = 0.0;
        report.reads_per_second = 0.0;
    }
    
    // Process a protein match from the pipeline
    void processProteinMatch(uint32_t read_id, const std::string& species_name, 
                           const std::string& gene_name, uint32_t gene_id, uint32_t species_id,
                           const std::vector<std::tuple<int, char, char>>& mutations,
                           float alignment_score, float identity, int match_length,
                           int ref_start, int ref_end) {
        
        if (mutations.empty()) return;
        
        // Update species summary
        auto& species_summary = report.species_summaries[species_name];
        species_summary.species_name = species_name;
        species_summary.total_reads++;
        species_summary.reads_with_mutations++;
        species_summary.max_identity_score = std::max(species_summary.max_identity_score, identity);
        
        bool has_fq_in_match = false;
        bool has_qrdr_in_match = false;
        
        // Process each mutation
        for (const auto& [pos, ref_aa, mut_aa] : mutations) {
            int global_position = ref_start + pos + 1;  // Convert to 1-based
            
            // Check if in QRDR region
            bool is_qrdr = false;
            auto qrdr_it = qrdr_regions.find(gene_name);
            if (qrdr_it != qrdr_regions.end()) {
                is_qrdr = (global_position >= qrdr_it->second.first && 
                          global_position <= qrdr_it->second.second);
            }
            
            // Check if known FQ resistance mutation
            bool is_fq_resistance = fq_mapper.isResistanceMutation(
                species_name, gene_name, global_position, ref_aa, mut_aa
            );
            
            if (is_qrdr) {
                has_qrdr_in_match = true;
                report.has_qrdr_mutations = true;
                species_summary.qrdr_mutations++;
            }
            
            if (is_fq_resistance) {
                has_fq_in_match = true;
                report.has_fq_resistance = true;
                species_summary.fq_resistance_mutations++;
            }
            
            // Create mutation record
            std::string mutation_code = std::string(1, ref_aa) + std::to_string(global_position) + mut_aa;
            species_summary.gene_mutations[gene_name].push_back(mutation_code);
            
            // Update mutation counts
            std::string full_mutation_key = species_name + "_" + gene_name + "_" + mutation_code;
            report.mutation_counts[full_mutation_key]++;
            
            // Add to detailed mutations list
            ClinicalFQReport::GeneMutation gene_mut;
            gene_mut.species = species_name;
            gene_mut.gene = gene_name;
            gene_mut.mutation_code = mutation_code;
            gene_mut.position = global_position;
            gene_mut.wildtype_aa = ref_aa;
            gene_mut.mutant_aa = mut_aa;
            gene_mut.occurrence_count = 1;  // Will be aggregated later
            gene_mut.avg_identity = identity;
            gene_mut.avg_alignment_score = alignment_score;
            gene_mut.is_known_fq_resistance = is_fq_resistance;
            gene_mut.is_in_qrdr = is_qrdr;
            
            // Assign clinical significance
            std::string sig_key = gene_name + "_" + mutation_code;
            auto sig_it = clinical_significance.find(sig_key);
            if (sig_it != clinical_significance.end()) {
                gene_mut.clinical_significance = sig_it->second;
            } else if (is_fq_resistance) {
                gene_mut.clinical_significance = "Known FQ resistance mutation";
            } else if (is_qrdr) {
                gene_mut.clinical_significance = "QRDR mutation of uncertain significance";
            } else {
                gene_mut.clinical_significance = "Non-QRDR mutation";
            }
            
            report.all_mutations.push_back(gene_mut);
        }
        
        if (has_fq_in_match) {
            report.reads_with_fq_resistance++;
            species_summary.likely_resistant = true;
        }
        
        report.reads_with_any_mutations++;
    }
    
    // Update statistics
    void updateReadStats(int total_reads, int reads_with_matches) {
        report.total_reads_analyzed = total_reads;
        report.reads_with_protein_matches = reads_with_matches;
    }
    
    // Update performance statistics
    void updatePerformanceStats(double processing_seconds, double reads_per_sec) {
        report.processing_time_seconds = processing_seconds;
        report.reads_per_second = reads_per_sec;
    }
    
    // Add allele frequency data
    void addAlleleFrequencyData(const std::string& species, const std::string& gene,
                                uint16_t position, uint32_t total_depth,
                                char wildtype_aa, uint32_t wildtype_count, float wildtype_frequency,
                                char dominant_mutant_aa, uint32_t dominant_mutant_count, 
                                float dominant_mutant_frequency, float total_resistant_frequency,
                                bool has_resistance_mutation, const std::string& mutation_summary,
                                const std::map<char, uint32_t>& aa_counts) {
        ClinicalFQReport::AlleleFrequencyEntry entry;
        entry.species = species;
        entry.gene = gene;
        entry.position = position;
        entry.total_depth = total_depth;
        entry.wildtype_aa = wildtype_aa;
        entry.wildtype_count = wildtype_count;
        entry.wildtype_frequency = wildtype_frequency;
        entry.dominant_mutant_aa = dominant_mutant_aa;
        entry.dominant_mutant_count = dominant_mutant_count;
        entry.dominant_mutant_frequency = dominant_mutant_frequency;
        entry.total_resistant_frequency = total_resistant_frequency;
        entry.has_resistance_mutation = has_resistance_mutation;
        entry.mutation_summary = mutation_summary;
        entry.aa_counts = aa_counts;
        
        report.allele_frequencies.push_back(entry);
        report.has_allele_frequency_data = true;
    }
    
    // Generate clinical interpretation
    void generateClinicalInterpretation() {
        // Calculate confidence based on mutation quality and quantity
        if (report.has_fq_resistance) {
            int high_confidence_mutations = 0;
            for (const auto& mut : report.all_mutations) {
                if (mut.is_known_fq_resistance && mut.avg_identity > 0.95) {
                    high_confidence_mutations++;
                }
            }
            
            if (high_confidence_mutations >= 2) {
                report.resistance_confidence = 0.95f;
                report.has_high_confidence_resistance = true;
                report.overall_interpretation = "HIGH CONFIDENCE: Fluoroquinolone resistance detected";
            } else if (high_confidence_mutations == 1) {
                report.resistance_confidence = 0.80f;
                report.overall_interpretation = "MODERATE CONFIDENCE: Fluoroquinolone resistance likely";
            } else {
                report.resistance_confidence = 0.60f;
                report.overall_interpretation = "LOW CONFIDENCE: Possible fluoroquinolone resistance";
            }
            
            // Add specific clinical notes
            std::set<std::string> resistant_species;
            for (const auto& [species, summary] : report.species_summaries) {
                if (summary.likely_resistant) {
                    resistant_species.insert(species);
                }
            }
            
            if (!resistant_species.empty()) {
                std::string species_list;
                for (const auto& sp : resistant_species) {
                    if (!species_list.empty()) species_list += ", ";
                    species_list += sp;
                }
                report.clinical_notes.push_back("Resistance detected in: " + species_list);
            }
            
            // Check for multi-drug resistance patterns
            std::set<std::string> affected_genes;
            for (const auto& mut : report.all_mutations) {
                if (mut.is_known_fq_resistance) {
                    affected_genes.insert(mut.gene);
                }
            }
            
            if (affected_genes.size() > 1) {
                report.clinical_notes.push_back("Multiple resistance genes affected - consider alternative antibiotics");
            }
            
        } else if (report.has_qrdr_mutations) {
            report.resistance_confidence = 0.40f;
            report.overall_interpretation = "QRDR mutations detected - resistance possible but not confirmed";
            report.clinical_notes.push_back("Consider phenotypic susceptibility testing");
        } else {
            report.resistance_confidence = 0.05f;
            report.overall_interpretation = "No fluoroquinolone resistance markers detected";
            report.clinical_notes.push_back("Sample appears susceptible to fluoroquinolones");
        }
    }
    
    // Generate the full clinical report
    void generateReport() {
        generateClinicalInterpretation();
        
        // Generate multiple output formats
        generateJSONReport();
        generateHTMLReport();
        generateTextReport();
    }
    
private:
    void generateJSONReport() {
        std::ofstream json_file(output_path + "_clinical_fq_report.json");
        json_file << std::setprecision(3) << std::fixed;
        
        json_file << "{\n";
        json_file << "  \"report_type\": \"clinical_fluoroquinolone_resistance\",\n";
        json_file << "  \"generated\": \"" << getCurrentTimestamp() << "\",\n";
        json_file << "  \"summary\": {\n";
        json_file << "    \"overall_interpretation\": \"" << report.overall_interpretation << "\",\n";
        json_file << "    \"resistance_confidence\": " << report.resistance_confidence << ",\n";
        json_file << "    \"has_fq_resistance\": " << (report.has_fq_resistance ? "true" : "false") << ",\n";
        json_file << "    \"has_qrdr_mutations\": " << (report.has_qrdr_mutations ? "true" : "false") << ",\n";
        json_file << "    \"total_reads_analyzed\": " << report.total_reads_analyzed << ",\n";
        json_file << "    \"reads_with_protein_matches\": " << report.reads_with_protein_matches << ",\n";
        json_file << "    \"reads_with_mutations\": " << report.reads_with_any_mutations << ",\n";
        json_file << "    \"reads_with_fq_resistance\": " << report.reads_with_fq_resistance << ",\n";
        json_file << "    \"processing_time_seconds\": " << report.processing_time_seconds << ",\n";
        json_file << "    \"reads_per_second\": " << std::fixed << std::setprecision(0) << report.reads_per_second << "\n";
        json_file << "  },\n";
        
        // Clinical notes
        json_file << "  \"clinical_notes\": [\n";
        for (size_t i = 0; i < report.clinical_notes.size(); i++) {
            json_file << "    \"" << report.clinical_notes[i] << "\"";
            if (i < report.clinical_notes.size() - 1) json_file << ",";
            json_file << "\n";
        }
        json_file << "  ],\n";
        
        // Species breakdown
        json_file << "  \"species_analysis\": [\n";
        bool first_species = true;
        for (const auto& [species, summary] : report.species_summaries) {
            if (!first_species) json_file << ",\n";
            json_file << "    {\n";
            json_file << "      \"species\": \"" << species << "\",\n";
            json_file << "      \"likely_resistant\": " << (summary.likely_resistant ? "true" : "false") << ",\n";
            json_file << "      \"fq_resistance_mutations\": " << summary.fq_resistance_mutations << ",\n";
            json_file << "      \"qrdr_mutations\": " << summary.qrdr_mutations << ",\n";
            json_file << "      \"genes_affected\": [";
            
            bool first_gene = true;
            for (const auto& [gene, muts] : summary.gene_mutations) {
                if (!first_gene) json_file << ", ";
                json_file << "\"" << gene << "\"";
                first_gene = false;
            }
            json_file << "]\n";
            json_file << "    }";
            first_species = false;
        }
        json_file << "\n  ],\n";
        
        // Detailed mutations
        json_file << "  \"mutations_detected\": [\n";
        
        // Aggregate and sort mutations by clinical importance
        std::map<std::string, ClinicalFQReport::GeneMutation> aggregated_mutations;
        for (const auto& mut : report.all_mutations) {
            std::string key = mut.species + "_" + mut.gene + "_" + mut.mutation_code;
            if (aggregated_mutations.find(key) == aggregated_mutations.end()) {
                aggregated_mutations[key] = mut;
            } else {
                aggregated_mutations[key].occurrence_count++;
                aggregated_mutations[key].avg_identity = 
                    (aggregated_mutations[key].avg_identity + mut.avg_identity) / 2;
                aggregated_mutations[key].avg_alignment_score = 
                    (aggregated_mutations[key].avg_alignment_score + mut.avg_alignment_score) / 2;
            }
        }
        
        // Sort by importance: FQ resistance > QRDR > other
        std::vector<std::pair<std::string, ClinicalFQReport::GeneMutation>> sorted_mutations(
            aggregated_mutations.begin(), aggregated_mutations.end()
        );
        std::sort(sorted_mutations.begin(), sorted_mutations.end(),
            [](const auto& a, const auto& b) {
                if (a.second.is_known_fq_resistance != b.second.is_known_fq_resistance)
                    return a.second.is_known_fq_resistance;
                if (a.second.is_in_qrdr != b.second.is_in_qrdr)
                    return a.second.is_in_qrdr;
                return a.second.occurrence_count > b.second.occurrence_count;
            }
        );
        
        bool first_mut = true;
        for (const auto& [key, mut] : sorted_mutations) {
            if (!first_mut) json_file << ",\n";
            json_file << "    {\n";
            json_file << "      \"species\": \"" << mut.species << "\",\n";
            json_file << "      \"gene\": \"" << mut.gene << "\",\n";
            json_file << "      \"mutation\": \"" << mut.mutation_code << "\",\n";
            json_file << "      \"position\": " << mut.position << ",\n";
            json_file << "      \"change\": \"" << mut.wildtype_aa << " → " << mut.mutant_aa << "\",\n";
            json_file << "      \"occurrences\": " << mut.occurrence_count << ",\n";
            json_file << "      \"avg_identity\": " << mut.avg_identity << ",\n";
            json_file << "      \"is_fq_resistance\": " << (mut.is_known_fq_resistance ? "true" : "false") << ",\n";
            json_file << "      \"is_qrdr\": " << (mut.is_in_qrdr ? "true" : "false") << ",\n";
            json_file << "      \"clinical_significance\": \"" << mut.clinical_significance << "\"\n";
            json_file << "    }";
            first_mut = false;
        }
        json_file << "\n  ]\n";
        json_file << "}\n";
        json_file.close();
    }
    
    void generateHTMLReport() {
        std::ofstream html_file(output_path + "_clinical_fq_report.html");
        
        html_file << "<!DOCTYPE html>\n<html>\n<head>\n";
        html_file << "<title>Clinical Fluoroquinolone Resistance Report</title>\n";
        html_file << "<style>\n";
        html_file << "body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }\n";
        html_file << ".container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }\n";
        html_file << "h1, h2, h3 { color: #2c3e50; }\n";
        html_file << ".alert { padding: 15px; margin: 20px 0; border-radius: 5px; }\n";
        html_file << ".alert-danger { background-color: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }\n";
        html_file << ".alert-warning { background-color: #fff3cd; color: #856404; border: 1px solid #ffeeba; }\n";
        html_file << ".alert-success { background-color: #d4edda; color: #155724; border: 1px solid #c3e6cb; }\n";
        html_file << ".summary-box { background-color: #e9ecef; padding: 15px; border-radius: 5px; margin: 20px 0; }\n";
        html_file << "table { border-collapse: collapse; width: 100%; margin-top: 20px; }\n";
        html_file << "th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }\n";
        html_file << "th { background-color: #343a40; color: white; }\n";
        html_file << "tr:nth-child(even) { background-color: #f2f2f2; }\n";
        html_file << ".resistance { background-color: #ffcccc; font-weight: bold; }\n";
        html_file << ".qrdr { background-color: #ffe6cc; }\n";
        html_file << ".confidence-bar { width: 200px; height: 20px; background-color: #e0e0e0; border-radius: 10px; display: inline-block; }\n";
        html_file << ".confidence-fill { height: 100%; border-radius: 10px; }\n";
        html_file << ".high-conf { background-color: #dc3545; }\n";
        html_file << ".med-conf { background-color: #ffc107; }\n";
        html_file << ".low-conf { background-color: #28a745; }\n";
        html_file << "</style>\n</head>\n<body>\n";
        
        html_file << "<div class='container'>\n";
        html_file << "<h1>Clinical Fluoroquinolone Resistance Report</h1>\n";
        html_file << "<p><strong>Generated:</strong> " << getCurrentTimestamp() << "</p>\n";
        
        // Alert box based on results
        std::string alert_class = "alert-success";
        if (report.has_high_confidence_resistance) {
            alert_class = "alert-danger";
        } else if (report.has_fq_resistance || report.has_qrdr_mutations) {
            alert_class = "alert-warning";
        }
        
        html_file << "<div class='alert " << alert_class << "'>\n";
        html_file << "<h2>" << report.overall_interpretation << "</h2>\n";
        html_file << "<p>Confidence: ";
        html_file << "<div class='confidence-bar'>";
        html_file << "<div class='confidence-fill ";
        if (report.resistance_confidence > 0.8) html_file << "high-conf";
        else if (report.resistance_confidence > 0.5) html_file << "med-conf";
        else html_file << "low-conf";
        html_file << "' style='width: " << (report.resistance_confidence * 100) << "%'></div>";
        html_file << "</div> " << std::fixed << std::setprecision(0) << (report.resistance_confidence * 100) << "%</p>\n";
        html_file << "</div>\n";
        
        // Clinical notes
        if (!report.clinical_notes.empty()) {
            html_file << "<div class='summary-box'>\n";
            html_file << "<h3>Clinical Recommendations</h3>\n";
            html_file << "<ul>\n";
            for (const auto& note : report.clinical_notes) {
                html_file << "<li>" << note << "</li>\n";
            }
            html_file << "</ul>\n";
            html_file << "</div>\n";
        }
        
        // Summary statistics
        html_file << "<div class='summary-box'>\n";
        html_file << "<h3>Analysis Summary</h3>\n";
        html_file << "<p>Total reads analyzed: " << report.total_reads_analyzed << "</p>\n";
        html_file << "<p>Reads with protein matches: " << report.reads_with_protein_matches 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * report.reads_with_protein_matches / std::max(1, report.total_reads_analyzed)) 
                  << "%)</p>\n";
        html_file << "<p>Reads with mutations: " << report.reads_with_any_mutations << "</p>\n";
        html_file << "<p>Reads with FQ resistance mutations: " << report.reads_with_fq_resistance << "</p>\n";
        html_file << "<p><strong>Performance:</strong> " << std::fixed << std::setprecision(0) 
                  << report.reads_per_second << " reads/second";
        if (report.processing_time_seconds > 0) {
            html_file << " (" << report.total_reads_analyzed << " reads in " 
                      << std::fixed << std::setprecision(1) << report.processing_time_seconds << " seconds)";
        }
        html_file << "</p>\n";
        html_file << "</div>\n";
        
        // Species breakdown
        if (!report.species_summaries.empty()) {
            html_file << "<h3>Species Analysis</h3>\n";
            html_file << "<table>\n";
            html_file << "<tr><th>Species</th><th>Status</th><th>FQ Resistance Mutations</th>"
                      << "<th>QRDR Mutations</th><th>Genes Affected</th></tr>\n";
            
            for (const auto& [species, summary] : report.species_summaries) {
                html_file << "<tr";
                if (summary.likely_resistant) html_file << " class='resistance'";
                html_file << ">\n";
                html_file << "<td>" << species << "</td>\n";
                html_file << "<td>" << (summary.likely_resistant ? "RESISTANT" : "Susceptible") << "</td>\n";
                html_file << "<td>" << summary.fq_resistance_mutations << "</td>\n";
                html_file << "<td>" << summary.qrdr_mutations << "</td>\n";
                html_file << "<td>";
                bool first = true;
                for (const auto& [gene, muts] : summary.gene_mutations) {
                    if (!first) html_file << ", ";
                    html_file << gene;
                    first = false;
                }
                html_file << "</td>\n";
                html_file << "</tr>\n";
            }
            html_file << "</table>\n";
        }
        
        // Detailed mutations table
        html_file << "<h3>Mutations Detected</h3>\n";
        html_file << "<table>\n";
        html_file << "<tr><th>Species</th><th>Gene</th><th>Mutation</th><th>Position</th>"
                  << "<th>Occurrences</th><th>Avg Identity</th><th>Clinical Significance</th></tr>\n";
        
        // Use the same sorted mutations from JSON
        std::map<std::string, ClinicalFQReport::GeneMutation> aggregated_mutations;
        for (const auto& mut : report.all_mutations) {
            std::string key = mut.species + "_" + mut.gene + "_" + mut.mutation_code;
            if (aggregated_mutations.find(key) == aggregated_mutations.end()) {
                aggregated_mutations[key] = mut;
            } else {
                aggregated_mutations[key].occurrence_count++;
            }
        }
        
        std::vector<std::pair<std::string, ClinicalFQReport::GeneMutation>> sorted_mutations(
            aggregated_mutations.begin(), aggregated_mutations.end()
        );
        std::sort(sorted_mutations.begin(), sorted_mutations.end(),
            [](const auto& a, const auto& b) {
                if (a.second.is_known_fq_resistance != b.second.is_known_fq_resistance)
                    return a.second.is_known_fq_resistance;
                if (a.second.is_in_qrdr != b.second.is_in_qrdr)
                    return a.second.is_in_qrdr;
                return a.second.occurrence_count > b.second.occurrence_count;
            }
        );
        
        for (const auto& [key, mut] : sorted_mutations) {
            html_file << "<tr";
            if (mut.is_known_fq_resistance) html_file << " class='resistance'";
            else if (mut.is_in_qrdr) html_file << " class='qrdr'";
            html_file << ">\n";
            html_file << "<td>" << mut.species << "</td>\n";
            html_file << "<td>" << mut.gene << "</td>\n";
            html_file << "<td>" << mut.mutation_code << "</td>\n";
            html_file << "<td>" << mut.position << "</td>\n";
            html_file << "<td>" << mut.occurrence_count << "</td>\n";
            html_file << "<td>" << std::fixed << std::setprecision(1) << (mut.avg_identity * 100) << "%</td>\n";
            html_file << "<td>" << mut.clinical_significance << "</td>\n";
            html_file << "</tr>\n";
        }
        
        html_file << "</table>\n";
        
        // Allele Frequency Analysis section
        if (report.has_allele_frequency_data && !report.allele_frequencies.empty()) {
            html_file << "<h3>Allele Frequency Analysis</h3>\n";
            html_file << "<p>This section shows the frequency of amino acid variants at key resistance positions.</p>\n";
            
            // Sort allele frequencies by resistance frequency
            std::vector<ClinicalFQReport::AlleleFrequencyEntry> sorted_freqs = report.allele_frequencies;
            std::sort(sorted_freqs.begin(), sorted_freqs.end(),
                [](const auto& a, const auto& b) {
                    if (a.has_resistance_mutation != b.has_resistance_mutation)
                        return a.has_resistance_mutation;
                    return a.total_resistant_frequency > b.total_resistant_frequency;
                }
            );
            
            html_file << "<table>\n";
            html_file << "<tr><th>Species</th><th>Gene</th><th>Position</th><th>Depth</th>"
                      << "<th>Wildtype</th><th>WT Freq</th><th>Dominant Mutant</th><th>Mut Freq</th>"
                      << "<th>Resistance Freq</th><th>Mutation Summary</th></tr>\n";
            
            for (const auto& freq : sorted_freqs) {
                html_file << "<tr";
                if (freq.has_resistance_mutation && freq.total_resistant_frequency > 0.1) {
                    html_file << " class='resistance'";
                } else if (freq.has_resistance_mutation) {
                    html_file << " class='qrdr'";
                }
                html_file << ">\n";
                
                html_file << "<td>" << freq.species << "</td>\n";
                html_file << "<td>" << freq.gene << "</td>\n";
                html_file << "<td>" << freq.position << "</td>\n";
                html_file << "<td>" << freq.total_depth << "</td>\n";
                html_file << "<td>" << freq.wildtype_aa << "</td>\n";
                html_file << "<td>" << std::fixed << std::setprecision(1) 
                          << (freq.wildtype_frequency * 100) << "%</td>\n";
                html_file << "<td>" << freq.dominant_mutant_aa << "</td>\n";
                html_file << "<td>" << std::fixed << std::setprecision(1) 
                          << (freq.dominant_mutant_frequency * 100) << "%</td>\n";
                html_file << "<td";
                if (freq.total_resistant_frequency > 0.5) {
                    html_file << " style='color: #dc3545; font-weight: bold;'";
                } else if (freq.total_resistant_frequency > 0.1) {
                    html_file << " style='color: #ffc107; font-weight: bold;'";
                }
                html_file << ">" << std::fixed << std::setprecision(1) 
                          << (freq.total_resistant_frequency * 100) << "%</td>\n";
                html_file << "<td>" << freq.mutation_summary << "</td>\n";
                html_file << "</tr>\n";
            }
            
            html_file << "</table>\n";
            
            // Add interpretation for allele frequencies
            html_file << "<div class='summary-box'>\n";
            html_file << "<h4>Allele Frequency Interpretation</h4>\n";
            html_file << "<ul>\n";
            
            // Count high-frequency resistance
            int high_freq_resistance = 0;
            int moderate_freq_resistance = 0;
            for (const auto& freq : report.allele_frequencies) {
                if (freq.has_resistance_mutation) {
                    if (freq.total_resistant_frequency > 0.5) {
                        high_freq_resistance++;
                    } else if (freq.total_resistant_frequency > 0.1) {
                        moderate_freq_resistance++;
                    }
                }
            }
            
            if (high_freq_resistance > 0) {
                html_file << "<li><strong>High-frequency resistance detected:</strong> " 
                          << high_freq_resistance << " position(s) with >50% resistant alleles. "
                          << "This indicates established resistance in the population.</li>\n";
            }
            if (moderate_freq_resistance > 0) {
                html_file << "<li><strong>Moderate-frequency resistance detected:</strong> " 
                          << moderate_freq_resistance << " position(s) with 10-50% resistant alleles. "
                          << "This may indicate emerging resistance or mixed populations.</li>\n";
            }
            if (high_freq_resistance == 0 && moderate_freq_resistance == 0) {
                html_file << "<li>No significant resistance allele frequencies detected.</li>\n";
            }
            
            html_file << "</ul>\n";
            html_file << "</div>\n";
        }
        
        html_file << "</div>\n";
        html_file << "</body>\n</html>\n";
        html_file.close();
    }
    
    void generateTextReport() {
        std::ofstream text_file(output_path + "_clinical_fq_report.txt");
        
        text_file << "=====================================\n";
        text_file << "CLINICAL FLUOROQUINOLONE RESISTANCE REPORT\n";
        text_file << "=====================================\n\n";
        text_file << "Generated: " << getCurrentTimestamp() << "\n\n";
        
        text_file << "CLINICAL INTERPRETATION\n";
        text_file << "----------------------\n";
        text_file << report.overall_interpretation << "\n";
        text_file << "Confidence: " << std::fixed << std::setprecision(0) 
                  << (report.resistance_confidence * 100) << "%\n\n";
        
        if (!report.clinical_notes.empty()) {
            text_file << "RECOMMENDATIONS\n";
            text_file << "---------------\n";
            for (const auto& note : report.clinical_notes) {
                text_file << "• " << note << "\n";
            }
            text_file << "\n";
        }
        
        text_file << "SUMMARY STATISTICS\n";
        text_file << "-----------------\n";
        text_file << "Total reads analyzed: " << report.total_reads_analyzed << "\n";
        text_file << "Reads with protein matches: " << report.reads_with_protein_matches << "\n";
        text_file << "Reads with mutations: " << report.reads_with_any_mutations << "\n";
        text_file << "Reads with FQ resistance: " << report.reads_with_fq_resistance << "\n";
        text_file << "Performance: " << std::fixed << std::setprecision(0) 
                  << report.reads_per_second << " reads/second";
        if (report.processing_time_seconds > 0) {
            text_file << " (" << report.total_reads_analyzed << " reads in " 
                      << std::fixed << std::setprecision(1) << report.processing_time_seconds << " seconds)";
        }
        text_file << "\n\n";
        
        text_file << "SPECIES BREAKDOWN\n";
        text_file << "-----------------\n";
        for (const auto& [species, summary] : report.species_summaries) {
            text_file << species << ": " 
                      << (summary.likely_resistant ? "RESISTANT" : "Susceptible") << "\n";
            text_file << "  FQ resistance mutations: " << summary.fq_resistance_mutations << "\n";
            text_file << "  QRDR mutations: " << summary.qrdr_mutations << "\n";
            text_file << "  Genes affected: ";
            bool first = true;
            for (const auto& [gene, muts] : summary.gene_mutations) {
                if (!first) text_file << ", ";
                text_file << gene;
                first = false;
            }
            text_file << "\n\n";
        }
        
        text_file << "DETAILED MUTATIONS\n";
        text_file << "-----------------\n";
        
        // Group mutations by clinical significance
        std::vector<ClinicalFQReport::GeneMutation> fq_mutations, qrdr_mutations, other_mutations;
        for (const auto& mut : report.all_mutations) {
            if (mut.is_known_fq_resistance) {
                fq_mutations.push_back(mut);
            } else if (mut.is_in_qrdr) {
                qrdr_mutations.push_back(mut);
            } else {
                other_mutations.push_back(mut);
            }
        }
        
        if (!fq_mutations.empty()) {
            text_file << "\nKNOWN FQ RESISTANCE MUTATIONS:\n";
            for (const auto& mut : fq_mutations) {
                text_file << "  " << mut.species << " " << mut.gene << " " << mut.mutation_code 
                          << " - " << mut.clinical_significance << "\n";
            }
        }
        
        if (!qrdr_mutations.empty()) {
            text_file << "\nQRDR MUTATIONS (uncertain significance):\n";
            for (const auto& mut : qrdr_mutations) {
                text_file << "  " << mut.species << " " << mut.gene << " " << mut.mutation_code << "\n";
            }
        }
        
        text_file.close();
    }
    
    std::string getCurrentTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }
};

// C interface for integration
extern "C" {
    void* create_clinical_fq_report_generator(const char* output_path) {
        return new ClinicalFQReportGenerator(std::string(output_path));
    }
    
    void destroy_clinical_report_generator(void* generator) {
        if (generator) {
            delete static_cast<ClinicalFQReportGenerator*>(generator);
        }
    }
    
    void process_match_for_clinical_report(void* generator, uint32_t read_id,
                                         const char* species_name, const char* gene_name,
                                         uint32_t gene_id, uint32_t species_id,
                                         int num_mutations, const int* positions,
                                         const char* ref_aas, const char* mut_aas,
                                         float alignment_score, float identity,
                                         int match_length, int ref_start, int ref_end) {
        if (generator && num_mutations > 0) {
            ClinicalFQReportGenerator* gen = static_cast<ClinicalFQReportGenerator*>(generator);
            
            std::vector<std::tuple<int, char, char>> mutations;
            for (int i = 0; i < num_mutations; i++) {
                mutations.push_back(std::make_tuple(positions[i], ref_aas[i], mut_aas[i]));
            }
            
            gen->processProteinMatch(read_id, std::string(species_name), std::string(gene_name),
                                   gene_id, species_id, mutations, alignment_score, identity,
                                   match_length, ref_start, ref_end);
        }
    }
    
    void update_clinical_report_stats(void* generator, int total_reads, int reads_with_matches) {
        if (generator) {
            static_cast<ClinicalFQReportGenerator*>(generator)->updateReadStats(total_reads, reads_with_matches);
        }
    }
    
    void update_clinical_report_performance(void* generator, double processing_seconds, double reads_per_sec) {
        if (generator) {
            static_cast<ClinicalFQReportGenerator*>(generator)->updatePerformanceStats(processing_seconds, reads_per_sec);
        }
    }
    
    void generate_clinical_report(void* generator) {
        if (generator) {
            static_cast<ClinicalFQReportGenerator*>(generator)->generateReport();
        }
    }
    
    void add_allele_frequency_to_clinical_report(void* generator, 
                                                  const char* species, const char* gene,
                                                  uint16_t position, uint32_t total_depth,
                                                  char wildtype_aa, uint32_t wildtype_count, 
                                                  float wildtype_frequency,
                                                  char dominant_mutant_aa, uint32_t dominant_mutant_count, 
                                                  float dominant_mutant_frequency, 
                                                  float total_resistant_frequency,
                                                  bool has_resistance_mutation, 
                                                  const char* mutation_summary) {
        if (generator) {
            // For simplicity, we'll pass an empty aa_counts map
            // In production, you might want to pass the full counts
            std::map<char, uint32_t> aa_counts;
            if (wildtype_aa != 'X' && wildtype_count > 0) {
                aa_counts[wildtype_aa] = wildtype_count;
            }
            if (dominant_mutant_aa != 'X' && dominant_mutant_count > 0) {
                aa_counts[dominant_mutant_aa] = dominant_mutant_count;
            }
            
            static_cast<ClinicalFQReportGenerator*>(generator)->addAlleleFrequencyData(
                std::string(species), std::string(gene), position, total_depth,
                wildtype_aa, wildtype_count, wildtype_frequency,
                dominant_mutant_aa, dominant_mutant_count, dominant_mutant_frequency,
                total_resistant_frequency, has_resistance_mutation,
                std::string(mutation_summary), aa_counts);
        }
    }
}