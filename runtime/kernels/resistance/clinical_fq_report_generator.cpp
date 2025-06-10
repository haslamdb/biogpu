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
    
    // Generate shared CSS styles for all HTML reports
    std::string generateCommonCSS() {
        return R"(
        body { 
            font-family: 'Segoe UI', Arial, sans-serif; 
            margin: 0; 
            padding: 20px; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }
        .container { 
            max-width: 1200px; 
            margin: 0 auto; 
            background-color: white; 
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1); 
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        .header h1 { 
            margin: 0; 
            font-size: 2.5em; 
            font-weight: 300;
        }
        .header .subtitle { 
            font-size: 1.2em; 
            opacity: 0.9; 
            margin-top: 10px;
        }
        .content { 
            padding: 30px; 
        }
        h2, h3 { 
            color: #2c3e50; 
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 10px;
        }
        .alert { 
            padding: 20px; 
            margin: 20px 0; 
            border-radius: 8px; 
            border-left: 5px solid;
        }
        .alert-danger { 
            background-color: #fff5f5; 
            color: #c53030; 
            border-left-color: #e53e3e;
        }
        .alert-warning { 
            background-color: #fffaf0; 
            color: #d69e2e; 
            border-left-color: #ed8936;
        }
        .alert-success { 
            background-color: #f0fff4; 
            color: #38a169; 
            border-left-color: #48bb78;
        }
        .summary-box { 
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            padding: 25px; 
            border-radius: 8px; 
            margin: 20px 0; 
            border: 1px solid #dee2e6;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        .stat-card {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            border-left: 4px solid #3498db;
        }
        .stat-number {
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
        }
        .stat-label {
            color: #7f8c8d;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin-top: 20px; 
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        th, td { 
            padding: 12px 15px; 
            text-align: left; 
            border-bottom: 1px solid #ecf0f1;
        }
        th { 
            background: linear-gradient(135deg, #34495e 0%, #2c3e50 100%);
            color: white; 
            font-weight: 500;
            text-transform: uppercase;
            font-size: 0.85em;
            letter-spacing: 0.5px;
        }
        tr:hover { 
            background-color: #f8f9fa; 
        }
        .resistance { 
            background: linear-gradient(135deg, #ffe6e6 0%, #ffcccc 100%);
            font-weight: bold; 
        }
        .qrdr { 
            background: linear-gradient(135deg, #fff0e6 0%, #ffe6cc 100%);
        }
        .confidence-bar { 
            width: 200px; 
            height: 20px; 
            background-color: #ecf0f1; 
            border-radius: 10px; 
            display: inline-block; 
            overflow: hidden;
        }
        .confidence-fill { 
            height: 100%; 
            border-radius: 10px; 
            transition: width 0.3s ease;
        }
        .high-conf { 
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
        }
        .med-conf { 
            background: linear-gradient(135deg, #f39c12 0%, #e67e22 100%);
        }
        .low-conf { 
            background: linear-gradient(135deg, #27ae60 0%, #229954 100%);
        }
        .nav-links {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            text-align: center;
        }
        .nav-links a {
            display: inline-block;
            margin: 10px 20px;
            padding: 12px 24px;
            background: linear-gradient(135deg, #3498db 0%, #2980b9 100%);
            color: white;
            text-decoration: none;
            border-radius: 6px;
            font-weight: 500;
            transition: transform 0.2s ease;
        }
        .nav-links a:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(52, 152, 219, 0.3);
        }
        .back-link {
            display: inline-block;
            margin-bottom: 20px;
            padding: 8px 16px;
            background: #6c757d;
            color: white;
            text-decoration: none;
            border-radius: 4px;
            font-size: 0.9em;
        }
        .back-link:hover {
            background: #5a6268;
        }
        .footer {
            background: #f8f9fa;
            padding: 20px;
            text-align: center;
            color: #6c757d;
            font-size: 0.9em;
            border-top: 1px solid #dee2e6;
        }
        .freq-bar { 
            width: 100px; 
            height: 20px; 
            background-color: #ecf0f1; 
            border-radius: 10px; 
            display: inline-block; 
            overflow: hidden; 
        }
        .freq-fill { 
            height: 100%; 
            border-radius: 10px; 
        }
        .freq-high { 
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%); 
        }
        .freq-med { 
            background: linear-gradient(135deg, #f39c12 0%, #e67e22 100%); 
        }
        .freq-low { 
            background: linear-gradient(135deg, #3498db 0%, #2980b9 100%); 
        }
        )";
    }

    void generateCoverPageHTML() {
        std::ofstream html_file(output_path + "_clinical_fq_report.html");
        
        // Extract sample name from output path
        std::string sample_name = output_path;
        size_t last_slash = sample_name.find_last_of("/\\");
        if (last_slash != std::string::npos) {
            sample_name = sample_name.substr(last_slash + 1);
        }
        
        html_file << "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n";
        html_file << "<meta charset=\"UTF-8\">\n";
        html_file << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
        html_file << "<title>Clinical FQ Resistance Report - " << sample_name << "</title>\n";
        html_file << "<style>\n" << generateCommonCSS() << "</style>\n";
        html_file << "</head>\n<body>\n";
        
        html_file << "<div class='container'>\n";
        
        // Header
        html_file << "<div class='header'>\n";
        html_file << "<h1>Clinical Fluoroquinolone Resistance Report</h1>\n";
        html_file << "<div class='subtitle'>Sample: " << sample_name << "</div>\n";
        html_file << "<div class='subtitle'>Generated: " << getCurrentTimestamp() << "</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='content'>\n";
        
        // Main alert box
        std::string alert_class = "alert-success";
        if (report.has_high_confidence_resistance) {
            alert_class = "alert-danger";
        } else if (report.has_fq_resistance || report.has_qrdr_mutations) {
            alert_class = "alert-warning";
        }
        
        html_file << "<div class='alert " << alert_class << "'>\n";
        html_file << "<h2 style='margin-top: 0; border: none; padding: 0;'>" << report.overall_interpretation << "</h2>\n";
        html_file << "<p style='margin-bottom: 0;'>Confidence: ";
        html_file << "<div class='confidence-bar' style='margin-left: 10px;'>";
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
        
        // Statistics grid
        html_file << "<h3>Analysis Summary</h3>\n";
        html_file << "<div class='stats-grid'>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << report.total_reads_analyzed << "</div>\n";
        html_file << "<div class='stat-label'>Total Reads Analyzed</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << report.reads_with_protein_matches << "</div>\n";
        html_file << "<div class='stat-label'>Protein Matches (" 
                  << std::fixed << std::setprecision(1) 
                  << (100.0 * report.reads_with_protein_matches / std::max(1, report.total_reads_analyzed)) 
                  << "%)</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << report.reads_with_any_mutations << "</div>\n";
        html_file << "<div class='stat-label'>Reads with Mutations</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << report.reads_with_fq_resistance << "</div>\n";
        html_file << "<div class='stat-label'>FQ Resistance Reads</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << std::fixed << std::setprecision(0) << report.reads_per_second << "</div>\n";
        html_file << "<div class='stat-label'>Reads/Second";
        if (report.processing_time_seconds > 0) {
            html_file << " (" << std::fixed << std::setprecision(1) << report.processing_time_seconds << "s total)";
        }
        html_file << "</div>\n";
        html_file << "</div>\n";
        
        // Count allele frequencies for summary
        int total_positions = 0;
        int resistant_positions = 0;
        if (report.has_allele_frequency_data) {
            for (const auto& freq : report.allele_frequencies) {
                total_positions++;
                if (freq.has_resistance_mutation && freq.total_resistant_frequency > 0.1) {
                    resistant_positions++;
                }
            }
        }
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << total_positions << "</div>\n";
        html_file << "<div class='stat-label'>Resistance Positions Analyzed</div>\n";
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End stats-grid
        
        // Species breakdown summary
        if (!report.species_summaries.empty()) {
            html_file << "<div class='summary-box'>\n";
            html_file << "<h3>Species Summary</h3>\n";
            for (const auto& [species, summary] : report.species_summaries) {
                html_file << "<div style='margin: 10px 0; padding: 10px; border-left: 4px solid " 
                          << (summary.likely_resistant ? "#e74c3c" : "#27ae60") << ";'>\n";
                html_file << "<strong>" << species << ":</strong> " 
                          << (summary.likely_resistant ? "RESISTANT" : "Susceptible");
                if (summary.fq_resistance_mutations > 0) {
                    html_file << " (" << summary.fq_resistance_mutations << " FQ resistance mutations)";
                }
                html_file << "</div>\n";
            }
            html_file << "</div>\n";
        }
        
        // Navigation links
        html_file << "<div class='nav-links'>\n";
        html_file << "<h3>Detailed Reports</h3>\n";
        html_file << "<a href='" << sample_name << "_mutations_report.html'>📊 Mutations Analysis</a>\n";
        if (report.has_allele_frequency_data && !report.allele_frequencies.empty()) {
            html_file << "<a href='" << sample_name << "_allele_frequencies_report.html'>🧬 Allele Frequencies</a>\n";
        }
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End content
        
        // Footer
        html_file << "<div class='footer'>\n";
        html_file << "Report generated by Clinical FQ Resistance Pipeline v1.0<br>\n";
        html_file << "For research use only. Clinical decisions should be based on validated diagnostic methods.\n";
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End container
        html_file << "</body>\n</html>\n";
        html_file.close();
    }
    
    void generateMutationsReportHTML() {
        // Extract sample name from output path
        std::string sample_name = output_path;
        size_t last_slash = sample_name.find_last_of("/\\");
        if (last_slash != std::string::npos) {
            sample_name = sample_name.substr(last_slash + 1);
        }
        
        std::ofstream html_file(output_path + "_mutations_report.html");
        
        html_file << "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n";
        html_file << "<meta charset=\"UTF-8\">\n";
        html_file << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
        html_file << "<title>Mutations Analysis - " << sample_name << "</title>\n";
        html_file << "<style>\n" << generateCommonCSS() << "</style>\n";
        html_file << "</head>\n<body>\n";
        
        html_file << "<div class='container'>\n";
        
        // Header
        html_file << "<div class='header'>\n";
        html_file << "<h1>Mutations Analysis</h1>\n";
        html_file << "<div class='subtitle'>Sample: " << sample_name << "</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='content'>\n";
        
        // Back link
        html_file << "<a href='" << sample_name << "_clinical_fq_report.html' class='back-link'>← Back to Summary</a>\n";
        
        // Species breakdown table
        if (!report.species_summaries.empty()) {
            html_file << "<h3>Species Analysis</h3>\n";
            html_file << "<table>\n";
            html_file << "<tr><th>Species</th><th>Status</th><th>FQ Resistance Mutations</th>"
                      << "<th>QRDR Mutations</th><th>Genes Affected</th><th>Max Identity</th></tr>\n";
            
            for (const auto& [species, summary] : report.species_summaries) {
                html_file << "<tr";
                if (summary.likely_resistant) html_file << " class='resistance'";
                html_file << ">\n";
                html_file << "<td><strong>" << species << "</strong></td>\n";
                html_file << "<td>" << (summary.likely_resistant ? "🔴 RESISTANT" : "🟢 Susceptible") << "</td>\n";
                html_file << "<td>" << summary.fq_resistance_mutations << "</td>\n";
                html_file << "<td>" << summary.qrdr_mutations << "</td>\n";
                html_file << "<td>";
                bool first = true;
                for (const auto& [gene, muts] : summary.gene_mutations) {
                    if (!first) html_file << ", ";
                    html_file << "<em>" << gene << "</em>";
                    first = false;
                }
                html_file << "</td>\n";
                html_file << "<td>" << std::fixed << std::setprecision(1) << (summary.max_identity_score * 100) << "%</td>\n";
                html_file << "</tr>\n";
            }
            html_file << "</table>\n";
        }
        
        // Detailed mutations table
        html_file << "<h3>Detailed Mutations</h3>\n";
        
        // Add filter legend
        html_file << "<div class='summary-box'>\n";
        html_file << "<h4>Legend</h4>\n";
        html_file << "<div style='display: flex; gap: 20px; flex-wrap: wrap;'>\n";
        html_file << "<div style='display: flex; align-items: center;'><div style='width: 20px; height: 20px; background: linear-gradient(135deg, #ffe6e6 0%, #ffcccc 100%); border-radius: 4px; margin-right: 8px;'></div>Known FQ Resistance</div>\n";
        html_file << "<div style='display: flex; align-items: center;'><div style='width: 20px; height: 20px; background: linear-gradient(135deg, #fff0e6 0%, #ffe6cc 100%); border-radius: 4px; margin-right: 8px;'></div>QRDR Region</div>\n";
        html_file << "<div style='display: flex; align-items: center;'><div style='width: 20px; height: 20px; background: white; border: 1px solid #ddd; border-radius: 4px; margin-right: 8px;'></div>Other Mutations</div>\n";
        html_file << "</div>\n";
        html_file << "</div>\n";
        
        html_file << "<table>\n";
        html_file << "<tr><th>Species</th><th>Gene</th><th>Mutation</th><th>Position</th>"
                  << "<th>Change</th><th>Occurrences</th><th>Avg Identity</th><th>Clinical Significance</th></tr>\n";
        
        // Use the same sorted mutations logic
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
            html_file << "<td><strong>" << mut.species << "</strong></td>\n";
            html_file << "<td><em>" << mut.gene << "</em></td>\n";
            html_file << "<td><code>" << mut.mutation_code << "</code></td>\n";
            html_file << "<td>" << mut.position << "</td>\n";
            html_file << "<td>" << mut.wildtype_aa << " → " << mut.mutant_aa << "</td>\n";
            html_file << "<td>" << mut.occurrence_count << "</td>\n";
            html_file << "<td>" << std::fixed << std::setprecision(1) << (mut.avg_identity * 100) << "%</td>\n";
            html_file << "<td>" << mut.clinical_significance << "</td>\n";
            html_file << "</tr>\n";
        }
        
        html_file << "</table>\n";
        
        // Summary statistics
        html_file << "<div class='summary-box'>\n";
        html_file << "<h4>Mutation Summary</h4>\n";
        int fq_count = 0, qrdr_count = 0, other_count = 0;
        for (const auto& [key, mut] : sorted_mutations) {
            if (mut.is_known_fq_resistance) fq_count++;
            else if (mut.is_in_qrdr) qrdr_count++;
            else other_count++;
        }
        html_file << "<p><strong>Known FQ Resistance Mutations:</strong> " << fq_count << "</p>\n";
        html_file << "<p><strong>QRDR Mutations:</strong> " << qrdr_count << "</p>\n";
        html_file << "<p><strong>Other Mutations:</strong> " << other_count << "</p>\n";
        html_file << "<p><strong>Total Unique Mutations:</strong> " << sorted_mutations.size() << "</p>\n";
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End content
        
        // Footer
        html_file << "<div class='footer'>\n";
        html_file << "Mutations Report | Generated: " << getCurrentTimestamp() << "\n";
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End container
        html_file << "</body>\n</html>\n";
        html_file.close();
    }
    
    void generateAlleleFrequenciesReportHTML() {
        if (!report.has_allele_frequency_data || report.allele_frequencies.empty()) {
            return; // Don't generate if no data
        }
        
        // Extract sample name from output path
        std::string sample_name = output_path;
        size_t last_slash = sample_name.find_last_of("/\\");
        if (last_slash != std::string::npos) {
            sample_name = sample_name.substr(last_slash + 1);
        }
        
        std::ofstream html_file(output_path + "_allele_frequencies_report.html");
        
        html_file << "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n";
        html_file << "<meta charset=\"UTF-8\">\n";
        html_file << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
        html_file << "<title>Allele Frequencies - " << sample_name << "</title>\n";
        html_file << "<style>\n" << generateCommonCSS() << "</style>\n";
        html_file << "</head>\n<body>\n";
        
        html_file << "<div class='container'>\n";
        
        // Header
        html_file << "<div class='header'>\n";
        html_file << "<h1>Allele Frequency Analysis</h1>\n";
        html_file << "<div class='subtitle'>Sample: " << sample_name << "</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='content'>\n";
        
        // Back link
        html_file << "<a href='" << sample_name << "_clinical_fq_report.html' class='back-link'>← Back to Summary</a>\n";
        
        html_file << "<div class='summary-box'>\n";
        html_file << "<h3>About Allele Frequency Analysis</h3>\n";
        html_file << "<p>This analysis shows the frequency of amino acid variants at key fluoroquinolone resistance positions. "
                  << "High frequencies of resistance alleles (>50%) indicate established resistance, while moderate frequencies "
                  << "(10-50%) may indicate emerging resistance or mixed populations.</p>\n";
        html_file << "</div>\n";
        
        // Sort allele frequencies by resistance frequency
        std::vector<ClinicalFQReport::AlleleFrequencyEntry> sorted_freqs = report.allele_frequencies;
        std::sort(sorted_freqs.begin(), sorted_freqs.end(),
            [](const auto& a, const auto& b) {
                if (a.has_resistance_mutation != b.has_resistance_mutation)
                    return a.has_resistance_mutation;
                return a.total_resistant_frequency > b.total_resistant_frequency;
            }
        );
        
        // Summary statistics
        int high_freq_resistance = 0;
        int moderate_freq_resistance = 0;
        int total_positions = sorted_freqs.size();
        float max_resistance_freq = 0.0f;
        
        for (const auto& freq : sorted_freqs) {
            if (freq.has_resistance_mutation) {
                if (freq.total_resistant_frequency > 0.5) {
                    high_freq_resistance++;
                } else if (freq.total_resistant_frequency > 0.1) {
                    moderate_freq_resistance++;
                }
            }
            max_resistance_freq = std::max(max_resistance_freq, freq.total_resistant_frequency);
        }
        
        html_file << "<div class='stats-grid'>\n";
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << total_positions << "</div>\n";
        html_file << "<div class='stat-label'>Total Positions</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << high_freq_resistance << "</div>\n";
        html_file << "<div class='stat-label'>High-Freq Resistance (>50%)</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << moderate_freq_resistance << "</div>\n";
        html_file << "<div class='stat-label'>Moderate-Freq Resistance (10-50%)</div>\n";
        html_file << "</div>\n";
        
        html_file << "<div class='stat-card'>\n";
        html_file << "<div class='stat-number'>" << std::fixed << std::setprecision(1) << (max_resistance_freq * 100) << "%</div>\n";
        html_file << "<div class='stat-label'>Maximum Resistance Frequency</div>\n";
        html_file << "</div>\n";
        html_file << "</div>\n";
        
        // Legend
        html_file << "<div class='summary-box'>\n";
        html_file << "<h4>Frequency Legend</h4>\n";
        html_file << "<div style='display: flex; gap: 20px; flex-wrap: wrap;'>\n";
        html_file << "<div style='display: flex; align-items: center;'><div style='width: 20px; height: 20px; background: linear-gradient(135deg, #ffe6e6 0%, #ffcccc 100%); border-radius: 4px; margin-right: 8px;'></div>Resistance Detected</div>\n";
        html_file << "<div style='display: flex; align-items: center;'><div style='width: 20px; height: 20px; background: linear-gradient(135deg, #fff0e6 0%, #ffe6cc 100%); border-radius: 4px; margin-right: 8px;'></div>Low-Level Resistance</div>\n";
        html_file << "<div style='display: flex; align-items: center;'><div style='width: 20px; height: 20px; background: white; border: 1px solid #ddd; border-radius: 4px; margin-right: 8px;'></div>No Resistance</div>\n";
        html_file << "</div>\n";
        html_file << "</div>\n";
        
        html_file << "<h3>Allele Frequencies by Position</h3>\n";
        html_file << "<table>\n";
        html_file << "<tr><th>Species</th><th>Gene</th><th>Position</th><th>Depth</th>"
                  << "<th>Wildtype</th><th>WT Freq</th><th>Dominant Mutant</th><th>Mut Freq</th>"
                  << "<th>Resistance Freq</th><th>Mutation Summary</th></tr>\n";
        
        for (const auto& freq : sorted_freqs) {
            html_file << "<tr";
            if (freq.has_resistance_mutation && freq.total_resistant_frequency > 0.5) {
                html_file << " class='resistance'";
            } else if (freq.has_resistance_mutation && freq.total_resistant_frequency > 0.1) {
                html_file << " class='qrdr'";
            }
            html_file << ">\n";
            
            html_file << "<td><strong>" << freq.species << "</strong></td>\n";
            html_file << "<td><em>" << freq.gene << "</em></td>\n";
            html_file << "<td>" << freq.position << "</td>\n";
            html_file << "<td>" << freq.total_depth << "</td>\n";
            html_file << "<td><code>" << freq.wildtype_aa << "</code></td>\n";
            
            // Wildtype frequency with bar
            html_file << "<td>";
            html_file << "<div class='freq-bar'>";
            html_file << "<div class='freq-fill freq-low' style='width: " << (freq.wildtype_frequency * 100) << "%'></div>";
            html_file << "</div>";
            html_file << " " << std::fixed << std::setprecision(1) << (freq.wildtype_frequency * 100) << "%</td>\n";
            
            html_file << "<td><code>" << freq.dominant_mutant_aa << "</code></td>\n";
            
            // Dominant mutant frequency with bar
            html_file << "<td>";
            html_file << "<div class='freq-bar'>";
            html_file << "<div class='freq-fill freq-med' style='width: " << (freq.dominant_mutant_frequency * 100) << "%'></div>";
            html_file << "</div>";
            html_file << " " << std::fixed << std::setprecision(1) << (freq.dominant_mutant_frequency * 100) << "%</td>\n";
            
            // Resistance frequency with bar and color coding
            html_file << "<td>";
            if (freq.total_resistant_frequency > 0) {
                html_file << "<div class='freq-bar'>";
                std::string freq_class = "freq-low";
                if (freq.total_resistant_frequency > 0.5) freq_class = "freq-high";
                else if (freq.total_resistant_frequency > 0.1) freq_class = "freq-med";
                html_file << "<div class='freq-fill " << freq_class << "' style='width: " << (freq.total_resistant_frequency * 100) << "%'></div>";
                html_file << "</div> ";
                if (freq.total_resistant_frequency > 0.5) {
                    html_file << "<strong style='color: #dc3545;'>";
                } else if (freq.total_resistant_frequency > 0.1) {
                    html_file << "<strong style='color: #ffc107;'>";
                }
                html_file << std::fixed << std::setprecision(1) << (freq.total_resistant_frequency * 100) << "%";
                if (freq.total_resistant_frequency > 0.1) {
                    html_file << "</strong>";
                }
            } else {
                html_file << "0%";
            }
            html_file << "</td>\n";
            
            html_file << "<td>" << freq.mutation_summary << "</td>\n";
            html_file << "</tr>\n";
        }
        
        html_file << "</table>\n";
        
        // Clinical interpretation
        html_file << "<div class='summary-box'>\n";
        html_file << "<h4>Clinical Interpretation</h4>\n";
        html_file << "<ul>\n";
        
        if (high_freq_resistance > 0) {
            html_file << "<li><strong style='color: #dc3545;'>High-frequency resistance detected:</strong> " 
                      << high_freq_resistance << " position(s) with >50% resistant alleles. "
                      << "This indicates established resistance in the population.</li>\n";
        }
        if (moderate_freq_resistance > 0) {
            html_file << "<li><strong style='color: #ffc107;'>Moderate-frequency resistance detected:</strong> " 
                      << moderate_freq_resistance << " position(s) with 10-50% resistant alleles. "
                      << "This may indicate emerging resistance or mixed populations.</li>\n";
        }
        if (high_freq_resistance == 0 && moderate_freq_resistance == 0) {
            html_file << "<li><strong style='color: #28a745;'>No significant resistance allele frequencies detected.</strong> "
                      << "The sample appears susceptible to fluoroquinolones based on allele frequency analysis.</li>\n";
        }
        
        html_file << "</ul>\n";
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End content
        
        // Footer
        html_file << "<div class='footer'>\n";
        html_file << "Allele Frequency Report | Generated: " << getCurrentTimestamp() << "\n";
        html_file << "</div>\n";
        
        html_file << "</div>\n"; // End container
        html_file << "</body>\n</html>\n";
        html_file.close();
    }

    std::string getCurrentTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }
    
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
        generateHTMLReport(); // This now generates the cover page + linked reports
        generateTextReport();
    }
    
    // Updated generateHTMLReport method
    void generateHTMLReport() {
        generateCoverPageHTML();
        generateMutationsReportHTML();
        generateAlleleFrequenciesReportHTML();
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