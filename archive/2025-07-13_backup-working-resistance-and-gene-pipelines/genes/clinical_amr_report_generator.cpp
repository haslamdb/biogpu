// runtime/kernels/genes/clinical_amr_report_generator.cpp
#include "clinical_amr_report_generator.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <chrono>

ClinicalAMRReportGenerator::ClinicalAMRReportGenerator(const std::string& output, const std::string& sample)
    : output_path(output), sample_name(sample), total_reads_processed(0), 
      reads_with_amr_hits(0), total_genes_detected(0), high_confidence_genes(0) {
    initializeDrugClassInterpretations();
    initializeCriticalGeneThresholds();
}

void ClinicalAMRReportGenerator::initializeDrugClassInterpretations() {
    drug_class_interpretations = {
        {"BETA_LACTAM", "Beta-lactam resistance detected. Consider alternative antibiotics such as fluoroquinolones or aminoglycosides if susceptible."},
        {"AMINOGLYCOSIDE", "Aminoglycoside resistance detected. Consider beta-lactams or fluoroquinolones based on susceptibility."},
        {"FLUOROQUINOLONE", "Fluoroquinolone resistance detected. Consider beta-lactams or aminoglycosides based on susceptibility."},
        {"TETRACYCLINE", "Tetracycline resistance detected. Alternative antibiotics should be considered."},
        {"MACROLIDE", "Macrolide resistance detected. Consider alternative antibiotic classes."},
        {"GLYCOPEPTIDE", "Glycopeptide resistance detected. This is a serious concern - consult infectious disease specialist."},
        {"COLISTIN", "Colistin resistance detected. Critical resistance - immediate consultation with infectious disease specialist recommended."},
        {"SULFONAMIDE", "Sulfonamide resistance detected. Consider alternative antibiotics."},
        {"UNKNOWN", "Resistance genes of unknown class detected. Further investigation recommended."}
    };
}

void ClinicalAMRReportGenerator::initializeCriticalGeneThresholds() {
    // Critical genes requiring higher identity thresholds
    critical_gene_thresholds = {
        {"mcr", 0.98f},      // Colistin resistance - extremely critical
        {"vanA", 0.97f},     // Vancomycin resistance
        {"vanB", 0.97f},     // Vancomycin resistance
        {"blaKPC", 0.96f},   // Carbapenem resistance
        {"blaNDM", 0.96f},   // Carbapenem resistance
        {"blaOXA-48", 0.96f} // Carbapenem resistance
    };
}

void ClinicalAMRReportGenerator::processAMRResults(const std::vector<AMRHit>& hits,
                                                   const std::vector<AMRCoverageStats>& coverage_stats,
                                                   const std::vector<AMRGeneEntry>& gene_entries,
                                                   int total_reads) {
    total_reads_processed = total_reads;
    
    // Debug logging
    std::cout << "\n=== Clinical Report Generator Debug ===" << std::endl;
    std::cout << "Total AMR hits received: " << hits.size() << std::endl;
    std::cout << "Total coverage stats entries: " << coverage_stats.size() << std::endl;
    std::cout << "Total gene entries: " << gene_entries.size() << std::endl;
    std::cout << "Total reads processed: " << total_reads << std::endl;
    
    // Count reads with AMR hits
    std::set<uint32_t> reads_with_hits;
    for (const auto& hit : hits) {
        reads_with_hits.insert(hit.read_id);
    }
    reads_with_amr_hits = reads_with_hits.size();
    std::cout << "Unique reads with AMR hits: " << reads_with_amr_hits << std::endl;
    
    // Process coverage statistics and build gene summaries
    int genes_with_coverage = 0;
    for (size_t i = 0; i < coverage_stats.size(); i++) {
        const auto& stats = coverage_stats[i];
        if (stats.total_reads > 0) {  // Gene was detected
            genes_with_coverage++;
            const auto& gene_entry = gene_entries[i];
            
            // Debug first few genes
            if (genes_with_coverage <= 5) {
                std::cout << "\nGene " << i << " detected:" << std::endl;
                std::cout << "  Name: " << gene_entry.gene_name << std::endl;
                std::cout << "  Total reads: " << stats.total_reads << std::endl;
                std::cout << "  Coverage: " << stats.percent_coverage << "%" << std::endl;
                std::cout << "  Mean depth: " << stats.mean_depth << std::endl;
            }
            
            GeneSummary summary;
            summary.gene_name = std::string(gene_entry.gene_name);
            summary.gene_family = std::string(gene_entry.gene_family);
            summary.drug_class = std::string(gene_entry.class_);
            summary.read_count = stats.total_reads;
            summary.percent_coverage = stats.percent_coverage;
            summary.mean_depth = stats.mean_depth;
            summary.tpm = stats.tpm;
            summary.rpkm = stats.rpkm;
            summary.is_complete_gene = false;  // Will update based on hits
            
            // Calculate average identity from hits for this gene
            float avg_identity = 0.0f;
            int identity_count = 0;
            for (const auto& hit : hits) {
                if (hit.gene_id == i) {
                    avg_identity += hit.identity;
                    identity_count++;
                }
            }
            if (identity_count > 0) {
                avg_identity /= identity_count;
            }
            
            // Determine confidence level with identity and gene family
            summary.confidence_level = getConfidenceLevel(stats.percent_coverage, stats.mean_depth, avg_identity, summary.gene_family);
            
            // Count high confidence genes based on the confidence level
            if (summary.confidence_level == "HIGH") {
                high_confidence_genes++;
            }
            
            // Check if any hit shows complete gene
            for (const auto& hit : hits) {
                if (hit.gene_id == i && hit.is_complete_gene) {
                    summary.is_complete_gene = true;
                    break;
                }
            }
            
            gene_summaries[summary.gene_name] = summary;
            total_genes_detected++;
            
            // Debug logging for confidence levels
            if (genes_with_coverage <= 5) {
                std::cout << "  Average identity: " << avg_identity << std::endl;
                std::cout << "  Confidence level: " << summary.confidence_level << std::endl;
            }
            
            // Update drug class summary
            auto& drug_summary = drug_class_summaries[summary.drug_class];
            drug_summary.drug_class = summary.drug_class;
            drug_summary.genes_detected.push_back(summary.gene_name);
            drug_summary.total_reads += summary.read_count;
            drug_summary.max_tpm = std::max(drug_summary.max_tpm, summary.tpm);
            
            // Use the already-calculated confidence level for consistency
            if (summary.confidence_level == "HIGH") {
                drug_summary.high_confidence_genes.push_back(summary.gene_name);
            } else if (summary.confidence_level == "MODERATE") {
                drug_summary.moderate_confidence_genes.push_back(summary.gene_name);
            } else {
                drug_summary.low_confidence_genes.push_back(summary.gene_name);
            }
        }
    }
    
    // Aggregate by gene family
    gene_family_summaries.clear();
    std::map<std::string, std::vector<float>> family_identities; // Store all identities for averaging
    
    for (const auto& [gene_name, summary] : gene_summaries) {
        auto& family_summary = gene_family_summaries[summary.gene_family];
        family_summary.gene_family = summary.gene_family;
        family_summary.gene_variants.push_back(gene_name);
        family_summary.total_tpm += summary.tpm;
        family_summary.total_reads += summary.read_count;
        
        // Store identities for this gene family (need to calculate from hits)
        // We'll update this after processing all hits
    }
    
    // Calculate average identity per gene family from hits
    std::map<std::string, float> family_identity_sums;
    std::map<std::string, int> family_identity_counts;
    
    for (const auto& hit : hits) {
        if (hit.gene_id < gene_entries.size()) {
            std::string gene_family(gene_entries[hit.gene_id].gene_family);
            family_identity_sums[gene_family] += hit.identity;
            family_identity_counts[gene_family]++;
        }
    }
    
    // Update family summaries with average identity
    for (auto& [family_name, summary] : gene_family_summaries) {
        if (family_identity_counts[family_name] > 0) {
            summary.mean_identity = family_identity_sums[family_name] / family_identity_counts[family_name];
        } else {
            summary.mean_identity = 0.0f;
        }
    }
    
    // Generate clinical interpretations for each drug class
    for (auto& [drug_class, summary] : drug_class_summaries) {
        summary.clinical_interpretation = generateDrugClassInterpretation(summary);
    }
    
    // Final debug summary
    std::cout << "\n=== Clinical Report Summary ===" << std::endl;
    std::cout << "Total genes with coverage > 0: " << genes_with_coverage << std::endl;
    std::cout << "Total genes detected (in summaries): " << total_genes_detected << std::endl;
    std::cout << "Total drug classes affected: " << drug_class_summaries.size() << std::endl;
    std::cout << "High confidence genes: " << high_confidence_genes << std::endl;
    
    // List drug classes
    std::cout << "\nDrug classes detected:" << std::endl;
    for (const auto& [drug_class, summary] : drug_class_summaries) {
        std::cout << "  " << drug_class << ": " << summary.genes_detected.size() 
                  << " genes (High: " << summary.high_confidence_genes.size() << ")" << std::endl;
    }
    std::cout << "==============================\n" << std::endl;
}

std::string ClinicalAMRReportGenerator::getConfidenceLevel(float coverage, float depth, float identity, const std::string& gene_family) {
    // Check if this is a critical gene requiring higher thresholds
    float min_identity_high = 0.95f;
    auto it = critical_gene_thresholds.find(gene_family);
    if (it != critical_gene_thresholds.end()) {
        min_identity_high = it->second;
    }
    
    // Metagenomics-appropriate criteria focusing on identity and presence
    // Note: coverage is passed as percentage (0-100), not fraction
    if (identity >= min_identity_high && coverage >= 20.0f && depth >= 2.0f) {
        return "HIGH";
    } else if (identity >= 0.90f && coverage >= 10.0f && depth >= 1.0f) {
        return "MODERATE";
    } else {
        return "LOW";
    }
}

std::string ClinicalAMRReportGenerator::generateDrugClassInterpretation(const DrugClassSummary& summary) {
    if (summary.high_confidence_genes.empty() && 
        summary.moderate_confidence_genes.empty() && 
        summary.low_confidence_genes.empty()) {
        return "No resistance genes detected for this drug class.";
    }
    
    std::string interpretation;
    
    if (!summary.high_confidence_genes.empty()) {
        interpretation = "HIGH CONFIDENCE: ";
        auto it = drug_class_interpretations.find(summary.drug_class);
        if (it != drug_class_interpretations.end()) {
            interpretation += it->second;
        } else {
            interpretation += "Resistance detected. Consider alternative antibiotics.";
        }
        
        interpretation += " Genes detected: ";
        for (size_t i = 0; i < summary.high_confidence_genes.size(); i++) {
            if (i > 0) interpretation += ", ";
            interpretation += summary.high_confidence_genes[i];
        }
    } else if (!summary.moderate_confidence_genes.empty()) {
        interpretation = "MODERATE CONFIDENCE: Possible resistance. ";
        interpretation += "Consider susceptibility testing. Genes detected: ";
        for (size_t i = 0; i < summary.moderate_confidence_genes.size(); i++) {
            if (i > 0) interpretation += ", ";
            interpretation += summary.moderate_confidence_genes[i];
        }
    } else {
        interpretation = "LOW CONFIDENCE: Weak evidence of resistance. ";
        interpretation += "Susceptibility testing recommended.";
    }
    
    return interpretation;
}

void ClinicalAMRReportGenerator::generateReports() {
    generateHTMLReport();
    generateTextReport();
    generateJSONReport();
    generateTSVReports();
}

void ClinicalAMRReportGenerator::generateHTMLReport() {
    std::ofstream html(output_path + "_clinical_amr_report.html");
    
    html << "<!DOCTYPE html>\n<html>\n<head>\n";
    html << "<title>Clinical AMR Report - " << sample_name << "</title>\n";
    html << "<style>\n";
    html << "body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }\n";
    html << ".container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; "
         << "border-radius: 10px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }\n";
    html << "h1, h2, h3 { color: #2c3e50; }\n";
    html << ".summary { background-color: #ecf0f1; padding: 15px; border-radius: 5px; margin: 20px 0; }\n";
    html << ".alert { padding: 15px; border-radius: 5px; margin: 10px 0; }\n";
    html << ".alert-danger { background-color: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }\n";
    html << ".alert-warning { background-color: #fff3cd; color: #856404; border: 1px solid #ffeeba; }\n";
    html << ".alert-success { background-color: #d4edda; color: #155724; border: 1px solid #c3e6cb; }\n";
    html << "table { border-collapse: collapse; width: 100%; margin-top: 20px; }\n";
    html << "th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }\n";
    html << "th { background-color: #3498db; color: white; }\n";
    html << "tr:nth-child(even) { background-color: #f2f2f2; }\n";
    html << ".high-confidence { background-color: #d4edda; }\n";
    html << ".moderate-confidence { background-color: #fff3cd; }\n";
    html << ".low-confidence { background-color: #f8d7da; }\n";
    html << ".stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); "
         << "gap: 20px; margin: 20px 0; }\n";
    html << ".stat-box { background-color: #ecf0f1; padding: 20px; border-radius: 5px; text-align: center; }\n";
    html << ".stat-number { font-size: 2em; font-weight: bold; color: #3498db; }\n";
    html << ".stat-label { color: #7f8c8d; margin-top: 5px; }\n";
    html << "</style>\n</head>\n<body>\n";
    
    html << "<div class='container'>\n";
    
    // Header
    html << "<h1>Clinical Antimicrobial Resistance Report</h1>\n";
    html << "<div class='summary'>\n";
    html << "<p><strong>Sample:</strong> " << sample_name << "</p>\n";
    html << "<p><strong>Report Generated:</strong> " << getCurrentTimestamp() << "</p>\n";
    html << "<p><strong>Analysis Type:</strong> Metagenomic AMR Gene Detection</p>\n";
    html << "</div>\n";
    
    // Executive Summary
    html << "<h2>Executive Summary</h2>\n";
    
    if (total_genes_detected == 0) {
        html << "<div class='alert alert-success'>\n";
        html << "<strong>No antimicrobial resistance genes detected.</strong> ";
        html << "Standard antibiotic therapy may be appropriate based on clinical presentation.\n";
        html << "</div>\n";
    } else {
        std::string alert_class = "alert-warning";
        if (high_confidence_genes > 5) {
            alert_class = "alert-danger";
        }
        
        html << "<div class='alert " << alert_class << "'>\n";
        html << "<strong>" << total_genes_detected << " antimicrobial resistance genes detected</strong> ";
        html << "(" << high_confidence_genes << " with high confidence).<br>\n";
        html << "Antibiotic selection should be guided by the resistance profile below.\n";
        html << "</div>\n";
    }
    
    // Statistics Grid
    html << "<div class='stats-grid'>\n";
    
    html << "<div class='stat-box'>\n";
    html << "<div class='stat-number'>" << total_reads_processed << "</div>\n";
    html << "<div class='stat-label'>Total Reads</div>\n";
    html << "</div>\n";
    
    html << "<div class='stat-box'>\n";
    html << "<div class='stat-number'>" << reads_with_amr_hits << "</div>\n";
    html << "<div class='stat-label'>Reads with AMR Genes</div>\n";
    html << "</div>\n";
    
    html << "<div class='stat-box'>\n";
    html << "<div class='stat-number'>" << total_genes_detected << "</div>\n";
    html << "<div class='stat-label'>AMR Genes Detected</div>\n";
    html << "</div>\n";
    
    html << "<div class='stat-box'>\n";
    html << "<div class='stat-number'>" << drug_class_summaries.size() << "</div>\n";
    html << "<div class='stat-label'>Drug Classes Affected</div>\n";
    html << "</div>\n";
    
    html << "</div>\n";
    
    // Drug Class Analysis
    html << "<h2>Resistance by Drug Class</h2>\n";
    
    for (const auto& [drug_class, summary] : drug_class_summaries) {
        html << "<div class='summary'>\n";
        html << "<h3>" << drug_class << "</h3>\n";
        html << "<p><strong>Clinical Interpretation:</strong> " << summary.clinical_interpretation << "</p>\n";
        html << "<p><strong>Genes Detected:</strong> " << summary.genes_detected.size() 
             << " (High confidence: " << summary.high_confidence_genes.size() << ")</p>\n";
        html << "<p><strong>Total Reads:</strong> " << summary.total_reads << "</p>\n";
        html << "<p><strong>Maximum TPM:</strong> " << std::fixed << std::setprecision(2) 
             << summary.max_tpm << "</p>\n";
        html << "</div>\n";
    }
    
    // Gene Family Analysis
    html << "<h2>Gene Families Detected</h2>\n";
    html << "<div class='summary'>\n";
    html << "<p>Gene families group related resistance genes that may have similar functions or evolutionary origins.</p>\n";
    html << "</div>\n";
    
    // Sort gene families by total TPM
    std::vector<std::pair<std::string, GeneFamilySummary>> sorted_families(
        gene_family_summaries.begin(), gene_family_summaries.end()
    );
    std::sort(sorted_families.begin(), sorted_families.end(),
        [](const auto& a, const auto& b) {
            return a.second.total_tpm > b.second.total_tpm;
        }
    );
    
    html << "<table>\n";
    html << "<tr><th>Gene Family</th><th>Variants Detected</th><th>Total Reads</th>"
         << "<th>Total TPM</th><th>Mean Identity</th></tr>\n";
    
    for (const auto& [family_name, family_summary] : sorted_families) {
        html << "<tr>\n";
        html << "<td><strong>" << family_name << "</strong></td>\n";
        html << "<td>";
        for (size_t i = 0; i < family_summary.gene_variants.size(); i++) {
            if (i > 0) html << ", ";
            html << family_summary.gene_variants[i];
        }
        html << "</td>\n";
        html << "<td>" << static_cast<int>(family_summary.total_reads) << "</td>\n";
        html << "<td>" << std::fixed << std::setprecision(2) << family_summary.total_tpm << "</td>\n";
        html << "<td>" << std::fixed << std::setprecision(1) << (family_summary.mean_identity * 100) << "%</td>\n";
        html << "</tr>\n";
    }
    
    html << "</table>\n";
    
    // Detailed Gene Table
    html << "<h2>Detected Resistance Genes</h2>\n";
    html << "<table>\n";
    html << "<tr><th>Gene</th><th>Gene Family</th><th>Drug Class</th><th>Coverage (%)</th><th>Depth</th>"
         << "<th>TPM</th><th>Confidence</th><th>Complete Gene</th></tr>\n";
    
    // Sort genes by TPM for display
    std::vector<std::pair<std::string, GeneSummary>> sorted_genes(
        gene_summaries.begin(), gene_summaries.end()
    );
    std::sort(sorted_genes.begin(), sorted_genes.end(),
        [](const auto& a, const auto& b) {
            return a.second.tpm > b.second.tpm;
        }
    );
    
    for (const auto& [gene_name, summary] : sorted_genes) {
        std::string row_class = "";
        if (summary.confidence_level == "HIGH") row_class = "high-confidence";
        else if (summary.confidence_level == "MODERATE") row_class = "moderate-confidence";
        else row_class = "low-confidence";
        
        html << "<tr class='" << row_class << "'>\n";
        html << "<td><strong>" << gene_name << "</strong></td>\n";
        html << "<td>" << summary.gene_family << "</td>\n";
        html << "<td>" << summary.drug_class << "</td>\n";
        html << "<td>" << std::fixed << std::setprecision(1) << summary.percent_coverage << "</td>\n";
        html << "<td>" << std::fixed << std::setprecision(1) << summary.mean_depth << "</td>\n";
        html << "<td>" << std::fixed << std::setprecision(2) << summary.tpm << "</td>\n";
        html << "<td>" << summary.confidence_level << "</td>\n";
        html << "<td>" << (summary.is_complete_gene ? "Yes" : "No") << "</td>\n";
        html << "</tr>\n";
    }
    
    html << "</table>\n";
    
    // Clinical Recommendations
    html << "<h2>Clinical Recommendations</h2>\n";
    html << "<div class='summary'>\n";
    html << "<ul>\n";
    
    if (high_confidence_genes > 0) {
        html << "<li>High-confidence resistance genes detected. Avoid antibiotics from affected drug classes.</li>\n";
    }
    
    if (!drug_class_summaries.empty()) {
        html << "<li>Consider antibiotic susceptibility testing to confirm resistance profile.</li>\n";
        html << "<li>Consult local antibiogram data for empiric therapy guidance.</li>\n";
    }
    
    if (drug_class_summaries.find("GLYCOPEPTIDE") != drug_class_summaries.end() ||
        drug_class_summaries.find("COLISTIN") != drug_class_summaries.end()) {
        html << "<li><strong>Critical resistance detected - infectious disease consultation recommended.</strong></li>\n";
    }
    
    html << "</ul>\n";
    html << "</div>\n";
    
    // Footer
    html << "<div class='summary' style='margin-top: 40px; text-align: center; font-size: 0.9em;'>\n";
    html << "<p>This report is for research use only. Clinical decisions should be based on "
         << "validated diagnostic methods and professional medical judgment.</p>\n";
    html << "</div>\n";
    
    html << "</div>\n"; // container
    html << "</body>\n</html>\n";
    
    html.close();
}

void ClinicalAMRReportGenerator::generateTextReport() {
    std::ofstream txt(output_path + "_clinical_amr_report.txt");
    
    txt << "CLINICAL ANTIMICROBIAL RESISTANCE REPORT\n";
    txt << "=====================================\n\n";
    txt << "Sample: " << sample_name << "\n";
    txt << "Report Generated: " << getCurrentTimestamp() << "\n\n";
    
    txt << "SUMMARY\n";
    txt << "-------\n";
    txt << "Total reads processed: " << total_reads_processed << "\n";
    txt << "Reads with AMR genes: " << reads_with_amr_hits << "\n";
    txt << "Total AMR genes detected: " << total_genes_detected << "\n";
    txt << "High confidence genes: " << high_confidence_genes << "\n";
    txt << "Drug classes affected: " << drug_class_summaries.size() << "\n\n";
    
    txt << "RESISTANCE BY DRUG CLASS\n";
    txt << "------------------------\n";
    
    for (const auto& [drug_class, summary] : drug_class_summaries) {
        txt << "\n" << drug_class << ":\n";
        txt << "  Genes detected: " << summary.genes_detected.size() << "\n";
        txt << "  High confidence: " << summary.high_confidence_genes.size() << "\n";
        txt << "  Clinical interpretation: " << summary.clinical_interpretation << "\n";
    }
    
    txt << "\nDETAILED GENE LIST\n";
    txt << "------------------\n";
    txt << std::left << std::setw(30) << "Gene" 
        << std::setw(20) << "Drug Class"
        << std::setw(12) << "Coverage(%)"
        << std::setw(10) << "Depth"
        << std::setw(10) << "TPM"
        << "Confidence\n";
    txt << std::string(92, '-') << "\n";
    
    // Sort by TPM
    std::vector<std::pair<std::string, GeneSummary>> sorted_genes(
        gene_summaries.begin(), gene_summaries.end()
    );
    std::sort(sorted_genes.begin(), sorted_genes.end(),
        [](const auto& a, const auto& b) {
            return a.second.tpm > b.second.tpm;
        }
    );
    
    for (const auto& [gene_name, summary] : sorted_genes) {
        txt << std::left << std::setw(30) << gene_name
            << std::setw(20) << summary.drug_class
            << std::setw(12) << std::fixed << std::setprecision(1) << summary.percent_coverage
            << std::setw(10) << std::fixed << std::setprecision(1) << summary.mean_depth
            << std::setw(10) << std::fixed << std::setprecision(2) << summary.tpm
            << summary.confidence_level << "\n";
    }
    
    txt.close();
}

void ClinicalAMRReportGenerator::generateJSONReport() {
    std::ofstream json(output_path + "_clinical_amr_report.json");
    
    json << "{\n";
    json << "  \"sample_name\": \"" << sample_name << "\",\n";
    json << "  \"report_timestamp\": \"" << getCurrentTimestamp() << "\",\n";
    json << "  \"summary\": {\n";
    json << "    \"total_reads\": " << total_reads_processed << ",\n";
    json << "    \"reads_with_amr\": " << reads_with_amr_hits << ",\n";
    json << "    \"total_genes_detected\": " << total_genes_detected << ",\n";
    json << "    \"high_confidence_genes\": " << high_confidence_genes << ",\n";
    json << "    \"drug_classes_affected\": " << drug_class_summaries.size() << "\n";
    json << "  },\n";
    
    json << "  \"drug_class_analysis\": [\n";
    bool first_class = true;
    for (const auto& [drug_class, summary] : drug_class_summaries) {
        if (!first_class) json << ",\n";
        json << "    {\n";
        json << "      \"drug_class\": \"" << drug_class << "\",\n";
        json << "      \"total_genes\": " << summary.genes_detected.size() << ",\n";
        json << "      \"high_confidence_genes\": " << summary.high_confidence_genes.size() << ",\n";
        json << "      \"total_reads\": " << summary.total_reads << ",\n";
        json << "      \"max_tpm\": " << summary.max_tpm << ",\n";
        json << "      \"clinical_interpretation\": \"" << summary.clinical_interpretation << "\"\n";
        json << "    }";
        first_class = false;
    }
    json << "\n  ],\n";
    
    json << "  \"genes_detected\": [\n";
    bool first_gene = true;
    for (const auto& [gene_name, summary] : gene_summaries) {
        if (!first_gene) json << ",\n";
        json << "    {\n";
        json << "      \"gene_name\": \"" << gene_name << "\",\n";
        json << "      \"drug_class\": \"" << summary.drug_class << "\",\n";
        json << "      \"read_count\": " << summary.read_count << ",\n";
        json << "      \"percent_coverage\": " << summary.percent_coverage << ",\n";
        json << "      \"mean_depth\": " << summary.mean_depth << ",\n";
        json << "      \"tpm\": " << summary.tpm << ",\n";
        json << "      \"rpkm\": " << summary.rpkm << ",\n";
        json << "      \"confidence_level\": \"" << summary.confidence_level << "\",\n";
        json << "      \"is_complete_gene\": " << (summary.is_complete_gene ? "true" : "false") << "\n";
        json << "    }";
        first_gene = false;
    }
    json << "\n  ]\n";
    
    json << "}\n";
    json.close();
}

void ClinicalAMRReportGenerator::generateTSVReports() {
    // Gene abundance table
    std::ofstream abundance(output_path + "_amr_abundance.tsv");
    abundance << "gene_name\tgene_family\tdrug_class\tread_count\tpercent_coverage\tmean_depth\ttpkm\trpm\tconfidence\n";
    
    for (const auto& [gene_name, summary] : gene_summaries) {
        abundance << gene_name << "\t"
                 << summary.gene_family << "\t"
                 << summary.drug_class << "\t"
                 << summary.read_count << "\t"
                 << summary.percent_coverage << "\t"
                 << summary.mean_depth << "\t"
                 << summary.tpm << "\t"
                 << summary.rpkm << "\t"
                 << summary.confidence_level << "\n";
    }
    abundance.close();
    
    // Drug class summary
    std::ofstream drug_summary(output_path + "_drug_class_summary.tsv");
    drug_summary << "drug_class\ttotal_genes\thigh_confidence\tmoderate_confidence\tlow_confidence\ttotal_reads\tmax_tpm\n";
    
    for (const auto& [drug_class, summary] : drug_class_summaries) {
        drug_summary << drug_class << "\t"
                    << summary.genes_detected.size() << "\t"
                    << summary.high_confidence_genes.size() << "\t"
                    << summary.moderate_confidence_genes.size() << "\t"
                    << summary.low_confidence_genes.size() << "\t"
                    << summary.total_reads << "\t"
                    << summary.max_tpm << "\n";
    }
    drug_summary.close();
    
    // Gene family summary
    std::ofstream family_summary(output_path + "_gene_family_summary.tsv");
    family_summary << "gene_family\tnum_variants\tvariant_names\ttotal_reads\ttotal_tpm\tmean_identity\n";
    
    for (const auto& [family_name, summary] : gene_family_summaries) {
        family_summary << family_name << "\t"
                      << summary.gene_variants.size() << "\t";
        
        // Join variant names with semicolons
        for (size_t i = 0; i < summary.gene_variants.size(); i++) {
            if (i > 0) family_summary << ";";
            family_summary << summary.gene_variants[i];
        }
        
        family_summary << "\t"
                      << static_cast<int>(summary.total_reads) << "\t"
                      << summary.total_tpm << "\t"
                      << summary.mean_identity << "\n";
    }
    family_summary.close();
}

std::string ClinicalAMRReportGenerator::getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}
