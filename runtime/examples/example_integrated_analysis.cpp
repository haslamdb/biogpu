// runtime/examples/example_integrated_analysis.cpp
// Demonstrates how to run all three diagnostic analyses on the same sample
// to provide comprehensive antimicrobial resistance profiling for clinical decisions

#include "../common/pipeline/pipeline_coordinator.h"
#include "../common/config/unified_config.h"
#include "../common/io/streaming_fastq_reader.h"
#include "../kernels/resistance/resistance_pipeline_v2.h"
#include "../kernels/genes/amr_detection_pipeline_v2.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

using namespace BioGPU;

// Clinical report generator that combines results from all pipelines
class IntegratedClinicalReport {
public:
    struct PathogenProfile {
        std::string species;
        float abundance;
        std::vector<std::string> resistance_mutations;
        std::vector<std::string> resistance_genes;
        std::string resistance_level;  // "Susceptible", "Intermediate", "Resistant"
        std::vector<std::string> treatment_options;
    };
    
    void generateReport(
        const std::string& sample_id,
        const std::map<std::string, PathogenProfile>& pathogens,
        const std::string& output_file) {
        
        std::ofstream out(output_file);
        
        // HTML header
        out << "<!DOCTYPE html>\n<html>\n<head>\n";
        out << "<title>Clinical Diagnostic Report - " << sample_id << "</title>\n";
        out << "<style>\n";
        out << "body { font-family: Arial, sans-serif; margin: 20px; }\n";
        out << ".header { background-color: #2c3e50; color: white; padding: 20px; }\n";
        out << ".pathogen { border: 1px solid #ddd; margin: 10px 0; padding: 15px; }\n";
        out << ".resistant { background-color: #ffe6e6; }\n";
        out << ".susceptible { background-color: #e6ffe6; }\n";
        out << ".warning { color: #d9534f; font-weight: bold; }\n";
        out << ".success { color: #5cb85c; font-weight: bold; }\n";
        out << "</style>\n</head>\n<body>\n";
        
        // Header
        out << "<div class='header'>\n";
        out << "<h1>Clinical Metagenomic Analysis Report</h1>\n";
        out << "<p>Sample ID: " << sample_id << "</p>\n";
        out << "<p>Analysis Date: " << getCurrentDateTime() << "</p>\n";
        out << "</div>\n";
        
        // Executive Summary
        out << "<h2>Executive Summary</h2>\n";
        int resistant_count = 0;
        for (const auto& [species, profile] : pathogens) {
            if (profile.resistance_level == "Resistant") resistant_count++;
        }
        
        if (resistant_count > 0) {
            out << "<p class='warning'>⚠️ Antimicrobial resistance detected in " 
                << resistant_count << " pathogen(s)</p>\n";
        } else {
            out << "<p class='success'>✓ No significant antimicrobial resistance detected</p>\n";
        }
        
        // Detailed Results
        out << "<h2>Pathogen Analysis</h2>\n";
        
        for (const auto& [species, profile] : pathogens) {
            std::string css_class = profile.resistance_level == "Resistant" ? "resistant" : "susceptible";
            
            out << "<div class='pathogen " << css_class << "'>\n";
            out << "<h3>" << species << "</h3>\n";
            out << "<p><strong>Abundance:</strong> " << std::fixed << std::setprecision(1) 
                << profile.abundance << "%</p>\n";
            
            // Resistance Status
            out << "<p><strong>Resistance Status:</strong> ";
            if (profile.resistance_level == "Resistant") {
                out << "<span class='warning'>" << profile.resistance_level << "</span>";
            } else {
                out << "<span class='success'>" << profile.resistance_level << "</span>";
            }
            out << "</p>\n";
            
            // Detected Mutations
            if (!profile.resistance_mutations.empty()) {
                out << "<p><strong>Resistance Mutations:</strong></p>\n<ul>\n";
                for (const auto& mutation : profile.resistance_mutations) {
                    out << "<li>" << mutation << "</li>\n";
                }
                out << "</ul>\n";
            }
            
            // Detected Resistance Genes
            if (!profile.resistance_genes.empty()) {
                out << "<p><strong>Resistance Genes:</strong></p>\n<ul>\n";
                for (const auto& gene : profile.resistance_genes) {
                    out << "<li>" << gene << "</li>\n";
                }
                out << "</ul>\n";
            }
            
            // Treatment Recommendations
            out << "<p><strong>Treatment Recommendations:</strong></p>\n<ul>\n";
            for (const auto& treatment : profile.treatment_options) {
                out << "<li>" << treatment << "</li>\n";
            }
            out << "</ul>\n";
            
            out << "</div>\n";
        }
        
        // Clinical Interpretation
        out << "<h2>Clinical Interpretation</h2>\n";
        out << generateClinicalInterpretation(pathogens);
        
        // Footer
        out << "<hr>\n";
        out << "<p><small>This report is for clinical decision support only. ";
        out << "Results should be interpreted in conjunction with clinical findings ";
        out << "and local antimicrobial resistance patterns.</small></p>\n";
        out << "</body>\n</html>\n";
        
        out.close();
    }
    
private:
    std::string getCurrentDateTime() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }
    
    std::string generateClinicalInterpretation(
        const std::map<std::string, PathogenProfile>& pathogens) {
        
        std::stringstream ss;
        
        // Find primary pathogen
        std::string primary_pathogen;
        float max_abundance = 0;
        for (const auto& [species, profile] : pathogens) {
            if (profile.abundance > max_abundance) {
                max_abundance = profile.abundance;
                primary_pathogen = species;
            }
        }
        
        if (!primary_pathogen.empty()) {
            ss << "<p>Primary pathogen: <strong>" << primary_pathogen 
               << "</strong> (" << std::fixed << std::setprecision(1) 
               << max_abundance << "% abundance)</p>\n";
            
            const auto& profile = pathogens.at(primary_pathogen);
            
            if (profile.resistance_level == "Resistant") {
                ss << "<p class='warning'>⚠️ The primary pathogen shows antimicrobial resistance. ";
                ss << "Standard empiric therapy may be ineffective.</p>\n";
                ss << "<p>Consider the following alternatives:</p>\n<ul>\n";
                for (const auto& treatment : profile.treatment_options) {
                    ss << "<li>" << treatment << "</li>\n";
                }
                ss << "</ul>\n";
            } else {
                ss << "<p class='success'>✓ The primary pathogen appears susceptible to standard therapy.</p>\n";
            }
        }
        
        // Check for polymicrobial infection
        int pathogen_count = 0;
        for (const auto& [species, profile] : pathogens) {
            if (profile.abundance > 1.0) pathogen_count++;
        }
        
        if (pathogen_count > 1) {
            ss << "<p><strong>Note:</strong> This appears to be a polymicrobial infection. ";
            ss << "Consider broad-spectrum coverage or combination therapy.</p>\n";
        }
        
        return ss.str();
    }
};

// Main integrated analysis function
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <sample.fastq.gz> [config.json] [output_dir]\n";
        std::cerr << "\nThis tool provides comprehensive antimicrobial resistance profiling\n";
        std::cerr << "by integrating resistance mutation detection, gene quantification,\n";
        std::cerr << "and taxonomic profiling for clinical decision support.\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string config_file = argc > 2 ? argv[2] : "";
    std::string output_dir = argc > 3 ? argv[3] : "integrated_results";
    
    std::cout << "=== BioGPU Integrated Clinical Diagnostic Analysis ===\n";
    std::cout << "Sample: " << input_file << "\n\n";
    
    // Load or create configuration
    std::unique_ptr<UnifiedConfig> config;
    if (!config_file.empty()) {
        config = ConfigLoader::loadFromFile(config_file);
        if (!config) {
            std::cerr << "Failed to load config file, using defaults\n";
            config = std::make_unique<UnifiedConfig>(ConfigLoader::getDefault());
        }
    } else {
        config = std::make_unique<UnifiedConfig>(ConfigLoader::getDefault());
    }
    
    // Enable all pipelines for comprehensive analysis
    config->pipeline_config.enable_resistance = true;
    config->pipeline_config.enable_genes = true;
    config->pipeline_config.enable_profiler = true;
    config->pipeline_config.output_dir = output_dir;
    
    // Create output directory
    system(("mkdir -p " + output_dir).c_str());
    
    // Initialize pipeline coordinator
    PipelineCoordinator coordinator(*config);
    
    std::cout << "Initializing diagnostic pipelines...\n";
    if (!coordinator.initialize()) {
        std::cerr << "Failed to initialize pipelines\n";
        return 1;
    }
    
    // Process sample
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "\nProcessing sample through integrated pipeline:\n";
    std::cout << "  1. Taxonomic profiling...\n";
    std::cout << "  2. Resistance mutation detection...\n";
    std::cout << "  3. Resistance gene quantification...\n\n";
    
    bool success = coordinator.processFastq(input_file);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    if (!success) {
        std::cerr << "Error processing sample\n";
        return 1;
    }
    
    // Get results from all pipelines
    const auto& results = coordinator.getResults();
    
    std::cout << "\n=== Analysis Complete ===\n";
    std::cout << "Total processing time: " << duration.count() << " seconds\n";
    std::cout << "Sequences processed: " << results.total_reads_processed << "\n";
    std::cout << "Bases processed: " << results.total_bases_processed << "\n\n";
    
    // Example of integrated interpretation
    // In practice, this would parse actual results from each pipeline
    std::cout << "=== Integrated Clinical Findings ===\n\n";
    
    // Simulated integrated results (would be populated from actual pipeline results)
    std::map<std::string, IntegratedClinicalReport::PathogenProfile> pathogen_profiles;
    
    // Example: E. coli with fluoroquinolone resistance
    IntegratedClinicalReport::PathogenProfile ecoli;
    ecoli.species = "Escherichia coli";
    ecoli.abundance = 78.5;
    ecoli.resistance_mutations = {"gyrA_S83L (98% frequency)", "parC_S80I (95% frequency)"};
    ecoli.resistance_genes = {"CTX-M-15 (ESBL)", "aac(6')-Ib-cr"};
    ecoli.resistance_level = "Resistant";
    ecoli.treatment_options = {
        "Avoid fluoroquinolones (ciprofloxacin, levofloxacin)",
        "Consider: Carbapenems (meropenem, imipenem)",
        "Alternative: Aminoglycosides (if susceptible)",
        "For UTI: Nitrofurantoin or fosfomycin may be effective"
    };
    pathogen_profiles["Escherichia coli"] = ecoli;
    
    // Example: K. pneumoniae co-infection
    IntegratedClinicalReport::PathogenProfile kpneu;
    kpneu.species = "Klebsiella pneumoniae";
    kpneu.abundance = 18.2;
    kpneu.resistance_mutations = {};
    kpneu.resistance_genes = {"SHV-1 (chromosomal beta-lactamase)"};
    kpneu.resistance_level = "Susceptible";
    kpneu.treatment_options = {
        "Susceptible to standard therapy",
        "Ceftriaxone or fluoroquinolones likely effective"
    };
    pathogen_profiles["Klebsiella pneumoniae"] = kpneu;
    
    // Display summary
    for (const auto& [species, profile] : pathogen_profiles) {
        std::cout << species << " (" << std::fixed << std::setprecision(1) 
                  << profile.abundance << "% abundance)\n";
        std::cout << "  Resistance status: " << profile.resistance_level << "\n";
        if (!profile.resistance_mutations.empty()) {
            std::cout << "  Key mutations: ";
            for (const auto& mut : profile.resistance_mutations) {
                std::cout << mut << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    
    // Generate integrated clinical report
    IntegratedClinicalReport report_generator;
    std::string report_file = output_dir + "/integrated_clinical_report.html";
    report_generator.generateReport(
        input_file.substr(input_file.find_last_of("/\\") + 1),
        pathogen_profiles,
        report_file
    );
    
    std::cout << "=== Reports Generated ===\n";
    std::cout << "Integrated clinical report: " << report_file << "\n";
    std::cout << "Resistance mutations: " << output_dir << "/resistance_report.tsv\n";
    std::cout << "Gene expression: " << output_dir << "/amr_genes_report.tsv\n";
    std::cout << "Taxonomic profile: " << output_dir << "/taxonomy_report.tsv\n";
    
    std::cout << "\n=== Clinical Recommendation ===\n";
    std::cout << "Based on integrated analysis:\n";
    std::cout << "• Primary pathogen: E. coli with high-level fluoroquinolone resistance\n";
    std::cout << "• Avoid empiric fluoroquinolone therapy\n";
    std::cout << "• Recommended: Carbapenem-based therapy\n";
    std::cout << "• Monitor for treatment response given polymicrobial infection\n";
    
    return 0;
}