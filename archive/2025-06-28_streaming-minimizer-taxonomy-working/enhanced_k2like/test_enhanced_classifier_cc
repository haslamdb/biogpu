// test_enhanced_classifier.cu
// Test program for the Phase 1 Enhanced Classifier

#include "phase1_enhanced_classification.cu"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

// Helper function to load reads from FASTQ file
std::vector<std::string> load_fastq_reads(const std::string& filename, int max_reads = -1) {
    std::vector<std::string> reads;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return reads;
    }
    
    std::string line;
    int count = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '@') {  // Header line
            // Get sequence line
            if (std::getline(file, line)) {
                reads.push_back(line);
                count++;
                
                // Skip quality header and quality line
                std::getline(file, line);
                std::getline(file, line);
                
                if (max_reads > 0 && count >= max_reads) break;
            }
        }
    }
    
    file.close();
    return reads;
}

// Helper function to print classification results
void print_enhanced_results(const std::vector<Phase1ClassificationResult>& results, int limit = 10) {
    std::cout << "\n=== SAMPLE ENHANCED CLASSIFICATION RESULTS ===" << std::endl;
    std::cout << "Showing first " << std::min(limit, (int)results.size()) << " results:" << std::endl;
    
    for (int i = 0; i < std::min(limit, (int)results.size()); i++) {
        const auto& r = results[i];
        std::cout << "\nRead " << i << ":" << std::endl;
        std::cout << "  Taxon ID: " << r.taxon_id << std::endl;
        std::cout << "  Stage 1: " << (r.passed_stage1 ? "PASS" : "FAIL") 
                  << " (conf: " << r.primary_confidence << ")" << std::endl;
        std::cout << "  Stage 2: " << (r.passed_stage2 ? "PASS" : "FAIL")
                  << " (conf: " << r.secondary_confidence << ")" << std::endl;
        std::cout << "  High confidence: " << (r.is_high_confidence_call ? "YES" : "NO") << std::endl;
        std::cout << "  K-mers: " << r.classified_kmers << "/" << r.total_kmers << std::endl;
        std::cout << "  Distinct minimizers: " << r.distinct_minimizers_found << std::endl;
        std::cout << "  Phylogenetic score: " << r.phylogenetic_consistency_score << std::endl;
        std::cout << "  Has conflicts: " << (r.has_taxonomic_conflicts ? "YES" : "NO") << std::endl;
        
        if (r.kmer_tracking.num_tracked_kmers > 0) {
            std::cout << "  Tracked k-mers: " << r.kmer_tracking.num_tracked_kmers 
                      << (r.kmer_tracking.tracking_overflow ? " (overflow)" : "") << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <database_dir> <reads.fastq> [options]" << std::endl;
        std::cout << "\nOptions:" << std::endl;
        std::cout << "  --confidence <float>     Primary confidence threshold (default: 0.1)" << std::endl;
        std::cout << "  --confidence2 <float>    Secondary confidence threshold (default: 0.3)" << std::endl;
        std::cout << "  --min-kmers <int>        Minimum k-mers for classification (default: 10)" << std::endl;
        std::cout << "  --min-minimizers <int>   Minimum distinct minimizers (default: 3)" << std::endl;
        std::cout << "  --no-phylo               Disable phylogenetic validation" << std::endl;
        std::cout << "  --no-tracking            Disable k-mer position tracking" << std::endl;
        std::cout << "  --max-reads <int>        Maximum reads to process (default: all)" << std::endl;
        std::cout << "  --compare                Compare with base classifier" << std::endl;
        return 1;
    }
    
    std::string database_dir = argv[1];
    std::string reads_file = argv[2];
    
    // Parse command line options
    Phase1EnhancedParams params;
    int max_reads_to_process = -1;
    bool compare_with_base = false;
    
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--confidence" && i + 1 < argc) {
            params.primary_confidence_threshold = std::stof(argv[++i]);
        } else if (arg == "--confidence2" && i + 1 < argc) {
            params.secondary_confidence_threshold = std::stof(argv[++i]);
        } else if (arg == "--min-kmers" && i + 1 < argc) {
            params.min_kmers_for_classification = std::stoi(argv[++i]);
        } else if (arg == "--min-minimizers" && i + 1 < argc) {
            params.min_distinct_minimizers = std::stoi(argv[++i]);
        } else if (arg == "--no-phylo") {
            params.enable_phylogenetic_validation = false;
        } else if (arg == "--no-tracking") {
            params.forward_compat.track_kmer_positions = false;
        } else if (arg == "--max-reads" && i + 1 < argc) {
            max_reads_to_process = std::stoi(argv[++i]);
        } else if (arg == "--compare") {
            compare_with_base = true;
        }
    }
    
    std::cout << "Enhanced Kraken2 Classifier Test" << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << "Database: " << database_dir << std::endl;
    std::cout << "Reads: " << reads_file << std::endl;
    std::cout << "\nEnhanced Parameters:" << std::endl;
    std::cout << "  Primary confidence: " << params.primary_confidence_threshold << std::endl;
    std::cout << "  Secondary confidence: " << params.secondary_confidence_threshold << std::endl;
    std::cout << "  Min k-mers: " << params.min_kmers_for_classification << std::endl;
    std::cout << "  Min minimizers: " << params.min_distinct_minimizers << std::endl;
    std::cout << "  Phylogenetic validation: " << (params.enable_phylogenetic_validation ? "ON" : "OFF") << std::endl;
    std::cout << "  K-mer tracking: " << (params.forward_compat.track_kmer_positions ? "ON" : "OFF") << std::endl;
    
    // Load reads
    std::cout << "\nLoading reads..." << std::endl;
    auto reads = load_fastq_reads(reads_file, max_reads_to_process);
    if (reads.empty()) {
        std::cerr << "Error: No reads loaded" << std::endl;
        return 1;
    }
    std::cout << "Loaded " << reads.size() << " reads" << std::endl;
    
    // Create enhanced classifier
    std::cout << "\nInitializing enhanced classifier..." << std::endl;
    Phase1EnhancedClassifier classifier(params);
    
    // Load database
    std::cout << "Loading database..." << std::endl;
    if (!classifier.load_database(database_dir)) {
        std::cerr << "Error: Failed to load database" << std::endl;
        return 1;
    }
    
    // Perform enhanced classification
    std::cout << "\nPerforming enhanced classification..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    auto enhanced_results = classifier.classify_enhanced(reads);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "Enhanced classification completed in " << duration.count() << " ms" << std::endl;
    std::cout << "Speed: " << (reads.size() * 1000.0 / duration.count()) << " reads/second" << std::endl;
    
    // Print statistics
    classifier.print_enhanced_statistics(enhanced_results);
    
    // Print sample results
    print_enhanced_results(enhanced_results, 5);
    
    // Optional: Compare with base classifier
    if (compare_with_base) {
        std::cout << "\n\n=== COMPARISON WITH BASE CLASSIFIER ===" << std::endl;
        
        // Create base classifier
        ClassificationParams base_params;
        base_params.confidence_threshold = params.primary_confidence_threshold;
        PairedEndGPUKrakenClassifier base_classifier(base_params);
        
        if (!base_classifier.load_database(database_dir)) {
            std::cerr << "Error: Failed to load database for base classifier" << std::endl;
            return 1;
        }
        
        // Convert reads to paired format (single-end)
        std::vector<PairedRead> paired_reads;
        for (size_t i = 0; i < reads.size(); i++) {
            paired_reads.emplace_back(reads[i], "", std::to_string(i));
        }
        
        std::cout << "Running base classifier..." << std::endl;
        start_time = std::chrono::high_resolution_clock::now();
        
        auto base_results = base_classifier.classify_paired_reads(paired_reads);
        
        end_time = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Base classification completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Speed: " << (reads.size() * 1000.0 / duration.count()) << " reads/second" << std::endl;
        
        // Compare results
        int agreement_count = 0;
        int enhanced_only = 0;
        int base_only = 0;
        int both_unclassified = 0;
        
        for (size_t i = 0; i < enhanced_results.size() && i < base_results.size(); i++) {
            bool enhanced_classified = enhanced_results[i].taxon_id != 0;
            bool base_classified = base_results[i].taxon_id != 0;
            
            if (enhanced_results[i].taxon_id == base_results[i].taxon_id) {
                agreement_count++;
                if (!enhanced_classified) both_unclassified++;
            } else if (enhanced_classified && !base_classified) {
                enhanced_only++;
            } else if (!enhanced_classified && base_classified) {
                base_only++;
            }
        }
        
        std::cout << "\nComparison Results:" << std::endl;
        std::cout << "  Agreement: " << agreement_count << "/" << enhanced_results.size() 
                  << " (" << (100.0 * agreement_count / enhanced_results.size()) << "%)" << std::endl;
        std::cout << "  Enhanced only: " << enhanced_only << std::endl;
        std::cout << "  Base only: " << base_only << std::endl;
        std::cout << "  Both unclassified: " << both_unclassified << std::endl;
        
        // Print classification statistics
        base_classifier.print_paired_classification_stats(base_results);
    }
    
    std::cout << "\nTest completed successfully!" << std::endl;
    
    return 0;
}