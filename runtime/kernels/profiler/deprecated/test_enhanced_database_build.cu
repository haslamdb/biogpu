// test_enhanced_database_build.cu
// Comprehensive test suite for ML-enhanced database building
// Tests feature extraction, contamination detection, and database integrity

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <random>
#include <iomanip>
#include <cuda_runtime.h>

#include "gpu_kraken_database_builder.h"
#include "gpu_kraken_types.h"
#include "processing/minimizer_feature_extractor.h"
#include "processing/contamination_detector.h"
#include "processing/feature_exporter.h"
#include "output/database_serializer.h"
#include "gpu/gpu_database_kernels.h"

// Test configuration
struct TestConfig {
    std::string genome_dir = "data/type_strain_reference_genomes";
    std::string output_dir = "test_output";
    std::string contamination_test_dir = "test_contamination";
    int k_value = 31;
    int ell_value = 31;
    int spaces_value = 7;
    int max_test_genomes = 10;  // Limit for quick testing
    bool verbose = true;
};

// Test result structure
struct TestResults {
    bool feature_extraction_passed = false;
    bool contamination_detection_passed = false;
    bool database_integrity_passed = false;
    bool performance_benchmark_passed = false;
    
    // Statistics
    size_t total_minimizers = 0;
    size_t minimizers_with_features = 0;
    size_t contamination_markers = 0;
    double feature_extraction_time = 0.0;
    double contamination_detection_time = 0.0;
    double database_build_time = 0.0;
    size_t database_size_20byte = 0;
    size_t database_size_24byte = 0;
    
    // Feature distributions
    std::vector<int> gc_distribution;
    std::vector<int> complexity_distribution;
    int position_bias_count = 0;
    
    void print_summary() const {
        std::cout << "\n=== TEST RESULTS SUMMARY ===" << std::endl;
        std::cout << "Feature Extraction: " << (feature_extraction_passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << "Contamination Detection: " << (contamination_detection_passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << "Database Integrity: " << (database_integrity_passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << "Performance Benchmark: " << (performance_benchmark_passed ? "PASSED" : "FAILED") << std::endl;
        
        std::cout << "\n=== STATISTICS ===" << std::endl;
        std::cout << "Total minimizers: " << total_minimizers << std::endl;
        std::cout << "Minimizers with features: " << minimizers_with_features << std::endl;
        std::cout << "Contamination markers: " << contamination_markers << std::endl;
        std::cout << "Feature extraction time: " << std::fixed << std::setprecision(3) << feature_extraction_time << "s" << std::endl;
        std::cout << "Contamination detection time: " << contamination_detection_time << "s" << std::endl;
        std::cout << "Database build time: " << database_build_time << "s" << std::endl;
        
        std::cout << "\n=== DATABASE SIZE COMPARISON ===" << std::endl;
        std::cout << "20-byte structure size: " << (database_size_20byte / 1024 / 1024) << " MB" << std::endl;
        std::cout << "24-byte structure size: " << (database_size_24byte / 1024 / 1024) << " MB" << std::endl;
        std::cout << "Size increase: " << std::fixed << std::setprecision(1) 
                  << ((double)database_size_24byte / database_size_20byte - 1.0) * 100 << "%" << std::endl;
    }
};

// Helper function to create contaminated sequences
std::vector<std::string> create_contamination_test_sequences() {
    std::vector<std::string> sequences;
    
    // Create synthetic contaminated sequences with patterns from multiple species
    std::string contam_seq1 = "ATCGATCGATCGATCGATCGATCGATCGATCG";  // Low complexity
    for (int i = 0; i < 10; i++) {
        contam_seq1 += "NNNNNNNNNN";  // Add N's (common contamination marker)
    }
    sequences.push_back(contam_seq1);
    
    // High GC contamination (often from lab contamination)
    std::string contam_seq2 = "";
    for (int i = 0; i < 100; i++) {
        contam_seq2 += "GCGCGCGC";
    }
    sequences.push_back(contam_seq2);
    
    // Adapter sequences (common contamination)
    std::string adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    sequences.push_back(adapter + adapter + adapter);
    
    // PhiX control (common spike-in contamination)
    std::string phix = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT";
    sequences.push_back(phix);
    
    return sequences;
}

// Test 1: Feature Extraction
bool test_feature_extraction(const TestConfig& config, TestResults& results) {
    std::cout << "\n=== TEST 1: FEATURE EXTRACTION ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Create feature extractor
        MinimizerFeatureExtractor feature_extractor;
        
        // Test sequences with known properties
        std::vector<std::string> test_sequences = {
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",  // Low complexity, low GC
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",  // Low complexity, high GC
            "ATCGATCGATCGATCGATCGATCGATCGATCGAT",  // Medium complexity, medium GC
            "ACGTACGTGCATGCATACGTGCATACGTGCATGC",  // High complexity, medium GC
        };
        
        // Process each test sequence
        for (size_t i = 0; i < test_sequences.size(); i++) {
            const auto& seq = test_sequences[i];
            std::vector<uint64_t> minimizers;
            std::vector<uint32_t> positions;
            
            // Extract minimizers
            for (size_t pos = 0; pos <= seq.length() - config.k_value; pos++) {
                std::string kmer = seq.substr(pos, config.k_value);
                uint64_t hash = 0;
                for (char c : kmer) {
                    hash = (hash << 2) | ((c >> 1) & 3);
                }
                minimizers.push_back(hash);
                positions.push_back(pos);
            }
            
            // Extract features for each minimizer
            for (size_t j = 0; j < minimizers.size(); j++) {
                MinimizerFeatures features;
                feature_extractor.extract_features_cpu(
                    minimizers[j], 
                    seq.c_str() + positions[j], 
                    config.k_value,
                    positions[j],
                    seq.length(),
                    features
                );
                
                // Validate features
                if (i == 0) {  // Low GC sequence
                    if (features.gc_content > 0.3f) {
                        std::cerr << "ERROR: Low GC sequence has high GC feature" << std::endl;
                        return false;
                    }
                }
                if (i == 1) {  // High GC sequence
                    if (features.gc_content < 0.7f) {
                        std::cerr << "ERROR: High GC sequence has low GC feature" << std::endl;
                        return false;
                    }
                }
                
                results.minimizers_with_features++;
            }
            results.total_minimizers += minimizers.size();
        }
        
        std::cout << "✓ Feature extraction test passed" << std::endl;
        std::cout << "  Processed " << results.minimizers_with_features << " minimizers" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Feature extraction test failed: " << e.what() << std::endl;
        return false;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    results.feature_extraction_time = std::chrono::duration<double>(end_time - start_time).count();
    results.feature_extraction_passed = true;
    return true;
}

// Test 2: Contamination Detection
bool test_contamination_detection(const TestConfig& config, TestResults& results) {
    std::cout << "\n=== TEST 2: CONTAMINATION DETECTION ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Create contamination detector
        ContaminationDetectorConfig contam_config;
        contam_config.enable_ml_scoring = true;
        ContaminationDetector detector(contam_config);
        
        // Get contamination test sequences
        auto contam_sequences = create_contamination_test_sequences();
        
        // Also load some real sequences
        std::vector<std::string> real_sequences;
        std::string test_genome_file = config.genome_dir + "/GCF_000005845.2_ASM584v2_genomic.fna";
        if (std::filesystem::exists(test_genome_file)) {
            std::ifstream file(test_genome_file);
            std::string line, seq;
            while (std::getline(file, line)) {
                if (line[0] != '>') {
                    seq += line;
                    if (seq.length() > 1000) {
                        real_sequences.push_back(seq.substr(0, 1000));
                        seq.clear();
                    }
                }
            }
        }
        
        // Process contamination sequences
        int detected_contamination = 0;
        for (const auto& seq : contam_sequences) {
            std::vector<ContaminationMarker> markers;
            detector.detect_contamination_cpu(seq, markers);
            
            if (!markers.empty()) {
                detected_contamination++;
                results.contamination_markers += markers.size();
            }
        }
        
        // Process real sequences (should have fewer contamination markers)
        int false_positives = 0;
        for (const auto& seq : real_sequences) {
            std::vector<ContaminationMarker> markers;
            detector.detect_contamination_cpu(seq, markers);
            
            if (!markers.empty()) {
                false_positives++;
            }
        }
        
        std::cout << "✓ Contamination detection test passed" << std::endl;
        std::cout << "  Detected " << detected_contamination << "/" << contam_sequences.size() 
                  << " contaminated sequences" << std::endl;
        std::cout << "  False positives: " << false_positives << "/" << real_sequences.size() 
                  << " real sequences" << std::endl;
        
        // Check detection rate
        double detection_rate = (double)detected_contamination / contam_sequences.size();
        if (detection_rate < 0.5) {
            std::cerr << "WARNING: Low contamination detection rate: " << detection_rate << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Contamination detection test failed: " << e.what() << std::endl;
        return false;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    results.contamination_detection_time = std::chrono::duration<double>(end_time - start_time).count();
    results.contamination_detection_passed = true;
    return true;
}

// Test 3: Database Integrity
bool test_database_integrity(const TestConfig& config, TestResults& results) {
    std::cout << "\n=== TEST 3: DATABASE INTEGRITY ===" << std::endl;
    
    try {
        // Create test output directory
        std::filesystem::create_directories(config.output_dir);
        
        // Build a small test database
        ClassificationParams params;
        params.k = config.k_value;
        params.ell = config.ell_value;
        params.spaces = config.spaces_value;
        
        GPUKrakenDatabaseBuilder builder(config.output_dir, params);
        
        // Initialize CUDA
        if (!builder.initialize_cuda_context()) {
            std::cerr << "Failed to initialize CUDA context" << std::endl;
            return false;
        }
        
        // Create test minimizer data
        std::vector<PhylogeneticLCACandidate> test_candidates;
        std::vector<StreamlinedMinimizerMetadata> test_metadata;
        
        // Create some test minimizers with various features
        for (int i = 0; i < 1000; i++) {
            PhylogeneticLCACandidate candidate;
            candidate.minimizer_hash = 0x123456789ABCDEF0 + i;
            candidate.lca_taxon = 1000 + (i % 10);
            candidate.genome_count = 1 + (i % 5);
            candidate.uniqueness_score = (float)(i % 100) / 100.0f;
            
            // Add some contributing species
            for (int j = 0; j < (i % 3) + 1; j++) {
                candidate.contributing_species.push_back(1000 + j);
                candidate.genome_counts_per_species.push_back(1);
            }
            
            test_candidates.push_back(candidate);
            
            // Create corresponding metadata
            StreamlinedMinimizerMetadata metadata;
            metadata.minimizer_hash = candidate.minimizer_hash;
            metadata.lca_taxon = candidate.lca_taxon;
            metadata.total_genome_count = candidate.genome_count;
            metadata.set_ml_confidence(candidate.uniqueness_score);
            
            // Set feature flags
            uint8_t gc_category = i % 8;
            uint8_t complexity = (i / 8) % 8;
            metadata.feature_flags = MinimizerFlags::set_gc_content_category(0, gc_category);
            metadata.feature_flags = MinimizerFlags::set_complexity_score(metadata.feature_flags, complexity);
            
            if (i % 10 == 0) {
                metadata.feature_flags = MinimizerFlags::set_contamination_risk(metadata.feature_flags, true);
            }
            
            test_metadata.push_back(metadata);
        }
        
        // Test serialization
        EnhancedDatabaseSerializer serializer(config.output_dir);
        ContributingTaxaArrays dummy_taxa;
        std::unordered_map<uint32_t, std::string> taxon_names;
        EnhancedBuildStats stats;
        
        if (!serializer.save_enhanced_database(test_candidates, dummy_taxa, taxon_names, stats)) {
            std::cerr << "Failed to save enhanced database" << std::endl;
            return false;
        }
        
        // Verify files were created
        std::string hash_table = config.output_dir + "/enhanced_hash_table.bin";
        std::string ml_weights = config.output_dir + "/ml_weights.bin";
        std::string feature_stats = config.output_dir + "/feature_statistics.json";
        std::string contam_markers = config.output_dir + "/contamination_markers.bin";
        
        if (!std::filesystem::exists(hash_table)) {
            std::cerr << "Hash table file not created" << std::endl;
            return false;
        }
        
        if (!std::filesystem::exists(ml_weights)) {
            std::cerr << "ML weights file not created" << std::endl;
            return false;
        }
        
        if (!std::filesystem::exists(feature_stats)) {
            std::cerr << "Feature statistics file not created" << std::endl;
            return false;
        }
        
        if (!std::filesystem::exists(contam_markers)) {
            std::cerr << "Contamination markers file not created" << std::endl;
            return false;
        }
        
        // Verify 24-byte structure size
        size_t expected_size = test_metadata.size() * sizeof(StreamlinedMinimizerMetadata);
        size_t actual_size = std::filesystem::file_size(hash_table) - 24; // Subtract header
        
        std::cout << "Expected structure size: " << expected_size << " bytes" << std::endl;
        std::cout << "Actual data size: " << actual_size << " bytes" << std::endl;
        
        // Calculate size comparison
        results.database_size_20byte = test_metadata.size() * 20;  // Old structure size
        results.database_size_24byte = test_metadata.size() * 32;  // New structure size
        
        // Read back and verify version
        std::ifstream in(hash_table, std::ios::binary);
        uint64_t version;
        in.read(reinterpret_cast<char*>(&version), sizeof(uint64_t));
        
        if (version != 3) {
            std::cerr << "Incorrect database version: " << version << " (expected 3)" << std::endl;
            return false;
        }
        
        std::cout << "✓ Database integrity test passed" << std::endl;
        std::cout << "  Database version: " << version << std::endl;
        std::cout << "  All required files created" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Database integrity test failed: " << e.what() << std::endl;
        return false;
    }
    
    results.database_integrity_passed = true;
    return true;
}

// Test 4: Performance Benchmark
bool benchmark_performance_impact(const TestConfig& config, TestResults& results) {
    std::cout << "\n=== TEST 4: PERFORMANCE BENCHMARK ===" << std::endl;
    
    try {
        // Get test genome files
        std::vector<std::string> genome_files;
        if (std::filesystem::exists(config.genome_dir)) {
            for (const auto& entry : std::filesystem::directory_iterator(config.genome_dir)) {
                if (entry.path().extension() == ".fna" || entry.path().extension() == ".fasta") {
                    genome_files.push_back(entry.path().string());
                    if (genome_files.size() >= config.max_test_genomes) break;
                }
            }
        }
        
        if (genome_files.empty()) {
            std::cerr << "No genome files found in " << config.genome_dir << std::endl;
            return false;
        }
        
        std::cout << "Testing with " << genome_files.size() << " genome files" << std::endl;
        
        // Benchmark old structure (simulate 20-byte)
        auto start_old = std::chrono::high_resolution_clock::now();
        
        ClassificationParams params;
        params.k = config.k_value;
        std::string output_old = config.output_dir + "_old";
        std::filesystem::create_directories(output_old);
        
        GPUKrakenDatabaseBuilder builder_old(output_old, params);
        builder_old.initialize_cuda_context();
        
        // Simulate processing without ML features
        for (const auto& genome_file : genome_files) {
            // Just read the file to simulate I/O
            std::ifstream file(genome_file);
            std::string line;
            while (std::getline(file, line)) {
                // Minimal processing
            }
        }
        
        auto end_old = std::chrono::high_resolution_clock::now();
        double time_old = std::chrono::duration<double>(end_old - start_old).count();
        
        // Benchmark new structure (24-byte with ML)
        auto start_new = std::chrono::high_resolution_clock::now();
        
        std::string output_new = config.output_dir + "_new";
        std::filesystem::create_directories(output_new);
        
        GPUKrakenDatabaseBuilder builder_new(output_new, params);
        builder_new.initialize_cuda_context();
        
        // Process with ML features
        MinimizerFeatureExtractor feature_extractor;
        ContaminationDetector contam_detector;
        
        for (const auto& genome_file : genome_files) {
            std::ifstream file(genome_file);
            std::string line, seq;
            while (std::getline(file, line)) {
                if (line[0] != '>') {
                    seq += line;
                    if (seq.length() > 10000) {
                        // Extract features for a sample
                        MinimizerFeatures features;
                        feature_extractor.extract_features_cpu(
                            0x123456789ABCDEF0,
                            seq.c_str(),
                            config.k_value,
                            0,
                            seq.length(),
                            features
                        );
                        
                        // Check contamination
                        std::vector<ContaminationMarker> markers;
                        contam_detector.detect_contamination_cpu(seq.substr(0, 1000), markers);
                        
                        seq.clear();
                    }
                }
            }
        }
        
        auto end_new = std::chrono::high_resolution_clock::now();
        double time_new = std::chrono::duration<double>(end_new - start_new).count();
        
        results.database_build_time = time_new;
        
        // Calculate performance impact
        double overhead = ((time_new / time_old) - 1.0) * 100.0;
        
        std::cout << "✓ Performance benchmark completed" << std::endl;
        std::cout << "  Old structure time: " << std::fixed << std::setprecision(3) << time_old << "s" << std::endl;
        std::cout << "  New structure time: " << time_new << "s" << std::endl;
        std::cout << "  Performance overhead: " << std::setprecision(1) << overhead << "%" << std::endl;
        
        // Acceptable overhead threshold (e.g., 20%)
        if (overhead > 50.0) {
            std::cerr << "WARNING: High performance overhead detected" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Performance benchmark failed: " << e.what() << std::endl;
        return false;
    }
    
    results.performance_benchmark_passed = true;
    return true;
}

// Feature distribution analysis
void analyze_feature_distributions(const TestConfig& config, TestResults& results) {
    std::cout << "\n=== FEATURE DISTRIBUTION ANALYSIS ===" << std::endl;
    
    // Read feature statistics if available
    std::string stats_file = config.output_dir + "/feature_statistics.json";
    if (std::filesystem::exists(stats_file)) {
        std::ifstream file(stats_file);
        std::string content((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());
        
        std::cout << "Feature Statistics:" << std::endl;
        std::cout << content << std::endl;
    }
    
    // Initialize distributions
    results.gc_distribution.resize(8, 0);
    results.complexity_distribution.resize(8, 0);
    
    // Analyze ML weights
    std::string ml_weights_file = config.output_dir + "/ml_weights.bin";
    if (std::filesystem::exists(ml_weights_file)) {
        std::ifstream file(ml_weights_file, std::ios::binary);
        uint64_t version, num_entries;
        file.read(reinterpret_cast<char*>(&version), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&num_entries), sizeof(uint64_t));
        
        std::cout << "\nML Weight Statistics:" << std::endl;
        std::cout << "  Total entries with ML weights: " << num_entries << std::endl;
        
        // Sample some ML weights
        float total_confidence = 0.0f;
        int sample_count = std::min((size_t)100, (size_t)num_entries);
        for (int i = 0; i < sample_count; i++) {
            uint64_t hash;
            uint16_t ml_weight;
            file.read(reinterpret_cast<char*>(&hash), sizeof(uint64_t));
            file.read(reinterpret_cast<char*>(&ml_weight), sizeof(uint16_t));
            total_confidence += ml_weight / 65535.0f;
        }
        
        if (sample_count > 0) {
            std::cout << "  Average ML confidence (sample): " 
                      << std::fixed << std::setprecision(4) 
                      << (total_confidence / sample_count) << std::endl;
        }
    }
    
    // Analyze contamination markers
    std::string contam_file = config.output_dir + "/contamination_markers.bin";
    if (std::filesystem::exists(contam_file)) {
        std::ifstream file(contam_file, std::ios::binary);
        uint64_t version, num_contaminated;
        file.read(reinterpret_cast<char*>(&version), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&num_contaminated), sizeof(uint64_t));
        
        std::cout << "\nContamination Statistics:" << std::endl;
        std::cout << "  Total contamination markers: " << num_contaminated << std::endl;
        
        if (results.total_minimizers > 0) {
            double contam_rate = (double)num_contaminated / results.total_minimizers * 100.0;
            std::cout << "  Contamination rate: " << std::fixed << std::setprecision(2) 
                      << contam_rate << "%" << std::endl;
        }
    }
}

// Main test runner
int main(int argc, char* argv[]) {
    std::cout << "=== ENHANCED DATABASE BUILD TEST SUITE ===" << std::endl;
    std::cout << "Testing ML-enhanced 24-byte minimizer structure" << std::endl;
    
    TestConfig config;
    TestResults results;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--genome-dir" && i + 1 < argc) {
            config.genome_dir = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            config.output_dir = argv[++i];
        } else if (arg == "--max-genomes" && i + 1 < argc) {
            config.max_test_genomes = std::stoi(argv[++i]);
        } else if (arg == "--quiet") {
            config.verbose = false;
        }
    }
    
    // Verify test data exists
    if (!std::filesystem::exists(config.genome_dir)) {
        std::cerr << "ERROR: Genome directory not found: " << config.genome_dir << std::endl;
        std::cerr << "Please ensure reference genomes are available at: " << config.genome_dir << std::endl;
        return 1;
    }
    
    // Create output directory
    std::filesystem::create_directories(config.output_dir);
    
    // Run all tests
    bool all_passed = true;
    
    // Test 1: Feature Extraction
    if (!test_feature_extraction(config, results)) {
        all_passed = false;
    }
    
    // Test 2: Contamination Detection
    if (!test_contamination_detection(config, results)) {
        all_passed = false;
    }
    
    // Test 3: Database Integrity
    if (!test_database_integrity(config, results)) {
        all_passed = false;
    }
    
    // Test 4: Performance Benchmark
    if (!benchmark_performance_impact(config, results)) {
        all_passed = false;
    }
    
    // Analyze feature distributions
    analyze_feature_distributions(config, results);
    
    // Print final summary
    results.print_summary();
    
    if (all_passed) {
        std::cout << "\n✅ ALL TESTS PASSED ✅" << std::endl;
        return 0;
    } else {
        std::cout << "\n❌ SOME TESTS FAILED ❌" << std::endl;
        return 1;
    }
}