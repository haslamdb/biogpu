// runtime/kernels/resistance/test_resistance_pipeline_v2.cpp
// Comprehensive test suite for validating the refactored diagnostic pipeline
// Tests known fluoroquinolone resistance mutations with clinical significance

#include "resistance_pipeline_v2.h"
#include "../../common/io/streaming_fastq_reader.h"
#include "../../common/gpu/gpu_sequence_buffer.h"
#include "../../common/config/unified_config.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <iomanip>
#include <cstring>

using namespace BioGPU;

// Test sequence generator for known mutations
class TestSequenceGenerator {
private:
    std::mt19937 rng;
    
    // Reference sequences around mutation sites
    struct MutationContext {
        std::string gene;
        std::string species;
        std::string upstream;
        std::string wildtype_codon;
        std::string mutant_codon;
        std::string downstream;
        int position;
        char wildtype_aa;
        char mutant_aa;
        std::string clinical_significance;
    };
    
    std::vector<MutationContext> known_mutations = {
        // E. coli gyrA S83L - High-level ciprofloxacin resistance
        {"gyrA", "Escherichia coli", 
         "ATGAGCGACCTTGCGAGAGAAATTACACCGGTCAACATTGAGGAAGAGCTGAAGAGCTCC",
         "TCG", "TTG",  // S -> L
         "GGTTACGCGCAGACCGTAGTTAAACTTGTACCGCAGGTGCTATGACCATCGATGTTTAC",
         83, 'S', 'L',
         "High-level fluoroquinolone resistance - avoid ciprofloxacin/levofloxacin"},
        
        // K. pneumoniae parC S80I
        {"parC", "Klebsiella pneumoniae",
         "ATGGAAACCTACCCGCATAAAGGTAAGAGCATTACGCTGCGCGTGATGTACGCGGCTGAA",
         "AGC", "ATC",  // S -> I
         "ATTGAGCGGCTTGTTTACGACAACCTGATCAAGCGCAGCCACCTGGCGGTAGAACGTAAA",
         80, 'S', 'I',
         "Moderate fluoroquinolone resistance - consider higher dose or alternative"},
        
        // S. aureus grlA S80F
        {"grlA", "Staphylococcus aureus",
         "ATGAGTGAAATCGATAATCGTGAAATCTATCAAGAAGAAACAATAGAAGATGAAATTGAT",
         "TCA", "TTC",  // S -> F
         "AAAGGGTTAACTCATCGTTTAGTTAAAGTTTATCCAAAAGTTCTTTTTACTATTGATACA",
         80, 'S', 'F',
         "Fluoroquinolone resistance in MRSA - use alternative antibiotics"},
        
        // E. coli gyrA D87N
        {"gyrA", "Escherichia coli",
         "AGAGAAATTACACCGGTCAACATTGAGGAAGAGCTGAAGAGCTCCTATCTGGATTATGCG",
         "GAC", "AAC",  // D -> N
         "ATGTATCGTGTCCTTTATGTGCTGGGAGCCGTACACCGTCATAACGGTTCAGTTCGCCTA",
         87, 'D', 'N',
         "High-level fluoroquinolone resistance - contraindicated"},
        
        // Control: Susceptible E. coli (no mutations)
        {"gyrA", "Escherichia coli",
         "ATGAGCGACCTTGCGAGAGAAATTACACCGGTCAACATTGAGGAAGAGCTGAAGAGCTCC",
         "TCG", "TCG",  // S -> S (no change)
         "GGTTACGCGCAGACCGTAGTTAAACTTGTACCGCAGGTGCTATGACCATCGATGTTTAC",
         83, 'S', 'S',
         "Susceptible - fluoroquinolones effective"}
    };
    
public:
    TestSequenceGenerator() : rng(42) {}  // Fixed seed for reproducibility
    
    // Generate a read containing a specific mutation
    std::string generateMutantRead(int mutation_idx, bool add_errors = false) {
        if (mutation_idx >= known_mutations.size()) return "";
        
        const auto& mut = known_mutations[mutation_idx];
        std::string sequence = mut.upstream + mut.mutant_codon + mut.downstream;
        
        // Add sequencing errors if requested
        if (add_errors) {
            std::uniform_real_distribution<> error_dist(0.0, 1.0);
            std::uniform_int_distribution<> base_dist(0, 3);
            const char bases[] = "ACGT";
            
            for (size_t i = 0; i < sequence.length(); i++) {
                if (error_dist(rng) < 0.01) {  // 1% error rate
                    sequence[i] = bases[base_dist(rng)];
                }
            }
        }
        
        return sequence;
    }
    
    // Generate wildtype (susceptible) read
    std::string generateWildtypeRead(int mutation_idx) {
        if (mutation_idx >= known_mutations.size()) return "";
        
        const auto& mut = known_mutations[mutation_idx];
        return mut.upstream + mut.wildtype_codon + mut.downstream;
    }
    
    // Generate quality string
    std::string generateQuality(size_t length) {
        std::string qual(length, 'I');  // High quality (Q40)
        return qual;
    }
    
    // Get mutation info
    const MutationContext& getMutation(int idx) const {
        return known_mutations[idx];
    }
    
    size_t getMutationCount() const { return known_mutations.size(); }
};

// Test result validator
class DiagnosticValidator {
private:
    struct ExpectedResult {
        std::string gene;
        int position;
        char wildtype_aa;
        char mutant_aa;
        float min_frequency;
        std::string clinical_interpretation;
    };
    
public:
    bool validateResults(const ResistanceResults& results,
                        const std::vector<ExpectedResult>& expected) {
        bool all_found = true;
        
        std::cout << "\n=== Diagnostic Validation Results ===\n";
        
        for (const auto& exp : expected) {
            bool found = false;
            
            for (const auto& hit : results.all_hits) {
                if (hit.gene_name == exp.gene && 
                    hit.position == exp.position &&
                    hit.wildtype_aa == exp.wildtype_aa &&
                    hit.mutant_aa == exp.mutant_aa) {
                    found = true;
                    
                    std::cout << "✓ DETECTED: " << hit.mutation_name 
                              << " (confidence: " << std::fixed << std::setprecision(2) 
                              << hit.confidence * 100 << "%)\n";
                    std::cout << "  Clinical significance: " << exp.clinical_interpretation << "\n";
                    
                    // Check if marked as QRDR
                    if (hit.is_qrdr) {
                        std::cout << "  Located in QRDR region - critical for resistance\n";
                    }
                    break;
                }
            }
            
            if (!found) {
                std::cout << "✗ MISSED: " << exp.gene << "_" << exp.wildtype_aa 
                          << exp.position << exp.mutant_aa << "\n";
                all_found = false;
            }
        }
        
        return all_found;
    }
    
    void printAlleleFrequencies(const ResistanceResults& results) {
        std::cout << "\n=== Allele Frequency Analysis ===\n";
        
        for (const auto& [sample, profile] : results.sample_profiles) {
            std::cout << "Sample: " << sample << "\n";
            
            for (const auto& [gene, frequencies] : profile.allele_frequencies) {
                std::cout << "  Gene: " << gene << "\n";
                
                for (const auto& freq : frequencies) {
                    std::cout << "    Position " << freq.position << ": ";
                    
                    // Sort alleles by frequency
                    std::vector<std::pair<char, float>> sorted_alleles(
                        freq.amino_acid_frequencies.begin(), 
                        freq.amino_acid_frequencies.end()
                    );
                    std::sort(sorted_alleles.begin(), sorted_alleles.end(),
                             [](const auto& a, const auto& b) { return a.second > b.second; });
                    
                    for (const auto& [aa, frequency] : sorted_alleles) {
                        if (frequency > 0.01) {  // Show alleles >1%
                            std::cout << aa << ":" << std::fixed << std::setprecision(1) 
                                     << frequency * 100 << "% ";
                            if (aa == freq.wildtype_aa) std::cout << "(WT) ";
                        }
                    }
                    std::cout << " [coverage: " << freq.total_coverage << "]\n";
                }
            }
        }
    }
};

// Performance benchmarking
class PerformanceBenchmark {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    size_t total_reads = 0;
    size_t total_bases = 0;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    void recordBatch(size_t reads, size_t bases) {
        total_reads += reads;
        total_bases += bases;
    }
    
    void printReport() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        double seconds = duration.count();
        
        std::cout << "\n=== Performance Metrics ===\n";
        std::cout << "Total processing time: " << seconds << " seconds\n";
        std::cout << "Total reads processed: " << total_reads << "\n";
        std::cout << "Total bases processed: " << total_bases << "\n";
        std::cout << "Throughput: " << std::fixed << std::setprecision(0) 
                  << total_reads / seconds << " reads/second\n";
        std::cout << "Throughput: " << std::fixed << std::setprecision(2) 
                  << total_bases / seconds / 1e6 << " Mbases/second\n";
        
        // GPU memory usage
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        size_t used_mem = total_mem - free_mem;
        std::cout << "GPU memory usage: " << used_mem / 1024 / 1024 << " MB / " 
                  << total_mem / 1024 / 1024 << " MB\n";
    }
};

// Main test function
int main(int argc, char* argv[]) {
    std::cout << "=== Fluoroquinolone Resistance Diagnostic Pipeline Test ===\n";
    std::cout << "Testing clinical detection of antibiotic resistance mutations\n\n";
    
    // Load configuration
    UnifiedConfig config = ConfigLoader::getDefault();
    
    // Configure for diagnostic testing
    config.resistance_config.batch_size = 10000;
    config.resistance_config.min_identity = 0.95f;
    config.resistance_config.min_coverage = 0.80f;
    config.resistance_config.calculate_allele_frequencies = true;
    config.resistance_config.generate_clinical_report = true;
    config.resistance_config.verbose = true;
    
    // Create resistance pipeline
    ResistancePipeline pipeline(config.resistance_config);
    
    // Initialize pipeline
    std::cout << "Initializing diagnostic pipeline...\n";
    if (!pipeline.initialize()) {
        std::cerr << "Failed to initialize resistance pipeline\n";
        return 1;
    }
    
    // Test 1: Known mutation detection
    std::cout << "\n=== Test 1: Known Mutation Detection ===\n";
    testKnownMutations(pipeline, config);
    
    // Test 2: Mixed sample with multiple species
    std::cout << "\n=== Test 2: Mixed Sample Analysis ===\n";
    testMixedSample(pipeline, config);
    
    // Test 3: Low-frequency mutation detection
    std::cout << "\n=== Test 3: Low-Frequency Mutation Detection ===\n";
    testLowFrequencyMutations(pipeline, config);
    
    // Test 4: Performance benchmarking
    if (argc > 1) {
        std::cout << "\n=== Test 4: Performance Benchmark ===\n";
        testPerformance(pipeline, config, argv[1]);
    }
    
    // Test 5: Clinical report generation
    std::cout << "\n=== Test 5: Clinical Report Generation ===\n";
    testClinicalReports(pipeline);
    
    return 0;
}

// Test known mutations
void testKnownMutations(ResistancePipeline& pipeline, const UnifiedConfig& config) {
    TestSequenceGenerator generator;
    DiagnosticValidator validator;
    
    // Create GPU buffer
    GPUSequenceBuffer gpu_buffer(10000, 3000000);
    gpu_buffer.allocate();
    
    // Expected results
    std::vector<DiagnosticValidator::ExpectedResult> expected = {
        {"gyrA", 83, 'S', 'L', 0.95f, "High-level fluoroquinolone resistance"},
        {"parC", 80, 'S', 'I', 0.95f, "Moderate fluoroquinolone resistance"},
        {"grlA", 80, 'S', 'F', 0.95f, "Fluoroquinolone resistance in S. aureus"},
        {"gyrA", 87, 'D', 'N', 0.95f, "High-level fluoroquinolone resistance"}
    };
    
    // Generate test batch
    SequenceBatch batch;
    batch.sample_name = "diagnostic_test_known_mutations";
    
    // Add mutant reads (100 copies of each mutation for strong signal)
    for (size_t i = 0; i < generator.getMutationCount() - 1; i++) {  // Skip control
        for (int copy = 0; copy < 100; copy++) {
            std::string seq = generator.generateMutantRead(i, copy % 10 == 0);  // 10% with errors
            std::string qual = generator.generateQuality(seq.length());
            std::string header = "@read_" + std::to_string(i) + "_" + std::to_string(copy);
            batch.addRead(header, seq, qual);
        }
    }
    
    // Add some wildtype reads
    for (int i = 0; i < 50; i++) {
        std::string seq = generator.generateWildtypeRead(0);  // E. coli wildtype
        std::string qual = generator.generateQuality(seq.length());
        std::string header = "@wildtype_" + std::to_string(i);
        batch.addRead(header, seq, qual);
    }
    
    std::cout << "Generated " << batch.size() << " test reads\n";
    
    // Process batch
    pipeline.setCurrentSample(batch.sample_name);
    gpu_buffer.transferBatch(batch);
    pipeline.processBatch(&gpu_buffer);
    
    // Get results
    auto results = pipeline.getResults();
    auto* resistance_results = dynamic_cast<ResistanceResults*>(results.get());
    
    // Validate
    bool all_detected = validator.validateResults(*resistance_results, expected);
    
    if (all_detected) {
        std::cout << "\n✅ All known mutations correctly detected!\n";
    } else {
        std::cout << "\n❌ Some mutations were missed\n";
    }
}

// Test mixed sample
void testMixedSample(ResistancePipeline& pipeline, const UnifiedConfig& config) {
    TestSequenceGenerator generator;
    
    // Create GPU buffer
    GPUSequenceBuffer gpu_buffer(10000, 3000000);
    gpu_buffer.allocate();
    
    // Generate mixed species batch
    SequenceBatch batch;
    batch.sample_name = "mixed_species_sample";
    
    // E. coli reads (60%)
    for (int i = 0; i < 600; i++) {
        std::string seq = i < 300 ? generator.generateMutantRead(0) : 
                                   generator.generateWildtypeRead(0);
        std::string qual = generator.generateQuality(seq.length());
        batch.addRead("@ecoli_" + std::to_string(i), seq, qual);
    }
    
    // K. pneumoniae reads (30%)
    for (int i = 0; i < 300; i++) {
        std::string seq = i < 200 ? generator.generateMutantRead(1) : 
                                   generator.generateWildtypeRead(1);
        std::string qual = generator.generateQuality(seq.length());
        batch.addRead("@kpneu_" + std::to_string(i), seq, qual);
    }
    
    // S. aureus reads (10%)
    for (int i = 0; i < 100; i++) {
        std::string seq = generator.generateMutantRead(2);
        std::string qual = generator.generateQuality(seq.length());
        batch.addRead("@saureus_" + std::to_string(i), seq, qual);
    }
    
    std::cout << "Processing mixed sample with " << batch.size() << " reads\n";
    std::cout << "Expected: E. coli (60%), K. pneumoniae (30%), S. aureus (10%)\n";
    
    // Process
    pipeline.setCurrentSample(batch.sample_name);
    gpu_buffer.transferBatch(batch);
    pipeline.processBatch(&gpu_buffer);
    
    // Calculate allele frequencies
    pipeline.calculateAlleleFrequencies();
    
    // Get results
    auto results = pipeline.getResults();
    auto* resistance_results = dynamic_cast<ResistanceResults*>(results.get());
    
    // Print species distribution
    std::map<std::string, int> species_counts;
    for (const auto& hit : resistance_results->all_hits) {
        species_counts[hit.species]++;
    }
    
    std::cout << "\nDetected species distribution:\n";
    for (const auto& [species, count] : species_counts) {
        std::cout << "  " << species << ": " << count << " resistance markers\n";
    }
}

// Test low-frequency mutations
void testLowFrequencyMutations(ResistancePipeline& pipeline, const UnifiedConfig& config) {
    TestSequenceGenerator generator;
    DiagnosticValidator validator;
    
    // Create GPU buffer
    GPUSequenceBuffer gpu_buffer(10000, 3000000);
    gpu_buffer.allocate();
    
    // Generate batch with varying mutation frequencies
    SequenceBatch batch;
    batch.sample_name = "low_frequency_test";
    
    // 95% wildtype
    for (int i = 0; i < 950; i++) {
        std::string seq = generator.generateWildtypeRead(0);
        std::string qual = generator.generateQuality(seq.length());
        batch.addRead("@wt_" + std::to_string(i), seq, qual);
    }
    
    // 5% mutant
    for (int i = 0; i < 50; i++) {
        std::string seq = generator.generateMutantRead(0);  // gyrA S83L
        std::string qual = generator.generateQuality(seq.length());
        batch.addRead("@mut_" + std::to_string(i), seq, qual);
    }
    
    std::cout << "Testing detection of 5% mutation frequency (50/1000 reads)\n";
    
    // Process
    pipeline.setCurrentSample(batch.sample_name);
    gpu_buffer.transferBatch(batch);
    pipeline.processBatch(&gpu_buffer);
    pipeline.calculateAlleleFrequencies();
    
    // Get results
    auto results = pipeline.getResults();
    auto* resistance_results = dynamic_cast<ResistanceResults*>(results.get());
    
    // Check allele frequencies
    validator.printAlleleFrequencies(*resistance_results);
    
    // Verify low-frequency detection
    bool detected_low_freq = false;
    for (const auto& hit : resistance_results->all_hits) {
        if (hit.gene_name == "gyrA" && hit.position == 83) {
            std::cout << "\n✅ Low-frequency mutation detected: " << hit.mutation_name << "\n";
            std::cout << "   Read support: " << hit.read_support << " reads\n";
            detected_low_freq = true;
        }
    }
    
    if (!detected_low_freq) {
        std::cout << "\n❌ Failed to detect low-frequency mutation\n";
    }
}

// Test performance
void testPerformance(ResistancePipeline& pipeline, const UnifiedConfig& config, 
                    const std::string& fastq_file) {
    
    PerformanceBenchmark benchmark;
    
    // Create reader and buffer
    StreamingFastqReader reader(config.resistance_config.batch_size);
    GPUSequenceBuffer gpu_buffer(
        config.resistance_config.batch_size * 2,
        config.resistance_config.batch_size * 300
    );
    
    if (!reader.open(fastq_file)) {
        std::cerr << "Cannot open file: " << fastq_file << "\n";
        return;
    }
    
    if (!gpu_buffer.allocate()) {
        std::cerr << "Failed to allocate GPU buffer\n";
        return;
    }
    
    benchmark.start();
    
    // Process file
    size_t batch_count = 0;
    while (reader.hasNext()) {
        auto batch = reader.getNextBatch();
        if (!batch || batch->empty()) break;
        
        batch_count++;
        benchmark.recordBatch(batch->size(), batch->getTotalBases());
        
        // Process batch
        pipeline.setCurrentSample("performance_test");
        gpu_buffer.transferBatch(*batch);
        pipeline.processBatch(&gpu_buffer);
        
        if (batch_count % 10 == 0) {
            std::cout << "Processed " << batch_count << " batches...\r" << std::flush;
        }
    }
    
    std::cout << "\nProcessed " << batch_count << " total batches\n";
    benchmark.printReport();
}

// Test clinical report generation
void testClinicalReports(ResistancePipeline& pipeline) {
    // Get accumulated results
    auto results = pipeline.getResults();
    auto* resistance_results = dynamic_cast<ResistanceResults*>(results.get());
    
    // Generate reports in different formats
    std::string output_dir = "diagnostic_test_output";
    system(("mkdir -p " + output_dir).c_str());
    
    // HTML report
    resistance_results->writeClinicalReport(output_dir + "/clinical_report.html", "html");
    std::cout << "Generated HTML clinical report: " << output_dir << "/clinical_report.html\n";
    
    // JSON report
    resistance_results->writeJSONReport(output_dir + "/resistance_data.json");
    std::cout << "Generated JSON report: " << output_dir << "/resistance_data.json\n";
    
    // TSV report
    resistance_results->writeReport(output_dir + "/resistance_mutations.tsv");
    std::cout << "Generated TSV report: " << output_dir << "/resistance_mutations.tsv\n";
    
    // Allele frequency CSV
    resistance_results->writeAlleleFrequencyReport(output_dir + "/allele_frequencies.csv");
    std::cout << "Generated allele frequency report: " << output_dir << "/allele_frequencies.csv\n";
    
    // Print sample clinical interpretation
    std::cout << "\n=== Sample Clinical Interpretation ===\n";
    std::cout << "Based on detected mutations:\n\n";
    
    for (const auto& [sample, profile] : resistance_results->sample_profiles) {
        std::cout << "Sample: " << sample << "\n";
        
        if (profile.has_fluoroquinolone_resistance) {
            std::cout << "⚠️  FLUOROQUINOLONE RESISTANCE DETECTED\n";
            std::cout << "Resistance mechanisms:\n";
            for (const auto& mechanism : profile.resistance_mechanisms) {
                std::cout << "  • " << mechanism << "\n";
            }
            std::cout << "\nTreatment recommendation:\n";
            std::cout << "  - Avoid fluoroquinolones (ciprofloxacin, levofloxacin)\n";
            std::cout << "  - Consider alternative antibiotics based on susceptibility testing\n";
            std::cout << "  - For UTI: nitrofurantoin, fosfomycin, or beta-lactams\n";
            std::cout << "  - For systemic infections: carbapenems or aminoglycosides\n";
        } else {
            std::cout << "✅ No fluoroquinolone resistance detected\n";
            std::cout << "  - Fluoroquinolones likely effective at standard doses\n";
        }
        std::cout << "\n";
    }
}