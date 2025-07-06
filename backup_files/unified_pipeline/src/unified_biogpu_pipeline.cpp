// runtime/unified_pipeline/unified_biogpu_pipeline.cpp

#include <iostream>
#include <memory>
#include <vector>
#include <chrono>
#include "sample_csv_parser.h"
#include "streaming_fastq_reader.h"
#include "sequence_batch.h"
#include "gpu_sequence_buffer.h"

namespace BioGPU {

// Shared k-mer/minimizer generator
class SharedKmerGenerator {
private:
    // GPU memory for k-mers
    uint64_t* d_kmers;
    uint32_t* d_kmer_counts;
    void* bloom_filter;
    
    // Configuration
    struct Config {
        int kmer_size_nucleotide = 31;  // For FQ resistance & taxonomy
        int kmer_size_protein = 15;      // For AMR genes
        int minimizer_k = 15;
        int minimizer_w = 10;
        bool generate_reverse_complement = true;
    } config;
    
public:
    SharedKmerGenerator() {
        // Initialize GPU memory
        size_t max_kmers = 100000000;  // 100M k-mers
        cudaMalloc(&d_kmers, max_kmers * sizeof(uint64_t));
        cudaMalloc(&d_kmer_counts, sizeof(uint32_t));
        
        // Create shared bloom filter
        bloom_filter = create_bloom_filter(config.kmer_size_nucleotide);
    }
    
    ~SharedKmerGenerator() {
        cudaFree(d_kmers);
        cudaFree(d_kmer_counts);
        if (bloom_filter) destroy_bloom_filter(bloom_filter);
    }
    
    // Generate k-mers once, use for all pipelines
    void generateKmers(GPUSequenceBuffer* buffer, cudaStream_t stream) {
        // This generates k-mers that can be used by all three pipelines
        launch_kmer_generation_kernel(
            buffer->getSequences(),
            buffer->getOffsets(), 
            buffer->getLengths(),
            buffer->getCurrentSequences(),
            d_kmers,
            d_kmer_counts,
            config.kmer_size_nucleotide,
            config.generate_reverse_complement,
            stream
        );
    }
    
    uint64_t* getKmers() { return d_kmers; }
    uint32_t* getKmerCounts() { return d_kmer_counts; }
    void* getBloomFilter() { return bloom_filter; }
};

// Base class for analysis modules
class AnalysisModule {
public:
    virtual ~AnalysisModule() = default;
    virtual bool initialize(const std::string& db_path) = 0;
    virtual void processSequences(GPUSequenceBuffer* buffer, SharedKmerGenerator* kmer_gen, cudaStream_t stream) = 0;
    virtual void finalizeSample(const std::string& sample_name) = 0;
    virtual void generateReport(const std::string& output_path) = 0;
};

// FQ Resistance Module (adapted from clean_resistance_pipeline)
class FQResistanceModule : public AnalysisModule {
private:
    // Existing clean_resistance_pipeline components
    FQMutationDetectorCUDA detector;
    void* translated_search_engine;
    void* clinical_report_generator;
    AlleleFrequencyAnalyzer* allele_analyzer;
    
    // Configuration
    bool use_bloom_filter = true;
    bool use_smith_waterman = true;
    uint32_t min_allele_depth = 5;
    
public:
    bool initialize(const std::string& db_path) override {
        // Load FQ resistance database
        detector.loadIndex((db_path + "/nucleotide").c_str());
        
        // Initialize translated search
        translated_search_engine = create_translated_search_engine_with_sw(10000, use_smith_waterman);
        load_protein_database(translated_search_engine, (db_path + "/protein").c_str());
        
        allele_analyzer = new AlleleFrequencyAnalyzer(gene_id_to_name, species_id_to_name);
        return true;
    }
    
    void processSequences(GPUSequenceBuffer* buffer, SharedKmerGenerator* kmer_gen, cudaStream_t stream) override {
        // Use pre-generated k-mers from shared generator
        uint64_t* d_kmers = kmer_gen->getKmers();
        void* bloom_filter = kmer_gen->getBloomFilter();
        
        // Run FQ resistance detection pipeline stages
        // This reuses the logic from clean_resistance_pipeline
        runFQResistanceDetection(buffer, d_kmers, bloom_filter, stream);
    }
    
    void finalizeSample(const std::string& sample_name) override {
        // Generate allele frequencies
        auto frequencies = allele_analyzer->generateAlleleFrequencies(min_allele_depth);
        writeAlleleFrequenciesToCSV(frequencies, sample_name + "_fq_alleles.csv", sample_name);
    }
    
    void generateReport(const std::string& output_path) override {
        generate_clinical_report(clinical_report_generator);
    }
};

// AMR Gene Module
class AMRGeneModule : public AnalysisModule {
private:
    // AMR gene database
    std::vector<std::string> amr_genes;
    std::map<std::string, std::vector<uint64_t>> gene_kmers;
    
    // Coverage tracking
    std::map<std::string, std::vector<uint32_t>> gene_coverage;
    std::map<std::string, float> gene_rpkm;
    
public:
    bool initialize(const std::string& db_path) override {
        // Load AMR gene database (CARD, ResFinder, etc.)
        loadAMRGeneDatabase(db_path + "/amr_genes.fasta");
        buildGeneKmerIndex();
        return true;
    }
    
    void processSequences(GPUSequenceBuffer* buffer, SharedKmerGenerator* kmer_gen, cudaStream_t stream) override {
        // Use shared k-mers to calculate gene coverage
        calculateGeneCoverage(buffer, kmer_gen->getKmers(), stream);
    }
    
    void finalizeSample(const std::string& sample_name) override {
        // Calculate RPKM/TPM values
        calculateGeneAbundance();
        
        // Write gene quantification
        writeGeneQuantification(sample_name + "_amr_genes.tsv");
    }
    
    void generateReport(const std::string& output_path) override {
        // Generate AMR gene report with coverage plots
        generateAMRGeneReport(output_path);
    }
};

// Taxonomy Module  
class TaxonomyModule : public AnalysisModule {
private:
    // Taxonomy database (Kraken2-style)
    void* taxonomy_db;
    std::map<uint32_t, std::string> taxid_to_name;
    std::map<uint32_t, uint32_t> taxid_counts;
    
public:
    bool initialize(const std::string& db_path) override {
        // Load taxonomy database
        loadTaxonomyDB(db_path + "/taxonomy");
        return true;
    }
    
    void processSequences(GPUSequenceBuffer* buffer, SharedKmerGenerator* kmer_gen, cudaStream_t stream) override {
        // Classify reads using k-mers
        classifyReads(buffer, kmer_gen->getKmers(), stream);
    }
    
    void finalizeSample(const std::string& sample_name) override {
        // Generate taxonomic profile
        writeTaxonomicProfile(sample_name + "_taxonomy.tsv");
    }
    
    void generateReport(const std::string& output_path) override {
        // Generate Kraken-style report
        generateKrakenReport(output_path);
    }
};

// Main Unified Pipeline
class UnifiedBioGPUPipeline {
private:
    // Modules
    std::unique_ptr<FQResistanceModule> fq_module;
    std::unique_ptr<AMRGeneModule> amr_module;
    std::unique_ptr<TaxonomyModule> tax_module;
    
    // Shared components
    std::unique_ptr<SharedKmerGenerator> kmer_generator;
    std::unique_ptr<GPUSequenceBuffer> gpu_buffer;
    
    // Configuration
    struct Config {
        bool enable_fq_resistance = true;
        bool enable_amr_genes = true;
        bool enable_taxonomy = true;
        int batch_size = 100000;
        std::string output_dir = "results";
        bool save_intermediate = false;
    } config;
    
    // CUDA streams for overlap
    cudaStream_t stream_io;
    cudaStream_t stream_compute;
    
public:
    UnifiedBioGPUPipeline() {
        // Create CUDA streams
        cudaStreamCreate(&stream_io);
        cudaStreamCreate(&stream_compute);
        
        // Initialize shared components
        kmer_generator = std::make_unique<SharedKmerGenerator>();
        
        size_t max_sequences = config.batch_size;
        size_t max_bases = config.batch_size * 300;  // Assuming ~300bp reads
        gpu_buffer = std::make_unique<GPUSequenceBuffer>(max_sequences, max_bases);
        gpu_buffer->allocate();
    }
    
    ~UnifiedBioGPUPipeline() {
        cudaStreamDestroy(stream_io);
        cudaStreamDestroy(stream_compute);
    }
    
    void initialize(const std::string& db_base_path) {
        std::cout << "Initializing Unified BioGPU Pipeline..." << std::endl;
        
        if (config.enable_fq_resistance) {
            fq_module = std::make_unique<FQResistanceModule>();
            fq_module->initialize(db_base_path);
            std::cout << "  ✓ FQ Resistance module initialized" << std::endl;
        }
        
        if (config.enable_amr_genes) {
            amr_module = std::make_unique<AMRGeneModule>();
            amr_module->initialize(db_base_path);
            std::cout << "  ✓ AMR Gene module initialized" << std::endl;
        }
        
        if (config.enable_taxonomy) {
            tax_module = std::make_unique<TaxonomyModule>();
            tax_module->initialize(db_base_path);
            std::cout << "  ✓ Taxonomy module initialized" << std::endl;
        }
    }
    
    void processSample(const BioGPU::SampleInfo& sample) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::cout << "\nProcessing sample: " << sample.sample_name << std::endl;
        
        // Create output directory for sample
        std::string sample_output_dir = config.output_dir + "/" + sample.sample_name;
        std::string mkdir_cmd = "mkdir -p " + sample_output_dir;
        system(mkdir_cmd.c_str());
        
        // Open FASTQ files
        StreamingFastqReader reader(config.batch_size);
        bool is_paired = !sample.read2_path.empty();
        
        if (is_paired) {
            reader.openPaired(sample.read1_path, sample.read2_path);
        } else {
            reader.open(sample.read1_path);
        }
        
        // Process batches
        int batch_num = 0;
        while (reader.hasNext()) {
            auto batch = reader.getNextBatch();
            if (!batch || batch->empty()) break;
            
            // Transfer sequences to GPU
            gpu_buffer->transferBatchAsync(*batch, stream_io);
            
            // Wait for transfer to complete
            cudaStreamSynchronize(stream_io);
            
            // Generate k-mers once for all modules
            kmer_generator->generateKmers(gpu_buffer.get(), stream_compute);
            
            // Process with each enabled module
            if (config.enable_fq_resistance) {
                fq_module->processSequences(gpu_buffer.get(), kmer_generator.get(), stream_compute);
            }
            
            if (config.enable_amr_genes) {
                amr_module->processSequences(gpu_buffer.get(), kmer_generator.get(), stream_compute);
            }
            
            if (config.enable_taxonomy) {
                tax_module->processSequences(gpu_buffer.get(), kmer_generator.get(), stream_compute);
            }
            
            // Synchronize before next batch
            cudaStreamSynchronize(stream_compute);
            
            batch_num++;
            if (batch_num % 10 == 0) {
                std::cout << "  Processed " << (batch_num * config.batch_size) << " reads..." << std::endl;
            }
        }
        
        // Finalize sample processing
        if (config.enable_fq_resistance) {
            fq_module->finalizeSample(sample_output_dir + "/" + sample.sample_name);
        }
        
        if (config.enable_amr_genes) {
            amr_module->finalizeSample(sample_output_dir + "/" + sample.sample_name);
        }
        
        if (config.enable_taxonomy) {
            tax_module->finalizeSample(sample_output_dir + "/" + sample.sample_name);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "Sample " << sample.sample_name << " completed in " << duration.count() << " seconds" << std::endl;
    }
    
    void processBatch(const std::string& csv_path) {
        // Load sample list
        SampleCSVParser parser;
        if (!parser.parseFile(csv_path)) {
            std::cerr << "Failed to parse CSV file: " << csv_path << std::endl;
            return;
        }
        
        std::cout << "Processing " << parser.getSampleCount() << " samples..." << std::endl;
        
        // Process each sample
        for (size_t i = 0; i < parser.getSampleCount(); i++) {
            const SampleInfo* sample = parser.getSample(i);
            if (sample) {
                processSample(*sample);
            }
        }
        
        // Generate combined reports
        generateCombinedReports();
    }
    
    void generateCombinedReports() {
        std::cout << "\nGenerating combined reports..." << std::endl;
        
        if (config.enable_fq_resistance) {
            fq_module->generateReport(config.output_dir + "/combined_fq_resistance_report.html");
        }
        
        if (config.enable_amr_genes) {
            amr_module->generateReport(config.output_dir + "/combined_amr_genes_report.html");
        }
        
        if (config.enable_taxonomy) {
            tax_module->generateReport(config.output_dir + "/combined_taxonomy_report.html");
        }
        
        // Generate master summary
        generateMasterSummary();
    }
    
    void generateMasterSummary() {
        std::string summary_path = config.output_dir + "/analysis_summary.html";
        std::ofstream summary(summary_path);
        
        summary << R"(
<!DOCTYPE html>
<html>
<head>
    <title>BioGPU Analysis Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .module { border: 1px solid #ddd; padding: 15px; margin: 10px 0; }
        .complete { background-color: #d4edda; }
        .warning { background-color: #fff3cd; }
        h2 { color: #0056b3; }
        .timestamp { color: #666; font-size: 0.9em; }
    </style>
</head>
<body>
    <h1>BioGPU Unified Analysis Summary</h1>
    <p class="timestamp">Generated: )" << getCurrentTimestamp() << R"(</p>
    
    <div class="module complete">
        <h2>Fluoroquinolone Resistance Detection</h2>
        <p>Status: Complete</p>
        <p><a href="combined_fq_resistance_report.html">View Full Report</a></p>
        <ul>
            <li>Samples analyzed: )" << sample_count << R"(</li>
            <li>Resistance mutations detected: )" << total_mutations << R"(</li>
            <li>Samples with resistance: )" << resistant_samples << R"(</li>
        </ul>
    </div>
    
    <div class="module complete">
        <h2>AMR Gene Quantification</h2>
        <p>Status: Complete</p>
        <p><a href="combined_amr_genes_report.html">View Full Report</a></p>
        <ul>
            <li>AMR genes detected: )" << amr_genes_found << R"(</li>
            <li>Average coverage: )" << avg_coverage << R"(x</li>
        </ul>
    </div>
    
    <div class="module complete">
        <h2>Taxonomic Classification</h2>
        <p>Status: Complete</p>
        <p><a href="combined_taxonomy_report.html">View Full Report</a></p>
        <ul>
            <li>Species identified: )" << species_count << R"(</li>
            <li>Dominant organism: )" << dominant_species << R"(</li>
        </ul>
    </div>
</body>
</html>
        )";
        
        summary.close();
        std::cout << "Master summary written to: " << summary_path << std::endl;
    }
    
    // Configuration methods
    void setConfig(const std::string& key, const std::string& value) {
        if (key == "enable_fq_resistance") config.enable_fq_resistance = (value == "true");
        else if (key == "enable_amr_genes") config.enable_amr_genes = (value == "true");
        else if (key == "enable_taxonomy") config.enable_taxonomy = (value == "true");
        else if (key == "batch_size") config.batch_size = std::stoi(value);
        else if (key == "output_dir") config.output_dir = value;
        else if (key == "save_intermediate") config.save_intermediate = (value == "true");
    }
};

} // namespace BioGPU

// Main program
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <database_path> --csv <samples.csv> [options]\n";
        std::cerr << "\nOptions:\n";
        std::cerr << "  --enable-fq          Enable FQ resistance detection (default: true)\n";
        std::cerr << "  --enable-amr         Enable AMR gene quantification (default: true)\n";
        std::cerr << "  --enable-taxonomy    Enable taxonomic classification (default: true)\n";
        std::cerr << "  --output-dir <dir>   Output directory (default: results)\n";
        std::cerr << "  --batch-size <n>     Batch size (default: 100000)\n";
        std::cerr << "\nRun specific modules only:\n";
        std::cerr << "  --only-fq            Run only FQ resistance detection\n";
        std::cerr << "  --only-amr           Run only AMR gene quantification\n";
        std::cerr << "  --only-taxonomy      Run only taxonomic classification\n";
        return 1;
    }
    
    std::string db_path = argv[1];
    std::string csv_path;
    
    // Create pipeline
    BioGPU::UnifiedBioGPUPipeline pipeline;
    
    // Parse arguments
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--csv" && i + 1 < argc) {
            csv_path = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            pipeline.setConfig("output_dir", argv[++i]);
        } else if (arg == "--batch-size" && i + 1 < argc) {
            pipeline.setConfig("batch_size", argv[++i]);
        } else if (arg == "--only-fq") {
            pipeline.setConfig("enable_fq_resistance", "true");
            pipeline.setConfig("enable_amr_genes", "false");
            pipeline.setConfig("enable_taxonomy", "false");
        } else if (arg == "--only-amr") {
            pipeline.setConfig("enable_fq_resistance", "false");
            pipeline.setConfig("enable_amr_genes", "true");
            pipeline.setConfig("enable_taxonomy", "false");
        } else if (arg == "--only-taxonomy") {
            pipeline.setConfig("enable_fq_resistance", "false");
            pipeline.setConfig("enable_amr_genes", "false");
            pipeline.setConfig("enable_taxonomy", "true");
        }
    }
    
    if (csv_path.empty()) {
        std::cerr << "Error: No CSV file specified\n";
        return 1;
    }
    
    std::cout << "=== BioGPU Unified Pipeline ===" << std::endl;
    std::cout << "Database: " << db_path << std::endl;
    std::cout << "Samples: " << csv_path << std::endl;
    
    try {
        // Initialize pipeline
        pipeline.initialize(db_path);
        
        // Process samples
        pipeline.processBatch(csv_path);
        
        std::cout << "\nPipeline completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
