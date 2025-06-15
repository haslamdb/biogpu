// kraken_pipeline_main.cpp
// Complete front-to-back GPU-accelerated Kraken2-style pipeline
// Build database from genomes + classify reads

// Forward declarations for the classes we need
class ClassificationParams;
class GPUKrakenDatabaseBuilder;
class PairedEndGPUKrakenClassifier;

// Include the actual headers with class definitions
#define GPU_KRAKEN_CLASSIFIER_HEADER_ONLY
#define GPU_KRAKEN_DATABASE_BUILDER_HEADER_ONLY
#include "gpu_kraken_classifier.cu"
#include "gpu_kraken_database_builder.cu"
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono>

void print_usage(const char* program_name) {
    std::cout << "GPU-Accelerated Kraken2-Style Pipeline" << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir <dir> --output <db_dir> [options]" << std::endl;
    std::cout << "  " << program_name << " classify --database <db_dir> --reads <fastq> [options]" << std::endl;
    std::cout << "  " << program_name << " pipeline --genome-dir <dir> --reads <fastq> --output <dir> [options]" << std::endl;
    std::cout << "\nCommands:" << std::endl;
    std::cout << "  build      Build database from genome files" << std::endl;
    std::cout << "  classify   Classify reads using existing database" << std::endl;
    std::cout << "  pipeline   Complete pipeline: build database + classify reads" << std::endl;
    std::cout << "\nRequired Arguments:" << std::endl;
    std::cout << "  --genome-dir <dir>    Directory containing genome FASTA files" << std::endl;
    std::cout << "  --database <db_dir>   Database directory (for classify mode)" << std::endl;
    std::cout << "  --reads <fastq>       FASTQ file to classify (R1 for paired-end)" << std::endl;
    std::cout << "  --output <dir>        Output directory" << std::endl;
    std::cout << "\nOptional Arguments:" << std::endl;
    std::cout << "  --reads2 <fastq>      R2 FASTQ file for paired-end reads" << std::endl;
    std::cout << "  --k <int>             k-mer length (default: 35)" << std::endl;
    std::cout << "  --minimizer-len <int> Minimizer length (default: 31)" << std::endl;
    std::cout << "  --spaces <int>        Spaced seed spacing (default: 7)" << std::endl;
    std::cout << "  --confidence <float>  Confidence threshold 0-1 (default: 0.1)" << std::endl;
    std::cout << "  --batch-size <int>    Reads per batch (default: 10000)" << std::endl;
    std::cout << "  --gpu-batch <int>     GPU batch size for building (default: 1000)" << std::endl;
    std::cout << "  --taxonomy <dir>      NCBI taxonomy directory (optional)" << std::endl;
    std::cout << "  --threads <int>       CPU threads for I/O (default: 4)" << std::endl;
    std::cout << "  --quick               Quick mode: stop at first hit" << std::endl;
    std::cout << "  --report <file>       Generate summary report" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  # Build database from RefSeq genomes" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir ./refseq_genomes --output ./kraken_db" << std::endl;
    std::cout << "\n  # Classify single-end reads" << std::endl;
    std::cout << "  " << program_name << " classify --database ./kraken_db --reads sample.fastq.gz --output results.txt" << std::endl;
    std::cout << "\n  # Classify paired-end reads" << std::endl;
    std::cout << "  " << program_name << " classify --database ./kraken_db --reads sample_R1.fastq.gz --reads2 sample_R2.fastq.gz --output results.txt" << std::endl;
    std::cout << "\n  # Complete pipeline" << std::endl;
    std::cout << "  " << program_name << " pipeline --genome-dir ./genomes --reads sample.fastq.gz --output ./results" << std::endl;
}

struct PipelineConfig {
    std::string command;
    std::string genome_dir;
    std::string database_dir;
    std::string reads_file;
    std::string reads_file2;  // For paired-end reads
    std::string output_path;
    std::string taxonomy_dir;
    std::string report_file;
    
    ClassificationParams classifier_params;
    int batch_size = 10000;
    int gpu_batch_size = 1000;  // GPU batch size for database building
    int cpu_threads = 4;
    bool quick_mode = false;
    bool verbose = false;
};

bool parse_arguments(int argc, char* argv[], PipelineConfig& config) {
    if (argc < 2) {
        return false;
    }
    
    config.command = argv[1];
    
    if (config.command != "build" && config.command != "classify" && config.command != "pipeline") {
        std::cerr << "Error: Unknown command '" << config.command << "'" << std::endl;
        return false;
    }
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--genome-dir" && i + 1 < argc) {
            config.genome_dir = argv[++i];
        } else if (arg == "--database" && i + 1 < argc) {
            config.database_dir = argv[++i];
        } else if (arg == "--reads" && i + 1 < argc) {
            config.reads_file = argv[++i];
        } else if (arg == "--reads2" && i + 1 < argc) {
            config.reads_file2 = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            config.output_path = argv[++i];
        } else if (arg == "--taxonomy" && i + 1 < argc) {
            config.taxonomy_dir = argv[++i];
        } else if (arg == "--report" && i + 1 < argc) {
            config.report_file = argv[++i];
        } else if (arg == "--k" && i + 1 < argc) {
            config.classifier_params.k = std::stoi(argv[++i]);
        } else if (arg == "--minimizer-len" && i + 1 < argc) {
            config.classifier_params.ell = std::stoi(argv[++i]);
        } else if (arg == "--spaces" && i + 1 < argc) {
            config.classifier_params.spaces = std::stoi(argv[++i]);
        } else if (arg == "--confidence" && i + 1 < argc) {
            config.classifier_params.confidence_threshold = std::stof(argv[++i]);
        } else if (arg == "--batch-size" && i + 1 < argc) {
            config.batch_size = std::stoi(argv[++i]);
        } else if (arg == "--gpu-batch" && i + 1 < argc) {
            config.gpu_batch_size = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            config.cpu_threads = std::stoi(argv[++i]);
        } else if (arg == "--quick") {
            config.quick_mode = true;
        } else if (arg == "--verbose") {
            config.verbose = true;
        } else {
            std::cerr << "Warning: Unknown argument '" << arg << "'" << std::endl;
        }
    }
    
    return true;
}

bool validate_config(const PipelineConfig& config) {
    if (config.command == "build" || config.command == "pipeline") {
        if (config.genome_dir.empty()) {
            std::cerr << "Error: --genome-dir is required for build/pipeline mode" << std::endl;
            return false;
        }
        if (!std::filesystem::exists(config.genome_dir)) {
            std::cerr << "Error: Genome directory does not exist: " << config.genome_dir << std::endl;
            return false;
        }
    }
    
    if (config.command == "classify") {
        if (config.database_dir.empty()) {
            std::cerr << "Error: --database is required for classify mode" << std::endl;
            return false;
        }
        if (!std::filesystem::exists(config.database_dir)) {
            std::cerr << "Error: Database directory does not exist: " << config.database_dir << std::endl;
            return false;
        }
    }
    
    if (config.command == "classify" || config.command == "pipeline") {
        if (config.reads_file.empty()) {
            std::cerr << "Error: --reads is required for classify/pipeline mode" << std::endl;
            return false;
        }
        if (!std::filesystem::exists(config.reads_file)) {
            std::cerr << "Error: Reads file does not exist: " << config.reads_file << std::endl;
            return false;
        }
    }
    
    if (config.output_path.empty()) {
        std::cerr << "Error: --output is required" << std::endl;
        return false;
    }
    
    // Validate parameters
    if (config.classifier_params.k < 10 || config.classifier_params.k > 50) {
        std::cerr << "Error: k-mer length must be between 10 and 50" << std::endl;
        return false;
    }
    
    if (config.classifier_params.ell >= config.classifier_params.k) {
        std::cerr << "Error: Minimizer length must be less than k-mer length" << std::endl;
        return false;
    }
    
    if (config.classifier_params.confidence_threshold < 0.0f || 
        config.classifier_params.confidence_threshold > 1.0f) {
        std::cerr << "Error: Confidence threshold must be between 0.0 and 1.0" << std::endl;
        return false;
    }
    
    return true;
}

bool build_database_command(const PipelineConfig& config) {
    std::cout << "\n=== BUILDING KRAKEN DATABASE ===" << std::endl;
    std::cout << "Genome directory: " << config.genome_dir << std::endl;
    std::cout << "Output directory: " << config.output_path << std::endl;
    std::cout << "Parameters: k=" << config.classifier_params.k 
              << ", ell=" << config.classifier_params.ell 
              << ", spaces=" << config.classifier_params.spaces << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        GPUKrakenDatabaseBuilder builder(config.output_path, config.classifier_params);
        
        // Set batch size if specified
        if (config.gpu_batch_size > 0) {
            builder.set_batch_size(config.gpu_batch_size);
        }
        
        bool success = builder.build_database_from_genomes(
            config.genome_dir,
            config.taxonomy_dir
        );
        
        if (!success) {
            std::cerr << "Database build failed!" << std::endl;
            return false;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\n✓ Database build completed successfully in " 
                  << duration.count() << " seconds" << std::endl;
        std::cout << "Database saved to: " << config.output_path << std::endl;
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error building database: " << e.what() << std::endl;
        return false;
    }
}

bool classify_reads_command(const PipelineConfig& config) {
    std::cout << "\n=== CLASSIFYING READS ===" << std::endl;
    std::cout << "Database: " << config.database_dir << std::endl;
    std::cout << "Reads file: " << config.reads_file << std::endl;
    if (!config.reads_file2.empty()) {
        std::cout << "Reads file 2: " << config.reads_file2 << std::endl;
        std::cout << "Mode: Paired-end" << std::endl;
    } else {
        std::cout << "Mode: Single-end" << std::endl;
    }
    std::cout << "Output: " << config.output_path << std::endl;
    std::cout << "Batch size: " << config.batch_size << std::endl;
    std::cout << "Confidence threshold: " << config.classifier_params.confidence_threshold << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize classifier
        PairedEndGPUKrakenClassifier classifier(config.classifier_params);
        
        // Load database
        if (!classifier.load_database(config.database_dir)) {
            std::cerr << "Failed to load database from " << config.database_dir << std::endl;
            return false;
        }
        
        // Check if paired-end mode
        bool is_paired = !config.reads_file2.empty();
        
        if (is_paired) {
            // Load paired-end reads
            std::vector<PairedRead> paired_reads;
            std::ifstream fastq_file1(config.reads_file);
            std::ifstream fastq_file2(config.reads_file2);
            
            if (!fastq_file1.is_open()) {
                std::cerr << "Cannot open reads file: " << config.reads_file << std::endl;
                return false;
            }
            if (!fastq_file2.is_open()) {
                std::cerr << "Cannot open reads file 2: " << config.reads_file2 << std::endl;
                return false;
            }
            
            std::string line1, line2;
            int line_count = 0;
            while (std::getline(fastq_file1, line1) && std::getline(fastq_file2, line2)) {
                line_count++;
                if (line_count % 4 == 2) {  // Sequence line
                    if (!line1.empty() && !line2.empty() && 
                        line1.length() >= config.classifier_params.k && 
                        line2.length() >= config.classifier_params.k) {
                        paired_reads.push_back(PairedRead(line1, line2, "read_" + std::to_string(paired_reads.size())));
                    }
                }
                
                // Limit for testing
                if (paired_reads.size() >= 1000000) {
                    std::cout << "Limited to " << paired_reads.size() << " read pairs for processing" << std::endl;
                    break;
                }
            }
            fastq_file1.close();
            fastq_file2.close();
            
            if (paired_reads.empty()) {
                std::cerr << "No valid read pairs found" << std::endl;
                return false;
            }
            
            std::cout << "Loaded " << paired_reads.size() << " read pairs for classification" << std::endl;
            
            // Classify paired-end reads
            auto results = classifier.classify_paired_reads(paired_reads);
            
            // Write results
            std::ofstream output(config.output_path);
            if (!output.is_open()) {
                std::cerr << "Cannot create output file: " << config.output_path << std::endl;
                return false;
            }
            
            // Kraken2-style output for paired-end
            output << "# GPU Kraken Classification Results (Paired-end)\n";
            for (size_t i = 0; i < results.size(); i++) {
                const auto& result = results[i];
                char status = (result.taxon_id > 0) ? 'C' : 'U';
                
                output << status << "\t"
                       << paired_reads[i].read_id << "\t"
                       << result.taxon_id << "\t"
                       << paired_reads[i].read1.length() << "|" << paired_reads[i].read2.length() << "\t"
                       << std::fixed << std::setprecision(3) << result.confidence_score << "\t"
                       << result.read1_votes << "|" << result.read2_votes << "\t" 
                       << result.read1_kmers << "|" << result.read2_kmers << "\n";
            }
            output.close();
            
            // Generate summary report if requested
            if (!config.report_file.empty()) {
                std::ofstream report(config.report_file);
                
                // Calculate statistics
                int classified = 0;
                int unclassified = 0;
                std::map<uint32_t, int> taxon_counts;
                
                for (const auto& result : results) {
                    if (result.taxon_id > 0) {
                        classified++;
                        taxon_counts[result.taxon_id]++;
                    } else {
                        unclassified++;
                    }
                }
                
                report << "Classification Summary Report (Paired-end)\n";
                report << "==========================================\n";
                report << "Total read pairs: " << results.size() << "\n";
                report << "Classified: " << classified << " (" 
                       << std::fixed << std::setprecision(1) 
                       << (100.0 * classified / results.size()) << "%)\n";
                report << "Unclassified: " << unclassified << " (" 
                       << std::fixed << std::setprecision(1) 
                       << (100.0 * unclassified / results.size()) << "%)\n";
                report << "Unique taxa detected: " << taxon_counts.size() << "\n\n";
                
                report << "Top Taxa by Read Count:\n";
                report << "Taxon ID\tRead Count\tPercentage\n";
                
                // Sort taxa by count
                std::vector<std::pair<uint32_t, int>> sorted_taxa;
                for (const auto& pair : taxon_counts) {
                    sorted_taxa.push_back(pair);
                }
                std::sort(sorted_taxa.begin(), sorted_taxa.end(),
                         [](const auto& a, const auto& b) { return a.second > b.second; });
                
                for (size_t i = 0; i < std::min(size_t(20), sorted_taxa.size()); i++) {
                    const auto& pair = sorted_taxa[i];
                    report << pair.first << "\t" << pair.second << "\t"
                           << std::fixed << std::setprecision(2)
                           << (100.0 * pair.second / classified) << "%\n";
                }
                
                report.close();
                std::cout << "Summary report written to: " << config.report_file << std::endl;
            }
        } else {
            // Single-end mode - original code
            std::vector<std::string> reads;
            std::ifstream fastq_file(config.reads_file);
            if (!fastq_file.is_open()) {
                std::cerr << "Cannot open reads file: " << config.reads_file << std::endl;
                return false;
            }
            
            std::string line;
            int line_count = 0;
            while (std::getline(fastq_file, line)) {
                line_count++;
                if (line_count % 4 == 2) {  // Sequence line
                    if (!line.empty() && line.length() >= config.classifier_params.k) {
                        reads.push_back(line);
                    }
                }
                
                // Limit for testing
                if (reads.size() >= 1000000) {
                    std::cout << "Limited to " << reads.size() << " reads for processing" << std::endl;
                    break;
                }
            }
            fastq_file.close();
            
            if (reads.empty()) {
                std::cerr << "No valid reads found in " << config.reads_file << std::endl;
                return false;
            }
            
            std::cout << "Loaded " << reads.size() << " reads for classification" << std::endl;
            
            // Classify reads (single-end)
            auto results = classifier.classify_reads(reads);
            
            // Write results
            std::ofstream output(config.output_path);
            if (!output.is_open()) {
                std::cerr << "Cannot create output file: " << config.output_path << std::endl;
                return false;
            }
            
            // Kraken2-style output
            output << "# GPU Kraken Classification Results\n";
            for (size_t i = 0; i < results.size(); i++) {
                const auto& result = results[i];
                char status = (result.taxon_id > 0) ? 'C' : 'U';
                
                output << status << "\t"
                       << "read_" << i << "\t"
                       << result.taxon_id << "\t"
                       << reads[i].length() << "\t"
                       << std::fixed << std::setprecision(3) << result.confidence_score << "\t"
                       << result.read1_votes << "\t" << result.read1_kmers << "\n";
            }
            output.close();
            
            // Generate summary report if requested
            if (!config.report_file.empty()) {
                std::ofstream report(config.report_file);
                
                // Calculate statistics
                int classified = 0;
                int unclassified = 0;
                std::map<uint32_t, int> taxon_counts;
                
                for (const auto& result : results) {
                    if (result.taxon_id > 0) {
                        classified++;
                        taxon_counts[result.taxon_id]++;
                    } else {
                        unclassified++;
                    }
                }
                
                report << "Classification Summary Report\n";
                report << "============================\n";
                report << "Total reads: " << results.size() << "\n";
                report << "Classified: " << classified << " (" 
                       << std::fixed << std::setprecision(1) 
                       << (100.0 * classified / results.size()) << "%)\n";
                report << "Unclassified: " << unclassified << " (" 
                       << std::fixed << std::setprecision(1) 
                       << (100.0 * unclassified / results.size()) << "%)\n";
                report << "Unique taxa detected: " << taxon_counts.size() << "\n\n";
                
                report << "Top Taxa by Read Count:\n";
                report << "Taxon ID\tRead Count\tPercentage\n";
                
                // Sort taxa by count
                std::vector<std::pair<uint32_t, int>> sorted_taxa;
                for (const auto& pair : taxon_counts) {
                    sorted_taxa.push_back(pair);
                }
                std::sort(sorted_taxa.begin(), sorted_taxa.end(),
                         [](const auto& a, const auto& b) { return a.second > b.second; });
                
                for (size_t i = 0; i < std::min(size_t(20), sorted_taxa.size()); i++) {
                    const auto& pair = sorted_taxa[i];
                    report << pair.first << "\t" << pair.second << "\t"
                           << std::fixed << std::setprecision(2)
                           << (100.0 * pair.second / classified) << "%\n";
                }
                
                report.close();
                std::cout << "Summary report written to: " << config.report_file << std::endl;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\n✓ Classification completed successfully in " 
                  << duration.count() << " seconds" << std::endl;
        std::cout << "Results written to: " << config.output_path << std::endl;
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during classification: " << e.what() << std::endl;
        return false;
    }
}

bool pipeline_command(const PipelineConfig& config) {
    std::cout << "\n=== RUNNING COMPLETE PIPELINE ===" << std::endl;
    
    // Create output directory structure
    std::filesystem::create_directories(config.output_path);
    std::string db_dir = config.output_path + "/database";
    std::string results_file = config.output_path + "/classification_results.txt";
    std::string report_file = config.output_path + "/summary_report.txt";
    
    // Step 1: Build database
    PipelineConfig build_config = config;
    build_config.command = "build";
    build_config.output_path = db_dir;
    
    if (!build_database_command(build_config)) {
        std::cerr << "Pipeline failed at database build step" << std::endl;
        return false;
    }
    
    // Step 2: Classify reads
    PipelineConfig classify_config = config;
    classify_config.command = "classify";
    classify_config.database_dir = db_dir;
    classify_config.output_path = results_file;
    classify_config.report_file = report_file;
    
    if (!classify_reads_command(classify_config)) {
        std::cerr << "Pipeline failed at classification step" << std::endl;
        return false;
    }
    
    std::cout << "\n✓ Complete pipeline finished successfully!" << std::endl;
    std::cout << "Results available in: " << config.output_path << std::endl;
    
    return true;
}

int main(int argc, char* argv[]) {
    std::cout << "GPU-Accelerated Kraken2-Style Taxonomic Classifier" << std::endl;
    std::cout << "==================================================" << std::endl;
    
    PipelineConfig config;
    
    if (!parse_arguments(argc, argv, config)) {
        print_usage(argv[0]);
        return 1;
    }
    
    if (!validate_config(config)) {
        return 1;
    }
    
    // Check CUDA availability
    int device_count;
    cudaError_t cuda_status = cudaGetDeviceCount(&device_count);
    if (cuda_status != cudaSuccess || device_count == 0) {
        std::cerr << "Error: No CUDA-capable GPU found!" << std::endl;
        return 1;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using GPU: " << prop.name << " with " 
              << prop.totalGlobalMem / (1024*1024*1024) << " GB memory" << std::endl;
    
    bool success = false;
    
    try {
        if (config.command == "build") {
            success = build_database_command(config);
        } else if (config.command == "classify") {
            success = classify_reads_command(config);
        } else if (config.command == "pipeline") {
            success = pipeline_command(config);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
    
    if (success) {
        std::cout << "\n🎉 Operation completed successfully!" << std::endl;
        return 0;
    } else {
        std::cerr << "\n❌ Operation failed!" << std::endl;
        return 1;
    }
}

// Compilation instructions:
/*
To compile this complete pipeline:

1. Create a Makefile:

```makefile
NVCC = nvcc
CXX = g++
CUDA_FLAGS = -std=c++14 -O3 -arch=sm_70 -lineinfo
CXX_FLAGS = -std=c++14 -O3

INCLUDES = -I/usr/local/cuda/include
LIBS = -lcuda -lcudart -lz

SOURCES = kraken_pipeline_main.cpp gpu_kraken_classifier.cu gpu_kraken_database_builder.cu
TARGET = gpu_kraken_pipeline

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(NVCC) $(CUDA_FLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)

install:
	cp $(TARGET) /usr/local/bin/

.PHONY: all clean install
```

2. Compile:
```bash
make
```

3. Run examples:
```bash
# Build database from RefSeq genomes
./gpu_kraken_pipeline build --genome-dir ./refseq_genomes --output ./my_kraken_db

# Classify reads
./gpu_kraken_pipeline classify --database ./my_kraken_db --reads sample.fastq.gz --output results.txt --report summary.txt

# Complete pipeline
./gpu_kraken_pipeline pipeline --genome-dir ./genomes --reads sample.fastq.gz --output ./results --confidence 0.1
```
*/