// kraken_pipeline_main.cpp
// Complete front-to-back GPU-accelerated Kraken2-style pipeline
// Build database from genomes + classify reads

// Include the headers
#include "gpu_kraken_classifier.h"
#include "gpu_kraken_database_builder.h"
#include "sample_csv_parser.h"
#include "classification_report_generator.h"

// ADD THESE NEW INCLUDES:
#include "fast_enhanced_classifier.h"              // The fast enhanced classifier
#include "phase1_enhanced_classifier_with_phylo.h" // The phylo-aware classifier
#include "ncbi_taxonomy_loader.h"                  // NCBI taxonomy support

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono>
#include <zlib.h>
#include <sys/stat.h>  // for stat()
#include <cstdio>      // for fopen()
#include <cstring>     // for strcasecmp()

// Comprehensive heap corruption debugging
#include <signal.h>
#include <execinfo.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

// Signal handler for catching crashes
void crash_handler(int sig) {
    void *array[10];
    size_t size;
    
    // Get void*'s for all entries on the stack
    size = backtrace(array, 10);
    
    // Print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

// Memory debugging helper class
class HeapValidator {
private:
    static bool heap_check_enabled;
    static size_t allocation_count;
    
public:
    static void enable_validation() {
        heap_check_enabled = true;
        allocation_count = 0;
        std::cout << "HEAP DEBUG: Validation enabled" << std::endl;
    }
    
    static bool validate_heap(const std::string& location) {
        if (!heap_check_enabled) return true;
        
        std::cout << "HEAP DEBUG: Checking at " << location << " (alloc count: " << allocation_count << ")" << std::endl;
        
        // Try a small allocation and deallocation
        try {
            void* test_ptr = malloc(16);
            if (!test_ptr) {
                std::cerr << "HEAP ERROR: malloc(16) failed at " << location << std::endl;
                return false;
            }
            
            // Write some data to ensure it's valid
            memset(test_ptr, 0xAA, 16);
            
            // Free it
            free(test_ptr);
            
            allocation_count++;
            
            if (allocation_count % 100 == 0) {
                std::cout << "HEAP DEBUG: " << allocation_count << " successful heap operations" << std::endl;
            }
            
            return true;
            
        } catch (...) {
            std::cerr << "HEAP ERROR: Exception during heap test at " << location << std::endl;
            return false;
        }
    }
    
    static void check_string_safety(const std::string& str, const std::string& name) {
        std::cout << "STRING DEBUG: Checking '" << name << "' (length: " << str.length() << ")" << std::endl;
        
        if (str.length() > 10000) {
            std::cerr << "STRING WARNING: Very long string '" << name << "': " << str.length() << " chars" << std::endl;
        }
        
        // Check for embedded nulls
        for (size_t i = 0; i < str.length(); i++) {
            if (str[i] == '\0') {
                std::cerr << "STRING ERROR: Embedded null in '" << name << "' at position " << i << std::endl;
                break;
            }
        }
        
        // Check string capacity vs size
        if (str.capacity() > str.size() * 10) {
            std::cerr << "STRING WARNING: Large capacity waste in '" << name << "': " 
                      << str.capacity() << " vs " << str.size() << std::endl;
        }
    }
};

bool HeapValidator::heap_check_enabled = false;
size_t HeapValidator::allocation_count = 0;

// Enhanced classification mode selection
enum class ClassificationMode {
    STANDARD,           // Original paired-end GPU classifier
    FAST_ENHANCED,      // Fast enhanced with compact taxonomy
    PHYLO_ENHANCED      // Full phylogenetic validation
};

void print_usage(const char* program_name) {
    std::cout << "GPU-Accelerated Kraken2-Style Pipeline (High Capacity)" << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir <dir> --output <db_dir> [options]" << std::endl;
    std::cout << "  " << program_name << " classify --database <db_dir> --reads <fastq> [options]" << std::endl;
    std::cout << "  " << program_name << " pipeline --genome-dir <dir> --reads <fastq> --output <dir> [options]" << std::endl;
    std::cout << "  " << program_name << " batch --csv <samples.csv> --database <db_dir> --batch-output <dir> [options]" << std::endl;
    std::cout << "\nCommands:" << std::endl;
    std::cout << "  build      Build database from genome files" << std::endl;
    std::cout << "  classify   Classify reads using existing database" << std::endl;
    std::cout << "  pipeline   Complete pipeline: build database + classify reads" << std::endl;
    std::cout << "  batch      Process multiple samples from CSV file" << std::endl;
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
    
    // NEW: Memory and capacity options
    std::cout << "\nMemory and Capacity Options:" << std::endl;
    std::cout << "  --minimizer-capacity <int>  Max minimizers per batch (default: 5000000)" << std::endl;
    std::cout << "  --auto-memory               Enable automatic memory scaling (default)" << std::endl;
    std::cout << "  --no-auto-memory            Disable automatic memory scaling" << std::endl;
    std::cout << "  --memory-fraction <int>     GPU memory usage % when auto-scaling (default: 80)" << std::endl;
    
    std::cout << "\nOther Options:" << std::endl;
    std::cout << "  --taxonomy <dir>      NCBI taxonomy directory (optional)" << std::endl;
    std::cout << "  --threads <int>       CPU threads for I/O (default: 4)" << std::endl;
    std::cout << "  --quick               Quick mode: stop at first hit" << std::endl;
    std::cout << "  --report <file>       Generate summary report" << std::endl;
    
    std::cout << "\nStreaming Options (for large datasets):" << std::endl;
    std::cout << "  --streaming-batch <int>    Read pairs per batch (default: 100000)" << std::endl;
    std::cout << "  --force-streaming          Force streaming mode regardless of file size" << std::endl;
    std::cout << "  --disable-streaming        Disable automatic streaming" << std::endl;
    std::cout << "  --gpu-memory-limit <mb>    GPU memory limit for reads (default: 8192)" << std::endl;
    std::cout << "\nStreaming Examples:" << std::endl;
    std::cout << "  # Large paired-end dataset with custom batch size" << std::endl;
    std::cout << "  " << program_name << " classify --database ./db --reads large_R1.fastq.gz --reads2 large_R2.fastq.gz --output results.txt --streaming-batch 25000" << std::endl;
    
    std::cout << "\nBatch Processing Options:" << std::endl;
    std::cout << "  --csv <file>              CSV file with sample information" << std::endl;
    std::cout << "  --batch-output <dir>      Output directory for batch results (default: batch_results)" << std::endl;
    std::cout << "  --no-sample-dirs          Don't create individual sample directories" << std::endl;
    std::cout << "  --stop-on-error           Stop batch processing on first error" << std::endl;
    std::cout << "  --no-validate-paths       Skip file existence validation on CSV load" << std::endl;
    
    std::cout << "\nCSV File Format:" << std::endl;
    std::cout << "Required columns (case-insensitive, flexible names):" << std::endl;
    std::cout << "  - Sample Name/ID: sample_name, Sample, ID, etc." << std::endl;
    std::cout << "  - File Path: file_path, path, directory, etc." << std::endl;
    std::cout << "  - R1 File: R1, read1, fastq1, etc." << std::endl;
    std::cout << "  - R2 File: R2, read2, fastq2, etc. (optional for single-end)" << std::endl;
    
    std::cout << "\nBatch Processing Examples:" << std::endl;
    std::cout << "  # Process multiple samples from CSV" << std::endl;
    std::cout << "  " << program_name << " batch --csv samples.csv --database ./db --batch-output ./results" << std::endl;
    std::cout << "\n  # Example CSV content:" << std::endl;
    std::cout << "  Sample Name,File Path,R1 file,R2 file" << std::endl;
    std::cout << "  Sample001,/data/fastq/,sample001_R1.fastq.gz,sample001_R2.fastq.gz" << std::endl;
    std::cout << "  Sample002,/data/fastq/,sample002_R1.fastq.gz,sample002_R2.fastq.gz" << std::endl;
    
    std::cout << "\nMemory Usage Examples:" << std::endl;
    std::cout << "  # Use 10M minimizers per batch (high memory, comprehensive)" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir ./genomes --output ./db --minimizer-capacity 10000000" << std::endl;
    std::cout << "\n  # Conservative memory usage (2M minimizers)" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir ./genomes --output ./db --minimizer-capacity 2000000" << std::endl;
    std::cout << "\n  # Let system auto-scale to use 90% of GPU memory" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir ./genomes --output ./db --auto-memory --memory-fraction 90" << std::endl;
    std::cout << "\n  # Manual control for specific hardware" << std::endl;
    std::cout << "  " << program_name << " build --genome-dir ./genomes --output ./db --no-auto-memory --minimizer-capacity 8000000" << std::endl;
    
    std::cout << "\nEnhanced Classification Options:" << std::endl;
    std::cout << "  --enhanced                     Enable enhanced classification features" << std::endl;
    std::cout << "  --classification-mode <mode>   Classification mode:" << std::endl;
    std::cout << "                                   standard        - Original GPU classifier" << std::endl;
    std::cout << "                                   fast-enhanced   - Fast enhanced with compact taxonomy" << std::endl;
    std::cout << "                                   phylo-enhanced  - Full phylogenetic validation" << std::endl;
    std::cout << "  --compact-taxonomy <file>      Compact taxonomy file (for fast-enhanced)" << std::endl;
    std::cout << "  --taxonomy-nodes <file>        NCBI nodes.dmp file (for phylo-enhanced)" << std::endl;
    std::cout << "  --taxonomy-names <file>        NCBI names.dmp file (for phylo-enhanced)" << std::endl;
    std::cout << "  --phylo-quality <level>        Phylogenetic quality: fast|balanced|accurate" << std::endl;
    
    std::cout << "\nReport Generation Options:" << std::endl;
    std::cout << "  --no-reports                   Disable automatic report generation" << std::endl;
    std::cout << "  --reports-mpa-only             Generate only MetaPhlAn-style profile" << std::endl;
    std::cout << "  --reports-species-only         Generate only species count table" << std::endl;
    std::cout << "  --reports-genus-only           Generate only genus count table" << std::endl;
    std::cout << "  --reports-min-abundance <f>    Minimum abundance threshold % (default: 0.01)" << std::endl;
    std::cout << "  --reports-include-unclassified Include unclassified reads in reports" << std::endl;
    std::cout << "  --reports-sample-name <name>   Custom sample name (default: auto-detect)" << std::endl;
    
    std::cout << "\nEnhanced Classification Examples:" << std::endl;
    std::cout << "  # Fast enhanced mode with compact taxonomy" << std::endl;
    std::cout << "  " << program_name << " classify --database ./db --reads sample.fastq.gz \\" << std::endl;
    std::cout << "    --classification-mode fast-enhanced --compact-taxonomy compact_tax.bin --output results.txt" << std::endl;
    std::cout << "\n  # Full phylogenetic validation mode" << std::endl;
    std::cout << "  " << program_name << " classify --database ./db --reads sample.fastq.gz \\" << std::endl;
    std::cout << "    --classification-mode phylo-enhanced --taxonomy-nodes nodes.dmp \\" << std::endl;
    std::cout << "    --taxonomy-names names.dmp --output results.txt" << std::endl;
    std::cout << "\n  # Build compact taxonomy for fast enhanced mode" << std::endl;
    std::cout << "  build_compact_taxonomy --nodes nodes.dmp --names names.dmp --output compact_taxonomy.bin" << std::endl;
}

struct StreamingConfig {
    size_t read_batch_size = 100000;          // Process 100K read pairs per batch
    size_t max_gpu_memory_for_reads_mb = 8192; // Reserve 8GB for read processing
    bool enable_streaming = true;
    bool show_batch_progress = true;
    size_t max_reads_in_memory = 1000000;     // Trigger streaming if more reads
};

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
    int gpu_batch_size = 100;  // GPU batch size for database building
    
    // NEW: Minimizer capacity controls
    int minimizer_capacity = 5000000;     // Default to 5M (up from 1M)
    bool auto_memory_scaling = true;      // Enable auto-scaling by default
    int memory_fraction = 80;             // Use 80% of GPU memory
    
    int cpu_threads = 4;
    bool quick_mode = false;
    bool verbose = false;
    
    // Streaming configuration
    StreamingConfig streaming_config;
    bool force_streaming = false;
    size_t batch_size_override = 0;
    
    // CSV batch processing fields
    std::string csv_file;
    std::string batch_output_dir = "batch_results";
    bool create_sample_dirs = true;
    bool stop_on_error = false;
    bool validate_paths_on_load = true;
    
    // NEW: Report generation options
    bool generate_reports = true;           // Generate comprehensive reports
    bool reports_only_mpa = false;         // Generate only MetaPhlAn-style
    bool reports_only_species = false;     // Generate only species counts
    bool reports_only_genus = false;       // Generate only genus counts
    float reports_min_abundance = 0.01f;   // Minimum abundance threshold (%)
    bool reports_include_unclassified = false;
    std::string reports_sample_name = "";  // Auto-detect from input file if empty
    
    // ADD THESE NEW FIELDS:
    ClassificationMode classification_mode = ClassificationMode::STANDARD;
    std::string compact_taxonomy_path;                    // For fast enhanced mode
    std::string taxonomy_nodes_path;                      // For phylo enhanced mode  
    std::string taxonomy_names_path;                      // For phylo enhanced mode
    bool enable_enhanced_features = false;               // Enable enhanced classification
    FastEnhancedParams::PhyloQualityLevel phylo_quality = FastEnhancedParams::BALANCED;
};

bool parse_arguments(int argc, char* argv[], PipelineConfig& config) {
    if (argc < 2) {
        return false;
    }
    
    config.command = argv[1];
    
    if (config.command != "build" && config.command != "classify" && 
        config.command != "pipeline" && config.command != "batch") {
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
        // NEW: Minimizer capacity control
        } else if (arg == "--minimizer-capacity" && i + 1 < argc) {
            config.minimizer_capacity = std::stoi(argv[++i]);
        } else if (arg == "--auto-memory") {
            config.auto_memory_scaling = true;
        } else if (arg == "--no-auto-memory") {
            config.auto_memory_scaling = false;
        } else if (arg == "--memory-fraction" && i + 1 < argc) {
            config.memory_fraction = std::stoi(argv[++i]);
        } else if (arg == "--quick") {
            config.quick_mode = true;
        } else if (arg == "--verbose") {
            config.verbose = true;
        // Streaming options
        } else if (arg == "--streaming-batch" && i + 1 < argc) {
            config.streaming_config.read_batch_size = std::stoi(argv[++i]);
        } else if (arg == "--force-streaming") {
            config.force_streaming = true;
        } else if (arg == "--gpu-memory-limit" && i + 1 < argc) {
            config.streaming_config.max_gpu_memory_for_reads_mb = std::stoi(argv[++i]);
        } else if (arg == "--disable-streaming") {
            config.streaming_config.enable_streaming = false;
        // CSV batch processing arguments
        } else if (arg == "--csv" && i + 1 < argc) {
            config.csv_file = argv[++i];
        } else if (arg == "--batch-output" && i + 1 < argc) {
            config.batch_output_dir = argv[++i];
        } else if (arg == "--no-sample-dirs") {
            config.create_sample_dirs = false;
        } else if (arg == "--stop-on-error") {
            config.stop_on_error = true;
        } else if (arg == "--no-validate-paths") {
            config.validate_paths_on_load = false;
        } else if (arg == "--no-reports") {
            config.generate_reports = false;
        } else if (arg == "--reports-mpa-only") {
            config.reports_only_mpa = true;
        } else if (arg == "--reports-species-only") {
            config.reports_only_species = true;
        } else if (arg == "--reports-genus-only") {
            config.reports_only_genus = true;
        } else if (arg == "--reports-min-abundance" && i + 1 < argc) {
            config.reports_min_abundance = std::stof(argv[++i]);
        } else if (arg == "--reports-include-unclassified") {
            config.reports_include_unclassified = true;
        } else if (arg == "--reports-sample-name" && i + 1 < argc) {
            config.reports_sample_name = argv[++i];
        // Enhanced classification options
        } else if (arg == "--enhanced") {
            config.enable_enhanced_features = true;
        } else if (arg == "--classification-mode" && i + 1 < argc) {
            std::string mode = argv[++i];
            if (mode == "standard") {
                config.classification_mode = ClassificationMode::STANDARD;
            } else if (mode == "fast-enhanced") {
                config.classification_mode = ClassificationMode::FAST_ENHANCED;
                config.enable_enhanced_features = true;
            } else if (mode == "phylo-enhanced") {
                config.classification_mode = ClassificationMode::PHYLO_ENHANCED;
                config.enable_enhanced_features = true;
            } else {
                std::cerr << "Unknown classification mode: " << mode << std::endl;
                return false;
            }
        } else if (arg == "--compact-taxonomy" && i + 1 < argc) {
            config.compact_taxonomy_path = argv[++i];
        } else if (arg == "--taxonomy-nodes" && i + 1 < argc) {
            config.taxonomy_nodes_path = argv[++i];
        } else if (arg == "--taxonomy-names" && i + 1 < argc) {
            config.taxonomy_names_path = argv[++i];
        } else if (arg == "--phylo-quality" && i + 1 < argc) {
            std::string quality = argv[++i];
            if (quality == "fast") {
                config.phylo_quality = FastEnhancedParams::FAST_APPROXIMATE;
            } else if (quality == "balanced") {
                config.phylo_quality = FastEnhancedParams::BALANCED;
            } else if (quality == "accurate") {
                config.phylo_quality = FastEnhancedParams::HIGH_ACCURACY;
            }
        } else {
            std::cerr << "Warning: Unknown argument '" << arg << "'" << std::endl;
        }
    }
    
    return true;
}

// FIXED: Safe version of is_concatenated_fna_file function
// This prevents the malloc error by avoiding unsafe std::filesystem operations
bool is_concatenated_fna_file(const std::string& path) {
    std::cout << "DEBUG: Conservative check if path is concatenated FNA file: '" << path << "'" << std::endl;
    
    // Basic validation
    if (path.empty() || path.size() > 4096) {
        std::cout << "DEBUG: Path validation failed" << std::endl;
        return false;
    }
    
    // Use stat() to check if it's a regular file (most reliable cross-platform way)
    struct stat file_stat;
    if (stat(path.c_str(), &file_stat) != 0) {
        std::cout << "DEBUG: stat() failed - file does not exist" << std::endl;
        return false;
    }
    
    // Check if it's a regular file
    if (!S_ISREG(file_stat.st_mode)) {
        std::cout << "DEBUG: Path is not a regular file" << std::endl;
        return false;
    }
    
    std::cout << "DEBUG: Path is a regular file" << std::endl;
    
    // Manual extension extraction
    size_t last_dot = path.find_last_of('.');
    size_t last_slash = path.find_last_of('/');
    size_t last_backslash = path.find_last_of('\\');
    
    // Make sure the dot is after the last path separator
    size_t last_separator = std::max(
        last_slash == std::string::npos ? 0 : last_slash,
        last_backslash == std::string::npos ? 0 : last_backslash
    );
    
    if (last_dot == std::string::npos || last_dot <= last_separator) {
        std::cout << "DEBUG: No valid extension found" << std::endl;
        return false;
    }
    
    std::string ext = path.substr(last_dot);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    std::cout << "DEBUG: Extracted extension: '" << ext << "'" << std::endl;
    
    bool result = (ext == ".fna" || ext == ".fasta" || ext == ".fa");
    std::cout << "DEBUG: Is FNA file: " << (result ? "true" : "false") << std::endl;
    
    return result;
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
    
    if (config.command == "batch") {
        if (config.csv_file.empty()) {
            std::cerr << "Error: --csv is required for batch mode" << std::endl;
            return false;
        }
        if (!std::filesystem::exists(config.csv_file)) {
            std::cerr << "Error: CSV file does not exist: " << config.csv_file << std::endl;
            return false;
        }
        if (config.database_dir.empty()) {
            std::cerr << "Error: --database is required for batch mode" << std::endl;
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
    std::cout << "\n=== BUILDING KRAKEN DATABASE (High Capacity) ===" << std::endl;
    std::cout << "Genome directory: " << config.genome_dir << std::endl;
    std::cout << "Output directory: " << config.output_path << std::endl;
    std::cout << "Parameters: k=" << config.classifier_params.k 
              << ", ell=" << config.classifier_params.ell 
              << ", spaces=" << config.classifier_params.spaces << std::endl;
    
    // Validate input paths before proceeding
    if (config.genome_dir.empty()) {
        std::cerr << "ERROR: Genome directory path is empty!" << std::endl;
        return false;
    }
    
    if (config.output_path.empty()) {
        std::cerr << "ERROR: Output path is empty!" << std::endl;
        return false;
    }
    
    // Check if genome directory exists
    if (!std::filesystem::exists(config.genome_dir)) {
        std::cerr << "ERROR: Genome directory does not exist: " << config.genome_dir << std::endl;
        return false;
    }
    
    // NEW: Display capacity settings
    std::cout << "Capacity settings:" << std::endl;
    if (config.auto_memory_scaling) {
        std::cout << "  Auto memory scaling: ENABLED (" << config.memory_fraction << "% of GPU memory)" << std::endl;
    } else {
        std::cout << "  Auto memory scaling: DISABLED" << std::endl;
        std::cout << "  Manual minimizer capacity: " << config.minimizer_capacity << std::endl;
    }
    
    // Create output directory if it doesn't exist
    try {
        std::filesystem::create_directories(config.output_path);
        std::cout << "Created/verified output directory: " << config.output_path << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create output directory: " << e.what() << std::endl;
        return false;
    }
    
    // Initialize CUDA runtime early
    std::cout << "Initializing CUDA runtime..." << std::endl;
    cudaError_t init_status = cudaFree(0);
    if (init_status != cudaSuccess) {
        std::cerr << "WARNING: Early CUDA initialization returned: " << cudaGetErrorString(init_status) << std::endl;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        std::cout << "Creating GPUKrakenDatabaseBuilder..." << std::endl;
        
        // Create the builder with a copy of the classification params to avoid reference issues
        ClassificationParams params_copy = config.classifier_params;
        std::unique_ptr<GPUKrakenDatabaseBuilder> builder;
        
        try {
            builder = std::make_unique<GPUKrakenDatabaseBuilder>(config.output_path, params_copy);
        } catch (const std::exception& e) {
            std::cerr << "Failed to create database builder: " << e.what() << std::endl;
            return false;
        }
        
        std::cout << "Database builder created successfully!" << std::endl;
        
        // Configure memory and capacity settings
        if (config.auto_memory_scaling) {
            std::cout << "Enabling auto memory scaling..." << std::endl;
            builder->enable_auto_memory_scaling(true, config.memory_fraction);
        } else {
            std::cout << "Using manual memory settings..." << std::endl;
            builder->enable_auto_memory_scaling(false);
            builder->set_minimizer_capacity(config.minimizer_capacity);
        }
        
        // Set batch size if specified
        if (config.gpu_batch_size > 0) {
            std::cout << "Setting batch size to " << config.gpu_batch_size << std::endl;
            builder->set_batch_size(config.gpu_batch_size);
        }
        
        std::cout << "Configuration complete, starting database build..." << std::endl;
        std::cout << "Input path: '" << config.genome_dir << "'" << std::endl;
        std::cout << "Taxonomy path: '" << config.taxonomy_dir << "'" << std::endl;
        
        // Debug: Check string integrity before is_concatenated_fna_file
        std::cout << "DEBUG: Checking genome_dir string integrity..." << std::endl;
        std::cout << "DEBUG: genome_dir size = " << config.genome_dir.size() << std::endl;
        for (size_t i = 0; i < config.genome_dir.size(); i++) {
            if (config.genome_dir[i] == '\0') {
                std::cerr << "ERROR: Null character found in genome_dir at position " << i << std::endl;
                return false;
            }
        }
        std::cout << "DEBUG: genome_dir string integrity OK" << std::endl;
        
        bool success = false;
        
        std::cout << "DEBUG: About to call is_concatenated_fna_file..." << std::endl;
        
        // Check if input is a concatenated FNA file or a directory
        bool is_fna_file = false;
        try {
            is_fna_file = is_concatenated_fna_file(config.genome_dir);
            std::cout << "DEBUG: is_concatenated_fna_file returned " << is_fna_file << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "ERROR in is_concatenated_fna_file: " << e.what() << std::endl;
            return false;
        }
        
        if (is_fna_file) {
            std::cout << "Detected concatenated FNA file - using streaming processor" << std::endl;
            success = builder->build_database_from_streaming_fna(config.genome_dir, config.taxonomy_dir);
        } else {
            std::cout << "Processing genome directory with individual files" << std::endl;
            success = builder->build_database_from_genomes(config.genome_dir, config.taxonomy_dir);
        }
        
        if (success) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
            
            std::cout << "\n✓ Database build completed successfully in " 
                      << duration.count() << " seconds" << std::endl;
            std::cout << "Database saved to: " << config.output_path << std::endl;
            
            // Explicitly reset builder before returning
            builder.reset();
            return true;
        } else {
            std::cerr << "Database build failed!" << std::endl;
            builder.reset();
            return false;
        }
        
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation error: " << e.what() << std::endl;
        return false;
    } catch (const std::runtime_error& e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "Error during database build: " << e.what() << std::endl;
        return false;
    } catch (...) {
        std::cerr << "Unknown error during database build" << std::endl;
        return false;
    }
}

// Streaming FASTQ reader class
class StreamingFASTQReader {
private:
    std::ifstream file1, file2;
    gzFile gz_file1 = nullptr, gz_file2 = nullptr;
    bool is_paired;
    bool is_gzipped;
    size_t current_read_count = 0;
    
public:
    StreamingFASTQReader(const std::string& file1_path, const std::string& file2_path = "") {
        is_paired = !file2_path.empty();
        is_gzipped = file1_path.substr(file1_path.find_last_of(".") + 1) == "gz";
        
        if (is_gzipped) {
            gz_file1 = gzopen(file1_path.c_str(), "rb");
            if (is_paired) {
                gz_file2 = gzopen(file2_path.c_str(), "rb");
            }
        } else {
            file1.open(file1_path);
            if (is_paired) {
                file2.open(file2_path);
            }
        }
    }
    
    ~StreamingFASTQReader() {
        if (gz_file1) gzclose(gz_file1);
        if (gz_file2) gzclose(gz_file2);
    }
    
    bool readBatch(std::vector<PairedRead>& batch, size_t batch_size) {
        batch.clear();
        batch.reserve(batch_size);
        
        std::string line1, line2;
        int line_count = 0;
        
        while (batch.size() < batch_size) {
            bool has_line1 = false, has_line2 = false;
            
            if (is_gzipped) {
                char buffer1[4096], buffer2[4096];
                has_line1 = gzgets(gz_file1, buffer1, sizeof(buffer1)) != nullptr;
                if (has_line1) {
                    line1 = buffer1;
                    if (!line1.empty() && line1.back() == '\n') line1.pop_back();
                }
                
                if (is_paired && has_line1) {
                    has_line2 = gzgets(gz_file2, buffer2, sizeof(buffer2)) != nullptr;
                    if (has_line2) {
                        line2 = buffer2;
                        if (!line2.empty() && line2.back() == '\n') line2.pop_back();
                    }
                }
            } else {
                has_line1 = static_cast<bool>(std::getline(file1, line1));
                if (is_paired && has_line1) {
                    has_line2 = static_cast<bool>(std::getline(file2, line2));
                }
            }
            
            if (!has_line1) break;
            
            line_count++;
            if (line_count % 4 == 2) {  // Sequence line
                if (is_paired && has_line2) {
                    if (line1.length() >= 35 && line2.length() >= 35) {
                        batch.emplace_back(line1, line2, "pair_" + std::to_string(current_read_count++));
                    }
                } else {
                    if (line1.length() >= 35) {
                        batch.emplace_back(line1);
                        batch.back().read_id = "single_" + std::to_string(current_read_count++);
                    }
                }
            }
        }
        
        return !batch.empty();
    }
    
    size_t getTotalReadsProcessed() const { return current_read_count; }
};

// Forward declarations
bool generate_classification_reports(const PipelineConfig& config);

bool classify_reads_with_streaming_and_classifier(const PipelineConfig& config, PairedEndGPUKrakenClassifier& classifier) {
    try {
        // Setup streaming reader
        StreamingFASTQReader reader(config.reads_file, config.reads_file2);
        
        // Open output file
        std::ofstream output(config.output_path);
        if (!output.is_open()) {
            std::cerr << "Cannot create output file: " << config.output_path << std::endl;
            return false;
        }
        
        output << "# GPU Kraken Classification Results (Streaming)\n";
        
        // Streaming classification loop
        std::vector<PairedRead> batch;
        size_t total_reads_processed = 0;
        size_t batch_number = 0;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        while (reader.readBatch(batch, config.streaming_config.read_batch_size)) {
            batch_number++;
            
            if (config.streaming_config.show_batch_progress && batch_number % 10 == 0) {
                std::cout << "    Batch " << batch_number 
                          << " (" << batch.size() << " reads)" << std::endl;
            }
            
            // Classify this batch
            auto batch_results = classifier.classify_paired_reads(batch);
            
            // Write results immediately
            for (size_t i = 0; i < batch_results.size(); i++) {
                const auto& result = batch_results[i];
                const auto& read_pair = batch[i];
                
                char status = (result.taxon_id > 0) ? 'C' : 'U';
                
                if (read_pair.is_paired) {
                    output << status << "\t"
                           << read_pair.read_id << "\t"
                           << result.taxon_id << "\t"
                           << read_pair.read1.length() << "|" << read_pair.read2.length() << "\t"
                           << std::fixed << std::setprecision(3) << result.confidence_score << "\t"
                           << result.read1_votes << "|" << result.read2_votes << "\t" 
                           << result.read1_kmers << "|" << result.read2_kmers << "\n";
                } else {
                    output << status << "\t"
                           << read_pair.read_id << "\t"
                           << result.taxon_id << "\t"
                           << read_pair.read1.length() << "\t"
                           << std::fixed << std::setprecision(3) << result.confidence_score << "\t"
                           << result.read1_votes << "\t" << result.read1_kmers << "\n";
                }
            }
            
            total_reads_processed += batch.size();
        }
        
        output.close();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "    Completed: " << total_reads_processed 
                  << " reads in " << total_duration.count() << "s" << std::endl;
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during streaming classification: " << e.what() << std::endl;
        return false;
    }
}

bool classify_reads_with_streaming(const PipelineConfig& config) {
    std::cout << "\n=== STREAMING CLASSIFICATION ===" << std::endl;
    
    try {
        // Initialize classifier
        PairedEndGPUKrakenClassifier classifier(config.classifier_params);
        
        // Load database
        if (!classifier.load_database(config.database_dir)) {
            std::cerr << "Failed to load database from " << config.database_dir << std::endl;
            return false;
        }
        
        // Setup streaming reader
        StreamingFASTQReader reader(config.reads_file, config.reads_file2);
        
        // Open output file
        std::ofstream output(config.output_path);
        if (!output.is_open()) {
            std::cerr << "Cannot create output file: " << config.output_path << std::endl;
            return false;
        }
        
        output << "# GPU Kraken Classification Results (Streaming)\n";
        
        // Streaming classification loop
        std::vector<PairedRead> batch;
        size_t total_reads_processed = 0;
        size_t batch_number = 0;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        while (reader.readBatch(batch, config.streaming_config.read_batch_size)) {
            batch_number++;
            
            if (config.streaming_config.show_batch_progress) {
                std::cout << "Processing batch " << batch_number 
                          << " (" << batch.size() << " read pairs)..." << std::endl;
            }
            
            // Classify this batch
            auto batch_results = classifier.classify_paired_reads(batch);
            
            // Write results immediately
            for (size_t i = 0; i < batch_results.size(); i++) {
                const auto& result = batch_results[i];
                const auto& read_pair = batch[i];
                
                char status = (result.taxon_id > 0) ? 'C' : 'U';
                
                if (read_pair.is_paired) {
                    output << status << "\t"
                           << read_pair.read_id << "\t"
                           << result.taxon_id << "\t"
                           << read_pair.read1.length() << "|" << read_pair.read2.length() << "\t"
                           << std::fixed << std::setprecision(3) << result.confidence_score << "\t"
                           << result.read1_votes << "|" << result.read2_votes << "\t" 
                           << result.read1_kmers << "|" << result.read2_kmers << "\n";
                } else {
                    output << status << "\t"
                           << read_pair.read_id << "\t"
                           << result.taxon_id << "\t"
                           << read_pair.read1.length() << "\t"
                           << std::fixed << std::setprecision(3) << result.confidence_score << "\t"
                           << result.read1_votes << "\t" << result.read1_kmers << "\n";
                }
            }
            
            total_reads_processed += batch.size();
            
            // Optional: Print progress
            if (batch_number % 10 == 0) {
                auto current_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);
                double reads_per_second = total_reads_processed / (duration.count() + 1);
                
                std::cout << "Processed " << total_reads_processed 
                          << " reads (" << std::fixed << std::setprecision(0) 
                          << reads_per_second << " reads/sec)" << std::endl;
            }
        }
        
        output.close();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\n✓ Streaming classification completed!" << std::endl;
        std::cout << "Total reads processed: " << total_reads_processed << std::endl;
        std::cout << "Total batches: " << batch_number << std::endl;
        std::cout << "Time: " << total_duration.count() << " seconds" << std::endl;
        std::cout << "Results written to: " << config.output_path << std::endl;
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during streaming classification: " << e.what() << std::endl;
        return false;
    }
}

// Forward declaration
bool classify_reads_enhanced(const PipelineConfig& config);

bool classify_reads_command(const PipelineConfig& config) {
    // Check if enhanced classification is requested
    if (config.enable_enhanced_features) {
        return classify_reads_enhanced(config);
    }
    
    // Otherwise use existing standard classification
    // Determine if we should use streaming
    bool should_stream = config.streaming_config.enable_streaming || config.force_streaming;
    
    // Check file size to auto-enable streaming for large files
    if (!should_stream) {
        std::ifstream file(config.reads_file, std::ios::ate | std::ios::binary);
        if (file.is_open()) {
            size_t file_size_mb = file.tellg() / (1024 * 1024);
            file.close();
            
            // Auto-enable streaming for files > 1GB
            if (file_size_mb > 1024) {
                std::cout << "Large input file detected (" << file_size_mb 
                          << " MB), enabling streaming mode" << std::endl;
                should_stream = true;
            }
        }
    }
    
    if (should_stream) {
        return classify_reads_with_streaming(config);
    }
    
    // Keep existing non-streaming code as fallback for small files
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
            
            // Check if files are gzipped
            bool is_gzipped1 = config.reads_file.substr(config.reads_file.find_last_of(".") + 1) == "gz";
            bool is_gzipped2 = config.reads_file2.substr(config.reads_file2.find_last_of(".") + 1) == "gz";
            
            if (is_gzipped1 && is_gzipped2) {
                // Use zlib for gzipped files
                gzFile gz_file1 = gzopen(config.reads_file.c_str(), "rb");
                gzFile gz_file2 = gzopen(config.reads_file2.c_str(), "rb");
                
                if (!gz_file1) {
                    std::cerr << "Cannot open reads file: " << config.reads_file << std::endl;
                    return false;
                }
                if (!gz_file2) {
                    std::cerr << "Cannot open reads file 2: " << config.reads_file2 << std::endl;
                    gzclose(gz_file1);
                    return false;
                }
                
                const int buffer_size = 1024;
                char buffer1[buffer_size];
                char buffer2[buffer_size];
                std::string line1, line2;
                int line_count1 = 0, line_count2 = 0;
                
                while (gzgets(gz_file1, buffer1, buffer_size) && gzgets(gz_file2, buffer2, buffer_size)) {
                    line1 = buffer1;
                    line2 = buffer2;
                    
                    // Remove newline
                    if (!line1.empty() && line1.back() == '\n') line1.pop_back();
                    if (!line2.empty() && line2.back() == '\n') line2.pop_back();
                    
                    line_count1++;
                    if (line_count1 % 4 == 2) {  // Sequence line
                        if (!line1.empty() && !line2.empty() && 
                            line1.length() >= config.classifier_params.k && 
                            line2.length() >= config.classifier_params.k) {
                            paired_reads.push_back(PairedRead(line1, line2, "read_" + std::to_string(paired_reads.size())));
                        }
                    }
                }
                gzclose(gz_file1);
                gzclose(gz_file2);
            } else {
                // Use standard ifstream for uncompressed files
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
                }
                fastq_file1.close();
                fastq_file2.close();
            }
            
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
            // Single-end mode - with gzip support
            std::vector<std::string> reads;
            
            // Check if file is gzipped
            bool is_gzipped = config.reads_file.substr(config.reads_file.find_last_of(".") + 1) == "gz";
            
            if (is_gzipped) {
                // Use zlib for gzipped files
                gzFile gz_file = gzopen(config.reads_file.c_str(), "rb");
                if (!gz_file) {
                    std::cerr << "Cannot open gzipped reads file: " << config.reads_file << std::endl;
                    return false;
                }
                
                char buffer[4096];
                std::string line;
                int line_count = 0;
                
                while (gzgets(gz_file, buffer, sizeof(buffer))) {
                    line = buffer;
                    // Remove newline if present
                    if (!line.empty() && line.back() == '\n') {
                        line.pop_back();
                    }
                    
                    line_count++;
                    if (line_count % 4 == 2) {  // Sequence line
                        if (!line.empty() && line.length() >= config.classifier_params.k) {
                            reads.push_back(line);
                        }
                    }
                    
                    // Remove limit - process all reads
                    // if (reads.size() >= 1000000) {
                    //     std::cout << "Limited to " << reads.size() << " reads for processing" << std::endl;
                    //     break;
                    // }
                }
                gzclose(gz_file);
            } else {
                // Use standard ifstream for uncompressed files
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
            }
            
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
        
        // Generate reports after successful classification
        bool success = true;
        if (config.generate_reports) {
            success = generate_classification_reports(config);
        }
        
        return success;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during classification: " << e.what() << std::endl;
        return false;
    }
}

// Forward declarations for enhanced classification functions
bool classify_reads_fast_enhanced(const PipelineConfig& config);
bool classify_reads_phylo_enhanced(const PipelineConfig& config);

// Enhanced classification dispatcher
bool classify_reads_enhanced(const PipelineConfig& config) {
    std::cout << "\n=== ENHANCED CLASSIFICATION ===" << std::endl;
    std::cout << "Mode: ";
    
    switch (config.classification_mode) {
        case ClassificationMode::FAST_ENHANCED:
            std::cout << "Fast Enhanced (with compact taxonomy)" << std::endl;
            return classify_reads_fast_enhanced(config);
            
        case ClassificationMode::PHYLO_ENHANCED:
            std::cout << "Phylogenetic Enhanced (with NCBI taxonomy)" << std::endl;
            return classify_reads_phylo_enhanced(config);
            
        default:
            std::cout << "Standard (falling back to original classifier)" << std::endl;
            return classify_reads_command(config);
    }
}

bool classify_reads_fast_enhanced(const PipelineConfig& config) {
    if (config.compact_taxonomy_path.empty()) {
        std::cerr << "Error: --compact-taxonomy required for fast enhanced mode" << std::endl;
        std::cerr << "Build compact taxonomy first with build_compact_taxonomy tool" << std::endl;
        return false;
    }
    
    try {
        // Configure fast enhanced parameters
        FastEnhancedParams params;
        params.k = config.classifier_params.k;
        params.ell = config.classifier_params.ell;
        params.spaces = config.classifier_params.spaces;
        params.primary_confidence_threshold = config.classifier_params.confidence_threshold;
        params.secondary_confidence_threshold = config.classifier_params.confidence_threshold * 1.5f;
        params.compact_taxonomy_path = config.compact_taxonomy_path;
        params.phylo_quality = config.phylo_quality;
        
        // Initialize fast enhanced classifier
        FastEnhancedClassifier classifier(params, 50000);
        
        // Load database
        if (!classifier.load_database(config.database_dir)) {
            std::cerr << "Failed to load database from " << config.database_dir << std::endl;
            return false;
        }
        
        // Load reads (similar to existing code but adapted for enhanced results)
        std::vector<std::string> reads;
        
        // Check if file is gzipped
        bool is_gzipped = config.reads_file.substr(config.reads_file.find_last_of(".") + 1) == "gz";
        
        if (is_gzipped) {
            // Use zlib for gzipped files
            gzFile gz_file = gzopen(config.reads_file.c_str(), "rb");
            if (!gz_file) {
                std::cerr << "Cannot open reads file: " << config.reads_file << std::endl;
                return false;
            }
            
            const int buffer_size = 1024;
            char buffer[buffer_size];
            std::string current_line;
            int line_count = 0;
            
            while (gzgets(gz_file, buffer, buffer_size)) {
                current_line += buffer;
                if (current_line.back() == '\n') {
                    current_line.pop_back();
                    line_count++;
                    if (line_count % 4 == 2) {  // Sequence line
                        if (!current_line.empty() && current_line.length() >= params.k) {
                            reads.push_back(current_line);
                        }
                    }
                    current_line.clear();
                }
            }
            gzclose(gz_file);
        } else {
            // Use standard ifstream for uncompressed files
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
                    if (!line.empty() && line.length() >= params.k) {
                        reads.push_back(line);
                    }
                }
            }
            fastq_file.close();
        }
        
        if (reads.empty()) {
            std::cerr << "No valid reads found" << std::endl;
            return false;
        }
        
        std::cout << "Loaded " << reads.size() << " reads for fast enhanced classification" << std::endl;
        
        // Classify with enhanced features
        auto results = classifier.classify_fast_enhanced(reads);
        
        // Write enhanced results
        std::ofstream output(config.output_path);
        output << "# GPU Kraken Enhanced Classification Results\n";
        
        for (size_t i = 0; i < results.size(); i++) {
            const auto& result = results[i];
            char status = (result.taxon_id > 0) ? 'C' : 'U';
            
            output << status << "\t"
                   << "read_" << i << "\t"
                   << result.taxon_id << "\t"
                   << reads[i].length() << "\t"
                   << std::fixed << std::setprecision(3) << result.primary_confidence << "\t"
                   << result.total_kmers << "\t" 
                   << result.classified_kmers << "\t"
                   << result.phylogenetic_consistency_score << "\n";
        }
        output.close();
        
        // Print enhanced statistics
        classifier.print_fast_enhanced_statistics(results);
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in fast enhanced classification: " << e.what() << std::endl;
        return false;
    }
}

bool classify_reads_phylo_enhanced(const PipelineConfig& config) {
    if (config.taxonomy_nodes_path.empty() || config.taxonomy_names_path.empty()) {
        std::cerr << "Error: --taxonomy-nodes and --taxonomy-names required for phylo enhanced mode" << std::endl;
        return false;
    }
    
    try {
        // Configure phylogenetic parameters
        Phase1EnhancedParams enhanced_params;
        enhanced_params.k = config.classifier_params.k;
        enhanced_params.ell = config.classifier_params.ell;
        enhanced_params.spaces = config.classifier_params.spaces;
        enhanced_params.primary_confidence_threshold = config.classifier_params.confidence_threshold;
        enhanced_params.secondary_confidence_threshold = config.classifier_params.confidence_threshold * 1.5f;
        enhanced_params.enable_phylogenetic_validation = true;
        
        PhylogeneticClassificationParams phylo_params;
        phylo_params.use_phylogenetic_validation = true;
        phylo_params.taxonomy_nodes_path = config.taxonomy_nodes_path;
        phylo_params.taxonomy_names_path = config.taxonomy_names_path;
        
        // Initialize phylogenetic enhanced classifier
        Phase1EnhancedClassifierWithPhylogeny classifier(enhanced_params, phylo_params, 50000);
        
        // Load database
        if (!classifier.load_database(config.database_dir)) {
            std::cerr << "Failed to load database from " << config.database_dir << std::endl;
            return false;
        }
        
        // Load reads
        std::vector<std::string> reads;
        
        // Check if file is gzipped
        bool is_gzipped = config.reads_file.substr(config.reads_file.find_last_of(".") + 1) == "gz";
        
        if (is_gzipped) {
            // Use zlib for gzipped files
            gzFile gz_file = gzopen(config.reads_file.c_str(), "rb");
            if (!gz_file) {
                std::cerr << "Cannot open reads file: " << config.reads_file << std::endl;
                return false;
            }
            
            const int buffer_size = 1024;
            char buffer[buffer_size];
            std::string current_line;
            int line_count = 0;
            
            while (gzgets(gz_file, buffer, buffer_size)) {
                current_line += buffer;
                if (current_line.back() == '\n') {
                    current_line.pop_back();
                    line_count++;
                    if (line_count % 4 == 2) {  // Sequence line
                        if (!current_line.empty() && current_line.length() >= enhanced_params.k) {
                            reads.push_back(current_line);
                        }
                    }
                    current_line.clear();
                }
            }
            gzclose(gz_file);
        } else {
            // Use standard ifstream for uncompressed files
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
                    if (!line.empty() && line.length() >= enhanced_params.k) {
                        reads.push_back(line);
                    }
                }
            }
            fastq_file.close();
        }
        
        if (reads.empty()) {
            std::cerr << "No valid reads found" << std::endl;
            return false;
        }
        
        std::cout << "Loaded " << reads.size() << " reads for phylogenetic enhanced classification" << std::endl;
        
        // Classify with phylogenetic validation
        auto results = classifier.classify_enhanced_with_phylogeny(reads);
        
        // Write results with lineage information
        std::ofstream output(config.output_path);
        output << "# GPU Kraken Phylogenetic Enhanced Classification Results\n";
        
        for (size_t i = 0; i < results.size(); i++) {
            const auto& result = results[i];
            char status = (result.taxon_id > 0) ? 'C' : 'U';
            
            output << status << "\t"
                   << "read_" << i << "\t"
                   << result.taxon_id << "\t"
                   << reads[i].length() << "\t"
                   << std::fixed << std::setprecision(3) << result.primary_confidence << "\t"
                   << result.secondary_confidence << "\t"
                   << result.phylogenetic_consistency_score << "\t"
                   << (result.passed_phylogenetic_filter ? "PHYLO_OK" : "PHYLO_FAIL") << "\t"
                   << result.classification_path << "\n";
        }
        output.close();
        
        // Print phylogenetic statistics
        classifier.print_enhanced_statistics_with_phylogeny(results);
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in phylogenetic enhanced classification: " << e.what() << std::endl;
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
    
    // Step 1: Build database (same as before)
    PipelineConfig build_config = config;
    build_config.command = "build";
    build_config.output_path = db_dir;
    
    if (!build_database_command(build_config)) {
        std::cerr << "Pipeline failed at database build step" << std::endl;
        return false;
    }
    
    // Step 2: Classify reads (enhanced or standard)
    PipelineConfig classify_config = config;
    classify_config.command = "classify";
    classify_config.database_dir = db_dir;
    classify_config.output_path = results_file;
    classify_config.report_file = report_file;
    
    bool success;
    if (config.enable_enhanced_features) {
        success = classify_reads_enhanced(classify_config);
    } else {
        success = classify_reads_command(classify_config);
    }
    
    if (!success) {
        std::cerr << "Pipeline failed at classification step" << std::endl;
        return false;
    }
    
    // Generate reports
    if (classify_config.generate_reports) {
        PipelineConfig report_config = config;
        report_config.output_path = results_file;
        report_config.database_dir = db_dir;
        generate_classification_reports(report_config);
    }
    
    std::cout << "\n✓ Complete pipeline finished successfully!" << std::endl;
    std::cout << "Results available in: " << config.output_path << std::endl;
    
    return true;
}

bool batch_classify_command(const PipelineConfig& config) {
    std::cout << "\n=== BATCH CLASSIFICATION FROM CSV ===" << std::endl;
    std::cout << "CSV file: " << config.csv_file << std::endl;
    std::cout << "Database: " << config.database_dir << std::endl;
    std::cout << "Batch output: " << config.batch_output_dir << std::endl;
    
    try {
        // Initialize CSV parser and batch processor
        BioGPU::BatchProcessor processor(config.batch_output_dir, config.create_sample_dirs);
        
        // Load samples from CSV
        if (!processor.loadSamples(config.csv_file)) {
            std::cerr << "Failed to load samples from CSV file" << std::endl;
            return false;
        }
        
        // Validate that database exists
        if (!std::filesystem::exists(config.database_dir)) {
            std::cerr << "Database directory not found: " << config.database_dir << std::endl;
            return false;
        }
        
        // Load database ONCE for all samples
        std::cout << "\nLoading database (once for all samples)..." << std::endl;
        PairedEndGPUKrakenClassifier classifier(config.classifier_params);
        if (!classifier.load_database(config.database_dir)) {
            std::cerr << "Failed to load database from " << config.database_dir << std::endl;
            return false;
        }
        std::cout << "Database loaded successfully - will be reused for all samples" << std::endl;
        
        // Print sample summary
        processor.getParser().printSummary();
        
        // Check for validation errors if enabled
        if (config.validate_paths_on_load) {
            auto errors = processor.getParser().getValidationErrors();
            if (!errors.empty()) {
                std::cerr << "\nValidation errors found:" << std::endl;
                for (const auto& error : errors) {
                    std::cerr << "  - " << error << std::endl;
                }
                if (config.stop_on_error) {
                    return false;
                }
            }
        }
        
        // Define the processing function for each sample
        auto process_sample = [&](const BioGPU::SampleInfo& sample, const std::string& output_path) -> int {
            try {
                // Create a temporary config for this sample
                PipelineConfig sample_config = config;
                sample_config.reads_file = sample.read1_path;
                sample_config.reads_file2 = sample.read2_path;
                sample_config.output_path = output_path + "_classification.txt";
                
                // Add report if not disabled
                if (!config.report_file.empty() || config.create_sample_dirs) {
                    sample_config.report_file = output_path + "_report.txt";
                }
                
                std::cout << "  Processing: " << sample.sample_name << std::endl;
                std::cout << "    R1: " << sample.read1_path << std::endl;
                if (!sample.read2_path.empty()) {
                    std::cout << "    R2: " << sample.read2_path << std::endl;
                }
                std::cout << "    Output: " << sample_config.output_path << std::endl;
                
                // Use streaming classification with pre-loaded classifier
                bool success = classify_reads_with_streaming_and_classifier(sample_config, classifier);
                
                // Generate reports for this sample
                if (success && config.generate_reports) {
                    PipelineConfig report_config = config;
                    report_config.output_path = sample_config.output_path;
                    report_config.reports_sample_name = sample.sample_name;
                    generate_classification_reports(report_config);
                }
                
                return success ? 0 : 1;
                
            } catch (const std::exception& e) {
                std::cerr << "    Error processing sample " << sample.sample_name 
                          << ": " << e.what() << std::endl;
                return 1;
            }
        };
        
        // Run batch processing
        int failed_count = processor.processBatch(process_sample, config.stop_on_error);
        
        if (failed_count == 0) {
            std::cout << "\n✓ All samples processed successfully!" << std::endl;
        } else {
            std::cout << "\n⚠ Batch completed with " << failed_count << " failures" << std::endl;
        }
        
        // Generate batch summary report
        std::string summary_file = config.batch_output_dir + "/batch_summary.txt";
        std::ofstream summary(summary_file);
        if (summary.is_open()) {
            summary << "Batch Classification Summary\n";
            summary << "============================\n";
            summary << "CSV file: " << config.csv_file << "\n";
            summary << "Database: " << config.database_dir << "\n";
            summary << "Total samples: " << processor.getParser().getSampleCount() << "\n";
            summary << "Failed samples: " << failed_count << "\n";
            summary << "Success rate: " << std::fixed << std::setprecision(1) 
                    << (100.0 * (processor.getParser().getSampleCount() - failed_count) / processor.getParser().getSampleCount()) << "%\n";
            
            summary << "\nProcessing Parameters:\n";
            summary << "k-mer length: " << config.classifier_params.k << "\n";
            summary << "Confidence threshold: " << config.classifier_params.confidence_threshold << "\n";
            summary << "Streaming batch size: " << config.streaming_config.read_batch_size << "\n";
            
            summary.close();
            std::cout << "Batch summary written to: " << summary_file << std::endl;
        }
        
        return failed_count == 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during batch processing: " << e.what() << std::endl;
        return false;
    }
}

// Quick memory check function for user guidance
void print_memory_recommendations() {
    int device_count;
    cudaError_t cuda_status = cudaGetDeviceCount(&device_count);
    if (cuda_status != cudaSuccess || device_count == 0) {
        std::cerr << "No CUDA-capable GPU found for memory analysis" << std::endl;
        return;
    }
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    
    size_t total_memory = prop.totalGlobalMem;
    size_t memory_gb = total_memory / (1024*1024*1024);
    
    std::cout << "\n=== MEMORY RECOMMENDATIONS ===" << std::endl;
    std::cout << "GPU: " << prop.name << " (" << memory_gb << " GB)" << std::endl;
    
    // Provide capacity recommendations based on GPU memory
    if (memory_gb >= 40) {
        std::cout << "Recommended minimizer capacity: 15,000,000 - 25,000,000" << std::endl;
        std::cout << "Your GPU has excellent memory - use high capacity for comprehensive databases" << std::endl;
    } else if (memory_gb >= 24) {
        std::cout << "Recommended minimizer capacity: 8,000,000 - 15,000,000" << std::endl;
        std::cout << "Your GPU has good memory - suitable for large microbial databases" << std::endl;
    } else if (memory_gb >= 16) {
        std::cout << "Recommended minimizer capacity: 5,000,000 - 10,000,000" << std::endl;
        std::cout << "Your GPU has moderate memory - good for standard databases" << std::endl;
    } else if (memory_gb >= 8) {
        std::cout << "Recommended minimizer capacity: 2,000,000 - 5,000,000" << std::endl;
        std::cout << "Your GPU has limited memory - consider smaller batches" << std::endl;
    } else {
        std::cout << "Recommended minimizer capacity: 1,000,000 - 2,000,000" << std::endl;
        std::cout << "Your GPU has very limited memory - use conservative settings" << std::endl;
    }
    
    std::cout << "Note: Enable --auto-memory for automatic optimization" << std::endl;
}

// NEW: Generate comprehensive reports after classification
bool generate_classification_reports(const PipelineConfig& config) {
    if (!config.generate_reports) {
        std::cout << "Report generation disabled - skipping reports" << std::endl;
        return true;
    }
    
    std::cout << "\n=== GENERATING CLASSIFICATION REPORTS ===" << std::endl;
    
    try {
        // Determine sample name
        std::string sample_name = config.reports_sample_name;
        if (sample_name.empty()) {
            // Auto-detect from reads file name
            std::filesystem::path reads_path(config.reads_file);
            sample_name = reads_path.stem().string();
            
            // Remove common suffixes
            if (sample_name.size() >= 3 && 
                (sample_name.substr(sample_name.size() - 3) == "_R1" || 
                 sample_name.substr(sample_name.size() - 3) == "_r1")) {
                sample_name = sample_name.substr(0, sample_name.length() - 3);
            }
            if ((sample_name.size() >= 6 && sample_name.substr(sample_name.size() - 6) == ".fastq") ||
                (sample_name.size() >= 3 && sample_name.substr(sample_name.size() - 3) == ".fq")) {
                sample_name = sample_name.substr(0, sample_name.find_last_of('.'));
            }
        }
        
        std::cout << "Sample name: " << sample_name << std::endl;
        
        // Create reports directory
        std::string reports_dir;
        std::filesystem::path output_path(config.output_path);
        
        if (output_path.extension() == ".txt") {
            // Output is a file, create reports directory next to it
            reports_dir = output_path.parent_path() / (sample_name + "_reports");
        } else {
            // Output is a directory, create reports subdirectory
            reports_dir = output_path / "reports";
        }
        
        std::filesystem::create_directories(reports_dir);
        std::cout << "Reports directory: " << reports_dir << std::endl;
        
        // Initialize report generator
        ClassificationReportGenerator generator(reports_dir);
        generator.set_abundance_threshold(config.reports_min_abundance);
        generator.include_unclassified_reads(config.reports_include_unclassified);
        
        // Load classification results
        std::string kraken_output = config.output_path;
        if (std::filesystem::is_directory(kraken_output)) {
            kraken_output = kraken_output + "/classification_results.txt";
        }
        
        if (!generator.load_kraken_results(kraken_output)) {
            std::cerr << "Failed to load classification results from: " << kraken_output << std::endl;
            return false;
        }
        
        // Load taxonomy
        std::string taxonomy_file;
        if (config.command == "classify" || config.command == "batch") {
            taxonomy_file = config.database_dir + "/taxonomy.tsv";
        } else if (config.command == "pipeline") {
            taxonomy_file = config.output_path + "/database/taxonomy.tsv";
        }
        
        if (!std::filesystem::exists(taxonomy_file)) {
            std::cerr << "Warning: Taxonomy file not found at: " << taxonomy_file << std::endl;
            std::cerr << "Reports will have limited taxonomic information" << std::endl;
        } else {
            if (!generator.load_taxonomy(taxonomy_file)) {
                std::cerr << "Warning: Failed to load taxonomy from: " << taxonomy_file << std::endl;
            }
        }
        
        // Generate reports based on configuration
        bool success = false;
        
        if (config.reports_only_mpa) {
            std::string mpa_file = reports_dir + "/" + sample_name + "_profile.mpa";
            generator.compute_statistics();
            success = generator.generate_mpa_report(mpa_file);
            generator.print_summary();
            
        } else if (config.reports_only_species) {
            std::string species_file = reports_dir + "/" + sample_name + "_species.tsv";
            generator.compute_statistics();
            success = generator.generate_species_counts(species_file);
            generator.print_summary();
            
        } else if (config.reports_only_genus) {
            std::string genus_file = reports_dir + "/" + sample_name + "_genus.tsv";
            generator.compute_statistics();
            success = generator.generate_genus_counts(genus_file);
            generator.print_summary();
            
        } else {
            // Generate all reports
            success = generator.generate_all_reports(sample_name);
        }
        
        if (success) {
            std::cout << "\n✓ Classification reports generated successfully!" << std::endl;
            std::cout << "Reports available in: " << reports_dir << std::endl;
        }
        
        return success;
        
    } catch (const std::exception& e) {
        std::cerr << "Error generating reports: " << e.what() << std::endl;
        return false;
    }
}

// Safe parse_arguments function with heap validation
bool parse_arguments_debug(int argc, char* argv[], PipelineConfig& config) {
    HeapValidator::validate_heap("parse_arguments start");
    
    if (argc < 2) {
        return false;
    }
    
    std::cout << "DEBUG: Setting command from argv[1]: " << argv[1] << std::endl;
    
    try {
        config.command = std::string(argv[1]);  // Explicit string construction
        HeapValidator::validate_heap("after command assignment");
        HeapValidator::check_string_safety(config.command, "config.command");
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Exception setting command: " << e.what() << std::endl;
        return false;
    }
    
    if (config.command != "build" && config.command != "classify" && 
        config.command != "pipeline" && config.command != "batch") {
        std::cerr << "Error: Unknown command '" << config.command << "'" << std::endl;
        return false;
    }
    
    std::cout << "DEBUG: Processing " << (argc - 2) << " additional arguments..." << std::endl;
    
    for (int i = 2; i < argc; i++) {
        std::cout << "DEBUG: Processing arg " << i << ": " << argv[i] << std::endl;
        HeapValidator::validate_heap("before argument processing");
        
        try {
            std::string arg(argv[i]);  // Explicit string construction
            HeapValidator::validate_heap("after arg string creation");
            
            if (arg == "--genome-dir" && i + 1 < argc) {
                std::cout << "DEBUG: Setting genome_dir from argv[" << (i+1) << "]: " << argv[i+1] << std::endl;
                config.genome_dir = std::string(argv[++i]);
                HeapValidator::validate_heap("after genome_dir assignment");
                HeapValidator::check_string_safety(config.genome_dir, "config.genome_dir");
            } else if (arg == "--database" && i + 1 < argc) {
                config.database_dir = std::string(argv[++i]);
                HeapValidator::validate_heap("after database_dir assignment");
                HeapValidator::check_string_safety(config.database_dir, "config.database_dir");
            } else if (arg == "--reads" && i + 1 < argc) {
                config.reads_file = std::string(argv[++i]);
                HeapValidator::validate_heap("after reads_file assignment");
                HeapValidator::check_string_safety(config.reads_file, "config.reads_file");
            } else if (arg == "--output" && i + 1 < argc) {
                config.output_path = std::string(argv[++i]);
                HeapValidator::validate_heap("after output_path assignment");
                HeapValidator::check_string_safety(config.output_path, "config.output_path");
            } else if (arg == "--taxonomy" && i + 1 < argc) {
                config.taxonomy_dir = std::string(argv[++i]);
                HeapValidator::validate_heap("after taxonomy_dir assignment");
                HeapValidator::check_string_safety(config.taxonomy_dir, "config.taxonomy_dir");
            }
            // Add other arguments as needed...
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR: Exception processing argument " << i << ": " << e.what() << std::endl;
            HeapValidator::validate_heap("after argument processing error");
            return false;
        }
        
        HeapValidator::validate_heap("after argument processing");
    }
    
    std::cout << "DEBUG: Argument parsing completed successfully" << std::endl;
    HeapValidator::validate_heap("parse_arguments end");
    
    return true;
}

// Replace the original main with a minimal version for debugging
int main(int argc, char* argv[]) {
    // Add heap corruption debugging at the very start
    
    std::cout << "GPU-Accelerated Kraken2-Style Taxonomic Classifier (High Capacity)" << std::endl;
    std::cout << "=================================================================" << std::endl;
    
    
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
    
    // Print memory recommendations
    print_memory_recommendations();
    
    bool success = false;
    
    try {
        if (config.command == "build") {
            success = build_database_command(config);
        } else if (config.command == "classify") {
            success = classify_reads_command(config);
        } else if (config.command == "pipeline") {
            success = pipeline_command(config);
        } else if (config.command == "batch") {
            success = batch_classify_command(config);
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