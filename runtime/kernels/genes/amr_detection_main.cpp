// amr_detection_main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <zlib.h>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include "amr_detection_pipeline.h"
#include "sample_csv_parser.h"  // Reuse from your FQ pipeline
#include "hdf5_amr_writer.h"
#include "clinical_amr_report_generator.h"

namespace fs = std::filesystem;

// Data structures from resistance pipeline for batch reading
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};

struct ReadBatch {
    char* sequences;
    int* lengths;
    int* offsets;
    int num_reads;
    int total_bases;
};

// Function to check if file is gzipped
bool isGzipped(const std::string& filename) {
    return filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz";
}


// Batch reading function from resistance pipeline that reads two files efficiently
bool readBatch(gzFile gz_r1, gzFile gz_r2, std::vector<FastqRecord>& batch_r1,
               std::vector<FastqRecord>& batch_r2, int max_size) {
    char buffer[1024];
    
    for (int i = 0; i < max_size; i++) {
        FastqRecord rec1, rec2;
        
        // Read R1
        if (gzgets(gz_r1, buffer, 1024) == NULL) return i > 0;
        rec1.header = std::string(buffer);
        if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
        rec1.sequence = std::string(buffer);
        rec1.sequence.pop_back(); // Remove newline
        if (gzgets(gz_r1, buffer, 1024) == NULL) return false; // +
        if (gzgets(gz_r1, buffer, 1024) == NULL) return false;
        rec1.quality = std::string(buffer);
        rec1.quality.pop_back();
        
        // Read R2
        if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
        rec2.header = std::string(buffer);
        if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
        rec2.sequence = std::string(buffer);
        rec2.sequence.pop_back();
        if (gzgets(gz_r2, buffer, 1024) == NULL) return false; // +
        if (gzgets(gz_r2, buffer, 1024) == NULL) return false;
        rec2.quality = std::string(buffer);
        rec2.quality.pop_back();
        
        batch_r1.push_back(rec1);
        batch_r2.push_back(rec2);
    }
    
    return true;
}

// Single-file batch reading function for single-end reads
bool readBatchSingle(gzFile gz_file, std::vector<FastqRecord>& batch, int max_size) {
    char buffer[1024];
    
    for (int i = 0; i < max_size; i++) {
        FastqRecord rec;
        
        // Read header
        if (gzgets(gz_file, buffer, 1024) == NULL) return i > 0;
        rec.header = std::string(buffer);
        // Read sequence
        if (gzgets(gz_file, buffer, 1024) == NULL) return false;
        rec.sequence = std::string(buffer);
        rec.sequence.pop_back(); // Remove newline
        // Read + line
        if (gzgets(gz_file, buffer, 1024) == NULL) return false;
        // Read quality
        if (gzgets(gz_file, buffer, 1024) == NULL) return false;
        rec.quality = std::string(buffer);
        rec.quality.pop_back();
        
        batch.push_back(rec);
    }
    
    return true;
}

// Prepare batch for GPU processing
ReadBatch prepareBatch(const std::vector<FastqRecord>& records) {
    ReadBatch batch;
    batch.num_reads = records.size();
    batch.total_bases = 0;
    
    for (const auto& rec : records) {
        batch.total_bases += rec.sequence.length();
    }
    
    batch.sequences = new char[batch.total_bases];
    batch.lengths = new int[batch.num_reads];
    batch.offsets = new int[batch.num_reads];
    
    int offset = 0;
    for (int i = 0; i < batch.num_reads; i++) {
        const std::string& seq = records[i].sequence;
        memcpy(batch.sequences + offset, seq.c_str(), seq.length());
        batch.lengths[i] = seq.length();
        batch.offsets[i] = offset;
        offset += seq.length();
    }
    
    return batch;
}


// Function to process paired-end sample
void processSamplePaired(AMRDetectionPipeline& pipeline, 
                        const std::string& sample_name,
                        const std::string& fastq_r1,
                        const std::string& fastq_r2,
                        const std::string& output_dir,
                        const AMRDetectionConfig& config,  // Add config parameter
                        bool merge_reads = true) {
    
    std::cout << "\n=== Processing paired-end sample: " << sample_name << " ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Create output directory if it doesn't exist
    fs::create_directories(output_dir);
    
    // Create HDF5 writer
    std::string hdf5_path = output_dir + "/" + sample_name + "_amr_results.h5";
    HDF5AMRWriter hdf5_writer(hdf5_path);
    
    // Get the database paths from config
    std::string dna_db_path = config.protein_db_path;  // Assuming this contains the database info
    std::string protein_db_path = config.protein_db_path;
    hdf5_writer.initialize(sample_name, dna_db_path, protein_db_path);
    
    // Disable auto-flush during batch processing for better performance
    hdf5_writer.setAutoFlush(false);
    
    // Create clinical report generator
    std::string report_prefix = output_dir + "/" + sample_name;
    ClinicalAMRReportGenerator report_generator(report_prefix, sample_name);
    
    // Open the gzipped FASTQ files
    gzFile gz_r1 = gzopen(fastq_r1.c_str(), "r");
    gzFile gz_r2 = gzopen(fastq_r2.c_str(), "r");
    if (!gz_r1 || !gz_r2) {
        std::cerr << "Error: Cannot open one or both FASTQ files." << std::endl;
        return;
    }
    
    std::cout << "Processing paired-end reads with enhanced paired-end handling..." << std::endl;
    
    // Collect all AMR hits across batches
    std::vector<AMRHit> all_amr_hits;
    
    // Initialize accumulated hits if EM is enabled
    if (config.use_em) {
        std::vector<AMRHit> empty_hits;
        pipeline.setAccumulatedHits(empty_hits);
        std::cout << "EM enabled - will accumulate hits across all batches" << std::endl;
    }
    
    const int batch_size = config.reads_per_batch;
    int batch_num = 0;
    int total_pairs_processed = 0;
    
    // Loop through the file by reading one batch at a time
    while (true) {
        std::vector<FastqRecord> batch_r1_records, batch_r2_records;
        if (!readBatch(gz_r1, gz_r2, batch_r1_records, batch_r2_records, batch_size)) {
            break; // End of file
        }
        
        if (batch_r1_records.empty()) {
            break;
        }
        
        std::cout << "\nProcessing batch " << ++batch_num << std::endl;
        
        // Build read pair data for new approach
        std::vector<ReadPairData> batch_pairs;
        batch_pairs.reserve(batch_r1_records.size());
        
        for (size_t i = 0; i < batch_r1_records.size(); ++i) {
            // Extract read ID by removing the '@' and potential newline
            std::string header1 = batch_r1_records[i].header;
            std::string header2 = batch_r2_records[i].header;
            header1.erase(std::remove(header1.begin(), header1.end(), '\n'), header1.end());
            header2.erase(std::remove(header2.begin(), header2.end(), '\n'), header2.end());
            
            // Skip pairs where either read is empty
            if (!batch_r1_records[i].sequence.empty() && !batch_r2_records[i].sequence.empty()) {
                ReadPairData pair_data;
                pair_data.read1_seq = batch_r1_records[i].sequence;
                pair_data.read2_seq = batch_r2_records[i].sequence;
                pair_data.read1_id = header1.substr(1) + "_R1";  // Remove '@'
                pair_data.read2_id = header2.substr(1) + "_R2";  // Remove '@'
                pair_data.pair_index = i;
                batch_pairs.push_back(pair_data);
            }
        }
        
        // Process batch of paired reads with new approach
        pipeline.processPairedBatch(batch_pairs);
        
        // Update total pairs processed
        total_pairs_processed += batch_pairs.size();
        
        // Get results from this batch
        auto batch_hits = pipeline.getAMRHits();
        
        std::cout << "Batch " << batch_num << " hits: " << batch_hits.size() << std::endl;
        
        // Accumulate hits locally for reporting
        all_amr_hits.insert(all_amr_hits.end(), batch_hits.begin(), batch_hits.end());
        
        std::cout << "Total accumulated hits so far: " << all_amr_hits.size() << std::endl;
        
        // Write batch hits to HDF5
        hdf5_writer.addAMRHits(batch_hits);
    }
    
    // Close the files
    gzclose(gz_r1);
    gzclose(gz_r2);
    
    // Finalize coverage statistics after all batches but BEFORE EM
    std::cout << "\nFinalizing coverage statistics from GPU..." << std::endl;
    pipeline.finalizeCoverageStats();
    
    // NOW run EM algorithm AFTER finalizeCoverageStats
    std::cout << "\n=== CHECKING IF EM SHOULD RUN ===" << std::endl;
    std::cout << "config.use_em = " << (config.use_em ? "true" : "false") << std::endl;
    std::cout << "all_amr_hits.size() = " << all_amr_hits.size() << std::endl;
    
    if (config.use_em) {
        std::cout << "\n=== CALLING EM ALGORITHM ===" << std::endl;
        std::cout << "\n=== Running EM Algorithm on All Accumulated Hits ===" << std::endl;
        std::cout << "Total hits before EM: " << all_amr_hits.size() << std::endl;
        
        // Set accumulated hits in pipeline explicitly
        pipeline.setAccumulatedHits(all_amr_hits);
        std::cout << "Set " << all_amr_hits.size() << " hits for EM processing" << std::endl;
        
        // Run the EM algorithm
        pipeline.resolveAmbiguousAssignmentsEM();
        
        // Get updated coverage stats after EM (the EM updates these internally)
        std::cout << "EM algorithm completed - coverage stats updated on host" << std::endl;
    }
    
    // Get final coverage statistics and gene entries
    auto coverage_stats = pipeline.getCoverageStats();
    auto gene_entries = pipeline.getGeneEntries();
    
    std::cout << "Total AMR hits accumulated: " << all_amr_hits.size() << std::endl;
    std::cout << "Gene entries available: " << gene_entries.size() << std::endl;
    
    // Write coverage stats to HDF5
    hdf5_writer.addCoverageStats(coverage_stats, gene_entries);
    
    // Manually flush all buffered data before finalizing
    std::cout << "Flushing all buffered HDF5 data..." << std::endl;
    hdf5_writer.manualFlush();
    
    // Finalize HDF5
    std::string json_summary = output_dir + "/" + sample_name + "_hdf5_summary.json";
    hdf5_writer.finalize(json_summary);
    
    // Generate clinical report with all accumulated hits
    std::cout << "\nGenerating clinical report..." << std::endl;
    std::cout << "Passing " << all_amr_hits.size() << " hits to report generator" << std::endl;
    report_generator.processAMRResults(all_amr_hits, coverage_stats, 
                                      gene_entries, total_pairs_processed);
    report_generator.generateReports();
    
    // Write basic TSV results (backward compatibility)
    std::string output_prefix = output_dir + "/" + sample_name;
    std::cout << "\nWriting TSV results..." << std::endl;
    pipeline.writeResults(output_prefix);
    pipeline.exportAbundanceTable(output_prefix + "_abundance.tsv");
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\nSample processing completed in " << duration.count() << " seconds" << std::endl;
    std::cout << "Results written to:" << std::endl;
    std::cout << "  - HDF5: " << hdf5_path << std::endl;
    std::cout << "  - Clinical report: " << report_prefix << "_clinical_amr_report.html" << std::endl;
    std::cout << "  - Abundance table: " << output_prefix << "_abundance.tsv" << std::endl;
    std::cout << "  - Coverage stats: " << output_prefix << "_coverage.tsv" << std::endl;
    std::cout << "  - Raw hits: " << output_prefix << "_hits.tsv" << std::endl;
}

// Function to process a single-end sample
void processSample(AMRDetectionPipeline& pipeline, 
                  const std::string& sample_name,
                  const std::string& fastq_file,
                  const std::string& output_dir,
                  const AMRDetectionConfig& config) {
    
    std::cout << "\n=== Processing single-end sample: " << sample_name << " ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Create output directory if it doesn't exist
    fs::create_directories(output_dir);
    
    // Create HDF5 writer
    std::string hdf5_path = output_dir + "/" + sample_name + "_amr_results.h5";
    HDF5AMRWriter hdf5_writer(hdf5_path);
    
    std::string dna_db_path = config.protein_db_path;
    std::string protein_db_path = config.protein_db_path;
    hdf5_writer.initialize(sample_name, dna_db_path, protein_db_path);
    
    // Disable auto-flush during batch processing for better performance
    hdf5_writer.setAutoFlush(false);
    
    // Create clinical report generator
    std::string report_prefix = output_dir + "/" + sample_name;
    ClinicalAMRReportGenerator report_generator(report_prefix, sample_name);
    
    // Open FASTQ file
    std::cout << "Opening FASTQ file: " << fastq_file << std::endl;
    gzFile gz_file = gzopen(fastq_file.c_str(), "r");
    if (!gz_file) {
        std::cerr << "Error: Cannot open FASTQ file: " << fastq_file << std::endl;
        return;
    }
    
    // Collect all AMR hits
    std::vector<AMRHit> all_amr_hits;
    
    const int batch_size = config.reads_per_batch;
    int batch_num = 0;
    int total_reads_processed = 0;
    
    // Loop through file by reading one batch at a time
    while (true) {
        std::vector<FastqRecord> batch_records;
        if (!readBatchSingle(gz_file, batch_records, batch_size)) {
            break; // End of file
        }
        
        if (batch_records.empty()) {
            break;
        }
        
        std::cout << "\nProcessing batch " << ++batch_num << std::endl;
        
        // Extract reads and IDs for this batch
        std::vector<std::string> batch_reads;
        std::vector<std::string> batch_ids;
        
        for (size_t i = 0; i < batch_records.size(); i++) {
            // Extract read ID by removing the '@' and potential newline
            std::string header = batch_records[i].header;
            header.erase(std::remove(header.begin(), header.end(), '\n'), header.end());
            
            // Skip empty reads
            if (!batch_records[i].sequence.empty()) {
                batch_ids.push_back(header.substr(1)); // Remove '@'
                batch_reads.push_back(batch_records[i].sequence);
            }
        }
        
        // Update total reads processed
        total_reads_processed += batch_reads.size();
        
        // Process batch
        pipeline.processBatch(batch_reads, batch_ids);
        
        // Get results from this batch
        auto batch_hits = pipeline.getAMRHits();
        
        std::cout << "Batch " << batch_num << " hits: " << batch_hits.size() << std::endl;
        
        // Accumulate hits locally for reporting
        all_amr_hits.insert(all_amr_hits.end(), batch_hits.begin(), batch_hits.end());
        
        std::cout << "Total accumulated hits so far: " << all_amr_hits.size() << std::endl;
        
        // Write batch hits to HDF5
        hdf5_writer.addAMRHits(batch_hits);
    }
    
    // Close the file
    gzclose(gz_file);
    
    // Finalize coverage statistics after all batches but BEFORE EM
    std::cout << "\nFinalizing coverage statistics from GPU..." << std::endl;
    pipeline.finalizeCoverageStats();
    
    // Check if EM should run for single-end (if supported)
    if (config.use_em && all_amr_hits.size() > 0) {
        std::cout << "\n=== Running EM Algorithm on Single-End Data ===" << std::endl;
        std::cout << "Total hits before EM: " << all_amr_hits.size() << std::endl;
        
        // Set accumulated hits in pipeline
        pipeline.setAccumulatedHits(all_amr_hits);
        
        // Run the EM algorithm
        pipeline.resolveAmbiguousAssignmentsEM();
        
        std::cout << "EM algorithm completed - coverage stats updated on host" << std::endl;
    }
    
    // Get final coverage statistics and gene entries
    auto coverage_stats = pipeline.getCoverageStats();
    auto gene_entries = pipeline.getGeneEntries();
    
    // Write coverage stats to HDF5
    hdf5_writer.addCoverageStats(coverage_stats, gene_entries);
    
    // Manually flush all buffered data before finalizing
    std::cout << "Flushing all buffered HDF5 data..." << std::endl;
    hdf5_writer.manualFlush();
    
    // Finalize HDF5
    std::string json_summary = output_dir + "/" + sample_name + "_hdf5_summary.json";
    hdf5_writer.finalize(json_summary);
    
    // Generate clinical report
    report_generator.processAMRResults(all_amr_hits, coverage_stats, 
                                      gene_entries, total_reads_processed);
    report_generator.generateReports();
    
    // Write basic TSV results (backward compatibility)
    std::string output_prefix = output_dir + "/" + sample_name;
    std::cout << "\nWriting TSV results..." << std::endl;
    pipeline.writeResults(output_prefix);
    pipeline.exportAbundanceTable(output_prefix + "_abundance.tsv");
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\nSample processing completed in " << duration.count() << " seconds" << std::endl;
    std::cout << "Results written to:" << std::endl;
    std::cout << "  - HDF5: " << hdf5_path << std::endl;
    std::cout << "  - Clinical report: " << report_prefix << "_clinical_amr_report.html" << std::endl;
    std::cout << "  - Abundance table: " << output_prefix << "_abundance.tsv" << std::endl;
    std::cout << "  - Coverage stats: " << output_prefix << "_coverage.tsv" << std::endl;
}

// Main function
int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <amr_db_path> <input_csv> <output_dir> [options]" << std::endl;
        std::cerr << "  amr_db_path: Path to AMR database (format: dna.fasta,protein.fasta)" << std::endl;
        std::cerr << "  input_csv: CSV file with columns: sample_name,fastq_path" << std::endl;
        std::cerr << "  output_dir: Directory for output files" << std::endl;
        std::cerr << "\nOptions:" << std::endl;
        std::cerr << "  --use-bloom-filter    Use bloom filter for pre-screening (default: disabled)" << std::endl;
        std::cerr << "  --no-merge           Don't merge paired-end reads" << std::endl;
        std::cerr << "  --min-identity <f>   Minimum identity threshold (default: 0.90)" << std::endl;
        std::cerr << "  --min-coverage <f>   Minimum coverage threshold (default: 0.80)" << std::endl;
        std::cerr << "  --em                 Enable EM algorithm for multi-mapping reads" << std::endl;
        std::cerr << "  --em-iterations <n>  Number of EM iterations (default: 10)" << std::endl;
        std::cerr << "  --min-hit-coverage <f> Minimum hit coverage for EM (default: none)" << std::endl;
        return 1;
    }
    
    std::string amr_db_path = argv[1];
    std::string input_csv = argv[2];
    std::string output_dir = argv[3];
    
    // Configuration
    AMRDetectionConfig config;
    config.bloom_filter_size = 1ULL << 30;  // 1GB bloom filter
    config.kmer_length = 31;
    config.minimizer_k = 15;
    config.minimizer_w = 10;
    config.min_identity = 0.90f;
    config.min_coverage = 0.80f;
    config.min_alignment_length = 50;
    config.band_width = 15;
    
    // Batch size for paired-end reads (will be doubled when both R1 and R2 are processed)
    config.reads_per_batch = 100000;  // This means 200000 total reads when paired
    
    config.max_read_length = 300;
    config.output_prefix = output_dir + "/amr_results";
    
    // Additional configuration for paired-end processing
    bool merge_paired_reads = true;
    
    // Debug: Show all arguments
    std::cout << "Total arguments: " << argc << std::endl;
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "] = '" << argv[i] << "'" << std::endl;
    }
    
    // Optional: Parse additional config from command line
    int i = 4;
    while (i < argc) {
        std::string arg = argv[i];
        std::cout << "Processing argument " << i << ": '" << arg << "'" << std::endl;
        
        // Handle flags without values first
        if (arg == "--no-merge") {
            merge_paired_reads = false;
            i++; // Move to next argument
            continue;
        } else if (arg == "--use-bloom-filter") {
            config.use_bloom_filter = true;
            i++; // Move to next argument
            continue;
        } else if (arg == "--em") {
            config.use_em = true;
            std::cout << "EM algorithm ENABLED via command line" << std::endl;
            i++; // Move to next argument
            continue;
        }
        
        // Handle arguments with values
        if (i + 1 < argc) {
            if (arg == "--min-identity") {
                config.min_identity = std::stof(argv[i + 1]);
                i += 2; // Skip the value
            } else if (arg == "--min-coverage") {
                config.min_coverage = std::stof(argv[i + 1]);
                i += 2; // Skip the value
            } else if (arg == "--kmer-length") {
                config.kmer_length = std::stoi(argv[i + 1]);
                i += 2; // Skip the value
            } else if (arg == "--batch-size") {
                config.reads_per_batch = std::stoi(argv[i + 1]);
                i += 2; // Skip the value
            } else if (arg == "--protein-db") {
                config.protein_db_path = argv[i + 1];
                i += 2; // Skip the value
            } else if (arg == "--em-iterations") {
                config.em_iterations = std::stoi(argv[i + 1]);
                i += 2; // Skip the value
            } else if (arg == "--min-hit-coverage") {
                config.min_hit_coverage = std::stof(argv[i + 1]);
                i += 2; // Skip the value
            } else {
                // Unknown argument with no handler
                std::cerr << "Warning: Unknown argument '" << arg << "'" << std::endl;
                i++; // Skip it
            }
        } else {
            // No value available for this argument
            std::cerr << "Warning: Unknown argument '" << arg << "'" << std::endl;
            i++; // Skip it
        }
    }
    
    std::cout << "Final configuration:" << std::endl;
    std::cout << "  use_em: " << (config.use_em ? "true" : "false") << std::endl;
    
    std::cout << "=== Clinical Diagnostic Pipeline: Antibiotic Resistance Gene Detection ===" << std::endl;
    std::cout << "Purpose: Detect resistance genes in patient microbiome samples for treatment guidance" << std::endl;
    std::cout << "\nReference Database: " << amr_db_path << std::endl;
    std::cout << "Patient Samples CSV: " << input_csv << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "\nDetection Parameters:" << std::endl;
    std::cout << "  Min identity: " << config.min_identity << std::endl;
    std::cout << "  Min coverage: " << config.min_coverage << std::endl;
    std::cout << "  DNA k-mer length: " << config.kmer_length << std::endl;
    std::cout << "  Protein k-mer size: " << config.protein_kmer_size << std::endl;
    std::cout << "  Batch size: " << config.reads_per_batch << std::endl;
    if (config.use_em) {
        std::cout << "  EM algorithm: ENABLED (" << config.em_iterations << " iterations)" << std::endl;
    } else {
        std::cout << "  EM algorithm: DISABLED" << std::endl;
    }
    
    // Initialize pipeline
    AMRDetectionPipeline pipeline(config);
    
    std::cout << "\nLoading reference database of known resistance genes..." << std::endl;
    if (!pipeline.initialize(amr_db_path)) {
        std::cerr << "Failed to initialize diagnostic pipeline" << std::endl;
        return 1;
    }
    
    // Parse input CSV using the proper parser
    std::cout << "\nParsing input CSV..." << std::endl;
    BioGPU::SampleCSVParser csv_parser;
    
    if (!csv_parser.parseFile(input_csv, true)) {
        std::cerr << "Failed to parse input CSV: " << input_csv << std::endl;
        return 1;
    }
    
    const auto& samples = csv_parser.getSamples();
    std::cout << "Found " << samples.size() << " patient samples to analyze" << std::endl;
    
    // Print sample information
    for (const auto& sample : samples) {
        std::cout << "  " << sample.sample_name << ": ";
        if (sample.isPairedEnd()) {
            std::cout << "paired-end (R1 + R2)" << std::endl;
        } else {
            std::cout << "single-end" << std::endl;
        }
    }
    
    // Process each patient sample
    auto total_start = std::chrono::high_resolution_clock::now();
    
    for (size_t i = 0; i < samples.size(); i++) {
        const auto& sample = samples[i];
        std::cout << "\n[" << (i + 1) << "/" << samples.size() << "] ";
        
        // Clear any accumulated results from previous samples
        pipeline.clearResults();
        
        // Clear accumulated hits for new sample when EM is enabled
        if (config.use_em) {
            std::vector<AMRHit> empty_hits;
            pipeline.setAccumulatedHits(empty_hits);
        }
        
        if (sample.isPairedEnd()) {
            // Process paired-end reads with enhanced handling (no merging)
            processSamplePaired(pipeline, sample.sample_name, 
                               sample.read1_path, sample.read2_path, output_dir, config, false);
        } else {
            // Process single-end reads
            processSample(pipeline, sample.sample_name, sample.read1_path, output_dir, config);
        }
    }
    
    // Generate diagnostic summary across all patient samples
    std::cout << "\n=== Generating clinical diagnostic summary ===" << std::endl;
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::minutes>(total_end - total_start);
    
    std::cout << "\n=== Clinical Diagnostic Analysis Complete ===" << std::endl;
    std::cout << "Total analysis time: " << total_duration.count() << " minutes" << std::endl;
    std::cout << "Clinical reports and results written to: " << output_dir << std::endl;
    std::cout << "These results can guide antibiotic selection for patient treatment" << std::endl;
    
    
    return 0;
}
