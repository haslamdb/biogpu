// amr_detection_main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <zlib.h>
#include <filesystem>
#include "amr_detection_pipeline.h"
#include "sample_csv_parser.h"  // Reuse from your FQ pipeline
#include "hdf5_amr_writer.h"
#include "clinical_amr_report_generator.h"

namespace fs = std::filesystem;

// Function to check if file is gzipped
bool isGzipped(const std::string& filename) {
    return filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz";
}

// Function to read FASTQ file (supports gzipped)
std::vector<std::pair<std::string, std::string>> readFastq(const std::string& filename, 
                                                           int max_reads = -1) {
    std::vector<std::pair<std::string, std::string>> reads;
    
    if (isGzipped(filename)) {
        // Read gzipped file
        gzFile gz_file = gzopen(filename.c_str(), "rb");
        if (!gz_file) {
            std::cerr << "Error: Cannot open gzipped file " << filename << std::endl;
            return reads;
        }
        
        char buffer[4096];
        int count = 0;
        
        // Read FASTQ in groups of 4 lines (like resistance pipeline)
        while (gzgets(gz_file, buffer, sizeof(buffer))) {
            std::string header = buffer;
            header.pop_back(); // Remove newline
            
            // Read sequence line
            if (!gzgets(gz_file, buffer, sizeof(buffer))) break;
            std::string seq = buffer;
            seq.pop_back(); // Remove newline
            
            // Read + line
            if (!gzgets(gz_file, buffer, sizeof(buffer))) break;
            
            // Read quality line
            if (!gzgets(gz_file, buffer, sizeof(buffer))) break;
            
            // Check for empty sequence
            if (seq.empty()) {
                std::cerr << "WARNING: Empty sequence at line group starting with: " << header << std::endl;
                continue;  // Skip this read
            }
            
            // Add read (extract ID by removing @)
            if (!seq.empty() && !header.empty()) {
                std::string id = header.substr(1);  // Remove '@'
                reads.push_back({id, seq});
                
                // Debug output for first 10 reads
                if (reads.size() < 10) {
                    std::cout << "DEBUG: Read " << reads.size() << " - ID: '" << id << "', Length: " << seq.length() << std::endl;
                }
                
                count++;
                
                if (max_reads > 0 && count >= max_reads) {
                    break;
                }
            }
        }
        
        gzclose(gz_file);
    } else {
        // Read uncompressed file
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return reads;
        }
        
        std::string header, seq, plus, quality;
        int count = 0;
        
        // Read FASTQ in groups of 4 lines (like resistance pipeline)
        while (std::getline(file, header)) {
            // Read the next 3 lines
            if (std::getline(file, seq) && 
                std::getline(file, plus) && 
                std::getline(file, quality)) {
                
                // Check for empty sequence
                if (seq.empty()) {
                    std::cerr << "WARNING: Empty sequence at line group starting with: " << header << std::endl;
                    continue;  // Skip this read
                }
                
                // Add read if sequence is not empty
                if (!seq.empty() && !header.empty()) {
                    std::string id = header.substr(1);  // Remove '@'
                    reads.push_back({id, seq});
                    
                    // Debug output for first 10 reads
                    if (reads.size() < 10) {
                        std::cout << "DEBUG: Read " << reads.size() << " - ID: '" << id << "', Length: " << seq.length() << std::endl;
                    }
                    
                    count++;
                    
                    if (max_reads > 0 && count >= max_reads) {
                        break;
                    }
                }
            }
        }
        
        file.close();
    }
    
    return reads;
}

// Function to merge paired-end reads into single longer reads
std::vector<std::pair<std::string, std::string>> mergePairedReads(
    const std::vector<std::pair<std::string, std::string>>& reads1,
    const std::vector<std::pair<std::string, std::string>>& reads2,
    int merge_gap = 10) {
    
    std::vector<std::pair<std::string, std::string>> merged;
    size_t min_size = std::min(reads1.size(), reads2.size());
    
    // Simple merge strategy: concatenate R1 and R2 with N's in between
    std::string gap_seq(merge_gap, 'N');
    
    for (size_t i = 0; i < min_size; i++) {
        // Extract read ID without /1 or /2 suffix
        std::string id1 = reads1[i].first;
        std::string id2 = reads2[i].first;
        
        // Remove /1 or /2 suffixes if present (old format)
        size_t pos1 = id1.find('/');
        if (pos1 != std::string::npos) id1 = id1.substr(0, pos1);
        size_t pos2 = id2.find('/');
        if (pos2 != std::string::npos) id2 = id2.substr(0, pos2);
        
        // Handle new Illumina format (e.g., "1:N:0" vs "2:N:0")
        // Find the first space and check if it's followed by read pair info
        size_t space1 = id1.find(' ');
        size_t space2 = id2.find(' ');
        if (space1 != std::string::npos && space2 != std::string::npos) {
            // Check if the part after space starts with "1:" or "2:"
            std::string suffix1 = id1.substr(space1 + 1);
            std::string suffix2 = id2.substr(space2 + 1);
            if ((suffix1.length() >= 2 && suffix1[0] == '1' && suffix1[1] == ':') ||
                (suffix2.length() >= 2 && suffix2[0] == '2' && suffix2[1] == ':')) {
                // Remove the read pair suffix for comparison
                id1 = id1.substr(0, space1);
                id2 = id2.substr(0, space2);
            }
        }
        
        // Verify paired reads match
        if (id1 != id2) {
            std::cerr << "Warning: Read IDs don't match at position " << i 
                      << ": " << reads1[i].first << " vs " << reads2[i].first << std::endl;
            continue;
        }
        
        // Merge sequences
        std::string merged_seq = reads1[i].second + gap_seq + reads2[i].second;
        merged.push_back({id1 + "_merged", merged_seq});
    }
    
    return merged;
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
    
    // Create clinical report generator
    std::string report_prefix = output_dir + "/" + sample_name;
    ClinicalAMRReportGenerator report_generator(report_prefix, sample_name);
    
    // Read FASTQ files
    std::cout << "Reading R1 file: " << fastq_r1 << std::endl;
    auto reads_r1 = readFastq(fastq_r1);
    
    std::cout << "Reading R2 file: " << fastq_r2 << std::endl;
    auto reads_r2 = readFastq(fastq_r2);
    
    if (reads_r1.empty() || reads_r2.empty()) {
        std::cerr << "No reads found in one or both files" << std::endl;
        return;
    }
    
    std::cout << "R1 reads: " << reads_r1.size() << ", R2 reads: " << reads_r2.size() << std::endl;
    std::cout << "Processing paired-end reads separately..." << std::endl;
    
    // Collect all AMR hits across batches
    std::vector<AMRHit> all_amr_hits;
    
    // Process in batches
    const int batch_size = config.reads_per_batch;
    int num_pairs = std::min(reads_r1.size(), reads_r2.size());
    int num_batches = (num_pairs + batch_size - 1) / batch_size;
    
    for (int batch = 0; batch < num_batches; batch++) {
        int start_idx = batch * batch_size;
        int end_idx = std::min(start_idx + batch_size, num_pairs);
        
        std::cout << "\n=== BATCH TRANSITION DEBUG ===" << std::endl;
        std::cout << "Starting batch " << (batch + 1) << "/" << num_batches << std::endl;
        std::cout << "Batch range: pairs " << start_idx << " to " << end_idx << std::endl;
        
        // Extract reads and IDs for this batch
        std::vector<std::string> batch_reads1;
        std::vector<std::string> batch_reads2;
        std::vector<std::string> batch_ids;
        
        std::cout << "Extracting reads for batch..." << std::endl;
        
        for (int i = start_idx; i < end_idx; i++) {
            // Skip pairs where either read is empty
            if (!reads_r1[i].second.empty() && !reads_r2[i].second.empty()) {
                batch_ids.push_back(reads_r1[i].first);
                batch_reads1.push_back(reads_r1[i].second);
                batch_reads2.push_back(reads_r2[i].second);
            }
        }
        
        std::cout << "Batch " << (batch + 1) << " extracted: " 
                  << batch_reads1.size() << " R1, " 
                  << batch_reads2.size() << " R2, " 
                  << batch_ids.size() << " IDs" << std::endl;
        
        // Validate first few reads in batch
        for (int i = 0; i < std::min(3, (int)batch_reads1.size()); i++) {
            std::cout << "  R1[" << i << "]: length=" << batch_reads1[i].length() << std::endl;
            std::cout << "  R2[" << i << "]: length=" << batch_reads2[i].length() << std::endl;
        }
        
        std::cout << "Calling processBatchPaired..." << std::endl;
        
        // Process batch of paired reads
        pipeline.processBatchPaired(batch_reads1, batch_reads2, batch_ids);
        
        std::cout << "Batch " << (batch + 1) << " completed successfully" << std::endl;
        std::cout << "=== END BATCH TRANSITION DEBUG ===" << std::endl;
        
        // Get results from this batch
        auto batch_hits = pipeline.getAMRHits();
        
        // Accumulate hits
        all_amr_hits.insert(all_amr_hits.end(), batch_hits.begin(), batch_hits.end());
        
        // Write batch hits to HDF5
        hdf5_writer.addAMRHits(batch_hits);
    }
    
    // Get final coverage statistics and gene entries
    auto coverage_stats = pipeline.getCoverageStats();
    auto gene_entries = pipeline.getGeneEntries();  // Now using the new getter method
    
    // Write coverage stats to HDF5
    hdf5_writer.addCoverageStats(coverage_stats, gene_entries);
    
    // Finalize HDF5
    std::string json_summary = output_dir + "/" + sample_name + "_hdf5_summary.json";
    hdf5_writer.finalize(json_summary);
    
    // Generate clinical report with all accumulated hits
    report_generator.processAMRResults(all_amr_hits, coverage_stats, 
                                      gene_entries, num_pairs);
    report_generator.generateReports();
    
    // Write basic TSV results (backward compatibility)
    std::string output_prefix = output_dir + "/" + sample_name;
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
    
    // Create clinical report generator
    std::string report_prefix = output_dir + "/" + sample_name;
    ClinicalAMRReportGenerator report_generator(report_prefix, sample_name);
    
    // Read FASTQ file
    std::cout << "Reading FASTQ file: " << fastq_file << std::endl;
    auto fastq_data = readFastq(fastq_file);
    
    if (fastq_data.empty()) {
        std::cerr << "No reads found in " << fastq_file << std::endl;
        return;
    }
    
    std::cout << "Total reads: " << fastq_data.size() << std::endl;
    
    // Collect all AMR hits
    std::vector<AMRHit> all_amr_hits;
    
    // Process in batches
    const int batch_size = config.reads_per_batch;
    int num_batches = (fastq_data.size() + batch_size - 1) / batch_size;
    
    for (int batch = 0; batch < num_batches; batch++) {
        int start_idx = batch * batch_size;
        int end_idx = std::min(start_idx + batch_size, (int)fastq_data.size());
        
        std::cout << "\nProcessing batch " << (batch + 1) << "/" << num_batches 
                  << " (reads " << start_idx << "-" << end_idx << ")" << std::endl;
        
        // Extract reads and IDs for this batch
        std::vector<std::string> batch_reads;
        std::vector<std::string> batch_ids;
        
        for (int i = start_idx; i < end_idx; i++) {
            // Skip empty reads
            if (!fastq_data[i].second.empty()) {
                batch_ids.push_back(fastq_data[i].first);
                batch_reads.push_back(fastq_data[i].second);
            }
        }
        
        // Process batch
        pipeline.processBatch(batch_reads, batch_ids);
        
        // Get results from this batch
        auto batch_hits = pipeline.getAMRHits();
        
        // Accumulate hits
        all_amr_hits.insert(all_amr_hits.end(), batch_hits.begin(), batch_hits.end());
        
        // Write batch hits to HDF5
        hdf5_writer.addAMRHits(batch_hits);
    }
    
    // Get final coverage statistics and gene entries
    auto coverage_stats = pipeline.getCoverageStats();
    auto gene_entries = pipeline.getGeneEntries();
    
    // Write coverage stats to HDF5
    hdf5_writer.addCoverageStats(coverage_stats, gene_entries);
    
    // Finalize HDF5
    std::string json_summary = output_dir + "/" + sample_name + "_hdf5_summary.json";
    hdf5_writer.finalize(json_summary);
    
    // Generate clinical report
    report_generator.processAMRResults(all_amr_hits, coverage_stats, 
                                      gene_entries, fastq_data.size());
    report_generator.generateReports();
    
    // Write basic TSV results (backward compatibility)
    std::string output_prefix = output_dir + "/" + sample_name;
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
    
    // REDUCE BATCH SIZE TO AVOID ISSUES
    config.reads_per_batch = 50000;  // Reduced from 100000
    
    config.max_read_length = 300;
    config.output_prefix = output_dir + "/amr_results";
    
    // Additional configuration for paired-end processing
    bool merge_paired_reads = true;
    
    // Optional: Parse additional config from command line
    for (int i = 4; i < argc; i += 2) {
        std::string arg = argv[i];
        if (i + 1 < argc) {
            if (arg == "--min-identity") {
                config.min_identity = std::stof(argv[i + 1]);
            } else if (arg == "--min-coverage") {
                config.min_coverage = std::stof(argv[i + 1]);
            } else if (arg == "--kmer-length") {
                config.kmer_length = std::stoi(argv[i + 1]);
            } else if (arg == "--batch-size") {
                config.reads_per_batch = std::stoi(argv[i + 1]);
            } else if (arg == "--protein-db") {
                config.protein_db_path = argv[i + 1];
            } else if (arg == "--no-merge") {
                merge_paired_reads = false;
                i--; // This flag doesn't take a value
            } else if (arg == "--use-bloom-filter") {
                config.use_bloom_filter = true;
                i--; // This flag doesn't take a value
            }
        }
    }
    
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
        
        if (sample.isPairedEnd()) {
            // Process paired-end reads
            processSamplePaired(pipeline, sample.sample_name, 
                               sample.read1_path, sample.read2_path, output_dir, config, merge_paired_reads);
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
    
    // Check for validation mode
    bool validate_mode = false;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--validate") {
            validate_mode = true;
            break;
        }
    }
    
    if (validate_mode) {
        std::cout << "\n=== VALIDATION MODE ===" << std::endl;
        std::cout << "Please verify the following in the output files:" << std::endl;
        std::cout << "1. Gene families are correctly extracted (e.g., blaKPC-2 -> blaKPC)" << std::endl;
        std::cout << "2. Gene IDs match the FASTA headers" << std::endl;
        std::cout << "3. Multiple variants of same family are grouped together" << std::endl;
        std::cout << "4. TSV files contain gene_family column" << std::endl;
        std::cout << "5. HDF5 files contain gene_family dataset" << std::endl;
        std::cout << "\nCheck the following output files:" << std::endl;
        std::cout << "  - " << output_dir << "_amr_abundance.tsv" << std::endl;
        std::cout << "  - " << output_dir << "_gene_family_summary.tsv" << std::endl;
        std::cout << "  - " << output_dir << "_clinical_amr_report.html" << std::endl;
        std::cout << "  - " << output_dir << "_clinical_amr_report.json" << std::endl;
        std::cout << "  - " << output_dir << ".h5 (HDF5 output)" << std::endl;
    }
    
    return 0;
}
