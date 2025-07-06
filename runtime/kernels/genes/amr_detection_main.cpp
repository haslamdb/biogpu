// amr_detection_main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include "amr_detection_pipeline.h"
#include "sample_csv_parser.h"  // Reuse from your FQ pipeline

// Function to read FASTQ file
std::vector<std::pair<std::string, std::string>> readFastq(const std::string& filename, 
                                                           int max_reads = -1) {
    std::vector<std::pair<std::string, std::string>> reads;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return reads;
    }
    
    std::string line;
    int count = 0;
    
    while (std::getline(file, line)) {
        if (line[0] == '@') {
            std::string id = line.substr(1);  // Remove '@'
            std::string seq;
            std::getline(file, seq);
            
            // Skip quality header and quality string
            std::getline(file, line);  // +
            std::getline(file, line);  // quality
            
            reads.push_back({id, seq});
            count++;
            
            if (max_reads > 0 && count >= max_reads) {
                break;
            }
        }
    }
    
    file.close();
    return reads;
}

// Function to process a single sample
void processSample(AMRDetectionPipeline& pipeline, 
                  const std::string& sample_name,
                  const std::string& fastq_file,
                  const std::string& output_dir) {
    
    std::cout << "\n=== Processing sample: " << sample_name << " ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Read FASTQ file
    std::cout << "Reading FASTQ file: " << fastq_file << std::endl;
    auto fastq_data = readFastq(fastq_file);
    
    if (fastq_data.empty()) {
        std::cerr << "No reads found in " << fastq_file << std::endl;
        return;
    }
    
    std::cout << "Total reads: " << fastq_data.size() << std::endl;
    
    // Process in batches
    const int batch_size = 100000;  // From config
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
            batch_ids.push_back(fastq_data[i].first);
            batch_reads.push_back(fastq_data[i].second);
        }
        
        // Process batch
        pipeline.processBatch(batch_reads, batch_ids);
    }
    
    // Write results
    std::string output_prefix = output_dir + "/" + sample_name;
    pipeline.writeResults(output_prefix);
    pipeline.generateClinicalReport(output_prefix + "_clinical_report.txt");
    pipeline.exportAbundanceTable(output_prefix + "_abundance.tsv");
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "\nSample processing completed in " << duration.count() << " seconds" << std::endl;
}

// Main function
int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <amr_db_path> <input_csv> <output_dir>" << std::endl;
        std::cerr << "  amr_db_path: Path to AMR database (format: dna.fasta,protein.fasta)" << std::endl;
        std::cerr << "  input_csv: CSV file with columns: sample_name,fastq_path" << std::endl;
        std::cerr << "  output_dir: Directory for output files" << std::endl;
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
    config.min_identity = 0.85f;
    config.min_coverage = 0.80f;
    config.min_alignment_length = 50;
    config.band_width = 15;
    config.reads_per_batch = 100000;
    config.max_read_length = 300;
    config.output_prefix = output_dir + "/amr_results";
    
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
    
    // Parse input CSV
    std::cout << "\nParsing input CSV..." << std::endl;
    std::vector<std::pair<std::string, std::string>> samples;
    
    std::ifstream csv_file(input_csv);
    if (!csv_file.is_open()) {
        std::cerr << "Cannot open input CSV: " << input_csv << std::endl;
        return 1;
    }
    
    std::string line;
    // Skip header if present
    std::getline(csv_file, line);
    if (line.find("sample_name") == std::string::npos) {
        // No header, process first line
        size_t comma_pos = line.find(',');
        if (comma_pos != std::string::npos) {
            std::string sample_name = line.substr(0, comma_pos);
            std::string fastq_path = line.substr(comma_pos + 1);
            samples.push_back({sample_name, fastq_path});
        }
    }
    
    // Read remaining lines
    while (std::getline(csv_file, line)) {
        size_t comma_pos = line.find(',');
        if (comma_pos != std::string::npos) {
            std::string sample_name = line.substr(0, comma_pos);
            std::string fastq_path = line.substr(comma_pos + 1);
            
            // Remove quotes if present
            if (fastq_path.front() == '"') fastq_path = fastq_path.substr(1);
            if (fastq_path.back() == '"') fastq_path.pop_back();
            
            samples.push_back({sample_name, fastq_path});
        }
    }
    csv_file.close();
    
    std::cout << "Found " << samples.size() << " patient samples to analyze" << std::endl;
    
    // Process each patient sample
    auto total_start = std::chrono::high_resolution_clock::now();
    
    for (size_t i = 0; i < samples.size(); i++) {
        std::cout << "\n[" << (i + 1) << "/" << samples.size() << "] ";
        processSample(pipeline, samples[i].first, samples[i].second, output_dir);
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
