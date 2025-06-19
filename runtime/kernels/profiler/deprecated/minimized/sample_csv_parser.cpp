// sample_csv_parser.cpp
// Implementation of CSV parser for batch sample processing

#include "sample_csv_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <functional>
#include <iomanip>
#include <chrono>

namespace BioGPU {

// Constructor
SampleCSVParser::SampleCSVParser() : has_header(true), delimiter(',') {
}

// Trim whitespace from string
std::string SampleCSVParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
        return "";
    }
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// Split string by delimiter
std::vector<std::string> SampleCSVParser::split(const std::string& line, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(line);
    
    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(trim(token));
    }
    
    return tokens;
}

// Check if column name matches any variant
bool SampleCSVParser::matchesColumnName(const std::string& column, 
                                       const std::vector<std::string>& variants) {
    std::string col_lower = column;
    std::transform(col_lower.begin(), col_lower.end(), col_lower.begin(), ::tolower);
    
    for (const auto& variant : variants) {
        std::string var_lower = variant;
        std::transform(var_lower.begin(), var_lower.end(), var_lower.begin(), ::tolower);
        if (col_lower == var_lower) {
            return true;
        }
    }
    return false;
}

// Expand ~ in path
std::string SampleCSVParser::expandPath(const std::string& path) const {
    if (path.empty() || path[0] != '~') {
        return path;
    }
    
    const char* home = getenv("HOME");
    if (!home) {
        return path;
    }
    
    return std::string(home) + path.substr(1);
}

// Strip path from filename
std::string SampleCSVParser::stripPath(const std::string& filepath) const {
    size_t pos = filepath.find_last_of("/\\");
    if (pos != std::string::npos) {
        return filepath.substr(pos + 1);
    }
    return filepath;
}

// Combine directory and filename
std::string SampleCSVParser::combinePath(const std::string& dir, const std::string& filename) const {
    if (dir.empty()) {
        return filename;
    }
    
    std::string expanded_dir = expandPath(dir);
    
    // Ensure directory ends with separator
    if (expanded_dir.back() != '/' && expanded_dir.back() != '\\') {
        expanded_dir += '/';
    }
    
    // Strip any leading path from filename
    std::string clean_filename = stripPath(filename);
    
    return expanded_dir + clean_filename;
}

// Check if file exists
bool SampleCSVParser::fileExists(const std::string& path) const {
    struct stat buffer;
    std::string expanded = expandPath(path);
    return (stat(expanded.c_str(), &buffer) == 0);
}

// Auto-detect delimiter
void SampleCSVParser::detectDelimiter(const std::string& first_line) {
    // Count occurrences of common delimiters
    int comma_count = std::count(first_line.begin(), first_line.end(), ',');
    int tab_count = std::count(first_line.begin(), first_line.end(), '\t');
    int semicolon_count = std::count(first_line.begin(), first_line.end(), ';');
    
    // Choose the most common one
    if (tab_count > comma_count && tab_count > semicolon_count) {
        delimiter = '\t';
    } else if (semicolon_count > comma_count) {
        delimiter = ';';
    } else {
        delimiter = ',';
    }
}

// Parse CSV file
bool SampleCSVParser::parseFile(const std::string& csv_path, bool validate_paths) {
    samples.clear();
    header_columns.clear();
    
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open CSV file: " << csv_path << std::endl;
        return false;
    }
    
    std::string line;
    bool first_line = true;
    int line_number = 0;
    
    // Column indices
    int sample_name_idx = -1;
    int file_path_idx = -1;
    int read1_idx = -1;
    int read2_idx = -1;
    
    while (std::getline(file, line)) {
        line_number++;
        
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\n\r") == std::string::npos) {
            continue;
        }
        
        // Auto-detect delimiter on first line
        if (first_line) {
            detectDelimiter(line);
        }
        
        std::vector<std::string> fields = split(line, delimiter);
        
        // Process header
        if (first_line && has_header) {
            first_line = false;
            header_columns = fields;
            
            // Find column indices
            for (size_t i = 0; i < fields.size(); i++) {
                if (matchesColumnName(fields[i], SAMPLE_NAME_VARIANTS)) {
                    sample_name_idx = i;
                } else if (matchesColumnName(fields[i], FILE_PATH_VARIANTS)) {
                    file_path_idx = i;
                } else if (matchesColumnName(fields[i], READ1_VARIANTS)) {
                    read1_idx = i;
                } else if (matchesColumnName(fields[i], READ2_VARIANTS)) {
                    read2_idx = i;
                }
            }
            
            // Validate required columns
            if (sample_name_idx == -1) {
                std::cerr << "Error: Could not find sample name column. Expected one of: ";
                for (const auto& var : SAMPLE_NAME_VARIANTS) {
                    std::cerr << "'" << var << "' ";
                }
                std::cerr << std::endl;
                return false;
            }
            
            if (file_path_idx == -1) {
                std::cerr << "Error: Could not find file path column. Expected one of: ";
                for (const auto& var : FILE_PATH_VARIANTS) {
                    std::cerr << "'" << var << "' ";
                }
                std::cerr << std::endl;
                return false;
            }
            
            if (read1_idx == -1) {
                std::cerr << "Error: Could not find R1 file column. Expected one of: ";
                for (const auto& var : READ1_VARIANTS) {
                    std::cerr << "'" << var << "' ";
                }
                std::cerr << std::endl;
                return false;
            }
            
            std::cout << "Detected delimiter: '" << delimiter << "'" << std::endl;
            std::cout << "Column mapping:" << std::endl;
            std::cout << "  Sample name: column " << (sample_name_idx + 1) 
                      << " ('" << fields[sample_name_idx] << "')" << std::endl;
            std::cout << "  File path: column " << (file_path_idx + 1) 
                      << " ('" << fields[file_path_idx] << "')" << std::endl;
            std::cout << "  R1 file: column " << (read1_idx + 1) 
                      << " ('" << fields[read1_idx] << "')" << std::endl;
            if (read2_idx != -1) {
                std::cout << "  R2 file: column " << (read2_idx + 1) 
                          << " ('" << fields[read2_idx] << "')" << std::endl;
            }
            
            continue;
        }
        
        // If no header, assume first four columns are sample, path, r1, r2
        if (!has_header && first_line) {
            first_line = false;
            sample_name_idx = 0;
            file_path_idx = 1;
            read1_idx = 2;
            read2_idx = (fields.size() > 3) ? 3 : -1;
        }
        
        // Parse sample data
        if (fields.size() > sample_name_idx && fields.size() > file_path_idx && fields.size() > read1_idx) {
            SampleInfo sample;
            sample.sample_name = fields[sample_name_idx];
            sample.file_path = fields[file_path_idx];
            sample.read1_filename = stripPath(fields[read1_idx]);
            
            // Construct full paths
            sample.read1_path = combinePath(sample.file_path, sample.read1_filename);
            
            if (read2_idx != -1 && fields.size() > read2_idx && !fields[read2_idx].empty()) {
                sample.read2_filename = stripPath(fields[read2_idx]);
                sample.read2_path = combinePath(sample.file_path, sample.read2_filename);
            }
            
            // Store any additional columns as metadata
            for (size_t i = 0; i < fields.size(); i++) {
                if (i != sample_name_idx && i != file_path_idx && i != read1_idx && i != read2_idx && !header_columns.empty()) {
                    if (i < header_columns.size()) {
                        sample.metadata[header_columns[i]] = fields[i];
                    }
                }
            }
            
            // Validate paths if requested
            if (validate_paths) {
                if (!fileExists(sample.read1_path)) {
                    std::cerr << "Warning: R1 file not found for sample '" << sample.sample_name 
                              << "': " << sample.read1_path << std::endl;
                }
                if (!sample.read2_path.empty() && !fileExists(sample.read2_path)) {
                    std::cerr << "Warning: R2 file not found for sample '" << sample.sample_name 
                              << "': " << sample.read2_path << std::endl;
                }
            }
            
            samples.push_back(sample);
        } else {
            std::cerr << "Warning: Skipping incomplete line " << line_number 
                      << " (expected at least " << std::max({sample_name_idx, file_path_idx, read1_idx}) + 1 
                      << " fields, got " << fields.size() << ")" << std::endl;
        }
    }
    
    file.close();
    
    std::cout << "Parsed " << samples.size() << " samples from CSV" << std::endl;
    return !samples.empty();
}

// Find sample by name
const SampleInfo* SampleCSVParser::findSampleByName(const std::string& name) const {
    for (const auto& sample : samples) {
        if (sample.sample_name == name) {
            return &sample;
        }
    }
    return nullptr;
}

// Validate sample
bool SampleCSVParser::validateSample(const SampleInfo& sample) const {
    if (!sample.isValid()) {
        return false;
    }
    
    if (!fileExists(sample.read1_path)) {
        return false;
    }
    
    if (!sample.read2_path.empty() && !fileExists(sample.read2_path)) {
        return false;
    }
    
    return true;
}

// Get validation errors
std::vector<std::string> SampleCSVParser::getValidationErrors() const {
    std::vector<std::string> errors;
    
    for (const auto& sample : samples) {
        if (!sample.isValid()) {
            errors.push_back("Sample '" + sample.sample_name + "' is missing required fields");
        } else {
            if (!fileExists(sample.read1_path)) {
                errors.push_back("Sample '" + sample.sample_name + "': Read1 file not found: " + sample.read1_path);
            }
            if (!sample.read2_path.empty() && !fileExists(sample.read2_path)) {
                errors.push_back("Sample '" + sample.sample_name + "': Read2 file not found: " + sample.read2_path);
            }
        }
    }
    
    return errors;
}

// Print summary
void SampleCSVParser::printSummary() const {
    std::cout << "\n=== Sample Summary ===" << std::endl;
    std::cout << "Total samples: " << samples.size() << std::endl;
    
    int paired_end = 0;
    int single_end = 0;
    int valid = 0;
    
    for (const auto& sample : samples) {
        if (sample.isPairedEnd()) {
            paired_end++;
        } else {
            single_end++;
        }
        if (validateSample(sample)) {
            valid++;
        }
    }
    
    std::cout << "Paired-end: " << paired_end << std::endl;
    std::cout << "Single-end: " << single_end << std::endl;
    std::cout << "Valid (files exist): " << valid << std::endl;
    
    if (valid < samples.size()) {
        std::cout << "Invalid: " << (samples.size() - valid) << std::endl;
    }
}

// Print detailed information
void SampleCSVParser::printDetailed() const {
    std::cout << "\n=== Detailed Sample Information ===" << std::endl;
    
    for (size_t i = 0; i < samples.size(); i++) {
        const auto& sample = samples[i];
        std::cout << "\nSample " << (i + 1) << ":" << std::endl;
        std::cout << "  Name: " << sample.sample_name << std::endl;
        std::cout << "  Path: " << sample.file_path << std::endl;
        std::cout << "  R1 file: " << sample.read1_filename 
                  << (fileExists(sample.read1_path) ? " [EXISTS]" : " [NOT FOUND]") << std::endl;
        std::cout << "    Full path: " << sample.read1_path << std::endl;
        
        if (!sample.read2_filename.empty()) {
            std::cout << "  R2 file: " << sample.read2_filename 
                      << (fileExists(sample.read2_path) ? " [EXISTS]" : " [NOT FOUND]") << std::endl;
            std::cout << "    Full path: " << sample.read2_path << std::endl;
        }
        
        if (!sample.metadata.empty()) {
            std::cout << "  Additional metadata:" << std::endl;
            for (const auto& kv : sample.metadata) {
                std::cout << "    " << kv.first << ": " << kv.second << std::endl;
            }
        }
    }
}

// BatchProcessor implementation
BatchProcessor::BatchProcessor(const std::string& output_dir, bool create_dirs) 
    : output_base_dir(output_dir), create_sample_dirs(create_dirs) {
}

bool BatchProcessor::loadSamples(const std::string& csv_path) {
    return parser.parseFile(csv_path);
}

std::string BatchProcessor::getOutputPath(const SampleInfo& sample, const std::string& suffix) const {
    std::string path = output_base_dir;
    
    if (create_sample_dirs) {
        path += "/" + sample.sample_name;
        
        // Create directory if needed
        std::string mkdir_cmd = "mkdir -p " + path;
        system(mkdir_cmd.c_str());
    }
    
    path += "/" + sample.sample_name;
    if (!suffix.empty()) {
        path += "_" + suffix;
    }
    
    return path;
}

int BatchProcessor::processBatch(ProcessFunc process_func, bool stop_on_error) {
    int successful = 0;
    int failed = 0;
    
    std::cout << "\n=== Starting Batch Processing ===" << std::endl;
    std::cout << "Processing " << parser.getSampleCount() << " samples" << std::endl;
    
    for (size_t i = 0; i < parser.getSampleCount(); i++) {
        const SampleInfo* sample = parser.getSample(i);
        if (!sample) continue;
        
        std::cout << "\n[" << (i + 1) << "/" << parser.getSampleCount() << "] "
                  << "Processing sample: " << sample->sample_name << std::endl;
        
        if (!parser.validateSample(*sample)) {
            std::cerr << "  ERROR: Invalid sample (files missing)" << std::endl;
            failed++;
            if (stop_on_error) {
                break;
            }
            continue;
        }
        
        std::string output_path = getOutputPath(*sample);
        
        auto start_time = std::chrono::high_resolution_clock::now();
        int result = process_func(*sample, output_path);
        auto end_time = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        if (result == 0) {
            std::cout << "  SUCCESS - Completed in " << duration.count() << " seconds" << std::endl;
            successful++;
        } else {
            std::cerr << "  FAILED - Error code: " << result << std::endl;
            failed++;
            if (stop_on_error) {
                break;
            }
        }
    }
    
    std::cout << "\n=== Batch Processing Complete ===" << std::endl;
    std::cout << "Successful: " << successful << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    return failed;
}

} // namespace BioGPU