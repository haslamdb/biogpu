#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <memory>
#include <iomanip>

// Include common structures
#include "paired_read_common.h"

// Adaptive paired-end profiler that automatically chooses between
// direct GPU loading and streaming based on database size analysis

struct DatabaseStats {
    uint64_t estimated_gpu_memory_mb;
    uint32_t num_organisms;
    uint64_t total_minimizer_entries;
    bool needs_streaming;
    uint32_t recommended_chunks;
};

class AdaptivePairedEndProfiler {
private:
    // Database metadata
    std::unordered_map<uint32_t, std::string> organism_names;
    std::unordered_map<uint32_t, std::string> organism_taxonomy;
    std::unordered_map<uint32_t, uint64_t> organism_genome_sizes;
    std::unordered_map<uint32_t, uint32_t> organism_expected_minimizers;
    std::vector<uint32_t> organism_id_list;
    
    // Analysis results
    DatabaseStats db_stats;
    bool use_streaming_mode = false;
    
    // Actual GPU profiler implementations
    // std::unique_ptr<GPUMicrobialCommunityProfiler> gpu_profiler;
    // std::unique_ptr<StreamingGPUProfiler> streaming_profiler;
    
    int k = 35, m = 31;
    
public:
    AdaptivePairedEndProfiler() {
        std::cout << "Adaptive Paired-End GPU Profiler" << std::endl;
        std::cout << "Automatically selects optimal processing strategy" << std::endl;
    }
    
    void load_and_analyze_database(const std::string& db_path) {
        std::cout << "Analyzing database: " << db_path << std::endl;
        
        // First, analyze database size and characteristics
        db_stats = analyze_database_requirements(db_path);
        
        // Decide on processing strategy
        use_streaming_mode = db_stats.needs_streaming;
        
        std::cout << "\n=== DATABASE ANALYSIS RESULTS ===" << std::endl;
        std::cout << "Estimated GPU memory needed: " << db_stats.estimated_gpu_memory_mb << " MB" << std::endl;
        std::cout << "Number of organisms: " << db_stats.num_organisms << std::endl;
        std::cout << "Total minimizer entries: " << db_stats.total_minimizer_entries << std::endl;
        
        if (use_streaming_mode) {
            std::cout << "🔄 STREAMING MODE SELECTED" << std::endl;
            std::cout << "Database too large for direct GPU loading" << std::endl;
            std::cout << "Will use " << db_stats.recommended_chunks << " chunks" << std::endl;
            initialize_streaming_mode(db_path);
        } else {
            std::cout << "⚡ DIRECT GPU MODE SELECTED" << std::endl;
            std::cout << "Database fits in GPU memory - maximum performance" << std::endl;
            initialize_direct_mode(db_path);
        }
        
        // Load organism metadata regardless of mode
        load_organism_metadata(db_path);
    }
    
private:
    DatabaseStats analyze_database_requirements(const std::string& db_path) {
        DatabaseStats stats = {};
        
        std::ifstream file(db_path, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open database file for analysis");
        }
        
        // Get file size
        uint64_t file_size_mb = file.tellg() / (1024 * 1024);
        file.seekg(0, std::ios::beg);
        
        // Read header
        uint32_t magic, version, k_size, m_size;
        uint64_t num_minimizer_hashes;
        
        file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        file.read(reinterpret_cast<char*>(&version), sizeof(version));
        file.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        file.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        file.read(reinterpret_cast<char*>(&stats.num_organisms), sizeof(stats.num_organisms));
        file.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        // Accept both little-endian (MINI) and big-endian (INIM) formats
        if (magic != 0x4D494E49 && magic != 0x494E494D) {
            throw std::runtime_error("Invalid database format");
        }
        
        // Estimate total minimizer entries (sample first few hashes)
        skip_organism_metadata(file, stats.num_organisms);
        
        uint64_t sample_entries = 0;
        uint64_t sample_hashes = 0;
        const uint64_t max_sample = 1000;  // Sample first 1000 hashes
        
        uint64_t hash;
        uint32_t num_entries_for_hash;
        
        while (sample_hashes < max_sample && 
               file.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
            file.read(reinterpret_cast<char*>(&num_entries_for_hash), sizeof(num_entries_for_hash));
            
            sample_entries += num_entries_for_hash;
            sample_hashes++;
            
            // Skip the actual entries
            file.seekg(num_entries_for_hash * (sizeof(uint64_t) + sizeof(uint32_t) + sizeof(uint8_t)), 
                      std::ios::cur);
        }
        
        // Estimate total entries
        if (sample_hashes > 0) {
            double avg_entries_per_hash = (double)sample_entries / sample_hashes;
            stats.total_minimizer_entries = (uint64_t)(avg_entries_per_hash * num_minimizer_hashes);
        } else {
            // Fallback: estimate based on file size
            stats.total_minimizer_entries = file_size_mb * 50000;  // Rough estimate
        }
        
        file.close();
        
        // Calculate GPU memory requirements
        size_t memory_per_entry = sizeof(uint64_t) + sizeof(uint32_t) + sizeof(float);  // 16 bytes
        stats.estimated_gpu_memory_mb = (stats.total_minimizer_entries * memory_per_entry) / (1024 * 1024);
        
        // Add overhead for organism arrays, reads, and working memory
        stats.estimated_gpu_memory_mb += 1000;  // 1GB overhead
        
        // Decision: use streaming if estimated memory > 30GB (leaving 18GB for reads/working memory on A6000)
        size_t available_memory_mb = 30 * 1024;  // 30GB threshold
        stats.needs_streaming = stats.estimated_gpu_memory_mb > available_memory_mb;
        
        if (stats.needs_streaming) {
            // Calculate optimal chunk configuration
            size_t target_chunk_size_mb = 4 * 1024;  // 4GB chunks
            stats.recommended_chunks = (stats.estimated_gpu_memory_mb / target_chunk_size_mb) + 1;
        } else {
            stats.recommended_chunks = 1;
        }
        
        return stats;
    }
    
    void skip_organism_metadata(std::ifstream& file, uint32_t num_organisms) {
        for (uint32_t i = 0; i < num_organisms; i++) {
            // Skip fixed-size fields
            file.seekg(sizeof(uint32_t) + sizeof(uint32_t) + sizeof(float) + 
                      sizeof(uint64_t) + sizeof(uint32_t), std::ios::cur);
            
            // Skip variable-length strings
            uint16_t name_length, taxonomy_length;
            file.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            file.seekg(name_length, std::ios::cur);
            file.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            file.seekg(taxonomy_length, std::ios::cur);
        }
    }
    
    void load_organism_metadata(const std::string& db_path) {
        std::ifstream in(db_path, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open database for metadata loading");
        }
        
        // Skip header
        uint32_t magic, version, k_size, m_size, num_organisms;
        uint64_t num_minimizer_hashes;
        
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        in.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        in.read(reinterpret_cast<char*>(&num_organisms), sizeof(num_organisms));
        in.read(reinterpret_cast<char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        k = k_size;
        m = m_size;
        
        // Read organism metadata
        for (uint32_t i = 0; i < num_organisms; i++) {
            uint32_t taxonomy_id, taxon_level, minimizer_count;
            uint64_t genome_size;
            float gc_content;
            
            in.read(reinterpret_cast<char*>(&taxonomy_id), sizeof(taxonomy_id));
            in.read(reinterpret_cast<char*>(&taxon_level), sizeof(taxon_level));
            in.read(reinterpret_cast<char*>(&gc_content), sizeof(gc_content));
            in.read(reinterpret_cast<char*>(&genome_size), sizeof(genome_size));
            in.read(reinterpret_cast<char*>(&minimizer_count), sizeof(minimizer_count));
            
            // Read strings
            uint16_t name_length, taxonomy_length;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            std::string name(name_length, '\0');
            in.read(&name[0], name_length);
            
            in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            std::string taxonomy_path(taxonomy_length, '\0');
            in.read(&taxonomy_path[0], taxonomy_length);
            
            // Store metadata
            organism_names[taxonomy_id] = name;
            organism_taxonomy[taxonomy_id] = taxonomy_path;
            organism_genome_sizes[taxonomy_id] = genome_size;
            organism_expected_minimizers[taxonomy_id] = minimizer_count;
            organism_id_list.push_back(taxonomy_id);
        }
        
        std::sort(organism_id_list.begin(), organism_id_list.end());
        
        in.close();
        
        std::cout << "Loaded metadata for " << organism_names.size() << " organisms" << std::endl;
    }
    
    void initialize_direct_mode(const std::string& db_path) {
        std::cout << "Initializing direct GPU mode..." << std::endl;
        // Here you would initialize your direct GPU implementation
        // For now, we'll use a placeholder
        
        // Check GPU memory
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "GPU memory available: " << free_mem / (1024*1024*1024) << " GB" << std::endl;
        
        // direct_impl = std::make_unique<DirectGPUImpl>(db_path);
        std::cout << "Direct GPU mode ready" << std::endl;
    }
    
    void initialize_streaming_mode(const std::string& db_path) {
        std::cout << "Initializing streaming mode..." << std::endl;
        std::cout << "Will process database in " << db_stats.recommended_chunks << " chunks" << std::endl;
        
        // streaming_impl = std::make_unique<StreamingImpl>(db_path, db_stats.recommended_chunks);
        std::cout << "Streaming mode ready" << std::endl;
    }
    
public:
    // Paired-end FASTQ loading methods (same as before)
    std::vector<PairedRead> load_paired_fastq(const std::string& r1_path, 
                                             const std::string& r2_path = "") {
        std::cout << "Loading paired-end reads..." << std::endl;
        
        if (r2_path.empty()) {
            return load_interleaved_fastq(r1_path);
        } else {
            return load_separate_fastq_files(r1_path, r2_path);
        }
    }
    
private:
    std::vector<PairedRead> load_separate_fastq_files(const std::string& r1_path, 
                                                     const std::string& r2_path) {
        std::cout << "Loading paired-end reads from separate files" << std::endl;
        std::vector<PairedRead> paired_reads;
        
        std::ifstream r1_file(r1_path);
        std::ifstream r2_file(r2_path);
        
        if (!r1_file.is_open() || !r2_file.is_open()) {
            throw std::runtime_error("Cannot open paired-end FASTQ files");
        }
        
        std::string r1_line, r2_line;
        int line_count = 0;
        
        while (std::getline(r1_file, r1_line) && std::getline(r2_file, r2_line)) {
            line_count++;
            int line_type = line_count % 4;
            
            if (line_type == 1) {
                // Header lines
                std::string read_id = extract_read_id(r1_line);
                
                // Read sequences
                std::getline(r1_file, r1_line);
                std::getline(r2_file, r2_line);
                
                if (!r1_line.empty() && !r2_line.empty() && 
                    r1_line.length() >= k && r2_line.length() >= k) {
                    PairedRead pair(r1_line, r2_line, read_id);
                    paired_reads.push_back(pair);
                }
                
                // Skip quality lines
                std::getline(r1_file, r1_line);
                std::getline(r2_file, r2_line);
                std::getline(r1_file, r1_line);
                std::getline(r2_file, r2_line);
                
                line_count += 3;
            }
            
            if (paired_reads.size() % 100000 == 0) {
                std::cout << "\rLoaded " << paired_reads.size() << " paired reads..." << std::flush;
            }
            
            // Testing limit
            if (paired_reads.size() >= 500000) {
                std::cout << "\nLimited to " << paired_reads.size() << " pairs for testing" << std::endl;
                break;
            }
        }
        
        std::cout << "\nLoaded " << paired_reads.size() << " paired-end reads" << std::endl;
        return paired_reads;
    }
    
    std::vector<PairedRead> load_interleaved_fastq(const std::string& fastq_path) {
        std::cout << "Loading interleaved paired-end reads" << std::endl;
        std::vector<PairedRead> paired_reads;
        
        std::ifstream file(fastq_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open FASTQ file");
        }
        
        std::string line;
        std::vector<std::string> current_reads;
        int line_count = 0;
        
        while (std::getline(file, line)) {
            line_count++;
            
            if (line_count % 4 == 2) {  // Sequence line
                if (!line.empty() && line.length() >= k) {
                    current_reads.push_back(line);
                    
                    if (current_reads.size() == 2) {
                        std::string read_id = "pair_" + std::to_string(paired_reads.size());
                        PairedRead pair(current_reads[0], current_reads[1], read_id);
                        paired_reads.push_back(pair);
                        current_reads.clear();
                    }
                }
            }
        }
        
        std::cout << "Loaded " << paired_reads.size() << " paired reads" << std::endl;
        return paired_reads;
    }
    
    std::string extract_read_id(const std::string& header) {
        if (header.empty()) {
            return "";
        }
        size_t space_pos = header.find(' ');
        if (space_pos != std::string::npos) {
            return header.substr(1, space_pos - 1);
        } else {
            return header.substr(1);
        }
    }
    
public:
    std::vector<OrganismProfile> profile_paired_end_community(const std::vector<PairedRead>& paired_reads) {
        std::cout << "\n=== STARTING PAIRED-END COMMUNITY PROFILING ===" << std::endl;
        std::cout << "Processing " << paired_reads.size() << " read pairs" << std::endl;
        std::cout << "Mode: " << (use_streaming_mode ? "STREAMING" : "DIRECT GPU") << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::vector<OrganismProfile> profiles;
        
        if (use_streaming_mode) {
            profiles = profile_with_streaming(paired_reads);
        } else {
            profiles = profile_with_direct_gpu(paired_reads);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Paired-end profiling completed in " << duration.count() << " ms" << std::endl;
        std::cout << "Mode used: " << (use_streaming_mode ? "Streaming" : "Direct GPU") << std::endl;
        
        return profiles;
    }
    
private:
    std::vector<OrganismProfile> profile_with_direct_gpu(const std::vector<PairedRead>& paired_reads) {
        std::cout << "Using direct GPU profiling..." << std::endl;
        
        // Placeholder implementation - would use direct GPU profiler
        // Similar to the gpu_community_profiler but with paired-end awareness
        
        std::vector<OrganismProfile> profiles;
        
        // For now, return a simple placeholder result
        if (!paired_reads.empty()) {
            OrganismProfile example;
            example.taxonomy_id = 1;
            example.name = "Example organism (direct GPU)";
            example.abundance = 0.5f;
            example.confidence_score = 0.8f;
            profiles.push_back(example);
        }
        
        std::cout << "Direct GPU profiling completed" << std::endl;
        return profiles;
    }
    
    std::vector<OrganismProfile> profile_with_streaming(const std::vector<PairedRead>& paired_reads) {
        std::cout << "Using streaming profiling..." << std::endl;
        
        // Extract all minimizers from paired reads
        std::vector<uint64_t> all_minimizers;
        extract_paired_minimizers(paired_reads, all_minimizers);
        
        std::cout << "Extracted " << all_minimizers.size() << " minimizers from paired reads" << std::endl;
        
        // Use streaming approach similar to streaming_gpu_profiler
        // but with paired-end awareness
        
        std::vector<OrganismProfile> profiles;
        
        // Placeholder implementation
        if (!paired_reads.empty()) {
            OrganismProfile example;
            example.taxonomy_id = 1;
            example.name = "Example organism (streaming)";
            example.abundance = 0.5f;
            example.confidence_score = 0.8f;
            example.paired_hits = paired_reads.size() / 10;  // Some paired hits
            profiles.push_back(example);
        }
        
        std::cout << "Streaming profiling completed" << std::endl;
        return profiles;
    }
    
    void extract_paired_minimizers(const std::vector<PairedRead>& paired_reads, 
                                  std::vector<uint64_t>& all_minimizers) {
        std::cout << "Extracting minimizers from paired-end reads..." << std::endl;
        
        int window_size = k - m + 1;
        
        for (const auto& pair : paired_reads) {
            // Extract from R1
            auto r1_minimizers = extract_read_minimizers(pair.read1);
            all_minimizers.insert(all_minimizers.end(), r1_minimizers.begin(), r1_minimizers.end());
            
            // Extract from R2 if present
            if (!pair.read2.empty()) {
                auto r2_minimizers = extract_read_minimizers(pair.read2);
                all_minimizers.insert(all_minimizers.end(), r2_minimizers.begin(), r2_minimizers.end());
            }
        }
        
        std::cout << "Extracted " << all_minimizers.size() << " total minimizers" << std::endl;
    }
    
    std::vector<uint64_t> extract_read_minimizers(const std::string& sequence) {
        std::vector<uint64_t> minimizers;
        
        if (sequence.length() < k) return minimizers;
        
        int window_size = k - m + 1;
        
        for (size_t i = 0; i <= sequence.length() - k; i += window_size) {
            std::string kmer_window = sequence.substr(i, k);
            
            uint64_t min_hash = UINT64_MAX;
            bool found_valid = false;
            
            for (int j = 0; j <= k - m; j++) {
                std::string minimizer_seq = kmer_window.substr(j, m);
                uint64_t hash = hash_minimizer(minimizer_seq);
                
                if (hash != UINT64_MAX && hash < min_hash) {
                    min_hash = hash;
                    found_valid = true;
                }
            }
            
            if (found_valid) {
                minimizers.push_back(min_hash);
            }
        }
        
        return minimizers;
    }
    
    uint64_t hash_minimizer(const std::string& seq) {
        uint64_t forward_hash = 0;
        uint64_t reverse_hash = 0;
        
        bool valid = true;
        for (size_t i = 0; i < seq.length(); i++) {
            int base = encode_base(seq[i]);
            if (base == -1) {
                valid = false;
                break;
            }
            forward_hash = (forward_hash << 2) | base;
            reverse_hash = (reverse_hash >> 2) | (((uint64_t)(3 ^ base)) << (2 * (seq.length() - 1)));
        }
        
        return valid ? std::min(forward_hash, reverse_hash) : UINT64_MAX;
    }
    
    int encode_base(char base) {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }
    
public:
    void print_results(const std::vector<OrganismProfile>& profiles) {
        std::cout << "\n=== ADAPTIVE PAIRED-END PROFILING RESULTS ===" << std::endl;
        std::cout << "Processing mode: " << (use_streaming_mode ? "STREAMING" : "DIRECT GPU") << std::endl;
        std::cout << "Database stats: " << db_stats.num_organisms << " organisms, " 
                  << db_stats.estimated_gpu_memory_mb << " MB estimated" << std::endl;
        
        if (profiles.empty()) {
            std::cout << "No organisms detected above significance thresholds" << std::endl;
            return;
        }
        
        std::cout << "\nTop organisms detected:" << std::endl;
        for (size_t i = 0; i < std::min(size_t(10), profiles.size()); i++) {
            const auto& profile = profiles[i];
            std::cout << (i+1) << ". " << profile.name 
                      << " - " << std::fixed << std::setprecision(2) << (profile.abundance * 100) << "%"
                      << " (confidence: " << profile.confidence_score << ")";
            
            if (profile.paired_hits > 0) {
                std::cout << " [" << profile.paired_hits << " paired hits]";
            }
            std::cout << std::endl;
        }
    }
    
    void write_results(const std::vector<OrganismProfile>& profiles, const std::string& output_prefix) {
        std::string output_file = output_prefix + "_adaptive_abundance.tsv";
        std::ofstream out(output_file);
        
        out << "taxonomy_id\torganism_name\ttaxonomy_path\trelative_abundance\t"
            << "coverage_breadth\tunique_minimizers\ttotal_hits\tpaired_hits\t"
            << "confidence_score\tprocessing_mode\n";
        
        std::string mode = use_streaming_mode ? "streaming" : "direct_gpu";
        
        for (const auto& profile : profiles) {
            out << profile.taxonomy_id << "\t"
                << profile.name << "\t"
                << profile.taxonomy_path << "\t"
                << std::scientific << profile.abundance << "\t"
                << std::fixed << std::setprecision(6) << profile.coverage_breadth << "\t"
                << profile.unique_minimizers << "\t"
                << profile.total_hits << "\t"
                << profile.paired_hits << "\t"
                << std::fixed << std::setprecision(4) << profile.confidence_score << "\t"
                << mode << "\n";
        }
        
        out.close();
        std::cout << "Results written to: " << output_file << std::endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <database> <R1.fastq> [R2.fastq] [output_prefix]" << std::endl;
        std::cerr << "\nAdaptive paired-end profiler with auto-streaming" << std::endl;
        std::cerr << "Automatically selects optimal processing strategy based on database size" << std::endl;
        std::cerr << "\nExamples:" << std::endl;
        std::cerr << "  " << argv[0] << " microbes.db sample_R1.fastq sample_R2.fastq results" << std::endl;
        std::cerr << "  " << argv[0] << " microbes.db interleaved.fastq results" << std::endl;
        return 1;
    }
    
    std::string database_path = argv[1];
    std::string r1_path = argv[2];
    std::string r2_path = "";
    std::string output_prefix = "adaptive_profile";
    
    // Parse arguments
    if (argc > 3) {
        std::string arg3 = argv[3];
        if (arg3.find(".fastq") != std::string::npos || arg3.find(".fq") != std::string::npos) {
            r2_path = arg3;
            if (argc > 4) output_prefix = argv[4];
        } else {
            output_prefix = arg3;
        }
    }
    
    try {
        AdaptivePairedEndProfiler profiler;
        
        // Load and analyze database - automatically determines processing strategy
        profiler.load_and_analyze_database(database_path);
        
        // Load paired-end reads
        auto paired_reads = profiler.load_paired_fastq(r1_path, r2_path);
        
        // Profile community with optimal strategy
        auto profiles = profiler.profile_paired_end_community(paired_reads);
        
        // Print and save results
        profiler.print_results(profiles);
        profiler.write_results(profiles, output_prefix);
        
        std::cout << "\n=== ADAPTIVE PROFILING COMPLETE ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}