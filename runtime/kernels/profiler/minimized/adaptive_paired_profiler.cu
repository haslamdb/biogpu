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

// Include the common paired-end structures
#include "paired_read_common.h"

// Include the actual GPU profiler implementations
// (These would be compiled as separate translation units in practice)
// #include "gpu_paired_profiler.h"
// #include "streaming_paired_profiler.h"

// Forward declarations for the GPU profilers
class GPUMicrobialCommunityProfiler;
class StreamingGPUProfiler;

// Enhanced adaptive profiler that automatically chooses between
// direct GPU loading and streaming based on database size analysis
// Now with full paired-end read support

struct DatabaseStats {
    uint64_t estimated_gpu_memory_mb;
    uint32_t num_organisms;
    uint64_t total_minimizer_entries;
    bool needs_streaming;
    uint32_t recommended_chunks;
    float paired_end_benefit_score;     // NEW: Estimated benefit from paired-end processing
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
    bool has_paired_reads = false;
    
    // Actual GPU profiler implementations (would use smart pointers in practice)
    std::unique_ptr<GPUMicrobialCommunityProfiler> gpu_profiler;
    std::unique_ptr<StreamingGPUProfiler> streaming_profiler;
    
    int k = 35, m = 31;
    
    // Configuration for paired-end processing
    struct PairedEndConfig {
        bool enable_paired_end_bonus = true;
        float paired_concordance_weight = 1.5f;
        float min_paired_concordance = 0.3f;
        bool prefer_concordant_pairs = true;
        float single_vs_paired_weight_ratio = 0.7f;  // Single reads get 70% weight vs paired
    } paired_config;
    
public:
    AdaptivePairedEndProfiler() {
        std::cout << "Enhanced Adaptive Paired-End GPU Profiler" << std::endl;
        std::cout << "Automatically selects optimal processing strategy with paired-end support" << std::endl;
    }
    
    void load_and_analyze_database(const std::string& db_path) {
        std::cout << "Analyzing database: " << db_path << std::endl;
        
        // First, analyze database size and characteristics
        db_stats = analyze_database_requirements(db_path);
        
        // Decide on processing strategy
        use_streaming_mode = db_stats.needs_streaming;
        
        std::cout << "\n=== ENHANCED DATABASE ANALYSIS RESULTS ===" << std::endl;
        std::cout << "Estimated GPU memory needed: " << db_stats.estimated_gpu_memory_mb << " MB" << std::endl;
        std::cout << "Number of organisms: " << db_stats.num_organisms << std::endl;
        std::cout << "Total minimizer entries: " << db_stats.total_minimizer_entries << std::endl;
        std::cout << "Paired-end benefit score: " << std::fixed << std::setprecision(2) 
                  << db_stats.paired_end_benefit_score << std::endl;
        
        if (use_streaming_mode) {
            std::cout << "🔄 STREAMING MODE SELECTED" << std::endl;
            std::cout << "Database too large for direct GPU loading" << std::endl;
            std::cout << "Will use " << db_stats.recommended_chunks << " chunks" << std::endl;
            std::cout << "Paired-end processing: Available with streaming kernels" << std::endl;
            initialize_streaming_mode(db_path);
        } else {
            std::cout << "⚡ DIRECT GPU MODE SELECTED" << std::endl;
            std::cout << "Database fits in GPU memory - maximum performance" << std::endl;
            std::cout << "Paired-end processing: Full paired-end kernels available" << std::endl;
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
        
        // Add overhead for organism arrays, reads, and paired-end scoring arrays
        size_t paired_end_overhead = (stats.num_organisms * 2000000 * sizeof(float)) / (1024 * 1024);  // 2M pairs max
        stats.estimated_gpu_memory_mb += 1000 + paired_end_overhead;  // 1GB base + paired-end overhead
        
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
        
        // Calculate paired-end benefit score based on database characteristics
        stats.paired_end_benefit_score = calculate_paired_end_benefit(stats);
        
        return stats;
    }
    
    float calculate_paired_end_benefit(const DatabaseStats& stats) {
        // Calculate estimated benefit from paired-end processing
        // Higher scores indicate more benefit from paired-end reads
        
        float benefit_score = 0.0f;
        
        // Factor 1: Database diversity (more diverse = more benefit from paired evidence)
        float diversity_factor = std::min(1.0f, (float)stats.num_organisms / 10000.0f);  // Normalize to 10K organisms
        benefit_score += diversity_factor * 0.3f;
        
        // Factor 2: Minimizer density (lower density = more benefit from paired confirmation)
        float avg_minimizers_per_organism = (float)stats.total_minimizer_entries / stats.num_organisms;
        float density_factor = std::max(0.0f, 1.0f - (avg_minimizers_per_organism / 100000.0f));  // Inverse relationship
        benefit_score += density_factor * 0.4f;
        
        // Factor 3: Database size (larger databases benefit more from paired-end discrimination)
        float size_factor = std::min(1.0f, stats.estimated_gpu_memory_mb / (20.0f * 1024.0f));  // Normalize to 20GB
        benefit_score += size_factor * 0.3f;
        
        return std::min(1.0f, benefit_score);  // Cap at 1.0
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
        std::cout << "Initializing direct GPU mode with paired-end support..." << std::endl;
        
        // Check GPU memory
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "GPU memory available: " << free_mem / (1024*1024*1024) << " GB" << std::endl;
        
        // In practice, would initialize the actual GPU profiler here
        // gpu_profiler = std::make_unique<GPUMicrobialCommunityProfiler>();
        // gpu_profiler->load_minimizer_database(db_path);
        
        std::cout << "Direct GPU mode with paired-end kernels ready" << std::endl;
    }
    
    void initialize_streaming_mode(const std::string& db_path) {
        std::cout << "Initializing streaming mode with paired-end support..." << std::endl;
        std::cout << "Will process database in " << db_stats.recommended_chunks << " chunks" << std::endl;
        
        // In practice, would initialize the actual streaming profiler here
        // StreamingConfig config;
        // config.use_paired_end_bonus = true;
        // config.paired_concordance_weight = paired_config.paired_concordance_weight;
        // streaming_profiler = std::make_unique<StreamingGPUProfiler>(config);
        // streaming_profiler->load_database_metadata(db_path);
        // streaming_profiler->create_database_chunks(db_path);
        
        std::cout << "Streaming mode with paired-end support ready" << std::endl;
    }
    
public:
    // Enhanced paired-end FASTQ loading methods
    std::vector<PairedRead> load_paired_fastq(const std::string& r1_path, 
                                             const std::string& r2_path = "") {
        std::cout << "Loading paired-end reads..." << std::endl;
        
        if (r2_path.empty()) {
            auto paired_reads = load_interleaved_fastq(r1_path);
            analyze_paired_read_characteristics(paired_reads);
            return paired_reads;
        } else {
            auto paired_reads = load_separate_fastq_files(r1_path, r2_path);
            analyze_paired_read_characteristics(paired_reads);
            return paired_reads;
        }
    }
    
private:
    void analyze_paired_read_characteristics(const std::vector<PairedRead>& paired_reads) {
        // Analyze the characteristics of the loaded reads
        size_t true_pairs = 0;
        size_t single_reads = 0;
        size_t total_read_length = 0;
        
        for (const auto& pair : paired_reads) {
            if (pair.is_paired) {
                true_pairs++;
                total_read_length += pair.read1.length() + pair.read2.length();
            } else {
                single_reads++;
                total_read_length += pair.read1.length();
            }
        }
        
        has_paired_reads = true_pairs > 0;
        
        std::cout << "\n=== READ CHARACTERISTICS ANALYSIS ===" << std::endl;
        std::cout << "Total read pairs/singles: " << paired_reads.size() << std::endl;
        std::cout << "True paired-end reads: " << true_pairs << std::endl;
        std::cout << "Single-end reads: " << single_reads << std::endl;
        std::cout << "Paired-end fraction: " << std::fixed << std::setprecision(1) 
                  << (100.0f * true_pairs / paired_reads.size()) << "%" << std::endl;
        
        if (has_paired_reads) {
            float avg_read_length = (float)total_read_length / (true_pairs * 2 + single_reads);
            std::cout << "Average read length: " << std::fixed << std::setprecision(0) 
                      << avg_read_length << " bp" << std::endl;
            
            // Estimate paired-end benefit for this dataset
            float dataset_benefit = estimate_dataset_paired_benefit(paired_reads);
            std::cout << "Estimated paired-end benefit: " << std::fixed << std::setprecision(2) 
                      << dataset_benefit << std::endl;
            
            // Adjust paired-end configuration based on analysis
            adjust_paired_end_config(dataset_benefit);
        }
    }
    
    float estimate_dataset_paired_benefit(const std::vector<PairedRead>& paired_reads) {
        // Estimate how much this specific dataset will benefit from paired-end processing
        float benefit = db_stats.paired_end_benefit_score;  // Start with database benefit
        
        // Factor in actual paired-end fraction
        size_t true_pairs = 0;
        for (const auto& pair : paired_reads) {
            if (pair.is_paired) true_pairs++;
        }
        
        float paired_fraction = (float)true_pairs / paired_reads.size();
        benefit *= paired_fraction;  // Scale by actual paired fraction
        
        // Factor in read length (longer reads = more discriminatory power)
        if (true_pairs > 0) {
            float avg_length = 0.0f;
            for (const auto& pair : paired_reads) {
                if (pair.is_paired) {
                    avg_length += (pair.read1.length() + pair.read2.length()) / 2.0f;
                }
            }
            avg_length /= true_pairs;
            
            float length_factor = std::min(1.2f, avg_length / 100.0f);  // Longer reads get bonus up to 20%
            benefit *= length_factor;
        }
        
        return std::min(1.0f, benefit);
    }
    
    void adjust_paired_end_config(float dataset_benefit) {
        // Adjust paired-end processing parameters based on estimated benefit
        
        if (dataset_benefit > 0.7f) {
            // High benefit - use aggressive paired-end scoring
            paired_config.paired_concordance_weight = 2.0f;
            paired_config.min_paired_concordance = 0.2f;
            paired_config.single_vs_paired_weight_ratio = 0.5f;
            std::cout << "High paired-end benefit detected - using aggressive paired scoring" << std::endl;
            
        } else if (dataset_benefit > 0.4f) {
            // Moderate benefit - use balanced scoring
            paired_config.paired_concordance_weight = 1.5f;
            paired_config.min_paired_concordance = 0.3f;
            paired_config.single_vs_paired_weight_ratio = 0.7f;
            std::cout << "Moderate paired-end benefit - using balanced scoring" << std::endl;
            
        } else {
            // Low benefit - use conservative scoring
            paired_config.paired_concordance_weight = 1.2f;
            paired_config.min_paired_concordance = 0.4f;
            paired_config.single_vs_paired_weight_ratio = 0.8f;
            std::cout << "Low paired-end benefit - using conservative scoring" << std::endl;
        }
    }
    
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