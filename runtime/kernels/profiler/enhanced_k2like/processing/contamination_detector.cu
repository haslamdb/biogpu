// processing/contamination_detector.cu
// Contamination detection for minimizers (human, adapter, low-complexity)
// Identifies and flags potentially problematic sequences

#include "contamination_detector.h"
#include "../gpu_kraken_types.h"
#include "../gpu/gpu_database_kernels.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <cmath>

// ===========================
// Device Functions
// ===========================

__device__ bool is_low_complexity_kmer(uint64_t kmer_hash, int k) {
    // Check for low complexity patterns
    // Simple entropy-based check on the hash bits
    
    int bit_changes = 0;
    uint64_t prev_bit = kmer_hash & 1;
    
    for (int i = 1; i < 64; i++) {
        uint64_t current_bit = (kmer_hash >> i) & 1;
        if (current_bit != prev_bit) {
            bit_changes++;
        }
        prev_bit = current_bit;
    }
    
    // Low complexity if very few bit transitions
    return bit_changes < 8;
}

__device__ bool check_homopolymer_run(uint64_t kmer_hash, int k) {
    // Check if k-mer represents a homopolymer or near-homopolymer
    // Simplified check based on hash pattern
    
    // Extract 2-bit encoded bases (assuming standard encoding)
    uint64_t mask = 0x3; // 2 bits
    uint64_t first_base = kmer_hash & mask;
    int same_base_count = 1;
    int max_run = 1;
    
    for (int i = 1; i < k && i < 32; i++) {
        uint64_t current_base = (kmer_hash >> (i * 2)) & mask;
        if (current_base == first_base) {
            same_base_count++;
            max_run = max(max_run, same_base_count);
        } else {
            same_base_count = 1;
            first_base = current_base;
        }
    }
    
    // Flag if >70% same base
    return max_run > (k * 0.7);
}

__device__ uint8_t calculate_contamination_score(
    bool is_human, 
    bool is_adapter, 
    bool is_low_complexity,
    bool is_homopolymer) {
    
    uint8_t score = 0;
    
    if (is_human) score += 4;        // High weight for human
    if (is_adapter) score += 3;      // High weight for adapter
    if (is_low_complexity) score += 2; // Medium weight
    if (is_homopolymer) score += 1;   // Low weight
    
    // Map to 3-bit score (0-7)
    return min(score, (uint8_t)7);
}

// ===========================
// CUDA Kernels
// ===========================

__global__ void mark_contamination_kernel(
    GPUMinimizerHit* minimizer_hits,
    const uint64_t* d_human_hashes,
    const uint64_t* d_adapter_hashes,
    int num_hits,
    int num_human_hashes,
    int num_adapter_hashes,
    int k_value,
    ContaminationRiskLevel risk_threshold) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    GPUMinimizerHit& hit = minimizer_hits[tid];
    uint64_t hash = hit.minimizer_hash;
    
    // Check human contamination (binary search in sorted array)
    bool is_human = false;
    if (num_human_hashes > 0) {
        int left = 0, right = num_human_hashes - 1;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (d_human_hashes[mid] == hash) {
                is_human = true;
                break;
            } else if (d_human_hashes[mid] < hash) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
    }
    
    // Check adapter contamination
    bool is_adapter = false;
    if (num_adapter_hashes > 0) {
        int left = 0, right = num_adapter_hashes - 1;
        while (left <= right) {
            int mid = (left + right) / 2;
            if (d_adapter_hashes[mid] == hash) {
                is_adapter = true;
                break;
            } else if (d_adapter_hashes[mid] < hash) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
    }
    
    // Check low complexity
    bool is_low_complexity = is_low_complexity_kmer(hash, k_value);
    bool is_homopolymer = check_homopolymer_run(hash, k_value);
    
    // Calculate contamination score
    uint8_t contamination_score = calculate_contamination_score(
        is_human, is_adapter, is_low_complexity, is_homopolymer
    );
    
    // Set contamination flags based on risk level
    if (contamination_score >= (uint8_t)risk_threshold) {
        // Set contamination risk bit (bit 7)
        hit.feature_flags |= (1 << 7);
        
        // Store contamination type in bits 13-15
        uint16_t contam_type = 0;
        if (is_human) contam_type |= 1;
        if (is_adapter) contam_type |= 2;
        if (is_low_complexity || is_homopolymer) contam_type |= 4;
        
        hit.feature_flags |= (contam_type << 13);
    }
}

__global__ void batch_contamination_check_kernel(
    const uint64_t* minimizer_hashes,
    uint8_t* contamination_flags,
    const uint64_t* d_contamination_db,
    int num_minimizers,
    int db_size) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_minimizers) return;
    
    uint64_t hash = minimizer_hashes[tid];
    
    // Binary search in contamination database
    int left = 0, right = db_size - 1;
    bool found = false;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        if (d_contamination_db[mid] == hash) {
            found = true;
            break;
        } else if (d_contamination_db[mid] < hash) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    contamination_flags[tid] = found ? 1 : 0;
}

// ===========================
// ContaminationDetector Implementation
// ===========================

ContaminationDetector::ContaminationDetector(const ContaminationConfig& config)
    : config_(config), initialized_(false),
      d_human_hashes_(nullptr), d_adapter_hashes_(nullptr),
      d_combined_contamination_db_(nullptr),
      num_human_hashes_(0), num_adapter_hashes_(0),
      total_contamination_hashes_(0) {
    
    // Pre-allocate GPU memory for contamination databases
    size_t max_db_size = config_.max_contamination_patterns;
    cudaMalloc(&d_human_hashes_, max_db_size * sizeof(uint64_t));
    cudaMalloc(&d_adapter_hashes_, max_db_size * sizeof(uint64_t));
    cudaMalloc(&d_combined_contamination_db_, max_db_size * 2 * sizeof(uint64_t));
    
    // Initialize with built-in patterns
    initialize_builtin_patterns();
}

ContaminationDetector::~ContaminationDetector() {
    if (d_human_hashes_) cudaFree(d_human_hashes_);
    if (d_adapter_hashes_) cudaFree(d_adapter_hashes_);
    if (d_combined_contamination_db_) cudaFree(d_combined_contamination_db_);
}

bool ContaminationDetector::load_contamination_database(const std::string& database_path) {
    std::cout << "Loading contamination database from: " << database_path << std::endl;
    
    std::ifstream db_file(database_path);
    if (!db_file.is_open()) {
        std::cerr << "Cannot open contamination database: " << database_path << std::endl;
        return false;
    }
    
    // Clear existing patterns
    human_minimizer_hashes_.clear();
    adapter_minimizer_hashes_.clear();
    
    std::string line;
    int patterns_loaded = 0;
    
    while (std::getline(db_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string type;
        uint64_t hash;
        
        if (iss >> type >> std::hex >> hash) {
            if (type == "HUMAN") {
                human_minimizer_hashes_.insert(hash);
            } else if (type == "ADAPTER") {
                adapter_minimizer_hashes_.insert(hash);
            } else if (type == "LOWCOMPLEXITY") {
                low_complexity_patterns_.insert(hash);
            }
            patterns_loaded++;
        }
    }
    
    db_file.close();
    
    std::cout << "Loaded " << patterns_loaded << " contamination patterns" << std::endl;
    std::cout << "  Human: " << human_minimizer_hashes_.size() << std::endl;
    std::cout << "  Adapter: " << adapter_minimizer_hashes_.size() << std::endl;
    std::cout << "  Low complexity: " << low_complexity_patterns_.size() << std::endl;
    
    // Update GPU memory
    return update_gpu_databases();
}

ContaminationRiskLevel ContaminationDetector::check_minimizer_contamination(
    uint64_t minimizer_hash) const {
    
    // Quick host-side check
    bool is_human = human_minimizer_hashes_.find(minimizer_hash) != human_minimizer_hashes_.end();
    bool is_adapter = adapter_minimizer_hashes_.find(minimizer_hash) != adapter_minimizer_hashes_.end();
    bool is_low_complexity = low_complexity_patterns_.find(minimizer_hash) != low_complexity_patterns_.end();
    
    if (is_human) return ContaminationRiskLevel::HIGH_RISK;
    if (is_adapter) return ContaminationRiskLevel::MEDIUM_RISK;
    if (is_low_complexity) return ContaminationRiskLevel::LOW_RISK;
    
    return ContaminationRiskLevel::NO_RISK;
}

bool ContaminationDetector::mark_contamination_in_batch(
    GPUMinimizerHit* d_minimizer_hits,
    size_t num_hits,
    ContaminationRiskLevel risk_threshold) {
    
    if (!initialized_ || num_hits == 0) return false;
    
    int block_size = 256;
    int grid_size = (num_hits + block_size - 1) / block_size;
    
    mark_contamination_kernel<<<grid_size, block_size>>>(
        d_minimizer_hits,
        d_human_hashes_,
        d_adapter_hashes_,
        num_hits,
        num_human_hashes_,
        num_adapter_hashes_,
        config_.k_value,
        risk_threshold
    );
    
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        std::cerr << "Contamination marking kernel failed: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    cudaDeviceSynchronize();
    
    // Update statistics by analyzing results
    std::vector<GPUMinimizerHit> h_hits(num_hits);
    cudaMemcpy(h_hits.data(), d_minimizer_hits, num_hits * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
    
    for (const auto& hit : h_hits) {
        if (hit.feature_flags & (1 << 7)) {
            // Extract contamination type from bits 13-15
            uint16_t contam_type = (hit.feature_flags >> 13) & 0x7;
            if (contam_type & 1) stats_.human_contamination_count++;
            if (contam_type & 2) stats_.adapter_contamination_count++;
            if (contam_type & 4) stats_.low_complexity_count++;
        }
    }
    
    stats_.batches_processed++;
    stats_.total_minimizers_checked += num_hits;
    
    return true;
}

bool ContaminationDetector::quick_contamination_check(
    const uint64_t* d_minimizer_hashes,
    uint8_t* d_contamination_flags,
    size_t num_minimizers) {
    
    if (!initialized_ || num_minimizers == 0) return false;
    
    int block_size = 256;
    int grid_size = (num_minimizers + block_size - 1) / block_size;
    
    batch_contamination_check_kernel<<<grid_size, block_size>>>(
        d_minimizer_hashes,
        d_contamination_flags,
        d_combined_contamination_db_,
        num_minimizers,
        total_contamination_hashes_
    );
    
    cudaDeviceSynchronize();
    return true;
}

void ContaminationDetector::add_custom_contamination_pattern(
    uint64_t minimizer_hash,
    ContaminationType type) {
    
    switch (type) {
        case ContaminationType::HUMAN:
            human_minimizer_hashes_.insert(minimizer_hash);
            break;
        case ContaminationType::ADAPTER:
            adapter_minimizer_hashes_.insert(minimizer_hash);
            break;
        case ContaminationType::LOW_COMPLEXITY:
            low_complexity_patterns_.insert(minimizer_hash);
            break;
    }
}

bool ContaminationDetector::export_contamination_report(const std::string& output_file) const {
    std::ofstream out(output_file);
    if (!out.is_open()) return false;
    
    out << "# Contamination Detection Report\n";
    out << "Total minimizers checked: " << stats_.total_minimizers_checked << "\n";
    out << "Human contamination detected: " << stats_.human_contamination_count << "\n";
    out << "Adapter contamination detected: " << stats_.adapter_contamination_count << "\n";
    out << "Low complexity detected: " << stats_.low_complexity_count << "\n";
    out << "\n";
    
    out << "# Contamination Rate\n";
    if (stats_.total_minimizers_checked > 0) {
        double human_rate = (double)stats_.human_contamination_count / stats_.total_minimizers_checked * 100;
        double adapter_rate = (double)stats_.adapter_contamination_count / stats_.total_minimizers_checked * 100;
        double low_complexity_rate = (double)stats_.low_complexity_count / stats_.total_minimizers_checked * 100;
        
        out << "Human contamination rate: " << std::fixed << std::setprecision(2) << human_rate << "%\n";
        out << "Adapter contamination rate: " << adapter_rate << "%\n";
        out << "Low complexity rate: " << low_complexity_rate << "%\n";
    }
    
    out.close();
    return true;
}

// Private methods

void ContaminationDetector::initialize_builtin_patterns() {
    // Common Illumina adapter sequences (would compute all k-mers in production)
    // TruSeq Universal Adapter
    add_adapter_pattern("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    // TruSeq Index Adapter
    add_adapter_pattern("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
    
    // Nextera Transposase Sequence
    add_adapter_pattern("CTGTCTCTTATACACATCT");
    
    // Common low-complexity patterns
    // Poly-A/T runs
    for (int i = 0; i < 10; i++) {
        uint64_t polyA = 0x0000000000000000; // All A's
        uint64_t polyT = 0xFFFFFFFFFFFFFFFF; // All T's
        low_complexity_patterns_.insert(polyA + i);
        low_complexity_patterns_.insert(polyT - i);
    }
    
    std::cout << "Initialized built-in contamination patterns" << std::endl;
    
    // Update GPU databases with built-in patterns
    update_gpu_databases();
}

bool ContaminationDetector::update_gpu_databases() {
    // Convert sets to sorted vectors for GPU binary search
    std::vector<uint64_t> human_vec(human_minimizer_hashes_.begin(), human_minimizer_hashes_.end());
    std::vector<uint64_t> adapter_vec(adapter_minimizer_hashes_.begin(), adapter_minimizer_hashes_.end());
    
    std::sort(human_vec.begin(), human_vec.end());
    std::sort(adapter_vec.begin(), adapter_vec.end());
    
    num_human_hashes_ = human_vec.size();
    num_adapter_hashes_ = adapter_vec.size();
    
    // Copy to GPU
    if (num_human_hashes_ > 0) {
        cudaMemcpy(d_human_hashes_, human_vec.data(), 
                   num_human_hashes_ * sizeof(uint64_t), cudaMemcpyHostToDevice);
    }
    
    if (num_adapter_hashes_ > 0) {
        cudaMemcpy(d_adapter_hashes_, adapter_vec.data(), 
                   num_adapter_hashes_ * sizeof(uint64_t), cudaMemcpyHostToDevice);
    }
    
    // Create combined database for quick checks
    std::vector<uint64_t> combined;
    combined.insert(combined.end(), human_vec.begin(), human_vec.end());
    combined.insert(combined.end(), adapter_vec.begin(), adapter_vec.end());
    combined.insert(combined.end(), low_complexity_patterns_.begin(), low_complexity_patterns_.end());
    
    std::sort(combined.begin(), combined.end());
    combined.erase(std::unique(combined.begin(), combined.end()), combined.end());
    
    total_contamination_hashes_ = combined.size();
    
    if (total_contamination_hashes_ > 0) {
        cudaMemcpy(d_combined_contamination_db_, combined.data(),
                   total_contamination_hashes_ * sizeof(uint64_t), cudaMemcpyHostToDevice);
    }
    
    initialized_ = true;
    return true;
}

void ContaminationDetector::add_adapter_pattern(const std::string& adapter_sequence) {
    // In production, this would generate all k-mers from the adapter
    // For now, using a simple hash of the sequence
    std::hash<std::string> hasher;
    uint64_t adapter_hash = hasher(adapter_sequence);
    adapter_minimizer_hashes_.insert(adapter_hash);
}

uint64_t ContaminationDetector::compute_minimizer_from_sequence(
    const std::string& sequence, size_t pos) const {
    // Simplified - in production would use actual minimizer extraction
    if (pos + config_.k_value > sequence.length()) return 0;
    
    std::string kmer = sequence.substr(pos, config_.k_value);
    std::hash<std::string> hasher;
    return hasher(kmer);
}

bool ContaminationDetector::load_human_minimizers(const std::string& human_minimizer_file) {
    std::ifstream file(human_minimizer_file);
    if (!file.is_open()) {
        std::cerr << "Cannot open human minimizer file: " << human_minimizer_file << std::endl;
        return false;
    }
    
    std::string line;
    size_t count = 0;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        uint64_t hash;
        std::istringstream iss(line);
        if (iss >> std::hex >> hash) {
            human_minimizer_hashes_.insert(hash);
            count++;
        }
    }
    
    file.close();
    std::cout << "Loaded " << count << " human minimizers" << std::endl;
    
    return update_gpu_databases();
}

bool ContaminationDetector::load_adapter_sequences(const std::string& adapter_file) {
    std::ifstream file(adapter_file);
    if (!file.is_open()) {
        std::cerr << "Cannot open adapter file: " << adapter_file << std::endl;
        return false;
    }
    
    std::string line;
    std::string current_adapter;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // New adapter sequence
            current_adapter.clear();
        } else {
            // Sequence line
            current_adapter += line;
            
            // Generate k-mers from adapter
            if (current_adapter.length() >= config_.k_value) {
                for (size_t i = 0; i <= current_adapter.length() - config_.k_value; i++) {
                    uint64_t hash = compute_minimizer_from_sequence(current_adapter, i);
                    adapter_minimizer_hashes_.insert(hash);
                }
            }
        }
    }
    
    file.close();
    std::cout << "Loaded " << adapter_minimizer_hashes_.size() << " adapter k-mers" << std::endl;
    
    return update_gpu_databases();
}

void ContaminationDetector::clear_contamination_patterns() {
    human_minimizer_hashes_.clear();
    adapter_minimizer_hashes_.clear();
    low_complexity_patterns_.clear();
    
    num_human_hashes_ = 0;
    num_adapter_hashes_ = 0;
    total_contamination_hashes_ = 0;
    
    initialized_ = false;
}

// ===========================
// Utility Functions
// ===========================

namespace ContaminationDetectionUtils {
    
    bool create_contamination_database(
        const std::string& human_reference_path,
        const std::string& adapter_sequences_path,
        const std::string& output_database_path,
        int k_value) {
        
        std::cout << "Creating contamination database..." << std::endl;
        std::ofstream out(output_database_path);
        if (!out.is_open()) return false;
        
        out << "# Contamination Database\n";
        out << "# Format: TYPE HASH\n";
        out << "# Generated with k=" << k_value << "\n\n";
        
        // Process human reference (would extract all k-mers)
        // This is a placeholder - real implementation would process FASTA
        
        // Process adapter sequences
        std::ifstream adapters(adapter_sequences_path);
        if (adapters.is_open()) {
            std::string line;
            while (std::getline(adapters, line)) {
                if (line.empty() || line[0] == '>') continue;
                
                // Generate k-mers from adapter
                for (size_t i = 0; i <= line.length() - k_value; i++) {
                    std::string kmer = line.substr(i, k_value);
                    std::hash<std::string> hasher;
                    uint64_t hash = hasher(kmer);
                    out << "ADAPTER " << std::hex << hash << "\n";
                }
            }
            adapters.close();
        }
        
        out.close();
        std::cout << "Contamination database created: " << output_database_path << std::endl;
        return true;
    }
    
    std::vector<std::string> get_common_adapter_sequences() {
        std::vector<std::string> adapters;
        
        // TruSeq Universal Adapter
        adapters.push_back("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        
        // TruSeq Index Adapter
        adapters.push_back("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
        
        // Nextera Transposase
        adapters.push_back("CTGTCTCTTATACACATCT");
        
        // TruSeq3 Universal Adapter
        adapters.push_back("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
        
        // Illumina Small RNA 3' Adapter
        adapters.push_back("TGGAATTCTCGGGTGCCAAGG");
        
        // Illumina Small RNA 5' Adapter
        adapters.push_back("GTTCAGAGTTCTACAGTCCGACGATC");
        
        return adapters;
    }
    
    bool load_human_kmer_database(
        const std::string& database_path,
        std::unordered_set<uint64_t>& human_kmers) {
        
        std::ifstream db_file(database_path, std::ios::binary);
        if (!db_file.is_open()) {
            std::cerr << "Cannot open human k-mer database: " << database_path << std::endl;
            return false;
        }
        
        // Read header (if any)
        uint64_t num_kmers;
        db_file.read(reinterpret_cast<char*>(&num_kmers), sizeof(num_kmers));
        
        // Read k-mer hashes
        human_kmers.clear();
        for (uint64_t i = 0; i < num_kmers; i++) {
            uint64_t hash;
            db_file.read(reinterpret_cast<char*>(&hash), sizeof(hash));
            human_kmers.insert(hash);
        }
        
        db_file.close();
        std::cout << "Loaded " << human_kmers.size() << " human k-mers from database" << std::endl;
        
        return true;
    }
    
    void analyze_contamination_patterns(
        const std::vector<GPUMinimizerHit>& minimizer_hits,
        const ContaminationDetector& detector) {
        
        std::unordered_map<ContaminationType, int> type_counts;
        int total_contaminated = 0;
        
        for (const auto& hit : minimizer_hits) {
            if (hit.feature_flags & (1 << 7)) {
                total_contaminated++;
                
                // Extract contamination type from bits 13-15
                uint16_t contam_type = (hit.feature_flags >> 13) & 0x7;
                if (contam_type & 1) type_counts[ContaminationType::HUMAN]++;
                if (contam_type & 2) type_counts[ContaminationType::ADAPTER]++;
                if (contam_type & 4) type_counts[ContaminationType::LOW_COMPLEXITY]++;
            }
        }
        
        std::cout << "\nContamination Analysis:" << std::endl;
        std::cout << "Total contaminated minimizers: " << total_contaminated 
                  << " / " << minimizer_hits.size() << std::endl;
        
        for (const auto& [type, count] : type_counts) {
            std::string type_name;
            switch (type) {
                case ContaminationType::HUMAN: type_name = "Human"; break;
                case ContaminationType::ADAPTER: type_name = "Adapter"; break;
                case ContaminationType::LOW_COMPLEXITY: type_name = "Low Complexity"; break;
            }
            std::cout << "  " << type_name << ": " << count << std::endl;
        }
    }
    
    bool filter_contaminated_minimizers(
        GPUMinimizerHit* d_minimizer_hits,
        size_t& num_hits,
        ContaminationRiskLevel threshold) {
        
        // This would remove contaminated minimizers above threshold
        // Implementation would use thrust::remove_if
        return true;
    }
}