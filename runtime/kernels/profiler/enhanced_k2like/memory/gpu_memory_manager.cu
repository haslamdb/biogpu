// memory/gpu_memory_manager.cu
// GPU memory management implementation
// Handles allocation, auto-scaling, and memory optimization

#include "gpu_memory_manager.h"
#include "../gpu_kraken_types.h"
#include "../gpu/gpu_database_kernels.h"
#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

// ===========================
// GPU Memory Pool Implementation
// ===========================

GPUMemoryPool::GPUMemoryPool(size_t pool_size_mb) 
    : pool_base_(nullptr), pool_size_(pool_size_mb * 1024 * 1024), 
      current_offset_(0), initialized_(false) {
    
    cudaError_t error = cudaMalloc(&pool_base_, pool_size_);
    if (error == cudaSuccess) {
        initialized_ = true;
        std::cout << "GPU memory pool created: " << pool_size_mb << " MB" << std::endl;
    } else {
        std::cerr << "Failed to create GPU memory pool: " << cudaGetErrorString(error) << std::endl;
    }
}

GPUMemoryPool::~GPUMemoryPool() {
    if (pool_base_) {
        cudaFree(pool_base_);
        pool_base_ = nullptr;
    }
}

void* GPUMemoryPool::allocate(size_t size, size_t alignment) {
    if (!initialized_ || size == 0) return nullptr;
    
    // Align the current offset
    size_t aligned_offset = (current_offset_ + alignment - 1) & ~(alignment - 1);
    
    // Check if we have enough space
    if (aligned_offset + size > pool_size_) {
        return nullptr;  // Out of memory
    }
    
    void* ptr = static_cast<char*>(pool_base_) + aligned_offset;
    current_offset_ = aligned_offset + size;
    
    return ptr;
}

void GPUMemoryPool::deallocate(void* ptr) {
    // Simple pool doesn't support individual deallocation
    // Could be extended with a free list
}

void GPUMemoryPool::reset() {
    current_offset_ = 0;
}

size_t GPUMemoryPool::get_available_space() const {
    return pool_size_ - current_offset_;
}

// ===========================
// GPU Memory Manager Implementation
// ===========================

GPUMemoryManager::GPUMemoryManager(const MemoryConfig& config)
    : config_(config), initialized_(false), allocations_active_(false), total_allocated_(0),
      d_sequence_data_(nullptr), d_genome_info_(nullptr), d_minimizer_hits_(nullptr),
      d_lca_candidates_(nullptr), d_minimizer_counts_(nullptr) {
    
    // Initialize statistics
    memset(&stats_, 0, sizeof(stats_));
}

GPUMemoryManager::~GPUMemoryManager() {
    free_all_allocations();
}

bool GPUMemoryManager::initialize() {
    if (initialized_) return true;
    
    std::cout << "Initializing GPU Memory Manager..." << std::endl;
    
    // Query GPU memory information
    if (!query_gpu_memory_info()) {
        std::cerr << "Failed to query GPU memory information" << std::endl;
        return false;
    }
    
    // Validate GPU context
    if (!validate_gpu_context()) {
        std::cerr << "Invalid GPU context" << std::endl;
        return false;
    }
    
    // Create memory pool if enabled
    if (config_.enable_memory_pooling) {
        size_t pool_size_mb = (stats_.available_memory * config_.max_memory_fraction / 100) / (1024 * 1024);
        pool_size_mb = std::min(pool_size_mb, size_t(8192));  // Cap at 8GB
        
        memory_pool_ = std::make_unique<GPUMemoryPool>(pool_size_mb);
        if (!memory_pool_) {
            std::cerr << "Failed to create GPU memory pool" << std::endl;
            return false;
        }
    }
    
    // Calculate optimal configuration if auto-scaling is enabled
    if (config_.auto_scale_enabled) {
        if (!calculate_optimal_batch_sizes()) {
            std::cerr << "Failed to calculate optimal batch sizes" << std::endl;
            return false;
        }
    }
    
    initialized_ = true;
    std::cout << "✓ GPU Memory Manager initialized successfully" << std::endl;
    print_memory_usage();
    
    return true;
}

bool GPUMemoryManager::configure_auto_scaling(bool enable, size_t memory_fraction) {
    config_.auto_scale_enabled = enable;
    config_.max_memory_fraction = std::min(memory_fraction, size_t(95));  // Cap at 95%
    
    if (enable && initialized_) {
        return calculate_optimal_batch_sizes();
    }
    
    return true;
}

bool GPUMemoryManager::set_minimizer_capacity(int capacity) {
    if (capacity <= 0 || capacity > 50000000) {  // Cap at 50M
        std::cerr << "Invalid minimizer capacity: " << capacity << std::endl;
        return false;
    }
    
    config_.minimizer_capacity = capacity;
    std::cout << "Set minimizer capacity to " << capacity << std::endl;
    
    return true;
}

bool GPUMemoryManager::set_batch_size(int batch_size) {
    if (batch_size <= 0 || batch_size > 100) {
        std::cerr << "Invalid batch size: " << batch_size << std::endl;
        return false;
    }
    
    config_.sequence_batch_size = batch_size;
    std::cout << "Set sequence batch size to " << batch_size << std::endl;
    
    return true;
}

bool GPUMemoryManager::allocate_sequence_memory(size_t max_sequences, size_t max_total_length) {
    if (allocations_active_) {
        std::cerr << "Memory already allocated. Free existing allocations first." << std::endl;
        return false;
    }
    
    // Validate request
    size_t required_memory = max_total_length + max_sequences * sizeof(GPUGenomeInfo);
    if (!validate_allocation_request(required_memory)) {
        return false;
    }
    
    std::cout << "Allocating sequence memory:" << std::endl;
    std::cout << "  Max sequences: " << max_sequences << std::endl;
    std::cout << "  Total length: " << (max_total_length / 1024 / 1024) << " MB" << std::endl;
    
    // Allocate sequence data buffer
    cudaError_t error = cudaMalloc(&d_sequence_data_, max_total_length);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate sequence data: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    // Allocate genome info buffer
    error = cudaMalloc(&d_genome_info_, max_sequences * sizeof(GPUGenomeInfo));
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate genome info: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_sequence_data_);
        d_sequence_data_ = nullptr;
        return false;
    }
    
    // Initialize memory
    launch_memory_initialization_kernel(
        d_sequence_data_, max_total_length,
        d_genome_info_, max_sequences,
        nullptr, 0
    );
    
    stats_.current_sequence_memory = max_total_length + max_sequences * sizeof(GPUGenomeInfo);
    total_allocated_ += stats_.current_sequence_memory;
    
    std::cout << "✓ Sequence memory allocated successfully" << std::endl;
    return true;
}

bool GPUMemoryManager::allocate_minimizer_memory(size_t max_minimizers) {
    size_t minimizer_memory = max_minimizers * sizeof(GPUMinimizerHit);
    size_t count_memory = config_.sequence_batch_size * sizeof(uint32_t);
    size_t total_memory = minimizer_memory + count_memory;
    
    if (!validate_allocation_request(total_memory)) {
        return false;
    }
    
    std::cout << "Allocating minimizer memory for " << max_minimizers << " minimizers" << std::endl;
    
    // Allocate minimizer hits buffer
    cudaError_t error = cudaMalloc(&d_minimizer_hits_, minimizer_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate minimizer hits: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    // Allocate count buffer
    error = cudaMalloc(&d_minimizer_counts_, count_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate minimizer counts: " << cudaGetErrorString(error) << std::endl;
        cudaFree(d_minimizer_hits_);
        d_minimizer_hits_ = nullptr;
        return false;
    }
    
    // Initialize memory
    launch_memory_initialization_kernel(
        nullptr, 0,
        nullptr, 0,
        d_minimizer_hits_, max_minimizers
    );
    
    stats_.current_minimizer_memory = total_memory;
    total_allocated_ += stats_.current_minimizer_memory;
    
    std::cout << "✓ Minimizer memory allocated successfully" << std::endl;
    return true;
}

bool GPUMemoryManager::allocate_metadata_memory(size_t max_genomes) {
    // This would allocate additional metadata buffers
    stats_.current_metadata_memory = max_genomes * sizeof(uint32_t) * 4;  // Example
    return true;
}

bool GPUMemoryManager::allocate_results_memory(size_t max_candidates) {
    size_t results_memory = max_candidates * sizeof(LCACandidate);
    
    if (!validate_allocation_request(results_memory)) {
        return false;
    }
    
    std::cout << "Allocating results memory for " << max_candidates << " candidates" << std::endl;
    
    cudaError_t error = cudaMalloc(&d_lca_candidates_, results_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to allocate LCA candidates: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    // Initialize results memory
    cudaMemset(d_lca_candidates_, 0, results_memory);
    
    total_allocated_ += results_memory;
    
    std::cout << "✓ Results memory allocated successfully" << std::endl;
    return true;
}

void GPUMemoryManager::free_all_allocations() {
    std::cout << "Freeing all GPU memory allocations..." << std::endl;
    
    if (d_sequence_data_) {
        cudaFree(d_sequence_data_);
        d_sequence_data_ = nullptr;
    }
    
    if (d_genome_info_) {
        cudaFree(d_genome_info_);
        d_genome_info_ = nullptr;
    }
    
    if (d_minimizer_hits_) {
        cudaFree(d_minimizer_hits_);
        d_minimizer_hits_ = nullptr;
    }
    
    if (d_lca_candidates_) {
        cudaFree(d_lca_candidates_);
        d_lca_candidates_ = nullptr;
    }
    
    if (d_minimizer_counts_) {
        cudaFree(d_minimizer_counts_);
        d_minimizer_counts_ = nullptr;
    }
    
    // Reset memory pool if it exists
    if (memory_pool_) {
        memory_pool_->reset();
    }
    
    // Reset statistics
    stats_.current_sequence_memory = 0;
    stats_.current_minimizer_memory = 0;
    stats_.current_metadata_memory = 0;
    total_allocated_ = 0;
    allocations_active_ = false;
    
    std::cout << "✓ All GPU memory freed" << std::endl;
}

bool GPUMemoryManager::validate_memory_integrity() {
    if (!allocations_active_ || !d_sequence_data_ || !d_genome_info_) {
        return true;  // No memory to validate
    }
    
    std::cout << "Validating GPU memory integrity..." << std::endl;
    
    bool result = launch_memory_validation_kernel(
        d_sequence_data_, d_genome_info_, config_.sequence_batch_size
    );
    
    if (result) {
        std::cout << "✓ Memory integrity validation passed" << std::endl;
    } else {
        std::cerr << "⚠ Memory integrity validation failed" << std::endl;
    }
    
    return result;
}

// Accessor methods
char* GPUMemoryManager::get_sequence_buffer() const {
    return d_sequence_data_;
}

GPUGenomeInfo* GPUMemoryManager::get_genome_info_buffer() const {
    return d_genome_info_;
}

GPUMinimizerHit* GPUMemoryManager::get_minimizer_buffer() const {
    return d_minimizer_hits_;
}

LCACandidate* GPUMemoryManager::get_candidate_buffer() const {
    return d_lca_candidates_;
}

uint32_t* GPUMemoryManager::get_count_buffer() const {
    return d_minimizer_counts_;
}

MemoryStats GPUMemoryManager::get_memory_statistics() const {
    // Update current statistics
    const_cast<GPUMemoryManager*>(this)->update_memory_statistics();
    return stats_;
}

void GPUMemoryManager::print_memory_usage() const {
    std::cout << "\n=== GPU MEMORY USAGE ===" << std::endl;
    std::cout << "Total GPU memory: " << (stats_.total_gpu_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "Available memory: " << (stats_.available_memory / 1024 / 1024) << " MB" << std::endl;
    std::cout << "Currently allocated: " << (total_allocated_ / 1024 / 1024) << " MB" << std::endl;
    
    if (stats_.current_sequence_memory > 0) {
        std::cout << "  Sequence data: " << (stats_.current_sequence_memory / 1024 / 1024) << " MB" << std::endl;
    }
    if (stats_.current_minimizer_memory > 0) {
        std::cout << "  Minimizer data: " << (stats_.current_minimizer_memory / 1024 / 1024) << " MB" << std::endl;
    }
    if (stats_.current_metadata_memory > 0) {
        std::cout << "  Metadata: " << (stats_.current_metadata_memory / 1024 / 1024) << " MB" << std::endl;
    }
    
    double usage_percent = total_allocated_ > 0 ? 
        (double)total_allocated_ / stats_.total_gpu_memory * 100.0 : 0.0;
    std::cout << "Memory usage: " << std::fixed << std::setprecision(1) 
              << usage_percent << "%" << std::endl;
    
    if (memory_pool_) {
        std::cout << "Memory pool available: " 
                  << (memory_pool_->get_available_space() / 1024 / 1024) << " MB" << std::endl;
    }
}

bool GPUMemoryManager::check_memory_pressure() const {
    double usage_ratio = (double)total_allocated_ / stats_.available_memory;
    return usage_ratio > 0.85;  // Consider high pressure at 85%
}

size_t GPUMemoryManager::estimate_memory_requirements(int num_sequences, int avg_sequence_length) const {
    size_t sequence_memory = num_sequences * avg_sequence_length;
    size_t genome_info_memory = num_sequences * sizeof(GPUGenomeInfo);
    size_t minimizer_memory = config_.minimizer_capacity * sizeof(GPUMinimizerHit);
    size_t overhead = (sequence_memory + genome_info_memory + minimizer_memory) * 0.1;  // 10% overhead
    
    return sequence_memory + genome_info_memory + minimizer_memory + overhead;
}

bool GPUMemoryManager::scale_for_workload(size_t estimated_minimizers, size_t estimated_sequences) {
    if (!config_.auto_scale_enabled) return true;
    
    std::cout << "Auto-scaling for workload:" << std::endl;
    std::cout << "  Estimated minimizers: " << estimated_minimizers << std::endl;
    std::cout << "  Estimated sequences: " << estimated_sequences << std::endl;
    
    // Calculate required capacity
    size_t safe_minimizer_capacity = calculate_safe_minimizer_capacity();
    size_t required_capacity = std::max(estimated_minimizers, safe_minimizer_capacity);
    
    if (required_capacity != config_.minimizer_capacity) {
        std::cout << "Adjusting minimizer capacity: " << config_.minimizer_capacity 
                  << " → " << required_capacity << std::endl;
        config_.minimizer_capacity = required_capacity;
    }
    
    // Adjust batch size if needed
    int optimal_batch_size = calculate_optimal_batch_size();
    if (optimal_batch_size != config_.sequence_batch_size) {
        std::cout << "Adjusting batch size: " << config_.sequence_batch_size 
                  << " → " << optimal_batch_size << std::endl;
        config_.sequence_batch_size = optimal_batch_size;
    }
    
    return true;
}

void GPUMemoryManager::suggest_optimal_configuration() const {
    std::cout << "\n=== MEMORY OPTIMIZATION SUGGESTIONS ===" << std::endl;
    
    size_t total_gb = stats_.total_gpu_memory / (1024 * 1024 * 1024);
    
    if (total_gb >= 24) {
        std::cout << "High-memory GPU detected (" << total_gb << " GB)" << std::endl;
        std::cout << "Suggested configuration:" << std::endl;
        std::cout << "  Minimizer capacity: 15,000,000 - 25,000,000" << std::endl;
        std::cout << "  Sequence batch size: 40-50" << std::endl;
    } else if (total_gb >= 16) {
        std::cout << "Medium-memory GPU detected (" << total_gb << " GB)" << std::endl;
        std::cout << "Suggested configuration:" << std::endl;
        std::cout << "  Minimizer capacity: 8,000,000 - 15,000,000" << std::endl;
        std::cout << "  Sequence batch size: 25-35" << std::endl;
    } else if (total_gb >= 8) {
        std::cout << "Standard-memory GPU detected (" << total_gb << " GB)" << std::endl;
        std::cout << "Suggested configuration:" << std::endl;
        std::cout << "  Minimizer capacity: 3,000,000 - 8,000,000" << std::endl;
        std::cout << "  Sequence batch size: 15-25" << std::endl;
    } else {
        std::cout << "Limited-memory GPU detected (" << total_gb << " GB)" << std::endl;
        std::cout << "Suggested configuration:" << std::endl;
        std::cout << "  Minimizer capacity: 1,000,000 - 3,000,000" << std::endl;
        std::cout << "  Sequence batch size: 5-15" << std::endl;
        std::cout << "  Consider using streaming mode for large datasets" << std::endl;
    }
}

// Private implementation methods

bool GPUMemoryManager::query_gpu_memory_info() {
    cudaError_t error = cudaMemGetInfo(&stats_.available_memory, &stats_.total_gpu_memory);
    if (error != cudaSuccess) {
        std::cerr << "Failed to query GPU memory: " << cudaGetErrorString(error) << std::endl;
        return false;
    }
    
    stats_.allocated_memory = stats_.total_gpu_memory - stats_.available_memory;
    return true;
}

bool GPUMemoryManager::calculate_optimal_batch_sizes() {
    // Calculate optimal minimizer capacity based on available memory
    size_t usable_memory = stats_.available_memory * config_.max_memory_fraction / 100;
    usable_memory -= config_.reserved_memory_mb * 1024 * 1024;  // Reserve memory for system
    
    if (usable_memory <= 0) {
        std::cerr << "Insufficient GPU memory available" << std::endl;
        return false;
    }
    
    // Estimate memory distribution
    size_t minimizer_memory = usable_memory * 0.6;  // 60% for minimizers
    size_t sequence_memory = usable_memory * 0.3;   // 30% for sequences
    size_t overhead_memory = usable_memory * 0.1;   // 10% overhead
    
    // Calculate optimal capacities
    config_.minimizer_capacity = minimizer_memory / sizeof(GPUMinimizerHit);
    config_.minimizer_capacity = std::min(config_.minimizer_capacity, 50000000);  // Cap at 50M
    config_.minimizer_capacity = std::max(config_.minimizer_capacity, 1000000);   // Min 1M
    
    // Calculate optimal batch size
    size_t sequence_per_mb = 1024 * 1024 / 5000;  // Assume 5KB per sequence on average
    config_.sequence_batch_size = (sequence_memory / 1024 / 1024) / sequence_per_mb;
    config_.sequence_batch_size = std::min(config_.sequence_batch_size, 100);  // Cap at 100
    config_.sequence_batch_size = std::max(config_.sequence_batch_size, 5);    // Min 5
    
    std::cout << "Calculated optimal configuration:" << std::endl;
    std::cout << "  Minimizer capacity: " << config_.minimizer_capacity << std::endl;
    std::cout << "  Sequence batch size: " << config_.sequence_batch_size << std::endl;
    
    return true;
}

bool GPUMemoryManager::validate_gpu_context() {
    int device_count;
    cudaError_t error = cudaGetDeviceCount(&device_count);
    if (error != cudaSuccess || device_count == 0) {
        return false;
    }
    
    int current_device;
    error = cudaGetDevice(&current_device);
    return error == cudaSuccess;
}

void GPUMemoryManager::update_memory_statistics() {
    query_gpu_memory_info();
    stats_.memory_efficiency = calculate_memory_efficiency();
    
    if (total_allocated_ > stats_.peak_usage) {
        stats_.peak_usage = total_allocated_;
    }
}

bool GPUMemoryManager::validate_allocation_request(size_t size) const {
    if (size == 0) return false;
    
    size_t available_for_allocation = stats_.available_memory * config_.max_memory_fraction / 100;
    available_for_allocation -= config_.reserved_memory_mb * 1024 * 1024;
    
    if (size > available_for_allocation) {
        std::cerr << "Allocation request (" << (size / 1024 / 1024) << " MB) "
                  << "exceeds available memory (" << (available_for_allocation / 1024 / 1024) << " MB)" << std::endl;
        return false;
    }
    
    return true;
}

size_t GPUMemoryManager::calculate_safe_minimizer_capacity() const {
    size_t memory_gb = stats_.total_gpu_memory / (1024 * 1024 * 1024);
    
    if (memory_gb >= 24) return 15000000;      // 15M for high-end GPUs
    else if (memory_gb >= 16) return 10000000; // 10M for mid-range GPUs
    else if (memory_gb >= 8) return 5000000;   // 5M for standard GPUs
    else return 2000000;                       // 2M for low-end GPUs
}

int GPUMemoryManager::calculate_optimal_batch_size() const {
    size_t memory_gb = stats_.total_gpu_memory / (1024 * 1024 * 1024);
    
    if (memory_gb >= 24) return 50;
    else if (memory_gb >= 16) return 35;
    else if (memory_gb >= 8) return 25;
    else return 15;
}

double GPUMemoryManager::calculate_memory_efficiency() const {
    if (stats_.peak_usage == 0) return 1.0;
    return (double)total_allocated_ / stats_.peak_usage;
}
