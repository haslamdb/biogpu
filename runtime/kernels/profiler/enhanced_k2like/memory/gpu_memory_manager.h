// memory/gpu_memory_manager.h
// GPU memory management for Kraken database builder
// Handles allocation, auto-scaling, and memory optimization

#ifndef GPU_MEMORY_MANAGER_H
#define GPU_MEMORY_MANAGER_H

#include <memory>
#include <cstdint>

// Include the main types header instead of redefining
#include "../gpu_kraken_types.h"

// Note: The following types are already defined in gpu_kraken_types.h:
// - MemoryConfig (lines 17-24)
// - MemoryStats (lines 27-36)
// - GPUGenomeInfo, GPUMinimizerHit, LCACandidate (forward declared at lines 22-24)

// GPU Memory Pool for efficient allocation
class GPUMemoryPool {
private:
    void* pool_base_;
    size_t pool_size_;
    size_t current_offset_;
    bool initialized_;

public:
    explicit GPUMemoryPool(size_t pool_size_mb);
    ~GPUMemoryPool();
    
    void* allocate(size_t size, size_t alignment = 256);
    void deallocate(void* ptr);
    void reset();
    size_t get_available_space() const;
    bool is_initialized() const;
};

// Main GPU Memory Manager
class GPUMemoryManager {
private:
    MemoryConfig config_;
    bool initialized_;
    bool allocations_active_;
    size_t total_allocated_;
    
    // Memory pool
    std::unique_ptr<GPUMemoryPool> memory_pool_;
    
    // GPU memory pointers
    char* d_sequence_data_;
    GPUGenomeInfo* d_genome_info_;
    GPUMinimizerHit* d_minimizer_hits_;
    LCACandidate* d_lca_candidates_;
    uint32_t* d_minimizer_counts_;
    
    // Statistics
    MemoryStats stats_;
    
public:
    explicit GPUMemoryManager(const MemoryConfig& config = MemoryConfig());
    ~GPUMemoryManager();
    
    // Initialization and configuration
    bool initialize();
    bool configure_auto_scaling(bool enable, size_t memory_fraction = 80);
    bool set_minimizer_capacity(int capacity);
    bool set_batch_size(int batch_size);
    
    // Memory allocation methods
    bool allocate_sequence_memory(size_t max_sequences, size_t max_total_length);
    bool allocate_minimizer_memory(size_t max_minimizers);
    bool allocate_metadata_memory(size_t max_genomes);
    bool allocate_results_memory(size_t max_candidates);
    void free_all_allocations();
    
    // Memory validation and integrity
    bool validate_memory_integrity();
    
    // Accessor methods for GPU pointers
    char* get_sequence_buffer() const;
    GPUGenomeInfo* get_genome_info_buffer() const;
    GPUMinimizerHit* get_minimizer_buffer() const;
    LCACandidate* get_candidate_buffer() const;
    uint32_t* get_count_buffer() const;
    
    // Statistics and monitoring
    MemoryStats get_memory_statistics() const;
    void print_memory_usage() const;
    bool check_memory_pressure() const;
    size_t estimate_memory_requirements(int num_sequences, int avg_sequence_length) const;
    
    // Auto-scaling methods
    bool scale_for_workload(size_t estimated_minimizers, size_t estimated_sequences);
    void suggest_optimal_configuration() const;
    
    // Additional methods for accessing internal state
    uint32_t* get_global_counter() const { return d_minimizer_counts_; }
    size_t get_minimizer_capacity() const { return config_.minimizer_capacity; }
    
private:
    // Internal methods
    bool query_gpu_memory_info();
    bool calculate_optimal_batch_sizes();
    bool validate_gpu_context();
    void update_memory_statistics();
    bool validate_allocation_request(size_t size) const;
    size_t calculate_safe_minimizer_capacity() const;
    int calculate_optimal_batch_size() const;
    double calculate_memory_efficiency() const;
};

// Utility functions for memory debugging and validation
namespace MemoryUtils {
    void check_cuda_memory_leaks();
    bool validate_gpu_pointer(void* ptr, size_t size);
    size_t get_gpu_memory_usage();
    void print_gpu_memory_info();
}

#endif // GPU_MEMORY_MANAGER_H
