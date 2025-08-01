#ifndef GPU_MEMORY_POOL_H
#define GPU_MEMORY_POOL_H

#include <cuda_runtime.h>
#include <mutex>
#include <memory>
#include <unordered_map>
#include <list>
#include <iostream>
#include <cassert>

namespace BioGPU {

// GPU Memory Pool for efficient memory management
class GPUMemoryPool {
private:
    struct MemoryBlock {
        void* ptr;
        size_t size;
        bool in_use;
        int device_id;
        
        MemoryBlock(void* p, size_t s, int dev) 
            : ptr(p), size(s), in_use(false), device_id(dev) {}
    };
    
    // Pool per GPU device
    std::unordered_map<int, std::list<std::unique_ptr<MemoryBlock>>> device_pools_;
    std::mutex pool_mutex_;
    
    // Statistics
    size_t total_allocated_ = 0;
    size_t total_in_use_ = 0;
    size_t peak_usage_ = 0;
    
    // Configuration
    size_t max_pool_size_ = 4ULL * 1024 * 1024 * 1024; // 4GB default
    bool enable_logging_ = false;
    
    static GPUMemoryPool* instance_;
    static std::mutex instance_mutex_;
    
    GPUMemoryPool() = default;
    
public:
    static GPUMemoryPool& getInstance() {
        std::lock_guard<std::mutex> lock(instance_mutex_);
        if (!instance_) {
            instance_ = new GPUMemoryPool();
        }
        return *instance_;
    }
    
    // Allocate memory from pool
    void* allocate(size_t size, int device_id = -1) {
        if (device_id < 0) {
            cudaGetDevice(&device_id);
        }
        
        std::lock_guard<std::mutex> lock(pool_mutex_);
        
        // Try to find a free block of sufficient size
        auto& pool = device_pools_[device_id];
        for (auto& block : pool) {
            if (!block->in_use && block->size >= size && block->size <= size * 1.5) {
                block->in_use = true;
                total_in_use_ += block->size;
                peak_usage_ = std::max(peak_usage_, total_in_use_);
                
                if (enable_logging_) {
                    std::cout << "[GPU Memory Pool] Reused block: " << block->size 
                              << " bytes on GPU " << device_id << std::endl;
                }
                return block->ptr;
            }
        }
        
        // No suitable block found, allocate new
        void* ptr = nullptr;
        cudaError_t err = cudaMalloc(&ptr, size);
        if (err != cudaSuccess || !ptr) {
            std::cerr << "[GPU Memory Pool] Allocation failed: " << size 
                      << " bytes on GPU " << device_id 
                      << " - " << cudaGetErrorString(err) << std::endl;
            return nullptr;
        }
        
        // Add to pool
        pool.emplace_back(std::make_unique<MemoryBlock>(ptr, size, device_id));
        pool.back()->in_use = true;
        total_allocated_ += size;
        total_in_use_ += size;
        peak_usage_ = std::max(peak_usage_, total_in_use_);
        
        if (enable_logging_) {
            std::cout << "[GPU Memory Pool] New allocation: " << size 
                      << " bytes on GPU " << device_id 
                      << " (Total: " << (total_allocated_ / (1024*1024)) << " MB)" << std::endl;
        }
        
        return ptr;
    }
    
    // Release memory back to pool
    void deallocate(void* ptr) {
        if (!ptr) return;
        
        std::lock_guard<std::mutex> lock(pool_mutex_);
        
        for (auto& [device_id, pool] : device_pools_) {
            for (auto& block : pool) {
                if (block->ptr == ptr) {
                    assert(block->in_use);
                    block->in_use = false;
                    total_in_use_ -= block->size;
                    
                    if (enable_logging_) {
                        std::cout << "[GPU Memory Pool] Released block: " << block->size 
                                  << " bytes on GPU " << device_id << std::endl;
                    }
                    return;
                }
            }
        }
        
        std::cerr << "[GPU Memory Pool] Warning: Attempting to deallocate unknown pointer" << std::endl;
    }
    
    // Clear unused blocks to free memory
    void cleanup(int device_id = -1) {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        
        if (device_id < 0) {
            // Clean all devices
            for (auto& [dev_id, pool] : device_pools_) {
                cleanupDevice(dev_id);
            }
        } else {
            cleanupDevice(device_id);
        }
    }
    
    // Get memory statistics
    void printStats() const {
        std::lock_guard<std::mutex> lock(const_cast<std::mutex&>(pool_mutex_));
        
        std::cout << "\n=== GPU Memory Pool Statistics ===" << std::endl;
        std::cout << "Total allocated: " << (total_allocated_ / (1024*1024)) << " MB" << std::endl;
        std::cout << "Currently in use: " << (total_in_use_ / (1024*1024)) << " MB" << std::endl;
        std::cout << "Peak usage: " << (peak_usage_ / (1024*1024)) << " MB" << std::endl;
        
        for (const auto& [device_id, pool] : device_pools_) {
            size_t device_total = 0;
            size_t device_in_use = 0;
            int free_blocks = 0;
            int used_blocks = 0;
            
            for (const auto& block : pool) {
                device_total += block->size;
                if (block->in_use) {
                    device_in_use += block->size;
                    used_blocks++;
                } else {
                    free_blocks++;
                }
            }
            
            std::cout << "\nGPU " << device_id << ":" << std::endl;
            std::cout << "  Blocks: " << used_blocks << " in use, " 
                      << free_blocks << " free" << std::endl;
            std::cout << "  Memory: " << (device_in_use / (1024*1024)) << " MB in use, " 
                      << ((device_total - device_in_use) / (1024*1024)) << " MB free" << std::endl;
        }
        std::cout << "=================================" << std::endl;
    }
    
    // Enable/disable logging
    void setLogging(bool enable) { enable_logging_ = enable; }
    
    // Set maximum pool size
    void setMaxPoolSize(size_t size) { max_pool_size_ = size; }
    
private:
    void cleanupDevice(int device_id) {
        auto it = device_pools_.find(device_id);
        if (it == device_pools_.end()) return;
        
        auto& pool = it->second;
        auto block_it = pool.begin();
        while (block_it != pool.end()) {
            if (!(*block_it)->in_use) {
                cudaFree((*block_it)->ptr);
                total_allocated_ -= (*block_it)->size;
                block_it = pool.erase(block_it);
            } else {
                ++block_it;
            }
        }
    }
    
    // Destructor - clean up all allocated memory
    ~GPUMemoryPool() {
        for (auto& [device_id, pool] : device_pools_) {
            for (auto& block : pool) {
                cudaFree(block->ptr);
            }
        }
    }
};

// Static members
GPUMemoryPool* GPUMemoryPool::instance_ = nullptr;
std::mutex GPUMemoryPool::instance_mutex_;

// Convenience wrapper class for RAII memory management
template<typename T>
class GPUMemoryHandle {
private:
    T* ptr_;
    size_t count_;
    
public:
    GPUMemoryHandle() : ptr_(nullptr), count_(0) {}
    
    explicit GPUMemoryHandle(size_t count) : count_(count) {
        if (count > 0) {
            ptr_ = static_cast<T*>(GPUMemoryPool::getInstance().allocate(count * sizeof(T)));
        } else {
            ptr_ = nullptr;
        }
    }
    
    ~GPUMemoryHandle() {
        if (ptr_) {
            GPUMemoryPool::getInstance().deallocate(ptr_);
        }
    }
    
    // Move semantics
    GPUMemoryHandle(GPUMemoryHandle&& other) noexcept 
        : ptr_(other.ptr_), count_(other.count_) {
        other.ptr_ = nullptr;
        other.count_ = 0;
    }
    
    GPUMemoryHandle& operator=(GPUMemoryHandle&& other) noexcept {
        if (this != &other) {
            if (ptr_) {
                GPUMemoryPool::getInstance().deallocate(ptr_);
            }
            ptr_ = other.ptr_;
            count_ = other.count_;
            other.ptr_ = nullptr;
            other.count_ = 0;
        }
        return *this;
    }
    
    // Delete copy constructor and assignment
    GPUMemoryHandle(const GPUMemoryHandle&) = delete;
    GPUMemoryHandle& operator=(const GPUMemoryHandle&) = delete;
    
    // Accessors
    T* get() { return ptr_; }
    const T* get() const { return ptr_; }
    size_t size() const { return count_; }
    bool valid() const { return ptr_ != nullptr; }
    
    // Release ownership
    T* release() {
        T* temp = ptr_;
        ptr_ = nullptr;
        count_ = 0;
        return temp;
    }
};

} // namespace BioGPU

#endif // GPU_MEMORY_POOL_H