// heap_corruption_detective.h
// Header for heap corruption detection tool

#ifndef HEAP_CORRUPTION_DETECTIVE_H
#define HEAP_CORRUPTION_DETECTIVE_H

#include <cstddef>

class HeapCorruptionDetective {
private:
    static bool monitoring_enabled;
    static size_t allocation_count;
    static size_t total_allocated;
    static void* last_allocation;
    static size_t last_size;
    
public:
    static void enable_monitoring();
    static void* safe_malloc(size_t size, const char* location);
    static void safe_free(void* ptr, const char* location);
    static void crash_handler(int sig);
    static void print_stack_trace();
    static void print_memory_info();
    static void checkpoint(const char* location);
    static void enable_advanced_debugging();
};

// Convenience macros
#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)
#define SAFE_MALLOC(size) HeapCorruptionDetective::safe_malloc(size, __FILE__ ":" STRINGIFY(__LINE__))
#define SAFE_FREE(ptr) HeapCorruptionDetective::safe_free(ptr, __FILE__ ":" STRINGIFY(__LINE__))
#define HEAP_CHECKPOINT(msg) HeapCorruptionDetective::checkpoint(msg)

#endif // HEAP_CORRUPTION_DETECTIVE_H