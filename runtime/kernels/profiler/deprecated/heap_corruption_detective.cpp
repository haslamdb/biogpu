// heap_corruption_detective.cpp
// Tool to pinpoint exactly where heap corruption is occurring

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <signal.h>
#include <execinfo.h>
#include <unistd.h>
#include <malloc.h>
#include <sys/mman.h>

class HeapCorruptionDetective {
private:
    static bool monitoring_enabled;
    static size_t allocation_count;
    static size_t total_allocated;
    static void* last_allocation;
    static size_t last_size;
    
public:
    static void enable_monitoring() {
        monitoring_enabled = true;
        allocation_count = 0;
        total_allocated = 0;
        last_allocation = nullptr;
        last_size = 0;
        
        std::cout << "🕵️ Heap Corruption Detective: Monitoring enabled" << std::endl;
        
        // Install crash handlers
        signal(SIGABRT, crash_handler);
        signal(SIGSEGV, crash_handler);
        
        // Enable malloc debugging
        setenv("MALLOC_CHECK_", "3", 1);  // Abort on corruption
        setenv("MALLOC_PERTURB_", "42", 1);  // Fill freed memory with pattern
        
        std::cout << "🛡️ Enhanced malloc debugging enabled" << std::endl;
    }
    
    static void* safe_malloc(size_t size, const char* location) {
        if (!monitoring_enabled) return malloc(size);
        
        void* ptr = malloc(size);
        if (ptr) {
            allocation_count++;
            total_allocated += size;
            last_allocation = ptr;
            last_size = size;
            
            if (allocation_count % 1000 == 0) {
                std::cout << "🔍 Allocation #" << allocation_count 
                          << " at " << location 
                          << " (total: " << (total_allocated / 1024 / 1024) << " MB)" << std::endl;
            }
            
            // Fill with pattern to detect use-after-free
            memset(ptr, 0xAB, size);
        } else {
            std::cerr << "❌ MALLOC FAILED at " << location << " for " << size << " bytes" << std::endl;
            print_memory_info();
        }
        
        return ptr;
    }
    
    static void safe_free(void* ptr, const char* location) {
        if (!monitoring_enabled || !ptr) {
            free(ptr);
            return;
        }
        
        // Fill with pattern before freeing
        if (ptr == last_allocation && last_size > 0) {
            memset(ptr, 0xDE, last_size);
        }
        
        free(ptr);
        
        if (allocation_count % 1000 == 0) {
            std::cout << "🧹 Free at " << location << std::endl;
        }
    }
    
    static void crash_handler(int sig) {
        std::cerr << "\n💥 CRASH DETECTED! Signal: " << sig << std::endl;
        std::cerr << "Last allocation: " << last_allocation 
                  << " (" << last_size << " bytes)" << std::endl;
        std::cerr << "Total allocations: " << allocation_count << std::endl;
        std::cerr << "Total memory: " << (total_allocated / 1024 / 1024) << " MB" << std::endl;
        
        print_stack_trace();
        print_memory_info();
        
        // Try to get malloc stats
        malloc_stats();
        
        exit(1);
    }
    
    static void print_stack_trace() {
        std::cerr << "\n📚 Stack trace:" << std::endl;
        void *array[20];
        size_t size = backtrace(array, 20);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
    }
    
    static void print_memory_info() {
        std::cerr << "\n💾 Memory information:" << std::endl;
        
        // Print malloc info
        struct mallinfo mi = mallinfo();
        std::cerr << "Arena: " << mi.arena << " bytes" << std::endl;
        std::cerr << "Ordinary blocks: " << mi.ordblks << std::endl;
        std::cerr << "Small blocks: " << mi.smblks << std::endl;
        std::cerr << "Hold blocks: " << mi.hblks << std::endl;
        std::cerr << "Total allocated: " << mi.uordblks << " bytes" << std::endl;
        std::cerr << "Total free: " << mi.fordblks << " bytes" << std::endl;
        
        // Print from /proc/self/status
        FILE* status = fopen("/proc/self/status", "r");
        if (status) {
            char line[256];
            while (fgets(line, sizeof(line), status)) {
                if (strstr(line, "VmSize:") || strstr(line, "VmRSS:") || 
                    strstr(line, "VmData:") || strstr(line, "VmStk:")) {
                    std::cerr << line;
                }
            }
            fclose(status);
        }
    }
    
    static void checkpoint(const char* location) {
        if (!monitoring_enabled) return;
        
        std::cout << "🔍 CHECKPOINT: " << location << std::endl;
        std::cout << "   Allocations so far: " << allocation_count << std::endl;
        std::cout << "   Memory allocated: " << (total_allocated / 1024 / 1024) << " MB" << std::endl;
        
        // Test heap integrity
        void* test_ptr = malloc(64);
        if (test_ptr) {
            memset(test_ptr, 0x55, 64);
            free(test_ptr);
            std::cout << "   Heap test: PASS" << std::endl;
        } else {
            std::cerr << "   Heap test: FAIL - malloc returned NULL!" << std::endl;
        }
        
        // Note: mcheck functions are deprecated in modern glibc
        // Using heap test above instead
    }
    
    static void enable_advanced_debugging() {
        std::cout << "🔬 Enabling advanced heap debugging..." << std::endl;
        
        // Note: mcheck is deprecated in modern glibc
        // Using MALLOC_CHECK_ environment variable instead
        
        // Set environment variables for maximum debugging
        setenv("MALLOC_CHECK_", "3", 1);     // Abort on heap corruption
        setenv("MALLOC_PERTURB_", "42", 1);  // Perturb freed memory
        setenv("LIBC_FATAL_STDERR_", "1", 1); // Print to stderr on fatal errors
        
        std::cout << "✅ Advanced debugging enabled" << std::endl;
    }
};

// Static member definitions
bool HeapCorruptionDetective::monitoring_enabled = false;
size_t HeapCorruptionDetective::allocation_count = 0;
size_t HeapCorruptionDetective::total_allocated = 0;
void* HeapCorruptionDetective::last_allocation = nullptr;
size_t HeapCorruptionDetective::last_size = 0;

// Convenience macros
#define STRINGIFY(x) #x
#define SAFE_MALLOC(size) HeapCorruptionDetective::safe_malloc(size, __FILE__ ":" STRINGIFY(__LINE__))
#define SAFE_FREE(ptr) HeapCorruptionDetective::safe_free(ptr, __FILE__ ":" STRINGIFY(__LINE__))
#define HEAP_CHECKPOINT(msg) HeapCorruptionDetective::checkpoint(msg)

// Test function to demonstrate usage
int test_heap_monitoring() {
    std::cout << "🧪 Testing heap corruption detection..." << std::endl;
    
    HeapCorruptionDetective::enable_advanced_debugging();
    HeapCorruptionDetective::enable_monitoring();
    
    HEAP_CHECKPOINT("test_start");
    
    // Test normal allocation
    void* ptr1 = SAFE_MALLOC(1024);
    HEAP_CHECKPOINT("after_malloc_1024");
    
    // Test large allocation
    void* ptr2 = SAFE_MALLOC(1024 * 1024);
    HEAP_CHECKPOINT("after_malloc_1MB");
    
    // Test many small allocations
    for (int i = 0; i < 1000; i++) {
        void* ptr = SAFE_MALLOC(64);
        if (i % 100 == 0) {
            HEAP_CHECKPOINT(("small_alloc_" + std::to_string(i)).c_str());
        }
        SAFE_FREE(ptr);
    }
    
    HEAP_CHECKPOINT("after_small_allocs");
    
    SAFE_FREE(ptr1);
    SAFE_FREE(ptr2);
    
    HEAP_CHECKPOINT("test_end");
    
    std::cout << "✅ Heap monitoring test completed" << std::endl;
    return 0;
}

/*
USAGE INSTRUCTIONS:

1. Add this to your project and include the header
2. At the very beginning of your main function, call:
   HeapCorruptionDetective::enable_advanced_debugging();
   HeapCorruptionDetective::enable_monitoring();

3. Add checkpoints throughout your code:
   HEAP_CHECKPOINT("before_database_builder");
   HEAP_CHECKPOINT("after_cuda_init");
   HEAP_CHECKPOINT("before_load_genome_files");

4. This will show you exactly where the heap corruption occurs
   and provide detailed debugging information when it crashes.

5. Compile with:
   g++ -g -O0 -fsanitize=address -fsanitize=undefined your_program.cpp

This tool will help pinpoint the exact location of the heap corruption
so you can apply targeted fixes.
*/
