// feature_tests.cpp
// Modular feature test implementations

#include "test_framework.h"
#include "features/uniqueness_score_implementation.cu"
#include "features/cooccurrence_scoring.cu"
#include "features/namespace_conflict_resolution.h"

// ===========================
// GC Content Test
// ===========================
DEFINE_FEATURE_TEST(GCContentTest, "GC Content Analysis", enable_gc_content)
    std::vector<int> gc_distribution_;
    
    bool initialize() override {
        gc_distribution_.resize(8, 0);
        return true;
    }
    
    bool processOnGPU(GPUMinimizerHit* d_hits, size_t num_hits, 
                     const std::vector<uint32_t>& taxon_ids) override {
        // GC content is calculated during extraction, just verify it's there
        return true;
    }
    
    void analyzeResults(const std::vector<GPUMinimizerHit>& hits) override {
        for (const auto& hit : hits) {
            uint8_t gc_cat = hit.feature_flags & 0x7;
            if (gc_cat < 8) gc_distribution_[gc_cat]++;
        }
    }
    
    void printSummary() const override {
        const char* gc_names[] = {
            "0-20%", "20-30%", "30-40%", "40-50%",
            "50-60%", "60-70%", "70-80%", "80-100%"
        };
        
        std::cout << "GC Content Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            std::cout << "  " << gc_names[i] << ": " << gc_distribution_[i] << std::endl;
        }
    }
};

// ===========================
// Uniqueness Test
// ===========================
DEFINE_FEATURE_TEST(UniquenessTest, "Uniqueness Scoring", enable_uniqueness)
    std::vector<int> uniqueness_distribution_;
    size_t unique_count_ = 0;
    size_t rare_count_ = 0;
    size_t reliable_count_ = 0;
    double avg_uniqueness_ = 0.0;
    
    bool initialize() override {
        uniqueness_distribution_.resize(8, 0);
        return true;
    }
    
    bool processOnGPU(GPUMinimizerHit* d_hits, size_t num_hits, 
                     const std::vector<uint32_t>& taxon_ids) override {
        // Call the uniqueness computation
        return compute_and_encode_uniqueness_scores(
            d_hits, num_hits, taxon_ids, taxon_ids.size()
        );
    }
    
    void analyzeResults(const std::vector<GPUMinimizerHit>& hits) override {
        double uniqueness_sum = 0.0;
        
        for (const auto& hit : hits) {
            uint8_t cat = MinimizerFlags::get_uniqueness_category_safe(hit.feature_flags);
            if (cat < 8) uniqueness_distribution_[cat]++;
            
            if (MinimizerFlags::is_unique_minimizer_safe(hit.feature_flags)) unique_count_++;
            if (MinimizerFlags::is_rare_minimizer_safe(hit.feature_flags)) rare_count_++;
            if (MinimizerFlags::is_reliable_minimizer_safe(hit.feature_flags)) reliable_count_++;
            
            uniqueness_sum += MinimizerFlags::ml_weight_to_uniqueness(hit.ml_weight);
        }
        
        if (!hits.empty()) {
            avg_uniqueness_ = uniqueness_sum / hits.size();
        }
    }
    
    void printSummary() const override {
        const char* names[] = {
            "Extremely common (0-10%)", "Very low (10-30%)",
            "Low (30-50%)", "Moderate (50-70%)",
            "High (70-80%)", "Very high (80-90%)",
            "Extremely high (90-95%)", "Singleton-like (95-100%)"
        };
        
        std::cout << "Uniqueness Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            std::cout << "  " << names[i] << ": " << uniqueness_distribution_[i] << std::endl;
        }
        
        std::cout << "\nQuality Metrics:" << std::endl;
        std::cout << "  Unique minimizers (≥90%): " << unique_count_ << std::endl;
        std::cout << "  Rare minimizers (≤3 occurrences): " << rare_count_ << std::endl;
        std::cout << "  Reliable minimizers: " << reliable_count_ << std::endl;
        std::cout << "  Average uniqueness score: " << std::fixed 
                  << std::setprecision(3) << avg_uniqueness_ << std::endl;
    }
};

// ===========================
// Co-occurrence Test
// ===========================
DEFINE_FEATURE_TEST(CooccurrenceTest, "Co-occurrence Analysis", enable_cooccurrence)
    std::vector<int> cooccurrence_distribution_;
    
    bool initialize() override {
        cooccurrence_distribution_.resize(8, 0);
        return true;
    }
    
    bool processOnGPU(GPUMinimizerHit* d_hits, size_t num_hits, 
                     const std::vector<uint32_t>& taxon_ids) override {
        // Build unique minimizer list
        std::vector<GPUMinimizerHit> h_hits(num_hits);
        cudaMemcpy(h_hits.data(), d_hits, num_hits * sizeof(GPUMinimizerHit), 
                   cudaMemcpyDeviceToHost);
        
        std::unordered_map<uint64_t, uint32_t> counts;
        for (const auto& hit : h_hits) {
            counts[hit.minimizer_hash]++;
        }
        
        std::vector<uint64_t> unique_minimizers;
        std::vector<uint32_t> occurrence_counts;
        for (const auto& [hash, count] : counts) {
            unique_minimizers.push_back(hash);
            occurrence_counts.push_back(count);
        }
        
        // Create dummy genome info
        std::vector<GPUGenomeInfo> genome_info;
        for (size_t i = 0; i < taxon_ids.size(); i++) {
            GPUGenomeInfo info = {i, taxon_ids[i], 100000, 0};
            genome_info.push_back(info);
        }
        
        return compute_and_encode_cooccurrence_scores(
            d_hits, num_hits, unique_minimizers, occurrence_counts, genome_info
        );
    }
    
    void analyzeResults(const std::vector<GPUMinimizerHit>& hits) override {
        for (const auto& hit : hits) {
            uint8_t score = MinimizerFlags::get_cooccurrence_score(hit.feature_flags);
            if (score < 8) cooccurrence_distribution_[score]++;
        }
    }
    
    void printSummary() const override {
        const char* names[] = {
            "No co-occurrence", "Very low", "Low", "Moderate",
            "High", "Very high", "Extremely high", "Perfect pattern"
        };
        
        std::cout << "Co-occurrence Distribution:" << std::endl;
        for (int i = 0; i < 8; i++) {
            std::cout << "  " << names[i] << ": " << cooccurrence_distribution_[i] << std::endl;
        }
    }
};

// ===========================
// Main test runner
// ===========================
int main(int argc, char** argv) {
    // Create test configuration
    TestConfig config;
    
    // Parse command line arguments to override defaults
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--fna" && i + 1 < argc) {
            config.fna_file = argv[++i];
        } else if (arg == "--batch-size" && i + 1 < argc) {
            config.batch_size = std::stoi(argv[++i]);
        } else if (arg == "--no-uniqueness") {
            config.enable_uniqueness = false;
        } else if (arg == "--no-cooccurrence") {
            config.enable_cooccurrence = false;
        } else if (arg == "--examples" && i + 1 < argc) {
            config.num_examples_to_print = std::stoi(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --fna FILE              Input FNA file" << std::endl;
            std::cout << "  --batch-size N          Process N genomes per batch" << std::endl;
            std::cout << "  --no-uniqueness         Disable uniqueness scoring" << std::endl;
            std::cout << "  --no-cooccurrence       Disable co-occurrence analysis" << std::endl;
            std::cout << "  --examples N            Print N example minimizers" << std::endl;
            return 0;
        }
    }
    
    // Create framework
    BioGPUTestFramework framework(config);
    
    // Register all feature tests
    framework.registerFeatureTest(std::make_unique<GCContentTest>());
    framework.registerFeatureTest(std::make_unique<UniquenessTest>());
    framework.registerFeatureTest(std::make_unique<CooccurrenceTest>());
    
    // Run all tests
    bool success = framework.runAllTests();
    
    return success ? 0 : 1;
}

// ===========================
// Easy extension example
// ===========================
// To add a new feature test, just create a new class:
/*
DEFINE_FEATURE_TEST(MyNewFeatureTest, "My New Feature", enable_my_feature)
    // Your data members
    
    bool initialize() override {
        // Setup code
        return true;
    }
    
    bool processOnGPU(GPUMinimizerHit* d_hits, size_t num_hits, 
                     const std::vector<uint32_t>& taxon_ids) override {
        // GPU processing
        return true;
    }
    
    void analyzeResults(const std::vector<GPUMinimizerHit>& hits) override {
        // Analyze results
    }
    
    void printSummary() const override {
        // Print summary
    }
};

// Then in main(), just add:
framework.registerFeatureTest(std::make_unique<MyNewFeatureTest>());
*/
