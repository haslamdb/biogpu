// output/feature_exporter.cu
// Feature export pipeline for ML training data generation
// Exports minimizer features in multiple formats for downstream ML applications

#include "feature_exporter.h"
#include "../gpu_kraken_types.h"
#include "../processing/minimizer_feature_extractor.h"
#include "../processing/contamination_detector.h"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <chrono>

#ifdef USE_HDF5
#include <H5Cpp.h>
#endif

// ===========================
// Device Functions
// ===========================

__device__ float calculate_position_std_dev(float mean, float sum_squares, int count) {
    if (count <= 1) return 0.0f;
    float variance = (sum_squares - mean * mean * count) / (count - 1);
    return sqrtf(fmaxf(0.0f, variance));
}

__device__ void extract_feature_values(
    const GPUMinimizerHit& hit,
    MinimizerFeatureValues& features) {
    
    // Extract encoded features from feature_flags
    features.gc_category = hit.feature_flags & 0x7;
    features.complexity_score = (hit.feature_flags >> 3) & 0x7;
    features.has_position_bias = (hit.feature_flags & (1 << 6)) != 0;
    features.is_contamination = (hit.feature_flags & (1 << 7)) != 0;
    features.conservation_category = (hit.feature_flags >> 8) & 0x7;
    features.is_unique = (hit.feature_flags & (1 << 11)) != 0;
    features.is_rare = (hit.feature_flags & (1 << 12)) != 0;
    
    // Decode ML weight to float
    features.uniqueness_score = (float)hit.ml_weight / 65535.0f;
}

// ===========================
// CUDA Kernels
// ===========================

__global__ void aggregate_minimizer_features_kernel(
    const GPUMinimizerHit* minimizer_hits,
    const uint64_t* unique_minimizers,
    AggregatedMinimizerFeatures* aggregated_features,
    int num_hits,
    int num_unique) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_unique) return;
    
    uint64_t target_hash = unique_minimizers[tid];
    AggregatedMinimizerFeatures& agg = aggregated_features[tid];
    
    // Initialize aggregated features
    agg.minimizer_hash = target_hash;
    agg.total_occurrences = 0;
    agg.mean_position = 0.0f;
    agg.position_std_dev = 0.0f;
    agg.mean_gc_content = 0.0f;
    agg.mean_complexity = 0.0f;
    agg.contamination_frequency = 0.0f;
    
    // Aggregate features from all occurrences
    float position_sum = 0.0f;
    float position_sum_squares = 0.0f;
    int gc_sum = 0;
    int complexity_sum = 0;
    int contamination_count = 0;
    
    for (int i = 0; i < num_hits; i++) {
        if (minimizer_hits[i].minimizer_hash == target_hash) {
            agg.total_occurrences++;
            
            float pos = (float)minimizer_hits[i].position;
            position_sum += pos;
            position_sum_squares += pos * pos;
            
            MinimizerFeatureValues features;
            extract_feature_values(minimizer_hits[i], features);
            
            gc_sum += features.gc_category;
            complexity_sum += features.complexity_score;
            if (features.is_contamination) contamination_count++;
            
            // Store first occurrence details
            if (agg.total_occurrences == 1) {
                agg.avg_uniqueness_score = features.uniqueness_score;
                agg.conservation_level = features.conservation_category;
            }
        }
    }
    
    // Calculate statistics
    if (agg.total_occurrences > 0) {
        agg.mean_position = position_sum / agg.total_occurrences;
        agg.position_std_dev = calculate_position_std_dev(
            agg.mean_position, position_sum_squares, agg.total_occurrences
        );
        agg.mean_gc_content = (float)gc_sum / agg.total_occurrences;
        agg.mean_complexity = (float)complexity_sum / agg.total_occurrences;
        agg.contamination_frequency = (float)contamination_count / agg.total_occurrences;
    }
}

__global__ void compute_cooccurrence_kernel(
    const GPUMinimizerHit* minimizer_hits,
    CooccurrenceMatrix* cooccurrence,
    int num_hits,
    int window_size) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_hits) return;
    
    uint64_t center_hash = minimizer_hits[tid].minimizer_hash;
    uint16_t center_genome = minimizer_hits[tid].genome_id;
    uint32_t center_pos = minimizer_hits[tid].position;
    
    // Look for neighboring minimizers within window
    for (int i = max(0, tid - 100); i < min(num_hits, tid + 100); i++) {
        if (i == tid) continue;
        
        const GPUMinimizerHit& neighbor = minimizer_hits[i];
        if (neighbor.genome_id == center_genome) {
            int pos_diff = abs((int)neighbor.position - (int)center_pos);
            if (pos_diff <= window_size) {
                // Record co-occurrence (simplified - would need atomic operations)
                // In production, use a more sophisticated approach
            }
        }
    }
}

// ===========================
// FeatureExporter Implementation
// ===========================

FeatureExporter::FeatureExporter(const FeatureExportConfig& config)
    : config_(config), total_minimizers_exported_(0), export_start_time_(0) {
    
    std::cout << "Initializing Feature Exporter with configuration:" << std::endl;
    std::cout << "  Output directory: " << config_.output_directory << std::endl;
    std::cout << "  Export formats: ";
    if (config_.export_tsv) std::cout << "TSV ";
    if (config_.export_hdf5) std::cout << "HDF5 ";
    if (config_.export_binary) std::cout << "Binary ";
    std::cout << std::endl;
}

bool FeatureExporter::export_training_features(
    const std::vector<GPUMinimizerHit>& minimizer_hits,
    const MinimizerFeatureExtractor& feature_extractor,
    const std::unordered_map<uint32_t, std::string>& taxon_names) {
    
    export_start_time_ = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n=== EXPORTING MINIMIZER FEATURES FOR ML TRAINING ===" << std::endl;
    std::cout << "Processing " << minimizer_hits.size() << " minimizer hits" << std::endl;
    
    // Collect and aggregate features
    if (!collect_and_aggregate_features(minimizer_hits, feature_extractor)) {
        std::cerr << "Failed to collect and aggregate features" << std::endl;
        return false;
    }
    
    // Calculate co-occurrence patterns if requested
    if (config_.include_cooccurrence) {
        calculate_cooccurrence_patterns(minimizer_hits);
    }
    
    // Export in requested formats
    bool success = true;
    
    if (config_.export_tsv) {
        success &= export_to_tsv(taxon_names);
    }
    
    if (config_.export_hdf5) {
        success &= export_to_hdf5(taxon_names);
    }
    
    if (config_.export_binary) {
        success &= export_to_binary();
    }
    
    // Generate summary report
    generate_export_summary();
    
    auto export_end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(export_end_time - export_start_time_).count();
    
    std::cout << "\nFeature export completed in " << std::fixed << std::setprecision(2) 
              << duration << " seconds" << std::endl;
    std::cout << "Total minimizers exported: " << total_minimizers_exported_ << std::endl;
    
    return success;
}

bool FeatureExporter::export_to_tsv(const std::unordered_map<uint32_t, std::string>& taxon_names) {
    std::string tsv_path = config_.output_directory + "/minimizer_features.tsv";
    std::cout << "Exporting features to TSV: " << tsv_path << std::endl;
    
    std::ofstream out(tsv_path);
    if (!out.is_open()) {
        std::cerr << "Cannot create TSV file: " << tsv_path << std::endl;
        return false;
    }
    
    // Write header
    out << "minimizer_hash\t"
        << "total_occurrences\t"
        << "unique_taxa\t"
        << "mean_position\t"
        << "position_std_dev\t"
        << "gc_category\t"
        << "mean_complexity\t"
        << "uniqueness_score\t"
        << "conservation_level\t"
        << "contamination_freq\t";
    
    if (config_.include_taxonomic_distribution) {
        out << "taxonomic_distribution\t";
    }
    
    if (config_.include_cooccurrence) {
        out << "top_cooccurring_minimizers";
    }
    
    out << "\n";
    
    // Write data
    for (const auto& agg_feature : aggregated_features_) {
        out << std::hex << agg_feature.minimizer_hash << std::dec << "\t"
            << agg_feature.total_occurrences << "\t"
            << agg_feature.unique_taxa << "\t"
            << std::fixed << std::setprecision(2) << agg_feature.mean_position << "\t"
            << std::setprecision(2) << agg_feature.position_std_dev << "\t"
            << std::setprecision(2) << agg_feature.mean_gc_content << "\t"
            << std::setprecision(2) << agg_feature.mean_complexity << "\t"
            << std::setprecision(4) << agg_feature.avg_uniqueness_score << "\t"
            << (int)agg_feature.conservation_level << "\t"
            << std::setprecision(4) << agg_feature.contamination_frequency;
        
        if (config_.include_taxonomic_distribution) {
            out << "\t";
            // Write taxonomic distribution
            auto it = taxonomic_distributions_.find(agg_feature.minimizer_hash);
            if (it != taxonomic_distributions_.end()) {
                bool first = true;
                for (const auto& [taxon_id, count] : it->second) {
                    if (!first) out << ";";
                    out << taxon_id << ":" << count;
                    first = false;
                }
            }
        }
        
        if (config_.include_cooccurrence) {
            out << "\t";
            // Write top co-occurring minimizers
            auto it = cooccurrence_patterns_.find(agg_feature.minimizer_hash);
            if (it != cooccurrence_patterns_.end()) {
                bool first = true;
                int count = 0;
                for (const auto& [neighbor_hash, score] : it->second) {
                    if (!first) out << ";";
                    out << std::hex << neighbor_hash << std::dec << ":" 
                        << std::fixed << std::setprecision(3) << score;
                    first = false;
                    if (++count >= config_.top_n_cooccurrences) break;
                }
            }
        }
        
        out << "\n";
    }
    
    out.close();
    
    std::cout << "✓ TSV export completed: " << aggregated_features_.size() << " entries" << std::endl;
    return true;
}

bool FeatureExporter::export_to_hdf5(const std::unordered_map<uint32_t, std::string>& taxon_names) {
#ifdef USE_HDF5
    std::string hdf5_path = config_.output_directory + "/minimizer_features.h5";
    std::cout << "Exporting features to HDF5: " << hdf5_path << std::endl;
    
    try {
        H5::H5File file(hdf5_path, H5F_ACC_TRUNC);
        
        // Create main feature dataset
        hsize_t dims[1] = {aggregated_features_.size()};
        H5::DataSpace dataspace(1, dims);
        
        // Define compound type for features
        H5::CompType feature_type(sizeof(AggregatedMinimizerFeatures));
        feature_type.insertMember("minimizer_hash", HOFFSET(AggregatedMinimizerFeatures, minimizer_hash), H5::PredType::NATIVE_UINT64);
        feature_type.insertMember("total_occurrences", HOFFSET(AggregatedMinimizerFeatures, total_occurrences), H5::PredType::NATIVE_UINT32);
        feature_type.insertMember("mean_position", HOFFSET(AggregatedMinimizerFeatures, mean_position), H5::PredType::NATIVE_FLOAT);
        feature_type.insertMember("position_std_dev", HOFFSET(AggregatedMinimizerFeatures, position_std_dev), H5::PredType::NATIVE_FLOAT);
        feature_type.insertMember("mean_gc_content", HOFFSET(AggregatedMinimizerFeatures, mean_gc_content), H5::PredType::NATIVE_FLOAT);
        feature_type.insertMember("mean_complexity", HOFFSET(AggregatedMinimizerFeatures, mean_complexity), H5::PredType::NATIVE_FLOAT);
        feature_type.insertMember("uniqueness_score", HOFFSET(AggregatedMinimizerFeatures, avg_uniqueness_score), H5::PredType::NATIVE_FLOAT);
        feature_type.insertMember("conservation_level", HOFFSET(AggregatedMinimizerFeatures, conservation_level), H5::PredType::NATIVE_UINT8);
        feature_type.insertMember("contamination_freq", HOFFSET(AggregatedMinimizerFeatures, contamination_frequency), H5::PredType::NATIVE_FLOAT);
        
        // Create and write dataset
        H5::DataSet dataset = file.createDataSet("minimizer_features", feature_type, dataspace);
        dataset.write(aggregated_features_.data(), feature_type);
        
        // Add metadata attributes
        H5::Attribute attr_k = dataset.createAttribute("k_value", H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
        attr_k.write(H5::PredType::NATIVE_INT, &config_.k_value);
        
        // Write taxonomic distribution if included
        if (config_.include_taxonomic_distribution) {
            // Create group for taxonomic data
            H5::Group tax_group = file.createGroup("/taxonomic_distribution");
            
            for (const auto& [minimizer_hash, distribution] : taxonomic_distributions_) {
                std::string dataset_name = std::to_string(minimizer_hash);
                
                std::vector<uint32_t> taxon_ids;
                std::vector<uint32_t> counts;
                
                for (const auto& [taxon_id, count] : distribution) {
                    taxon_ids.push_back(taxon_id);
                    counts.push_back(count);
                }
                
                if (!taxon_ids.empty()) {
                    hsize_t tax_dims[1] = {taxon_ids.size()};
                    H5::DataSpace tax_space(1, tax_dims);
                    
                    H5::DataSet tax_dataset = tax_group.createDataSet(
                        dataset_name + "_taxa", H5::PredType::NATIVE_UINT32, tax_space
                    );
                    tax_dataset.write(taxon_ids.data(), H5::PredType::NATIVE_UINT32);
                    
                    H5::DataSet count_dataset = tax_group.createDataSet(
                        dataset_name + "_counts", H5::PredType::NATIVE_UINT32, tax_space
                    );
                    count_dataset.write(counts.data(), H5::PredType::NATIVE_UINT32);
                }
            }
        }
        
        file.close();
        std::cout << "✓ HDF5 export completed" << std::endl;
        return true;
        
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getCDetailMsg() << std::endl;
        return false;
    }
#else
    std::cerr << "HDF5 support not compiled in. Rebuild with USE_HDF5=1" << std::endl;
    return false;
#endif
}

bool FeatureExporter::export_to_binary() {
    std::string binary_path = config_.output_directory + "/minimizer_features.bin";
    std::cout << "Exporting features to binary format: " << binary_path << std::endl;
    
    std::ofstream out(binary_path, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Cannot create binary file: " << binary_path << std::endl;
        return false;
    }
    
    // Write header
    uint32_t version = 1;
    uint32_t num_features = aggregated_features_.size();
    uint32_t feature_size = sizeof(AggregatedMinimizerFeatures);
    
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&num_features), sizeof(num_features));
    out.write(reinterpret_cast<const char*>(&feature_size), sizeof(feature_size));
    out.write(reinterpret_cast<const char*>(&config_.k_value), sizeof(config_.k_value));
    
    // Write features
    out.write(reinterpret_cast<const char*>(aggregated_features_.data()),
              aggregated_features_.size() * sizeof(AggregatedMinimizerFeatures));
    
    out.close();
    
    std::cout << "✓ Binary export completed: " << (num_features * feature_size / 1024 / 1024) << " MB" << std::endl;
    return true;
}

// Selection methods

void FeatureExporter::select_features_by_occurrence(uint32_t min_occurrences) {
    if (min_occurrences <= 1) return;
    
    std::cout << "Filtering features with < " << min_occurrences << " occurrences..." << std::endl;
    
    auto new_end = std::remove_if(aggregated_features_.begin(), aggregated_features_.end(),
                                  [min_occurrences](const AggregatedMinimizerFeatures& f) {
                                      return f.total_occurrences < min_occurrences;
                                  });
    
    size_t removed = aggregated_features_.end() - new_end;
    aggregated_features_.erase(new_end, aggregated_features_.end());
    
    std::cout << "Removed " << removed << " low-occurrence features" << std::endl;
    total_minimizers_exported_ = aggregated_features_.size();
}

void FeatureExporter::select_features_by_uniqueness(float min_uniqueness) {
    if (min_uniqueness <= 0.0f) return;
    
    std::cout << "Filtering features with uniqueness < " << min_uniqueness << "..." << std::endl;
    
    auto new_end = std::remove_if(aggregated_features_.begin(), aggregated_features_.end(),
                                  [min_uniqueness](const AggregatedMinimizerFeatures& f) {
                                      return f.avg_uniqueness_score < min_uniqueness;
                                  });
    
    size_t removed = aggregated_features_.end() - new_end;
    aggregated_features_.erase(new_end, aggregated_features_.end());
    
    std::cout << "Removed " << removed << " low-uniqueness features" << std::endl;
    total_minimizers_exported_ = aggregated_features_.size();
}

void FeatureExporter::filter_contaminated_features() {
    std::cout << "Filtering contaminated features..." << std::endl;
    
    auto new_end = std::remove_if(aggregated_features_.begin(), aggregated_features_.end(),
                                  [](const AggregatedMinimizerFeatures& f) {
                                      return f.contamination_frequency > 0.5f;
                                  });
    
    size_t removed = aggregated_features_.end() - new_end;
    aggregated_features_.erase(new_end, aggregated_features_.end());
    
    std::cout << "Removed " << removed << " contaminated features" << std::endl;
    total_minimizers_exported_ = aggregated_features_.size();
}

// Private methods

bool FeatureExporter::collect_and_aggregate_features(
    const std::vector<GPUMinimizerHit>& minimizer_hits,
    const MinimizerFeatureExtractor& feature_extractor) {
    
    std::cout << "Collecting and aggregating features..." << std::endl;
    
    // Get unique minimizers
    std::unordered_set<uint64_t> unique_hashes;
    for (const auto& hit : minimizer_hits) {
        unique_hashes.insert(hit.minimizer_hash);
    }
    
    std::vector<uint64_t> unique_minimizers(unique_hashes.begin(), unique_hashes.end());
    std::sort(unique_minimizers.begin(), unique_minimizers.end());
    
    std::cout << "Found " << unique_minimizers.size() << " unique minimizers" << std::endl;
    
    // Allocate aggregated features
    aggregated_features_.resize(unique_minimizers.size());
    
    // GPU aggregation if available
    if (config_.use_gpu_acceleration && minimizer_hits.size() > 10000) {
        // Use GPU kernel for aggregation
        return aggregate_features_gpu(minimizer_hits, unique_minimizers);
    } else {
        // CPU aggregation
        return aggregate_features_cpu(minimizer_hits, unique_minimizers, feature_extractor);
    }
}

bool FeatureExporter::aggregate_features_cpu(
    const std::vector<GPUMinimizerHit>& minimizer_hits,
    const std::vector<uint64_t>& unique_minimizers,
    const MinimizerFeatureExtractor& feature_extractor) {
    
    // Build hash to index map
    std::unordered_map<uint64_t, size_t> hash_to_index;
    for (size_t i = 0; i < unique_minimizers.size(); i++) {
        hash_to_index[unique_minimizers[i]] = i;
    }
    
    // Initialize aggregated features
    for (size_t i = 0; i < unique_minimizers.size(); i++) {
        aggregated_features_[i].minimizer_hash = unique_minimizers[i];
    }
    
    // Aggregate from hits
    for (const auto& hit : minimizer_hits) {
        auto it = hash_to_index.find(hit.minimizer_hash);
        if (it != hash_to_index.end()) {
            size_t idx = it->second;
            AggregatedMinimizerFeatures& agg = aggregated_features_[idx];
            
            agg.total_occurrences++;
            
            // Track taxonomic distribution
            if (config_.include_taxonomic_distribution) {
                taxonomic_distributions_[hit.minimizer_hash][hit.taxon_id]++;
            }
            
            // Extract and aggregate feature values
            MinimizerFeatureValues features;
            features.gc_category = hit.feature_flags & 0x7;
            features.complexity_score = (hit.feature_flags >> 3) & 0x7;
            features.is_contamination = (hit.feature_flags & (1 << 7)) != 0;
            features.conservation_category = (hit.feature_flags >> 8) & 0x7;
            features.uniqueness_score = (float)hit.ml_weight / 65535.0f;
            
            // Update running averages
            float n = agg.total_occurrences;
            agg.mean_gc_content = (agg.mean_gc_content * (n-1) + features.gc_category) / n;
            agg.mean_complexity = (agg.mean_complexity * (n-1) + features.complexity_score) / n;
            
            if (features.is_contamination) {
                agg.contamination_frequency = 
                    (agg.contamination_frequency * (n-1) + 1.0f) / n;
            } else {
                agg.contamination_frequency = 
                    (agg.contamination_frequency * (n-1)) / n;
            }
            
            // Position statistics
            float pos = hit.position;
            float old_mean = agg.mean_position;
            agg.mean_position = (agg.mean_position * (n-1) + pos) / n;
            
            // Online variance calculation
            if (n > 1) {
                float variance = ((n-2) * agg.position_std_dev * agg.position_std_dev + 
                                 (pos - old_mean) * (pos - agg.mean_position)) / (n-1);
                agg.position_std_dev = sqrt(variance);
            }
            
            // First occurrence sets some values
            if (agg.total_occurrences == 1) {
                agg.avg_uniqueness_score = features.uniqueness_score;
                agg.conservation_level = features.conservation_category;
            }
        }
    }
    
    // Calculate unique taxa counts
    for (auto& agg : aggregated_features_) {
        auto it = taxonomic_distributions_.find(agg.minimizer_hash);
        if (it != taxonomic_distributions_.end()) {
            agg.unique_taxa = it->second.size();
        }
    }
    
    total_minimizers_exported_ = aggregated_features_.size();
    return true;
}

bool FeatureExporter::aggregate_features_gpu(
    const std::vector<GPUMinimizerHit>& minimizer_hits,
    const std::vector<uint64_t>& unique_minimizers) {
    
    // Allocate GPU memory
    GPUMinimizerHit* d_hits;
    uint64_t* d_unique;
    AggregatedMinimizerFeatures* d_aggregated;
    
    cudaMalloc(&d_hits, minimizer_hits.size() * sizeof(GPUMinimizerHit));
    cudaMalloc(&d_unique, unique_minimizers.size() * sizeof(uint64_t));
    cudaMalloc(&d_aggregated, unique_minimizers.size() * sizeof(AggregatedMinimizerFeatures));
    
    // Copy data to GPU
    cudaMemcpy(d_hits, minimizer_hits.data(), 
               minimizer_hits.size() * sizeof(GPUMinimizerHit), cudaMemcpyHostToDevice);
    cudaMemcpy(d_unique, unique_minimizers.data(), 
               unique_minimizers.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
    
    // Launch aggregation kernel
    int block_size = 256;
    int grid_size = (unique_minimizers.size() + block_size - 1) / block_size;
    
    aggregate_minimizer_features_kernel<<<grid_size, block_size>>>(
        d_hits, d_unique, d_aggregated, 
        minimizer_hits.size(), unique_minimizers.size()
    );
    
    cudaDeviceSynchronize();
    
    // Copy results back
    cudaMemcpy(aggregated_features_.data(), d_aggregated,
               unique_minimizers.size() * sizeof(AggregatedMinimizerFeatures), 
               cudaMemcpyDeviceToHost);
    
    // Cleanup
    cudaFree(d_hits);
    cudaFree(d_unique);
    cudaFree(d_aggregated);
    
    total_minimizers_exported_ = aggregated_features_.size();
    return true;
}

void FeatureExporter::calculate_cooccurrence_patterns(
    const std::vector<GPUMinimizerHit>& minimizer_hits) {
    
    if (!config_.include_cooccurrence) return;
    
    std::cout << "Calculating co-occurrence patterns..." << std::endl;
    
    // Group hits by genome
    std::unordered_map<uint16_t, std::vector<const GPUMinimizerHit*>> genome_hits;
    for (const auto& hit : minimizer_hits) {
        genome_hits[hit.genome_id].push_back(&hit);
    }
    
    // Calculate co-occurrences within each genome
    for (const auto& [genome_id, hits] : genome_hits) {
        // Sort by position
        std::vector<const GPUMinimizerHit*> sorted_hits = hits;
        std::sort(sorted_hits.begin(), sorted_hits.end(),
                  [](const GPUMinimizerHit* a, const GPUMinimizerHit* b) {
                      return a->position < b->position;
                  });
        
        // Find co-occurrences within window
        for (size_t i = 0; i < sorted_hits.size(); i++) {
            uint64_t center_hash = sorted_hits[i]->minimizer_hash;
            
            // Look ahead within window
            for (size_t j = i + 1; j < sorted_hits.size(); j++) {
                int pos_diff = sorted_hits[j]->position - sorted_hits[i]->position;
                if (pos_diff > config_.cooccurrence_window) break;
                
                uint64_t neighbor_hash = sorted_hits[j]->minimizer_hash;
                float score = 1.0f / (1.0f + pos_diff); // Distance-based score
                
                cooccurrence_patterns_[center_hash][neighbor_hash] += score;
                cooccurrence_patterns_[neighbor_hash][center_hash] += score;
            }
        }
    }
    
    // Normalize and keep only top N
    for (auto& [center_hash, neighbors] : cooccurrence_patterns_) {
        // Sort by score
        std::vector<std::pair<uint64_t, float>> sorted_neighbors(
            neighbors.begin(), neighbors.end()
        );
        std::sort(sorted_neighbors.begin(), sorted_neighbors.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Keep only top N
        neighbors.clear();
        for (size_t i = 0; i < std::min((size_t)config_.top_n_cooccurrences, 
                                        sorted_neighbors.size()); i++) {
            neighbors[sorted_neighbors[i].first] = sorted_neighbors[i].second;
        }
    }
}

void FeatureExporter::generate_export_summary() {
    std::string summary_path = config_.output_directory + "/feature_export_summary.txt";
    std::ofstream out(summary_path);
    if (!out.is_open()) return;
    
    out << "Feature Export Summary\n";
    out << "=====================\n\n";
    
    out << "Export Configuration:\n";
    out << "  K-mer size: " << config_.k_value << "\n";
    out << "  Output formats: ";
    if (config_.export_tsv) out << "TSV ";
    if (config_.export_hdf5) out << "HDF5 ";
    if (config_.export_binary) out << "Binary ";
    out << "\n\n";
    
    out << "Export Statistics:\n";
    out << "  Total unique minimizers: " << total_minimizers_exported_ << "\n";
    out << "  Total occurrences: " << std::accumulate(
        aggregated_features_.begin(), aggregated_features_.end(), 0ULL,
        [](uint64_t sum, const auto& feat) { return sum + feat.total_occurrences; }
    ) << "\n";
    
    // Feature distribution statistics
    out << "\nFeature Distributions:\n";
    
    // GC content distribution
    std::vector<int> gc_dist(8, 0);
    for (const auto& feat : aggregated_features_) {
        int gc_cat = (int)feat.mean_gc_content;
        if (gc_cat >= 0 && gc_cat < 8) gc_dist[gc_cat]++;
    }
    out << "  GC Content Categories:\n";
    for (int i = 0; i < 8; i++) {
        out << "    Category " << i << ": " << gc_dist[i] << "\n";
    }
    
    // Contamination statistics
    int contaminated = 0;
    for (const auto& feat : aggregated_features_) {
        if (feat.contamination_frequency > 0) contaminated++;
    }
    out << "\n  Contaminated minimizers: " << contaminated << " ("
        << std::fixed << std::setprecision(2) 
        << (100.0 * contaminated / total_minimizers_exported_) << "%)\n";
    
    out.close();
}

// ===========================
// Integration Functions
// ===========================

namespace FeatureExportUtils {
    
    bool export_features_for_ml_training(
        const std::string& database_path,
        const std::string& output_directory,
        const FeatureExportConfig& config) {
        
        std::cout << "Loading minimizer data from database: " << database_path << std::endl;
        
        // Load minimizer hits from database
        // This would interface with the database serializer
        std::vector<GPUMinimizerHit> minimizer_hits;
        // ... loading code ...
        
        // Create feature extractor (would be passed from database builder)
        MinimizerFeatureExtractor feature_extractor;
        
        // Create exporter
        FeatureExporter exporter(config);
        
        // Export features
        std::unordered_map<uint32_t, std::string> taxon_names; // Would be loaded
        return exporter.export_training_features(minimizer_hits, feature_extractor, taxon_names);
    }
    
    void print_feature_summary(const std::vector<AggregatedMinimizerFeatures>& features) {
        if (features.empty()) return;
        
        std::cout << "\nFeature Summary:" << std::endl;
        std::cout << "Total features: " << features.size() << std::endl;
        
        // Calculate statistics
        float avg_occurrences = 0;
        float avg_uniqueness = 0;
        int unique_count = 0;
        int rare_count = 0;
        int contaminated_count = 0;
        
        for (const auto& feat : features) {
            avg_occurrences += feat.total_occurrences;
            avg_uniqueness += feat.avg_uniqueness_score;
            
            if (feat.total_occurrences == 1) unique_count++;
            if (feat.total_occurrences < 10) rare_count++;
            if (feat.contamination_frequency > 0) contaminated_count++;
        }
        
        avg_occurrences /= features.size();
        avg_uniqueness /= features.size();
        
        std::cout << "Average occurrences: " << std::fixed << std::setprecision(2) 
                  << avg_occurrences << std::endl;
        std::cout << "Average uniqueness score: " << std::setprecision(4) 
                  << avg_uniqueness << std::endl;
        std::cout << "Singleton minimizers: " << unique_count << " ("
                  << std::setprecision(1) << (100.0 * unique_count / features.size()) << "%)" << std::endl;
        std::cout << "Rare minimizers (<10 occurrences): " << rare_count << " ("
                  << (100.0 * rare_count / features.size()) << "%)" << std::endl;
        std::cout << "Contaminated minimizers: " << contaminated_count << " ("
                  << (100.0 * contaminated_count / features.size()) << "%)" << std::endl;
    }
}