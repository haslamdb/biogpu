// core/gpu_database_builder_core.cu
// Implementation of the main GPU Kraken database builder
// Orchestrates all specialized modules and provides unified API

#include "gpu_database_builder_core.h"
#include "../memory/gpu_memory_manager.h"
#include "../gpu/gpu_database_kernels.h"
#include "../processing/minimizer_feature_extractor.h"
#include "../processing/feature_exporter.h"
#include "../processing/genome_file_processor.h"
#include "../taxonomy/taxonomy_processor.h"
#include "../output/database_serializer.h"
#include "../processing/contamination_detector.h"
#include <cuda_runtime.h>
#include <cstring>  // for memset
#include <iostream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <set>
#include <unordered_map>

// Include uniqueness implementation early to avoid forward declaration issues
#include "../features/working_uniqueness_implementation.cu"

// ===========================
// MinimizerStatistics Method Implementations
// ===========================

void MinimizerStatistics::add_occurrence(uint32_t taxon_id, float position, 
                                        uint8_t gc_category, uint8_t complexity_score,
                                        bool is_clustered, bool is_contamination) {
    occurrence_count++;
    taxon_occurrences.push_back(taxon_id);
    
    // Update running averages
    average_position_in_genome = (average_position_in_genome * (occurrence_count - 1) + position) / occurrence_count;
    
    // Convert gc_category (0-7) to percentage (0-100)
    float gc_percent = (gc_category * 100.0f) / 7.0f;
    gc_content_sum += gc_percent;
    
    // Convert complexity_score (0-7) to normalized value (0-1)
    float complexity_normalized = complexity_score / 7.0f;
    complexity_sum += complexity_normalized;
    
    if (is_clustered) clustered_count++;
    if (is_contamination) contamination_risk_count++;
}

void MinimizerStatistics::add_neighbor(uint64_t neighbor_hash) {
    neighbor_minimizers.insert(neighbor_hash);
}

float MinimizerStatistics::get_average_gc_content() const {
    return occurrence_count > 0 ? gc_content_sum / occurrence_count : 0.0f;
}

float MinimizerStatistics::get_average_complexity() const {
    return occurrence_count > 0 ? complexity_sum / occurrence_count : 0.0f;
}

float MinimizerStatistics::get_clustering_ratio() const {
    return occurrence_count > 0 ? static_cast<float>(clustered_count) / occurrence_count : 0.0f;
}

float MinimizerStatistics::get_contamination_risk_ratio() const {
    return occurrence_count > 0 ? static_cast<float>(contamination_risk_count) / occurrence_count : 0.0f;
}

size_t MinimizerStatistics::get_taxonomic_diversity() const {
    std::unordered_set<uint32_t> unique_taxons(taxon_occurrences.begin(), taxon_occurrences.end());
    return unique_taxons.size();
}


// ===========================
// GPUKrakenDatabaseBuilder Implementation
// ===========================

GPUKrakenDatabaseBuilder::GPUKrakenDatabaseBuilder(const std::string& output_dir,
                                                 const DatabaseBuildConfig& config)
    : config_(config), output_directory_(output_dir), initialized_(false),
      cuda_context_ready_(false), modules_initialized_(false) {
    
    std::cout << "Initializing GPU Kraken Database Builder (Modular Architecture)" << std::endl;
    std::cout << "Output directory: " << output_directory_ << std::endl;
    
    // Validate configuration
    if (!validate_configuration()) {
        std::cerr << "Invalid configuration provided" << std::endl;
        return;
    }
    
    // Initialize minimizer parameters from config
    minimizer_params_.k = config_.k_value;
    minimizer_params_.ell = config_.ell_value;
    minimizer_params_.spaces = config_.spaces_value;
    minimizer_params_.xor_mask = 0x3c8bfbb395c60474ULL;
    
    // Initialize statistics
    memset(&build_stats_, 0, sizeof(build_stats_));
    enhanced_stats_ = EnhancedBuildStats();
    
    // Create output directory
    try {
        std::filesystem::create_directories(output_directory_);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not create output directory: " << e.what() << std::endl;
    }
    
    std::cout << "Core builder initialized with parameters: k=" << config_.k_value 
              << ", ell=" << config_.ell_value << ", spaces=" << config_.spaces_value << std::endl;
    
    initialized_ = true;
}

GPUKrakenDatabaseBuilder::~GPUKrakenDatabaseBuilder() {
    try {
        cleanup_on_failure();
    } catch (...) {
        // Suppress exceptions in destructor
    }
}

bool GPUKrakenDatabaseBuilder::build_database_from_genomes(
    const std::string& genome_library_path,
    const std::string& taxonomy_path) {
    
    std::cout << "\n=== BUILDING DATABASE FROM GENOME DIRECTORY ===" << std::endl;
    std::cout << "Genome library: " << genome_library_path << std::endl;
    std::cout << "Taxonomy: " << (taxonomy_path.empty() ? "None" : taxonomy_path) << std::endl;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Initialize all modules
    if (!initialize_all_modules()) {
        std::cerr << "Failed to initialize modules" << std::endl;
        return false;
    }
    
    // Load and validate inputs
    if (!load_and_validate_inputs(genome_library_path, taxonomy_path)) {
        std::cerr << "Failed to load and validate inputs" << std::endl;
        return false;
    }
    
    // Execute the main build pipeline
    if (!execute_build_pipeline()) {
        std::cerr << "Build pipeline execution failed" << std::endl;
        return false;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
    
    std::cout << "\n=== DATABASE BUILD COMPLETE ===" << std::endl;
    std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
    
    // Print final statistics
    aggregate_module_statistics();
    print_build_progress();
    
    if (config_.validate_output) {
        validate_database_output();
    }
    
    return true;
}

bool GPUKrakenDatabaseBuilder::build_database_from_concatenated_fna(
    const std::string& fna_file_path,
    const std::string& taxonomy_nodes_path,
    const std::string& taxonomy_names_path,
    const std::string& compact_taxonomy_path) {
    
    std::cout << "\n=== BUILDING DATABASE FROM CONCATENATED FNA ===" << std::endl;
    std::cout << "FNA file: " << fna_file_path << std::endl;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Initialize modules
    if (!initialize_all_modules()) {
        std::cerr << "Failed to initialize modules" << std::endl;
        return false;
    }
    
    // Load taxonomy data first
    if (!taxonomy_nodes_path.empty() && !taxonomy_names_path.empty()) {
        std::cout << "Loading NCBI taxonomy..." << std::endl;
        if (!taxonomy_processor_->load_ncbi_taxonomy(taxonomy_nodes_path, taxonomy_names_path)) {
            std::cerr << "Warning: Failed to load NCBI taxonomy" << std::endl;
        }
    } else if (!compact_taxonomy_path.empty()) {
        std::cout << "Loading compact taxonomy..." << std::endl;
        if (!taxonomy_processor_->load_from_compact_file(compact_taxonomy_path)) {
            std::cerr << "Warning: Failed to load compact taxonomy" << std::endl;
        }
    }
    
    // Process FNA file using the file processor
    std::cout << "Processing concatenated FNA file..." << std::endl;
    ConcatenatedFnaProcessor fna_processor(fna_file_path, config_.file_config);
    
    // Clear existing data
    genome_files_.clear();
    genome_taxon_ids_.clear();
    taxon_names_.clear();
    
    // Process FNA and get temporary genome files
    std::string temp_dir = output_directory_ + "/temp_genomes";
    if (!fna_processor.process_fna_file(genome_files_, genome_taxon_ids_, taxon_names_, temp_dir)) {
        std::cerr << "Failed to process FNA file" << std::endl;
        return false;
    }
    
    // Get species tracking data
    species_tracking_ = fna_processor.get_species_data();
    
    // Update enhanced statistics
    enhanced_stats_.species_represented = species_tracking_.total_species();
    enhanced_stats_.total_sequences = genome_files_.size();
    build_stats_.total_sequences = genome_files_.size();
    
    std::cout << "FNA processing complete:" << std::endl;
    std::cout << "  Genomes: " << genome_files_.size() << std::endl;
    std::cout << "  Species: " << species_tracking_.total_species() << std::endl;
    
    // Continue with standard pipeline
    if (!process_genomes_with_gpu()) {
        std::cerr << "Failed to process genomes on GPU" << std::endl;
        return false;
    }
    
    // Enhanced processing if enabled
    if (config_.enable_phylogenetic_analysis && taxonomy_processor_->is_loaded()) {
        if (!enhance_with_phylogenetic_data()) {
            std::cerr << "Failed to enhance with phylogenetic data" << std::endl;
            return false;
        }
    }
    
    // Save outputs
    if (!save_database_outputs()) {
        std::cerr << "Failed to save database outputs" << std::endl;
        return false;
    }
    
    // Cleanup temporary files
    try {
        std::filesystem::remove_all(temp_dir);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not clean up temporary files: " << e.what() << std::endl;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
    
    std::cout << "\n=== ENHANCED DATABASE BUILD COMPLETE ===" << std::endl;
    std::cout << "Total build time: " << total_duration.count() << " seconds" << std::endl;
    
    if (config_.enable_phylogenetic_analysis) {
        enhanced_stats_.print_enhanced_stats();
    }
    
    return true;
}

bool GPUKrakenDatabaseBuilder::build_database_from_streaming_fna(
    const std::string& fna_file_path,
    const std::string& taxonomy_path) {
    
    std::cout << "\n=== BUILDING DATABASE FROM STREAMING FNA ===" << std::endl;
    
    // Initialize modules
    if (!initialize_all_modules()) {
        return false;
    }
    
    // Load taxonomy if provided
    if (!taxonomy_path.empty()) {
        load_taxonomy_data(taxonomy_path);
    }
    
    // Use streaming processor
    std::string temp_dir = output_directory_ + "/temp_streaming";
    StreamingFnaProcessor streaming_processor(fna_file_path, temp_dir, 
                                            config_.memory_config.sequence_batch_size);
    
    // Process in streaming batches
    std::vector<std::string> batch_files;
    std::vector<uint32_t> batch_taxons;
    int batch_number = 0;
    
    while (streaming_processor.process_next_batch(batch_files, batch_taxons)) {
        batch_number++;
        std::cout << "\nProcessing streaming batch " << batch_number 
                  << " (" << batch_files.size() << " genomes)" << std::endl;
        
        // Update internal data for this batch
        genome_files_ = batch_files;
        genome_taxon_ids_ = batch_taxons;
        
        // Process this batch
        if (!process_genomes_with_gpu()) {
            std::cerr << "Failed to process batch " << batch_number << std::endl;
            return false;
        }
        
        // Update statistics
        build_stats_.total_sequences += batch_files.size();
        
        // Clean up batch files
        for (const auto& file : batch_files) {
            std::filesystem::remove(file);
        }
    }
    
    // Final save
    if (!save_database_outputs()) {
        return false;
    }
    
    // Cleanup
    std::filesystem::remove_all(temp_dir);
    
    std::cout << "\nStreaming processing complete:" << std::endl;
    std::cout << "  Total batches: " << batch_number << std::endl;
    std::cout << "  Total genomes: " << streaming_processor.get_total_genomes() << std::endl;
    
    return true;
}

// Configuration methods
bool GPUKrakenDatabaseBuilder::configure_memory_settings(const MemoryConfig& memory_config) {
    config_.memory_config = memory_config;
    
    if (memory_manager_ && memory_manager_->initialize()) {
        return memory_manager_->configure_auto_scaling(
            memory_config.auto_scale_enabled,
            memory_config.max_memory_fraction
        );
    }
    
    return true;
}

bool GPUKrakenDatabaseBuilder::configure_phylogenetic_analysis(bool enable) {
    config_.enable_phylogenetic_analysis = enable;
    
    if (enable) {
        config_.save_enhanced_format = true;
        std::cout << "Phylogenetic analysis enabled - enhanced format will be saved" << std::endl;
    }
    
    return true;
}

void GPUKrakenDatabaseBuilder::set_minimizer_parameters(uint32_t k, uint32_t ell, uint32_t spaces) {
    config_.k_value = k;
    config_.ell_value = ell;
    config_.spaces_value = spaces;
    
    minimizer_params_.k = k;
    minimizer_params_.ell = ell;
    minimizer_params_.spaces = spaces;
    
    std::cout << "Updated minimizer parameters: k=" << k << ", ell=" << ell << ", spaces=" << spaces << std::endl;
}

void GPUKrakenDatabaseBuilder::set_subsampling_rate(double rate) {
    if (rate > 0.0 && rate <= 1.0) {
        config_.subsampling_rate = rate;
        std::cout << "Set subsampling rate to " << rate << std::endl;
    } else {
        std::cerr << "Invalid subsampling rate: " << rate << " (must be 0.0 < rate <= 1.0)" << std::endl;
    }
}

// Private implementation methods

bool GPUKrakenDatabaseBuilder::initialize_all_modules() {
    std::cout << "Initializing all modules..." << std::endl;
    
    record_timing_checkpoint("module_initialization_start");
    
    // Initialize in dependency order
    if (!initialize_cuda_context()) {
        std::cerr << "Failed to initialize CUDA context" << std::endl;
        return false;
    }
    
    if (!initialize_memory_manager()) {
        std::cerr << "Failed to initialize memory manager" << std::endl;
        return false;
    }
    
    if (!initialize_file_processor()) {
        std::cerr << "Failed to initialize file processor" << std::endl;
        return false;
    }
    
    if (!initialize_taxonomy_processor()) {
        std::cerr << "Failed to initialize taxonomy processor" << std::endl;
        return false;
    }
    
    if (!initialize_serializers()) {
        std::cerr << "Failed to initialize serializers" << std::endl;
        return false;
    }
    
    // Initialize feature extractor
    try {
        size_t max_minimizers = config_.memory_config.minimizer_capacity;
        size_t max_genomes = config_.memory_config.sequence_batch_size * 100; // Estimate
        feature_extractor_ = std::make_unique<MinimizerFeatureExtractor>(max_minimizers, max_genomes);
        std::cout << "Feature extractor initialized with capacity for " 
                  << max_minimizers << " minimizers" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize feature extractor: " << e.what() << std::endl;
        // Non-fatal - continue without feature extraction
        feature_extractor_.reset();
    }
    
    modules_initialized_ = true;
    record_timing_checkpoint("module_initialization_end");
    
    std::cout << "✓ All modules initialized successfully" << std::endl;
    
    return true;
}

bool GPUKrakenDatabaseBuilder::initialize_cuda_context() {
    std::cout << "Initializing CUDA context..." << std::endl;
    
    // Force CUDA runtime initialization
    cudaError_t init_status = cudaFree(0);
    if (init_status != cudaSuccess) {
        std::cerr << "CUDA runtime initialization failed: " << cudaGetErrorString(init_status) << std::endl;
        return false;
    }
    
    // Check device count and capabilities
    int device_count;
    cudaError_t cuda_status = cudaGetDeviceCount(&device_count);
    if (cuda_status != cudaSuccess || device_count == 0) {
        std::cerr << "No CUDA devices available" << std::endl;
        return false;
    }
    
    cuda_status = cudaSetDevice(0);
    if (cuda_status != cudaSuccess) {
        std::cerr << "Failed to set CUDA device" << std::endl;
        return false;
    }
    
    // Query GPU memory
    size_t free_mem, total_mem;
    cuda_status = cudaMemGetInfo(&free_mem, &total_mem);
    if (cuda_status != cudaSuccess) {
        std::cerr << "Failed to query GPU memory" << std::endl;
        return false;
    }
    
    std::cout << "GPU Memory: " << (free_mem / 1024 / 1024) << " MB free / " 
              << (total_mem / 1024 / 1024) << " MB total" << std::endl;
    
    cuda_context_ready_ = true;
    return true;
}

bool GPUKrakenDatabaseBuilder::initialize_memory_manager() {
    std::cout << "Initializing memory manager..." << std::endl;
    
    memory_manager_ = std::make_unique<GPUMemoryManager>(config_.memory_config);
    
    if (!memory_manager_->initialize()) {
        std::cerr << "Memory manager initialization failed" << std::endl;
        return false;
    }
    
    std::cout << "✓ Memory manager initialized" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::initialize_file_processor() {
    std::cout << "Initializing file processor..." << std::endl;
    
    file_processor_ = std::make_unique<GenomeFileProcessor>(config_.file_config);
    
    std::cout << "✓ File processor initialized" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::initialize_taxonomy_processor() {
    std::cout << "Initializing taxonomy processor..." << std::endl;
    
    taxonomy_processor_ = std::make_unique<EnhancedNCBITaxonomyProcessor>();
    
    std::cout << "✓ Taxonomy processor initialized" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::initialize_serializers() {
    std::cout << "Initializing database serializers..." << std::endl;
    
    standard_serializer_ = std::make_unique<StandardDatabaseSerializer>(output_directory_);
    enhanced_serializer_ = std::make_unique<EnhancedDatabaseSerializer>(output_directory_);
    
    std::cout << "✓ Database serializers initialized" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::execute_build_pipeline() {
    std::cout << "\n=== EXECUTING BUILD PIPELINE ===" << std::endl;
    
    // Step 1: Process genomes on GPU
    update_build_progress("GPU Processing", 0.0);
    if (!process_genomes_with_gpu()) {
        return false;
    }
    update_build_progress("GPU Processing", 100.0);
    
    // Step 2: Enhanced phylogenetic analysis (if enabled)
    if (config_.enable_phylogenetic_analysis) {
        update_build_progress("Phylogenetic Analysis", 0.0);
        if (!enhance_with_phylogenetic_data()) {
            return false;
        }
        update_build_progress("Phylogenetic Analysis", 100.0);
    }
    
    // Step 3: Save database outputs
    update_build_progress("Database Serialization", 0.0);
    if (!save_database_outputs()) {
        return false;
    }
    update_build_progress("Database Serialization", 100.0);
    
    std::cout << "✓ Build pipeline completed successfully" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::load_and_validate_inputs(const std::string& input_path, const std::string& taxonomy_path) {
    std::cout << "Loading and validating inputs..." << std::endl;
    
    // Load genome files
    if (!load_genome_files(input_path)) {
        return false;
    }
    
    // Load taxonomy if provided
    if (!taxonomy_path.empty()) {
        if (!load_taxonomy_data(taxonomy_path)) {
            std::cerr << "Warning: Failed to load taxonomy data" << std::endl;
        }
    }
    
    // Validate consistency
    if (!validate_inter_module_consistency()) {
        return false;
    }
    
    return true;
}

bool GPUKrakenDatabaseBuilder::load_genome_files(const std::string& library_path) {
    std::cout << "Loading genome files from: " << library_path << std::endl;
    
    // Use file processor to find and validate files
    genome_files_ = file_processor_->find_genome_files(library_path);
    
    if (genome_files_.empty()) {
        std::cerr << "No genome files found in " << library_path << std::endl;
        return false;
    }
    
    // Extract taxon IDs from filenames
    genome_taxon_ids_.clear();
    genome_taxon_ids_.reserve(genome_files_.size());
    
    for (const auto& file_path : genome_files_) {
        uint32_t taxon_id = file_processor_->extract_taxon_from_filename(file_path);
        genome_taxon_ids_.push_back(taxon_id);
        
        // Add to taxon names if not present
        if (taxon_names_.find(taxon_id) == taxon_names_.end()) {
            std::filesystem::path p(file_path);
            std::string filename = p.stem().string();
            taxon_names_[taxon_id] = filename.empty() ? ("taxon_" + std::to_string(taxon_id)) : filename;
        }
    }
    
    build_stats_.total_sequences = genome_files_.size();
    
    std::cout << "✓ Loaded " << genome_files_.size() << " genome files" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::load_taxonomy_data(const std::string& taxonomy_path) {
    std::cout << "Loading taxonomy data from: " << taxonomy_path << std::endl;
    
    // Try to load with enhanced processor first
    bool loaded = false;
    
    // Check if it's a directory with nodes.dmp and names.dmp
    std::filesystem::path tax_path(taxonomy_path);
    
    if (std::filesystem::is_directory(tax_path)) {
        std::string nodes_file = tax_path / "nodes.dmp";
        std::string names_file = tax_path / "names.dmp";
        
        if (std::filesystem::exists(nodes_file) && std::filesystem::exists(names_file)) {
            loaded = taxonomy_processor_->load_ncbi_taxonomy(nodes_file, names_file);
        }
    } else if (std::filesystem::exists(taxonomy_path)) {
        // Try to load as compact taxonomy file
        loaded = taxonomy_processor_->load_from_compact_file(taxonomy_path);
    }
    
    if (loaded) {
        std::cout << "✓ Taxonomy data loaded successfully" << std::endl;
        
        // Extract basic taxonomy information for compatibility
        if (taxon_parents_.empty()) {
            // Build basic parent relationships from loaded taxonomy
            const auto& parent_lookup = taxonomy_processor_->get_parent_lookup();
            const auto& name_lookup = taxonomy_processor_->get_name_lookup();
            
            taxon_parents_ = parent_lookup;
            
            // Merge with existing taxon names
            for (const auto& [taxon_id, name] : name_lookup) {
                if (taxon_names_.find(taxon_id) == taxon_names_.end()) {
                    taxon_names_[taxon_id] = name;
                }
            }
        }
    } else {
        std::cerr << "Warning: Could not load taxonomy data from " << taxonomy_path << std::endl;
    }
    
    return loaded;
}

bool GPUKrakenDatabaseBuilder::process_genomes_with_gpu() {
    std::cout << "Processing genomes on GPU..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Coordinate memory allocation
    if (!coordinate_memory_allocation()) {
        return false;
    }
    
    // Initialize accumulation variables for batch processing
    accumulated_genome_info_.clear();
    accumulated_sequence_length_ = 0;
    accumulated_sequences_.clear();
    
    // Process in batches
    size_t batch_size = config_.memory_config.sequence_batch_size;
    
    for (size_t batch_start = 0; batch_start < genome_files_.size(); batch_start += batch_size) {
        size_t buffer_limit = memory_manager_->get_sequence_buffer_size() * 0.9;  // 90% of buffer size
        size_t batch_end = std::min(batch_start + batch_size, genome_files_.size());
        
        std::cout << "Processing batch " << (batch_start / batch_size + 1) 
                  << ": genomes " << batch_start << "-" << batch_end << std::endl;
        
        // Load sequences for this batch
        std::vector<std::string> batch_sequences;
        std::vector<uint32_t> batch_taxon_ids;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            auto sequences = file_processor_->load_sequences_from_fasta(genome_files_[i]);
            if (sequences.size() > 1) {
                std::cout << "  File " << i << " contains " << sequences.size() << " sequences" << std::endl;
            }
            for (const auto& seq : sequences) {
                batch_sequences.push_back(seq);
                batch_taxon_ids.push_back(genome_taxon_ids_[i]);
                build_stats_.total_bases += seq.length();
            }
        }
        
        std::cout << "  Total sequences in batch: " << batch_sequences.size() << std::endl;
        
        // Process this batch
        if (!process_sequence_batch(batch_sequences, batch_taxon_ids)) {
            std::cerr << "Failed to process batch starting at " << batch_start << std::endl;
            return false;
        }
        
        if (config_.enable_progress_reporting) {
            double progress = (double)batch_end / genome_files_.size() * 100.0;
            update_build_progress("Genome Processing", progress);
        }
        
        // Check if buffer is getting full
        if (accumulated_sequence_length_ > buffer_limit) {
            std::cout << "Buffer 90% full at " << accumulated_sequence_length_ << " bytes, processing accumulated sequences..." << std::endl;
            if (!process_accumulated_sequences()) {
                std::cerr << "Failed to process accumulated sequences" << std::endl;
                return false;
            }
            std::cout << "✓ Intermediate processing completed, continuing with more genomes..." << std::endl;
            std::cout << "Cleared accumulation buffers for next batch" << std::endl;
        }
    }
    
    // Now process all accumulated data on GPU
    std::cout << "\nProcessing all accumulated genomes on GPU..." << std::endl;
    std::cout << "Total genomes: " << accumulated_genome_info_.size() << std::endl;
    std::cout << "Total sequence data: " << accumulated_sequences_.length() << " bytes" << std::endl;
    std::cout << "Accumulated sequence length: " << accumulated_sequence_length_ << " bytes" << std::endl;
    
    // Debug: Print first 3 genome info entries
    std::cout << "\nFirst 3 genome info entries:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(3), accumulated_genome_info_.size()); i++) {
        const auto& info = accumulated_genome_info_[i];
        std::cout << "  Genome " << i << ": id=" << info.genome_id 
                  << ", offset=" << info.sequence_offset 
                  << ", length=" << info.sequence_length 
                  << ", taxon=" << info.taxon_id << std::endl;
    }
    
    // Process any remaining accumulated sequences
    if (!accumulated_genome_info_.empty()) {
        std::cout << "\nProcessing final " << accumulated_genome_info_.size() 
                  << " sequences (" << accumulated_sequence_length_ << " bytes)" << std::endl;
        if (!process_accumulated_sequences()) {
            std::cerr << "Failed to process final batch" << std::endl;
            return false;
        }
    }
    
    // Run second pass of feature extraction after all batches are processed
    if (feature_extractor_ && !all_lca_candidates_.empty()) {
        std::cout << "Running second pass of feature extraction..." << std::endl;
        
        // Prepare minimizer hits for second pass
        // Note: This is a simplified approach - in production, we'd maintain GPU memory across batches
        std::vector<GPUMinimizerHit> all_minimizer_hits;
        // Convert LCA candidates back to minimizer hits (simplified)
        for (const auto& candidate : all_lca_candidates_) {
            GPUMinimizerHit hit;
            hit.minimizer_hash = candidate.minimizer_hash;
            hit.taxon_id = static_cast<uint16_t>(candidate.lca_taxon);
            hit.ml_weight = MinimizerFlags::float_to_ml_weight(candidate.uniqueness_score);
            // Other fields would be preserved in a real implementation
            all_minimizer_hits.push_back(hit);
        }
        
        // Allocate GPU memory for all hits
        GPUMinimizerHit* d_all_hits = nullptr;
        cudaMalloc(&d_all_hits, all_minimizer_hits.size() * sizeof(GPUMinimizerHit));
        cudaMemcpy(d_all_hits, all_minimizer_hits.data(), 
                   all_minimizer_hits.size() * sizeof(GPUMinimizerHit), cudaMemcpyHostToDevice);
        
        // Run second pass
        if (feature_extractor_->process_second_pass(d_all_hits, all_minimizer_hits.size())) {
            // Copy updated hits back
            cudaMemcpy(all_minimizer_hits.data(), d_all_hits, 
                       all_minimizer_hits.size() * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
            
            // Update ML weights in candidates
            for (size_t i = 0; i < all_lca_candidates_.size() && i < all_minimizer_hits.size(); i++) {
                all_lca_candidates_[i].uniqueness_score = 
                    MinimizerFlags::ml_weight_to_float(all_minimizer_hits[i].ml_weight);
            }
            
            std::cout << "✓ Feature extraction second pass completed" << std::endl;
        }
        
        cudaFree(d_all_hits);
        
        // Export feature statistics if requested
        if (config_.enable_debug_mode && !config_.debug_output_dir.empty()) {
            feature_extractor_->export_feature_statistics(
                config_.debug_output_dir + "/minimizer_features.tsv");
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    build_stats_.sequence_processing_time = 
        std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "✓ GPU processing completed" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::process_sequence_batch(
    const std::vector<std::string>& sequences,
    const std::vector<uint32_t>& taxon_ids) {
    
    // This method coordinates with the GPU kernels module
    // Implementation would call the kernels from gpu_database_kernels.cu
    
    if (sequences.empty()) return true;
    
    // Use memory manager to get GPU buffers
    char* d_sequence_data = memory_manager_->get_sequence_buffer();
    GPUGenomeInfo* d_genome_info = memory_manager_->get_genome_info_buffer();
    GPUMinimizerHit* d_minimizer_hits = memory_manager_->get_minimizer_buffer();
    LCACandidate* d_lca_candidates = memory_manager_->get_candidate_buffer();
    
    if (!d_sequence_data || !d_genome_info || !d_minimizer_hits || !d_lca_candidates) {
        std::cerr << "GPU memory buffers not allocated" << std::endl;
        return false;
    }
    
    // Prepare batch data for GPU processing
    GPUBatchData batch_data;
    batch_data.d_sequence_data = d_sequence_data;
    batch_data.d_genome_info = d_genome_info;
    batch_data.d_minimizer_hits = d_minimizer_hits;
    batch_data.d_global_counter = memory_manager_->get_global_counter();
    batch_data.d_lca_candidates = d_lca_candidates;
    batch_data.max_genomes = sequences.size();
    batch_data.max_minimizers = memory_manager_->get_minimizer_capacity();
    batch_data.sequence_buffer_size = memory_manager_->get_sequence_buffer_size();
    
    // Prepare genome information and accumulate across batches
    std::vector<uint32_t> genome_lengths;
    size_t sequence_offset = accumulated_sequence_length_;
    
    for (size_t i = 0; i < sequences.size(); i++) {
        GPUGenomeInfo info;
        info.genome_id = accumulated_genome_info_.size();  // Use accumulated count as genome_id
        info.sequence_offset = sequence_offset;  // Offset in accumulated buffer
        info.sequence_length = sequences[i].length();
        info.minimizer_count = 0;
        info.taxon_id = taxon_ids[i];
        
        if (i < 3 || i >= sequences.size() - 3) { // Debug first and last few
            std::cout << "  Genome " << info.genome_id << ": offset=" << info.sequence_offset 
                      << ", length=" << info.sequence_length 
                      << ", taxon=" << info.taxon_id << std::endl;
        }
        
        accumulated_genome_info_.push_back(info);
        genome_lengths.push_back(sequences[i].length());
        sequence_offset += sequences[i].length() + 1; // +1 for null terminator
        accumulated_sequence_length_ += sequences[i].length() + 1;
    }
    
    // Append sequences to accumulated buffer
    for (const auto& seq : sequences) {
        accumulated_sequences_ += seq + '\0';
    }
    
    // Just accumulate data in this method, GPU operations moved to after batch loop
    std::cout << "  Batch accumulated. Total sequences so far: " << accumulated_genome_info_.size() << std::endl;
    std::cout << "  Total accumulated length: " << accumulated_sequences_.length() << " bytes" << std::endl;
    std::cout << "  accumulated_sequence_length_: " << accumulated_sequence_length_ << " bytes" << std::endl;
    
    // Check if we're exceeding buffer size
    if (accumulated_sequences_.length() > memory_manager_->get_sequence_buffer_size()) {
        std::cerr << "ERROR: Accumulated sequences exceed buffer size!" << std::endl;
        return false;
    }
    
    // Return early - GPU processing will happen after all batches are accumulated
    return true;
    
    // Reset global counter before kernel launch
    uint32_t zero = 0;
    cudaMemcpy(batch_data.d_global_counter, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    // Launch minimizer extraction kernel
    uint32_t total_hits = 0;
    if (!launch_minimizer_extraction_kernel(batch_data, minimizer_params_, &total_hits)) {
        std::cerr << "Failed to launch minimizer extraction kernel" << std::endl;
        return false;
    }
    
    std::cout << "  Extracted " << total_hits << " minimizers from batch" << std::endl;
    
    // Copy minimizer hits back to host
    std::vector<GPUMinimizerHit> minimizer_hits(total_hits);
    if (total_hits > 0) {
        cudaMemcpy(minimizer_hits.data(), d_minimizer_hits, 
                   total_hits * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
        
        // Collect statistics for the extracted minimizers
        collect_minimizer_statistics(minimizer_hits, genome_lengths);
        
        // Use two-pass feature extraction if available
        if (feature_extractor_) {
            // First pass: collect statistics
            if (!feature_extractor_->process_first_pass(d_minimizer_hits, total_hits, taxon_ids)) {
                std::cerr << "Failed to process first pass of feature extraction" << std::endl;
            }
            
            // Note: Second pass will be called after all batches are processed
        }
        
        // Convert minimizer hits to LCA candidates (simplified for now)
        for (const auto& hit : minimizer_hits) {
            LCACandidate candidate;
            candidate.minimizer_hash = hit.minimizer_hash;
            candidate.lca_taxon = hit.taxon_id;
            candidate.genome_count = 1;
            candidate.uniqueness_score = MinimizerFlags::ml_weight_to_float(hit.ml_weight);
            
            all_lca_candidates_.push_back(candidate);
        }
    }
    
    build_stats_.valid_minimizers_extracted += total_hits;
    build_stats_.lca_assignments += total_hits;
    
    return true;
}

void GPUKrakenDatabaseBuilder::collect_minimizer_statistics(
    const std::vector<GPUMinimizerHit>& minimizer_hits,
    const std::vector<uint32_t>& genome_lengths) {
    
    // Sort hits by genome_id and position to facilitate neighbor detection
    std::vector<GPUMinimizerHit> sorted_hits = minimizer_hits;
    std::sort(sorted_hits.begin(), sorted_hits.end(), 
              [](const GPUMinimizerHit& a, const GPUMinimizerHit& b) {
                  if (a.genome_id != b.genome_id) return a.genome_id < b.genome_id;
                  return a.position < b.position;
              });
    
    // Process each hit and update statistics
    for (size_t i = 0; i < sorted_hits.size(); i++) {
        const GPUMinimizerHit& hit = sorted_hits[i];
        
        // Calculate relative position in genome (0.0 to 1.0)
        float relative_position = 0.5f; // Default to middle
        if (hit.genome_id < genome_lengths.size() && genome_lengths[hit.genome_id] > 0) {
            relative_position = static_cast<float>(hit.position) / genome_lengths[hit.genome_id];
        }
        
        // Extract features from feature_flags
        uint8_t gc_category = MinimizerFlags::get_gc_content_category(hit.feature_flags);
        uint8_t complexity_score = MinimizerFlags::get_complexity_score(hit.feature_flags);
        bool is_clustered = MinimizerFlags::has_position_bias(hit.feature_flags);
        bool is_contamination = MinimizerFlags::has_contamination_risk(hit.feature_flags);
        
        // Update or create statistics for this minimizer
        auto& stats = minimizer_stats_[hit.minimizer_hash];
        stats.add_occurrence(hit.taxon_id, relative_position, gc_category, 
                           complexity_score, is_clustered, is_contamination);
        
        // Track neighboring minimizers (within same genome)
        const int neighbor_window = 5; // Look for neighbors within +/- 5 positions
        
        // Look backward for neighbors
        for (int j = i - 1; j >= 0 && sorted_hits[j].genome_id == hit.genome_id; j--) {
            if (i - j > neighbor_window) break;
            stats.add_neighbor(sorted_hits[j].minimizer_hash);
        }
        
        // Look forward for neighbors
        for (size_t j = i + 1; j < sorted_hits.size() && sorted_hits[j].genome_id == hit.genome_id; j++) {
            if (j - i > neighbor_window) break;
            stats.add_neighbor(sorted_hits[j].minimizer_hash);
        }
    }
    
    // Update build statistics
    build_stats_.unique_minimizers = minimizer_stats_.size();
}

bool GPUKrakenDatabaseBuilder::process_accumulated_sequences() {
    std::cout << "Processing " << accumulated_genome_info_.size() << " accumulated sequences..." << std::endl;
    
    // Debug: Show first few genome info entries
    for (size_t i = 0; i < std::min(size_t(3), accumulated_genome_info_.size()); i++) {
        std::cout << "  Genome " << i << ": offset=" << accumulated_genome_info_[i].sequence_offset 
                  << ", length=" << accumulated_genome_info_[i].sequence_length 
                  << ", taxon=" << accumulated_genome_info_[i].taxon_id << std::endl;
    }
    
    // Copy GPU buffers
    char* d_sequence_data = memory_manager_->get_sequence_buffer();
    GPUGenomeInfo* d_genome_info = memory_manager_->get_genome_info_buffer();
    GPUMinimizerHit* d_minimizer_hits = memory_manager_->get_minimizer_buffer();
    LCACandidate* d_lca_candidates = memory_manager_->get_candidate_buffer();
    
    // Copy accumulated data to GPU
    cudaMemcpy(d_sequence_data, accumulated_sequences_.data(), accumulated_sequences_.length(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_genome_info, accumulated_genome_info_.data(), 
               accumulated_genome_info_.size() * sizeof(GPUGenomeInfo), cudaMemcpyHostToDevice);
    
    // Prepare batch data for kernel
    GPUBatchData batch_data;
    batch_data.d_sequence_data = d_sequence_data;
    batch_data.d_genome_info = d_genome_info;
    batch_data.d_minimizer_hits = d_minimizer_hits;
    batch_data.d_lca_candidates = d_lca_candidates;
    batch_data.sequence_buffer_size = accumulated_sequences_.length();
    batch_data.max_genomes = accumulated_genome_info_.size();
    batch_data.max_minimizers = memory_manager_->get_minimizer_capacity();
    
    // Get hit counter from memory manager
    uint32_t* d_hit_counter = memory_manager_->get_global_counter();
    batch_data.d_global_counter = d_hit_counter;
    
    // Launch minimizer extraction kernel with improved work distribution
    uint32_t total_hits_extracted = 0;
    if (!launch_improved_minimizer_kernel(
            batch_data,
            minimizer_params_,
            0,  // min_clear_hash_value
            config_.toggle_mask,
            &total_hits_extracted)) {
        std::cerr << "Minimizer extraction kernel failed" << std::endl;
        return false;
    }
    
    std::cout << "Minimizers extracted: " << total_hits_extracted << std::endl;
    build_stats_.valid_minimizers_extracted += total_hits_extracted;
    
    if (total_hits_extracted > 0) {
        // Process minimizers with feature extraction
        if (feature_extractor_) {
            std::cout << "Running feature extraction on minimizers..." << std::endl;
            
            // Create taxon IDs vector from genome info
            std::vector<uint32_t> genome_taxon_ids;
            for (const auto& info : accumulated_genome_info_) {
                genome_taxon_ids.push_back(info.taxon_id);
            }
            
            // Include uniqueness calculation
            if (config_.enable_uniqueness_scoring) {
                if (!integrate_uniqueness_with_feature_extractor(
                        feature_extractor_.get(),
                        d_minimizer_hits, 
                        total_hits_extracted, 
                        genome_taxon_ids)) {
                    std::cerr << "Uniqueness integration failed" << std::endl;
                    return false;
                }
            } else {
                // Original feature extraction without uniqueness
                if (!feature_extractor_->process_first_pass(
                        d_minimizer_hits, total_hits_extracted, genome_taxon_ids)) {
                    std::cerr << "Feature extraction failed" << std::endl;
                    return false;
                }
            }
            std::cout << "✓ Feature extraction completed" << std::endl;
            
            // Collect and print uniqueness statistics
            if (config_.enable_debug_mode && config_.enable_uniqueness_scoring) {
                std::vector<GPUMinimizerHit> h_hits(total_hits_extracted);
                cudaMemcpy(h_hits.data(), d_minimizer_hits, 
                           total_hits_extracted * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
                
                UniquenessStats uniqueness_stats = collect_uniqueness_statistics(h_hits);
                uniqueness_stats.print();
                
                // Update enhanced stats
                enhanced_stats_.minimizers_with_uniqueness_scores = total_hits_extracted;
                enhanced_stats_.unique_minimizers_count = uniqueness_stats.unique_minimizers;
                enhanced_stats_.rare_minimizers_count = uniqueness_stats.rare_minimizers;
                enhanced_stats_.reliable_minimizers_count = uniqueness_stats.reliable_minimizers;
                enhanced_stats_.average_uniqueness_score = uniqueness_stats.average_uniqueness;
            }
        }
        
        // Add co-occurrence scoring after uniqueness computation
        if (config_.enable_cooccurrence_scoring) {
            std::cout << "Computing co-occurrence scores..." << std::endl;
            
            // Get unique minimizers for co-occurrence analysis
            std::vector<uint64_t> unique_minimizers;
            std::unordered_set<uint64_t> seen_hashes;
            
            // Extract unique minimizer hashes from hits
            std::vector<GPUMinimizerHit> h_hits(total_hits_extracted);
            cudaMemcpy(h_hits.data(), d_minimizer_hits, 
                       total_hits_extracted * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
            
            for (const auto& hit : h_hits) {
                if (seen_hashes.insert(hit.minimizer_hash).second) {
                    unique_minimizers.push_back(hit.minimizer_hash);
                }
            }
            
            size_t num_unique = unique_minimizers.size();
            std::cout << "Found " << num_unique << " unique minimizers for co-occurrence analysis" << std::endl;
            
            // Sort unique minimizers for binary search in update kernel
            std::sort(unique_minimizers.begin(), unique_minimizers.end());
            
            // Allocate GPU memory for unique minimizers
            uint64_t* d_unique_minimizers;
            cudaMalloc(&d_unique_minimizers, num_unique * sizeof(uint64_t));
            cudaMemcpy(d_unique_minimizers, unique_minimizers.data(), 
                       num_unique * sizeof(uint64_t), cudaMemcpyHostToDevice);
            
            // Allocate memory for co-occurrence scores
            float* d_cooccurrence_scores;
            cudaMalloc(&d_cooccurrence_scores, num_unique * sizeof(float));
            
            // Use the optimized kernel if we have enough minimizers
            bool use_optimized = (num_unique > 1000 && total_hits_extracted > 10000);
            
            if (use_optimized) {
                int block_size = 128;
                int grid_size = (num_unique + block_size - 1) / block_size;
                size_t shared_mem_size = block_size * sizeof(GPUMinimizerHit);
                
                compute_cooccurrence_scores_optimized_kernel<<<grid_size, block_size, shared_mem_size>>>(
                    d_minimizer_hits,
                    d_unique_minimizers,
                    d_cooccurrence_scores,
                    total_hits_extracted,
                    num_unique,
                    config_.cooccurrence_window_size
                );
            } else {
                int block_size = 256;
                int grid_size = (num_unique + block_size - 1) / block_size;
                
                compute_cooccurrence_scores_kernel<<<grid_size, block_size>>>(
                    d_minimizer_hits,
                    d_unique_minimizers,
                    d_cooccurrence_scores,
                    total_hits_extracted,
                    num_unique,
                    config_.cooccurrence_window_size
                );
            }
            
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "Co-occurrence kernel launch failed: " << cudaGetErrorString(err) << std::endl;
                cudaFree(d_unique_minimizers);
                cudaFree(d_cooccurrence_scores);
                return false;
            }
            
            cudaDeviceSynchronize();
            
            // Update feature flags with scores
            int update_block_size = 256;
            int update_grid_size = (total_hits_extracted + update_block_size - 1) / update_block_size;
            
            update_cooccurrence_flags_kernel<<<update_grid_size, update_block_size>>>(
                d_minimizer_hits,
                d_unique_minimizers,
                d_cooccurrence_scores,
                total_hits_extracted,
                num_unique
            );
            
            cudaDeviceSynchronize();
            
            // Check for errors
            err = cudaGetLastError();
            if (err != cudaSuccess) {
                std::cerr << "Co-occurrence update kernel error: " << cudaGetErrorString(err) << std::endl;
                cudaFree(d_unique_minimizers);
                cudaFree(d_cooccurrence_scores);
                return false;
            }
            
            // Collect co-occurrence statistics if in debug mode
            if (config_.enable_debug_mode) {
                std::vector<float> cooccurrence_scores(num_unique);
                cudaMemcpy(cooccurrence_scores.data(), d_cooccurrence_scores, 
                           num_unique * sizeof(float), cudaMemcpyDeviceToHost);
                
                // Copy updated hits back for statistics
                cudaMemcpy(h_hits.data(), d_minimizer_hits, 
                           total_hits_extracted * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
                
                // Collect statistics
                size_t high_cooccurrence = 0;
                size_t low_cooccurrence = 0;
                double score_sum = 0.0;
                
                for (const auto& hit : h_hits) {
                    uint8_t category = MinimizerFlags::get_cooccurrence_score(hit.feature_flags);
                    if (category >= 5) high_cooccurrence++;
                    if (category <= 2) low_cooccurrence++;
                }
                
                for (const auto& score : cooccurrence_scores) {
                    score_sum += score;
                }
                
                // Update enhanced stats
                enhanced_stats_.minimizers_with_cooccurrence_scores = total_hits_extracted;
                enhanced_stats_.high_cooccurrence_minimizers = high_cooccurrence;
                enhanced_stats_.low_cooccurrence_minimizers = low_cooccurrence;
                enhanced_stats_.average_cooccurrence_score = score_sum / num_unique;
                
                std::cout << "Co-occurrence analysis results:" << std::endl;
                std::cout << "  High co-occurrence minimizers: " << high_cooccurrence 
                          << " (" << (100.0 * high_cooccurrence / total_hits_extracted) << "%)" << std::endl;
                std::cout << "  Low co-occurrence minimizers: " << low_cooccurrence 
                          << " (" << (100.0 * low_cooccurrence / total_hits_extracted) << "%)" << std::endl;
                std::cout << "  Average co-occurrence score: " << enhanced_stats_.average_cooccurrence_score << std::endl;
            }
            
            // Clean up
            cudaFree(d_unique_minimizers);
            cudaFree(d_cooccurrence_scores);
            
            std::cout << "✓ Co-occurrence scoring completed" << std::endl;
        }
        
        // Copy minimizer hits back to host
        std::vector<GPUMinimizerHit> minimizer_hits(total_hits_extracted);
        cudaMemcpy(minimizer_hits.data(), d_minimizer_hits, 
                   total_hits_extracted * sizeof(GPUMinimizerHit), cudaMemcpyDeviceToHost);
        
        // Collect statistics
        std::vector<uint32_t> genome_lengths;
        for (const auto& info : accumulated_genome_info_) {
            genome_lengths.push_back(info.sequence_length);
        }
        collect_minimizer_statistics(minimizer_hits, genome_lengths);
        
        // Convert to LCA candidates
        std::cout << "Computing LCA assignments for " << total_hits_extracted << " minimizers..." << std::endl;
        
        size_t new_candidates_start = all_lca_candidates_.size();
        for (const auto& hit : minimizer_hits) {
            LCACandidate candidate;
            candidate.minimizer_hash = hit.minimizer_hash;
            candidate.lca_taxon = hit.taxon_id;
            candidate.genome_count = 1;
            candidate.uniqueness_score = MinimizerFlags::ml_weight_to_float(hit.ml_weight);
            all_lca_candidates_.push_back(candidate);
        }
        
        // Merge and deduplicate new candidates
        if (new_candidates_start > 0) {
            std::sort(all_lca_candidates_.begin() + new_candidates_start, all_lca_candidates_.end(),
                     [](const LCACandidate& a, const LCACandidate& b) {
                         return a.minimizer_hash < b.minimizer_hash;
                     });
            
            std::inplace_merge(all_lca_candidates_.begin(),
                              all_lca_candidates_.begin() + new_candidates_start,
                              all_lca_candidates_.end(),
                              [](const LCACandidate& a, const LCACandidate& b) {
                                  return a.minimizer_hash < b.minimizer_hash;
                              });
        }
        
        std::cout << "✓ LCA assignment completed. Total candidates: " << all_lca_candidates_.size() << std::endl;
    }
    
    // Clear accumulation variables after processing
    accumulated_genome_info_.clear();
    accumulated_sequence_length_ = 0;
    accumulated_sequences_.clear();
    
    return true;
}

bool GPUKrakenDatabaseBuilder::enhance_with_phylogenetic_data() {
    if (!taxonomy_processor_->is_loaded()) {
        std::cerr << "Cannot enhance without loaded taxonomy" << std::endl;
        return false;
    }
    
    std::cout << "Enhancing candidates with phylogenetic data..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Group candidates by minimizer hash
    std::unordered_map<uint64_t, std::vector<LCACandidate*>> minimizer_groups;
    
    for (auto& candidate : all_lca_candidates_) {
        minimizer_groups[candidate.minimizer_hash].push_back(&candidate);
    }
    
    // Process each group
    phylogenetic_candidates_.clear();
    phylogenetic_candidates_.reserve(minimizer_groups.size());
    
    for (auto& [minimizer_hash, candidates] : minimizer_groups) {
        // Collect species from candidates
        std::vector<uint32_t> contributing_species;
        std::vector<uint16_t> genome_counts_per_species;
        
        std::unordered_map<uint32_t, uint16_t> species_counts;
        for (LCACandidate* candidate : candidates) {
            species_counts[candidate->lca_taxon] += candidate->genome_count;
        }
        
        for (const auto& [species, count] : species_counts) {
            contributing_species.push_back(species);
            genome_counts_per_species.push_back(count);
        }
        
        // Compute phylogenetic data
        uint32_t lca = taxonomy_processor_->compute_lca_of_species(contributing_species);
        uint8_t phylo_spread = taxonomy_processor_->calculate_phylogenetic_spread(contributing_species, lca);
        
        uint8_t max_distance = 0;
        for (uint32_t species : contributing_species) {
            uint8_t distance = taxonomy_processor_->calculate_distance_to_lca(species, lca);
            max_distance = std::max(max_distance, distance);
        }
        
        // Create enhanced candidate
        PhylogeneticLCACandidate phylo_candidate;
        phylo_candidate.minimizer_hash = minimizer_hash;
        phylo_candidate.lca_taxon = lca;
        phylo_candidate.contributing_species = contributing_species;
        phylo_candidate.genome_counts_per_species = genome_counts_per_species;
        phylo_candidate.phylogenetic_spread = phylo_spread;
        phylo_candidate.max_phylogenetic_distance = max_distance;
        
        phylo_candidate.genome_count = 0;
        for (uint16_t count : genome_counts_per_species) {
            phylo_candidate.genome_count += count;
        }
        phylo_candidate.uniqueness_score = 1.0f / phylo_candidate.genome_count;
        
        phylogenetic_candidates_.push_back(phylo_candidate);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    enhanced_stats_.phylogenetic_processing_time = 
        std::chrono::duration<double>(end_time - start_time).count();
    
    enhanced_stats_.minimizers_with_phylo_data = phylogenetic_candidates_.size();
    enhanced_stats_.phylogenetic_lca_computations = phylogenetic_candidates_.size();
    
    std::cout << "✓ Phylogenetic enhancement completed" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::save_database_outputs() {
    std::cout << "Saving database outputs..." << std::endl;
    
    bool success = true;
    
    // Always save standard format
    if (!standard_serializer_->save_standard_database(all_lca_candidates_, taxon_names_, 
                                                     taxon_parents_, build_stats_)) {
        std::cerr << "Failed to save standard database format" << std::endl;
        success = false;
    }
    
    // Save enhanced format if requested and available
    if (config_.save_enhanced_format && !phylogenetic_candidates_.empty()) {
        if (!enhanced_serializer_->save_enhanced_database(phylogenetic_candidates_, 
                                                         contributing_taxa_arrays_,
                                                         taxon_names_, enhanced_stats_)) {
            std::cerr << "Failed to save enhanced database format" << std::endl;
            success = false;
        }
    }
    
    return success;
}

// Validation and monitoring methods

bool GPUKrakenDatabaseBuilder::validate_configuration() const {
    if (output_directory_.empty()) {
        std::cerr << "Output directory not specified" << std::endl;
        return false;
    }
    
    if (config_.k_value < 15 || config_.k_value > 31) {
        std::cerr << "Invalid k value: " << config_.k_value << std::endl;
        return false;
    }
    
    if (config_.subsampling_rate <= 0.0 || config_.subsampling_rate > 1.0) {
        std::cerr << "Invalid subsampling rate: " << config_.subsampling_rate << std::endl;
        return false;
    }
    
    return true;
}

bool GPUKrakenDatabaseBuilder::coordinate_memory_allocation() {
    std::cout << "Coordinating memory allocation..." << std::endl;
    
    // Estimate memory requirements
    size_t max_sequences = config_.memory_config.sequence_batch_size;
    size_t max_sequence_length = max_sequences * 5000000; // 5MB per sequence estimate
    
    // Allocate through memory manager
    if (!memory_manager_->allocate_sequence_memory(max_sequences, max_sequence_length)) {
        return false;
    }
    
    if (!memory_manager_->allocate_minimizer_memory(config_.memory_config.minimizer_capacity)) {
        return false;
    }
    
    if (!memory_manager_->allocate_results_memory(config_.memory_config.minimizer_capacity)) {
        return false;
    }
    
    std::cout << "✓ Memory allocation coordinated successfully" << std::endl;
    return true;
}

bool GPUKrakenDatabaseBuilder::validate_inter_module_consistency() {
    if (genome_files_.size() != genome_taxon_ids_.size()) {
        std::cerr << "Inconsistent genome file and taxon ID counts" << std::endl;
        return false;
    }
    
    if (!modules_initialized_) {
        std::cerr << "Modules not properly initialized" << std::endl;
        return false;
    }
    
    return true;
}

void GPUKrakenDatabaseBuilder::update_build_progress(const std::string& stage, double progress_percent) {
    if (config_.enable_progress_reporting) {
        std::cout << "[" << stage << "] " << std::fixed << std::setprecision(1) 
                  << progress_percent << "% complete" << std::endl;
    }
}

void GPUKrakenDatabaseBuilder::record_timing_checkpoint(const std::string& stage_name) {
    if (config_.enable_detailed_timing) {
        auto now = std::chrono::high_resolution_clock::now();
        auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()).count();
        
        std::cout << "TIMING: " << stage_name << " at " << timestamp << "ms" << std::endl;
    }
}

void GPUKrakenDatabaseBuilder::print_build_progress() const {
    std::cout << "\n=== BUILD STATISTICS ===" << std::endl;
    std::cout << "Total sequences: " << build_stats_.total_sequences << std::endl;
    std::cout << "Total bases: " << build_stats_.total_bases << std::endl;
    std::cout << "Minimizers extracted: " << build_stats_.valid_minimizers_extracted << std::endl;
    std::cout << "LCA assignments: " << build_stats_.lca_assignments << std::endl;
    std::cout << "Processing time: " << std::fixed << std::setprecision(2) 
              << build_stats_.sequence_processing_time << "s" << std::endl;
    
    // Print enhanced statistics if available
    if (config_.enable_phylogenetic_analysis) {
        enhanced_stats_.print_enhanced_stats();
    }
    
    // NEW: Print uniqueness statistics
    if (config_.enable_uniqueness_scoring) {
        enhanced_stats_.print_uniqueness_stats();
    }
}

void GPUKrakenDatabaseBuilder::print_memory_usage() const {
    if (memory_manager_) {
        memory_manager_->print_memory_usage();
    }
}

void GPUKrakenDatabaseBuilder::aggregate_module_statistics() {
    // Aggregate statistics from all modules
    if (file_processor_) {
        const auto& file_stats = file_processor_->get_statistics();
        build_stats_.sequence_processing_time += file_stats.processing_time;
    }
    
    enhanced_stats_.unique_minimizers = all_lca_candidates_.size();
    build_stats_.unique_minimizers = all_lca_candidates_.size();
    
    if (config_.enable_phylogenetic_analysis) {
        enhanced_stats_.minimizers_with_phylo_data = phylogenetic_candidates_.size();
    }
}

void GPUKrakenDatabaseBuilder::cleanup_on_failure() {
    // Clean up resources on failure
    if (memory_manager_) {
        memory_manager_->free_all_allocations();
    }
    
    // Clear data structures
    genome_files_.clear();
    genome_taxon_ids_.clear();
    all_lca_candidates_.clear();
    phylogenetic_candidates_.clear();
}

bool GPUKrakenDatabaseBuilder::validate_database_output() {
    std::cout << "Validating database output..." << std::endl;
    
    // Check that output files exist
    std::string hash_table_file = output_directory_ + "/hash_table.k2d";
    std::string taxonomy_file = output_directory_ + "/taxonomy.tsv";
    std::string config_file = output_directory_ + "/config.txt";
    
    bool valid = true;
    
    if (!std::filesystem::exists(hash_table_file)) {
        std::cerr << "Missing hash table file: " << hash_table_file << std::endl;
        valid = false;
    }
    
    if (!std::filesystem::exists(taxonomy_file)) {
        std::cerr << "Missing taxonomy file: " << taxonomy_file << std::endl;
        valid = false;
    }
    
    if (!std::filesystem::exists(config_file)) {
        std::cerr << "Missing config file: " << config_file << std::endl;
        valid = false;
    }
    
    if (valid) {
        std::cout << "✓ Database output validation passed" << std::endl;
    } else {
        std::cerr << "✗ Database output validation failed" << std::endl;
    }
    
    return valid;
}

// ===========================
// Factory Implementation
// ===========================

std::unique_ptr<GPUKrakenDatabaseBuilder> DatabaseBuilderFactory::create_standard_builder(
    const std::string& output_dir) {
    
    DatabaseBuildConfig config;
    config.memory_config.minimizer_capacity = 5000000;
    config.memory_config.sequence_batch_size = 25;
    config.output_format = DatabaseFormat::STANDARD_KRAKEN2;
    
    return std::make_unique<GPUKrakenDatabaseBuilder>(output_dir, config);
}

std::unique_ptr<GPUKrakenDatabaseBuilder> DatabaseBuilderFactory::create_phylogenetic_builder(
    const std::string& output_dir) {
    
    DatabaseBuildConfig config;
    config.memory_config.minimizer_capacity = 10000000;
    config.memory_config.sequence_batch_size = 30;
    config.enable_phylogenetic_analysis = true;
    config.save_enhanced_format = true;
    config.output_format = DatabaseFormat::ENHANCED_PHYLO;
    
    return std::make_unique<GPUKrakenDatabaseBuilder>(output_dir, config);
}

// Static configuration creators
DatabaseBuildConfig GPUKrakenDatabaseBuilder::create_default_config() {
    return DatabaseBuildConfig();
}

DatabaseBuildConfig GPUKrakenDatabaseBuilder::create_high_memory_config() {
    DatabaseBuildConfig config;
    config.memory_config.minimizer_capacity = 20000000;
    config.memory_config.sequence_batch_size = 50;
    config.memory_config.max_memory_fraction = 90;
    return config;
}

DatabaseBuildConfig GPUKrakenDatabaseBuilder::create_enhanced_phylo_config() {
    DatabaseBuildConfig config;
    config.enable_phylogenetic_analysis = true;
    config.save_enhanced_format = true;
    config.memory_config.minimizer_capacity = 15000000;
    config.output_format = DatabaseFormat::ENHANCED_PHYLO;
    return config;
}

// Feature export method for ML training
bool GPUKrakenDatabaseBuilder::export_features_for_training(
    const FeatureExportConfig& export_config) {
    
    std::cout << "\n=== EXPORTING FEATURES FOR ML TRAINING ===" << std::endl;
    
    // Check if we have processed data
    if (all_lca_candidates_.empty()) {
        std::cerr << "No minimizer data available. Build database first." << std::endl;
        return false;
    }
    
    // Create feature exporter
    FeatureExporter exporter(export_config);
    
    // Convert LCA candidates to minimizer hits for export
    std::vector<GPUMinimizerHit> minimizer_hits;
    minimizer_hits.reserve(all_lca_candidates_.size());
    
    for (const auto& candidate : all_lca_candidates_) {
        GPUMinimizerHit hit;
        hit.minimizer_hash = candidate.minimizer_hash;
        hit.taxon_id = candidate.lca_taxon;
        hit.genome_id = 0; // Would need to track this
        hit.position = 0; // Would need to track this
        
        // Encode features if available from minimizer_stats_
        auto stats_it = minimizer_stats_.find(candidate.minimizer_hash);
        if (stats_it != minimizer_stats_.end()) {
            const auto& stats = stats_it->second;
            
            // Encode GC content category (0-7)
            uint8_t gc_category = static_cast<uint8_t>(
                std::min(7, static_cast<int>(stats.get_average_gc_content() / 12.5f))
            );
            
            // Encode complexity score (0-7)
            uint8_t complexity_score = static_cast<uint8_t>(
                std::min(7, static_cast<int>(stats.get_average_complexity() * 8))
            );
            
            // Build feature flags
            hit.feature_flags = gc_category | (complexity_score << 3);
            
            // Add contamination flag if high risk
            if (stats.get_contamination_risk_ratio() > 0.5f) {
                hit.feature_flags |= (1 << 7);
            }
            
            // Encode uniqueness as ML weight
            hit.ml_weight = static_cast<uint16_t>(candidate.uniqueness_score * 65535);
        }
        
        minimizer_hits.push_back(hit);
    }
    
    // Export features
    bool success = exporter.export_training_features(
        minimizer_hits, 
        *feature_extractor_,
        taxon_names_
    );
    
    if (success) {
        std::cout << "Feature export completed successfully!" << std::endl;
        std::cout << "Exported " << exporter.get_total_features_exported() 
                  << " unique minimizer features" << std::endl;
    }
    
    return success;
}

// ===========================
// Uniqueness Integration Implementation
// ===========================
// Implementation is included at the top of the file to avoid forward declaration issues
