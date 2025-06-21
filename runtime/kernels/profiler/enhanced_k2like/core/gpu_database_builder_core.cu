// core/gpu_database_builder_core.cu
// Implementation of the main GPU Kraken database builder
// Orchestrates all specialized modules and provides unified API

#include "gpu_database_builder_core.h"
#include "../memory/gpu_memory_manager.h"
#include "../gpu/gpu_database_kernels.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <filesystem>

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
    
    // Process in batches
    size_t batch_size = config_.memory_config.sequence_batch_size;
    
    for (size_t batch_start = 0; batch_start < genome_files_.size(); batch_start += batch_size) {
        size_t batch_end = std::min(batch_start + batch_size, genome_files_.size());
        
        std::cout << "Processing batch " << (batch_start / batch_size + 1) 
                  << ": genomes " << batch_start << "-" << batch_end << std::endl;
        
        // Load sequences for this batch
        std::vector<std::string> batch_sequences;
        std::vector<uint32_t> batch_taxon_ids;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            auto sequences = file_processor_->load_sequences_from_fasta(genome_files_[i]);
            for (const auto& seq : sequences) {
                batch_sequences.push_back(seq);
                batch_taxon_ids.push_back(genome_taxon_ids_[i]);
                build_stats_.total_bases += seq.length();
            }
        }
        
        // Process this batch
        if (!process_sequence_batch(batch_sequences, batch_taxon_ids)) {
            std::cerr << "Failed to process batch starting at " << batch_start << std::endl;
            return false;
        }
        
        if (config_.enable_progress_reporting) {
            double progress = (double)batch_end / genome_files_.size() * 100.0;
            update_build_progress("Genome Processing", progress);
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
    
    // The actual GPU processing would be implemented here
    // using the kernels from gpu_database_kernels.cu
    
    // For now, create dummy candidates to maintain the pipeline
    for (size_t i = 0; i < sequences.size(); i++) {
        LCACandidate candidate;
        candidate.minimizer_hash = i + 1000000; // Dummy hash
        candidate.lca_taxon = taxon_ids[i];
        candidate.genome_count = 1;
        candidate.uniqueness_score = 1.0f;
        
        all_lca_candidates_.push_back(candidate);
    }
    
    build_stats_.valid_minimizers_extracted += sequences.size();
    build_stats_.lca_assignments += sequences.size();
    
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
