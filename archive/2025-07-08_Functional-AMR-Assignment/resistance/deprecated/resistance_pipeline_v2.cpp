// runtime/kernels/resistance/resistance_pipeline_v2.cpp
#include "resistance_pipeline_v2.h"
#include "fq_mutation_detector.cuh"
#include "../../common/gpu/bloom_filter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cuda_runtime.h>

// External functions for report generation (from existing implementation)
extern "C" {
    void generate_html_report(const char* output_path, const char* sample_name, 
                             const ResistanceMutation* mutations, int num_mutations);
    void generate_json_report(const char* output_path, const char* sample_name,
                             const ResistanceMutation* mutations, int num_mutations);
    void generate_text_report(const char* output_path, const char* sample_name,
                             const ResistanceMutation* mutations, int num_mutations);
}

namespace BioGPU {

// AlleleFrequencyAnalyzer implementation
void AlleleFrequencyAnalyzer::addObservation(const std::string& gene, int position, char amino_acid) {
    std::lock_guard<std::mutex> lock(mutex);
    position_counts[gene][position][amino_acid]++;
    total_counts[gene][position]++;
}

void AlleleFrequencyAnalyzer::addMultipleObservations(const std::string& gene, int position, 
                                                      char amino_acid, int count) {
    std::lock_guard<std::mutex> lock(mutex);
    position_counts[gene][position][amino_acid] += count;
    total_counts[gene][position] += count;
}

ResistanceResults::AlleleFrequency AlleleFrequencyAnalyzer::calculateFrequency(
    const std::string& gene, int position, char wildtype_aa) const {
    
    std::lock_guard<std::mutex> lock(mutex);
    ResistanceResults::AlleleFrequency freq;
    freq.gene_name = gene;
    freq.position = position;
    freq.wildtype_aa = wildtype_aa;
    
    auto gene_it = total_counts.find(gene);
    if (gene_it != total_counts.end()) {
        auto pos_it = gene_it->second.find(position);
        if (pos_it != gene_it->second.end()) {
            freq.total_coverage = pos_it->second.load();
            
            auto aa_it = position_counts.find(gene);
            if (aa_it != position_counts.end()) {
                auto pos_aa_it = aa_it->second.find(position);
                if (pos_aa_it != aa_it->second.end()) {
                    for (const auto& [aa, count] : pos_aa_it->second) {
                        float frequency = count.load() / freq.total_coverage;
                        freq.amino_acid_frequencies[aa] = frequency;
                        if (aa == wildtype_aa) {
                            freq.wildtype_frequency = frequency;
                        }
                    }
                }
            }
        }
    }
    
    return freq;
}

std::vector<ResistanceResults::AlleleFrequency> AlleleFrequencyAnalyzer::getAllFrequencies(
    const std::string& gene) const {
    
    std::lock_guard<std::mutex> lock(mutex);
    std::vector<ResistanceResults::AlleleFrequency> frequencies;
    
    auto gene_it = position_counts.find(gene);
    if (gene_it != position_counts.end()) {
        for (const auto& [position, aa_counts] : gene_it->second) {
            // Get wildtype AA for this position (would come from mutation database)
            char wildtype_aa = 'X';  // Default, should be looked up
            frequencies.push_back(calculateFrequency(gene, position, wildtype_aa));
        }
    }
    
    return frequencies;
}

void AlleleFrequencyAnalyzer::clear() {
    std::lock_guard<std::mutex> lock(mutex);
    position_counts.clear();
    total_counts.clear();
}

// ResistanceResults implementation
void ResistanceResults::writeReport(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open output file " << filename << std::endl;
        return;
    }
    
    out << "# Fluoroquinolone Resistance Detection Report\n";
    out << "# Total sequences processed: " << sequences_processed << "\n";
    out << "# Total bases processed: " << bases_processed << "\n";
    out << "# Processing time: " << processing_time_seconds << " seconds\n";
    out << "# Total resistance mutations detected: " << all_hits.size() << "\n\n";
    
    // Header
    out << "sample\tgene\tspecies\tposition\twildtype\tmutant\tmutation\t"
        << "confidence\tread_support\tis_qrdr\tdrug_class\n";
    
    // Write all hits
    for (const auto& hit : all_hits) {
        out << hit.sample_name << "\t"
            << hit.gene_name << "\t"
            << hit.species << "\t"
            << hit.position << "\t"
            << hit.wildtype_aa << "\t"
            << hit.mutant_aa << "\t"
            << hit.mutation_name << "\t"
            << std::fixed << std::setprecision(3) << hit.confidence << "\t"
            << hit.read_support << "\t"
            << (hit.is_qrdr ? "Yes" : "No") << "\t"
            << "Fluoroquinolone\n";
    }
}

void ResistanceResults::writeSummary(std::ostream& out) const {
    out << "Fluoroquinolone Resistance Detection Summary:\n";
    out << "  Total sequences: " << sequences_processed << "\n";
    out << "  Samples analyzed: " << sample_profiles.size() << "\n";
    
    int resistant_samples = 0;
    std::map<std::string, int> mutation_counts;
    
    for (const auto& [sample, profile] : sample_profiles) {
        if (profile.has_fluoroquinolone_resistance) {
            resistant_samples++;
        }
        for (const auto& mutation : profile.mutations) {
            mutation_counts[mutation.mutation_name]++;
        }
    }
    
    out << "  Resistant samples: " << resistant_samples << "/" << sample_profiles.size() << "\n";
    out << "  Most common mutations:\n";
    
    // Sort mutations by frequency
    std::vector<std::pair<std::string, int>> sorted_mutations(mutation_counts.begin(), mutation_counts.end());
    std::sort(sorted_mutations.begin(), sorted_mutations.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    for (size_t i = 0; i < std::min(size_t(5), sorted_mutations.size()); i++) {
        out << "    " << sorted_mutations[i].first << ": " 
            << sorted_mutations[i].second << " samples\n";
    }
}

void ResistanceResults::writeClinicalReport(const std::string& filename, const std::string& format) const {
    if (format == "html") {
        // Convert mutations to C-style array for external function
        std::vector<ResistanceMutation> mutations;
        for (const auto& hit : all_hits) {
            ResistanceMutation mut;
            mut.id = mutations.size();
            strncpy(mut.gene_name, hit.gene_name.c_str(), sizeof(mut.gene_name) - 1);
            strncpy(mut.species, hit.species.c_str(), sizeof(mut.species) - 1);
            mut.position = hit.position;
            mut.wildtype_aa = hit.wildtype_aa;
            mut.mutant_aa = hit.mutant_aa;
            strncpy(mut.mutation_name, hit.mutation_name.c_str(), sizeof(mut.mutation_name) - 1);
            mutations.push_back(mut);
        }
        
        generate_html_report(filename.c_str(), "Combined", mutations.data(), mutations.size());
    } else if (format == "json") {
        writeJSONReport(filename);
    } else {
        writeReport(filename);  // Default TSV format
    }
}

void ResistanceResults::writeAlleleFrequencyReport(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open allele frequency file " << filename << std::endl;
        return;
    }
    
    out << "sample,gene,position,amino_acid,frequency,is_wildtype,total_coverage\n";
    
    for (const auto& [sample, profile] : sample_profiles) {
        for (const auto& [gene, frequencies] : profile.allele_frequencies) {
            for (const auto& freq : frequencies) {
                for (const auto& [aa, frequency] : freq.amino_acid_frequencies) {
                    out << sample << ","
                        << gene << ","
                        << freq.position << ","
                        << aa << ","
                        << std::fixed << std::setprecision(4) << frequency << ","
                        << (aa == freq.wildtype_aa ? "Yes" : "No") << ","
                        << freq.total_coverage << "\n";
                }
            }
        }
    }
}

void ResistanceResults::writeJSONReport(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open JSON file " << filename << std::endl;
        return;
    }
    
    out << "{\n";
    out << "  \"metadata\": {\n";
    out << "    \"sequences_processed\": " << sequences_processed << ",\n";
    out << "    \"bases_processed\": " << bases_processed << ",\n";
    out << "    \"processing_time_seconds\": " << processing_time_seconds << "\n";
    out << "  },\n";
    out << "  \"samples\": [\n";
    
    bool first_sample = true;
    for (const auto& [sample, profile] : sample_profiles) {
        if (!first_sample) out << ",\n";
        first_sample = false;
        
        out << "    {\n";
        out << "      \"sample_name\": \"" << sample << "\",\n";
        out << "      \"has_resistance\": " << (profile.has_fluoroquinolone_resistance ? "true" : "false") << ",\n";
        out << "      \"mutations\": [\n";
        
        bool first_mut = true;
        for (const auto& mutation : profile.mutations) {
            if (!first_mut) out << ",\n";
            first_mut = false;
            
            out << "        {\n";
            out << "          \"gene\": \"" << mutation.gene_name << "\",\n";
            out << "          \"position\": " << mutation.position << ",\n";
            out << "          \"wildtype\": \"" << mutation.wildtype_aa << "\",\n";
            out << "          \"mutant\": \"" << mutation.mutant_aa << "\",\n";
            out << "          \"name\": \"" << mutation.mutation_name << "\",\n";
            out << "          \"confidence\": " << mutation.confidence << ",\n";
            out << "          \"read_support\": " << mutation.read_support << ",\n";
            out << "          \"is_qrdr\": " << (mutation.is_qrdr ? "true" : "false") << "\n";
            out << "        }";
        }
        
        out << "\n      ]\n";
        out << "    }";
    }
    
    out << "\n  ]\n";
    out << "}\n";
}

// ResistancePipeline implementation
ResistancePipeline::ResistancePipeline(const ResistanceConfig& config)
    : PipelineBase("Fluoroquinolone Resistance", config), resistance_config(config) {
    accumulated_results = std::make_unique<ResistanceResults>();
    fq_mapper = std::make_unique<GlobalFQResistanceMapper>();
}

bool ResistancePipeline::initializePipeline() {
    log("Initializing fluoroquinolone resistance detection pipeline");
    
    // Load FQ mutation database
    if (!loadFQMutationDatabase()) {
        logError("Failed to load FQ mutation database");
        return false;
    }
    
    // Initialize QRDR regions
    initializeQRDRRegions();
    
    // Build k-mer table for nucleotide matching
    if (!buildKmerTable()) {
        logError("Failed to build k-mer table");
        return false;
    }
    
    // Prepare protein references
    if (!prepareProteinReferences()) {
        logError("Failed to prepare protein references");
        return false;
    }
    
    // Initialize bloom filter
    bloom_filter = std::make_unique<BloomFilter>(
        resistance_config.bloom_filter_size,
        resistance_config.bloom_filter_hashes,
        resistance_config.kmer_size
    );
    
    if (!bloom_filter->initialize()) {
        logError("Failed to initialize bloom filter");
        return false;
    }
    
    // Build bloom filter from mutations
    if (!buildBloomFilterFromMutations()) {
        logError("Failed to build bloom filter from mutations");
        return false;
    }
    
    // Estimate and allocate GPU memory
    size_t required_memory = estimateMemoryRequirements(BATCH_SIZE);
    if (!allocateWorkspace(required_memory)) {
        logError("Failed to allocate GPU workspace");
        return false;
    }
    
    // Set up GPU memory pointers within workspace
    size_t offset = 0;
    
    // Stage 1: Bloom filter matches
    d_bloom_matches = (BloomFilterMatch*)((char*)d_workspace + offset);
    offset += BATCH_SIZE * sizeof(BloomFilterMatch);
    
    // Stage 2: Nucleotide matches
    d_nucleotide_matches = (NucleotideMatch*)((char*)d_workspace + offset);
    offset += BATCH_SIZE * 100 * sizeof(NucleotideMatch);  // Up to 100 matches per read
    
    d_nucleotide_match_counts = (uint32_t*)((char*)d_workspace + offset);
    offset += BATCH_SIZE * sizeof(uint32_t);
    
    // Stage 3: Protein matches
    d_protein_matches = (ProteinMatch*)((char*)d_workspace + offset);
    offset += BATCH_SIZE * 50 * sizeof(ProteinMatch);  // Up to 50 protein matches per read
    
    d_protein_match_counts = (uint32_t*)((char*)d_workspace + offset);
    offset += BATCH_SIZE * sizeof(uint32_t);
    
    log("Fluoroquinolone resistance pipeline initialized successfully");
    log("Loaded " + std::to_string(mutation_database.size()) + " resistance mutations");
    
    return true;
}

void ResistancePipeline::cleanupPipeline() {
    // Free k-mer table
    if (d_kmer_table) {
        cudaFree(d_kmer_table);
        d_kmer_table = nullptr;
    }
    
    // Free mutation sequences
    if (d_mutation_sequences) {
        cudaFree(d_mutation_sequences);
        d_mutation_sequences = nullptr;
    }
    
    // Clear results
    allele_analyzer.clear();
    accumulated_results.reset();
    
    // Base class cleanup handles workspace and streams
}

bool ResistancePipeline::processBatchImpl(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) {
    if (!gpu_buffer || gpu_buffer->getCurrentSequences() == 0) {
        return true;  // Nothing to process
    }
    
    log("Processing batch with " + std::to_string(gpu_buffer->getCurrentSequences()) + 
        " sequences for sample: " + current_sample_name);
    
    // Clear previous batch results
    clearBatchResults();
    
    // Three-stage pipeline execution
    stage1_bloomFilterScreening(gpu_buffer, stream);
    stage2_nucleotideKmerMatching(gpu_buffer, stream);
    stage3_translatedProteinSearch(gpu_buffer, stream);
    
    // Process results
    processProteinMatches(stream);
    
    // Synchronize to ensure all operations complete
    cudaStreamSynchronize(stream);
    
    return true;
}

void ResistancePipeline::stage1_bloomFilterScreening(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) {
    log("Stage 1: Bloom filter screening");
    
    // Use shared bloom filter component
    std::vector<bool> passed_filter = bloom_filter->screenReads(
        gpu_buffer->getSequences(),
        gpu_buffer->getOffsets(),
        gpu_buffer->getLengths(),
        gpu_buffer->getCurrentSequences(),
        stream
    );
    
    // Convert to GPU format
    std::vector<BloomFilterMatch> bloom_matches(gpu_buffer->getCurrentSequences());
    size_t passed_count = 0;
    
    for (size_t i = 0; i < gpu_buffer->getCurrentSequences(); i++) {
        bloom_matches[i].read_id = i;
        bloom_matches[i].passed_r1 = passed_filter[i];
        bloom_matches[i].passed_r2 = false;  // Handle R2 separately if paired-end
        if (passed_filter[i]) passed_count++;
    }
    
    // Copy to GPU
    cudaMemcpyAsync(d_bloom_matches, bloom_matches.data(), 
                    bloom_matches.size() * sizeof(BloomFilterMatch),
                    cudaMemcpyHostToDevice, stream);
    
    stage1_passed = passed_count;
    log("Stage 1 complete: " + std::to_string(passed_count) + "/" + 
        std::to_string(gpu_buffer->getCurrentSequences()) + " reads passed");
}

void ResistancePipeline::stage2_nucleotideKmerMatching(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) {
    log("Stage 2: Nucleotide k-mer matching");
    
    // Launch nucleotide matching kernel
    launchNucleotideMatchingKernel(
        gpu_buffer->getSequences(),
        gpu_buffer->getOffsets(),
        gpu_buffer->getLengths(),
        d_bloom_matches,
        gpu_buffer->getCurrentSequences(),
        d_kmer_table,
        d_nucleotide_matches,
        d_nucleotide_match_counts,
        stream
    );
    
    // Count matches
    std::vector<uint32_t> match_counts(gpu_buffer->getCurrentSequences());
    cudaMemcpyAsync(match_counts.data(), d_nucleotide_match_counts,
                    match_counts.size() * sizeof(uint32_t),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    stage2_passed = 0;
    for (auto count : match_counts) {
        if (count > 0) stage2_passed++;
    }
    
    log("Stage 2 complete: " + std::to_string(stage2_passed) + " reads with nucleotide matches");
}

void ResistancePipeline::stage3_translatedProteinSearch(GPUSequenceBuffer* gpu_buffer, cudaStream_t stream) {
    log("Stage 3: Translated protein search");
    
    // Launch protein translation and matching kernel
    launchProteinTranslationKernel(
        gpu_buffer->getSequences(),
        gpu_buffer->getOffsets(),
        gpu_buffer->getLengths(),
        d_nucleotide_matches,
        d_nucleotide_match_counts,
        gpu_buffer->getCurrentSequences(),
        d_mutation_sequences,
        d_protein_matches,
        d_protein_match_counts,
        stream
    );
    
    // Count matches
    std::vector<uint32_t> protein_counts(gpu_buffer->getCurrentSequences());
    cudaMemcpyAsync(protein_counts.data(), d_protein_match_counts,
                    protein_counts.size() * sizeof(uint32_t),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    stage3_passed = 0;
    for (auto count : protein_counts) {
        if (count > 0) stage3_passed++;
    }
    
    log("Stage 3 complete: " + std::to_string(stage3_passed) + " reads with protein matches");
}

bool ResistancePipeline::loadFQMutationDatabase() {
    // Load hardcoded mutations from fq_mutations_hardcoded.h
    const FQMutation* hardcoded_mutations = get_fq_mutations();
    int num_mutations = get_num_fq_mutations();
    
    mutation_database.clear();
    mutation_database.reserve(num_mutations);
    
    for (int i = 0; i < num_mutations; i++) {
        ResistanceMutation mut;
        mut.id = i;
        mut.gene_name = hardcoded_mutations[i].gene;
        mut.species = hardcoded_mutations[i].species;
        mut.position = hardcoded_mutations[i].position;
        mut.wildtype_aa = hardcoded_mutations[i].wildtype;
        mut.mutant_aa = hardcoded_mutations[i].mutant;
        mut.mutation_name = std::string(mut.gene_name) + "_" + 
                           mut.wildtype_aa + std::to_string(mut.position) + mut.mutant_aa;
        mut.drug_class = "Fluoroquinolone";
        mut.resistance_level = 1.0;  // High resistance
        mut.is_primary_mutation = true;
        mut.nucleotide_sequence = hardcoded_mutations[i].nucleotide_context;
        mut.protein_context = hardcoded_mutations[i].protein_context;
        
        mutation_database.push_back(mut);
    }
    
    return true;
}

bool ResistancePipeline::buildBloomFilterFromMutations() {
    // Extract k-mers from mutation nucleotide sequences
    std::vector<uint64_t> kmers;
    
    for (const auto& mutation : mutation_database) {
        if (mutation.nucleotide_sequence.length() >= resistance_config.kmer_size) {
            // Extract k-mers from the nucleotide context
            for (size_t i = 0; i <= mutation.nucleotide_sequence.length() - resistance_config.kmer_size; i++) {
                std::string kmer = mutation.nucleotide_sequence.substr(i, resistance_config.kmer_size);
                // Convert to canonical form and add to bloom filter
                uint64_t kmer_hash = 0;  // Would compute hash here
                kmers.push_back(kmer_hash);
            }
        }
    }
    
    // Add k-mers to bloom filter
    return bloom_filter->addKmers(kmers);
}

bool ResistancePipeline::buildKmerTable() {
    // Build k-mer lookup table for Stage 2
    size_t table_size = 1ULL << (2 * KMER_LENGTH);  // 4^k entries
    
    if (cudaMalloc(&d_kmer_table, table_size * sizeof(uint64_t)) != cudaSuccess) {
        logError("Failed to allocate k-mer table");
        return false;
    }
    
    // Initialize table with mutation k-mers
    // This would populate the table with k-mers from mutation sequences
    cudaMemset(d_kmer_table, 0, table_size * sizeof(uint64_t));
    
    return true;
}

bool ResistancePipeline::prepareProteinReferences() {
    // Prepare protein sequences for Stage 3
    size_t total_protein_length = 0;
    for (const auto& mutation : mutation_database) {
        total_protein_length += mutation.protein_context.length();
    }
    
    if (cudaMalloc(&d_mutation_sequences, total_protein_length + mutation_database.size()) != cudaSuccess) {
        logError("Failed to allocate mutation sequences");
        return false;
    }
    
    // Copy protein sequences to GPU
    std::string concatenated_proteins;
    for (const auto& mutation : mutation_database) {
        concatenated_proteins += mutation.protein_context + '\0';
    }
    
    cudaMemcpy(d_mutation_sequences, concatenated_proteins.c_str(), 
               concatenated_proteins.length(), cudaMemcpyHostToDevice);
    
    return true;
}

void ResistancePipeline::initializeQRDRRegions() {
    // Initialize QRDR regions for different species and genes
    qrdr_regions.clear();
    
    // E. coli GyrA QRDR
    QRDRRegion ecoli_gyra;
    ecoli_gyra.gene = "gyrA";
    ecoli_gyra.species = "Escherichia coli";
    ecoli_gyra.start_position = 67;
    ecoli_gyra.end_position = 106;
    ecoli_gyra.critical_positions = {83, 87};
    qrdr_regions["gyrA_ecoli"] = ecoli_gyra;
    
    // E. coli ParC QRDR
    QRDRRegion ecoli_parc;
    ecoli_parc.gene = "parC";
    ecoli_parc.species = "Escherichia coli";
    ecoli_parc.start_position = 64;
    ecoli_parc.end_position = 102;
    ecoli_parc.critical_positions = {80, 84};
    qrdr_regions["parC_ecoli"] = ecoli_parc;
    
    // Add more species/genes as needed
}

bool ResistancePipeline::isInQRDR(const std::string& gene, int position) const {
    for (const auto& [key, region] : qrdr_regions) {
        if (region.gene == gene && 
            position >= region.start_position && 
            position <= region.end_position) {
            return true;
        }
    }
    return false;
}

void ResistancePipeline::processProteinMatches(cudaStream_t stream) {
    // Get total number of protein matches
    std::vector<uint32_t> match_counts(BATCH_SIZE);
    cudaMemcpyAsync(match_counts.data(), d_protein_match_counts,
                    BATCH_SIZE * sizeof(uint32_t),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    size_t total_matches = 0;
    for (auto count : match_counts) {
        total_matches += count;
    }
    
    if (total_matches == 0) return;
    
    // Copy protein matches from GPU
    std::vector<ProteinMatch> protein_matches(total_matches);
    cudaMemcpyAsync(protein_matches.data(), d_protein_matches,
                    total_matches * sizeof(ProteinMatch),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    
    // Process matches
    aggregateMutationResults(protein_matches);
    updateAlleleFrequencies(protein_matches);
}

void ResistancePipeline::aggregateMutationResults(const std::vector<ProteinMatch>& matches) {
    // Group matches by mutation
    std::map<std::string, std::vector<const ProteinMatch*>> mutation_groups;
    
    for (const auto& match : matches) {
        auto mutation = getMutationInfo(match.gene_id);
        std::string key = mutation.mutation_name;
        mutation_groups[key].push_back(&match);
    }
    
    // Create resistance hits
    for (const auto& [mutation_name, match_group] : mutation_groups) {
        if (match_group.empty()) continue;
        
        auto mutation = getMutationInfo(match_group[0]->gene_id);
        
        ResistanceResults::ResistanceHit hit;
        hit.read_id = match_group[0]->read_id;
        hit.sample_name = current_sample_name;
        hit.gene_name = mutation.gene_name;
        hit.species = mutation.species;
        hit.position = mutation.position;
        hit.wildtype_aa = mutation.wildtype_aa;
        hit.mutant_aa = mutation.mutant_aa;
        hit.confidence = 0.0f;
        hit.is_qrdr = isInQRDR(mutation.gene_name, mutation.position);
        hit.read_support = match_group.size();
        hit.mutation_name = mutation_name;
        
        // Calculate confidence based on match quality
        for (const auto* match : match_group) {
            hit.confidence = std::max(hit.confidence, match->score);
        }
        
        accumulated_results->all_hits.push_back(hit);
        
        // Update sample profile
        auto& profile = accumulated_results->sample_profiles[current_sample_name];
        profile.sample_name = current_sample_name;
        profile.mutations.push_back(hit);
        
        // Mark as resistant if high-confidence QRDR mutation
        if (hit.is_qrdr && hit.confidence > 0.9f) {
            profile.has_fluoroquinolone_resistance = true;
            profile.resistance_mechanisms.push_back(mutation_name);
        }
    }
}

void ResistancePipeline::updateAlleleFrequencies(const std::vector<ProteinMatch>& matches) {
    for (const auto& match : matches) {
        auto mutation = getMutationInfo(match.gene_id);
        allele_analyzer.addObservation(mutation.gene_name, match.protein_position, match.amino_acid);
    }
}

void ResistancePipeline::calculateAlleleFrequencies() {
    // Calculate frequencies for all observed positions
    for (auto& [sample_name, profile] : accumulated_results->sample_profiles) {
        for (const auto& mutation : mutation_database) {
            auto frequencies = allele_analyzer.getAllFrequencies(mutation.gene_name);
            if (!frequencies.empty()) {
                profile.allele_frequencies[mutation.gene_name] = frequencies;
            }
        }
    }
}

void ResistancePipeline::generateClinicalReport(const std::string& output_path, const std::string& format) {
    accumulated_results->writeClinicalReport(output_path, format);
}

// Kernel launcher implementations (stubs - actual implementations would call CUDA kernels)
void ResistancePipeline::launchBloomFilterKernel(
    const char* d_sequences,
    const int* d_offsets,
    const int* d_lengths,
    int num_sequences,
    BloomFilterMatch* d_matches,
    cudaStream_t stream) {
    
    // This would launch the actual CUDA kernel from fq_mutation_detector.cu
    // For now, using the bloom filter component handles this
}

void ResistancePipeline::launchNucleotideMatchingKernel(
    const char* d_sequences,
    const int* d_offsets,
    const int* d_lengths,
    const BloomFilterMatch* d_bloom_matches,
    int num_sequences,
    const uint64_t* d_kmer_table,
    NucleotideMatch* d_matches,
    uint32_t* d_match_counts,
    cudaStream_t stream) {
    
    // Launch the nucleotide k-mer matching kernel
    // This adapts kmer_resistance_mapper.cu to work with GPUSequenceBuffer format
    extern void launch_nucleotide_kmer_matching(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        const BloomFilterMatch* d_bloom_matches,
        int num_sequences,
        const uint64_t* d_kmer_table,
        NucleotideMatch* d_matches,
        uint32_t* d_match_counts,
        int kmer_length,
        cudaStream_t stream);
    
    launch_nucleotide_kmer_matching(
        d_sequences, d_offsets, d_lengths, d_bloom_matches,
        num_sequences, d_kmer_table, d_matches, d_match_counts,
        KMER_LENGTH, stream);
}

void ResistancePipeline::launchProteinTranslationKernel(
    const char* d_sequences,
    const int* d_offsets,
    const int* d_lengths,
    const NucleotideMatch* d_nucleotide_matches,
    const uint32_t* d_match_counts,
    int num_sequences,
    const char* d_mutation_sequences,
    ProteinMatch* d_protein_matches,
    uint32_t* d_protein_match_counts,
    cudaStream_t stream) {
    
    // Launch the translated protein search kernel
    // This adapts translated_resistance_search.cu for the new format
    extern void launch_translated_protein_search(
        const char* d_sequences,
        const int* d_offsets,
        const int* d_lengths,
        const NucleotideMatch* d_nucleotide_matches,
        const uint32_t* d_match_counts,
        int num_sequences,
        const char* d_mutation_sequences,
        ProteinMatch* d_protein_matches,
        uint32_t* d_protein_match_counts,
        cudaStream_t stream);
    
    launch_translated_protein_search(
        d_sequences, d_offsets, d_lengths, d_nucleotide_matches, d_match_counts,
        num_sequences, d_mutation_sequences, d_protein_matches, d_protein_match_counts,
        stream);
}

size_t ResistancePipeline::estimateMemoryRequirements(size_t num_sequences) {
    size_t memory_needed = 0;
    
    // Stage 1: Bloom filter matches
    memory_needed += num_sequences * sizeof(BloomFilterMatch);
    
    // Stage 2: Nucleotide matches (up to 100 per read)
    memory_needed += num_sequences * 100 * sizeof(NucleotideMatch);
    memory_needed += num_sequences * sizeof(uint32_t);  // match counts
    
    // Stage 3: Protein matches (up to 50 per read)
    memory_needed += num_sequences * 50 * sizeof(ProteinMatch);
    memory_needed += num_sequences * sizeof(uint32_t);  // match counts
    
    return memory_needed;
}

void ResistancePipeline::clearBatchResults() {
    cudaMemset(d_nucleotide_match_counts, 0, BATCH_SIZE * sizeof(uint32_t));
    cudaMemset(d_protein_match_counts, 0, BATCH_SIZE * sizeof(uint32_t));
}

std::string ResistancePipeline::normalizeSpeciesName(const std::string& species) const {
    return fq_mapper->normalizeSpeciesName(species);
}

std::string ResistancePipeline::normalizeGeneName(const std::string& gene, const std::string& species) const {
    return fq_mapper->normalizeGeneName(gene, species);
}

ResistanceMutation ResistancePipeline::getMutationInfo(uint32_t gene_id) const {
    if (gene_id < mutation_database.size()) {
        return mutation_database[gene_id];
    }
    return ResistanceMutation();
}

} // namespace BioGPU