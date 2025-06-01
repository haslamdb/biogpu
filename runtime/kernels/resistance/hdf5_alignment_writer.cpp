// hdf5_alignment_writer.cpp
// HDF5 writer implementation with translated search support

#include "hdf5_alignment_writer.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>

// Include the ProteinMatch structure definition
// This should match the structure in translated_search.cu
struct ProteinMatch {
    uint32_t read_id;
    int8_t frame;
    uint32_t protein_id;
    uint32_t gene_id;
    uint32_t species_id;
    uint16_t query_start;
    uint16_t ref_start;
    uint16_t match_length;
    float alignment_score;
    float identity;
    uint8_t num_mutations;
    uint8_t mutation_positions[10];
    char ref_aas[10];
    char query_aas[10];
    float blosum_scores[10];
};

HDF5AlignmentWriter::HDF5AlignmentWriter(const std::string& output_path) 
    : filename(output_path), total_alignments_written(0), total_reads_processed(0),
      total_translated_alignments_written(0) {
    
    try {
        // Create new HDF5 file (overwrite if exists)
        file = new H5::H5File(output_path, H5F_ACC_TRUNC);
        
        // Set up compression
        H5::DSetCreatPropList plist;
        hsize_t chunk_dims[1] = {10000};  // Chunk size for compression
        plist.setChunk(1, chunk_dims);
        plist.setDeflate(6);  // Compression level 6
        
    } catch (H5::Exception& e) {
        std::cerr << "Error creating HDF5 file: " << e.getCDetailMsg() << std::endl;
        throw;
    }
}

HDF5AlignmentWriter::~HDF5AlignmentWriter() {
    if (file) {
        flush();  // Ensure all data is written
        flushTranslatedResults();  // Flush translated results too
        delete file;
    }
}

void HDF5AlignmentWriter::initialize(const std::string& index_path, 
                                    const std::string& r1_path, 
                                    const std::string& r2_path) {
    // Create group structure
    file->createGroup("/metadata");
    file->createGroup("/alignments");
    file->createGroup("/kmer_screening");
    file->createGroup("/ml_features");
    file->createGroup("/translated_search");  // NEW
    
    // Write metadata
    H5::Group metadata = file->openGroup("/metadata");
    
    // Write input file paths
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace scalar_space(H5S_SCALAR);
    
    H5::Attribute index_attr = metadata.createAttribute("index_path", str_type, scalar_space);
    index_attr.write(str_type, index_path);
    
    H5::Attribute r1_attr = metadata.createAttribute("r1_path", str_type, scalar_space);
    r1_attr.write(str_type, r1_path);
    
    H5::Attribute r2_attr = metadata.createAttribute("r2_path", str_type, scalar_space);
    r2_attr.write(str_type, r2_path);
    
    // Write timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    
    H5::Attribute time_attr = metadata.createAttribute("analysis_timestamp", str_type, scalar_space);
    time_attr.write(str_type, ss.str());
    
    // Write pipeline version (updated to indicate translated search support)
    H5::Attribute version_attr = metadata.createAttribute("pipeline_version", str_type, scalar_space);
    version_attr.write(str_type, std::string("0.4.0-translated"));
    
    createDatasets();
    createTranslatedDatasets();  // NEW: Create translated search datasets
}

void HDF5AlignmentWriter::createDatasets() {
    // Create extensible datasets for alignments
    hsize_t initial_dims[1] = {0};
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    H5::DataSpace dataspace(1, initial_dims, max_dims);
    
    // Set up chunking and compression
    H5::DSetCreatPropList plist;
    hsize_t chunk_dims[1] = {10000};
    plist.setChunk(1, chunk_dims);
    plist.setDeflate(6);
    
    // Create alignment datasets
    H5::Group align_group = file->openGroup("/alignments");
    
    align_group.createDataSet("read_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    align_group.createDataSet("gene_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    align_group.createDataSet("species_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    align_group.createDataSet("seq_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    align_group.createDataSet("alignment_score", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    align_group.createDataSet("identity", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    align_group.createDataSet("start_pos", H5::PredType::NATIVE_UINT16, dataspace, plist);
    align_group.createDataSet("matches", H5::PredType::NATIVE_UINT16, dataspace, plist);
    align_group.createDataSet("num_mutations", H5::PredType::NATIVE_UINT8, dataspace, plist);
    
    // Create datasets for mutation positions (variable length)
    align_group.createDataSet("mutation_positions", H5::PredType::NATIVE_UINT32, dataspace, plist);
    align_group.createDataSet("mutation_offsets", H5::PredType::NATIVE_UINT32, dataspace, plist);
    
    // Create k-mer screening datasets
    H5::Group kmer_group = file->openGroup("/kmer_screening");
    kmer_group.createDataSet("read_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    kmer_group.createDataSet("num_candidates", H5::PredType::NATIVE_UINT32, dataspace, plist);
    kmer_group.createDataSet("total_kmers_matched", H5::PredType::NATIVE_UINT32, dataspace, plist);
}

void HDF5AlignmentWriter::createTranslatedDatasets() {
    // Create translated search group
    H5::Group trans_group = file->openGroup("/translated_search");
    
    // Create extensible datasets
    hsize_t initial_dims[1] = {0};
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    H5::DataSpace dataspace(1, initial_dims, max_dims);
    
    // Set up chunking and compression
    H5::DSetCreatPropList plist;
    hsize_t chunk_dims[1] = {10000};
    plist.setChunk(1, chunk_dims);
    plist.setDeflate(6);
    
    // Basic match information
    trans_group.createDataSet("read_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    trans_group.createDataSet("frame", H5::PredType::NATIVE_INT8, dataspace, plist);
    trans_group.createDataSet("protein_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    trans_group.createDataSet("gene_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    trans_group.createDataSet("species_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    
    // Alignment positions
    trans_group.createDataSet("query_start", H5::PredType::NATIVE_UINT16, dataspace, plist);
    trans_group.createDataSet("ref_start", H5::PredType::NATIVE_UINT16, dataspace, plist);
    trans_group.createDataSet("match_length", H5::PredType::NATIVE_UINT16, dataspace, plist);
    trans_group.createDataSet("codon_start_in_read", H5::PredType::NATIVE_UINT32, dataspace, plist);
    
    // Scoring
    trans_group.createDataSet("alignment_score", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    trans_group.createDataSet("identity", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    trans_group.createDataSet("num_mutations", H5::PredType::NATIVE_UINT8, dataspace, plist);
    
    // Mutation details (variable length)
    trans_group.createDataSet("mutation_positions", H5::PredType::NATIVE_UINT16, dataspace, plist);
    trans_group.createDataSet("ref_aas", H5::PredType::NATIVE_CHAR, dataspace, plist);
    trans_group.createDataSet("query_aas", H5::PredType::NATIVE_CHAR, dataspace, plist);
    trans_group.createDataSet("blosum_scores", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    trans_group.createDataSet("mutation_offsets", H5::PredType::NATIVE_UINT32, dataspace, plist);
    
    // Confidence scores
    trans_group.createDataSet("confidence_score", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    trans_group.createDataSet("is_resistance_mutation", H5::PredType::NATIVE_HBOOL, dataspace, plist);
    
    // Also create a summary group
    file->createGroup("/translated_search/summary");
    H5::Group summary_group = file->openGroup("/translated_search/summary");
    
    // Per-gene resistance counts
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    summary_group.createDataSet("gene_names", str_type, dataspace, plist);
    summary_group.createDataSet("resistance_counts", H5::PredType::NATIVE_UINT32, dataspace, plist);
    summary_group.createDataSet("total_matches", H5::PredType::NATIVE_UINT32, dataspace, plist);
}

void HDF5AlignmentWriter::addAlignmentBatch(const AlignmentResult* h_results, 
                                           size_t num_results,
                                           size_t batch_offset) {
    // Buffer the results
    for (size_t i = 0; i < num_results; i++) {
        const AlignmentResult& result = h_results[i];
        
        read_ids.push_back(result.read_id + batch_offset);
        gene_ids.push_back(result.gene_id);
        species_ids.push_back(result.species_id);
        seq_ids.push_back(result.seq_id);
        alignment_scores.push_back(result.alignment_score);
        identities.push_back(result.identity);
        start_positions.push_back(result.start_pos);
        match_counts.push_back(result.matches);
        num_mutations.push_back(result.num_mutations_detected);
        
        // Handle mutation positions
        uint32_t mutation_offset = mutation_positions.size();
        mutation_offsets.push_back(mutation_offset);
        
        for (int j = 0; j < result.num_mutations_detected; j++) {
            mutation_positions.push_back(result.mutations_detected[j]);
        }
    }
    
    total_alignments_written += num_results;
    
    // Flush if buffer is getting large
    if (read_ids.size() > 100000) {
        flush();
    }
}

void HDF5AlignmentWriter::addTranslatedResults(
    const ProteinMatch* h_matches,
    const uint32_t* h_match_counts,
    size_t num_reads,
    size_t batch_offset
) {
    // Buffer the results for batch writing
    for (size_t read_idx = 0; read_idx < num_reads; read_idx++) {
        uint32_t num_matches = h_match_counts[read_idx];
        const ProteinMatch* read_matches = &h_matches[read_idx * 32];  // max_matches_per_read
        
        for (uint32_t i = 0; i < num_matches; i++) {
            const ProteinMatch& match = read_matches[i];
            
            // Calculate codon start position in original read
            uint32_t codon_start = 0;
            if (match.frame > 0) {
                codon_start = match.query_start * 3 + (match.frame - 1);
            } else {
                // Reverse frame - calculate from end
                // TODO: Implement proper calculation based on read length
                codon_start = 0;
            }
            
            // Calculate confidence score
            float confidence = calculateResistanceConfidence(match);
            bool is_resistance = checkIfResistanceMutation(match);
            
            // Buffer the data
            trans_read_ids.push_back(match.read_id + batch_offset);
            trans_frames.push_back(match.frame);
            trans_protein_ids.push_back(match.protein_id);
            trans_gene_ids.push_back(match.gene_id);
            trans_species_ids.push_back(match.species_id);
            trans_query_starts.push_back(match.query_start);
            trans_ref_starts.push_back(match.ref_start);
            trans_match_lengths.push_back(match.match_length);
            trans_codon_starts.push_back(codon_start);
            trans_alignment_scores.push_back(match.alignment_score);
            trans_identities.push_back(match.identity);
            trans_num_mutations.push_back(match.num_mutations);
            trans_confidence_scores.push_back(confidence);
            trans_is_resistance.push_back(is_resistance);
            
            // Handle mutations
            uint32_t mutation_offset = trans_mutation_positions.size();
            trans_mutation_offsets.push_back(mutation_offset);
            
            for (int m = 0; m < match.num_mutations; m++) {
                trans_mutation_positions.push_back(match.mutation_positions[m]);
                trans_ref_aas.push_back(match.ref_aas[m]);
                trans_query_aas.push_back(match.query_aas[m]);
                trans_blosum_scores.push_back(match.blosum_scores[m]);
            }
        }
    }
    
    total_translated_alignments_written += num_reads;
    
    // Flush if buffer is large
    if (trans_read_ids.size() > 100000) {
        flushTranslatedResults();
    }
}

float HDF5AlignmentWriter::calculateResistanceConfidence(const ProteinMatch& match) {
    // Calculate confidence based on multiple factors
    float confidence = 0.0f;
    
    // Check each mutation
    for (int i = 0; i < match.num_mutations; i++) {
        uint16_t pos = match.mutation_positions[i];
        
        // Check known resistance positions
        bool is_known_position = false;
        if (match.gene_id == 0) {  // gyrA
            if (pos == 83 || pos == 87) is_known_position = true;
        } else if (match.gene_id == 1) {  // parC
            if (pos == 80 || pos == 84) is_known_position = true;
        } else if (match.gene_id == 2) {  // gyrB
            if (pos == 426 || pos == 447) is_known_position = true;
        } else if (match.gene_id == 3) {  // parE
            if (pos == 416 || pos == 420) is_known_position = true;
        }
        
        if (is_known_position) {
            // Calculate confidence components
            float position_weight = 2.0f;  // Known positions weighted higher
            float blosum_component = (match.blosum_scores[i] + 4.0f) / 8.0f;
            float identity_component = match.identity;
            float length_component = std::min(1.0f, match.match_length / 30.0f);
            
            float mutation_confidence = position_weight * blosum_component * 
                                       identity_component * length_component;
            
            confidence = std::max(confidence, mutation_confidence);
        }
    }
    
    return confidence;
}

bool HDF5AlignmentWriter::checkIfResistanceMutation(const ProteinMatch& match) {
    // Check for specific known resistance mutations
    for (int i = 0; i < match.num_mutations; i++) {
        uint16_t pos = match.mutation_positions[i];
        char ref_aa = match.ref_aas[i];
        char query_aa = match.query_aas[i];
        
        // Check specific known mutations
        if (match.gene_id == 0) {  // gyrA
            if (pos == 83 && ref_aa == 'S' && (query_aa == 'L' || query_aa == 'F' || query_aa == 'W')) {
                return true;  // S83L/F/W
            }
            if (pos == 87 && ref_aa == 'D' && (query_aa == 'N' || query_aa == 'G' || query_aa == 'Y')) {
                return true;  // D87N/G/Y
            }
        } else if (match.gene_id == 1) {  // parC
            if (pos == 80 && ref_aa == 'S' && (query_aa == 'I' || query_aa == 'R')) {
                return true;  // S80I/R
            }
            if (pos == 84 && ref_aa == 'E' && (query_aa == 'V' || query_aa == 'K' || query_aa == 'G')) {
                return true;  // E84V/K/G
            }
        }
        // Add other genes as needed
    }
    
    return false;
}

void HDF5AlignmentWriter::flush() {
    if (read_ids.empty()) return;
    
    try {
        H5::Group align_group = file->openGroup("/alignments");
        
        // Helper function to extend and write dataset
        auto extendAndWrite = [&](const std::string& name, 
                                 const void* data, 
                                 const H5::DataType& dtype,
                                 size_t count) {
            H5::DataSet dataset = align_group.openDataSet(name);
            
            // Get current dimensions
            H5::DataSpace filespace = dataset.getSpace();
            hsize_t current_dims[1];
            filespace.getSimpleExtentDims(current_dims);
            
            // Extend dataset
            hsize_t new_dims[1] = {current_dims[0] + count};
            dataset.extend(new_dims);
            
            // Select hyperslab in extended dataset
            filespace = dataset.getSpace();
            hsize_t offset[1] = {current_dims[0]};
            hsize_t write_dims[1] = {count};
            filespace.selectHyperslab(H5S_SELECT_SET, write_dims, offset);
            
            // Create memory dataspace
            H5::DataSpace memspace(1, write_dims);
            
            // Write data
            dataset.write(data, dtype, memspace, filespace);
        };
        
        // Write all arrays
        extendAndWrite("read_id", read_ids.data(), H5::PredType::NATIVE_UINT32, read_ids.size());
        extendAndWrite("gene_id", gene_ids.data(), H5::PredType::NATIVE_UINT32, gene_ids.size());
        extendAndWrite("species_id", species_ids.data(), H5::PredType::NATIVE_UINT32, species_ids.size());
        extendAndWrite("seq_id", seq_ids.data(), H5::PredType::NATIVE_UINT32, seq_ids.size());
        extendAndWrite("alignment_score", alignment_scores.data(), H5::PredType::NATIVE_FLOAT, alignment_scores.size());
        extendAndWrite("identity", identities.data(), H5::PredType::NATIVE_FLOAT, identities.size());
        extendAndWrite("start_pos", start_positions.data(), H5::PredType::NATIVE_UINT16, start_positions.size());
        extendAndWrite("matches", match_counts.data(), H5::PredType::NATIVE_UINT16, match_counts.size());
        extendAndWrite("num_mutations", num_mutations.data(), H5::PredType::NATIVE_UINT8, num_mutations.size());
        
        // Write mutation data
        if (!mutation_positions.empty()) {
            extendAndWrite("mutation_positions", mutation_positions.data(), 
                          H5::PredType::NATIVE_UINT32, mutation_positions.size());
            extendAndWrite("mutation_offsets", mutation_offsets.data(), 
                          H5::PredType::NATIVE_UINT32, mutation_offsets.size());
        }
        
        // Clear buffers
        read_ids.clear();
        gene_ids.clear();
        species_ids.clear();
        seq_ids.clear();
        alignment_scores.clear();
        identities.clear();
        start_positions.clear();
        match_counts.clear();
        num_mutations.clear();
        mutation_positions.clear();
        mutation_offsets.clear();
        
    } catch (H5::Exception& e) {
        std::cerr << "Error writing to HDF5: " << e.getCDetailMsg() << std::endl;
        throw;
    }
}

void HDF5AlignmentWriter::flushTranslatedResults() {
    if (trans_read_ids.empty()) return;
    
    try {
        H5::Group trans_group = file->openGroup("/translated_search");
        
        // Helper function to extend and write dataset
        auto extendAndWrite = [&](const std::string& name, 
                                 const void* data, 
                                 const H5::DataType& dtype,
                                 size_t count) {
            H5::DataSet dataset = trans_group.openDataSet(name);
            
            // Get current dimensions
            H5::DataSpace filespace = dataset.getSpace();
            hsize_t current_dims[1];
            filespace.getSimpleExtentDims(current_dims);
            
            // Extend dataset
            hsize_t new_dims[1] = {current_dims[0] + count};
            dataset.extend(new_dims);
            
            // Select hyperslab in extended dataset
            filespace = dataset.getSpace();
            hsize_t offset[1] = {current_dims[0]};
            hsize_t write_dims[1] = {count};
            filespace.selectHyperslab(H5S_SELECT_SET, write_dims, offset);
            
            // Create memory dataspace
            H5::DataSpace memspace(1, write_dims);
            
            // Write data
            dataset.write(data, dtype, memspace, filespace);
        };
        
        // Write all translated search arrays
        extendAndWrite("read_id", trans_read_ids.data(), H5::PredType::NATIVE_UINT32, trans_read_ids.size());
        extendAndWrite("frame", trans_frames.data(), H5::PredType::NATIVE_INT8, trans_frames.size());
        extendAndWrite("protein_id", trans_protein_ids.data(), H5::PredType::NATIVE_UINT32, trans_protein_ids.size());
        extendAndWrite("gene_id", trans_gene_ids.data(), H5::PredType::NATIVE_UINT32, trans_gene_ids.size());
        extendAndWrite("species_id", trans_species_ids.data(), H5::PredType::NATIVE_UINT32, trans_species_ids.size());
        extendAndWrite("query_start", trans_query_starts.data(), H5::PredType::NATIVE_UINT16, trans_query_starts.size());
        extendAndWrite("ref_start", trans_ref_starts.data(), H5::PredType::NATIVE_UINT16, trans_ref_starts.size());
        extendAndWrite("match_length", trans_match_lengths.data(), H5::PredType::NATIVE_UINT16, trans_match_lengths.size());
        extendAndWrite("codon_start_in_read", trans_codon_starts.data(), H5::PredType::NATIVE_UINT32, trans_codon_starts.size());
        extendAndWrite("alignment_score", trans_alignment_scores.data(), H5::PredType::NATIVE_FLOAT, trans_alignment_scores.size());
        extendAndWrite("identity", trans_identities.data(), H5::PredType::NATIVE_FLOAT, trans_identities.size());
        extendAndWrite("num_mutations", trans_num_mutations.data(), H5::PredType::NATIVE_UINT8, trans_num_mutations.size());
        extendAndWrite("confidence_score", trans_confidence_scores.data(), H5::PredType::NATIVE_FLOAT, trans_confidence_scores.size());
        
        // Write boolean array
        std::vector<uint8_t> resistance_flags(trans_is_resistance.begin(), trans_is_resistance.end());
        extendAndWrite("is_resistance_mutation", resistance_flags.data(), H5::PredType::NATIVE_UINT8, resistance_flags.size());
        
        // Write mutation details
        if (!trans_mutation_positions.empty()) {
            extendAndWrite("mutation_positions", trans_mutation_positions.data(), 
                          H5::PredType::NATIVE_UINT16, trans_mutation_positions.size());
            extendAndWrite("ref_aas", trans_ref_aas.data(), 
                          H5::PredType::NATIVE_CHAR, trans_ref_aas.size());
            extendAndWrite("query_aas", trans_query_aas.data(), 
                          H5::PredType::NATIVE_CHAR, trans_query_aas.size());
            extendAndWrite("blosum_scores", trans_blosum_scores.data(), 
                          H5::PredType::NATIVE_FLOAT, trans_blosum_scores.size());
            extendAndWrite("mutation_offsets", trans_mutation_offsets.data(), 
                          H5::PredType::NATIVE_UINT32, trans_mutation_offsets.size());
        }
        
        // Clear buffers
        trans_read_ids.clear();
        trans_frames.clear();
        trans_protein_ids.clear();
        trans_gene_ids.clear();
        trans_species_ids.clear();
        trans_query_starts.clear();
        trans_ref_starts.clear();
        trans_match_lengths.clear();
        trans_codon_starts.clear();
        trans_alignment_scores.clear();
        trans_identities.clear();
        trans_num_mutations.clear();
        trans_confidence_scores.clear();
        trans_is_resistance.clear();
        trans_mutation_positions.clear();
        trans_ref_aas.clear();
        trans_query_aas.clear();
        trans_blosum_scores.clear();
        trans_mutation_offsets.clear();
        
    } catch (H5::Exception& e) {
        std::cerr << "Error writing translated results to HDF5: " << e.getCDetailMsg() << std::endl;
        throw;
    }
}

void HDF5AlignmentWriter::finalize(const std::string& json_summary_path) {
    flush();  // Write any remaining nucleotide data
    flushTranslatedResults();  // Write any remaining translated data
    
    // Write summary statistics
    H5::Group metadata = file->openGroup("/metadata");
    H5::DataSpace scalar_space(H5S_SCALAR);
    
    H5::Attribute total_attr = metadata.createAttribute("total_alignments", 
                                                       H5::PredType::NATIVE_UINT64, 
                                                       scalar_space);
    total_attr.write(H5::PredType::NATIVE_UINT64, &total_alignments_written);
    
    H5::Attribute trans_attr = metadata.createAttribute("total_translated_alignments", 
                                                       H5::PredType::NATIVE_UINT64, 
                                                       scalar_space);
    trans_attr.write(H5::PredType::NATIVE_UINT64, &total_translated_alignments_written);
    
    // Also write JSON summary for compatibility
    std::ofstream json_file(json_summary_path);
    json_file << "{\n";
    json_file << "  \"hdf5_output\": \"" << filename << "\",\n";
    json_file << "  \"total_alignments\": " << total_alignments_written << ",\n";
    json_file << "  \"total_translated_alignments\": " << total_translated_alignments_written << ",\n";
    json_file << "  \"total_reads_processed\": " << total_reads_processed << "\n";
    json_file << "}\n";
    json_file.close();
}

// Placeholder for optional K-mer screening results method
void HDF5AlignmentWriter::addKmerScreeningResults(const CandidateMatch* h_candidates,
                                                  const uint32_t* h_candidate_counts,
                                                  size_t num_reads,
                                                  size_t batch_offset) {
    // TODO: Implement if needed for ML features
}

// Placeholder for optional SAM export
void HDF5AlignmentWriter::exportToSAM(const std::string& sam_path) {
    // TODO: Implement if needed
}