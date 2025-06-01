#include "hdf5_alignment_writer.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>

HDF5AlignmentWriter::HDF5AlignmentWriter(const std::string& output_path) 
    : filename(output_path), total_alignments_written(0), total_reads_processed(0) {
    
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
    
    // Write pipeline version
    H5::Attribute version_attr = metadata.createAttribute("pipeline_version", str_type, scalar_space);
    version_attr.write(str_type, std::string("0.3.0"));
    
    createDatasets();
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

void HDF5AlignmentWriter::finalize(const std::string& json_summary_path) {
    flush();  // Write any remaining data
    
    // Write summary statistics
    H5::Group metadata = file->openGroup("/metadata");
    H5::DataSpace scalar_space(H5S_SCALAR);
    
    H5::Attribute total_attr = metadata.createAttribute("total_alignments", 
                                                       H5::PredType::NATIVE_UINT64, 
                                                       scalar_space);
    total_attr.write(H5::PredType::NATIVE_UINT64, &total_alignments_written);
    
    // Also write JSON summary for compatibility
    std::ofstream json_file(json_summary_path);
    json_file << "{\n";
    json_file << "  \"hdf5_output\": \"" << filename << "\",\n";
    json_file << "  \"total_alignments\": " << total_alignments_written << ",\n";
    json_file << "  \"total_reads_processed\": " << total_reads_processed << "\n";
    json_file << "}\n";
    json_file.close();
}