// runtime/kernels/genes/hdf5_amr_writer.cpp
#include "hdf5_amr_writer.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>

HDF5AMRWriter::HDF5AMRWriter(const std::string& output_path) 
    : filename(output_path), total_hits_written(0), total_genes_detected(0) {
    
    try {
        file = new H5::H5File(output_path, H5F_ACC_TRUNC);
        
        // Set up compression
        H5::DSetCreatPropList plist;
        hsize_t chunk_dims[1] = {10000};
        plist.setChunk(1, chunk_dims);
        plist.setDeflate(6);
        
    } catch (H5::Exception& e) {
        std::cerr << "Error creating HDF5 file: " << e.getCDetailMsg() << std::endl;
        throw;
    }
}

HDF5AMRWriter::~HDF5AMRWriter() {
    if (file) {
        flush();
        delete file;
    }
}

void HDF5AMRWriter::initialize(const std::string& sample_name,
                              const std::string& dna_fasta_path,
                              const std::string& protein_fasta_path) {
    // Create group structure
    file->createGroup("/metadata");
    file->createGroup("/amr_hits");
    file->createGroup("/coverage_stats");
    file->createGroup("/positional_coverage");
    file->createGroup("/gene_abundance");
    
    // Write metadata
    H5::Group metadata = file->openGroup("/metadata");
    
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace scalar_space(H5S_SCALAR);
    
    H5::Attribute sample_attr = metadata.createAttribute("sample_name", str_type, scalar_space);
    sample_attr.write(str_type, sample_name);
    
    H5::Attribute dna_attr = metadata.createAttribute("dna_database", str_type, scalar_space);
    dna_attr.write(str_type, dna_fasta_path);
    
    H5::Attribute protein_attr = metadata.createAttribute("protein_database", str_type, scalar_space);
    protein_attr.write(str_type, protein_fasta_path);
    
    // Write timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    
    H5::Attribute time_attr = metadata.createAttribute("analysis_timestamp", str_type, scalar_space);
    time_attr.write(str_type, ss.str());
    
    H5::Attribute version_attr = metadata.createAttribute("pipeline_version", str_type, scalar_space);
    version_attr.write(str_type, std::string("AMR-Detection-1.0"));
    
    createDatasets();
}

void HDF5AMRWriter::createDatasets() {
    // Create extensible datasets for AMR hits
    hsize_t initial_dims[1] = {0};
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    H5::DataSpace dataspace(1, initial_dims, max_dims);
    
    H5::DSetCreatPropList plist;
    hsize_t chunk_dims[1] = {10000};
    plist.setChunk(1, chunk_dims);
    plist.setDeflate(6);
    
    // AMR hits datasets
    H5::Group hits_group = file->openGroup("/amr_hits");
    
    hits_group.createDataSet("read_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    hits_group.createDataSet("gene_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    hits_group.createDataSet("gene_name", str_type, dataspace, plist);
    hits_group.createDataSet("gene_family", str_type, dataspace, plist);
    hits_group.createDataSet("drug_class", str_type, dataspace, plist);
    
    hits_group.createDataSet("identity", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    hits_group.createDataSet("coverage", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    hits_group.createDataSet("ref_start", H5::PredType::NATIVE_UINT16, dataspace, plist);
    hits_group.createDataSet("ref_end", H5::PredType::NATIVE_UINT16, dataspace, plist);
    hits_group.createDataSet("frame", H5::PredType::NATIVE_INT8, dataspace, plist);
    hits_group.createDataSet("is_complete_gene", H5::PredType::NATIVE_HBOOL, dataspace, plist);
    hits_group.createDataSet("concordant", H5::PredType::NATIVE_HBOOL, dataspace, plist);
    
    createCoverageDatasets();
}

void HDF5AMRWriter::createCoverageDatasets() {
    hsize_t initial_dims[1] = {0};
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    H5::DataSpace dataspace(1, initial_dims, max_dims);
    
    H5::DSetCreatPropList plist;
    hsize_t chunk_dims[1] = {1000};
    plist.setChunk(1, chunk_dims);
    plist.setDeflate(6);
    
    // Coverage statistics datasets
    H5::Group cov_group = file->openGroup("/coverage_stats");
    
    cov_group.createDataSet("gene_id", H5::PredType::NATIVE_UINT32, dataspace, plist);
    
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    cov_group.createDataSet("gene_name", str_type, dataspace, plist);
    
    cov_group.createDataSet("total_reads", H5::PredType::NATIVE_UINT32, dataspace, plist);
    cov_group.createDataSet("total_bases_mapped", H5::PredType::NATIVE_UINT32, dataspace, plist);
    cov_group.createDataSet("covered_positions", H5::PredType::NATIVE_UINT16, dataspace, plist);
    cov_group.createDataSet("gene_length", H5::PredType::NATIVE_UINT16, dataspace, plist);
    cov_group.createDataSet("percent_coverage", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    cov_group.createDataSet("mean_depth", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    cov_group.createDataSet("rpkm", H5::PredType::NATIVE_FLOAT, dataspace, plist);
    cov_group.createDataSet("tpm", H5::PredType::NATIVE_FLOAT, dataspace, plist);
}

void HDF5AMRWriter::addAMRHits(const std::vector<AMRHit>& hits) {
    // Buffer the hits
    for (const auto& hit : hits) {
        read_ids.push_back(hit.read_id);
        gene_ids.push_back(hit.gene_id);
        gene_names.push_back(std::string(hit.gene_name));
        gene_families.push_back(std::string(hit.gene_family));
        drug_classes.push_back(std::string(hit.drug_class));
        identities.push_back(hit.identity);
        coverages.push_back(hit.coverage);
        ref_starts.push_back(hit.ref_start);
        ref_ends.push_back(hit.ref_end);
        frames.push_back(hit.frame);
        is_complete_genes.push_back(hit.is_complete_gene);
        concordant_flags.push_back(hit.concordant);
    }
    
    total_hits_written += hits.size();
    
    // Flush if buffer is large
    if (read_ids.size() > 100000) {
        flush();
    }
}

void HDF5AMRWriter::addCoverageStats(const std::vector<AMRCoverageStats>& coverage_stats,
                                    const std::vector<AMRGeneEntry>& gene_entries) {
    for (size_t i = 0; i < coverage_stats.size(); i++) {
        const auto& stats = coverage_stats[i];
        if (stats.total_reads > 0) {  // Only store genes with coverage
            const auto& gene = gene_entries[i];
            
            cov_gene_ids.push_back(i);
            cov_gene_names.push_back(std::string(gene.gene_name));
            total_reads.push_back(stats.total_reads);
            total_bases_mapped.push_back(stats.total_bases_mapped);
            covered_positions.push_back(stats.covered_positions);
            gene_lengths.push_back(stats.gene_length);
            percent_coverages.push_back(stats.percent_coverage);
            mean_depths.push_back(stats.mean_depth);
            rpkms.push_back(stats.rpkm);
            tpms.push_back(stats.tpm);
            
            // Store positional coverage if needed
            if (stats.position_counts && stats.gene_length > 0) {
                std::vector<uint32_t> positions(stats.gene_length);
                cudaMemcpy(positions.data(), stats.position_counts, 
                          stats.gene_length * sizeof(uint32_t), cudaMemcpyDeviceToHost);
                position_coverage_map[i] = positions;
            }
            
            total_genes_detected++;
        }
    }
}

void HDF5AMRWriter::flush() {
    if (read_ids.empty() && cov_gene_ids.empty()) return;
    
    try {
        // Helper lambda for extending and writing datasets
        auto extendAndWrite = [&](H5::Group& group, const std::string& name, 
                                 const void* data, const H5::DataType& dtype, size_t count) {
            H5::DataSet dataset = group.openDataSet(name);
            
            H5::DataSpace filespace = dataset.getSpace();
            hsize_t current_dims[1];
            filespace.getSimpleExtentDims(current_dims);
            
            hsize_t new_dims[1] = {current_dims[0] + count};
            dataset.extend(new_dims);
            
            filespace = dataset.getSpace();
            hsize_t offset[1] = {current_dims[0]};
            hsize_t write_dims[1] = {count};
            filespace.selectHyperslab(H5S_SELECT_SET, write_dims, offset);
            
            H5::DataSpace memspace(1, write_dims);
            dataset.write(data, dtype, memspace, filespace);
        };
        
        // Write AMR hits
        if (!read_ids.empty()) {
            H5::Group hits_group = file->openGroup("/amr_hits");
            
            extendAndWrite(hits_group, "read_id", read_ids.data(), 
                          H5::PredType::NATIVE_UINT32, read_ids.size());
            extendAndWrite(hits_group, "gene_id", gene_ids.data(), 
                          H5::PredType::NATIVE_UINT32, gene_ids.size());
            
            // Write string data
            H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
            
            // Convert string vectors to char* arrays for HDF5
            std::vector<const char*> gene_name_ptrs;
            std::vector<const char*> gene_family_ptrs;
            std::vector<const char*> drug_class_ptrs;
            for (size_t i = 0; i < gene_names.size(); i++) {
                gene_name_ptrs.push_back(gene_names[i].c_str());
                gene_family_ptrs.push_back(gene_families[i].c_str());
                drug_class_ptrs.push_back(drug_classes[i].c_str());
            }
            
            extendAndWrite(hits_group, "gene_name", gene_name_ptrs.data(), str_type, gene_names.size());
            extendAndWrite(hits_group, "gene_family", gene_family_ptrs.data(), str_type, gene_families.size());
            extendAndWrite(hits_group, "drug_class", drug_class_ptrs.data(), str_type, drug_classes.size());
            
            extendAndWrite(hits_group, "identity", identities.data(), 
                          H5::PredType::NATIVE_FLOAT, identities.size());
            extendAndWrite(hits_group, "coverage", coverages.data(), 
                          H5::PredType::NATIVE_FLOAT, coverages.size());
            extendAndWrite(hits_group, "ref_start", ref_starts.data(), 
                          H5::PredType::NATIVE_UINT16, ref_starts.size());
            extendAndWrite(hits_group, "ref_end", ref_ends.data(), 
                          H5::PredType::NATIVE_UINT16, ref_ends.size());
            extendAndWrite(hits_group, "frame", frames.data(), 
                          H5::PredType::NATIVE_INT8, frames.size());
            
            // Convert bool to uint8 for HDF5
            std::vector<uint8_t> complete_flags(is_complete_genes.begin(), is_complete_genes.end());
            std::vector<uint8_t> concordant_uint(concordant_flags.begin(), concordant_flags.end());
            
            extendAndWrite(hits_group, "is_complete_gene", complete_flags.data(), 
                          H5::PredType::NATIVE_UINT8, complete_flags.size());
            extendAndWrite(hits_group, "concordant", concordant_uint.data(), 
                          H5::PredType::NATIVE_UINT8, concordant_uint.size());
            
            // Clear buffers
            read_ids.clear();
            gene_ids.clear();
            gene_names.clear();
            gene_families.clear();
            drug_classes.clear();
            identities.clear();
            coverages.clear();
            ref_starts.clear();
            ref_ends.clear();
            frames.clear();
            is_complete_genes.clear();
            concordant_flags.clear();
        }
        
        // Write coverage statistics
        if (!cov_gene_ids.empty()) {
            H5::Group cov_group = file->openGroup("/coverage_stats");
            
            extendAndWrite(cov_group, "gene_id", cov_gene_ids.data(), 
                          H5::PredType::NATIVE_UINT32, cov_gene_ids.size());
            
            // String data
            H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
            std::vector<const char*> cov_gene_name_ptrs;
            for (const auto& name : cov_gene_names) {
                cov_gene_name_ptrs.push_back(name.c_str());
            }
            extendAndWrite(cov_group, "gene_name", cov_gene_name_ptrs.data(), 
                          str_type, cov_gene_names.size());
            
            extendAndWrite(cov_group, "total_reads", total_reads.data(), 
                          H5::PredType::NATIVE_UINT32, total_reads.size());
            extendAndWrite(cov_group, "total_bases_mapped", total_bases_mapped.data(), 
                          H5::PredType::NATIVE_UINT32, total_bases_mapped.size());
            extendAndWrite(cov_group, "covered_positions", covered_positions.data(), 
                          H5::PredType::NATIVE_UINT16, covered_positions.size());
            extendAndWrite(cov_group, "gene_length", gene_lengths.data(), 
                          H5::PredType::NATIVE_UINT16, gene_lengths.size());
            extendAndWrite(cov_group, "percent_coverage", percent_coverages.data(), 
                          H5::PredType::NATIVE_FLOAT, percent_coverages.size());
            extendAndWrite(cov_group, "mean_depth", mean_depths.data(), 
                          H5::PredType::NATIVE_FLOAT, mean_depths.size());
            extendAndWrite(cov_group, "rpkm", rpkms.data(), 
                          H5::PredType::NATIVE_FLOAT, rpkms.size());
            extendAndWrite(cov_group, "tpm", tpms.data(), 
                          H5::PredType::NATIVE_FLOAT, tpms.size());
            
            // Write positional coverage data
            if (!position_coverage_map.empty()) {
                H5::Group pos_group = file->openGroup("/positional_coverage");
                
                for (const auto& [gene_id, positions] : position_coverage_map) {
                    std::string dataset_name = "gene_" + std::to_string(gene_id);
                    
                    hsize_t dims[1] = {positions.size()};
                    H5::DataSpace dataspace(1, dims);
                    
                    H5::DSetCreatPropList plist;
                    hsize_t chunk_dims[1] = {std::min((size_t)1000, positions.size())};
                    plist.setChunk(1, chunk_dims);
                    plist.setDeflate(6);
                    
                    H5::DataSet dataset = pos_group.createDataSet(dataset_name, 
                                                                 H5::PredType::NATIVE_UINT32, 
                                                                 dataspace, plist);
                    dataset.write(positions.data(), H5::PredType::NATIVE_UINT32);
                }
                position_coverage_map.clear();
            }
            
            // Clear coverage buffers
            cov_gene_ids.clear();
            cov_gene_names.clear();
            total_reads.clear();
            total_bases_mapped.clear();
            covered_positions.clear();
            gene_lengths.clear();
            percent_coverages.clear();
            mean_depths.clear();
            rpkms.clear();
            tpms.clear();
        }
        
    } catch (H5::Exception& e) {
        std::cerr << "Error writing to HDF5: " << e.getCDetailMsg() << std::endl;
        throw;
    }
}

void HDF5AMRWriter::finalize(const std::string& json_summary_path) {
    flush();
    
    // Write summary statistics
    H5::Group metadata = file->openGroup("/metadata");
    H5::DataSpace scalar_space(H5S_SCALAR);
    
    H5::Attribute hits_attr = metadata.createAttribute("total_amr_hits", 
                                                       H5::PredType::NATIVE_UINT64, 
                                                       scalar_space);
    hits_attr.write(H5::PredType::NATIVE_UINT64, &total_hits_written);
    
    H5::Attribute genes_attr = metadata.createAttribute("total_genes_detected", 
                                                        H5::PredType::NATIVE_UINT64, 
                                                        scalar_space);
    genes_attr.write(H5::PredType::NATIVE_UINT64, &total_genes_detected);
    
    // Write JSON summary
    std::ofstream json_file(json_summary_path);
    json_file << "{\n";
    json_file << "  \"hdf5_output\": \"" << filename << "\",\n";
    json_file << "  \"total_amr_hits\": " << total_hits_written << ",\n";
    json_file << "  \"total_genes_detected\": " << total_genes_detected << "\n";
    json_file << "}\n";
    json_file.close();
}
