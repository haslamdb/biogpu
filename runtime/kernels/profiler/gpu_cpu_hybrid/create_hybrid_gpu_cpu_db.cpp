#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <algorithm>
#include <cctype>

// Database preparation tool for hybrid profiler
// Converts FASTA genome collections to memory-mapped binary format

struct OrganismRecord {
    uint32_t taxonomy_id;
    std::string name;
    std::string sequence;
    std::string taxonomy_path;  // Full taxonomic lineage
    uint32_t taxon_level;       // 0=strain, 1=species, 2=genus, etc.
    float gc_content;
};

class DatabaseBuilder {
private:
    std::vector<OrganismRecord> organisms;
    std::unordered_map<std::string, uint32_t> species_to_taxid;
    
public:
    DatabaseBuilder() {
    }
    
    // Load taxonomy mapping from file if it exists
    void load_taxonomy_mapping(const std::string& mapping_file) {
        std::ifstream file(mapping_file);
        if (!file.is_open()) {
            std::cout << "No taxonomy mapping file found at: " << mapping_file << std::endl;
            std::cout << "Will use auto-generated taxonomy IDs" << std::endl;
            return;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            size_t tab_pos = line.find('\t');
            if (tab_pos != std::string::npos) {
                std::string taxid_str = line.substr(0, tab_pos);
                std::string species_name = line.substr(tab_pos + 1);
                
                try {
                    uint32_t taxid = std::stoul(taxid_str);
                    species_to_taxid[species_name] = taxid;
                } catch (...) {
                    std::cerr << "Invalid taxonomy ID in line: " << line << std::endl;
                }
            }
        }
        
        std::cout << "Loaded " << species_to_taxid.size() << " taxonomy mappings" << std::endl;
        file.close();
    }
    
private:
    
public:
    void load_fasta_directory(const std::string& fasta_dir) {
        std::cout << "Loading FASTA files from: " << fasta_dir << std::endl;
        
        for (const auto& entry : std::filesystem::directory_iterator(fasta_dir)) {
            if (entry.path().extension() == ".fasta" || 
                entry.path().extension() == ".fa" ||
                entry.path().extension() == ".fna") {
                
                std::cout << "Processing: " << entry.path().filename() << std::endl;
                load_fasta_file(entry.path().string());
            }
        }
        
        std::cout << "Loaded " << organisms.size() << " organisms" << std::endl;
    }
    
    void load_fasta_file(const std::string& fasta_path) {
        std::ifstream file(fasta_path);
        if (!file.is_open()) {
            std::cerr << "Cannot open: " << fasta_path << std::endl;
            return;
        }
        
        OrganismRecord current_organism;
        std::string line;
        bool in_sequence = false;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Save previous organism if we have one
                if (in_sequence && !current_organism.sequence.empty()) {
                    finalize_organism(current_organism);
                    organisms.push_back(current_organism);
                }
                
                // Parse header with filename context
                parse_fasta_header(line, current_organism, fasta_path);
                current_organism.sequence.clear();
                in_sequence = true;
                
            } else if (in_sequence) {
                // Add to sequence, removing whitespace
                for (char c : line) {
                    if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
                        current_organism.sequence += std::toupper(c);
                    }
                }
            }
        }
        
        // Don't forget the last organism
        if (in_sequence && !current_organism.sequence.empty()) {
            finalize_organism(current_organism);
            organisms.push_back(current_organism);
        }
        
        file.close();
    }
    
    void parse_fasta_header(const std::string& header, OrganismRecord& organism, const std::string& fasta_path) {
        // Parse FASTA header to extract organism info
        // Format examples:
        // >NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
        // >GCF_000005825.2_ASM582v2_genomic
        
        std::string clean_header = header.substr(1);  // Remove '>'
        
        // First try to extract species name from filename
        std::filesystem::path filepath(fasta_path);
        std::string filename = filepath.stem().string();  // filename without extension
        
        // Extract species name from filename pattern: Species_name_ID.fasta
        std::string species_from_filename;
        size_t last_underscore = filename.find_last_of('_');
        if (last_underscore != std::string::npos) {
            // Check if the part after last underscore is numeric (ID)
            std::string potential_id = filename.substr(last_underscore + 1);
            bool is_numeric = !potential_id.empty() && 
                             std::all_of(potential_id.begin(), potential_id.end(), ::isdigit);
            
            if (is_numeric) {
                species_from_filename = filename.substr(0, last_underscore);
            }
        }
        
        // Try to get taxonomy ID from our mapping
        if (!species_from_filename.empty() && species_to_taxid.count(species_from_filename) > 0) {
            organism.taxonomy_id = species_to_taxid[species_from_filename];
        } else {
            // Extract taxonomy ID if present in header
            size_t taxid_pos = clean_header.find("taxid:");
            if (taxid_pos != std::string::npos) {
                std::string taxid_str = clean_header.substr(taxid_pos + 6);
                size_t space_pos = taxid_str.find(' ');
                if (space_pos != std::string::npos) {
                    taxid_str = taxid_str.substr(0, space_pos);
                }
                try {
                    organism.taxonomy_id = std::stoul(taxid_str);
                } catch (...) {
                    organism.taxonomy_id = 0;
                }
            } else {
                // Generate a hash-based ID if no taxonomy ID found
                organism.taxonomy_id = std::hash<std::string>{}(clean_header) % 1000000;
            }
        }
        
        // Extract organism name and build taxonomy path
        size_t first_space = clean_header.find(' ');
        if (first_space != std::string::npos) {
            organism.name = clean_header.substr(first_space + 1);
            
            // Clean up common suffixes
            size_t comma_pos = organism.name.find(',');
            if (comma_pos != std::string::npos) {
                organism.name = organism.name.substr(0, comma_pos);
            }
        } else {
            organism.name = clean_header;
        }
        
        // Build taxonomy path
        organism.taxonomy_path = build_taxonomy_path(organism.name);
        
        // Truncate very long names
        if (organism.name.length() > 200) {
            organism.name = organism.name.substr(0, 200);
        }
    }
    
    void finalize_organism(OrganismRecord& organism) {
        // Calculate GC content
        int gc_count = 0;
        int total_bases = 0;
        
        for (char base : organism.sequence) {
            if (base == 'G' || base == 'C' || base == 'g' || base == 'c') {
                gc_count++;
            }
            if (base != 'N' && base != 'n') {
                total_bases++;
            }
        }
        
        organism.gc_content = total_bases > 0 ? (float)gc_count / total_bases * 100.0f : 0.0f;
        
        // Determine taxonomic level (basic heuristic)
        organism.taxon_level = determine_taxonomic_level(organism.name);
        
        
        // Quality control
        if (organism.sequence.length() < 1000) {
            std::cout << "Warning: Very short sequence for " << organism.name 
                      << " (" << organism.sequence.length() << " bp)" << std::endl;
        }
        
        if (organism.sequence.length() > 50000000) {  // 50 Mbp
            std::cout << "Warning: Very large genome for " << organism.name 
                      << " (" << organism.sequence.length() << " bp)" << std::endl;
        }
    }
    
    uint32_t determine_taxonomic_level(const std::string& organism_name) {
        // Basic heuristic to determine taxonomic level from name
        // 0=strain, 1=species, 2=genus, 3=family, etc.
        
        if (organism_name.find("strain") != std::string::npos ||
            organism_name.find("str.") != std::string::npos) {
            return 0;  // Strain level
        }
        
        // Count words - species names typically have 2+ words
        int word_count = 1;
        for (char c : organism_name) {
            if (c == ' ') word_count++;
        }
        
        if (word_count >= 2) {
            return 1;  // Species level
        } else {
            return 2;  // Genus level
        }
    }
    
    std::string build_taxonomy_path(const std::string& organism_name) {
        // Build a basic taxonomy path from organism name
        // In a real implementation, this would query NCBI taxonomy
        
        std::string taxonomy = "Bacteria"; // Default to Bacteria
        
        // Add some basic classifications
        if (organism_name.find("Escherichia") != std::string::npos) {
            taxonomy = "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia";
        } else if (organism_name.find("Staphylococcus") != std::string::npos) {
            taxonomy = "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus";
        } else if (organism_name.find("Pseudomonas") != std::string::npos) {
            taxonomy = "Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas";
        } else if (organism_name.find("Klebsiella") != std::string::npos) {
            taxonomy = "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella";
        } else {
            // Extract genus from organism name
            size_t space_pos = organism_name.find(' ');
            if (space_pos != std::string::npos) {
                std::string genus = organism_name.substr(0, space_pos);
                taxonomy = "Bacteria;Unknown_phylum;Unknown_class;Unknown_order;Unknown_family;" + genus;
            }
        }
        
        return taxonomy;
    }
    
    
    void write_binary_database(const std::string& output_path) {
        std::cout << "Writing comprehensive binary database to: " << output_path << std::endl;
        
        std::ofstream out(output_path, std::ios::binary);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot create output file: " + output_path);
        }
        
        // Write header with explicit byte order
        const char magic_bytes[] = {'B', 'I', 'O', 'G'};  // Write as chars to ensure order
        out.write(magic_bytes, 4);
        
        uint32_t version = 2;         // Updated version for comprehensive format
        uint32_t num_organisms = organisms.size();
        
        out.write(reinterpret_cast<const char*>(&version), sizeof(version));
        out.write(reinterpret_cast<const char*>(&num_organisms), sizeof(num_organisms));
        
        // Calculate offsets for metadata and sequences
        uint64_t current_offset = 4 + sizeof(version) + sizeof(num_organisms);
        
        // First, write all organism metadata
        std::vector<uint64_t> sequence_offsets;
        uint64_t sequence_start_offset = current_offset;
        
        // Calculate size needed for metadata
        for (const auto& organism : organisms) {
            sequence_start_offset += sizeof(uint32_t);  // taxonomy_id
            sequence_start_offset += sizeof(uint64_t);  // genome_offset
            sequence_start_offset += sizeof(uint64_t);  // genome_size
            sequence_start_offset += sizeof(uint32_t);  // taxon_level
            sequence_start_offset += sizeof(float);     // gc_content
            sequence_start_offset += sizeof(uint16_t);  // name_length
            sequence_start_offset += organism.name.length();
            sequence_start_offset += sizeof(uint16_t);  // taxonomy_path_length
            sequence_start_offset += organism.taxonomy_path.length();
        }
        
        current_offset = sequence_start_offset;
        
        // Write organism metadata
        for (const auto& organism : organisms) {
            // Basic info
            out.write(reinterpret_cast<const char*>(&organism.taxonomy_id), sizeof(organism.taxonomy_id));
            
            sequence_offsets.push_back(current_offset);
            uint64_t sequence_length = organism.sequence.length();
            out.write(reinterpret_cast<const char*>(&current_offset), sizeof(current_offset));
            out.write(reinterpret_cast<const char*>(&sequence_length), sizeof(sequence_length));
            out.write(reinterpret_cast<const char*>(&organism.taxon_level), sizeof(organism.taxon_level));
            out.write(reinterpret_cast<const char*>(&organism.gc_content), sizeof(organism.gc_content));
            
            // Variable length strings
            uint16_t name_length = organism.name.length();
            out.write(reinterpret_cast<const char*>(&name_length), sizeof(name_length));
            out.write(organism.name.c_str(), name_length);
            
            uint16_t taxonomy_length = organism.taxonomy_path.length();
            out.write(reinterpret_cast<const char*>(&taxonomy_length), sizeof(taxonomy_length));
            out.write(organism.taxonomy_path.c_str(), taxonomy_length);
            
            current_offset += sequence_length;
        }
        
        // Write sequences
        for (const auto& organism : organisms) {
            out.write(organism.sequence.c_str(), organism.sequence.length());
        }
        
        out.close();
        
        // Write summary statistics
        write_comprehensive_database_summary(output_path + ".summary");
        
        std::cout << "Comprehensive database written successfully!" << std::endl;
        std::cout << "Total size: " << std::filesystem::file_size(output_path) / (1024*1024) << " MB" << std::endl;
    }
    
    void write_comprehensive_database_summary(const std::string& summary_path) {
        std::ofstream summary(summary_path);
        
        summary << "Comprehensive Metagenomic Database Summary\n";
        summary << "==========================================\n\n";
        
        summary << "Total organisms: " << organisms.size() << "\n";
        
        // Count by taxonomic level
        std::unordered_map<uint32_t, int> level_counts;
        std::unordered_map<uint32_t, std::string> level_names = {
            {0, "Strain"}, {1, "Species"}, {2, "Genus"}, {3, "Family"}, {4, "Order"}
        };
        
        uint64_t total_bases = 0;
        float total_gc = 0.0f;
        
        std::unordered_map<std::string, int> phylum_counts;
        
        for (const auto& org : organisms) {
            level_counts[org.taxon_level]++;
            total_bases += org.sequence.length();
            total_gc += org.gc_content;
            
            // Extract phylum from taxonomy path
            std::vector<std::string> taxa;
            std::stringstream ss(org.taxonomy_path);
            std::string taxon;
            while (std::getline(ss, taxon, ';')) {
                taxa.push_back(taxon);
            }
            if (taxa.size() > 1) {
                phylum_counts[taxa[1]]++;
            }
        }
        
        summary << "\nTaxonomic Level Distribution:\n";
        for (const auto& [level, name] : level_names) {
            if (level_counts[level] > 0) {
                summary << "- " << name << ": " << level_counts[level] << "\n";
            }
        }
        
        summary << "\nPhylum Distribution:\n";
        for (const auto& [phylum, count] : phylum_counts) {
            summary << "- " << phylum << ": " << count << "\n";
        }
        
        summary << "\nSequence Statistics:\n";
        summary << "- Total sequence data: " << total_bases / 1000000 << " Mbp\n";
        summary << "- Average GC content: " << std::fixed << std::setprecision(1) 
                << (total_gc / organisms.size()) << "%\n";
        
        summary << "\nLargest genomes:\n";
        auto sorted_organisms = organisms;
        std::sort(sorted_organisms.begin(), sorted_organisms.end(),
                 [](const auto& a, const auto& b) { return a.sequence.length() > b.sequence.length(); });
        
        for (int i = 0; i < std::min(20, (int)sorted_organisms.size()); i++) {
            summary << "- " << sorted_organisms[i].name 
                    << ": " << std::fixed << std::setprecision(2) 
                    << sorted_organisms[i].sequence.length() / 1000000.0 << " Mbp "
                    << "(" << sorted_organisms[i].gc_content << "% GC)\n";
        }
        
        summary.close();
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input_fasta_directory> <output_database>\n";
        std::cout << "\nPrepares a binary database for the hybrid metagenomic profiler\n";
        std::cout << "Input should be a directory containing FASTA files of microbial genomes\n";
        std::cout << "\nExample:\n";
        std::cout << "  " << argv[0] << " /data/refseq_genomes /data/clinical_microbes.db\n";
        return 1;
    }
    
    std::string input_dir = argv[1];
    std::string output_db = argv[2];
    
    try {
        std::cout << "BioGPU Database Preparation Tool\n";
        std::cout << "================================\n";
        
        DatabaseBuilder builder;
        
        // Check for taxonomy mapping file in the parent directory of the FASTA directory
        std::filesystem::path input_path(input_dir);
        std::filesystem::path parent_dir = input_path.parent_path();
        std::filesystem::path taxon_mapping_file = parent_dir / "taxon_mapping.txt";
        
        if (std::filesystem::exists(taxon_mapping_file)) {
            std::cout << "Found taxonomy mapping file: " << taxon_mapping_file << std::endl;
            builder.load_taxonomy_mapping(taxon_mapping_file.string());
        } else {
            // Also check in the input directory itself
            taxon_mapping_file = input_path / "taxon_mapping.txt";
            if (std::filesystem::exists(taxon_mapping_file)) {
                std::cout << "Found taxonomy mapping file: " << taxon_mapping_file << std::endl;
                builder.load_taxonomy_mapping(taxon_mapping_file.string());
            }
        }
        
        builder.load_fasta_directory(input_dir);
        builder.write_binary_database(output_db);
        
        std::cout << "\nDatabase preparation complete!\n";
        std::cout << "Database file: " << output_db << "\n";
        std::cout << "Summary file: " << output_db << ".summary\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}