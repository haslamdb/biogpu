#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include <map>
#include <chrono>

// Kraken2-inspired minimizer-based database for microbial profiling
// Clean, homogeneous approach without resistance detection
// Designed for dynamic updates and easy expansion

class MinimizerDatabaseBuilder {
private:
    struct MinimizerEntry {
        uint64_t minimizer_hash;
        uint32_t taxonomy_id;
        uint8_t uniqueness_score;  // 0-255, higher = more unique
        
        // Comparison operator for sorting
        bool operator<(const MinimizerEntry& other) const {
            if (minimizer_hash != other.minimizer_hash)
                return minimizer_hash < other.minimizer_hash;
            return taxonomy_id < other.taxonomy_id;
        }
    } __attribute__((packed));
    
    struct OrganismInfo {
        uint32_t taxonomy_id;
        std::string name;
        std::string taxonomy_path;
        uint32_t taxon_level;
        float gc_content;
        uint64_t genome_size;
        uint32_t minimizer_count;  // For tracking representation
    };
    
    // Kraken2-style parameters
    int k = 35;           // K-mer size (Kraken2 default)
    int m = 31;           // Minimizer length
    int window_size = 4;  // Minimizer window size
    int stride = 1;       // Process every position (no skipping)
    
    std::vector<OrganismInfo> organisms;
    std::map<uint64_t, std::vector<MinimizerEntry>> minimizer_index;
    std::unordered_map<uint64_t, int> minimizer_frequency;  // For uniqueness scoring
    
    // Dynamic build support
    std::unordered_set<std::string> processed_files;  // Track processed files
    bool incremental_mode = false;
    
public:
    MinimizerDatabaseBuilder(int k_size = 35, int min_size = 31, bool incremental = false) 
        : k(k_size), m(min_size), incremental_mode(incremental) {
        window_size = k - m + 1;
        std::cout << "Minimizer database builder: k=" << k 
                  << ", m=" << m << ", window=" << window_size 
                  << ", incremental=" << (incremental ? "yes" : "no") << std::endl;
    }
    
    // Load existing database for incremental updates
    void load_existing_database(const std::string& db_path) {
        if (!std::filesystem::exists(db_path)) {
            std::cout << "No existing database found, starting fresh" << std::endl;
            return;
        }
        
        std::cout << "Loading existing database for incremental update..." << std::endl;
        
        std::ifstream in(db_path, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open existing database file");
        }
        
        // Read header
        uint32_t magic, version, k_size, m_size, num_organisms;
        uint64_t num_minimizers;
        
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        in.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
        in.read(reinterpret_cast<char*>(&num_organisms), sizeof(num_organisms));
        in.read(reinterpret_cast<char*>(&num_minimizers), sizeof(num_minimizers));
        
        if (magic != 0x4D494E49) {  // "MINI" for minimizer database
            throw std::runtime_error("Invalid database format");
        }
        
        // Verify parameters match
        if (k_size != k || m_size != m) {
            std::cout << "Warning: Database parameters (k=" << k_size << ", m=" << m_size 
                      << ") differ from current (k=" << k << ", m=" << m << ")" << std::endl;
        }
        
        // Read organism metadata
        for (uint32_t i = 0; i < num_organisms; i++) {
            OrganismInfo org;
            
            in.read(reinterpret_cast<char*>(&org.taxonomy_id), sizeof(org.taxonomy_id));
            in.read(reinterpret_cast<char*>(&org.taxon_level), sizeof(org.taxon_level));
            in.read(reinterpret_cast<char*>(&org.gc_content), sizeof(org.gc_content));
            in.read(reinterpret_cast<char*>(&org.genome_size), sizeof(org.genome_size));
            in.read(reinterpret_cast<char*>(&org.minimizer_count), sizeof(org.minimizer_count));
            
            uint16_t name_length;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            org.name.resize(name_length);
            in.read(&org.name[0], name_length);
            
            uint16_t taxonomy_length;
            in.read(reinterpret_cast<char*>(&taxonomy_length), sizeof(taxonomy_length));
            org.taxonomy_path.resize(taxonomy_length);
            in.read(&org.taxonomy_path[0], taxonomy_length);
            
            organisms.push_back(org);
        }
        
        // Read minimizer index
        for (uint64_t i = 0; i < num_minimizers; i++) {
            uint64_t hash;
            uint32_t num_entries;
            
            in.read(reinterpret_cast<char*>(&hash), sizeof(hash));
            in.read(reinterpret_cast<char*>(&num_entries), sizeof(num_entries));
            
            std::vector<MinimizerEntry> entries;
            for (uint32_t j = 0; j < num_entries; j++) {
                MinimizerEntry entry;
                in.read(reinterpret_cast<char*>(&entry), sizeof(entry));
                entries.push_back(entry);
                
                // Update frequency count
                minimizer_frequency[hash]++;
            }
            
            minimizer_index[hash] = entries;
        }
        
        in.close();
        
        std::cout << "Loaded existing database: " << organisms.size() << " organisms, " 
                  << minimizer_index.size() << " unique minimizers" << std::endl;
    }
    
    void load_fasta_directory(const std::string& fasta_dir) {
        std::cout << "Scanning FASTA files from: " << fasta_dir << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        for (const auto& entry : std::filesystem::directory_iterator(fasta_dir)) {
            if (entry.path().extension() == ".fasta" || 
                entry.path().extension() == ".fa" ||
                entry.path().extension() == ".fna") {
                
                std::string file_path = entry.path().string();
                
                // Skip already processed files in incremental mode
                if (incremental_mode && processed_files.count(file_path)) {
                    std::cout << "Skipping already processed: " << entry.path().filename() << std::endl;
                    continue;
                }
                
                std::cout << "Processing: " << entry.path().filename() << std::endl;
                load_fasta_file(file_path);
                processed_files.insert(file_path);
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "Loaded " << organisms.size() << " total organisms in " 
                  << duration.count() << " seconds" << std::endl;
        
        // Calculate minimizer uniqueness scores
        calculate_uniqueness_scores();
        
        // Build final optimized index
        build_optimized_index();
    }
    
private:
    void load_fasta_file(const std::string& fasta_path) {
        std::ifstream file(fasta_path);
        if (!file.is_open()) {
            std::cerr << "Cannot open: " << fasta_path << std::endl;
            return;
        }
        
        OrganismInfo current_organism;
        std::string sequence;
        std::string line;
        bool in_sequence = false;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Process previous organism if we have one
                if (in_sequence && !sequence.empty()) {
                    process_organism_sequence(current_organism, sequence);
                    organisms.push_back(current_organism);
                }
                
                // Parse new header
                parse_fasta_header(line, current_organism, fasta_path);
                sequence.clear();
                in_sequence = true;
                
            } else if (in_sequence) {
                // Clean and add sequence
                for (char c : line) {
                    if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
                        sequence += std::toupper(c);
                    }
                }
            }
        }
        
        // Process last organism
        if (in_sequence && !sequence.empty()) {
            process_organism_sequence(current_organism, sequence);
            organisms.push_back(current_organism);
        }
        
        file.close();
    }
    
    void parse_fasta_header(const std::string& header, OrganismInfo& organism, 
                           const std::string& fasta_path) {
        std::string clean_header = header.substr(1);  // Remove '>'
        
        // Extract organism name
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
        
        // Generate consistent taxonomy ID from name hash
        organism.taxonomy_id = std::hash<std::string>{}(organism.name) % 1000000;
        
        // Build basic taxonomy path
        organism.taxonomy_path = build_taxonomy_path(organism.name);
        organism.taxon_level = determine_taxonomic_level(organism.name);
        
        // Truncate very long names
        if (organism.name.length() > 200) {
            organism.name = organism.name.substr(0, 200);
        }
    }
    
    void process_organism_sequence(OrganismInfo& organism, const std::string& sequence) {
        organism.genome_size = sequence.length();
        
        // Calculate GC content
        calculate_gc_content(organism, sequence);
        
        // Extract minimizers using consistent approach
        organism.minimizer_count = extract_minimizers_from_sequence(organism, sequence);
        
        std::cout << "  " << organism.name << ": " << sequence.length() 
                  << " bp, " << std::fixed << std::setprecision(1) 
                  << organism.gc_content << "% GC, " 
                  << organism.minimizer_count << " minimizers" << std::endl;
    }
    
    uint32_t extract_minimizers_from_sequence(const OrganismInfo& organism, 
                                            const std::string& sequence) {
        if (sequence.length() < k) return 0;
        
        uint32_t minimizer_count = 0;
        
        // Use consistent minimizer extraction - every position with stride=1
        for (size_t i = 0; i <= sequence.length() - k; i += stride) {
            // Extract k-mer window
            std::string kmer_window = sequence.substr(i, k);
            
            // Find minimizer in window (smallest hash)
            uint64_t min_hash = UINT64_MAX;
            bool valid_minimizer = false;
            
            for (int j = 0; j <= k - m; j++) {
                std::string minimizer = kmer_window.substr(j, m);
                uint64_t hash = hash_sequence(minimizer);
                
                if (hash != UINT64_MAX && hash < min_hash) {
                    min_hash = hash;
                    valid_minimizer = true;
                }
            }
            
            if (valid_minimizer) {
                // Count frequency for later uniqueness calculation
                minimizer_frequency[min_hash]++;
                
                // Store minimizer with taxonomy
                MinimizerEntry entry{
                    min_hash,
                    organism.taxonomy_id,
                    0   // Uniqueness score calculated later
                };
                
                minimizer_index[min_hash].push_back(entry);
                minimizer_count++;
            }
        }
        
        return minimizer_count;
    }
    
    uint64_t hash_sequence(const std::string& seq) {
        // Canonical k-mer hashing (forward and reverse complement)
        uint64_t forward_hash = 0;
        uint64_t reverse_hash = 0;
        
        bool valid = true;
        
        for (size_t i = 0; i < seq.length(); i++) {
            int base = encode_base(seq[i]);
            if (base == -1) {
                valid = false;
                break;
            }
            
            forward_hash = (forward_hash << 2) | base;
            reverse_hash = (reverse_hash >> 2) | (((uint64_t)(3 ^ base)) << (2 * (seq.length() - 1)));
        }
        
        return valid ? std::min(forward_hash, reverse_hash) : UINT64_MAX;
    }
    
    int encode_base(char base) {
        switch (base) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    }
    
    void calculate_gc_content(OrganismInfo& organism, const std::string& sequence) {
        int gc_count = 0;
        int total_bases = 0;
        
        for (char base : sequence) {
            if (base == 'G' || base == 'C') gc_count++;
            if (base != 'N') total_bases++;
        }
        
        organism.gc_content = total_bases > 0 ? (float)gc_count / total_bases * 100.0f : 0.0f;
    }
    
    void calculate_uniqueness_scores() {
        std::cout << "Calculating minimizer uniqueness scores..." << std::endl;
        
        // Calculate uniqueness score for each minimizer
        for (auto& [hash, entries] : minimizer_index) {
            int frequency = minimizer_frequency[hash];
            
            // Uniqueness score: higher for rarer minimizers
            uint8_t uniqueness = 0;
            if (frequency == 1) {
                uniqueness = 255;  // Completely unique
            } else if (frequency <= 5) {
                uniqueness = 200;  // Very unique
            } else if (frequency <= 20) {
                uniqueness = 150;  // Moderately unique
            } else if (frequency <= 100) {
                uniqueness = 100;  // Somewhat unique
            } else {
                uniqueness = 50;   // Common
            }
            
            // Update all entries for this minimizer
            for (auto& entry : entries) {
                entry.uniqueness_score = uniqueness;
            }
        }
        
        std::cout << "Uniqueness scores calculated for " << minimizer_index.size() 
                  << " unique minimizers" << std::endl;
    }
    
    void build_optimized_index() {
        std::cout << "Building optimized index..." << std::endl;
        
        // Filter out overly common minimizers to reduce database size
        auto it = minimizer_index.begin();
        size_t removed_common = 0;
        
        while (it != minimizer_index.end()) {
            int frequency = minimizer_frequency[it->first];
            
            // Remove minimizers that occur in >75% of organisms (too common to be informative)
            if (frequency > organisms.size() * 0.75) {
                it = minimizer_index.erase(it);
                removed_common++;
            } else {
                ++it;
            }
        }
        
        std::cout << "Removed " << removed_common << " overly common minimizers" << std::endl;
        std::cout << "Final index contains " << minimizer_index.size() << " unique minimizers" << std::endl;
        
        // Calculate average minimizers per organism
        size_t total_entries = 0;
        for (const auto& [hash, entries] : minimizer_index) {
            total_entries += entries.size();
        }
        
        std::cout << "Average minimizers per organism: " 
                  << std::fixed << std::setprecision(1) 
                  << (double)total_entries / organisms.size() << std::endl;
    }
    
    std::string build_taxonomy_path(const std::string& organism_name) {
        std::string taxonomy = "Bacteria";
        
        // Basic taxonomy assignment based on organism name patterns
        if (organism_name.find("Escherichia") != std::string::npos) {
            taxonomy = "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia";
        } else if (organism_name.find("Staphylococcus") != std::string::npos) {
            taxonomy = "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus";
        } else if (organism_name.find("Pseudomonas") != std::string::npos) {
            taxonomy = "Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas";
        } else if (organism_name.find("Klebsiella") != std::string::npos) {
            taxonomy = "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella";
        } else if (organism_name.find("Enterococcus") != std::string::npos) {
            taxonomy = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus";
        } else if (organism_name.find("Streptococcus") != std::string::npos) {
            taxonomy = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus";
        } else if (organism_name.find("Clostridioides") != std::string::npos || organism_name.find("Clostridium") != std::string::npos) {
            taxonomy = "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridioides";
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
    
    uint32_t determine_taxonomic_level(const std::string& organism_name) {
        if (organism_name.find("strain") != std::string::npos ||
            organism_name.find("str.") != std::string::npos) {
            return 0;  // Strain level
        }
        
        int word_count = 1;
        for (char c : organism_name) {
            if (c == ' ') word_count++;
        }
        
        return (word_count >= 2) ? 1 : 2;  // Species or genus level
    }
    
public:
    void write_minimizer_database(const std::string& output_path) {
        std::cout << "Writing minimizer database to: " << output_path << std::endl;
        
        std::ofstream out(output_path, std::ios::binary);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot create output file: " + output_path);
        }
        
        // Write header
        uint32_t magic = 0x4D494E49;  // "MINI" - minimizer database format
        uint32_t version = 1;
        uint32_t k_size = k;
        uint32_t m_size = m;
        uint32_t num_organisms = organisms.size();
        uint64_t num_minimizer_hashes = minimizer_index.size();
        
        out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
        out.write(reinterpret_cast<const char*>(&version), sizeof(version));
        out.write(reinterpret_cast<const char*>(&k_size), sizeof(k_size));
        out.write(reinterpret_cast<const char*>(&m_size), sizeof(m_size));
        out.write(reinterpret_cast<const char*>(&num_organisms), sizeof(num_organisms));
        out.write(reinterpret_cast<const char*>(&num_minimizer_hashes), sizeof(num_minimizer_hashes));
        
        // Write organism metadata
        for (const auto& org : organisms) {
            out.write(reinterpret_cast<const char*>(&org.taxonomy_id), sizeof(org.taxonomy_id));
            out.write(reinterpret_cast<const char*>(&org.taxon_level), sizeof(org.taxon_level));
            out.write(reinterpret_cast<const char*>(&org.gc_content), sizeof(org.gc_content));
            out.write(reinterpret_cast<const char*>(&org.genome_size), sizeof(org.genome_size));
            out.write(reinterpret_cast<const char*>(&org.minimizer_count), sizeof(org.minimizer_count));
            
            uint16_t name_length = org.name.length();
            out.write(reinterpret_cast<const char*>(&name_length), sizeof(name_length));
            out.write(org.name.c_str(), name_length);
            
            uint16_t taxonomy_length = org.taxonomy_path.length();
            out.write(reinterpret_cast<const char*>(&taxonomy_length), sizeof(taxonomy_length));
            out.write(org.taxonomy_path.c_str(), taxonomy_length);
        }
        
        // Write minimizer index (sorted for binary search)
        std::vector<std::pair<uint64_t, std::vector<MinimizerEntry>>> sorted_index(
            minimizer_index.begin(), minimizer_index.end());
        
        std::sort(sorted_index.begin(), sorted_index.end());
        
        for (const auto& [hash, entries] : sorted_index) {
            out.write(reinterpret_cast<const char*>(&hash), sizeof(hash));
            uint32_t num_entries = entries.size();
            out.write(reinterpret_cast<const char*>(&num_entries), sizeof(num_entries));
            
            for (const auto& entry : entries) {
                out.write(reinterpret_cast<const char*>(&entry), sizeof(entry));
            }
        }
        
        out.close();
        
        // Calculate database statistics
        size_t file_size = std::filesystem::file_size(output_path);
        size_t total_minimizer_entries = 0;
        for (const auto& [hash, entries] : minimizer_index) {
            total_minimizer_entries += entries.size();
        }
        
        std::cout << "Database written successfully!" << std::endl;
        std::cout << "File size: " << file_size / (1024*1024) << " MB" << std::endl;
        std::cout << "Total minimizer entries: " << total_minimizer_entries << std::endl;
        std::cout << "Bytes per entry: " << std::fixed << std::setprecision(1) 
                  << (double)file_size / total_minimizer_entries << std::endl;
        
        write_database_summary(output_path + ".summary");
    }
    
    void write_database_summary(const std::string& summary_path) {
        std::ofstream summary(summary_path);
        
        summary << "Minimizer Database Summary\n";
        summary << "=========================\n\n";
        
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        summary << "Generated: " << std::ctime(&time_t) << "\n";
        
        summary << "Parameters:\n";
        summary << "- K-mer size: " << k << "\n";
        summary << "- Minimizer size: " << m << "\n";
        summary << "- Window size: " << window_size << "\n";
        summary << "- Stride: " << stride << "\n";
        summary << "- Incremental mode: " << (incremental_mode ? "yes" : "no") << "\n\n";
        
        summary << "Database contents:\n";
        summary << "- Total organisms: " << organisms.size() << "\n";
        summary << "- Unique minimizers: " << minimizer_index.size() << "\n";
        
        size_t total_entries = 0;
        for (const auto& [hash, entries] : minimizer_index) {
            total_entries += entries.size();
        }
        summary << "- Total minimizer entries: " << total_entries << "\n";
        summary << "- Average entries per minimizer: " 
                << std::fixed << std::setprecision(1) 
                << (double)total_entries / minimizer_index.size() << "\n";
        summary << "- Average minimizers per organism: " 
                << (double)total_entries / organisms.size() << "\n\n";
        
        // Genome size distribution
        std::vector<uint64_t> genome_sizes;
        for (const auto& org : organisms) {
            genome_sizes.push_back(org.genome_size);
        }
        std::sort(genome_sizes.begin(), genome_sizes.end());
        
        summary << "Genome size distribution:\n";
        summary << "- Smallest: " << genome_sizes.front() / 1000 << " kb\n";
        summary << "- Largest: " << genome_sizes.back() / 1000000 << " Mb\n";
        summary << "- Median: " << genome_sizes[genome_sizes.size()/2] / 1000000 << " Mb\n\n";
        
        // Taxonomic distribution
        std::unordered_map<uint32_t, int> level_counts;
        for (const auto& org : organisms) {
            level_counts[org.taxon_level]++;
        }
        
        summary << "Taxonomic level distribution:\n";
        std::vector<std::string> level_names = {"Strain", "Species", "Genus", "Family", "Order"};
        for (const auto& [level, count] : level_counts) {
            if (level < level_names.size()) {
                summary << "- " << level_names[level] << ": " << count << "\n";
            }
        }
        summary << "\n";
        
        // Top organisms by minimizer count
        std::vector<OrganismInfo> sorted_orgs = organisms;
        std::sort(sorted_orgs.begin(), sorted_orgs.end(),
                 [](const OrganismInfo& a, const OrganismInfo& b) {
                     return a.minimizer_count > b.minimizer_count;
                 });
        
        summary << "Top organisms by minimizer representation:\n";
        for (int i = 0; i < std::min(15, (int)sorted_orgs.size()); i++) {
            const auto& org = sorted_orgs[i];
            summary << "- " << org.name << ": " << org.minimizer_count << " minimizers "
                    << "(" << std::fixed << std::setprecision(1) << org.gc_content << "% GC, "
                    << org.genome_size / 1000000 << " Mb)\n";
        }
        
        summary.close();
    }
    
    // Add organism from individual FASTA file (for dynamic updates)
    bool add_organism_from_fasta(const std::string& fasta_file) {
        if (processed_files.count(fasta_file)) {
            std::cout << "File already processed: " << fasta_file << std::endl;
            return false;
        }
        
        std::cout << "Adding organism from: " << fasta_file << std::endl;
        
        size_t initial_organism_count = organisms.size();
        load_fasta_file(fasta_file);
        processed_files.insert(fasta_file);
        
        if (organisms.size() > initial_organism_count) {
            // Recalculate uniqueness scores with new data
            calculate_uniqueness_scores();
            build_optimized_index();
            
            std::cout << "Successfully added " << (organisms.size() - initial_organism_count) 
                      << " new organisms" << std::endl;
            return true;
        }
        
        return false;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input_fasta_directory> <output_database> [options]\n";
        std::cout << "\nOptions:\n";
        std::cout << "  --incremental     Update existing database instead of rebuilding\n";
        std::cout << "  --k-size N        K-mer size (default: 35)\n";
        std::cout << "  --minimizer-size N Minimizer size (default: 31)\n";
        std::cout << "\nMinimizer-based database builder for microbial profiling\n";
        std::cout << "Designed for dynamic updates and easy expansion\n";
        std::cout << "\nExamples:\n";
        std::cout << "  # Build new database\n";
        std::cout << "  " << argv[0] << " /data/genomes /data/microbes.db\n";
        std::cout << "\n  # Update existing database with new genomes\n";
        std::cout << "  " << argv[0] << " /data/new_genomes /data/microbes.db --incremental\n";
        return 1;
    }
    
    std::string input_dir = argv[1];
    std::string output_db = argv[2];
    
    // Parse options
    bool incremental = false;
    int k_size = 35;
    int m_size = 31;
    
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--incremental") {
            incremental = true;
        } else if (arg == "--k-size" && i + 1 < argc) {
            k_size = std::stoi(argv[++i]);
        } else if (arg == "--minimizer-size" && i + 1 < argc) {
            m_size = std::stoi(argv[++i]);
        }
    }
    
    try {
        std::cout << "Minimizer Database Builder\n";
        std::cout << "=========================\n";
        
        MinimizerDatabaseBuilder builder(k_size, m_size, incremental);
        
        // Load existing database if in incremental mode
        if (incremental) {
            builder.load_existing_database(output_db);
        }
        
        // Process FASTA directory
        builder.load_fasta_directory(input_dir);
        
        // Write updated database
        builder.write_minimizer_database(output_db);
        
        std::cout << "\nDatabase preparation complete!\n";
        std::cout << "Database file: " << output_db << "\n";
        std::cout << "Summary file: " << output_db << ".summary\n";
        
        if (incremental) {
            std::cout << "\nTo add more organisms later, use:\n";
            std::cout << "  " << argv[0] << " <new_fasta_dir> " << output_db << " --incremental\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}